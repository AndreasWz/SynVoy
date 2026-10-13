#!/usr/bin/env python3
"""Stand-in for ``nextflow run main.nf`` in the tests of the web worker.

It takes the same command line the worker builds, writes console lines, a trace file
and an output folder in the shape of the real pipeline, and ends in the way the first
word of the query's FASTA header asks for:

    ok         succeed
    slow       succeed after FAKE_SLOW_SECONDS (default 30) inside the search step
    fail       a task fails                        -> exit 1
    nogene     the gene is not in the home genome  -> exit 1
    noreport   exit 0 without writing a report
    noparasail a task cannot load the pipeline's Python packages -> exit 1
    hang       never end; ignores TERM, and so does a child process
    crash      die together with the launcher, leaving no exit code
    crashonce  like crash on the first attempt, like ok when resumed
    crashfirst like crash on the first attempt, like ok on a later one (resumed or not)
    leak       succeed, leaving a detached process behind

A header such as ``slow_2`` behaves like ``slow`` (the suffix only makes the query,
and with it the cache key, different).
"""

import json
import os
import signal
import subprocess
import sys
import time


def arg(name, default=None):
    if name in sys.argv:
        return sys.argv[sys.argv.index(name) + 1]
    return default


def say(text):
    sys.stdout.write(text)
    sys.stdout.flush()


def model_line(stem, role, gene, block, confidence, identity, coverage, start):
    """One gene model as the search writes it (the attributes the result page reads)."""
    tag = "l1_exon_ann" if role == "goi" else "fl1_flank_ann"
    return "\t".join([
        "chr%s" % block, "flanking_annotation" if role == "flanking" else "exon_annotation",
        "mRNA", str(start), str(start + 90), "%.1f" % identity, "+", ".",
        "ID=%s|%s_fna_b%s_%s%d;Name=%s;SynVoy_Parent=%s;Identity=%.1f;SynVoyRole=%s;"
        "Confidence=%s;QueryCoverage=%.3f" % (gene, stem, block, tag, start, gene, gene,
                                              identity, role, confidence, coverage)]) + "\n"


# What the stand-in search "places" per genome, by its position in the list: the
# neighbours N0..N5 per search block, and the models of the gene itself.
PLACED = [
    {"blocks": {0: [0, 1, 2, 3, 4, 5]}, "goi": [("HIGH", 80.0, 1.0)]},
    {"blocks": {0: [0, 1, 2, 3]}, "goi": [("MEDIUM", 65.0, 0.8)]},
    {"blocks": {0: [1, 2, 3], 1: [5]}, "goi": [("AMBIGUOUS", 50.0, 0.6)]},
    {"blocks": {0: [0, 1], 1: [1, 4]}, "goi": [("LOW", 21.5, 0.43), ("LOW", 19.0, 0.2)]},
]


def write_outputs(out, stems, query_header):
    os.makedirs(os.path.join(out, "intermediate", "split_loci"), exist_ok=True)
    os.makedirs(os.path.join(out, "intermediate", "initial_db"), exist_ok=True)
    with open(os.path.join(out, "intermediate", "initial_db", "initial_db_locus_1.faa"), "w") as fh:
        for i in range(6):                           # the neighbours that have a protein
            fh.write(">gene-N%d|exon_1\nMKV\n>gene-N%d|exon_2\nLLA\n" % (i, i))
        fh.write(">GOI_chrH_1000|exon_1\nMKF\n")
    os.makedirs(os.path.join(out, "intermediate", "flanking"), exist_ok=True)
    os.makedirs(os.path.join(out, "intermediate", "locate_gene"), exist_ok=True)
    os.makedirs(os.path.join(out, "plot_inputs_synteny_block_locus_1"), exist_ok=True)
    with open(os.path.join(out, "intermediate", "split_loci", "locus_1.bed"), "w") as fh:
        fh.write("chrH\t1000\t2000\tgene_loc\t0.0\t+\t500\n")
    with open(os.path.join(out, "intermediate", "flanking", "synteny_block_locus_1.bed"), "w") as fh:
        for i in range(6):
            fh.write("chrH\t%d\t%d\tgene-N%d\t.\t+\tN%d\n" % (i * 100, i * 100 + 50, i, i))
        fh.write("chrH\t1000\t2000\tGOI_chrH_1000\t.\t+\tGOI_chrH_1000\n")
    with open(os.path.join(out, "intermediate", "locate_gene", "effective_params.json"), "w") as fh:
        json.dump({"preset_applied": "preset_single_copy"}, fh)

    plan = ["HIGH", "MEDIUM", "AMBIGUOUS"]
    classes = {"HIGH": "confident_goi", "MEDIUM": "probable_goi",
               "AMBIGUOUS": "syntenic_candidate_unconfirmed"}
    records = []
    for index, stem in enumerate(stems):
        with open(os.path.join(out, "plot_inputs_synteny_block_locus_1", stem + ".fna.gff"), "w") as fh:
            fh.write("##gff-version 3\n")
            placed = PLACED[index] if index < len(PLACED) else {"blocks": {}, "goi": []}
            for block, numbers in sorted(placed["blocks"].items()):
                for n in numbers:
                    fh.write(model_line(stem, "flanking", "gene-N%d" % n, block, "HIGH", 90.0, 1.0,
                                        1000 + 200 * n))
            for k, (confidence, identity, coverage) in enumerate(placed["goi"]):
                fh.write(model_line(stem, "goi", "GOI_chrH_1000", 0, confidence, identity,
                                    coverage, 5000 + 200 * k))
        if index < len(plan):
            confidence = plan[index]
            records.append({
                "genome": stem, "chrom": "chr%d" % (index + 1), "start": 100 * (index + 1),
                "end": 100 * (index + 1) + 300, "confidence": confidence,
                "identity": "%.1f" % (80 - 15 * index), "query_cov": "%.3f" % (1.0 - 0.2 * index),
                "goi_class": classes[confidence], "original_goi_class": classes[confidence],
                "model_status": "complete", "provenance": ["locus_1"], "owning_locus": "locus_1",
            })
    # The home annotation sits in the same folder in a real run; it must not be offered.
    with open(os.path.join(out, "plot_inputs_synteny_block_locus_1", "Homus_primus.gff"), "w") as fh:
        fh.write("##gff-version 3\n")
    report = {
        "genome_qc": [{"genome": stem + ".fna", "total_len": 1000000, "status": "PASS"}
                      for stem in stems],
        "qc_summary": {"failed_genomes": []},
        "goi_dedup": {"records": records},
        "summary": {
            "headline": "fake",
            "high_confidence_goi": sum(1 for r in records if r["confidence"] == "HIGH"),
            "medium_confidence_goi": sum(1 for r in records if r["confidence"] == "MEDIUM"),
            "low_confidence_goi": 2,
        },
        "query_header": query_header,
    }
    with open(os.path.join(out, "synvoy_report.json"), "w") as fh:
        json.dump(report, fh)
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="10" height="10"/></svg>'
    for name in ("anchor_grid.svg", "synteny_plot_view.svg", "gene_positions.svg"):
        with open(os.path.join(out, "synteny_block_locus_1_" + name), "w") as fh:
            fh.write(svg)
    for name in ("anchor_grid.html", "synteny_plot.html", "tree.html"):
        with open(os.path.join(out, "synteny_block_locus_1_" + name), "w") as fh:
            fh.write("<!doctype html><title>figure</title>" + svg)
    with open(os.path.join(out, "locus_1_tree.nwk"), "w") as fh:
        fh.write("(a,b);\n")
    with open(os.path.join(out, "not_listed.txt"), "w") as fh:
        fh.write("a file of the output folder that the result does not offer\n")


def main():
    query = arg("--query")
    out = arg("--outdir")
    trace = arg("-with-trace")
    work = arg("-work-dir")
    targets_dir = arg("--target_genomes")
    resumed = "-resume" in sys.argv
    with open(query) as fh:
        header = fh.readline().strip().lstrip(">")
    mode = header.split("_")[0]
    stems = sorted(os.path.splitext(name)[0] for name in os.listdir(targets_dir))

    os.makedirs(".nextflow", exist_ok=True)         # what a real run leaves in its start folder
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(work, "attempts.txt"), "a") as fh:
        fh.write("resumed\n" if resumed else "fresh\n")
    with open(os.path.join(work, "attempts.txt")) as fh:
        known = len(fh.read().splitlines())         # 1 = nothing of an earlier attempt is here
    digits = "".join(ch for ch in os.path.basename(trace or "") if ch.isdigit())
    attempt = int(digits) if digits else 1          # the worker numbers the trace files
    say("[PIPELINE] fake pipeline | mode=%s%s\n" % (mode, " (resumed)" if resumed else ""))
    seen = os.stat(query)                           # what Nextflow recognises an input by
    say("[PIPELINE] query file: %d bytes, modified %d\n" % (seen.st_size, seen.st_mtime_ns))
    say("[PIPELINE] attempts the working folder has seen: %d\n" % known)
    say("[PROCESS aa/111111] NORMALIZE_QUERY (normalize_query)\n")

    if mode == "nogene":
        say("ERROR ~ Error executing process > 'SPLIT_LOCI (split)'\n")
        say("  ERROR: Gene not found in home genome. Both BLAST and MMseqs2 returned no hits\n")
        return 1
    if mode == "noparasail":
        say("ERROR ~ Error executing process > 'ADJUDICATE_AND_REPORT:BUILD_HOME_PARALOG_PANEL'\n")
        say("  [ERROR] parasail is not importable, but the home-paralog panel needs it\n")
        return 1
    if mode == "crash" or (mode == "crashonce" and not resumed) \
            or (mode == "crashfirst" and attempt == 1):
        os.kill(os.getppid(), signal.SIGKILL)       # the launcher, so no exit code is written
        os.kill(os.getpid(), signal.SIGKILL)

    say("[OK  ] LOCATE_GENE              1 hit(s) -> chrH:1000-2000")
    say("[OK  ] SPLIT_LOCI               identified 1 locus/loci")
    say("[PROCESS bb/222222] EXTRACT_FLANKING (locus_1.bed)\n")
    say("[PROCESS cc/333333] ITERATIVE_SEARCH (iterative)\n")
    task = os.path.join(work, "cc", "333333abcdef")
    os.makedirs(task, exist_ok=True)
    total = len(stems)
    with open(os.path.join(task, ".command.log"), "w") as fh:
        fh.write("[INFO] === Wave 1/2 [genomes 1-1/%d] (1 genome(s)) ===\n" % total)
        fh.write("[INFO] Wave 1 complete: 1 new gene(s) in 1s.\n")
        fh.write("[INFO] === Wave 2/2 [genomes 2-%d/%d] ===\n" % (total, total))

    if mode == "fail":
        say("ERROR ~ Error executing process > 'ITERATIVE_SEARCH (iterative)'\n")
        say("Caused by: Process `ITERATIVE_SEARCH (iterative)` terminated with an error exit status (1)\n")
        return 1
    if mode == "hang":
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        subprocess.Popen(["bash", "-c", "trap '' TERM; sleep 600"])
        while True:
            time.sleep(1)
    if mode == "slow":
        deadline = time.time() + float(os.environ.get("FAKE_SLOW_SECONDS", "30"))
        while time.time() < deadline:
            time.sleep(0.05)

    with open(os.path.join(task, ".exitcode"), "w") as fh:
        fh.write("0\n")
    say("[OK  ] ITERATIVE_SEARCH         complete for locus_1")
    say("[PROCESS dd/444444] ADJUDICATE_AND_REPORT:GENERATE_REPORT (1)\n")
    if mode != "noreport":
        write_outputs(out, stems, header)
    if trace:
        with open(trace, "w") as fh:
            fh.write("task_id\thash\tnative_id\tname\tstatus\texit\n")
            fh.write("1\taa/111111\t1\tNORMALIZE_QUERY (normalize_query)\tCOMPLETED\t0\n")
            fh.write("2\tcc/333333\t2\tITERATIVE_SEARCH (iterative)\tCOMPLETED\t0\n")
            fh.write("3\tee/555555\t3\tRESCUE_GOI_HULL (locus_1/x)\tFAILED\t1\n")
            fh.write("4\tdd/444444\t4\tGENERATE_REPORT (1)\tCOMPLETED\t0\n")
    if mode == "leak":
        subprocess.Popen(["sleep", "600"], start_new_session=True)
    say("[OK  ] GENERATE_REPORT          analysis report generated")
    say("\n[SUCCESS] completed=4 failed=0 cached=0\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
