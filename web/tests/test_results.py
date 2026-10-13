"""Report -> result page, progress reading, examples, and the points where the web
layer depends on what the pipeline writes."""

import json
import os
import zipfile
from pathlib import Path

import pytest

from conftest import FAKE_PIPELINE, REPO_DIR, WEB_DIR
from synvoy_web import gallery, progress, results, runner


def record(genome, confidence, goi_class, identity="60.0", cov="0.900", **extra):
    data = {"genome": genome, "chrom": "chr1", "start": 100, "end": 400, "confidence": confidence,
            "goi_class": goi_class, "original_goi_class": goi_class, "identity": identity,
            "query_cov": cov, "model_status": "complete", "provenance": ["locus_1"]}
    data.update(extra)
    return data


def report(records, genomes, high=None, medium=None, **summary):
    summary.setdefault("low_confidence_goi", 0)
    if high is not None:
        summary["high_confidence_goi"] = high
    if medium is not None:
        summary["medium_confidence_goi"] = medium
    return {"genome_qc": [{"genome": g + ".fna", "total_len": 5000000} for g in genomes],
            "qc_summary": {"failed_genomes": []},
            "goi_dedup": {"records": records}, "summary": summary}


# ── wording and counts ──────────────────────────────────────────────────────

def test_labels_and_counts():
    data = report([
        record("A_a", "HIGH", "confident_goi", "84.3", "1.000"),
        record("B_b", "HIGH", "synteny_hull_rescue", "79.2"),
        record("B_b", "MEDIUM", "tandem_goi_copy", "56.0", "0.357"),
        record("C_c", "MEDIUM", "probable_goi", "55.9"),
        record("D_d", "AMBIGUOUS", "syntenic_candidate_unconfirmed", "57.6"),
        record("E_e", "HIGH", "paralog_not_goi", "70.0", owning_gene="gene-BGN"),
        record("F_f", "MEDIUM", "identity_coverage_decoupled", "90.0", "0.100"),
        record("G_g", "LOW", "whatever"),                       # not a call the page lists
    ], ["A_a", "B_b", "C_c", "D_d", "E_e", "F_f", "G_g", "H_h"], high=2, medium=2,
        low_confidence_goi=5)
    genomes, summary, check = results.summarise(data)
    by_id = {g["id"]: g for g in genomes}
    assert [by_id[k]["status"] for k in "A_a B_b C_c D_d E_e F_f G_g H_h".split()] == [
        "found", "found", "found", "candidate", "none", "none", "none", "none"]
    assert [c["label"] for c in by_id["B_b"]["calls"]] == ["found", "probable"]
    assert by_id["B_b"]["calls"][0]["by_rescue"] is True
    assert by_id["E_e"]["calls"][0]["label"] == "set_aside"
    assert by_id["E_e"]["calls"][0]["reason"] == "related_gene"
    assert by_id["E_e"]["calls"][0]["related_gene"] == "BGN"
    assert by_id["E_e"]["related_gene_only"] is True
    assert by_id["F_f"]["calls"][0]["reason"] == "short_match"
    assert by_id["G_g"]["calls"] == [] and by_id["H_h"]["species"] == "H h"
    assert summary == {
        "genomes": 8, "found": 3, "clear": 2, "candidate": 1, "none": 4, "not_searched": 0,
        "copies": 4, "copies_clear": 2, "copies_probable": 2, "candidates": 1, "set_aside": 2,
        "weak_leads": 5, "identity_min": 55.9, "identity_max": 84.3}
    assert check["agrees"] is True


def test_disagreement_with_the_report_is_flagged():
    data = report([record("A_a", "HIGH", "confident_goi")], ["A_a"], high=0, medium=0)
    assert results.summarise(data)[2] == {"agrees": False, "report_high": 0, "report_medium": 0}


def test_genome_left_out_for_quality_is_named():
    data = report([], ["A_a", "B_b"], high=0, medium=0)
    data["qc_summary"]["failed_genomes"] = [{"genome": "B_b.fna", "msg": "N50 too low"}]
    genomes, summary, _ = results.summarise(data)
    assert {g["id"]: g["status"] for g in genomes} == {"A_a": "none", "B_b": "not_searched"}
    assert summary["not_searched"] == 1 and summary["none"] == 1
    # Searched anyway (the visitor chose it): it counts, and carries the warning.
    data["summary"]["low_quality_genomes_searched"] = ["B_b"]
    genomes, summary, _ = results.summarise(data)
    assert {g["id"]: (g["status"], g["low_quality"]) for g in genomes}["B_b"] == ("none", True)


def test_names_and_order_come_from_the_request():
    data = report([record("B_b", "HIGH", "confident_goi")], ["A_a", "B_b"], high=1, medium=0)
    genomes, _, _ = results.summarise(
        data, {"B_b": {"species": "Bee bee", "common_name": "a bee"}}, ["B_b", "A_a"])
    assert [(g["id"], g["species"], g["common_name"]) for g in genomes] == [
        ("B_b", "Bee bee", "a bee"), ("A_a", "A a", "")]


def test_odd_values_do_not_break_the_page():
    data = report([record("A_a", "HIGH", "confident_goi", identity="", cov=None),
                   record("A_a", "MEDIUM", "probable_goi", identity="nan"),
                   {"confidence": "HIGH"}, "junk", record("", "HIGH", "confident_goi")], ["A_a"])
    genomes, summary, _ = results.summarise(data)
    assert [c["identity"] for c in genomes[0]["calls"]] == [None, None]
    assert summary["identity_min"] is None and summary["copies"] == 2
    assert results.summarise({})[1]["genomes"] == 0


# ── the bundle ──────────────────────────────────────────────────────────────

@pytest.fixture
def output(tmp_path):
    """An output folder as the pipeline leaves it (made by the stand-in pipeline)."""
    import subprocess
    import sys
    targets = tmp_path / "targets"
    targets.mkdir()
    for stem in ("Alpha_beta", "Gamma_delta", "Epsilon_zeta", "Eta_theta"):
        (targets / (stem + ".fna")).write_text(">c\nACGT\n")
    query = tmp_path / "query.faa"
    query.write_text(">ok\nMKF\n")
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(FAKE_PIPELINE), "--query", str(query), "--outdir", str(out),
                    "-work-dir", str(tmp_path / "work"), "--target_genomes", str(targets)],
                   check=True, cwd=str(tmp_path), capture_output=True)
    return out


def test_bundle_holds_what_the_page_links_to(output, tmp_path):
    bundle = tmp_path / "bundle"
    meta = {"id": "x", "query": {"name": "Gene", "length": 3, "fasta": ">ok\nMKF\n"},
            "home": {"species": "Homus primus"}, "pipeline": {"version": "2.0.0"},
            "timing": {"runtime_s": 60}, "failed_tasks": 2,
            "targets": [{"stem": s, "species": s.replace("_", " ")} for s in
                        ("Eta_theta", "Alpha_beta", "Gamma_delta", "Epsilon_zeta")]}
    result = results.build_bundle(output, bundle, meta)
    assert json.loads((bundle / "result.json").read_text()) == result
    assert result["schema"] == 1 and result["pipeline"]["preset"] == "preset_single_copy"
    assert [g["id"] for g in result["genomes"]][0] == "Eta_theta"          # request order
    assert result["summary"]["failed_tasks"] == 2 and result["summary"]["loci"] == 1
    assert "fasta" not in result["query"]
    locus = result["loci"][0]
    assert (locus["n"], locus["chrom"], locus["start"], locus["neighbours"]) == (1, "chrH", 1000, 6)
    assert set(locus["files"]) >= {"grid", "grid_page", "ribbons", "ribbons_page", "positions",
                                   "tree_page", "tree_newick"}
    listed = results.bundle_files(result)
    on_disk = {"files/" + p.name for p in (bundle / "files").iterdir()}
    assert listed == on_disk                                 # nothing unlisted, nothing missing
    assert "files/not_listed.txt" not in listed
    assert [d["key"] for d in result["downloads"]] == ["all", "report", "gene_models", "query", "table"]
    with zipfile.ZipFile(str(bundle / "files" / "gene_models.zip")) as archive:
        names = archive.namelist()
    assert "locus_1/Alpha_beta.fna.gff" in names
    assert not any("Homus_primus" in n for n in names)       # the home annotation is not ours to give
    assert str(tmp_path) not in (bundle / "result.json").read_text()


META = {"id": "x", "query": {"name": "Gene", "length": 3, "accession": "P12345",
                             "organism": "Homus primus", "fasta": ">ok\nMKF\n"},
        "home": {"species": "Homus primus", "assembly": "GCF_1.1"},
        "pipeline": {"label": "pre-release", "commit": "abc1234", "profile": "auto,low_mem"},
        "timing": {"runtime_s": 60, "finished_at": "2026-10-11T10:00:00Z"},
        "targets": [{"stem": s, "species": s.replace("_", " ")} for s in
                    ("Alpha_beta", "Epsilon_zeta", "Eta_theta", "Gamma_delta")]}


def test_every_genome_says_what_was_placed_there(output, tmp_path):
    """A genome without a call is a result too: the page can say how much of the
    neighbourhood the search found there, and what it saw of the gene itself."""
    result = results.build_bundle(output, tmp_path / "bundle", META)
    by_id = {g["id"]: g for g in result["genomes"]}
    assert [by_id[s]["status"] for s in ("Alpha_beta", "Epsilon_zeta", "Eta_theta", "Gamma_delta")] \
        == ["found", "found", "candidate", "none"]
    # all six neighbours side by side; four; three in one place and one elsewhere
    assert by_id["Alpha_beta"]["neighbours"] == {"found": 6, "together": 6, "places": 1,
                                                 "total": 6, "locus": 1}
    assert by_id["Epsilon_zeta"]["neighbours"]["together"] == 4
    assert (by_id["Eta_theta"]["neighbours"]["found"], by_id["Eta_theta"]["neighbours"]["together"],
            by_id["Eta_theta"]["neighbours"]["places"]) == (4, 3, 2)
    # no call: three different neighbours in two places, never more than two together,
    # and two weak matches of the gene, the better one 21.5 % identical
    nothing = by_id["Gamma_delta"]
    assert nothing["calls"] == [] and nothing["neighbours"] == {
        "found": 3, "together": 2, "places": 2, "total": 6, "locus": 1}
    assert nothing["leads"] == {"n": 2, "identity": 21.5, "coverage": 0.43}
    assert by_id["Alpha_beta"]["leads"] is None                  # only where there is no call


def test_genome_in_which_nothing_was_placed(output, tmp_path):
    (output / "plot_inputs_synteny_block_locus_1" / "Gamma_delta.fna.gff").write_text("##gff-version 3\n")
    result = results.build_bundle(output, tmp_path / "bundle", META)
    nothing = {g["id"]: g for g in result["genomes"]}["Gamma_delta"]
    assert nothing["neighbours"] == {"found": 0, "together": 0, "places": 0, "total": 6, "locus": 1}
    assert nothing["leads"] is None
    # A run folder without these files (an older run) still gives a result.
    for path in (output / "plot_inputs_synteny_block_locus_1").glob("*.fna.gff"):
        path.unlink()
    result = results.build_bundle(output, tmp_path / "bundle", META)
    assert all(g["neighbours"] is None for g in result["genomes"])


def test_neighbours_are_counted_among_those_that_were_searched(output, tmp_path):
    # The neighbourhood list also holds genes without a protein (they cannot be found).
    with open(str(output / "intermediate" / "flanking" / "synteny_block_locus_1.bed"), "a") as fh:
        fh.write("chrH\t9000\t9050\tgene-Mir1\t.\t+\tMir1\n")
    result = results.build_bundle(output, tmp_path / "bundle", META)
    assert result["loci"][0]["neighbours"] == 7
    assert result["genomes"][0]["neighbours"]["total"] == 6
    # Without the start database the list is all there is.
    (output / "intermediate" / "initial_db" / "initial_db_locus_1.faa").unlink()
    result = results.build_bundle(output, tmp_path / "bundle", META)
    assert result["genomes"][0]["neighbours"]["total"] == 7


def test_everything_comes_in_one_archive(output, tmp_path):
    result = results.build_bundle(output, tmp_path / "bundle", META)
    first = result["downloads"][0]
    assert first["key"] == "all" and first["file"] == "files/synvoy_Gene.zip"
    path = tmp_path / "bundle" / first["file"]
    assert first["bytes"] == path.stat().st_size
    with zipfile.ZipFile(str(path)) as archive:
        names = set(archive.namelist())
        read_me = archive.read("README.txt").decode()
        table = archive.read("genome_by_genome.tsv").decode()
    assert {"README.txt", "synvoy_report.json", "query.faa", "genome_by_genome.tsv",
            "figures/synteny_block_locus_1_anchor_grid.svg", "figures/locus_1_tree.nwk",
            "gene_models/locus_1/Alpha_beta.fna.gff"} <= names
    assert not any("Homus_primus" in n or n.endswith(".zip") for n in names)   # no archive inside
    assert "SynVoy result: Gene in 4 genomes" in read_me and "P12345" in read_me
    # With a gene symbol, the search goes by that, as on the page.
    named = results.build_bundle(output, tmp_path / "named", dict(
        META, query=dict(META["query"], gene="GEN1", name="A long protein name")))
    assert named["downloads"][0]["file"] == "files/synvoy_GEN1.zip"
    assert "SynVoy pre-release, commit abc1234" in read_me and "2026-10-11 10:00:00 UTC" in read_me
    assert "the gene found in 2, a candidate in 1, no call in 1" in read_me
    assert str(tmp_path) not in read_me + table
    # the table: a header, one row per call, one row for the genome without a call
    rows = [line.split("\t") for line in table.splitlines()]
    assert tuple(rows[0]) == results.TABLE_COLUMNS and len(rows) == 5
    none = dict(zip(rows[0], rows[4]))
    assert (none["species"], none["result"], none["neighbours_side_by_side"],
            none["neighbours_anywhere"], none["neighbours_searched"], none["weak_matches"]) == \
        ("Gamma delta", "not found", "2", "3", "6", "2")
    found = dict(zip(rows[0], rows[1]))
    assert (found["result"], found["confidence"], found["identity_pct"],
            found["query_coverage_pct"]) == ("found", "HIGH", "80.0", "100")


def test_archive_name_is_made_from_the_query_name_safely(output, tmp_path):
    for name, expected in (("../../etc/pass wd", "synvoy_etc_pass_wd.zip"), ("", "synvoy_result.zip"),
                           ("Dpp-4 (α)", "synvoy_Dpp-4.zip"), ("x" * 80, "synvoy_" + "x" * 40 + ".zip")):
        result = results.build_bundle(output, tmp_path / "bundle", dict(META, query={"name": name}))
        assert result["downloads"][0]["file"] == "files/" + expected
        assert (tmp_path / "bundle" / "files" / expected).is_file()


def test_bundle_is_the_same_every_time(output, tmp_path):
    first = results.build_bundle(output, tmp_path / "a", {"id": "x"})
    second = results.build_bundle(output, tmp_path / "b", {"id": "x"})
    assert first == second
    for name in ("gene_models.zip", "synvoy_result.zip", "genome_by_genome.tsv"):
        assert (tmp_path / "a" / "files" / name).read_bytes() == \
            (tmp_path / "b" / "files" / name).read_bytes()
    # Building again over an existing bundle replaces it completely.
    (tmp_path / "a" / "files" / "stale.txt").write_text("x")
    results.build_bundle(output, tmp_path / "a", {"id": "x"})
    assert not (tmp_path / "a" / "files" / "stale.txt").exists()
    assert not (tmp_path / "a.building").exists()


def test_stand_in_tree_is_not_offered(output, tmp_path):
    (output / "locus_1_tree.nwk").write_text("(GOI_placeholder:0.0);\n")
    result = results.build_bundle(output, tmp_path / "bundle", {})
    assert "tree_newick" not in result["loci"][0]["files"]
    assert not (tmp_path / "bundle" / "files" / "locus_1_tree.nwk").exists()


def test_no_bundle_without_a_report(output, tmp_path):
    (output / "synvoy_report.json").write_text("{broken")
    with pytest.raises(results.ResultError):
        results.build_bundle(output, tmp_path / "bundle", {})
    assert not (tmp_path / "bundle").exists()
    (output / "synvoy_report.json").write_text(json.dumps({"summary": {}}))
    with pytest.raises(results.ResultError):
        results.build_bundle(output, tmp_path / "bundle", {})


def test_files_outside_the_output_folder_are_not_copied(output, tmp_path):
    secret = tmp_path / "secret.svg"
    secret.write_text("<svg>secret</svg>")
    figure = output / "synteny_block_locus_1_anchor_grid.svg"
    figure.unlink()
    figure.symlink_to(secret)
    result = results.build_bundle(output, tmp_path / "bundle", {})
    assert "grid" not in result["loci"][0]["files"]
    assert not (tmp_path / "bundle" / "files" / figure.name).exists()


REAL_RUNS = {
    # run folder (kept out of git) -> genomes with the gene, with a candidate, without a call
    "local_runs/robust_yeast": (3, 0, 0),
    "local_runs/robust_new": (3, 0, 1),
}


@pytest.mark.parametrize("run, expected", sorted(REAL_RUNS.items()))
def test_counts_match_real_reports(run, expected):
    path = REPO_DIR / run / "synvoy_report.json"
    if not path.is_file():
        pytest.skip("%s is not on this machine" % run)
    _, summary, check = results.summarise(json.loads(path.read_text()))
    assert (summary["found"], summary["candidate"], summary["none"]) == expected
    assert check["agrees"] is True


# ── progress ────────────────────────────────────────────────────────────────

CONSOLE = (
    "[PIPELINE] SynVoy 2.0.0 | profile=auto,low_mem\n"
    "[SKIP] BORROW_ANNOT             home GFF found (x.gff)[OK  ] STAGE_GENOMES            staged 3 "
    "target genomes[PROCESS d5/6aa3b2] STAGE_GENOMES (staging)\n"
    "[PROCESS d6/236768] LOCATE_GENE (x.fna)\n"
    "[OK  ] LOCATE_GENE              1 hit(s) → chr1:82579-83872[PROCESS e3/99208d] SPLIT_LOCI (split)\n"
)
AFTER_SPLIT = (
    "[OK  ] SPLIT_LOCI               identified 2 locus/loci[PROCESS 34/56e569] "
    "ADJUDICATE_AND_REPORT:BUILD_HOME_PARALOG_PANEL (home_paralog_panel)\n"
    "[PROCESS cf/818ed6] EXTRACT_FLANKING (locus_1.bed)\n"
)
SEARCHING = (
    "[OK  ] QC_FILTER                3 target genome(s) passed QC filter for locus_1"
    "[PROCESS a3/8548f8] ITERATIVE_SEARCH (iterative)\n"
)


def states(text, **kwargs):
    made = progress.build_progress(text, **kwargs)
    return {s["key"]: s for s in made["steps"]}, made


def test_steps_follow_the_console():
    steps, _ = states("")
    assert [steps[k]["state"] for k in progress.STEP_KEYS] == [
        "active", "waiting", "waiting", "waiting", "waiting"]
    steps, _ = states(CONSOLE)
    assert steps["locate"]["state"] == "active"
    steps, _ = states(CONSOLE + AFTER_SPLIT)
    assert steps["locate"] == {"key": "locate", "state": "done", "detail": {"loci": 2}}
    assert steps["neighbours"]["state"] == "active"
    steps, _ = states(CONSOLE + AFTER_SPLIT + SEARCHING)
    assert steps["neighbours"]["state"] == "done" and steps["search"]["state"] == "active"
    assert steps["search"]["detail"] == {"loci_done": 0, "loci_total": 2}
    one = CONSOLE + AFTER_SPLIT + SEARCHING + "[OK  ] ITERATIVE_SEARCH         complete for locus_1"
    steps, _ = states(one)
    assert steps["search"]["state"] == "active" and steps["search"]["detail"]["loci_done"] == 1
    two = one + "[OK  ] ITERATIVE_SEARCH         complete for locus_2"
    steps, _ = states(two)
    assert steps["search"]["state"] == "done" and steps["check"]["state"] == "active"
    three = two + "[PROCESS c4/75507a] ADJUDICATE_AND_REPORT:GENERATE_REPORT (1)\n"
    steps, _ = states(three)
    assert steps["check"]["state"] == "done" and steps["report"]["state"] == "active"
    steps, _ = states(three + "\n[SUCCESS] completed=55 failed=0 cached=0\n")
    assert all(s["state"] == "done" for s in steps.values())


def test_steps_do_not_depend_on_how_nextflow_prints_tasks():
    """Only the pipeline's own stage lines, as an older or newer Nextflow would leave them."""
    stages = ("[OK  ] SPLIT_LOCI               identified 1 locus/loci"
              "[OK  ] QC_FILTER                3 target genome(s) passed QC filter for locus_1")
    steps, _ = states(stages)
    assert steps["neighbours"]["state"] == "done" and steps["search"]["state"] == "active"
    steps, _ = states(stages + "[OK  ] ITERATIVE_SEARCH         complete for locus_1")
    assert steps["search"]["state"] == "done" and steps["check"]["state"] == "active"
    steps, _ = states(stages + "[OK  ] GENERATE_REPORT          analysis report generated")
    assert steps["check"]["state"] == "done"
    old = "[a3/8548f8] Submitted process > ITERATIVE_SEARCH (iterative)\n"
    assert progress.parse_console(old)["submitted"] == [("a3", "8548f8", "ITERATIVE_SEARCH")]


def test_colour_codes_and_noise_are_ignored():
    coloured = "\x1b[32m[OK  ]\x1b[0m \x1b[37mSPLIT_LOCI              \x1b[0m identified 3 locus/loci"
    assert progress.parse_console(coloured)["loci"] == 3
    garbage = "\x00\xff[PROCESS zz/zzzzzz] [OK  ] [[[ ]]] identified many locus/loci" * 50
    assert progress.build_progress(garbage)["steps"][0]["state"] == "active"


def test_genome_counts_from_the_running_search(tmp_path):
    task = tmp_path / "a3" / "8548f8cafe"
    task.mkdir(parents=True)
    log = task / ".command.log"
    log.write_text("2026 [INFO] ═══ Wave 1/3 [genomes 1-1/10] (1 genome(s)) ═══\n")
    text = CONSOLE + AFTER_SPLIT + SEARCHING
    assert states(text, work_dir=tmp_path)[0]["search"]["detail"]["genomes_done"] == 0
    with open(str(log), "a") as fh:
        fh.write("[INFO] Wave 1 complete: 1 new gene(s) in 56s.\n"
                 "[INFO] === Wave 2/3 [genomes 2-5/10] (4 genome(s)) ===\n"
                 "[INFO] Wave 2 complete: no new genes (10s).\n"
                 "[INFO] === Wave 3/3 [genomes 6-10/10] ===\n")
    detail = states(text, work_dir=tmp_path)[0]["search"]["detail"]
    assert (detail["genomes_done"], detail["genomes_total"]) == (5, 10)
    (task / ".exitcode").write_text("0")                     # that task is over
    assert "genomes_done" not in states(text, work_dir=tmp_path)[0]["search"]["detail"]
    assert progress.parse_search_log("nothing useful") is None


def test_trace_counts_finished_and_failed_tasks(tmp_path):
    first, second = tmp_path / "trace_1.tsv", tmp_path / "trace_2.tsv"
    first.write_text("task_id\thash\tname\tstatus\texit\n1\tab/cd\tA\tCOMPLETED\t0\n"
                     "2\tab/ce\tB\tFAILED\t1\n3\tab/cf\tC\n")
    second.write_text("task_id\tname\tstatus\n4\tA\tCACHED\n5\tD\tCOMPLETED\n")
    (tmp_path / "trace_3.tsv").write_text("no header here\n")
    assert progress.parse_trace([first, second, tmp_path / "trace_3.tsv",
                                 tmp_path / "missing.tsv"]) == (3, 1)


def test_only_the_latest_attempt_is_read():
    text = ("\n===== attempt 1, t =====\ncmd\n[OK  ] SPLIT_LOCI   identified 1 locus/loci"
            "[OK  ] ITERATIVE_SEARCH   complete for locus_1"
            "\n===== attempt 2, t =====\ncmd -resume\n[PIPELINE] again\n")
    steps, _ = states(runner.last_attempt(text))
    assert steps["locate"]["state"] == "active"
    assert runner.last_attempt("no marker") == "no marker"


# ── examples ────────────────────────────────────────────────────────────────

def test_gallery_is_built_from_a_spec(output, tmp_path):
    site = tmp_path / "site"
    console = tmp_path / "run.log"
    console.write_text("...Run Summary...Duration:          1h 3m 56sTasks Completed:   55\n")
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"examples": [
        {"slug": "demo", "run_dir": "out", "name": "Demo gene", "about": "a test",
         "scope": "test animals", "console_log": "run.log",
         "home": {"species": "Homus primus"}, "label": "1.0.0", "commit": "abc1234"},
        {"slug": "second", "run_dir": str(output), "name": "Second", "runtime_s": 61}]}))
    stale = site / "data" / "examples" / "old"
    stale.mkdir(parents=True)
    messages = []
    cards = gallery.build_gallery(spec, site, log=messages.append)
    assert [c["slug"] for c in cards] == ["demo", "second"]
    assert cards[0]["runtime_s"] == 3836 and cards[1]["runtime_s"] == 61
    assert (cards[0]["genomes"], cards[0]["found"], cards[0]["candidate"], cards[0]["none"]) == (4, 2, 1, 1)
    listing = json.loads((site / "data" / "gallery.json").read_text())
    assert listing["examples"] == cards
    result = json.loads((site / "data" / "examples" / "demo" / "result.json").read_text())
    assert result["scope"] == "test animals" and result["query"]["name"] == "Demo gene"
    assert (result["pipeline"]["label"], result["pipeline"]["release"]) == ("1.0.0", True)
    second = json.loads((site / "data" / "examples" / "second" / "result.json").read_text())
    assert (second["pipeline"]["label"], second["pipeline"]["release"]) == ("pre-release", False)
    assert (site / "data" / "examples" / "demo" / result["loci"][0]["files"]["grid"]).is_file()
    assert not stale.exists()                                # no longer in the spec
    assert not any("local path" in m for m in messages)


def test_finished_search_can_become_an_example(output, tmp_path):
    job_id = "Aa0000000000000000000a"
    jobs_dir = tmp_path / "jobs"
    results.build_bundle(output, jobs_dir / job_id / "result",
                         {"id": job_id, "query": {"name": "Gene", "accession": "P01501"}})
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"examples": [
        {"slug": "mine", "job": job_id, "about": "from a search", "scope": "test animals"}]}))
    cards = gallery.build_gallery(spec, tmp_path / "site", log=lambda m: None, jobs_dir=jobs_dir)
    assert cards[0]["slug"] == "mine" and cards[0]["name"] == "Gene" and cards[0]["found"] == 2
    made = tmp_path / "site" / "data" / "examples" / "mine"
    text = (made / "result.json").read_text()
    assert job_id not in text                              # the private link is not published
    assert json.loads(text)["scope"] == "test animals"
    assert (made / "files" / "synvoy_report.json").is_file()
    for bad in ("../../etc", "short", 7):
        spec.write_text(json.dumps({"examples": [{"slug": "mine", "job": bad}]}))
        with pytest.raises(gallery.GalleryError):
            gallery.build_gallery(spec, tmp_path / "site", log=lambda m: None, jobs_dir=jobs_dir)
    spec.write_text(json.dumps({"examples": [{"slug": "mine", "job": "Bb0000000000000000000b"}]}))
    with pytest.raises(gallery.GalleryError):
        gallery.build_gallery(spec, tmp_path / "site", log=lambda m: None, jobs_dir=jobs_dir)


def test_gallery_warns_about_local_paths(output, tmp_path):
    figure = output / "synteny_block_locus_1_anchor_grid.svg"
    figure.write_text('<svg><!-- /home/someone/runs/secret_project/x.gff --></svg>')
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"examples": [{"slug": "demo", "run_dir": str(output), "name": "D"}]}))
    messages = []
    gallery.build_gallery(spec, tmp_path / "site", log=messages.append)
    assert any("local path" in m and "/home/someone" in m for m in messages)


@pytest.mark.parametrize("examples", [
    [], [{"slug": "Bad Slug", "run_dir": "out", "name": "x"}],
    [{"slug": "a", "run_dir": "out", "name": "x"}, {"slug": "a", "run_dir": "out", "name": "y"}],
    [{"slug": "a", "run_dir": "nowhere", "name": "x"}], [{"slug": "a", "run_dir": "out"}],
    [{"slug": "a", "run_dir": "out", "name": "x", "typo": 1}],
    [{"slug": "../a", "run_dir": "out", "name": "x"}],
])
def test_gallery_spec_mistakes_are_refused(output, tmp_path, examples):
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({"examples": examples}))
    with pytest.raises(gallery.GalleryError):
        gallery.build_gallery(spec, tmp_path / "site", log=lambda m: None)


# ── what the web layer expects of the pipeline ──────────────────────────────

def read(*parts):
    return (REPO_DIR.joinpath(*parts)).read_text(errors="replace")


def test_pipeline_still_writes_the_progress_lines():
    main = read("main.nf")
    for marker in progress.MAIN_NF_MARKERS:
        assert marker in main, "main.nf no longer contains %r (progress.py reads it)" % marker
    search = read("bin", "iterative_search_runner.py")
    for marker in progress.SEARCH_LOG_MARKERS:
        assert marker in search, "the search log no longer contains %r" % marker
    assert "[genomes {cumulative_genome_idx+1}-" in search


def test_pipeline_still_takes_the_options_the_worker_passes():
    config = read("nextflow.config")
    for name in ("mode", "query", "home_genome", "home_gff", "home_species", "target_genomes",
                 "target_gffs", "outdir"):
        assert ("\n    %s = " % name) in config, "nextflow.config has no parameter %r" % name
    assert "resolveTargetGenomeFiles" in read("main.nf")    # a folder of genomes is accepted
    assert "low_mem {" in config and "auto {" in config     # the default profile exists


def test_pipeline_still_fails_with_the_messages_that_are_translated():
    assert "Gene not found in home genome" in read("bin", "split_loci.py")
    assert "aa, below the" in read("bin", "normalize_query.py")
    steps = {step for _, step in runner._PROCESS_STEP}
    assert steps <= set(progress.STEP_KEYS) and set(runner._STEP_NAMES) == set(progress.STEP_KEYS)
    modules = "".join(p.read_text(errors="replace") for p in (REPO_DIR / "modules").glob("*.nf"))
    for process in ("ITERATIVE_SEARCH", "EXTRACT_FLANKING", "GENERATE_REPORT", "SPLIT_LOCI",
                    "PLOT_SYNTENY", "CLUSTER_REGIONS"):
        assert ("process %s" % process) in modules, "no process %s in modules/" % process


def test_report_still_has_the_fields_the_result_page_reads():
    source = read("bin", "generate_report.py")
    for field in ('"goi_dedup"', '"records"', "high_confidence_goi", "medium_confidence_goi",
                  "low_confidence_goi", "paralog_not_goi", "identity_coverage_decoupled",
                  "query_cov", "owning_gene", "model_status", "low_quality_assembly",
                  "low_quality_genomes_searched", "failed_qc_genomes_with_downstream_results"):
        assert field in source, "generate_report.py no longer mentions %s" % field
    output_doc = read("docs", "OUTPUT.md")
    for name in ("synvoy_report.json", "anchor_grid", "synteny_plot.html", "gene_positions"):
        assert name in output_doc


def test_pipeline_still_marks_the_gene_models_the_result_page_counts():
    """results.placed_in_genome reads these from the per-genome GFF of a search."""
    search = read("bin", "iterative_search_runner.py")
    for marker in ('attrs["SynVoyRole"] = "flanking"', 'attrs["SynVoyRole"] = "goi"',
                   "'SynVoy_Parent': parent_id", "_b{block_idx}_fl{locus_idx}_flank_ann",
                   "QueryCoverage", "Identity"):
        assert marker in search, "the search no longer writes %r" % marker
    # ... and the genes that were searched as neighbours from the start database.
    start_db = read("modules", "prepare_initial_db.nf")
    assert "intermediate/initial_db" in start_db and "initial_db_${locus_id}.faa" in start_db
    assert "'|exon_' in clean_id" in start_db
    assert "plot_inputs_" in read("modules", "plot_synteny.nf")


def test_launcher_script_is_valid_shell():
    import subprocess
    script = WEB_DIR / "synvoy_web" / "run_job.sh"
    assert os.access(str(script), os.X_OK)
    assert subprocess.run(["bash", "-n", str(script)]).returncode == 0


def test_example_is_published_under_its_own_name_not_under_the_runs_label(output, tmp_path):
    """A run is often started under an internal label. An example can say which words
    must not be published, and under which name the run is shown instead."""
    (output / "intermediate" / "query").mkdir(parents=True, exist_ok=True)
    (output / "intermediate" / "query" / "normalized_query.faa").write_text(">INTERNAL_77\nMKFLVNVALV\n")
    rescue = output / "rescue" / "locus_1"
    rescue.mkdir(parents=True)
    (rescue / "Alpha_beta.hull_rescue.gff").write_text(
        "chr1\tx\tmRNA\t1\t9\t.\t+\t.\tID=r1;Name=INTERNAL_77\n")
    site = tmp_path / "site"
    spec = tmp_path / "spec.json"
    example = {"slug": "demo", "run_dir": str(output), "name": "Demo gene",
               "must_not_contain": ["INTERNAL_77"]}

    def build(**changes):
        spec.write_text(json.dumps({"examples": [dict(example, **changes)]}))
        return gallery.build_gallery(spec, site, log=lambda message: None)

    # The query is published under the example's name whatever the run called it ...
    with pytest.raises(gallery.GalleryError) as refused:
        build()
    # ... but the label is also in the gene models of the rescue step: not published.
    assert "gene_models.zip:locus_1/rescue/Alpha_beta.hull_rescue.gff contains 'INTERNAL_77'" in str(refused.value)
    assert "query.faa" not in str(refused.value)
    assert not (site / "data" / "examples" / "demo").exists()

    build(rename={"INTERNAL_77": "GeneX"})
    files = site / "data" / "examples" / "demo" / "files"
    assert files.joinpath("query.faa").read_text().startswith(">Demo_gene\nMKFLVNVALV")
    with zipfile.ZipFile(str(files / "gene_models.zip")) as archive:
        assert b"Name=GeneX" in archive.read("locus_1/rescue/Alpha_beta.hull_rescue.gff")
    assert gallery.scan_for_words(site / "data", ["INTERNAL_77"]) == []
    assert "INTERNAL_77" in (rescue / "Alpha_beta.hull_rescue.gff").read_text()    # the run is untouched
    assert not list((site / "data").glob(".relabel-*"))                           # nothing left behind

    for bad in ({"INTERNAL_77": "a b"}, {"X": "Gene"}, {"INTERNAL_77": "../x"}, ["INTERNAL_77"]):
        with pytest.raises(gallery.GalleryError):
            build(rename=bad)
    with pytest.raises(gallery.GalleryError):
        build(must_not_contain="INTERNAL_77")
