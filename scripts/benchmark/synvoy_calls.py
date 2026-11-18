#!/usr/bin/env python3
"""SynVoy's calls of a finished run as the benchmark's calls.tsv (README.md, "Normalized tool
output contract"), read from the report's records (synvoy_report.json, goi_dedup.records).

  counted   records labelled HIGH or MEDIUM whose goi_class is not paralog_not_goi or
            identity_coverage_decoupled (the report's demoted classes), and no answer of the
            genome-wide path (genome_wide_strong_hit, never counted)
  curation  records labelled AMBIGUOUS whose goi_class is not paralog_not_goi ("needs curation")
  status    PRESENT with a counted record, AMBIGUOUS with a curation record, else ABSENT

LOW models (the class ambiguous_goi_family_member) are not calls, and regions are not calls:
neither sets a status. The locus is the first counted record of the species (else the first
curation record) in the report's order, best first; its coordinates are 1-based as in the record.
The strand comes from the run's GFF row of the record's model (plot_inputs_*/ and rescue/), else
the placeholder '-'. ``extra`` names the counts (counted, curation, the genome's LOW rows) and
the top region of the species, for information only. A species is the report's genome id
(report_genome_id), for the region BEDs too, without .fa / .fna (.gz).

Usage: synvoy_calls.py RESULTS_DIR CALLS_TSV
"""
import argparse
import glob
import json
import os
import re
import sys

NOT_COUNTED = frozenset({"paralog_not_goi", "identity_coverage_decoupled"})
COLUMNS = ["species", "accession", "called_status", "locus_chrom", "locus_start", "locus_end",
           "strand", "confidence", "extra"]
REGION_RANK = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}


# The report's genome id (bin/generate_report.py, canonical_genome_id): the per-locus staging
# prefix goes, a GCA_/GCF_ accession in the name is the id, else the known file suffixes are
# stripped. Kept in step with the report by tests/test_benchmark_calls_mapping.py.
_STAGING_PREFIX_RE = re.compile(r"^(locus_[A-Za-z0-9.]+|\d+)__")
_KNOWN_SUFFIXES = (".homology.tsv", ".scores.tsv", ".regions.bed", ".gff3", ".gff", ".faa", ".fna",
                   ".fasta", ".fa", ".m8", ".hull_rescue", ".rescue")


def report_genome_id(name):
    """The report's genome id of a file or genome name (generate_report.canonical_genome_id)."""
    name = _STAGING_PREFIX_RE.sub("", os.path.basename(name or ""), count=1)
    m = re.search(r"(GC[AF]_\d+\.\d+)", name)
    if m:
        return m.group(1)
    changed = True
    while changed:
        changed = False
        for suffix in _KNOWN_SUFFIXES:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                changed = True
    return re.sub(r"(_new_genes|\.regions|\.candidates)$", "", name)


def canon_species(name):
    """The species of a report genome or of a run file (a region BED): its report genome id,
    without .fa / .fna (.gz), so a BED and the report's records name a genome alike."""
    return re.sub(r"\.f(n)?a(\.gz)?$", "", report_genome_id(name))


def load_report(results_dir):
    for cand in (os.path.join(results_dir, "synvoy_report.json"),
                 os.path.join(results_dir, "logs", "GENERATE_REPORT__synvoy_report.json")):
        if os.path.isfile(cand):
            with open(cand) as fh:
                return json.load(fh)
    print(f"WARNING: synvoy_report.json not found under {results_dir}", file=sys.stderr)
    return {}


def _attrs(text):
    out = {}
    for item in text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def strand_index(results_dir):
    """{(model id, chrom, start, end): set of strands} of the GOI mRNA and gene rows of the run's
    GFFs (plot_inputs_*/ and rescue/)."""
    rows = {}
    files = sorted(glob.glob(os.path.join(results_dir, "plot_inputs_*", "*.gff"))
                   + glob.glob(os.path.join(results_dir, "rescue", "**", "*.gff*"), recursive=True))
    for path in files:
        if not os.path.isfile(path):
            continue
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("mRNA", "gene"):
                    continue
                a = _attrs(f[8])
                if a.get("SynVoyRole") == "goi" and "ID" in a:
                    try:
                        key = (a["ID"], f[0], int(f[3]), int(f[4]))
                    except ValueError:
                        continue
                    rows.setdefault(key, set()).add(f[6])
    return rows


def top_regions(results_dir):
    """{species: name} of the best region of each species' region BEDs (any home locus): label
    first, then score, then the name."""
    best = {}
    for bed in sorted(glob.glob(os.path.join(results_dir, "regions", "**", "*.regions.bed"), recursive=True)):
        sp = canon_species(os.path.basename(bed))
        with open(bed) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                m = re.search(r"_C([A-Z]+)_S([0-9.]+)", f[3])
                key = (REGION_RANK.get(m.group(1), 0) if m else 0, float(m.group(2)) if m else 0.0)
                if sp not in best or key > best[sp][0] or (key == best[sp][0] and f[3] < best[sp][1]):
                    best[sp] = (key, f[3])
    return {sp: name for sp, (_, name) in best.items()}


def is_counted(rec):
    return (rec.get("confidence") in ("HIGH", "MEDIUM") and rec.get("goi_class") not in NOT_COUNTED
            and not rec.get("genome_wide_strong_hit"))


def is_curation(rec):
    return (rec.get("confidence") == "AMBIGUOUS" and rec.get("goi_class") != "paralog_not_goi"
            and not rec.get("genome_wide_strong_hit"))


def calls(report, strands=None, regions=None):
    """The rows of calls.tsv (dicts with COLUMNS), one per species, sorted."""
    strands, regions = strands or {}, regions or {}
    records = (report.get("goi_dedup") or {}).get("records") or []
    ann = {canon_species(e.get("genome")): e for e in (report.get("annotations") or {}).get("per_genome") or []}
    species = set(ann) | {canon_species(e.get("genome")) for e in (report.get("regions") or {}).get("per_genome") or []}
    species |= {canon_species(r.get("genome")) for r in records}
    species |= {canon_species(g) for g in (report.get("target_policy") or {})}
    species |= set(regions)
    by_species = {}
    for r in records:
        by_species.setdefault(canon_species(r.get("genome")), []).append(r)
    out = []
    for sp in sorted(s for s in species if s):
        recs = by_species.get(sp, [])
        counted = [r for r in recs if is_counted(r)]
        curation = [r for r in recs if is_curation(r)]
        status = "PRESENT" if counted else "AMBIGUOUS" if curation else "ABSENT"
        rec = (counted or curation or [None])[0]
        if rec is not None:
            s = strands.get((rec.get("mrna_id", ""), rec.get("chrom"), int(rec["start"]), int(rec["end"])), set())
            strand = next(iter(s)) if len(s) == 1 and next(iter(s)) in ("+", "-") else "-"
            locus = (rec.get("chrom", "-"), str(rec["start"]), str(rec["end"]), strand, rec.get("confidence", "-"))
        else:
            locus = ("-", "-", "-", "-", "-")
        low = ((ann.get(sp) or {}).get("goi_confidence_counts") or {}).get("LOW", 0)
        extra = f"counted={len(counted)};curation={len(curation)};low={low}"
        if sp in regions:
            extra += f";region={regions[sp]}"
        out.append(dict(zip(COLUMNS, (sp, "-", status) + locus + (extra,))))
    return out


def write_calls(rows, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results_dir", help="The run's --outdir (holds synvoy_report.json)")
    ap.add_argument("calls_tsv", help="calls.tsv to write")
    args = ap.parse_args(argv)
    rows = calls(load_report(args.results_dir), strand_index(args.results_dir), top_regions(args.results_dir))
    if not rows:
        print("WARNING: no SynVoy data found in the report", file=sys.stderr)
    write_calls(rows, args.calls_tsv)
    print(f"wrote {args.calls_tsv}")


if __name__ == "__main__":
    main()
