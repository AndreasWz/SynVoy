"""Fix-phase item 20: the reciprocal-best paralog check of each home locus reaches the report.

RECIPROCAL_BEST_PARALOG wrote <genome>.paralog_check.tsv, without the locus. GENERATE_REPORT
stages every check into numbered folders (paralog_check/?/*) and copies them into one folder, so
with two or more home loci the checks of all but one locus were overwritten, and which one
survived depended on the order in which the tasks finished. The file is now named like
ASSIGN_LOCUS_OWNERSHIP's (<locus>__<genome>...), and every list the report stages is sorted.
"""
import csv
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402


def _read(rel):
    return open(os.path.join(ROOT, rel)).read()


def test_the_check_file_carries_the_locus():
    nf = _read("modules/reciprocal_best_paralog.nf")
    assert 'path("${locus_id}__${genome_name}.paralog_check.tsv"), emit: tsv' in nf
    assert "--output ${locus_id}__${genome_name}.paralog_check.tsv" in nf
    assert 'path("${genome_name}.paralog_check.tsv")' not in nf
    own = _read("modules/assign_locus_ownership.nf")
    assert 'path("${locus_id}__${genome_name}.locus_ownership.tsv")' in own


def test_every_list_the_report_stages_is_sorted():
    sw = _read("subworkflows/adjudicate_and_report.nf")
    block = sw[sw.index("    STAGE_REGION_GFF("):sw.index("    GENERATE_REPORT(")]
    assert block.count(".collect(sort: true)") == 8
    assert not re.search(r"\.collect\(\)", block)


def _check(path, locus, genome, best):
    cols = ["genome", "locus", "chrom", "start", "end", "target_id", "best_paralog", "best_bit",
            "second_paralog", "second_bit", "bitscore_gap"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerow({"genome": genome, "locus": locus, "chrom": "c1", "start": 100, "end": 400,
                    "target_id": f"m_{locus}", "best_paralog": best, "best_bit": 90,
                    "second_paralog": "p2", "second_bit": 40, "bitscore_gap": 50})


def test_two_loci_of_one_genome_both_reach_the_report(tmp_path):
    """Both files of one genome are read, one row per locus; with the old name the second
    copy into the staging folder replaced the first."""
    staged = tmp_path / "paralog_check"
    staged.mkdir()
    for locus in ("locus_1", "locus_2"):
        _check(staged / f"{locus}__G.fna.paralog_check.tsv", locus, "G.fna", f"p_{locus}")
    rows = gr.load_paralog_check_rows(str(staged))
    assert sorted((r["locus"], r["best_paralog"]) for r in rows) == [("locus_1", "p_locus_1"),
                                                                      ("locus_2", "p_locus_2")]
    old = tmp_path / "old"
    old.mkdir()
    for locus in ("locus_1", "locus_2"):            # what the flat copy of the old names left
        _check(old / "G.fna.paralog_check.tsv", locus, "G.fna", f"p_{locus}")
    assert [r["locus"] for r in gr.load_paralog_check_rows(str(old))] == ["locus_2"]
