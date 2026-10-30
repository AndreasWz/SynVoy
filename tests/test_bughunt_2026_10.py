"""Fixes from the bug hunt of 2026-10-13 that have no better home among the existing tests.

A2  report dedup: a compact model replaces a wider one only when it is itself a counted call
A9  a miniprot timeout skips one rescue window instead of killing the rescue task
A10 a genome the search could not search is recorded and listed by the report
A13 the ownership check also reads tandem copies (lone `gene` rows)
A14 the melittin coordinate scorer: AMBIGUOUS rank, calls on lost-gene loci
and the strong-synteny rescue writing identity as a fraction.
"""
import os
import subprocess
import sys
import types

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "scripts", "benchmark"))

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import reciprocal_best_paralog_check as rbp  # noqa: E402
import rescue_strong_synteny as rss  # noqa: E402
import score_coordinates as sc  # noqa: E402


# --------------------------------------------------------------------------- A2

def _rec(conf, start, end, cov, ident="60.0", source="locus_1"):
    return {"genome": "g1", "source": source, "chrom": "chr1", "start": start, "end": end,
            "confidence": conf, "goi_class": "x", "identity": ident, "query_cov": cov,
            "target_gene": "", "model_status": "partial", "mrna_id": f"m{start}"}


def test_ambiguous_compact_model_does_not_replace_a_high_one():
    wide = _rec("HIGH", 1000, 12000, "0.95")
    compact = _rec("AMBIGUOUS", 3000, 5000, "0.90", ident="40.0")
    recs = gr.dedupe_goi_annotations([wide, compact])["records"]
    assert len(recs) == 1 and recs[0]["confidence"] == "HIGH"
    assert not recs[0]["superseded_wide_model"]


def test_counted_compact_model_still_replaces_a_fused_one():
    """The case the rule was made for (2026-09-15) is unchanged."""
    wide = _rec("HIGH", 1000, 12000, "0.95")
    compact = _rec("MEDIUM", 3000, 5000, "0.90")
    recs = gr.dedupe_goi_annotations([wide, compact])["records"]
    assert len(recs) == 1 and recs[0]["confidence"] == "MEDIUM"
    assert recs[0]["superseded_wide_model"]


# --------------------------------------------------------------------------- A9

def test_miniprot_beyond_its_guard_stops_the_rescue(monkeypatch):
    """A9 (2026-10-13) made a timeout skip the window instead of killing the task, which
    errorStrategy 'ignore' hid; fix-phase item 18 (2026-10-29): the guard stops the rescue with
    StepFailed, the script exits STEP_FAILED_EXIT and the module ends the run."""
    def slow(*_a, **kw):
        raise subprocess.TimeoutExpired(cmd="miniprot", timeout=kw.get("timeout"))
    monkeypatch.setattr(rss.subprocess, "run", slow)
    with pytest.raises(isr.StepFailed, match="did not finish within 3 s"):
        rss._run_miniprot_relaxed("w.fa", "q.fa", 0.05, 3)


def test_missing_miniprot_still_fails_loud(monkeypatch):
    """Since fix-phase item 18 a missing miniprot stops the rescue with StepFailed (exit 3, which
    the module turns into the run's end); the FileNotFoundError before exited 1, which the module
    ignores."""
    def missing(*_a, **_kw):
        raise FileNotFoundError("miniprot")
    monkeypatch.setattr(rss.subprocess, "run", missing)
    with pytest.raises(isr.StepFailed, match="miniprot is not installed here"):
        rss._run_miniprot_relaxed("w.fa", "q.fa", 0.05, 3)


# --------------------------------------------------------------------------- A10

def test_skipped_genome_is_recorded_and_reported(tmp_path):
    out = tmp_path / "search"
    args = types.SimpleNamespace(output_dir=str(out))
    isr._write_not_searched(args, "locus_1_", "Apis_cerana.fna", "mmseqs_resource_failure",
                            "Prefilter died\n  Cannot fit databases into 2G")
    written = out / "hits" / "locus_1_Apis_cerana.fna.not_searched.tsv"
    assert written.exists()

    results = tmp_path / "results"
    (results / "regions").mkdir(parents=True)
    (results / "hits").mkdir()
    row = ("chr1\tx\tmRNA\t100\t400\t60\t+\t.\tID=GOI_q|g2_b0_l1;SynVoy_Parent=GOI_q;"
           "SynVoyRole=goi;Confidence=HIGH;GOIClass=confident_goi;EvidenceType=exon_annotation;"
           "Identity=60;QueryCoverage=0.9")
    (results / "regions" / "locus_1__Apis_florea.fna.gff").write_text("##gff-version 3\n" + row + "\n")
    (results / "hits" / written.name).write_text(written.read_text())
    report = gr.build_report(str(results))
    assert report["summary"]["genomes_not_searched"] == {"Apis_cerana": ["mmseqs_resource_failure"]}
    assert any("could not be searched" in w for w in report["warnings"])
    assert "Apis_cerana" not in report["summary"]["goi_absent_genomes"]


# --------------------------------------------------------------------------- A13

def test_ownership_check_reads_tandem_copies(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text("\n".join([
        "##gff-version 3",
        # an ordinary model: gene + mRNA (the gene row is the mRNA's parent, not a model)
        "chr1\tr\tgene\t100\t900\t.\t+\t.\tID=GOI_q|g_r1;SynVoyRole=goi;Confidence=HIGH",
        "chr1\tr\tmRNA\t100\t900\t.\t+\t.\tID=GOI_q|g_m1;Parent=GOI_q|g_r1;SynVoyRole=goi;Confidence=HIGH",
        # a tandem copy: a lone gene row
        "chr1\ttandem_copy\tgene\t5000\t5200\t55\t+\t.\tID=GOI_copy_1|g_b0_l1;SynVoyRole=goi;Confidence=MEDIUM",
        # a flanking gene: not a GOI model
        "chr1\tf\tgene\t9000\t9900\t.\t+\t.\tID=gene-X|g_fl1;SynVoyRole=flanking;Confidence=HIGH",
    ]) + "\n")
    ids = [m["mrna_id"] for m in rbp._find_goi_models(str(gff))]
    assert ids == ["GOI_q|g_m1", "GOI_copy_1|g_b0_l1"]


# --------------------------------------------------------------------------- A14

def _truth(kind, chrom, species):
    return {"species": species, "chrom": chrom, "chrom_keys": {chrom}, "start": 100, "end": 500,
            "truth_kind": kind, "model_id": f"{species}_{kind}", "confidence": "gold"}


def _call(chrom, conf="MEDIUM"):
    return {"chrom": chrom, "start": 150, "end": 450, "confidence": conf, "identity": "50",
            "goi_class": "probable_goi", "query_cov": "0.8"}


def test_call_on_a_lost_gene_locus_is_not_a_precision_hit():
    truth = [_truth("gene", "chr1", "A"), _truth("locus_gene_lost", "chr2", "B")]
    _models, calls, summary = sc.score(truth, [_call("chr1"), _call("chr2")])
    assert summary["locus_precision"] == 0.5
    assert summary["false_positive_calls_on_lost_loci"] == 1
    assert [c["on_truth_model"] for c in calls] == ["yes", "no"]


def test_ambiguous_calls_are_kept_at_a_low_floor_and_dropped_at_medium():
    truth = [_truth("gene", "chr1", "A")]
    assert sc.score(truth, [_call("chr1", "AMBIGUOUS")], min_confidence="LOW")[2]["calls_total"] == 1
    assert sc.score(truth, [_call("chr1", "AMBIGUOUS")], min_confidence="MEDIUM")[2]["calls_total"] == 0


# ------------------------------------------------------------ rescue identity scale

def test_strong_synteny_rescue_identity_is_percent():
    assert rss._identity_from_mrna_attrs("ID=x;Identity=0.962") == pytest.approx(96.2)
    assert rss._identity_from_mrna_attrs("ID=x;Identity=1") == 100.0
    assert rss._identity_from_mrna_attrs("ID=x") == 0.0
