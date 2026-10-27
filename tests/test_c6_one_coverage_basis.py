"""Constraint C6, one coverage basis (profile handoff, section 4): with core_completeness on
(the signal_anchor module, B3), every coverage test reads the core of the home protein, the
protein without its signal peptide and GPI tail, as the judge does (the larger of
QueryCoverage and CoreCoverage). A signal peptide must never count against a call in one
place and not in another.

Only the judge reads the core today. The fix is fix-phase item 11 of the profile work, not
stage 4: each test below is a strict xfail and flips when its place is fixed (then remove its
marker). Each comes with a control that passes today, showing that the setup reaches the
coverage test at all, so the xfail fails for the reason it names. Places:

* ModelStatus (iterative_search_runner._goi_feature_attrs) reads QueryCoverage only;
* the hull rescue's cost filter (CDS aa / query length) and the coverage its judge reads
  (the query span of the alignment), rescue_goi_hull.py;
* the report's coverage demotion (apply_coverage_demotion) and its identity_coverage_decoupled
  flag (build_self_consistency): generate_report.py never reads CoreCoverage.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402

C6 = pytest.mark.xfail(strict=True, reason="fix-phase item 11 (C6, one coverage basis): this place "
                                           "does not read the core under core_completeness yet")

# UniProt Q16553 (LY6E): signal peptide 1-20 and a GPI tail of 32 aa (query_features), so its
# core is residues 21-99 (79 aa) of 131.
LY6E = ("MKIFLPVLLAALLGVERASSLMCFSCLNQKSNLYCLKPTICSDQDNYCVTVSASAGIGNLVTFGHSLSKTCSPACPIPEGVNVG"
        "VASMGISCCQSFLCNFSAADGGLRASVTLLGAGLLLSLLPALLRFGP")


def test_the_core_of_ly6e():
    assert isr.query_core(LY6E) == (21, 99)
    assert isr.core_coverage([(21, 80)], LY6E) == pytest.approx(60 / 79)


# ----------------------------------------------------------------------------- ModelStatus
def _status(monkeypatch, core_on):
    monkeypatch.setattr(isr, "CORE_COMPLETENESS", core_on)
    attrs = isr._goi_feature_attrs({"ID": "m"}, evidence_type="exon_annotation", identity=90.0,
                                   exon_count=2, query_cov=0.35, core_cov=0.95)
    return attrs["ModelStatus"], attrs["CoreCoverage"]


def test_model_status_control(monkeypatch):
    assert _status(monkeypatch, False) == ("fragment", "0.950")


@C6
def test_model_status_reads_the_core(monkeypatch):
    assert _status(monkeypatch, True)[0] == "complete"


# ----------------------------------------------------------------------------- the hull rescue
def _hull(tmp_path, monkeypatch, settings):
    """rescue_goi_hull.main on one hull with the window, miniprot and the judge replaced:
    miniprot reports one model of LY6E residues 21-80 (60 aa: 0.46 of the precursor, 0.76 of
    its core). Returns the coverage each judged model was judged on (empty = dropped)."""
    rows = ["##gff-version 3"] + [
        f"chr1\tx\tmRNA\t{1_000_000 + i * 200_000}\t{1_010_000 + i * 200_000}\t.\t+\t.\t"
        f"ID=f{i};SynVoyRole=flanking;Confidence=HIGH" for i in range(4)]
    (tmp_path / "region.gff").write_text("\n".join(rows) + "\n")
    (tmp_path / "q.faa").write_text(f">LY6E\n{LY6E}\n")
    (tmp_path / "eff.json").write_text(json.dumps({"settings": settings}))
    model = ("chr1\tminiprot\tmRNA\t100\t279\t50\t+\t.\tID=MP1;Identity=0.9000;Target=LY6E 21 80\n"
             "chr1\tminiprot\tCDS\t100\t279\t50\t+\t0\tParent=MP1;Identity=0.9000;Target=LY6E 21 80\n")
    monkeypatch.setattr(hull, "_extract_window", lambda *a: "ATGGCC" * 1000)
    monkeypatch.setattr(hull, "_run_miniprot_relaxed", lambda *a: model)
    monkeypatch.setattr(hull, "_genome_length", lambda path: 10_000_000)
    judged = []

    def judge(query_id, query_seq, protein, cds_rows, genome_len, proteome, coverage=None, **kw):
        judged.append(coverage)
        return "LOW", (9,), "Basis=test", 0.0
    monkeypatch.setattr(hull, "judge_hull_model", judge)
    monkeypatch.setattr(sys, "argv", [
        "rescue_goi_hull.py", "--region_gff", str(tmp_path / "region.gff"),
        "--target_genome", str(tmp_path / "t.fna"), "--query", str(tmp_path / "q.faa"),
        "--genome_name", "t", "--output", str(tmp_path / "out.gff"),
        "--cluster_max_gap", "3000000", "--window_pad", "10000", "--max_window", "20000000",
        "--effective_params", str(tmp_path / "eff.json")])
    assert hull.main() == 0
    return judged


def test_hull_control(tmp_path, monkeypatch):
    """The cost filter at 0.5 drops the model; at 0.3 it is judged, on its span of the whole
    precursor."""
    assert _hull(tmp_path, monkeypatch, {"core_completeness": True, "goi_hull_min_coverage": 0.5}) == []
    judged = _hull(tmp_path, monkeypatch, {"core_completeness": True, "goi_hull_min_coverage": 0.3})
    assert judged == [pytest.approx(60 / 131)]


@C6
def test_hull_cost_filter_reads_the_core(tmp_path, monkeypatch):
    assert _hull(tmp_path, monkeypatch, {"core_completeness": True, "goi_hull_min_coverage": 0.5})


@C6
def test_hull_judge_reads_the_core(tmp_path, monkeypatch):
    judged = _hull(tmp_path, monkeypatch, {"core_completeness": True, "goi_hull_min_coverage": 0.3})
    assert judged and judged[0] >= 60 / 79 - 0.01


# ----------------------------------------------------------------------------- the report
def _report(tmp_path, core_on):
    """generate_report.py on one HIGH call at 90 % identity over 0.3 of the precursor and 0.9
    of its core, with core_completeness as given in the run's settings."""
    regions = tmp_path / "staged" / "regions"
    regions.mkdir(parents=True)
    attrs = ("ID=GOI_LY6E|chrA_1000;Name=GOI_LY6E;SynVoy_Parent=GOI_LY6E;SynVoyRole=goi;"
             "Confidence=HIGH;GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=90.0;"
             "ModelStatus=partial;QueryCoverage=0.300;CoreCoverage=0.900")
    (regions / "locus_1__GCA_1.1.fna.gff").write_text(
        "##gff-version 3\n" + "\t".join(["chrA", "exon_annotation", "mRNA", "1000", "1500", "90.0",
                                          "+", ".", attrs]) + "\n")
    eff = tmp_path / "effective_params.json"
    eff.write_text(json.dumps({"modules_mode": "auto", "settings": {"core_completeness": core_on}}))
    out = tmp_path / "report.json"
    proc = subprocess.run([sys.executable, str(ROOT / "bin" / "generate_report.py"),
                           "--results_dir", str(tmp_path / "staged"), "--effective_params", str(eff),
                           "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    report = json.loads(out.read_text())
    record = report["goi_dedup"]["records"][0]
    flagged = [f for f in report["self_consistency"]["flags"] if f["type"] == "identity_coverage_decoupled"]
    return record.get("goi_class"), flagged


def test_report_control(tmp_path):
    goi_class, flagged = _report(tmp_path, False)
    assert goi_class == "identity_coverage_decoupled" and len(flagged) == 1


@C6
def test_report_demotion_reads_the_core(tmp_path):
    assert _report(tmp_path, True)[0] != "identity_coverage_decoupled"


@C6
def test_report_flag_reads_the_core(tmp_path):
    assert _report(tmp_path, True)[1] == []
