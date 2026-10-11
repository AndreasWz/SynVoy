"""A low-quality genome the user chose is searched and flagged, not silently left out.

Before 2026-10-08 a target genome with an N50 below 5 kb or more than 500,000 sequences
was removed before the search in every mode, also when the user had supplied the FASTA
(pro mode) or named the species (``--target_species``). The run then showed no row for
it, which reads as "no gene there". The default policy is now ``auto``: a genome the
pipeline picked on its own is still dropped, a genome the user chose is kept and flagged
in the console, in ``synvoy_report.json`` and in the figures.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import fetch_related_genomes as frg  # noqa: E402
import generate_report as gr  # noqa: E402

PY = sys.executable


def _args(policy):
    return SimpleNamespace(bad_quality_policy=policy, bad_max_contigs=500000,
                           bad_max_scaffolds=500000, bad_min_n50=5000, bad_quality_timeout=0)


def _assemblies():
    good = {"species": "Good one", "accession": "GCF_1", "assembly_status": "Scaffold",
            "contig_n50": 2_000_000, "scaffold_n50": 9_000_000,
            "contig_count": None, "scaffold_count": None}
    bad = {"species": "Bad one", "accession": "GCA_2", "assembly_status": "Scaffold",
           "contig_n50": 900, "scaffold_n50": 1200,
           "contig_count": None, "scaffold_count": None}
    return [good, bad]


# ---------------------------------------------------------------- easy mode: the fetch
@pytest.mark.parametrize("policy,user_named,expected", [
    ("auto", False, "drop"), ("auto", True, "keep"),
    ("drop", True, "drop"), ("keep", False, "keep"), ("ask", True, "ask"),
])
def test_effective_policy(policy, user_named, expected):
    assert frg.effective_bad_quality_policy(policy, user_named) == expected


def test_auto_drops_a_bad_assembly_the_pipeline_picked(capsys):
    kept = frg.apply_bad_quality_policy(_assemblies(), _args("auto"), user_named=False)
    assert [a["accession"] for a in kept] == ["GCF_1"]
    assert "[drop]" in capsys.readouterr().out


def test_auto_keeps_a_bad_assembly_of_a_named_species_and_says_so(capsys):
    kept = frg.apply_bad_quality_policy(_assemblies(), _args("auto"), user_named=True)
    assert [a["accession"] for a in kept] == ["GCF_1", "GCA_2"]
    assert kept[1]["bad_quality"] is True and "N50" in kept[1]["bad_quality_reasons"]
    err = capsys.readouterr().err
    assert "LOW-QUALITY" in err and "Bad one" in err and "--target_species" in err


def test_explicit_drop_still_drops_a_named_species():
    kept = frg.apply_bad_quality_policy(_assemblies(), _args("drop"), user_named=True)
    assert [a["accession"] for a in kept] == ["GCF_1"]


# ---------------------------------------------------------------- all modes: the filter
def _run_filter(tmp_path, policy):
    sorted_list = tmp_path / "sorted.tsv"
    sorted_list.write_text("good.fna\t0.1\nbad.fna\t0.4\n")
    qc = tmp_path / "qc.json"
    qc.write_text(json.dumps([
        {"genome": "good.fna", "status": "PASS", "msg": ""},
        {"genome": "bad.fna", "status": "FAIL", "msg": "N50 (812) below threshold (5000)"},
    ]))
    out = tmp_path / "out.tsv"
    proc = subprocess.run(
        [PY, str(ROOT / "bin" / "filter_sorted_genomes.py"), "--sorted", str(sorted_list),
         "--qc_json", str(qc), "--policy", policy, "--output", str(out)],
        capture_output=True, text=True, check=True)
    return out.read_text().split(), proc.stderr


def test_filter_keeps_a_failed_genome_with_a_warning(tmp_path):
    kept, err = _run_filter(tmp_path, "keep")
    assert "bad.fna" in kept and "good.fna" in kept
    assert "WARNING" in err and "bad.fna" in err and "N50 (812)" in err
    assert "searched anyway" in err


def test_filter_drops_a_failed_genome_under_drop(tmp_path):
    kept, err = _run_filter(tmp_path, "drop")
    assert "bad.fna" not in kept and "good.fna" in kept
    assert "dropped bad.fna" in err and "WARNING" not in err


# ---------------------------------------------------------------- the workflow's choice
def test_workflow_resolves_auto_by_who_chose_the_genomes():
    main_nf = (ROOT / "main.nf").read_text()
    body = re.search(r"def qcFailPolicy\(\) \{(.*?)\n\}", main_nf, re.S).group(1)
    assert "params.mode == 'pro'" in body and "params.target_species" in body
    assert "'keep' : 'drop'" in body
    # the resolved value, not the raw parameter, reaches the filter
    assert re.search(r"FILTER_SORTED_GENOMES\(\s*PHYLO_SORT\.out\.sorted_list,\s*qc_summary_ch,\s*qc_policy\s*\)",
                     main_nf)
    assert "['auto', 'drop', 'keep']" in main_nf
    assert "['auto', 'ask', 'drop', 'keep']" in main_nf


def test_shipped_defaults_are_auto():
    config = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*qc_fail_policy\s*=\s*'auto'", config, re.M)
    assert re.search(r"^\s*bad_quality_policy\s*=\s*'auto'", config, re.M)


@pytest.mark.parametrize("script", ["fetch_related_genomes.py", "fetch_home_genome.py"])
def test_fetch_scripts_accept_auto(script):
    text = (ROOT / "bin" / script).read_text()
    assert re.search(r'"--bad-quality-policy",\s*choices=\["auto", "ask", "drop", "keep"\]', text)


# ---------------------------------------------------------------- the report
def test_report_names_low_quality_genomes_that_were_searched():
    failed = [
        {"genome": "bad_searched", "msg": "N50 (812) below threshold (5000)"},
        {"genome": "bad_dropped", "msg": "N50 (10) below threshold (5000)"},
    ]
    ids, sentences = gr.low_quality_warnings(failed, {"bad_searched", "good"})
    assert ids == ["bad_searched"]
    assert len(sentences) == 1
    assert "LOW-QUALITY ASSEMBLY: bad_searched" in sentences[0]
    assert "N50 (812)" in sentences[0] and "missing from the assembly" in sentences[0]


def test_report_has_no_warning_on_a_clean_run():
    assert gr.low_quality_warnings([], {"a", "b"}) == ([], [])
    assert gr.low_quality_warnings([{"genome": "x", "msg": "m"}], set()) == ([], [])


# ---------------------------------------------------------------- the figures
import plot_synteny as ps  # noqa: E402


@pytest.fixture
def flagged():
    """One genome (g1) marked as failed, for the duration of a test."""
    ps._LOW_QUALITY_GENOMES = {"g1"}
    yield
    ps._LOW_QUALITY_GENOMES = set()


def test_low_quality_genomes_are_read_from_the_quality_summary(tmp_path):
    qc = tmp_path / "genome_qc_summary.json"
    qc.write_text(json.dumps([
        {"genome": "/work/x/Bombus_impatiens.fna", "status": "FAIL", "msg": "N50 (812)"},
        {"genome": "GCF_1.2.fa.gz", "status": "fail", "msg": ""},
        {"genome": "fine.fna", "status": "PASS", "msg": ""},
    ]))
    assert ps.load_low_quality_genomes(str(qc)) == {"Bombus_impatiens", "GCF_1.2"}
    assert ps.load_low_quality_genomes(str(tmp_path / "missing.json")) == set()
    assert ps.load_low_quality_genomes(None) == set()
    (tmp_path / "broken.json").write_text("{not json")
    assert ps.load_low_quality_genomes(str(tmp_path / "broken.json")) == set()


def test_no_mark_without_flagged_genomes():
    assert ps._low_quality_mark({"genome_id": "g1"}) == ""


def test_a_flagged_row_carries_a_dagger_and_the_home_row_never_does(flagged):
    assert ps._low_quality_mark({"genome_id": "g1.fna"}) == " " + ps.LOW_QUALITY_MARK
    assert ps._low_quality_mark({"genome_id": "g0"}) == ""
    assert ps._low_quality_mark({"genome_id": "g1", "is_home": True}) == ""
    labels, _width = ps._figure_row_labels(
        [{"species": "Apis mellifera", "genome_id": "home", "is_home": True},
         {"species": "Bombus one", "genome_id": "g1"},
         {"species": "Bombus two", "genome_id": "g2"}], 11)
    assert labels == ["Apis mellifera", "Bombus one " + ps.LOW_QUALITY_MARK, "Bombus two"]


def test_ribbon_legend_explains_the_dagger_only_when_it_is_drawn():
    def symbols(groups):
        return [entry[0] for _title, entries in groups for entry in entries]
    assert ps.LOW_QUALITY_MARK not in symbols(ps._ribbon_legend_groups("GOI"))
    assert ps.LOW_QUALITY_MARK in symbols(ps._ribbon_legend_groups("GOI", low_quality=True))


def test_print_grid_marks_the_row_and_says_why(flagged):
    from test_anchor_grid_print import _print, _tracks
    svg, legend, source, info = _print(*_tracks())
    assert "Bombus species1 " + ps.LOW_QUALITY_MARK in svg
    assert "Bombus species0 " + ps.LOW_QUALITY_MARK not in svg
    assert info["low_quality_rows"] == ["Bombus species1"]
    assert ps.LOW_QUALITY_MARK in legend and "Bombus species1" in legend
    assert "missing from the assembly" in legend
    # the data table keeps the species name as it is
    assert ps.LOW_QUALITY_MARK not in str(source)
    assert "Bombus species1" in str(source)


def test_print_grid_is_unchanged_without_a_flagged_genome():
    from test_anchor_grid_print import _print, _tracks
    svg, legend, _source, info = _print(*_tracks())
    assert ps.LOW_QUALITY_MARK not in svg and ps.LOW_QUALITY_MARK not in legend
    assert info["low_quality_rows"] == []


def test_plot_module_passes_the_quality_summary():
    module = (ROOT / "modules" / "plot_synteny.nf").read_text()
    assert "path genome_qc" in module and "--genome_qc ${genome_qc}" in module
    assert "genome_qc_summary.json" in module          # kept with the plot inputs
    assert "--genome_qc" in (ROOT / "scripts" / "replot.sh").read_text()
