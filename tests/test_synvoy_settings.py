"""bin/synvoy_settings.py, the read-only loader of the run's settings (profile handoff, stage 4),
and the three scripts that read a rule key through it after the search:

* build_home_paralog_panel.py  locus_ownership_max_genes_per_locus, paralog_panel_max_family
* rescue_goi_hull.py           goi_hull_min_coverage
* plot_synteny.py              max_goi_per_genome (tested through scripts/replot.sh in
                               tests/test_replot_run.py)

A value on the script's own command line wins, then the settings value, then the script's
default (a call without the file: a direct call, a re-draw of an older run).
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import rescue_goi_hull as hull  # noqa: E402
import synvoy_settings as ss  # noqa: E402


def _eff(tmp_path, settings, name="effective_params.json"):
    path = tmp_path / name
    path.write_text(json.dumps({"modules_mode": "auto", "settings": settings, "source": {}}))
    return str(path)


# ----------------------------------------------------------------------------- the loader
def test_load_settings_returns_the_settings_map_read_only(tmp_path):
    settings = ss.load_settings(_eff(tmp_path, {"max_goi_per_genome": 11, "core_completeness": True}))
    assert dict(settings) == {"max_goi_per_genome": 11, "core_completeness": True}
    with pytest.raises(TypeError):
        settings["max_goi_per_genome"] = 3          # read-only
    assert ss.load_settings(_eff(tmp_path, {"a": 1}))["a"] == 1


def test_a_sentinel_or_a_missing_file_is_an_empty_mapping(tmp_path):
    sentinel = tmp_path / "NO_EFFECTIVE_PARAMS"
    sentinel.write_text((ROOT / "assets" / "sentinels" / "NO_EFFECTIVE_PARAMS").read_text())
    empty = tmp_path / "empty.json"
    empty.write_text("")
    for path in (None, "", str(tmp_path / "absent.json"), str(sentinel), str(empty)):
        settings = ss.load_settings(path)
        assert dict(settings) == {}, path
        with pytest.raises(TypeError):
            settings["x"] = 1


def test_a_file_that_is_no_settings_file_stops_the_step(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(ValueError, match="not a JSON file"):
        ss.load_settings(str(bad))
    report = tmp_path / "synvoy_report.json"
    report.write_text(json.dumps({"summary": {}}))
    with pytest.raises(ValueError, match="no settings map"):
        ss.load_settings(str(report))


@pytest.mark.parametrize("value,want", [(8, 8), (8.0, 8), ("8", 8), (" 8.0 ", 8)])
def test_whole_reads_a_count_as_nextflow_hands_it_over(value, want):
    assert ss.whole(value) == want


@pytest.mark.parametrize("value", ["8.5", "auto", True, None])
def test_whole_refuses_what_is_no_count(value):
    with pytest.raises((ValueError, TypeError)):
        ss.whole(value)


def test_pick_option_before_settings_before_default():
    settings = {"k": "7"}
    assert ss.pick(3, settings, "k", 10, ss.whole) == (3, "option")
    assert ss.pick(None, settings, "k", 10, ss.whole) == (7, "settings")
    assert ss.pick(None, {"k": None}, "k", 10, ss.whole) == (10, "default")
    assert ss.pick(None, {}, "k", 10, ss.whole) == (10, "default")
    assert ss.pick(None, {"k": "0.3"}, "k", 0.5, float) == (0.3, "settings")
    with pytest.raises(ValueError, match="settings value of k"):
        ss.pick(None, {"k": "many"}, "k", 10, ss.whole)


def test_the_loader_writes_into_no_namespace(tmp_path):
    import argparse
    args = argparse.Namespace(max_goi_per_genome=None)
    settings = ss.load_settings(_eff(tmp_path, {"max_goi_per_genome": 12}))
    assert ss.pick(args.max_goi_per_genome, settings, "max_goi_per_genome", 10, ss.whole) == (12, "settings")
    assert vars(args) == {"max_goi_per_genome": None}


# ----------------------------------------------------------------------------- the panel
def _panel_inputs(tmp_path, n_genes=8):
    """One home locus over n_genes annotated genes of different lengths (the panel ranks
    them by overlap with the locus), each with a protein in the home proteome."""
    gff, faa = ["##gff-version 3"], []
    for i in range(n_genes):
        s = 10_000 + i * 5_000
        gff.append(f"chr1\tt\tgene\t{s}\t{s + 1000 + i * 100}\t.\t+\t.\tID=gene-g{i}")
        faa.append(f">gene-g{i}\nMKT{'AILVG' * (10 + i)}")
    (tmp_path / "home.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "home_proteome.faa").write_text("\n".join(faa) + "\n")
    (tmp_path / "locus_1.bed").write_text(f"chr1\t9000\t{12_000 + n_genes * 5_000}\tGOI_x|exon_1\t100\t+\n")


def _panel(tmp_path, *extra):
    out = tmp_path / "panel.faa"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "build_home_paralog_panel.py"),
         "--home_gff", str(tmp_path / "home.gff"), "--home_proteome", str(tmp_path / "home_proteome.faa"),
         "--locus_beds", str(tmp_path / "locus_1.bed"), "--output", str(out), *extra],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return sum(1 for line in out.read_text().splitlines() if line.startswith(">locus_1|")), proc.stderr


def test_the_panel_takes_its_genes_per_locus_from_the_run_settings(tmp_path):
    _panel_inputs(tmp_path)
    assert _panel(tmp_path)[0] == 5                                       # no file: the default
    eff = _eff(tmp_path, {"locus_ownership_max_genes_per_locus": 7, "paralog_panel_max_family": 15})
    n, log = _panel(tmp_path, "--effective_params", eff)
    assert n == 7 and "7 genes per locus (settings)" in log and "15 family paralogs (settings)" in log
    assert _panel(tmp_path, "--effective_params", eff, "--max_genes_per_locus", "6")[0] == 6
    sentinel = ROOT / "assets" / "sentinels" / "NO_EFFECTIVE_PARAMS"
    assert _panel(tmp_path, "--effective_params", str(sentinel))[0] == 5


# ----------------------------------------------------------------------------- the hull rescue
QUERY = "MKTAYIAKQR" * 10        # 100 aa, no signal peptide, no GPI tail


def _hull(tmp_path, monkeypatch, *extra, query=QUERY, cds_aa=40, target=(1, 40)):
    """rescue_goi_hull.main on one hull (four HIGH flanking models on chr1), with the window,
    miniprot and the judge replaced: miniprot reports one model of cds_aa residues at 90 %
    identity, aligned to the query residues `target`. Returns the coverage each judged model
    was judged on (empty when the cost filters dropped the model) and the console."""
    rows = ["##gff-version 3"] + [
        f"chr1\tx\tmRNA\t{1_000_000 + i * 200_000}\t{1_010_000 + i * 200_000}\t.\t+\t.\t"
        f"ID=f{i};SynVoyRole=flanking;Confidence=HIGH" for i in range(4)]
    (tmp_path / "region.gff").write_text("\n".join(rows) + "\n")
    (tmp_path / "q.faa").write_text(f">GOI_X\n{query}\n")
    end = 100 + cds_aa * 3 - 1
    tgt = f"Target=GOI_X {target[0]} {target[1]}"
    model = (f"chr1\tminiprot\tmRNA\t100\t{end}\t50\t+\t.\tID=MP1;Identity=0.9000;{tgt}\n"
             f"chr1\tminiprot\tCDS\t100\t{end}\t50\t+\t0\tParent=MP1;Identity=0.9000;{tgt}\n")
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
        "--cluster_max_gap", "3000000", "--window_pad", "10000", "--max_window", "20000000", *extra])
    assert hull.main() == 0
    return judged


def test_the_hull_takes_its_coverage_filter_from_the_run_settings(tmp_path, monkeypatch, capsys):
    """A model over 0.4 of the query: dropped by the default 0.5, judged under B1's 0.3."""
    assert _hull(tmp_path, monkeypatch) == []
    assert "CDS >= 0.5 of the query (default)" in capsys.readouterr().err
    b1 = _eff(tmp_path, {"goi_hull_min_coverage": 0.3})
    assert _hull(tmp_path, monkeypatch, "--effective_params", b1) == [pytest.approx(0.4)]
    assert "CDS >= 0.3 of the query (settings)" in capsys.readouterr().err
    assert _hull(tmp_path, monkeypatch, "--effective_params", b1, "--min_coverage", "0.5") == []
    assert "CDS >= 0.5 of the query (option)" in capsys.readouterr().err
    legacy = _eff(tmp_path, {"goi_hull_min_coverage": "0.5"}, name="legacy.json")
    assert _hull(tmp_path, monkeypatch, "--effective_params", legacy) == []
