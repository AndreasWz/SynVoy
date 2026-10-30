"""Fix-phase item 18, the places after Smith-Waterman (b-k): every tool run either finishes or
stops its task loudly, through a guard set by a pipeline parameter; no timeout gives an empty
result any more.

Each test puts a fake tool first on PATH (a script that answers `--version` and otherwise sleeps
or fails) and sets a guard of 1 s: the task stops with a message that names the tool, the place
and the parameter, and no empty result is passed on. Until 2026-10-28:
  b  a GOI miniprot run beyond 120 s gave no model, and in the copy rounds ended them;
  c  a flanking miniprot run beyond 120 s gave a hit-chain model;
  d  ANNOTATE_GOI's home model came from the next chromosome with hits;
  e  BORROW_ANNOTATIONS left that target protein out;
  f  the strong-synteny rescue skipped the block (its GFF comment read like "no model");
  g  the hull rescue skipped the hull window (600 s);
  h  Augustus beyond 600 s (or failing) fell back to Prodigal, while the record named Augustus;
  i  Prodigal beyond 300 s failed PREPARE_HOME_PROTEOME but skipped a region of the flank walk;
  j  an `augustus --species=help` slower than 30 s gave the model 'fly', "not installed here";
  k  the miniprot probes had 5 s, after which every gene got a hit-chain model.
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import annotate_goi_exons as age  # noqa: E402
import borrow_annotations as ba  # noqa: E402
import gene_predictor as gp  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
import rescue_strong_synteny as rss  # noqa: E402
import sequence_utils as su  # noqa: E402

QUERY = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
DNA = "ATGAAATTTCTGGTTAATGTTGCTCTGGTTTTTATGGTTGTTTATATTTCTTATATTTATGCTGCT" * 40


def _fake(tmp_path, name, body="sleep 30", version="0.18-r281"):
    """A fake tool on a PATH of its own: `--version` answers at once, anything else runs ``body``."""
    d = tmp_path / "fakebin"
    d.mkdir(exist_ok=True)
    tool = d / name
    tool.write_text(f"#!/bin/bash\nif [ \"$1\" = \"--version\" ]; then echo {version}; exit 0; fi\n{body}\n")
    tool.chmod(0o755)
    return str(d)


@pytest.fixture
def on_path(monkeypatch):
    def put(folder):
        monkeypatch.setenv("PATH", folder + os.pathsep + os.environ.get("PATH", ""))
    return put


def _quick(fn, *a, **k):
    """fn() must stop within a few seconds of a 1-s guard (the fake sleeps 30 s)."""
    t0 = time.monotonic()
    try:
        return fn(*a, **k)
    finally:
        assert time.monotonic() - t0 < 15


# ------------------------------------------------------------------ b, c: the search's miniprot
@pytest.mark.parametrize("sensitive,param", [(True, "goi_miniprot_timeout"),
                                             (False, "flanking_miniprot_timeout")])
def test_b_c_a_miniprot_run_beyond_its_guard_stops_the_task(tmp_path, on_path, monkeypatch, sensitive, param):
    on_path(_fake(tmp_path, "miniprot"))
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age, "GOI_MINIPROT_TIMEOUT", 1)
    monkeypatch.setattr(age, "FLANKING_MINIPROT_TIMEOUT", 1)
    kind = "goi" if sensitive else "flanking"
    with pytest.raises(su.StepFailed, match=rf"miniprot \({kind} model: a {len(QUERY)}-aa protein against "
                                            rf"{len(DNA):,} bp of chr1\) did not finish within 1 s \(--{param}\)"):
        _quick(age._run_miniprot, QUERY, DNA, "chr1", sensitive=sensitive)
    # the copy rounds of a GOI locus pass it on (no round ends on time)
    with pytest.raises(su.StepFailed):
        _quick(age.miniprot_gene_models, QUERY, DNA, "chr1", sensitive=sensitive, max_models=5)
    assert not any(k.startswith("miniprot_timeout") or k.endswith("by_timeout") for k in age.LIMIT_COUNTS)


def test_b_the_guard_is_not_caught_by_the_searchs_gene_model_fallback(tmp_path, on_path, monkeypatch):
    """annotate_exons_from_hit_list catches Exception around the gene model; the guard passes."""
    on_path(_fake(tmp_path, "miniprot"))
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age, "GOI_MINIPROT_TIMEOUT", 1)
    hits = [{"qstart": 1, "qend": 60, "start": 0, "end": 180, "strand": "+", "pident": 90.0,
             "evalue": 1e-20, "bits": 100.0, "chrom": "chr1", "query": "GOI_q", "seq": QUERY[:60]}]
    with pytest.raises(su.StepFailed):
        _quick(age.annotate_exons_from_hit_list, hits, QUERY, DNA, "chr1", sensitive=True)


# ------------------------------------------------------------------ d: ANNOTATE_GOI's home model
def test_d_the_home_model_run_beyond_its_guard_stops_annotate_goi(tmp_path):
    fake = _fake(tmp_path, "miniprot")
    (tmp_path / "q.faa").write_text(f">q1\n{QUERY}\n")
    (tmp_path / "g.fna").write_text(f">chr1\n{DNA}\n>chr2\n{DNA}\n")
    with open(tmp_path / "hits.m8", "w") as fh:          # hits on two chromosomes
        for chrom in ("chr1", "chr2"):
            fh.write(f"q1\t{chrom}\t90.0\t60\t0\t0\t1\t60\t1\t180\t1e-20\t100.0\n")
    env = dict(os.environ, PATH=fake + os.pathsep + os.environ["PATH"])
    out = {k: tmp_path / f"out.{k}" for k in ("faa", "bed", "json")}
    t0 = time.monotonic()
    r = subprocess.run([sys.executable, str(BIN / "annotate_goi_exons.py"), "--query", str(tmp_path / "q.faa"),
                        "--genome", str(tmp_path / "g.fna"), "--gff", "NO_GFF",
                        "--blast_hits", str(tmp_path / "hits.m8"), "--query_id", "q1",
                        "--miniprot_timeout", "1", "--output_exons", str(out["faa"]),
                        "--output_bed", str(out["bed"]), "--output_info", str(out["json"])],
                       capture_output=True, text=True, env=env)
    assert time.monotonic() - t0 < 30                   # not one guard per chromosome in turn
    assert r.returncode != 0
    assert "did not finish within 1 s (--goi_miniprot_timeout)" in r.stderr
    assert not any(p.exists() for p in out.values())   # no home model from the next chromosome


# ------------------------------------------------------------------ e: BORROW_ANNOTATIONS
def test_e_a_transfer_beyond_its_guard_stops_the_borrowing(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "miniprot"))
    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age, "FLANKING_MINIPROT_TIMEOUT", 1)
    with pytest.raises(su.StepFailed, match=r"--flanking_miniprot_timeout"):
        _quick(ba.borrowed_models, [("g1", "n1", QUERY)], DNA, "chr1", 0, 2000)


def test_e_the_script_takes_its_guard_from_the_pipeline():
    nf = (ROOT / "modules" / "borrow_annotations.nf").read_text()
    assert "--miniprot_timeout ${params.flanking_miniprot_timeout}" in nf
    src = (BIN / "borrow_annotations.py").read_text()
    assert "ag.FLANKING_MINIPROT_TIMEOUT = args.miniprot_timeout" in src


# ------------------------------------------------------------------ f, g: the rescue passes
def _scores(path):
    fields = ["region_rank", "region_name", "species", "chrom", "start", "end", "strand", "score",
              "region_class", "goi_missing", "confidence", "high_flanking_count"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerow({"region_rank": 1, "region_name": "r1", "species": "X", "chrom": "chr1",
                    "start": 100, "end": 2000, "strand": "+", "score": 0.9,
                    "region_class": "goi_missing_but_strong_synteny", "goi_missing": "True",
                    "confidence": "MEDIUM", "high_flanking_count": 6})


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_f_the_strong_synteny_rescue_ends_the_run_with_exit_3(tmp_path):
    fake = _fake(tmp_path, "miniprot")
    _scores(tmp_path / "X.scores.tsv")
    (tmp_path / "t.fna").write_text(f">chr1\n{DNA}\n")
    (tmp_path / "q.faa").write_text(f">GOI_q\n{QUERY}\n")
    env = dict(os.environ, PATH=fake + os.pathsep + os.environ["PATH"])
    r = subprocess.run([sys.executable, str(BIN / "rescue_strong_synteny.py"), "--scores",
                        str(tmp_path / "X.scores.tsv"), "--target_genome", str(tmp_path / "t.fna"),
                        "--query", str(tmp_path / "q.faa"), "--genome_name", "X",
                        "--output", str(tmp_path / "out.gff"), "--miniprot_timeout", "1"],
                       capture_output=True, text=True, env=env, timeout=120)
    assert r.returncode == su.STEP_FAILED_EXIT == 3, r.stderr
    assert ("miniprot on X block 1 (chr1:1-" in r.stderr
            and "did not finish within 1 s (--strong_synteny_rescue_miniprot_timeout)" in r.stderr)
    assert not (tmp_path / "out.gff").exists()       # no GFF that reads like "no model"


def test_g_the_hull_rescue_stops_on_its_guard(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "miniprot"))
    rows = ["##gff-version 3"] + [
        f"chr1\tx\tmRNA\t{1_000_000 + i * 200_000}\t{1_010_000 + i * 200_000}\t.\t+\t.\t"
        f"ID=f{i};SynVoyRole=flanking;Confidence=HIGH" for i in range(4)]
    (tmp_path / "region.gff").write_text("\n".join(rows) + "\n")
    (tmp_path / "q.faa").write_text(f">GOI_q\n{QUERY}\n")
    (tmp_path / "eff.json").write_text(json.dumps({"settings": {}}))
    monkeypatch.setattr(hull, "_extract_window", lambda *a: DNA)
    monkeypatch.setattr(hull, "_genome_length", lambda path: 10_000_000)
    monkeypatch.setattr(sys, "argv", [
        "rescue_goi_hull.py", "--region_gff", str(tmp_path / "region.gff"),
        "--target_genome", str(tmp_path / "t.fna"), "--query", str(tmp_path / "q.faa"),
        "--genome_name", "t", "--output", str(tmp_path / "out.gff"),
        "--cluster_max_gap", "3000000", "--window_pad", "10000", "--max_window", "20000000",
        "--effective_params", str(tmp_path / "eff.json"), "--miniprot_timeout", "1"])
    with pytest.raises(su.StepFailed, match=r"miniprot on t hull chr1:\d+-\d+ \([\d,]+ bp\) did not "
                                            r"finish within 1 s \(--goi_miniprot_timeout\)"):
        _quick(hull.main)
    with pytest.raises(SystemExit) as exc:
        su.exit_step_failed(hull.main)
    assert exc.value.code == 3
    assert not (tmp_path / "out.gff").exists()


def test_f_g_a_rescue_without_miniprot_stops_the_run_too(tmp_path, monkeypatch):
    """A missing miniprot is an installation fault, not a window: exit 3 like a guard (until the
    review of 2026-10-30 it raised FileNotFoundError, exit 1, which the modules ignore)."""
    monkeypatch.setenv("PATH", str(tmp_path))              # no miniprot anywhere
    (tmp_path / "w.fa").write_text(f">w\n{DNA}\n")
    (tmp_path / "q.faa").write_text(f">q\n{QUERY}\n")
    with pytest.raises(su.StepFailed, match=r"\[rescue\] miniprot is not installed here"):
        rss._run_miniprot_relaxed(str(tmp_path / "w.fa"), str(tmp_path / "q.faa"), 0.05, 5)


def test_f_g_the_rescue_modules_end_the_run_on_exit_3_and_on_spent_kills():
    """Exit 3 (a guard) ends the run; a kill for memory or time (137/140/143) is retried twice and
    then ends the run too (until 2026-10-28 a third kill dropped the genome's rescue in silence);
    other failures are still ignored. Checked on Nextflow 26.04.4 (preview/errstrat)."""
    for name in ("rescue_strong_synteny.nf", "rescue_goi_hull.nf"):
        nf = (ROOT / "modules" / name).read_text()
        assert ("errorStrategy { task.exitStatus == 3 ? 'terminate' : (task.exitStatus in [137, 140, 143] "
                "? (task.attempt <= 2 ? 'retry' : 'terminate') : 'ignore') }") in nf
    assert "--miniprot_timeout ${params.goi_miniprot_timeout}" in (ROOT / "modules" / "rescue_goi_hull.nf").read_text()
    for script in ("rescue_strong_synteny.py", "rescue_goi_hull.py"):
        assert "    exit_step_failed(main)" in (BIN / script).read_text()


def test_exit_step_failed_passes_a_return_value_and_maps_the_guard_to_3():
    with pytest.raises(SystemExit) as exc:
        su.exit_step_failed(lambda: 0)
    assert exc.value.code == 0

    def boom():
        raise su.StepFailed("x did not finish")
    with pytest.raises(SystemExit) as exc:
        su.exit_step_failed(boom)
    assert exc.value.code == su.STEP_FAILED_EXIT


# ------------------------------------------------------------------ h, i, j: the gene prediction
def _augustus_config(tmp_path, monkeypatch, model="fly"):
    cfg = tmp_path / "augcfg"
    (cfg / "species" / model).mkdir(parents=True)
    (cfg / "species" / model / f"{model}_parameters.cfg").write_text("# fake\n")
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(cfg))


def test_h_augustus_beyond_its_guard_stops_the_prediction(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "augustus"))
    _augustus_config(tmp_path, monkeypatch)
    monkeypatch.setattr(gp, "GENE_PREDICTION_TIMEOUT", 1)
    (tmp_path / "r.fna").write_text(f">r1\n{DNA}\n")
    with pytest.raises(su.StepFailed, match=r"Augustus \(fly\) on .*did not finish within 1 s "
                                            r"\(--gene_prediction_timeout\)"):
        _quick(gp.run_augustus, str(tmp_path / "r.fna"), str(tmp_path / "o.faa"),
               str(tmp_path / "o.gff"), species="fly")


def test_h_a_failed_augustus_run_does_not_fall_back_to_prodigal(tmp_path, on_path, monkeypatch):
    folder = _fake(tmp_path, "augustus", body="echo broken >&2; exit 1")
    _fake(tmp_path, "prodigal", body="touch prodigal_ran; exit 0")
    on_path(folder)
    _augustus_config(tmp_path, monkeypatch)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "r.fna").write_text(f">r1\n{DNA}\n")
    with pytest.raises(su.StepFailed, match=r"Augustus \(fly\) failed on .*\(exit 1: broken\); the "
                                            r"prediction does not fall back to Prodigal"):
        gp.predict_genes(str(tmp_path / "r.fna"), str(tmp_path / "o.faa"), str(tmp_path / "o.gff"),
                         predictor="augustus", augustus_species="fly")
    with pytest.raises(su.StepFailed, match=r"Augustus \(fly\) failed .*exit 1: broken"):
        gp.predict_orfs(str(tmp_path / "r.fna"), str(tmp_path), predictor="augustus", augustus_species="fly")
    assert not (tmp_path / "prodigal_ran").exists()


def test_h_a_missing_augustus_stops_the_prediction_and_says_so(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))              # no augustus anywhere
    monkeypatch.setattr(gp, "resolve_augustus_species", lambda s: s)
    (tmp_path / "r.fna").write_text(f">r1\n{DNA}\n")
    with pytest.raises(su.StepFailed, match="Augustus is not installed here"):
        gp.predict_genes(str(tmp_path / "r.fna"), str(tmp_path / "o.faa"), str(tmp_path / "o.gff"),
                         predictor="augustus", augustus_species="fly")


def test_h_the_discovery_refuses_to_start_without_augustus(tmp_path, monkeypatch):
    """The structural discovery predicts with Augustus only; without it the search failed only
    after every wave (review of 2026-10-29), now it stops at its start."""
    import argparse
    import iterative_search_runner as runner
    monkeypatch.setenv("PATH", str(tmp_path))
    on = argparse.Namespace(enable_plm_search=True, enable_structural_search=False,
                            enable_structural_discovery=True, augustus_species="fly")
    assert "Augustus is not installed here" in runner.augustus_steps_model_error(on)
    off = argparse.Namespace(**{**vars(on), "enable_structural_discovery": False})
    assert runner.augustus_steps_model_error(off) is None
    src = (BIN / "iterative_search_runner.py").read_text()
    assert "_gp.GENE_PREDICTION_TIMEOUT = args.gene_prediction_timeout" in src
    assert "--gene_prediction_timeout ${params.gene_prediction_timeout}" in (
        ROOT / "modules" / "iterative_search.nf").read_text()


def test_i_prodigal_beyond_its_guard_stops_the_prediction(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "prodigal"))
    monkeypatch.setattr(gp, "GENE_PREDICTION_TIMEOUT", 1)
    (tmp_path / "r.fna").write_text(f">r1\n{DNA}\n")
    with pytest.raises(su.StepFailed, match=r"Prodigal on .*\(--gene_prediction_timeout\)"):
        _quick(gp.run_prodigal, str(tmp_path / "r.fna"), str(tmp_path / "o.faa"), str(tmp_path / "o.gff"))


def test_h_i_both_prediction_scripts_take_the_guard_from_the_pipeline():
    for script, module in (("prodigal_on_regions.py", "prepare_home.nf"),
                           ("extract_flanking_genes.py", "extract_flanking.nf")):
        src = (BIN / script).read_text()
        assert "gene_predictor.GENE_PREDICTION_TIMEOUT = args.gene_prediction_timeout" in src
        nf = (ROOT / "modules" / module).read_text()
        assert "--gene_prediction_timeout ${params.gene_prediction_timeout}" in nf
    assert (ROOT / "modules" / "prepare_home.nf").read_text().count("--gene_prediction_timeout") == 2


def test_j_an_augustus_whose_model_list_fails_stops_the_run(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "augustus", body="echo usage only >&2; exit 0"))
    _augustus_config(tmp_path, monkeypatch)
    with pytest.raises(su.StepFailed, match="listed no model"):
        gp.installed_augustus_models()
    src = (BIN / "gene_predictor.py").read_text()
    assert "list_augustus_species" not in src                  # the unused 10-s list is gone


@pytest.mark.parametrize("predictor,stops", [("auto", True), ("prodigal", False)])
def test_j_a_failed_model_list_stops_the_resolver_only_when_augustus_runs(tmp_path, predictor, stops):
    """--gene_predictor prodigal starts no Augustus, so a broken one does not stop that run."""
    fake = _fake(tmp_path, "augustus", body="echo usage only >&2; exit 0")
    defaults = tmp_path / "defaults.json"
    defaults.write_text(json.dumps({"goi_margin_genes": 2}))
    env = dict(os.environ, PATH=fake + os.pathsep + os.environ["PATH"])
    proc = subprocess.run([sys.executable, str(BIN / "resolve_effective_params.py"), "--defaults", str(defaults),
                           "--output", str(tmp_path / "ep.json"), "--modules", "auto", "--home_gff_present",
                           "false", "--home_species", "Homo sapiens", "--gene_predictor", predictor],
                          capture_output=True, text=True, env=env)
    assert "listed no model" in proc.stderr
    assert (proc.returncode != 0) == stops, proc.stderr
    if not stops:      # the record names the failed list, not a missing Augustus
        line = json.loads((tmp_path / "ep.json").read_text())["derivation"]["augustus_species"]["line"]
        assert "the installed Augustus gave no model list" in line and "not installed here" not in line
        assert "runs no Augustus" not in proc.stderr


def test_j_a_model_that_is_not_installed_is_named_also_when_the_list_fails(tmp_path, on_path, monkeypatch):
    """Item 16's message (the hint, the model, the file Augustus lacks) also when the model list
    fails, e.g. $AUGUSTUS_CONFIG_PATH names a folder without species (review of 2026-10-30: the
    list's StepFailed replaced the message)."""
    on_path(_fake(tmp_path, "augustus", body="echo 'config path broken' >&2; exit 1"))
    _augustus_config(tmp_path, monkeypatch)
    with pytest.raises(gp.AugustusModelNotInstalled) as err:
        gp.resolve_augustus_species("silkworm")
    text = str(err.value)
    assert "augustus_species 'silkworm' maps to the Augustus model 'bombyx_mori'" in text
    assert "The installed models could not be listed: `augustus --species=help` listed no model" in text


def test_j_the_model_list_has_no_time_limit(tmp_path, on_path, monkeypatch):
    on_path(_fake(tmp_path, "augustus", body="sleep 2; echo 'fly | Drosophila melanogaster' >&2"))
    _augustus_config(tmp_path, monkeypatch)
    assert gp.installed_augustus_models() == {"fly": "Drosophila melanogaster"}


# ------------------------------------------------------------------ k: the miniprot probes
def test_k_the_miniprot_probes_have_no_time_limit():
    import inspect
    for fn in (age._check_miniprot, age._miniprot_genetic_code_args):
        assert "timeout" not in inspect.getsource(fn).split('"""')[-1]


def test_k_the_search_refuses_to_start_without_miniprot():
    src = (BIN / "iterative_search_runner.py").read_text()
    main = src[src.index("    if args.flanking_miniprot_timeout < 1:"):]
    block = main[:main.index("sys.exit(1)", main.index("if not MINIPROT_AVAILABLE:")) + 12]
    assert "if not MINIPROT_AVAILABLE:" in block and "miniprot not found" in block


def test_every_guard_is_a_pipeline_parameter_validated_at_launch():
    cfg = (ROOT / "nextflow.config").read_text()
    for name in ("sw_timeout_seconds", "goi_miniprot_timeout", "flanking_miniprot_timeout",
                 "strong_synteny_rescue_miniprot_timeout", "gene_prediction_timeout"):
        assert f"    {name} = 3600" in cfg or f"    {name} = 3600 " in cfg or f"{name} = 3600" in cfg, name
    main_nf = (ROOT / "main.nf").read_text()
    assert ("['sw_timeout_seconds', 'goi_miniprot_timeout', 'flanking_miniprot_timeout',\n"
            "     'strong_synteny_rescue_miniprot_timeout', 'gene_prediction_timeout'].each") in main_nf
    # whole numbers only: the scripts read the guards as int (3600.0 stopped a rescue with exit 2,
    # which its module ignores)
    assert "if (!(params[guardName].toString() ==~ /[1-9][0-9]*/)) {" in main_nf
