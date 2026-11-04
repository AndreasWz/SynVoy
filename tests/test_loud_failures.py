"""Fix-phase d5 (2 Nov; items 21 and 35): a result that a task fails to produce is never silent.

- bin/run_step.py records a failed paralog check, ownership check, rescue or (with --require_paralog_panel false) panel in
  the step's own output, where errorStrategy 'ignore' left nothing; a kill or a guard (exit 3) keeps its code.
- generate_report.py names every recorded failure (warnings, failed_steps) and marks the calls a failed ownership check
  would have checked (ownership_not_checked).
- the search: a window whose tblastn fails, or whose GOI search fails otherwise, stops the genome's search (it kept the
  window without those hits, or without any), and a window's temporary files go also when it fails.
- the genomes are indexed once (bin/index_genomes.py) before the rescues read them side by side.
- the strong-synteny rescue reads a plain-gzip target (item 35): samtools faidx refuses it, and every block was skipped.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import sequence_utils as su  # noqa: E402


def _run_step(tmp_path, code_or_cmd, *extra):
    out = tmp_path / "x.tsv"
    out.write_text("stale\n")
    faa = tmp_path / "x.faa"
    faa.write_text(">a\nMK\n")
    cmd = code_or_cmd if isinstance(code_or_cmd, list) else \
        ["bash", "-c", f"echo 'line one' >&2; echo 'the last line' >&2; exit {code_or_cmd}"]
    p = subprocess.run([sys.executable, os.path.join(BIN, "run_step.py"), "--step", "ASSIGN_LOCUS_OWNERSHIP",
                        "--locus", "locus_1", "--genome", "Gx.fna", "--record", str(out), "--empty", str(faa),
                        *extra, "--", *cmd], capture_output=True, text=True)
    return p, out.read_text(), faa.read_text()


# --------------------------------------------------------------------------- run_step.py

def test_a_success_passes_through(tmp_path):
    p, out, faa = _run_step(tmp_path, ["bash", "-c", "echo hi; echo result > " + str(tmp_path / "x.tsv")])
    assert p.returncode == 0 and p.stdout == "hi\n" and out == "result\n" and faa == ">a\nMK\n"


def test_a_failure_is_recorded_in_the_steps_output(tmp_path):
    p, out, faa = _run_step(tmp_path, 1)
    assert p.returncode == 0
    assert out == ("# SYNVOY_STEP_FAILED\tstep=ASSIGN_LOCUS_OWNERSHIP\tlocus=locus_1\tgenome=Gx.fna\texit=1\t"
                   "reason=the last line\n")
    assert faa == "" and "the last line" in p.stderr and "recorded for the report" in p.stderr


def test_a_header_goes_before_the_line(tmp_path):
    _p, out, _f = _run_step(tmp_path, 2, "--header", "##gff-version 3")
    assert out.startswith("##gff-version 3\n# SYNVOY_STEP_FAILED\t")


@pytest.mark.parametrize("code", [3, 137, 140, 143])
def test_a_guard_or_a_kill_keeps_its_code(tmp_path, code):
    p, out, faa = _run_step(tmp_path, code)
    assert p.returncode == code and out == "stale\n" and faa == ">a\nMK\n"


def test_a_signal_reads_as_the_shells_code(tmp_path):
    p, _out, _f = _run_step(tmp_path, ["bash", "-c", "kill -9 $$"])
    assert p.returncode == 137


# --------------------------------------------------------------------------- the report

def _gff_row(start, end, mid, conf="HIGH"):
    attrs = (f"ID={mid};Name=GOI_G;SynVoy_Parent=GOI_G;SynVoyRole=goi;Confidence={conf};GOIClass=confident_goi;"
             f"EvidenceType=exon_annotation;Identity=80.0;QueryCoverage=0.950;ModelStatus=complete;"
             f"Basis=significant_specific")
    return "\t".join(["c1", "exon_annotation", "mRNA", str(start), str(end), "80", "+", ".", attrs]) + "\n"


MARK = "# SYNVOY_STEP_FAILED\tstep={step}\tlocus={locus}\tgenome={genome}\texit=1\treason=boom\n"


def _report(tmp_path, ownership_failed=(), panel_failed=False, rescue_failed=False):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "locus_ownership").mkdir()
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n" + _gff_row(1000, 2000, "GOI_G|m1"))
    (stage / "regions" / "locus_1__Gy.fna.gff").write_text("##gff-version 3\n" + _gff_row(5000, 6000, "GOI_G|m2"))
    for locus, genome in ownership_failed:
        (stage / "locus_ownership" / f"{locus}__{genome}.locus_ownership.tsv").write_text(
            MARK.format(step="ASSIGN_LOCUS_OWNERSHIP", locus=locus, genome=genome))
    if panel_failed:
        (stage / "locus_ownership" / "home_paralog_panel.faa.meta.tsv").write_text(
            MARK.format(step="BUILD_HOME_PARALOG_PANEL", locus="", genome=""))
    if rescue_failed:
        (stage / "regions" / "locus_1__Gy.hull_rescue.gff").write_text(
            "##gff-version 3\n" + MARK.format(step="RESCUE_GOI_HULL", locus="locus_1", genome="Gy"))
    out = tmp_path / "report.json"
    proc = subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                           "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(out.read_text())


def test_a_run_without_failures_has_no_new_keys(tmp_path):
    rep = _report(tmp_path)
    assert "failed_steps" not in rep and "ownership_not_checked_calls" not in rep["summary"]
    assert all("ownership_not_checked" not in r for r in rep["goi_dedup"]["records"])
    assert not any("failed_steps" in w for w in rep["warnings"])


def test_a_failed_ownership_check_marks_its_calls_only(tmp_path):
    rep = _report(tmp_path, ownership_failed=[("locus_1", "Gx.fna")])
    recs = {r["genome"]: r for r in rep["goi_dedup"]["records"]}
    assert recs["Gx"].get("ownership_not_checked") is True and "ownership_not_checked" not in recs["Gy"]
    assert rep["summary"]["ownership_not_checked_calls"] == 1
    assert rep["failed_steps"] == [{"step": "ASSIGN_LOCUS_OWNERSHIP", "locus": "locus_1", "genome": "Gx",
                                    "exit": "1", "reason": "boom"}]
    assert any("ASSIGN_LOCUS_OWNERSHIP failed for Gx, locus_1 (exit 1: boom)" in w for w in rep["warnings"])


def test_another_locus_does_not_mark_the_call(tmp_path):
    rep = _report(tmp_path, ownership_failed=[("locus_2", "Gx.fna")])
    assert all("ownership_not_checked" not in r for r in rep["goi_dedup"]["records"])
    assert rep["summary"]["ownership_not_checked_calls"] == 0


def test_a_failed_panel_marks_every_call(tmp_path):
    rep = _report(tmp_path, panel_failed=True)
    assert all(r.get("ownership_not_checked") for r in rep["goi_dedup"]["records"])
    assert any("BUILD_HOME_PARALOG_PANEL failed" in w for w in rep["warnings"])


def test_a_failed_rescue_is_named(tmp_path):
    rep = _report(tmp_path, rescue_failed=True)
    assert rep["failed_steps"][0]["step"] == "RESCUE_GOI_HULL" and rep["failed_steps"][0]["genome"] == "Gy"
    assert any("RESCUE_GOI_HULL failed for Gy, locus_1" in w and "hull rescue would model" in w
               for w in rep["warnings"])
    assert all("ownership_not_checked" not in r for r in rep["goi_dedup"]["records"])


# --------------------------------------------------------------------------- the modules

@pytest.mark.parametrize("module,step", [("assign_locus_ownership.nf", "ASSIGN_LOCUS_OWNERSHIP"),
                                         ("assign_locus_ownership.nf", "BUILD_HOME_PARALOG_PANEL"),
                                         ("reciprocal_best_paralog.nf", "RECIPROCAL_BEST_PARALOG"),
                                         ("rescue_goi_hull.nf", "RESCUE_GOI_HULL"),
                                         ("rescue_strong_synteny.nf", "RESCUE_STRONG_SYNTENY")])
def test_the_steps_record_and_never_ignore(module, step):
    text = open(os.path.join(ROOT, "modules", module)).read()
    assert f"run_step.py --step {step}" in text
    assert "errorStrategy 'ignore'" not in text and ": 'ignore'" not in text and "? 'terminate' : 'ignore'" not in text


def test_no_module_ignores_a_failure_any_more():
    for name in sorted(os.listdir(os.path.join(ROOT, "modules"))):
        text = open(os.path.join(ROOT, "modules", name)).read()
        code = "\n".join(l for l in text.splitlines() if not l.strip().startswith("//"))
        assert "'ignore'" not in code, name


@pytest.mark.parametrize("module,fallbacks", [
    ("locate_gene.nf", ["MMSeqs search failed", "hit profiling skipped", "preset advisor skipped"]),
    ("home_context.nf", ["|| echo"]),
    ("generate_report.nf", ["|| true", "2>/dev/null\n"]),
    ("plot_synteny.nf", ["|| true", "|| cp"])])
def test_no_silent_fallback_in(module, fallbacks):
    text = open(os.path.join(ROOT, "modules", module)).read()
    code = "\n".join(l for l in text.splitlines() if not l.strip().startswith(("//", "#")))
    for fallback in fallbacks:
        assert fallback not in code, (module, fallback)


# --------------------------------------------------------------------------- the search

def _m8(path, rows):
    with open(path, "w") as fh:
        for q, aln, e in rows:
            fh.write(f"{q}\tregion_seq\t60.0\t{aln}\t0\t0\t1\t{aln}\t101\t{100 + 3 * aln}\t{e}\t60.0\n")


def _args(tmp_path, **over):
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
             aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
             aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6, min_identity=10.0,
             aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=False, deterministic_goi_search=True,
             sw_min_score=20, sw_min_identity=15.0, sw_method="parasail",
             sw_max_alignments_per_piece=1, sw_timeout_seconds=3600, genome_ctx={})
    a.update(over)
    return types.SimpleNamespace(**a)


def _window(tmp_path, monkeypatch, fake_run, fake_mmseqs=None):
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")

    def mmseqs(**k):
        _m8(k["aug_hits_file"], [("GOI_q", 40, "1e-20")])
        return True, ""
    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", fake_mmseqs or mmseqs)
    return isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": "M" + "A" * 119}],
                                    "g1_b3", _args(tmp_path), "t1", 1, region_len=2_000)


def _leftovers(tmp_path):
    return sorted(p for p in os.listdir(tmp_path) if p.startswith("tmp_"))


def test_a_failed_tblastn_stops_the_window(tmp_path, monkeypatch):
    def fake_run(cmd, *a, **k):
        if cmd[0] == "makeblastdb":
            open(f"{cmd[cmd.index('-out') + 1]}.nhr", "w").close()
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        raise subprocess.CalledProcessError(2, cmd, stderr="BLAST Database error")
    monkeypatch.setattr(isr, "run_tblastn_tabular_and_json",
                        lambda cmd, out: (_ for _ in ()).throw(subprocess.CalledProcessError(2, cmd, stderr="db error")))
    with pytest.raises(su.StepFailed) as stop:
        _window(tmp_path, monkeypatch, fake_run)
    assert "tblastn of the GOI queries against a window failed (exit 2)" in str(stop.value)
    assert "db error" in str(stop.value)
    assert _leftovers(tmp_path) == []


def test_a_missing_tblastn_stops_the_window(tmp_path, monkeypatch):
    def fake_run(cmd, *a, **k):
        raise FileNotFoundError(cmd[0])
    with pytest.raises(su.StepFailed) as stop:
        _window(tmp_path, monkeypatch, fake_run)
    assert "Install BLAST+" in str(stop.value) and _leftovers(tmp_path) == []


def test_any_other_failure_of_the_window_stops_it(tmp_path, monkeypatch):
    def fake_run(cmd, *a, **k):
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    def broken_mmseqs(**k):
        open(k["aug_hits_file"], "w").write("not\ta\ttable\n")
        os.makedirs(os.path.join(k["aug_tmp_dir"], "sub"), exist_ok=True)
        return True, ""
    monkeypatch.setattr(isr, "parse_hits", lambda *a, **k: (_ for _ in ()).throw(ValueError("bad row")))
    with pytest.raises(su.StepFailed) as stop:
        _window(tmp_path, monkeypatch, fake_run, broken_mmseqs)
    assert "the GOI search of a window failed (ValueError: bad row" in str(stop.value)
    assert _leftovers(tmp_path) == []


def test_step_failed_reaches_the_genome_not_the_window_catch():
    """A StepFailed is a BaseException: the block's own catch (now a StepFailed too) does not hold it."""
    assert not issubclass(su.StepFailed, Exception)


# --------------------------------------------------------------------------- the genomes

def test_index_genomes_indexes_once_and_names_a_plain_gzip(tmp_path):
    if not shutil.which("samtools"):
        pytest.skip("samtools not installed")
    d = tmp_path / "g"
    d.mkdir()
    (d / "a.fna").write_text(">c1\nACGTACGT\n")
    with gzip.open(d / "b.fna.gz", "wt") as fh:
        fh.write(">c1\nACGTACGT\n")
    (d / "notes.txt").write_text("x")
    p = subprocess.run([sys.executable, os.path.join(BIN, "index_genomes.py"), "--dir", str(d)],
                       capture_output=True, text=True)
    assert p.returncode == 0
    assert (d / "a.fna.fai").exists() and not (d / "b.fna.gz.fai").exists() and not (d / "notes.txt.fai").exists()
    assert "b.fna.gz: no index (a plain-gzip file)" in p.stderr


def test_the_strong_synteny_rescue_reads_a_plain_gzip_target(tmp_path):
    import rescue_strong_synteny as rss
    g = tmp_path / "t.fna.gz"
    with gzip.open(g, "wt") as fh:
        fh.write(">c1\n" + "A" * 10 + "CCCGGG" + "T" * 10 + "\n")
    assert rss._read_fasta_subseq("samtools", str(g), "c1", 11, 16) == "CCCGGG"
    assert rss._read_fasta_subseq("no-such-samtools", str(g), "c1", 11, 16) == "CCCGGG"
    assert su.genome_window(str(g), "c9", 1, 5) == ""


def test_a_crash_by_another_signal_is_recorded(tmp_path):
    """Review of the session: a segfault (139) passed through and ended the run; only 3 and the kill codes do."""
    p, out, faa = _run_step(tmp_path, ["bash", "-c", "kill -SEGV $$"])
    assert p.returncode == 0 and "exit=139" in out and faa == ""


def test_stderr_that_is_not_utf8_is_recorded(tmp_path):
    p, out, _f = _run_step(tmp_path, [sys.executable, "-c", "import sys; sys.stderr.buffer.write(b'\\xff\\xfe bad\\n'); sys.exit(1)"])
    assert p.returncode == 0 and out.startswith("# SYNVOY_STEP_FAILED") and "bad" in out


def test_a_failed_hull_rescue_ownership_marks_only_the_hull_records(tmp_path):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "locus_ownership").mkdir()
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n" + _gff_row(1000, 2000, "GOI_G|m1"))
    hull = _gff_row(5000, 6000, "GOI_G|h1").replace("GOIClass=confident_goi", "GOIClass=synteny_hull_rescue")
    (stage / "regions" / "locus_1__Gx.hull_rescue.gff").write_text("##gff-version 3\n" + hull)
    (stage / "locus_ownership" / "locus_1__Gx.hull_rescue.locus_ownership.tsv").write_text(
        MARK.format(step="ASSIGN_LOCUS_OWNERSHIP", locus="locus_1", genome="Gx.hull_rescue"))
    out = tmp_path / "report.json"
    subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                    "--output", str(out), "--allow-empty"], check=True, capture_output=True)
    rep = json.loads(out.read_text())
    marked = {r["start"]: bool(r.get("ownership_not_checked")) for r in rep["goi_dedup"]["records"]}
    assert marked == {1000: False, 5000: True}
    assert rep["failed_steps"][0]["only"] == "hull_rescue"


# --------------------------------------------------------------------------- --preset_override

def test_a_mistyped_preset_override_stops_naming_the_presets(tmp_path):
    """Coordinator, 2 Nov: keep the stop (before, LOCATE_GENE ignored a mistyped override without a word), with the
    valid names in the message; main.nf stops at its start with the same names."""
    import auto_select_preset as asp
    (tmp_path / "hp.json").write_text('{"n_hits": 1}')
    p = subprocess.run([sys.executable, os.path.join(BIN, "auto_select_preset.py"), "--hit_profile", str(tmp_path / "hp.json"),
                        "--output", str(tmp_path / "a.json"), "--preset_override", "preset_singlecopy"],
                       capture_output=True, text=True)
    assert p.returncode == 2 and "'preset_singlecopy' is not a preset" in p.stderr
    assert all(n in p.stderr for n in asp.VALID_PRESETS)
    main = open(os.path.join(ROOT, "main.nf")).read()
    names = main[main.index("def presetNames() {"):].split("[", 1)[1].split("]", 1)[0]
    assert sorted(x.strip().strip("'") for x in names.split(",")) == sorted(asp.VALID_PRESETS)
    assert "checkPresetOverride(validationErrors)" in main
