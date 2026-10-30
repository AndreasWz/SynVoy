"""Fix-phase item 18, Smith-Waterman: the wall clock decides whether a result exists, never what
it contains.

Until 2026-10-28 the runner gave each window's Smith-Waterman child 300 s (sw_timeout_seconds);
on a timeout the window kept its MMseqs2 and tblastn hits only, so a slower or busier machine
got other models, bits and labels from the same inputs (fly case without its home GFF: one
window timed out in one of two back-to-back runs). Now:

- the child scores every piece without a traceback first and runs the traceback only where a
  hit the runner keeps can come out (min_reportable_score), so the kept hits are the same and
  the time is several times shorter;
- the timeout is a guard: a search that does not finish within it, a child that fails or one
  that writes no hits file stops the genome's search (StepFailed), and the wave ends the task;
- the cells of each window are logged, and the backend that ran is recorded in the
  target-policy tables and the report's run_settings.
"""
from __future__ import annotations

import os
import random
import re
import subprocess
import sys
import threading
import time
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import smith_waterman_search as sws  # noqa: E402
import sequence_utils as su  # noqa: E402

CODON = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
         'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
         'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'}
AA = "ACDEFGHIKLMNPQRSTVWY"


def _kept(rows):
    """The m8 lines the runner keeps from the child's file, in file order: E <= 1 as it reads
    the printed value (parse_hits: float(column 11) <= 1.0)."""
    return [r for r in rows if float(r.split("\t")[10]) <= 1.0]


# ------------------------------------------------------------------ the score threshold
@pytest.mark.parametrize("qlen,space", [(30, 1_000), (70, 60_000), (256, 2_000_000),
                                        (2048, 236_000_000), (5, 10)])
def test_min_reportable_score_is_one_below_the_first_score_with_e_at_most_1(qlen, space):
    s = sws.min_reportable_score(qlen, space, 11, 1)

    def e(score):
        return sws.karlin_altschul(score, qlen, space, 11, 1)[1]
    if s >= 0:
        assert e(s) > 1.0 >= e(s + 1)
    else:
        assert e(0) <= 1.0                      # a tiny search space: every score is kept
    # no score below the threshold can be printed as E <= 1 ('%.2e', as the child writes it)
    for score in range(max(0, s - 40), s):
        assert float(f"{e(score):.2e}") > 1.0


def test_a_hit_printed_as_1_00e_plus_00_is_above_the_threshold():
    """The runner reads the printed E: a score whose E is 1.004 prints as 1.00e+00 and is kept,
    so the threshold sits one raw score below the first score with E <= 1."""
    found = None
    for space in range(10_000, 400_000, 7):
        s = sws.min_reportable_score(100, space, 11, 1)
        e = sws.karlin_altschul(s, 100, space, 11, 1)[1]
        if 1.0 < e < 1.005:
            found = (space, s, e)
            break
    assert found, "no search space gives E just above 1 at the threshold"
    _space, s, e = found
    assert f"{e:.2e}" == "1.00e+00"           # kept by the runner, and at the threshold itself


def test_six_frame_length_is_the_childs_search_space():
    rng = random.Random(1)
    for n in (0, 1, 2, 3, 4, 5, 100, 1_001, 30_002):
        dna = ''.join(rng.choice("ACGTN") for _ in range(n))
        frames = sws.translate_in_six_frames(dna)
        assert sws.six_frame_length(n) == sum(len(p) for _f, p, _o in frames)
    assert sws.smith_waterman_cells(3_000, [100, 12_000, 50]) == 150 * sws.six_frame_length(3_000)


# ------------------------------------------------------------------ the two parasail kernels
parasail = pytest.importorskip("parasail")


def test_score_only_and_traceback_kernels_agree_on_score_and_end():
    """The prefilter skips a piece on the score-only kernel's score and end; the traceback kernel
    reports the same, also on ties (the alignment that ends first) and on masked residues."""
    rng = random.Random(7)
    for trial in range(600):
        qlen = rng.choice([5, 12, 40, 120])
        q = ''.join(rng.choice(AA) for _ in range(qlen))
        if trial % 3 == 0:
            unit = q[rng.randrange(qlen // 2):]
            t = ''.join(rng.choice(AA + "*") for _ in range(rng.randrange(60))) + unit \
                + ''.join(rng.choice(AA) for _ in range(rng.randrange(60))) + unit
        else:
            t = ''.join(rng.choice(AA + "**") for _ in range(rng.randrange(20, 1_500)))
        best = sws.best_local_score(q, t)
        full = sws.smith_waterman_parasail(q, t)
        assert (best.score, best.end_ref) == (full['score'], full['end_ref'])


def test_a_kernel_disagreement_stops_the_child():
    best = types.SimpleNamespace(score=60, end_ref=10)
    with pytest.raises(RuntimeError, match="disagree"):
        sws.check_same_alignment(best, {'score': 61, 'end_ref': 10}, "a piece")
    sws.check_same_alignment(best, {'score': 60, 'end_ref': 10}, "a piece")


# ------------------------------------------------------------------ equivalence of the hit sets
def _window(rng, queries, n=90_000, copies=3):
    """Random DNA with mutated copies of each query, some reverse-complemented."""
    dna = ''.join(rng.choice("ACGT") for _ in range(n))
    for q in queries:
        for k in range(copies):
            mut = ''.join(c if rng.random() > 0.3 + 0.15 * k else rng.choice(AA) for c in q)
            gene = ''.join(CODON[a] for a in mut)
            if k % 2:
                gene = gene[::-1].translate(str.maketrans("ACGT", "TGCA"))
            pos = rng.randrange(0, len(dna) - len(gene))
            dna = dna[:pos] + gene + dna[pos + len(gene):]
    return dna


def _run_main(tmp_path, monkeypatch, queries, window, name, extra=(), max_cells=None, prefilter=True):
    q, t, o = tmp_path / f"{name}.faa", tmp_path / f"{name}.fna", tmp_path / f"{name}.m8"
    q.write_text(''.join(f">GOI_q{i}\n{s}\n" for i, s in enumerate(queries)))
    t.write_text(f">region_seq\n{window}\n")
    if max_cells is not None:
        monkeypatch.setattr(sws, "SW_MAX_CELLS", max_cells)
    if not prefilter:     # the traceback on every piece, as until 2026-10-28
        monkeypatch.setattr(sws, "min_reportable_score", lambda *a, **k: -10 ** 9)
    monkeypatch.setattr(sys, "argv", ["smith_waterman_search.py", "--query", str(q), "--target", str(t),
                                      "--output", str(o), "--min_score", "20", "--min_identity", "10.0",
                                      "--method", "parasail", *extra])
    sws.main()
    monkeypatch.undo()
    return o.read_text().splitlines()


@pytest.mark.parametrize("seed", [3, 11])
@pytest.mark.parametrize("alignments", ["1", "3"])
@pytest.mark.parametrize("max_cells", [None, 70 * 1_400])
def test_the_prefilter_keeps_exactly_the_hits_the_runner_keeps(tmp_path, monkeypatch, seed, alignments,
                                                               max_cells):
    """Planted copies (also across piece borders with small pieces), one or three alignments
    per piece: the lines with E <= 1 are byte-identical to a traceback on every piece, in the
    same order; only lines with E > 1, which the runner drops, are no longer written."""
    rng = random.Random(seed)
    queries = [''.join(rng.choice(AA) for _ in range(n)) for n in (70, 40, 150)]
    window = _window(rng, queries)
    extra = ("--max_alignments_per_piece", alignments)
    new = _run_main(tmp_path, monkeypatch, queries, window, "new", extra, max_cells)
    old = _run_main(tmp_path, monkeypatch, queries, window, "old", extra, max_cells, prefilter=False)
    assert _kept(new) == _kept(old)
    assert _kept(new), "the planted copies give kept hits"
    assert len(new) <= len(old)
    if alignments == "1":
        # every line the prefilter writes is a kept one or one raw score below the threshold
        # (E at most exp(lambda) = 1.31 times 1); the chance hits of the other pieces are gone
        assert all(float(r.split("\t")[10]) < 1.31 for r in new)
        assert len(new) < len(old)


def test_a_window_without_a_homolog_writes_no_line(tmp_path, monkeypatch):
    """Chance alignments (E > 1) were written for every frame and dropped by the runner; the
    prefilter does not align them at all."""
    rng = random.Random(43)
    window = ''.join(rng.choice("ACGT") for _ in range(30_000))
    query = ''.join(rng.choice(AA) for _ in range(70))
    old = _run_main(tmp_path, monkeypatch, [query], window, "old", prefilter=False)
    new = _run_main(tmp_path, monkeypatch, [query], window, "new")
    assert old and not _kept(old)
    assert new == []


# ------------------------------------------------------------------ failures end the child
def _files(tmp_path, window="ACGT" * 300):
    q, t = tmp_path / "q.faa", tmp_path / "t.fna"
    q.write_text(">GOI_q\nMKVLAAGIVGLLLAGCSSNKEEAKQWLDE\n")
    t.write_text(f">region_seq\n{window}\n")
    return q, t


@pytest.mark.parametrize("alignments", ["1", "3"])
def test_a_failed_alignment_ends_the_child_instead_of_skipping_the_piece(tmp_path, monkeypatch, alignments):
    q, t = _files(tmp_path, ''.join(CODON[a] for a in "MKVLAAGIVGLLLAGCSSNKEEAKQWLDE") * 3)

    def boom(*_a, **_k):
        raise MemoryError("no room for the traceback")
    monkeypatch.setattr(sws, "smith_waterman_parasail", boom)
    monkeypatch.setattr(sys, "argv", ["smith_waterman_search.py", "--query", str(q), "--target", str(t),
                                      "--output", str(tmp_path / "o.m8"), "--min_score", "20",
                                      "--method", "parasail", "--max_alignments_per_piece", alignments])
    with pytest.raises(MemoryError):
        sws.main()
    assert not (tmp_path / "o.m8").exists()


# ------------------------------------------------------------------ the guard (sequence_utils)
def test_run_guarded_returns_the_completed_process():
    r = su.run_guarded([sys.executable, "-c", "import sys; print('out'); sys.exit(3)"], 30, "a test")
    assert (r.returncode, r.stdout.strip()) == (3, "out")


def _grandchildren(tmp_path, marker):
    """A command whose two children run a script named after ``marker``, each with a child of its
    own (as ssearch36 under smith_waterman_search.py): the marker is in the command line of the
    command and of both children, so pgrep -f finds them while they live."""
    script = tmp_path / f"{marker}.sh"
    script.write_text("#!/bin/bash\nsleep 300\n")
    return ["bash", "-c", f"bash {script} & bash {script}; wait"]


def _alive(marker):
    return subprocess.run(["pgrep", "-f", marker], capture_output=True, text=True).stdout.split()


def test_the_guard_stops_the_whole_process_group(tmp_path):
    """A child that does not finish within the guard is killed with its own children;
    StepFailed names the place and the parameter. Killing the child alone would leave both
    scripts running."""
    marker = f"synvoy_guard_test_{os.getpid()}_{time.time_ns()}"
    seen = []
    watcher = threading.Timer(0.5, lambda: seen.extend(_alive(marker)))
    watcher.start()
    t0 = time.monotonic()
    with pytest.raises(su.StepFailed) as err:
        su.run_guarded(_grandchildren(tmp_path, marker), 1,
                       "[g_b3] Smith-Waterman search of the window", param="sw_timeout_seconds")
    watcher.join()
    assert time.monotonic() - t0 < 30
    assert "[g_b3] Smith-Waterman search of the window did not finish within 1 s (--sw_timeout_seconds)" \
        in str(err.value)
    assert len(seen) == 3                 # the command and its two scripts, found while they ran
    time.sleep(0.3)
    assert _alive(marker) == []


def test_an_interrupt_stops_the_whole_process_group_too(tmp_path, monkeypatch):
    """The child has a session of its own: an interrupt of the caller (Ctrl-C) must stop it and
    its children too."""
    marker = f"synvoy_guard_int_{os.getpid()}_{time.time_ns()}"
    real = subprocess.Popen.communicate
    seen = []

    def interrupted(self, *a, **k):
        if k.get("timeout"):
            time.sleep(0.5)
            seen.extend(_alive(marker))
            raise KeyboardInterrupt
        return real(self, *a, **k)
    monkeypatch.setattr(subprocess.Popen, "communicate", interrupted)
    with pytest.raises(KeyboardInterrupt):
        su.run_guarded(_grandchildren(tmp_path, marker), 60, "a test")
    assert len(seen) == 3
    time.sleep(0.3)
    assert _alive(marker) == []


def test_step_failed_is_not_caught_by_except_exception():
    """The search's fallbacks catch Exception; the guard's failure must pass through them."""
    caught = []
    with pytest.raises(su.StepFailed):
        try:
            raise su.StepFailed("guard")
        except Exception:          # the fallback of a window, a block, a gene model
            caught.append(True)
    assert caught == []


# ------------------------------------------------------------------ the runner's window search
def test_backend_info_names_what_runs():
    info = sws.backend_info("auto")
    assert info == {'backend': 'parasail', 'parasail_version': info['parasail_version'], 'simd': info['simd']}
    # the binding's version and the C library's, which holds the alignment code
    assert info['parasail_version'] == (f"{parasail.__version__} (libparasail "
                                        f"{'.'.join(str(v) for v in parasail.version())})")
    assert info['simd'] in sws.SIMD_ORDER + ('none',)
    assert sws.backend_info("ssearch36")['backend'] == 'parasail_best_hit_rounds'


import csv  # noqa: E402

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402


def _m8(path, rows):
    with open(path, "w") as fh:
        for q, aln, e in rows:
            fh.write(f"{q}\tregion_seq\t60.0\t{aln}\t0\t0\t1\t{aln}\t101\t{100 + 3 * aln}\t{e}\t60.0\n")


def _args(tmp_path, **over):
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
             aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
             aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6, min_identity=10.0,
             aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=True, deterministic_goi_search=True,
             sw_min_score=20, sw_min_identity=15.0, sw_method="parasail",
             sw_max_alignments_per_piece=1, sw_timeout_seconds=3600, genome_ctx={})
    a.update(over)
    return types.SimpleNamespace(**a)


def _window_search(tmp_path, monkeypatch, child):
    """run_augmented_search on a 2-kb window with MMseqs2 and tblastn faked; ``child(cmd,
    guard, what, param)`` stands in for the Smith-Waterman child under run_guarded."""
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")

    def fake_mmseqs(**k):
        _m8(k["aug_hits_file"], [("GOI_q", 40, "1e-20")])
        return True, ""

    def fake_run(cmd, *a, **k):
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr.shutil, "which", lambda name: None)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", fake_mmseqs)
    monkeypatch.setattr(isr, "run_guarded", child)
    return isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": "M" + "A" * 119}],
                                    "g1_b3", _args(tmp_path), "t1", 1, region_len=2_000)


def test_the_window_search_names_its_cells_and_seconds(tmp_path, monkeypatch, capsys):
    def child(cmd, guard, what, param=None, **k):
        assert guard == 3600 and param == "sw_timeout_seconds"
        assert "[g1_b3] Smith-Waterman search of the window (2,000 bp, 1 queries of 120 aa: " \
               f"{sws.smith_waterman_cells(2_000, [120]):,} cells)" == what
        _m8(cmd[cmd.index("--output") + 1], [("GOI_q", 30, "0.01")])
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")
    _window_search(tmp_path, monkeypatch, child)
    out = capsys.readouterr().out
    assert re.search(r"\[g1_b3\] Smith-Waterman: 479,520 cells in [0-9.]+ s, 1 hits kept\.", out)   # 120 x 3,996
    assert "[g1_b3] Smith-Waterman found 1 additional hits." in out


@pytest.mark.parametrize("rc,write,expect", [(1, True, r"failed \(exit 1\): .*broken piece"),
                                             (0, False, r"failed \(exit 0, no hits file\)")])
def test_a_failed_child_stops_the_genomes_search(tmp_path, monkeypatch, rc, write, expect):
    """A child that fails, or writes no hits file, was a warning (the window kept its other
    hits; any other exception dropped them uncounted); now it stops the genome's search."""
    def child(cmd, guard, what, param=None, **k):
        if write:
            _m8(cmd[cmd.index("--output") + 1], [])
        return types.SimpleNamespace(returncode=rc, stdout="", stderr="Traceback ...\nValueError: broken piece\n")
    with pytest.raises(isr.StepFailed, match=expect):
        _window_search(tmp_path, monkeypatch, child)


@pytest.mark.parametrize("rc,stderr,why", [(-9, "", "killed by signal 9"),
                                            (1, "Traceback ...\nMemoryError\n", "MemoryError")])
def test_a_child_killed_for_memory_says_so(tmp_path, monkeypatch, rc, stderr, why):
    """The wave exits 137 when every failure is this one, so Nextflow retries with more memory."""
    def child(cmd, guard, what, param=None, **k):
        return types.SimpleNamespace(returncode=rc, stdout="", stderr=stderr)
    with pytest.raises(isr.StepFailed, match=rf"ran out of memory \({why}\)"):
        _window_search(tmp_path, monkeypatch, child)
    src = open(isr.__file__).read()
    wave = src[src.index("Aborting iterative search due to per-genome processing errors."):]
    assert ("if all('Smith-Waterman' in msg and 'ran out of memory' in msg for _, msg in wave_errors):\n"
            "                sys.exit(137)\n            sys.exit(1)") in wave


def test_a_backend_that_cannot_run_stops_the_search_at_its_start():
    src = open(isr.__file__).read()
    block = src[src.index("sw_missing = None"):]
    block = block[:block.index("sys.exit(1)") + 11]
    assert "args.sw_method == 'parasail' and not HAS_PARASAIL" in block
    assert "args.sw_record['backend'] == 'ssearch36' and not shutil.which('ssearch36')" in block


def test_a_child_that_cannot_start_stops_the_genomes_search(tmp_path, monkeypatch):
    def child(cmd, guard, what, param=None, **k):
        raise FileNotFoundError("python3")
    with pytest.raises(isr.StepFailed, match="could not start"):
        _window_search(tmp_path, monkeypatch, child)


def test_a_genome_whose_search_stopped_fails_its_wave():
    """process_single_genome turns StepFailed into the genome's error (before its catch-all of
    Exception, which StepFailed, a BaseException, would pass anyway); the wave then exits 1."""
    src = open(isr.__file__).read()
    body = src[src.index("def process_single_genome("):]
    body = body[:body.index("\ndef ", 10)]
    handler = body[body.index("    except StepFailed as e:\n"):]
    assert handler.index("    except Exception as e:\n        error_message = str(e)") < 600
    assert "error_message = str(e)" in handler[:handler.index("    except Exception as e:")]
    wave = src[src.index("if wave_errors:"):]
    assert "sys.exit(1)" in wave[:wave.index("# Update DB after Wave")]


def test_the_backend_reaches_the_target_policy_and_the_report(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path), sw_record=sws.backend_info("auto"))
    isr._write_target_policy(args, "", "G1.fna", 70.0, 10)
    path = tmp_path / "hits" / "G1.fna.target_policy.tsv"
    row = next(csv.DictReader(open(path), delimiter="\t"))
    assert (row["sw_backend"], row["sw_parasail_version"]) == ("parasail", sws.backend_info("auto")["parasail_version"])
    assert gr.collect_sw_backend([str(path)]) == sws.backend_info("auto")
    assert "sw_backend" not in gr.collect_target_policy([str(path)])["G1"]

    other = tmp_path / "G2.fna.target_policy.tsv"
    other.write_text(path.read_text().replace("\tparasail\t", "\tssearch36\t"))
    assert len(gr.collect_sw_backend([str(path), str(other)])["mixed"]) == 2
    old = tmp_path / "old.target_policy.tsv"
    old.write_text("genome\tmedian_flank_identity\nG3.fna\t70.0\n")
    assert gr.collect_sw_backend([str(old)]) is None
