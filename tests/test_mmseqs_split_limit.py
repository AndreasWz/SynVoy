"""A split limit too small for a genome's database (MMseqs2: "Cannot fit databases into ...").

The genome ladder keeps the split limit on every rung and lowers threads and sensitivity, which do
not change what has to fit, so on that failure it stops after the first attempt and the genome is
recorded as not searched with the reason mmseqs_split_memory_limit_too_small and the settings that
fix it. A memory kill still takes the lower rungs. A window's ladder goes on only to a rung with a
larger limit (its 2G rung, below a 2G limit). The searches run MMseqs2 at -v 1 or more, so the
cause is printed (at -v 0 MMseqs2 says only "Prefilter died"), and a failed attempt's folder is
removed before the next one.
"""
import json
import os
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import iterative_search_runner as isr  # noqa: E402
import generate_report as gr  # noqa: E402
import sequence_utils as su  # noqa: E402

# MMseqs2 18-8cc5c's streams (measured, -v 1 and -v 3): its own error on stderr, the workflow's
# lines on stdout. A memory kill: the shell's "Killed" on stderr.
SPLIT_ERR = "Cannot fit databases into 2G. Please use a computer with more main memory.\n"
KILLED_ERR = "Killed\n"
STEP_DIED = "Error: Prefilter died\nError: Search step died\n"
# The details text of such failures, as _format_mmseqs_failure joins the two streams.
SPLIT = SPLIT_ERR + STEP_DIED
KILLED = KILLED_ERR + STEP_DIED


def _proc(rc, stderr="", stdout=""):
    return types.SimpleNamespace(returncode=rc, stdout=stdout, stderr=stderr)


def _args(tmp_path, **kw):
    a = dict(output_dir=str(tmp_path), mmseqs_sens=9.5, mmseqs_split_memory_limit="2G", evalue=0.01,
             mmseqs_verbosity=0, aug_relaxed_evalue_cap=10.0, aug_relaxed_evalue_mult=1000.0,
             window_mmseqs_min_length=30, genome_ctx={})
    a.update(kw)
    return types.SimpleNamespace(**a)


def _fake_mmseqs(stderr_of_failure, fail_first, record, folders=None, stdout_of_failure=STEP_DIED):
    """MMseqs2 whose first ``fail_first`` searches fail with ``stderr_of_failure`` on stderr and
    the workflow's lines on stdout; records each search command and, with ``folders``, which
    attempt folders exist when a search starts."""
    calls = {"search": 0}

    def run(cmd, *a, **k):
        if cmd[1] == "search":
            calls["search"] += 1
            record.append(list(cmd))
            if folders is not None:
                parent = os.path.dirname(cmd[5])
                folders.append(sorted(p for p in os.listdir(parent) if p.startswith("attempt_")))
            # an attempt writes its temporary files
            with open(os.path.join(cmd[5], "prefilter.tmp"), "w") as fh:
                fh.write("x")
            if calls["search"] <= fail_first:
                return _proc(1, stderr_of_failure, stdout_of_failure)
        if cmd[1] == "convertalis":
            open(cmd[5], "w").close()
        return _proc(0)
    return run


def _genome_ladder(tmp_path, args):
    return isr.run_mmseqs_easy_search_with_retries("q.faa", "g.fna", str(tmp_path / "h.m8"),
                                                   str(tmp_path / "t"), args, 4, "G", used={})


# --------------------------------------------------------------------------- the genome ladder

def test_the_failure_kind():
    assert isr._is_split_limit_failure(SPLIT)
    assert isr._is_mmseqs_resource_failure(SPLIT)
    assert not isr._is_split_limit_failure(KILLED)
    assert not isr._is_split_limit_failure("")


def test_a_split_limit_failure_stops_the_genome_ladder_after_one_attempt(tmp_path, monkeypatch):
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(SPLIT_ERR, 9, cmds))
    ok, details, resource = _genome_ladder(tmp_path, _args(tmp_path))
    assert (ok, resource) == (False, True)
    assert len(cmds) == 1
    assert isr._is_split_limit_failure(details)


def test_a_memory_kill_still_takes_the_lower_rungs(tmp_path, monkeypatch):
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(KILLED_ERR, 1, cmds))
    ok, _details, _resource = _genome_ladder(tmp_path, _args(tmp_path))
    assert ok and len(cmds) == 2
    assert cmds[1][cmds[1].index("-s") + 1] == "7.0"


def test_a_failed_attempts_folder_goes_before_the_next_attempt(tmp_path, monkeypatch):
    cmds, folders = [], []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(KILLED_ERR, 2, cmds, folders))
    ok, _d, _r = _genome_ladder(tmp_path, _args(tmp_path))
    assert ok and len(cmds) == 3
    assert folders == [["attempt_1"], ["attempt_2"], ["attempt_3"]]
    # the attempt that ran keeps its folder (the genome's end removes it)
    assert sorted(p for p in os.listdir(tmp_path / "t") if p.startswith("attempt_")) == ["attempt_3"]


def test_the_searches_print_mmseqs2s_errors(tmp_path, monkeypatch):
    """-v 0 (accepted with a start warning) becomes 1 for the searches: errors only."""
    for level, expect in ((0, "1"), (1, "1"), (3, "3")):
        cmds = []
        monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(KILLED_ERR, 0, cmds))
        _genome_ladder(tmp_path, _args(tmp_path, mmseqs_verbosity=level))
        assert cmds[0][cmds[0].index("-v") + 1] == expect
        cmds.clear()
        isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"), str(tmp_path / "w"),
                                              _args(tmp_path, mmseqs_verbosity=level), 1, "G_b3")
        assert cmds[0][cmds[0].index("-v") + 1] == expect


# --------------------------------------------------------------------------- the genome is not searched

def _search_args(tmp_path):
    (tmp_path / "hits").mkdir(exist_ok=True)
    return types.SimpleNamespace(output_dir=str(tmp_path), min_identity=20.0, min_length=10, evalue=1.0,
                                 mmseqs_split_memory_limit="2G")


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return str(path)


def test_the_genome_is_recorded_as_not_searched_with_the_settings_that_fix_it(tmp_path, monkeypatch):
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", "ACGT" * 50)])
    db = _write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV"), ("GOI_x", "MKT")])
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (False, SPLIT.strip(), True))
    name, genes, extras, error = isr.process_single_genome(genome, db, _search_args(tmp_path), None, "locus_1_", 1)
    assert (name, genes, extras, error) == ("Target.fna", [], [], None)
    rows = (tmp_path / "hits" / "locus_1_Target.fna.not_searched.tsv").read_text().splitlines()
    assert rows[0].split("\t") == isr.NOT_SEARCHED_COLUMNS
    genome_col, reason, details = rows[1].split("\t")
    assert (genome_col, reason) == ("Target.fna", "mmseqs_split_memory_limit_too_small")
    assert "--mmseqs_split_memory_limit" in details and "--iterative_search_memory" in details
    assert "Cannot fit databases into 2G" in details


def test_a_memory_kill_keeps_its_reason(tmp_path, monkeypatch):
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", "ACGT" * 50)])
    db = _write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV")])
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (False, KILLED.strip(), True))
    isr.process_single_genome(genome, db, _search_args(tmp_path), None, "locus_1_", 1)
    row = (tmp_path / "hits" / "locus_1_Target.fna.not_searched.tsv").read_text().splitlines()[1]
    assert row.split("\t")[1] == "mmseqs_resource_failure"


def test_the_report_names_the_settings(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path / "search"))
    isr._write_not_searched(args, "locus_1_", "Big_genome.fna", "mmseqs_split_memory_limit_too_small", SPLIT)
    results = tmp_path / "results"
    (results / "regions").mkdir(parents=True)
    (results / "hits").mkdir()
    name = "locus_1_Big_genome.fna.not_searched.tsv"
    (results / "hits" / name).write_text((tmp_path / "search" / "hits" / name).read_text())
    (results / "regions" / "locus_1__Other.fna.gff").write_text("##gff-version 3\n")
    report = gr.build_report(str(results))
    assert report["summary"]["genomes_not_searched"] == {"Big_genome": ["mmseqs_split_memory_limit_too_small"]}
    advice = [w for w in report["warnings"] if "split limit" in w]
    assert len(advice) == 1
    assert "--mmseqs_split_memory_limit" in advice[0] and "--iterative_search_memory" in advice[0]


# --------------------------------------------------------------------------- the windows

def test_a_split_limit_failure_stops_the_window_ladder(tmp_path, monkeypatch):
    """At 2G every rung of the window ladder has the same limit (its 2G rung included): one attempt."""
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(SPLIT_ERR, 9, cmds))
    ok, details = isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"),
                                                        str(tmp_path / "w"), _args(tmp_path), 1, "G_b3")
    assert not ok and len(cmds) == 1 and isr._is_split_limit_failure(details)
    assert not [p for p in os.listdir(tmp_path / "w") if p.startswith("attempt_")]


def test_the_window_ladder_still_tries_a_larger_split_limit(tmp_path, monkeypatch):
    """Below 2G the window ladder's 2G rung is a larger limit: the rungs with the failed limit are
    skipped, that one runs."""
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(SPLIT_ERR, 1, cmds))
    args = _args(tmp_path, mmseqs_split_memory_limit="1G")
    ok, _details = isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"),
                                                         str(tmp_path / "w"), args, 1, "G_b3")
    assert ok
    assert [c[c.index("--split-memory-limit") + 1] for c in cmds] == ["1G", "2G"]
    assert args.genome_ctx["limits"]["mmseqs_window_lower_rung"] == 1


def test_split_limits_in_bytes():
    assert isr._split_limit_bytes("2G") == 2 * 1024 ** 3
    assert isr._split_limit_bytes("1500M") < isr._split_limit_bytes("2G") < isr._split_limit_bytes("8g")
    assert isr._split_limit_bytes("0") == float("inf")


def test_a_window_that_cannot_fit_stops_the_search_with_the_settings(tmp_path, monkeypatch):
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", lambda **k: (False, SPLIT.strip()))
    monkeypatch.setattr(isr.subprocess, "run", lambda *a, **k: _proc(0))
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6,
             min_identity=10.0, aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=False, deterministic_goi_search=True, genome_ctx={})
    with pytest.raises(su.StepFailed) as stop:
        isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": "M" + "A" * 119}], "G_b3",
                                 types.SimpleNamespace(**a), "t1", 1, region_len=2000)
    assert "split limit is too small" in str(stop.value)
    assert "--mmseqs_split_memory_limit" in str(stop.value) and "--iterative_search_memory" in str(stop.value)


# --------------------------------------------------------------------------- both streams (review round 2)

def test_both_streams_reach_the_details_and_the_kind():
    """MMseqs2 names the cause on stderr and the step on stdout: both are in the details, and
    the cause decides the kind."""
    details, resource = isr._mmseqs_search_failure(_proc(1, SPLIT_ERR, STEP_DIED))
    assert resource and isr._is_split_limit_failure(details)
    assert "Cannot fit databases into 2G" in details and "Prefilter died" in details
    # the cause on stdout (another MMseqs2 build) is found too
    details, resource = isr._mmseqs_search_failure(_proc(1, "", SPLIT_ERR + STEP_DIED))
    assert resource and isr._is_split_limit_failure(details)


@pytest.mark.parametrize("message", [
    "Not enough memory to keep dbreader/write in memory!",
    "Can not allocate entries memory in IndexTable",
    "Allocation error : not enough memory",
    "Failed to mmap memory dataSize=1234 File=x",
    "mem_align could not allocate memory.",
    "Write buffer could not be allocated (bufferSize=1234)",
    "Killed",
])
def test_mmseqs2s_memory_messages_are_resource_failures(message):
    details, resource = isr._mmseqs_search_failure(_proc(1, message + "\n", STEP_DIED))
    assert resource and not isr._is_split_limit_failure(details)


# MMseqs2 18-8cc5c under an address-space limit (measured at 400M, 600M and 1300M): its message on
# stderr ends in a NUL byte.
MEM_ALIGN_ERR = "mem_align could not allocate memory.\n\x00"


def test_a_memory_message_with_a_nul_byte_records_the_genome_without_it(tmp_path, monkeypatch):
    """Review round 3: "mem_align could not allocate memory." stopped the run after one attempt
    (no listed message, stderr not empty). It takes the three rungs now and the genome is recorded
    as not searched, and the NUL byte that ends the message is not written to the record: Python
    3.10's csv reader refuses a line that holds one, and the report would then leave the genome out
    of genomes_not_searched."""
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", "ACGT" * 50)])
    db = _write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV"), ("GOI_x", "MKT")])
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(MEM_ALIGN_ERR, 9, cmds))
    args = _search_args(tmp_path)
    for k, v in dict(mmseqs_sens=7.5, mmseqs_verbosity=1).items():
        setattr(args, k, v)
    assert isr.process_single_genome(genome, db, args, None, "locus_1_", 1) == ("Target.fna", [], [], None)
    assert len(cmds) == 3
    path = tmp_path / "hits" / "locus_1_Target.fna.not_searched.tsv"
    text = path.read_text()
    assert "\x00" not in text and "mem_align could not allocate memory." in text
    assert text.splitlines()[1].split("\t")[1] == "mmseqs_resource_failure"
    assert gr.collect_not_searched([str(path)]) == {"Target": ["mmseqs_resource_failure"]}
    assert "\x00" not in isr._tail_lines("a\x00\nb\n\x00")


def test_a_step_that_died_without_a_message_is_a_resource_failure():
    _details, resource = isr._mmseqs_search_failure(_proc(1, "", STEP_DIED))
    assert resource


def test_a_step_that_died_after_another_error_is_not_one():
    """A file that cannot be written is no case for the lower rungs: the search stops with it."""
    details, resource = isr._mmseqs_search_failure(
        _proc(1, "Can not write to data file tmp/pref_0.0\n", STEP_DIED))
    assert not resource and "Can not write to data file" in details


def test_the_genome_is_recorded_as_not_searched_through_the_ladder(tmp_path, monkeypatch):
    """The ladder's own return, from MMseqs2's real stream layout, reaches the record."""
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", "ACGT" * 50)])
    db = _write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", "MKV"), ("GOI_x", "MKT")])
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(SPLIT_ERR, 9, cmds))
    args = _search_args(tmp_path)
    for k, v in dict(mmseqs_sens=7.5, mmseqs_verbosity=1).items():
        setattr(args, k, v)
    assert isr.process_single_genome(genome, db, args, None, "locus_1_", 1) == ("Target.fna", [], [], None)
    assert len(cmds) == 1
    row = (tmp_path / "hits" / "locus_1_Target.fna.not_searched.tsv").read_text().splitlines()[1]
    assert row.split("\t")[1] == "mmseqs_split_memory_limit_too_small"


def test_the_wave_hint_says_raise_for_a_split_limit(tmp_path):
    """The hint of a failed wave does not tell a split limit that is too small to go lower."""
    import inspect
    src = inspect.getsource(isr)
    hint = src[src.index("split_errors = ["):src.index('"    --mmseqs_split_memory_limit <lower value>')]
    assert "_is_split_limit_failure(m)" in hint and "not _is_split_limit_failure(m)" in hint
    assert "Raise\\n" in hint


# --------------------------------------------------------------------------- the real MMseqs2

def _random_dna(n, seed):
    import random
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _random_protein(n, seed):
    import random
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


@pytest.mark.skipif(not __import__("shutil").which("mmseqs"), reason="mmseqs not on PATH")
def test_the_real_mmseqs2_split_limit_failure_is_recorded(tmp_path):
    """MMseqs2 cannot fit a 200 kb target into 100M (measured with 18-8cc5c: under 1 s). The
    genome ladder stops after one attempt and the genome is recorded with the split-limit reason."""
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", _random_dna(200_000, 1))])
    db = _write_fasta(tmp_path / "db.faa", [("gene-a|exon_1", _random_protein(150, 2)),
                                            ("GOI_x", _random_protein(120, 3))])
    args = types.SimpleNamespace(mmseqs_sens=7.5, mmseqs_split_memory_limit="100M", evalue=0.01,
                                 mmseqs_verbosity=1, mmseqs_max_seqs=300)
    ok, details, resource = isr.run_mmseqs_easy_search_with_retries(
        db, genome, str(tmp_path / "h.m8"), str(tmp_path / "t"), args, 1, "Target.fna", used={})
    if ok:
        pytest.skip("this MMseqs2 fits the target into 100M")
    assert resource and isr._is_split_limit_failure(details)

    s_args = _search_args(tmp_path)
    for k, v in vars(args).items():
        setattr(s_args, k, v)
    assert isr.process_single_genome(genome, db, s_args, None, "locus_1_", 1) == ("Target.fna", [], [], None)
    row = (tmp_path / "hits" / "locus_1_Target.fna.not_searched.tsv").read_text().splitlines()[1]
    assert row.split("\t")[1] == "mmseqs_split_memory_limit_too_small"


@pytest.mark.skipif(not __import__("shutil").which("mmseqs"), reason="mmseqs not on PATH")
def test_the_real_mmseqs2_window_ladder_goes_on_to_its_2g_rung(tmp_path, monkeypatch):
    """A window that does not fit into 1G (a 200 kb window, MMseqs2 18-8cc5c) is searched on the
    window ladder's 2G rung."""
    window = _write_fasta(tmp_path / "window.fna", [("w1", _random_dna(200_000, 4))])
    query = _write_fasta(tmp_path / "v.faa", [("GOI_x", _random_protein(120, 5))])
    args = _args(tmp_path, mmseqs_split_memory_limit="1G", mmseqs_sens=7.5, mmseqs_verbosity=1,
                 mmseqs_max_seqs=300)
    first = []
    real_run = subprocess.run

    def run(cmd, *a, **k):
        proc = real_run(cmd, *a, **k)
        if cmd[1] == "search":
            first.append((cmd[cmd.index("--split-memory-limit") + 1], proc.returncode, proc.stderr))
        return proc
    monkeypatch.setattr(isr.subprocess, "run", run)
    ok, details = isr.run_augmented_mmseqs_with_retries(query, window, str(tmp_path / "a.m8"),
                                                        str(tmp_path / "w"), args, 1, "G_b0")
    if first and first[0][1] == 0:
        pytest.skip("this MMseqs2 fits the window into 1G")
    assert ok, details
    assert [f[0] for f in first] == ["1G", "2G"]
    assert isr._is_split_limit_failure(first[0][2])
