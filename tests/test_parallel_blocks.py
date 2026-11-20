"""The search blocks of one genome side by side (bin/synvoy_parallel.py, --block_workers): the
results, the genome context, the counters and the task's stdout/stderr must be the ones the
blocks give one after another, a failing block must fail the genome as it did, and a write the
merge has no rule for must stop the search instead of being lost.

The scenario runs in fresh interpreters (the blocks fork, and the output is compared as the task
writes it); its blocks use the runner's own writers (record_judge_numbers, count_limit, the
miniprot counter) with a model id two blocks share, so the repeat count crosses blocks."""
import ast
import json
import os
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import synvoy_parallel as sp  # noqa: E402

_SCENARIO = r"""
import json, os, sys
sys.path.insert(0, sys.argv[1])
work, workers = sys.argv[2], int(sys.argv[3])
fail_at, kill_at, rogue_at = (int(x) if x != '-' else None for x in sys.argv[4:7])
from types import SimpleNamespace
import iterative_search_runner as isr
import annotate_goi_exons as age
import synvoy_parallel as sp

args = SimpleNamespace(genome_ctx={'places': {'g1': 2}, 'genome_len': 1000, 'decoy_lines': [],
                                   'block_esyn': {}})
age.LIMIT_COUNTS.clear()
blocks = [{'chrom': 'c', 'start': i * 100, 'end': i * 100 + 50} for i in range(7)]


def run_block(i, block):
    ctx = args.genome_ctx
    print(f'stdout of block {i}', flush=(i % 2 == 0))
    isr.logger.info(f'log of block {i}')
    ctx.setdefault('block_context', {})[i] = (i % 3 == 0)
    ctx.setdefault('block_windows', []).append(['c', block['start'], block['end'], i])
    ctx.setdefault('goi_cores', []).append(['c', block['start'], block['end']])
    ctx.setdefault('block_esyn', {})[i] = 0.01 * i
    if i % 2:
        ctx.setdefault('decoy_lines', []).extend([f'decoy {i} a', f'decoy {i} b'])
    if i == 4:
        ctx.setdefault('reduced_search', {}).setdefault('windows', []).append({'window': f'b{i}'})
    isr.count_limit(args, 'goi_core_filter', i)
    age._count('miniprot_no_alignment_goi', 1)
    cands = [{'gff': [f'c\tsrc\tmRNA\t1\t2\t.\t+\t.\tID=M{i};SynVoy_Parent=P'], 'judge_numbers': {'E': i}}]
    if i in (2, 5):
        cands.append({'gff': [f'c\tsrc\tmRNA\t1\t2\t.\t+\t.\tID=SHARED;SynVoy_Parent=P{i}'],
                      'judge_numbers': {'E': -i}})
    ctx['judge_repeats'] = ctx.get('judge_repeats', 0) + isr.record_judge_numbers(
        ctx.setdefault('judge_numbers', {}), cands, ctx.setdefault('judge_writers', {}))
    if i == fail_at:
        raise isr.StepFailed(f'block {i} failed')
    if i == kill_at:
        os.kill(os.getpid(), 9)
    if i == rogue_at:
        ctx['places']['g1'] = 99
    return [{'id': f'gene{i}'}], [f'gff line {i}']


out = []
try:
    for res in sp.search_blocks(blocks, run_block, args, workers, work_dir=work,
                                **isr.block_worker_state()):
        out.append(res)
        isr.logger.info('consumer took ' + res[0][0]['id'])
except BaseException as exc:
    out.append(['error', type(exc).__name__, str(exc)])
print('STATE ' + json.dumps({'ctx': args.genome_ctx, 'limits': age.LIMIT_COUNTS, 'out': out,
                             'work_dir_left': os.path.exists(work)}, sort_keys=True), flush=True)
"""

_STAMP = re.compile(r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d ", re.M)


def _run(tmp_path, workers, fail_at=None, kill_at=None, rogue_at=None):
    script = tmp_path / "scenario.py"
    script.write_text(_SCENARIO)
    work = tmp_path / f"work_{workers}_{fail_at}_{kill_at}_{rogue_at}"
    env = dict(os.environ, PYTHONHASHSEED="0")
    proc = subprocess.run(
        [sys.executable, str(script), BIN, str(work), str(workers)]
        + ['-' if x is None else str(x) for x in (fail_at, kill_at, rogue_at)],
        capture_output=True, text=True, env=env, timeout=120)
    assert proc.returncode == 0, proc.stderr
    lines = proc.stdout.splitlines()
    state = json.loads(lines[-1][len("STATE "):])
    return state, "\n".join(lines[:-1]), _STAMP.sub("", proc.stderr)


def test_side_by_side_gives_what_one_after_another_gives(tmp_path):
    serial = _run(tmp_path, 1)
    for workers in (2, 3, 7):
        parallel = _run(tmp_path, workers)
        assert parallel[0] == serial[0]          # results, genome context, counters
        assert parallel[1] == serial[1]          # stdout, line for line
        assert parallel[2] == serial[2]          # stderr (the log) without its timestamps
    state = serial[0]
    assert state["ctx"]["judge_repeats"] == 1    # SHARED written by blocks 2 and 5
    assert state["ctx"]["judge_numbers"]["SHARED"] == {"E": -5}
    assert state["ctx"]["judge_writers"]["SHARED"] == ["P2", "P5"]
    assert state["ctx"]["limits"] == {"goi_core_filter": 21}
    assert state["limits"] == {"miniprot_no_alignment_goi": 7}
    assert state["ctx"]["reduced_search"] == {"windows": [{"window": "b4"}]}
    assert [r[0][0]["id"] for r in state["out"]] == [f"gene{i}" for i in range(7)]
    assert not state["work_dir_left"]
    assert "log of block 6" in serial[2] and "stdout of block 6" in serial[1]


def test_a_failing_block_fails_the_genome_as_before(tmp_path):
    serial = _run(tmp_path, 1, fail_at=3)
    parallel = _run(tmp_path, 3, fail_at=3)
    assert serial[0]["out"][-1] == ["error", "StepFailed", "block 3 failed"]
    assert parallel[0] == serial[0]              # blocks 0-3 taken back, 4-6 dropped
    assert parallel[1] == serial[1]
    assert parallel[2] == serial[2]
    assert not parallel[0]["work_dir_left"]


def test_a_killed_worker_is_a_memory_kill(tmp_path):
    state, _out, _err = _run(tmp_path, 3, kill_at=2)
    kind, message = state["out"][-1][1:]
    assert kind == "MemoryKilled" and "ran out of memory" in message
    # The pool breaks for every block not finished when the worker dies, so how many of blocks 0
    # and 1 are taken back first depends on timing; never one out of order, never a later one.
    taken = [r[0][0]["id"] for r in state["out"][:-1]]
    assert taken == ["gene0", "gene1"][:len(taken)]
    assert not state["work_dir_left"]


def test_a_write_without_a_rule_stops_the_search(tmp_path):
    state, _out, _err = _run(tmp_path, 3, rogue_at=1)
    kind, message = state["out"][-1][1:]
    assert kind == "UnmergedBlockChange" and "places" in message


def test_worker_count():
    assert sp.block_worker_count("auto", 4, 10) == (4, "auto: 4 core(s) for this genome")
    assert sp.block_worker_count("auto", 4, 3)[0] == 3
    assert sp.block_worker_count("auto", 4, 1)[0] == 1
    assert sp.block_worker_count("2", 4, 10)[0] == 2
    assert sp.block_worker_count(8, 4, 10)[0] == 4          # never more than the genome's cores
    assert sp.block_worker_count("1", 4, 10)[0] == 1
    assert sp.block_worker_count("auto", 4, 10, serial_reason="plm") == (1, "plm")
    assert sp.block_worker_count("auto", 4, 10, auto_serial_reason="threads") == (1, "threads")
    assert sp.block_worker_count("3", 4, 10, auto_serial_reason="threads")[0] == 3
    for bad in ("0", "-1", "many", "2.5"):
        with pytest.raises(ValueError):
            sp.block_worker_count(bad, 4, 10)


def test_no_workers_off_linux(monkeypatch):
    monkeypatch.setattr(sp.sys, "platform", "darwin")
    workers, why = sp.block_worker_count("auto", 4, 10)
    assert workers == 1 and "Linux only" in why


def test_serial_reasons_of_the_runner():
    args = SimpleNamespace(enable_plm_search=False, enable_structural_search=False,
                           enable_structural_discovery=False, deterministic_goi_search=True,
                           enable_smith_waterman=True, sw_record={"backend": "parasail"})
    assert isr.block_workers_serial_reason(args) is None
    assert isr.block_workers_auto_serial_reason(args) is None
    assert "plm" in isr.block_workers_serial_reason(SimpleNamespace(**dict(vars(args), enable_plm_search=True)))
    assert isr.block_workers_auto_serial_reason(
        SimpleNamespace(**dict(vars(args), deterministic_goi_search=False)))
    assert isr.block_workers_auto_serial_reason(
        SimpleNamespace(**dict(vars(args), sw_record={"backend": "ssearch36"})))


def test_the_runner_passes_the_setting():
    parser_src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    assert '"--block_workers", type=str, default="auto"' in parser_src
    nf = open(os.path.join(BIN, "..", "modules", "iterative_search.nf")).read()
    # the setting on the first attempt; a retried attempt (after a memory kill) one after another
    assert "def block_workers = task.attempt > 1 ? 1 : params.iterative_search_block_workers" in nf
    assert "--block_workers ${block_workers} \\\\" in nf
    cfg = open(os.path.join(BIN, "..", "nextflow.config")).read()
    assert re.search(r"(?m)^\s*iterative_search_block_workers\s*=\s*'auto'", cfg)


def test_delta_rules():
    s0 = {"places": {"a": 1}, "decoy_lines": ["x"], "limits": {"k": 2}}
    ctx = {"places": {"a": 1}, "decoy_lines": ["x", "y"], "limits": {"k": 5, "n": 1},
           "block_slot": {3: (None, None)}}
    delta = sp.ctx_delta(s0, ctx)
    assert delta == {"decoy_lines": (False, ["y"]), "limits": (False, [("k", 3), ("n", 1)]),
                     "block_slot": (True, [(3, (None, None))])}
    merged = {"places": {"a": 1}, "decoy_lines": ["x", "w"], "limits": {"k": 4}}
    sp.apply_ctx_delta(merged, delta, s0)
    assert merged == {"places": {"a": 1}, "decoy_lines": ["x", "w", "y"],
                      "limits": {"k": 7, "n": 1}, "block_slot": {3: (None, None)}}
    with pytest.raises(sp.UnmergedBlockChange):
        sp.ctx_delta(s0, dict(ctx, decoy_lines=["z", "y"]))       # rewrote, not appended
    with pytest.raises(sp.UnmergedBlockChange):
        sp.ctx_delta(s0, dict(ctx, limits={}))                    # removed an entry
    with pytest.raises(sp.UnmergedBlockChange):
        sp.ctx_delta(s0, dict(ctx, new_key=1))                    # no rule


# Keys of args.genome_ctx that only the genome's own code writes, before its blocks or after
# them (process_single_genome, apply_rescue, genome_wide_goi_pass): a block never changes them.
_OUTSIDE_BLOCKS = {
    "axis", "decoy_counts", "decoy_gw_counts", "decoy_gw_parents", "flank_median", "genome_len",
    "places", "in_zone_hits", "pre_rescue_zones", "rescue_dropped",
    "rescue_dropped_outside_zones", "rescue_record", "rescue_refused",
}


def _written_ctx_keys():
    """Every constant key the runner stores into a genome context (ctx[...] = ...,
    ctx.setdefault(...), and the two keys of the conditional setdefault of a block's zone)."""
    tree = ast.parse(open(os.path.join(BIN, "iterative_search_runner.py")).read())

    def is_ctx(n):
        if isinstance(n, ast.Name):
            return n.id in ("ctx", "_ctx")
        if isinstance(n, ast.Attribute):
            return n.attr == "genome_ctx"
        if isinstance(n, ast.BoolOp):
            return any(is_ctx(v) for v in n.values)
        return (isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "getattr"
                and len(n.args) >= 2 and isinstance(n.args[1], ast.Constant)
                and n.args[1].value == "genome_ctx")

    keys = set()
    for n in ast.walk(tree):
        if (isinstance(n, ast.Subscript) and isinstance(n.ctx, ast.Store) and is_ctx(n.value)
                and isinstance(n.slice, ast.Constant)):
            keys.add(n.slice.value)
        if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                and n.func.attr in ("setdefault", "update", "pop") and is_ctx(n.func.value) and n.args):
            arg = n.args[0]
            if isinstance(arg, ast.Constant):
                keys.add(arg.value)
            elif isinstance(arg, ast.IfExp):
                keys.update(x.value for x in (arg.body, arg.orelse) if isinstance(x, ast.Constant))
    return keys


def test_every_key_the_runner_writes_is_classified():
    keys = _written_ctx_keys()
    assert {"block_context", "goi_cores", "genome_wide_windows", "judge_numbers"} <= keys
    unclassified = keys - set(sp.CTX_RULES) - _OUTSIDE_BLOCKS
    assert not unclassified, (
        f"genome_ctx keys {sorted(unclassified)} are new: a key a block writes needs a rule in "
        "synvoy_parallel.CTX_RULES, one written only before or after the blocks goes to "
        "_OUTSIDE_BLOCKS here")
