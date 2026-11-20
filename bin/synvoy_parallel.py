#!/usr/bin/env python3
"""The search blocks of one genome searched side by side (import-only helper of
iterative_search_runner.py). A speed-up only: the calls, records, decoys and counts are the
ones the blocks give one after another.

process_single_genome searches the synteny blocks of a genome one after another, and a block
uses about one core (its MMseqs2 window search runs on one thread, tblastn, the parasail
Smith-Waterman search and miniprot too), while the genome's worker holds ``threads_per_job``
cores. A block's search reads only what was there before the first block ran: the hits, the
genome, the query set, ``args`` and what the genome set up in ``args.genome_ctx`` (its places,
its length, the rescue zones, the axis). It writes only entries of its own there: its number
in block_context, block_esyn, block_rows and block_slot, its window and GOI zone at the end of
block_windows and goi_cores, its decoy rows, the judge numbers of its models, its limit
counts and the windows that ran on a lower MMseqs2 rung (checked by instrumenting every block
of the Candida gate genes, 2026-11-19, and by reading the code).

So the blocks run in worker processes forked from the genome's worker (they share the genome
and the hits without a copy). Each block starts from the state before the first block: the
genome context and the counters below are reset before each one. Their results and changes
are taken back in block order, as if the blocks had run one after another:

  - lists are extended (block_windows, goi_cores, genome_wide_windows, decoy_lines and
    reduced_search['windows']),
  - per-block entries are set (block_context, block_esyn, block_rows, block_slot,
    genome_wide_outcome),
  - counts are added (limits, genome_wide_dropped, and the module counters the runner
    names: annotate_goi_exons.LIMIT_COUNTS, the blastp and family-test failures),
  - the judge-number table is replayed: a model id an earlier block wrote counts as a repeat
    and the later write wins (record_judge_numbers), the writers of an id are kept in order.

A change to any other key of the genome context, or a removed entry, stops the genome's
search (UnmergedBlockChange): a new kind of write gets its rule here before it can run side
by side. Each block's stdout and stderr are kept in files and copied out in block order, so
the task's .command.out and .command.err hold the lines of the one-after-another run in the
same order (timestamps aside; a warning a process gives once can come once per worker). A
block that fails raises its own exception after the blocks before it are taken back, as it
would have one after another; the blocks after it are dropped. A worker killed by a signal
raises MemoryKilled (exit 137, which Nextflow retries with more memory; a retried attempt
passes --block_workers 1, modules/iterative_search.nf, so it searches one block at a time).
"""
import concurrent.futures
import contextlib
import copy
import gc
import multiprocessing
import os
import pickle
import shutil
import sys
import threading
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple

try:
    from sequence_utils import MemoryKilled
except ImportError:  # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import MemoryKilled

# How each key of args.genome_ctx that a block may change is taken back (module docstring).
APPEND, PER_BLOCK, COUNT, JUDGE, WINDOWS = 'append', 'per_block', 'count', 'judge', 'windows'
CTX_RULES = {
    'block_windows': APPEND, 'goi_cores': APPEND, 'genome_wide_windows': APPEND,
    'decoy_lines': APPEND,
    'block_context': PER_BLOCK, 'block_esyn': PER_BLOCK, 'block_rows': PER_BLOCK,
    'block_slot': PER_BLOCK, 'genome_wide_outcome': PER_BLOCK,
    'limits': COUNT, 'genome_wide_dropped': COUNT,
    'judge_numbers': JUDGE, 'judge_writers': JUDGE, 'judge_repeats': JUDGE,
    'reduced_search': WINDOWS,
}

_MISSING = object()


class UnmergedBlockChange(RuntimeError):
    """A block changed the genome context in a way no rule of CTX_RULES takes back."""


def block_worker_count(setting: Any, threads_per_job: int, n_blocks: int,
                       serial_reason: Optional[str] = None,
                       auto_serial_reason: Optional[str] = None) -> Tuple[int, str]:
    """(workers, why) for the blocks of one genome. ``setting``: 'auto' (the genome's cores,
    ``threads_per_job``) or a number, the most workers (never more than the genome's cores); 1
    searches the blocks one after another, as before. ``serial_reason``: why this run keeps
    them one after another whatever the setting; ``auto_serial_reason``: why 'auto' does (the
    caller's checks)."""
    cores = max(1, int(threads_per_job or 1))
    text = str(setting if setting is not None else 'auto').strip().lower()
    if text in ('', 'auto'):
        if auto_serial_reason:
            return 1, auto_serial_reason
        wanted, why = cores, f'auto: {cores} core(s) for this genome'
    else:
        if not text.isdigit() or int(text) < 1:
            raise ValueError(f"--block_workers must be 'auto' or a whole number >= 1, not {setting!r}")
        wanted = int(text)
        why = f'--block_workers {wanted}, {cores} core(s) for this genome'
        wanted = min(wanted, cores)
    if serial_reason:
        return 1, serial_reason
    if n_blocks < 2:
        return 1, f'{n_blocks} block(s)'
    if threading.active_count() > 1:
        # fork() copies only the calling thread; a lock held by another one would stay held.
        return 1, 'other threads run in this process'
    if not sys.platform.startswith('linux') or 'fork' not in multiprocessing.get_all_start_methods():
        # On macOS a forked process can crash in the system's libraries, which is why Python
        # starts its workers there without fork(); the blocks then go one after another.
        return 1, f'workers are forked on Linux only, not on {sys.platform}'
    return max(1, min(wanted, n_blocks)), why


# ------------------------------------------------------------------ genome context deltas
def _same(a: Any, b: Any) -> bool:
    return a is b or a == b


def _empty_like(value: Any) -> Any:
    if isinstance(value, list):
        return []
    if isinstance(value, dict):
        return {}
    return 0


def _list_tail(key: str, before: list, after: Any) -> list:
    if not isinstance(after, list) or after[:len(before)] != before:
        raise UnmergedBlockChange(f"genome_ctx[{key!r}]: a block rewrote entries of the list "
                                  "instead of adding to its end")
    return after[len(before):]


def _dict_changes(key: str, before: dict, after: Any) -> List[Tuple[Any, Any]]:
    if not isinstance(after, dict):
        raise UnmergedBlockChange(f"genome_ctx[{key!r}]: a block replaced the table")
    gone = [k for k in before if k not in after]
    if gone:
        raise UnmergedBlockChange(f"genome_ctx[{key!r}]: a block removed the entries {gone[:3]}")
    return [(k, v) for k, v in after.items() if k not in before or not _same(before[k], v)]


def ctx_delta(s0: Dict[str, Any], ctx: Dict[str, Any]) -> Dict[str, Tuple[bool, Any]]:
    """What one block changed in the genome context ``ctx``, which held ``s0`` (a deep copy)
    when the block started: {key: (created, payload)} by CTX_RULES. Raises UnmergedBlockChange
    for a change no rule takes back."""
    out: Dict[str, Tuple[bool, Any]] = {}
    for key in list(s0) + [k for k in ctx if k not in s0]:
        before, after = s0.get(key, _MISSING), ctx.get(key, _MISSING)
        if after is _MISSING:
            raise UnmergedBlockChange(f"genome_ctx[{key!r}] was removed by a block")
        if before is not _MISSING and _same(before, after):
            continue
        rule = CTX_RULES.get(key)
        if rule is None:
            raise UnmergedBlockChange(
                f"a block changed genome_ctx[{key!r}], which has no rule for blocks searched side "
                "by side (synvoy_parallel.CTX_RULES); search with --block_workers 1 or add one")
        created = before is _MISSING
        base = _empty_like(after) if created else before
        if rule == APPEND:
            payload: Any = _list_tail(key, base, after)
        elif rule == PER_BLOCK:
            payload = _dict_changes(key, base, after)
        elif rule == COUNT:
            payload = [(k, v - base.get(k, 0)) for k, v in _dict_changes(key, base, after)]
        elif key == 'judge_numbers':
            payload = _dict_changes(key, base, after)
        elif key == 'judge_writers':
            changes = _dict_changes(key, base, after)
            payload = [(k, _list_tail(f'{key}.{k}', base.get(k, []), v)) for k, v in changes]
        elif key == 'judge_repeats':
            payload = int(after) - int(base or 0)
        else:  # WINDOWS: reduced_search, whose 'windows' list blocks extend
            if not isinstance(after, dict):
                raise UnmergedBlockChange(f"genome_ctx[{key!r}]: a block replaced it")
            other = {k: v for k, v in after.items() if k != 'windows'}
            if other != {k: v for k, v in base.items() if k != 'windows'}:
                raise UnmergedBlockChange(f"genome_ctx[{key!r}]: a block changed more than its windows")
            payload = _list_tail(f'{key}.windows', base.get('windows', []), after.get('windows', []))
        out[key] = (created, payload)
    return out


def apply_ctx_delta(ctx: Dict[str, Any], delta: Dict[str, Tuple[bool, Any]],
                    s0: Dict[str, Any]) -> None:
    """Take one block's changes (ctx_delta) back into the genome context ``ctx``, which holds
    the changes of the blocks before it; ``s0`` is the state every block started from."""
    repeats = 0
    for key, (created, payload) in delta.items():
        rule = CTX_RULES[key]
        if rule == APPEND:
            ctx.setdefault(key, []).extend(payload)
        elif rule == PER_BLOCK:
            table = ctx.setdefault(key, {})
            for k, v in payload:
                table[k] = v
        elif rule == COUNT:
            table = ctx.setdefault(key, {})
            for k, n in payload:
                table[k] = table.get(k, 0) + n
        elif key == 'judge_numbers':
            table = ctx.setdefault(key, {})
            held_before = s0.get(key) or {}
            for k, v in payload:
                # The block saw only what was there before the first block: an id an earlier
                # block wrote is a repeat it could not count (record_judge_numbers).
                if k not in held_before and k in table:
                    repeats += 1
                table[k] = v
        elif key == 'judge_writers':
            table = ctx.setdefault(key, {})
            for k, parents in payload:
                table.setdefault(k, []).extend(parents)
        elif key == 'judge_repeats':
            ctx[key] = ctx.get(key, 0) + payload
        else:  # WINDOWS
            ctx.setdefault(key, {}).setdefault('windows', []).extend(payload)
    if repeats:
        ctx['judge_repeats'] = ctx.get('judge_repeats', 0) + repeats


# ------------------------------------------------------------------ module state
def _counter_delta(before: Dict[Any, Any], after: Dict[Any, Any]) -> List[Tuple[Any, Any]]:
    return [(k, v - before.get(k, 0)) for k, v in after.items()
            if isinstance(v, (int, float)) and v != before.get(k, 0)]


class _Job:
    """What the workers of one genome's blocks inherit when they are forked."""

    def __init__(self, blocks, run_block, args, counters, once_sets, flags, work_dir):
        self.blocks = blocks
        self.run_block = run_block
        self.args = args
        self.counters = counters
        self.once_sets = once_sets
        self.flags = flags
        self.work_dir = work_dir
        self.s0_ctx = copy.deepcopy(args.genome_ctx)
        self.s0_counters = [copy.deepcopy(c) for c in counters]


_JOB: Optional[_Job] = None


@contextlib.contextmanager
def _captured(out_path: str, err_path: str):
    """stdout and stderr of this process (Python's and every tool's) into two files."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = (os.dup(1), os.dup(2))
    out_fd = os.open(out_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    err_fd = os.open(err_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        os.dup2(out_fd, 1)
        os.dup2(err_fd, 2)
        yield
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        finally:
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            for fd in (out_fd, err_fd) + saved:
                os.close(fd)


def _portable(exc: BaseException) -> BaseException:
    """``exc`` as it can travel back to the genome's worker (pickled)."""
    try:
        pickle.loads(pickle.dumps(exc))
        return exc
    except Exception:
        return RuntimeError(f"{type(exc).__name__}: {exc}")


def _run_block_task(i: int):
    """One block in a worker: from the state before the first block, its output into files.
    Returns (result, error, ctx delta, counter deltas, new once-keys, flags)."""
    job = _JOB
    job.args.genome_ctx = copy.deepcopy(job.s0_ctx)
    for counter, start in zip(job.counters, job.s0_counters):
        counter.clear()
        counter.update(copy.deepcopy(start))
    once_start = [set(s) for s in job.once_sets]
    result, error = None, None
    with _captured(os.path.join(job.work_dir, f'{i}.out'), os.path.join(job.work_dir, f'{i}.err')):
        try:
            result = job.run_block(i, job.blocks[i])
        except BaseException as exc:      # StepFailed is a BaseException
            error = _portable(exc)
    try:
        delta = ctx_delta(job.s0_ctx, job.args.genome_ctx)
    except UnmergedBlockChange as exc:
        delta, error = {}, error or exc
    counter_deltas = [_counter_delta(start, counter)
                      for counter, start in zip(job.counters, job.s0_counters)]
    once_new = [sorted(s - start, key=repr) for s, start in zip(job.once_sets, once_start)]
    flags = [bool(f[0]) for f in job.flags]
    return result, error, delta, counter_deltas, once_new, flags


def _emit(work_dir: str, i: int) -> None:
    """Block ``i``'s stdout and stderr, as the block would have written them."""
    for name, stream in (('out', sys.stdout), ('err', sys.stderr)):
        path = os.path.join(work_dir, f'{i}.{name}')
        if os.path.exists(path) and os.path.getsize(path):
            with open(path, 'r', encoding='utf-8', errors='replace') as fh:
                stream.write(fh.read())
            stream.flush()


def search_blocks(blocks: Sequence[Any], run_block: Callable[[int, Any], Any], args,
                  workers: int, work_dir: str, counters: Sequence[Dict[Any, Any]] = (),
                  once_sets: Sequence[set] = (), flags: Sequence[list] = ()) -> Iterator[Any]:
    """The result of ``run_block(i, blocks[i])`` for each block, in block order, with
    ``args.genome_ctx`` and the ``counters`` (dicts of counts), ``once_sets`` (keys of
    warnings given once) and ``flags`` (one-element lists of a bool) changed as the blocks
    would have changed them one after another. ``workers`` <= 1 runs them one after another
    in this process, as before; more forks that many workers (module docstring). ``work_dir``
    (created, then removed) keeps the workers' output files."""
    if workers <= 1 or len(blocks) < 2:
        for i, block in enumerate(blocks):
            yield run_block(i, block)
        return
    global _JOB
    os.makedirs(work_dir, exist_ok=True)
    _JOB = _Job(blocks, run_block, args, list(counters), list(once_sets), list(flags), work_dir)
    s0_ctx = _JOB.s0_ctx
    sys.stdout.flush()
    sys.stderr.flush()
    gc.freeze()            # the workers share these objects' pages; the collector leaves them
    pool = concurrent.futures.ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context('fork'))
    done = False
    try:
        futures = [pool.submit(_run_block_task, i) for i in range(len(blocks))]
        for i, future in enumerate(futures):
            try:
                result, error, delta, counter_deltas, once_new, flag_values = future.result()
            except concurrent.futures.process.BrokenProcessPool as exc:
                raise MemoryKilled(
                    f"a block worker died abruptly ({exc}): killed by a signal, most often because "
                    "it ran out of memory (Nextflow retries the task with more, and the retry "
                    "searches the blocks one after another: --block_workers 1)") from None
            _emit(work_dir, i)
            apply_ctx_delta(args.genome_ctx, delta, s0_ctx)
            for counter, changes in zip(counters, counter_deltas):
                for k, n in changes:
                    counter[k] = counter.get(k, 0) + n
            for once, keys in zip(once_sets, once_new):
                once.update(keys)
            for flag, value in zip(flags, flag_values):
                flag[0] = bool(flag[0]) or value
            if error is not None:
                raise error
            yield result
        done = True
    finally:
        pool.shutdown(wait=True, cancel_futures=not done)
        gc.unfreeze()
        _JOB = None
        shutil.rmtree(work_dir, ignore_errors=True)
