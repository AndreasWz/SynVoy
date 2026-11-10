"""The tree's alignment does not depend on the threads: MAFFT's iterative refinement runs in one thread
(--threadit 0), while the threads (compute_tree_cpus) run the pairwise or distance stage and the progressive
step. Refined by several threads, MAFFT 7.526's L-INS-i alignment depended on the threads' timing: six GOI
proteins of one Candida gene gave two alignments in 16 runs at --thread 2, and so two trees from one search."""
import os
import random
import shutil
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))

import compute_tree as ct  # noqa: E402


def test_mafft_refines_in_one_thread(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(ct.subprocess, "run", lambda cmd, *a, **k: seen.append(cmd))
    ct.run_mafft("goi_only.faa", str(tmp_path / "out.aln"), threads=4)
    cmd = seen[0]
    assert cmd[cmd.index("--thread") + 1] == "4"
    assert cmd[cmd.index("--threadit") + 1] == "0"
    assert cmd[-1] == "goi_only.faa"


def _family(n=8, length=120, seed=7):
    """n diverged copies of one random protein (substitutions and short indels): made up, no real data."""
    rnd = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    root = [rnd.choice(aa) for _ in range(length)]
    out = []
    for i in range(n):
        s = list(root)
        for _ in range(length // 3):
            s[rnd.randrange(len(s))] = rnd.choice(aa)
        for _ in range(3):
            p = rnd.randrange(len(s) - 5)
            if rnd.random() < 0.5:
                del s[p:p + rnd.randint(1, 4)]
            else:
                s[p:p] = [rnd.choice(aa) for _ in range(rnd.randint(1, 4))]
        out.append((f"GOI_x|s{i}", "".join(s)))
    return out


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft is not on PATH")
def test_the_alignment_is_the_same_at_any_thread_count(tmp_path):
    """One alignment at 1, 2 and 4 threads. On these 25 made-up proteins the previous command gave 2 and 3
    alignments in two rounds of these 8 runs; the command test fails on it every time."""
    faa = tmp_path / "goi_only.faa"
    faa.write_text("".join(f">{name}\n{seq}\n" for name, seq in _family(n=25)))
    alignments = set()
    for i, threads in enumerate((1, 2, 4, 2, 4, 2, 4, 2)):
        aln = tmp_path / f"run{i}.aln"
        ct.run_mafft(str(faa), str(aln), threads=threads)
        alignments.add(aln.read_text())
    assert len(alignments) == 1
    assert next(iter(alignments)).count(">") == 25
