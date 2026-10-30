"""Fix-phase item 14: the report never follows a directory listing's order.

staging_diagnostics.dirs.*.sample_matches came from glob.glob, which returns files in the file
system's listing order, so two runs on LRZ scratch never gave byte-identical reports.
"""
import os
import random
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402


def test_sample_matches_are_sorted_whatever_the_listing_order(tmp_path, monkeypatch):
    names = [f"locus_{i}__G{j}.fna.gff" for i in (1, 2) for j in range(8)]
    for n in names:
        (tmp_path / n).write_text("##gff-version 3\n")
    real_glob = gr.glob.glob
    seen = []
    for seed in range(5):
        def shuffled(pattern, _seed=seed):
            out = real_glob(pattern)
            random.Random(_seed).shuffle(out)
            return out
        monkeypatch.setattr(gr.glob, "glob", shuffled)
        seen.append(gr._dir_diagnostics(str(tmp_path), ["*.gff"]))
    assert all(d == seen[0] for d in seen)
    assert seen[0]["sample_matches"] == sorted(names)[:5]
    assert seen[0]["pattern_match_count"] == len(names)


def test_every_listing_of_the_report_is_sorted():
    """No glob or listdir in generate_report.py whose order reaches the report unsorted: every
    call sits inside sorted(...); every scandir loop sorts its entries."""
    import ast
    src = open(os.path.join(ROOT, "bin", "generate_report.py")).read()
    tree = ast.parse(src)
    parents = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node

    def listing(call):
        f = call.func
        return (isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name)
                and (f.value.id, f.attr) in (("glob", "glob"), ("os", "listdir")))

    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and listing(n)]
    assert calls
    for call in calls:
        up = parents.get(call)
        assert isinstance(up, ast.Call) and getattr(up.func, "id", "") == "sorted", ast.unparse(call)
    assert src.count("with os.scandir(") == src.count("for entry in sorted(it, key=lambda e: e.name)") + 1
    assert "result[bucket].sort()" in src            # scan_dir_by_suffix's one scandir
