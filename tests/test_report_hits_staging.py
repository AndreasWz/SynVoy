"""The report reads the search's per-genome tables from hits/ (generate_report.py, hits_scan),
and modules/generate_report.nf stages them by name. Every table the report reads (every suffix
of hits_scan but the raw .m8 hits) must be in that find: the places of the genome-wide path
(.genome_wide.tsv) were not, so a run with --goi_genome_wide_path had neither the report's
genome_wide_path block nor its warning (found 2026-10-19)."""
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import generate_report as gr  # noqa: E402


def _hits_scan_suffixes():
    tree = ast.parse((ROOT / "bin" / "generate_report.py").read_text())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
                and any(isinstance(t, ast.Name) and t.id == "hits_scan" for t in node.targets)):
            buckets = node.value.args[1]
            return [e.value if isinstance(e, ast.Constant) else getattr(gr, e.id)
                    for v in buckets.values for e in v.elts]
    raise AssertionError("no hits_scan = scan_dir_by_suffix(...) in generate_report.py")


def test_every_table_the_report_reads_from_hits_is_staged():
    suffixes = [s for s in _hits_scan_suffixes() if s != ".m8"]
    assert {".nonsyntenic.tsv", ".not_searched.tsv", ".target_policy.tsv",
            ".genome_wide.tsv"} <= set(suffixes)
    module = (ROOT / "modules" / "generate_report.nf").read_text()
    start = module.index("find -L hits")
    find = module[start:module.index("staged_results/hits/", start)]
    missing = [s for s in suffixes if f"-name '*{s}'" not in find]
    assert not missing, f"read by the report but not staged: {missing}"
