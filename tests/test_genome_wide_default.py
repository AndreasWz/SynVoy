"""The strong hit without synteny (--goi_genome_wide_path) is on by default: in the config, in
the runner's own default, and in the docs that state the default."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_the_config_default_is_true():
    text = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*goi_genome_wide_path\s*=\s*true\b", text, re.M)
    assert not re.search(r"^\s*goi_genome_wide_path\s*=\s*false\b", text, re.M)


def test_the_runner_default_is_true():
    text = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    assert 'parser.add_argument("--goi_genome_wide_path", type=str2bool, default=True,' in text


def test_the_module_passes_the_parameter():
    nf = (ROOT / "modules" / "iterative_search.nf").read_text()
    assert "--goi_genome_wide_path ${params.goi_genome_wide_path}" in nf


def test_the_docs_state_the_default():
    params = (ROOT / "docs" / "PARAMETERS.md").read_text()
    assert "### `goi_genome_wide_path`\n**Default:** `true`" in params
    assert "| `--goi_genome_wide_path false` |" in params          # the comparison switch
    algorithm = (ROOT / "docs" / "ALGORITHM.md").read_text()
    assert "`goi_genome_wide_path` (off by default)" not in algorithm
    assert "`goi_genome_wide_path` (on by default; `false` switches it off)" in algorithm
    output = (ROOT / "docs" / "OUTPUT.md").read_text()
    assert "only with `--goi_genome_wide_path`" not in output.lower()
    assert "0 and empty without `--goi_genome_wide_path`" not in output
    runner = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    config = (ROOT / "nextflow.config").read_text()
    assert "with --goi_genome_wide_path also in the whole genome" not in runner   # --decoy_mode help
    assert "with goi_genome_wide_path also in the whole genome" not in config
    for doc in ("README.md", "docs/USAGE.md"):
        assert "genome_wide_strong_hit" in (ROOT / doc).read_text(), doc   # the answer is named
