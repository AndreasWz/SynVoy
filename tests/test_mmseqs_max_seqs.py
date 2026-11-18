"""MMseqs2's --max-seqs is a parameter (mmseqs_max_seqs, default 300 = MMseqs2's own default).

The genome-wide search, the window searches and LOCATE_GENE pass it; the counters of the queries
that reached it (mmseqs_maxseqs_genome, mmseqs_maxseqs_window) use the value; the report warns
when one is not 0 and says that the count is a lower bound; effective_params.json records it in
run_resources. The report's warning for the model cap of a GOI locus (goi_model_cap_stopped) is
tested here too, with the same staged report.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import iterative_search_runner as isr  # noqa: E402


def _proc(rc, stderr=""):
    return types.SimpleNamespace(returncode=rc, stdout="", stderr=stderr)


def _args(tmp_path, **kw):
    a = dict(output_dir=str(tmp_path), mmseqs_sens=9.5, mmseqs_split_memory_limit="8G", evalue=0.01,
             mmseqs_verbosity=1, aug_relaxed_evalue_cap=10.0, aug_relaxed_evalue_mult=1000.0,
             window_mmseqs_min_length=30, genome_ctx={})
    a.update(kw)
    return types.SimpleNamespace(**a)


def _recording(record):
    def run(cmd, *a, **k):
        if cmd[1] == "search":
            record.append(list(cmd))
        if cmd[1] == "convertalis":
            open(cmd[5], "w").close()
        return _proc(0)
    return run


def _max_seqs(cmd):
    return cmd[cmd.index("--max-seqs") + 1]


def test_the_default_is_mmseqs2s_own():
    assert isr.MMSEQS_DEFAULT_MAX_SEQS == 300
    config = open(os.path.join(ROOT, "nextflow.config")).read()
    assert re.search(r"^\s*mmseqs_max_seqs\s*=\s*300\b", config, re.M)


def test_the_genome_and_the_window_searches_pass_it(tmp_path, monkeypatch):
    for args, expect in ((_args(tmp_path), "300"), (_args(tmp_path, mmseqs_max_seqs=1000), "1000")):
        cmds = []
        monkeypatch.setattr(isr.subprocess, "run", _recording(cmds))
        assert isr.run_mmseqs_easy_search_with_retries("q.faa", "g.fna", str(tmp_path / "h.m8"),
                                                       str(tmp_path / "t"), args, 2, "G", used={})[0]
        assert isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"),
                                                     str(tmp_path / "w"), args, 1, "G_b1")[0]
        assert [_max_seqs(c) for c in cmds] == [expect, expect]


def test_the_counter_uses_the_value(tmp_path):
    p = tmp_path / "h.m8"
    p.write_text("".join(f"GOI_x\tc1\t{i}\n" for i in range(400)) + "gene-a\tc1\t1\n" * 50)
    assert isr.queries_at_max_seqs(str(p), 300) == 1
    assert isr.queries_at_max_seqs(str(p), 500) == 0
    assert isr.queries_at_max_seqs(str(p), 50) == 2


def test_the_runner_accepts_the_option():
    p = subprocess.run([sys.executable, os.path.join(BIN, "iterative_search_runner.py"), "--help"],
                       capture_output=True, text=True)
    assert p.returncode == 0 and "--mmseqs_max_seqs" in p.stdout


def test_the_modules_pass_it():
    search = open(os.path.join(ROOT, "modules", "iterative_search.nf")).read()
    assert "--mmseqs_max_seqs ${params.mmseqs_max_seqs}" in search
    locate = open(os.path.join(ROOT, "modules", "locate_gene.nf")).read()
    assert "--max-seqs ${params.mmseqs_max_seqs}" in locate
    resolver = open(os.path.join(ROOT, "modules", "resolve_effective_params.nf")).read()
    assert "mmseqs_max_seqs: params.mmseqs_max_seqs" in resolver
    main = open(os.path.join(ROOT, "main.nf")).read()
    assert "Invalid --mmseqs_max_seqs" in main


def test_locate_gene_prints_mmseqs2s_errors():
    locate = open(os.path.join(ROOT, "modules", "locate_gene.nf")).read()
    assert "Math.max(1, params.mmseqs_verbosity.toString().toInteger())" in locate
    assert "-v ${mmseqs_v}" in locate


# --------------------------------------------------------------------------- the report

def _report(tmp_path, extra_cols=None):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "hits").mkdir()
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    cols = {"genome": "Gx", "median_flank_identity": "80.0", "n_flanking_genes": "20", "blocks_searched": "3",
            "significant_blocks": "1", "best_block_esyn": "1e-9", "best_block_context": "complete",
            "truncated_blocks": "0", "genetic_code": "1"}
    cols.update({f"limit_{k}": "0" for k in isr.LIMIT_KEYS})
    cols.update(extra_cols or {})
    (stage / "hits" / "locus_1_Gx.fna.target_policy.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(cols.values()) + "\n")
    out = tmp_path / "r.json"
    p = subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                        "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return json.loads(out.read_text())


def test_the_report_warns_when_a_query_reached_the_limit(tmp_path):
    rep = _report(tmp_path, {"limit_mmseqs_maxseqs_genome": "2"})
    warn = [w for w in rep["warnings"] if "--mmseqs_max_seqs" in w]
    assert len(warn) == 1
    assert "mmseqs_maxseqs_genome 2 (Gx)" in warn[0] and "lower bound" in warn[0]
    assert "mmseqs_maxseqs_window" not in warn[0]


def test_no_count_no_warning(tmp_path):
    rep = _report(tmp_path)
    assert not any("--mmseqs_max_seqs" in w for w in rep["warnings"])
    assert not any("goi_max_models_per_locus" in w for w in rep["warnings"])


def test_the_report_warns_when_the_model_cap_stopped_the_rounds(tmp_path):
    rep = _report(tmp_path, {"limit_goi_model_cap_stopped": "3"})
    warn = [w for w in rep["warnings"] if "--goi_max_models_per_locus" in w]
    assert len(warn) == 1
    assert "3 time(s) (Gx)" in warn[0] and "limits.goi_model_cap_stopped" in warn[0]


# --------------------------------------------------------------------------- the start's check

def _preview(folder, *extra):
    sys.path.insert(0, os.path.join(ROOT, "tests"))
    from test_removed_params_preview import _preview as preview
    folder.mkdir()
    return preview(folder, *extra)


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")
def test_the_start_checks_the_value(tmp_path):
    code, text = _preview(tmp_path / "bad", "--mmseqs_max_seqs", "0")
    assert code != 0 and "Invalid --mmseqs_max_seqs (0). Must be a whole number >= 1." in text
    code, text = _preview(tmp_path / "good", "--mmseqs_max_seqs", "1000", "--mmseqs_verbosity", "0")
    assert code == 0, text[-2000:]
