"""B1 short_query, engine settings (preset plan of 2026-10-15, PR-7): finding only.

* window tblastn: word size 3 -> 2, composition statistics 2 -> 0; a run with other values
  than tblastn's defaults is used for finding, and the run with the defaults still gives
  EeffTblastn (the decision between the label's two E sources must not move);
* window MMseqs2: --min-length (ORF codons) 30 -> 10, the genome-wide search unchanged;
* miniprot: queries below 100 aa (60 today) get -L 10 in the windows;
* LOCATE_GENE: a query < 100 aa without any home hit is searched once more with tblastn
  -seg no, when --modules selects short_query (never under legacy or none).
The judge's blastp is untouched; legacy applies the old preset only.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as age  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import profile_hits  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402

B1_KEYS = {"window_tblastn_word_size": 2, "window_tblastn_comp_stats": 0,
           "window_mmseqs_min_length": 10, "miniprot_short_orf_max_len": 100}
DEFAULTS = {"classify_complete_min_qcov": 0.7, "classify_fragment_max_qcov": 0.4,
            "window_tblastn_word_size": 3, "window_tblastn_comp_stats": 2,
            "window_mmseqs_min_length": 30, "miniprot_short_orf_max_len": 60}


def test_b1_sets_the_engine_keys_and_legacy_does_not():
    f = {"query_len": 50}
    auto = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="auto", features=f)
    for k, v in B1_KEYS.items():
        assert auto["settings"][k] == v and auto["source"][k] == "module:B1", k
    legacy = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"},
                         modules_spec="legacy", features=f)
    for k in B1_KEYS:
        assert legacy["settings"][k] == DEFAULTS[k], k


def test_the_most_permissive_value_wins_when_two_modules_set_one():
    a = sm.Module("Z1", "z1", "", "", lambda f: (True, ""), lambda f: {"window_tblastn_word_size": 2,
                                                                    "miniprot_short_orf_max_len": 80})
    b = sm.Module("Z2", "z2", "", "", lambda f: (True, ""), lambda f: {"window_tblastn_word_size": 3,
                                                                    "miniprot_short_orf_max_len": 100})
    settings, _ = sm.merge([a, b], {})
    assert settings == {"window_tblastn_word_size": 2, "miniprot_short_orf_max_len": 100}


def test_tblastn_finding_options_only_away_from_the_defaults():
    class A:
        window_tblastn_word_size = 3
        window_tblastn_comp_stats = 2
    assert isr.window_tblastn_finding_options(A) == []
    A.window_tblastn_word_size, A.window_tblastn_comp_stats = 2, 0
    assert isr.window_tblastn_finding_options(A) == ["-word_size", "2", "-comp_based_stats", "0"]


def test_the_runner_wires_the_engine_keys():
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    # EeffTblastn keeps the run with tblastn's defaults; the other run finds
    i = src.index("finding_opts = window_tblastn_finding_options(args)")
    block = src[i - 1500:src.index("# ========== 3. Smith-Waterman Search", i)]
    assert "parse_hits(\n                tblastn_finding_file" in block
    assert "parse_hits(tblastn_hits_file, 0.0, 1, float('inf')" in block
    assert 'search_cmd += ["--min-length", str(int(args.window_mmseqs_min_length))]' in src
    assert "_age.MINIPROT_SHORT_ORF_MAX_LEN = int(getattr(args, 'miniprot_short_orf_max_len', 60))" in src
    # the genome-wide MMseqs2 search has no ORF option
    start = src.index("def run_mmseqs_easy_search_with_retries")
    gw = src[start:src.index("\ndef ", start + 10)]
    assert '"mmseqs", "easy-search"' in gw or "easy-search" in gw
    assert "--min-length" not in gw


def test_miniprot_short_peptide_mode_follows_the_module_setting(monkeypatch):
    seen = []

    class Done:
        returncode, stdout, stderr = 0, "", ""

    monkeypatch.setattr(age, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(age.subprocess, "run", lambda cmd, **kw: seen.append(cmd) or Done())
    query = "M" + "ACDEFGHIKLMNPQRSTVWY" * 4          # 81 aa
    age._run_miniprot(query, "ACGT" * 50, "c1", sensitive=True)
    assert "-L" not in seen[-1]
    monkeypatch.setattr(age, "MINIPROT_SHORT_ORF_MAX_LEN", 100)
    age._run_miniprot(query, "ACGT" * 50, "c1", sensitive=True)
    assert seen[-1][seen[-1].index("-L") + 1] == "10" and "-S" not in seen[-1]


# ============================================================ LOCATE_GENE retry
def test_the_hit_profile_records_the_retry(tmp_path, monkeypatch):
    (tmp_path / "q.faa").write_text(">q\nMKTAYIAKQRQISFVKSHFSRQ\n")
    (tmp_path / "b.txt").write_text("q\tc1\t90.0\t20\t2\t0\t1\t20\t100\t159\t1e-5\t40.0\n")
    (tmp_path / "m.m8").write_text("")
    for retry, want in (("tblastn_seg_no", "tblastn_seg_no"), ("", None)):
        out = tmp_path / f"hp{retry}.json"
        monkeypatch.setattr(sys, "argv", ["profile_hits.py", "--blast", str(tmp_path / "b.txt"),
                                          "--mmseqs", str(tmp_path / "m.m8"), "--query", str(tmp_path / "q.faa"),
                                          "--locus_gap", "40000", "--output", str(out)]
                            + (["--locate_retry", retry] if retry else []))
        profile_hits.main()
        hp = json.loads(out.read_text())
        assert hp.get("locate_retry") == want
        assert sm.flat_features(hp, None, None)["locate_retry"] == want


def _retry_block():
    text = (ROOT / "modules" / "locate_gene.nf").read_text()
    block = text[text.index("    LOCATE_RETRY=\"\""):text.index("    # 3. Combine and convert to BED")]
    # Groovy -> shell: \$ is $, \\ is \
    return block.replace("\\$", "$").replace("\\\\", "\\")


@pytest.mark.parametrize("selected,query,blast_rows,retried", [
    ("true", "MKTAYIAKQRQISF", "", True),                              # no hit at all
    ("true", "MKTAYIAKQRQISF", "q\tc1\t50\t10\t5\t0\t1\t10\t1\t30\t5.0\t20\n", True),   # E 5 only
    ("true", "MKTAYIAKQRQISF", "q\tc1\t90\t14\t1\t0\t1\t14\t1\t42\t1e-4\t30\n", False),  # a hit
    ("false", "MKTAYIAKQRQISF", "", False),                            # module not selected
    ("true", "M" * 120, "", False),                                    # not a short query
])
def test_locate_retry_shell_block(tmp_path, selected, query, blast_rows, retried):
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "tblastn").write_text("#!/bin/bash\necho \"$@\" > \"$(dirname \"$0\")/../tblastn_args\"\n"
                                  "printf 'q\\tc1\\t80\\t14\\t2\\t0\\t1\\t14\\t1\\t42\\t1e-3\\t28\\n' > hits_blast.txt\n")
    (stub / "tblastn").chmod(0o755)
    (tmp_path / "q.faa").write_text(f">q\n{query}\n")
    (tmp_path / "hits_blast.txt").write_text(blast_rows)
    (tmp_path / "hits_mmseqs.m8").write_text("")
    script = ("BLAST_CMD=tblastn\n" + _retry_block()
              .replace("${seg_no_retry}", selected).replace("${params.search_evalue}", "0.01")
              .replace("$gene", "q.faa").replace("$home_genome", "h.fna").replace("${tbn_gc}", "")
              + 'echo "RETRY=[$LOCATE_RETRY]"\n')
    proc = subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True,
                          env=dict(os.environ, PATH=f"{stub}:{os.environ['PATH']}"))
    assert proc.returncode == 0, proc.stderr
    if retried:
        assert "RETRY=[--locate_retry tblastn_seg_no]" in proc.stdout
        assert "-seg no" in (tmp_path / "tblastn_args").read_text()
    else:
        assert "RETRY=[]" in proc.stdout and not (tmp_path / "tblastn_args").exists()


NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_modules_select_in_groovy(tmp_path):
    text = (ROOT / "main.nf").read_text()

    def func(name):
        start = text.index(f"def {name}(")
        return text[start:text.index("\n}\n", start) + 3]

    (tmp_path / "t.nf").write_text(
        func("presetModules") + "\ndef modulesSpec() { return params.modules.toString().trim() }\n"
        + func("modulesSelect")
        + "\nworkflow { println \"SELECT=\" + modulesSelect('short_query') + \",\" + modulesSelect('B1') }\n")
    cases = {"legacy": False, "auto": True, "auto,-short_query": False, "auto,-B1": False,
             "none": False, "none,+short_query": True, "short_query,A1": True,
             "local_family": False, "auto,-local_family": True}
    for spec, want in cases.items():
        proc = subprocess.run([NEXTFLOW, "run", str(tmp_path / "t.nf"), "-ansi-log", "false",
                               "--modules", spec], cwd=tmp_path, capture_output=True, text=True, timeout=300)
        m = re.search(r"SELECT=(\w+),(\w+)", proc.stdout)
        assert m, proc.stdout + proc.stderr
        assert m.group(1) == m.group(2) == str(want).lower(), spec
