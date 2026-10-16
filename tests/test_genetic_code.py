"""Genetic codes (preset plan of 2026-10-15, PR-3): every translation can use an NCBI table.

Table 1 everywhere by default, and then no engine gets a genetic-code option, so a default
run is byte-identical. A table is set per process (sequence_utils.set_genetic_code): the
search sets the target's table at the start of each genome, the home-side scripts the home
genome's; translate(), stop_codons(), the frame guard and the engines' options follow it.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
import sys  # noqa: E402

sys.path.insert(0, str(ROOT / "bin"))

import sequence_utils as su  # noqa: E402
import synvoy_label  # noqa: E402


@pytest.fixture(autouse=True)
def standard_code():
    """The process default leaks between tests otherwise."""
    su.set_genetic_code(1)
    yield
    su.set_genetic_code(1)


# ============================================================ translate
@pytest.mark.parametrize("table,want", [
    (1, "LL***W"),      # standard
    (12, "SL***W"),     # CUG-Ser1 yeasts: CTG -> S
    (26, "AL***W"),     # Pachysolen: CTG -> A
    (6, "LLQQ*W"),      # ciliates: TAA/TAG -> Q, TGA a stop
    (4, "LL**WW"),      # TGA -> W
])
def test_translate_with_ncbi_tables(table, want):
    assert su.translate("CTGTTGTAATAGTGATGG", table) == want


def test_table_1_is_the_old_standard_code():
    assert su.translate("ATGNNNTGAuuu") == "MX*F"
    assert su.stop_codons(1) == {"TAA", "TAG", "TGA"}


def test_the_process_default_is_what_unnamed_calls_use():
    assert su.translate("CTG") == "L"
    su.set_genetic_code(12)
    assert su.get_genetic_code() == 12 and su.translate("CTG") == "S"
    assert su.stop_codons() == {"TAA", "TAG", "TGA"}
    su.set_genetic_code(6)
    assert su.stop_codons() == {"TGA"}


@pytest.mark.parametrize("bad", [7, 8, 17, 20, 32, 33, 0, "x", None])
def test_tables_not_every_engine_supports_are_refused(bad):
    with pytest.raises(ValueError):
        su.check_genetic_code(bad)


def test_frame_guard_follows_the_genetic_code():
    """A model with TAG inside is a broken frame under the standard code and a clean one
    under table 6, where TAG is glutamine."""
    cds = "ATGTAGTTGAATTACTCTAAATGA"       # every frame has a stop under table 1
    assert synvoy_label.frame_guard(cds)[1] == "broken"
    su.set_genetic_code(6)
    protein, status = synvoy_label.frame_guard(cds)
    assert (protein, status) == ("MQLNYSK", "ok")


# ============================================================ per target genome
def test_genetic_code_of_a_genome_by_file_name_or_stem():
    spec = su.parse_genetic_code_spec("Candida_dubliniensis.fna=12,Pachysolen=26")
    assert su.genetic_code_for("/data/Candida_dubliniensis.fna", spec) == 12
    assert su.genetic_code_for("Pachysolen.fa.gz", spec) == 26
    assert su.genetic_code_for("Other.fna", spec, 1) == 1
    assert su.genetic_code_for("Other.fna", spec, 4) == 4
    assert su.parse_genetic_code_spec("") == {}
    with pytest.raises(ValueError):
        su.parse_genetic_code_spec("x=33")


def test_engine_options_only_for_a_non_standard_code():
    for tool in ("mmseqs", "tblastn", "miniprot", "synvoy"):
        assert su.genetic_code_flag(tool) == []
    su.set_genetic_code(12)
    assert su.genetic_code_flag("mmseqs") == ["--translation-table", "12"]
    assert su.genetic_code_flag("tblastn") == ["-db_gencode", "12"]
    assert su.genetic_code_flag("miniprot") == ["-T", "12"]
    assert su.genetic_code_flag("synvoy") == ["--genetic_code", "12"]


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
def test_miniprot_gets_the_table_of_the_process():
    import annotate_goi_exons as age
    assert age._miniprot_genetic_code_args() == []
    su.set_genetic_code(12)
    assert age._miniprot_genetic_code_args() == ["-T", "12"]


# ============================================================ wiring
def test_every_engine_call_of_the_search_carries_the_code_option():
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    assert src.count('+ genetic_code_flag("mmseqs")') == 2           # window and genome-wide
    assert src.count('+ genetic_code_flag("tblastn")') == 2           # + the finding run (B1)
    assert src.count('+ genetic_code_flag("synvoy")') == 1           # Smith-Waterman
    assert "set_genetic_code(genetic_code_for(genome_path" in src
    rescue = (ROOT / "bin" / "rescue_strong_synteny.py").read_text()
    assert 'genetic_code_flag("miniprot")' in rescue
    age = (ROOT / "bin" / "annotate_goi_exons.py").read_text()
    assert "cmd.extend(_miniprot_genetic_code_args())" in age
    assert not re.search(r"\('TAA', 'TAG', 'TGA'\)", age)            # stops follow the table


@pytest.mark.parametrize("module,helper", [
    ("normalize_query.nf", "homeArg"), ("annotate_goi.nf", "homeArg"),
    ("prepare_home.nf", "homeArg"), ("extract_flanking.nf", "homeArg"),
    ("iterative_search.nf", "targetArgs"), ("rescue_goi_hull.nf", "targetArgs"),
    ("rescue_strong_synteny.nf", "targetArgs"), ("borrow_annotations.nf", "targetArgs"),
])
def test_modules_hand_the_code_over(module, helper):
    assert f"GeneticCodes.{helper}(params)" in (ROOT / "modules" / module).read_text()


def test_locate_gene_passes_the_home_code_to_translated_searches_only():
    text = (ROOT / "modules" / "locate_gene.nf").read_text()
    assert "--translation-table ${gc}" in text and "-db_gencode ${gc}" in text
    assert 'MM_GC=""' in text                     # nucleotide query: no table
    assert "blastn -query $gene -db $home_genome -out hits_blast.txt -outfmt 6\n" in text


def test_groovy_and_python_accept_the_same_tables():
    text = (ROOT / "lib" / "GeneticCodes.groovy").read_text()
    m = re.search(r"SUPPORTED = \(1\.\.6\) \+ \(9\.\.16\) \+ \(21\.\.31\)", text)
    assert m, "lib/GeneticCodes.groovy: SUPPORTED changed, update this test and sequence_utils"
    assert su.SUPPORTED_GENETIC_CODES == set(range(1, 7)) | set(range(9, 17)) | set(range(21, 32))


NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_an_unsupported_table_stops_the_run_at_the_start(tmp_path):
    (tmp_path / "q.faa").write_text(">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKV\n")
    (tmp_path / "h.fna").write_text(">c1\nACGTACGTACGTACGTACGTACGTACGT\n")
    (tmp_path / "codes.tsv").write_text("genome\ttable\nA.fna\t12\nB.fna\t33\n")
    (tmp_path / "noconda.config").write_text("conda.enabled = false\nprocess.conda = null\n")
    proc = subprocess.run(
        [NEXTFLOW, "run", str(ROOT / "main.nf"), "-preview", "-ansi-log", "false",
         "-c", str(tmp_path / "noconda.config"), "--mode", "pro", "--query", "q.faa",
         "--home_genome", "h.fna", "--home_species", "Testus_testus", "--outdir", "out",
         "--target_genetic_codes", str(tmp_path / "codes.tsv"), "--home_genetic_code", "12"],
        cwd=tmp_path, capture_output=True, text=True, timeout=300)
    text = re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout + proc.stderr)
    assert proc.returncode != 0
    assert "genome B.fna has table 33" in text
    assert "home_genetic_code" not in text.split("Input validation failed")[-1]
