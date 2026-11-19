"""The read-through rule of the genome-wide path's retrocopy test, and where it is described.

A genome-wide gene model with one exon across an exon junction of the home gene is a retrocopy
candidate (RetrocopyCandidate=1, no answer), unless the exon's CDS is READTHROUGH_INTRON_MIN_NT
(30) nt or more per spanned junction longer than 3 x its aligned query span: miniprot reads a short
intron that is a multiple of 3 and free of stops as coding. The excess is counted, not located, so
an exon across two junctions is excused by 60 nt wherever they lie, also when they are one
read-through intron next to a junction whose intron was lost. Every text that describes the test
says so, and every text that says when it runs names the miniprot home model with canonical
splice sites and the candidate it did not run on (retrocopy_untested). All synthetic."""
import ast
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

THREE_EXONS = isr.home_exon_junctions([120, 100, 80])            # junctions after 120 and 220
LIMIT = 2 * isr.READTHROUGH_INTRON_MIN_NT                        # two junctions: 60 nt


def _flag(cds_nt, walk_nt=0):
    """RetrocopyCandidate of a one-exon model over the whole 300-aa query (both junctions)."""
    exon = {"qstart": 1, "qend": 300, "gstart": 1_000, "gend": 1_000 + cds_nt}
    if walk_nt:
        exon["terminal_walk_nt"] = walk_nt
    cand = {"start": 1001, "end": 1_000 + cds_nt, "model_parent": "GOI_q", "model_exons": [exon],
            "gff": ["c1\tx\tmRNA\t1001\t2000\t.\t+\t.\tID=m;SynVoy_Parent=GOI_q"]}
    isr.mark_genome_wide_candidates([cand], THREE_EXONS, {})
    return isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8]).get("RetrocopyCandidate")


def test_one_read_through_intron_of_60_nt_next_to_a_lost_intron_is_excused():
    """The target lost intron 1 and kept intron 2 (in frame, stop-free), which miniprot read as
    coding: one exon across both junctions. 45 nt of excess hold one intron, not two, and stay
    flagged; from 60 nt on the excess is read as an intron at each junction, so the model is no
    retrocopy candidate although one junction has no intron."""
    assert LIMIT == 60
    assert _flag(900 + 45) == "1"
    assert _flag(900 + LIMIT - 1) == "1"
    assert _flag(900 + LIMIT) == "0"
    assert _flag(900 + 75) == "0"
    assert _flag(900 + LIMIT + 9, walk_nt=9) == "0"             # the terminal walk is no intron
    assert _flag(900 + LIMIT, walk_nt=9) == "1"
    assert _flag(903) == "1"                                     # an intron-less copy, its stop


def _norm(text):
    return " ".join(re.sub(r"(?m)^\s*(?://|#)\s?", " ", text).split())


def _between(text, start, end):
    return text[text.index(start):text.index(end, text.index(start) + len(start))]


def _help(flag):
    tree = ast.parse(open(os.path.join(ROOT, "bin", "iterative_search_runner.py"), encoding="utf-8").read())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument"
                and node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == flag):
            return next(ast.literal_eval(k.value) for k in node.keywords if k.arg == "help")
    raise AssertionError(flag)


def _texts():
    def read(*parts):
        return open(os.path.join(ROOT, *parts), encoding="utf-8").read()
    runner = read("bin", "iterative_search_runner.py")
    return {
        "ALGORITHM.md 7.11": _between(read("docs", "ALGORITHM.md"), "### 7.11", "### 7.12"),
        "OUTPUT.md RetrocopyCandidate": next(line for line in read("docs", "OUTPUT.md").splitlines()
                                             if line.startswith("| `RetrocopyCandidate` |")),
        "PARAMETERS.md goi_genome_wide_path": _between(read("docs", "PARAMETERS.md"),
                                                       "### `goi_genome_wide_path`", "\n### "),
        "runner comment": _between(runner, "RETROCOPY_JUNCTION_FLANK_AA = 10",
                                   "READTHROUGH_INTRON_MIN_NT = 30"),
        "nextflow.config comment": _between(read("nextflow.config"), "goi_join_leftover_hits = ",
                                            "goi_genome_wide_path = "),
        "--goi_genome_wide_path help": _help("--goi_genome_wide_path"),
        "GENOME_WIDE_NOTE": gr.GENOME_WIDE_NOTE,
    }


def test_every_description_of_the_retrocopy_test_names_the_read_through_rule_and_its_limit():
    limit = re.compile(r"read-through intron of (?:>= ?)?%d nt\b[^.;]*\blost\b" % LIMIT)
    for where, text in _texts().items():
        text = _norm(text)
        assert "retrocop" in text, where
        assert re.search(r"read as coding|read-through", text), where
        assert limit.search(text), (where, text)


def test_every_statement_of_when_the_retrocopy_test_runs_names_the_untested_rule():
    """F1 (laptop gate of 18 Nov 2026): the test ran only with an annotated home gene, and a
    candidate it did not run on could be an answer. The texts that say when it runs name the
    miniprot home model with canonical splice sites, retrocopy_untested and the one-exon case; the
    removed parameter's row no longer says that the test is always on."""
    def read(*parts):
        return open(os.path.join(ROOT, *parts), encoding="utf-8").read()
    texts = _texts()
    runner = read("bin", "iterative_search_runner.py")
    statements = {
        "ALGORITHM.md 7.11": texts["ALGORITHM.md 7.11"],
        "OUTPUT.md RetrocopyCandidate": texts["OUTPUT.md RetrocopyCandidate"],
        "PARAMETERS.md goi_genome_wide_path": texts["PARAMETERS.md goi_genome_wide_path"],
        "runner comment": _between(runner, "# FINAL step 5 (--goi_genome_wide_path)", "GENOME_WIDE_COLUMNS = "),
        "CHANGELOG": next(line for line in read("CHANGELOG.md").splitlines()
                          if line.startswith("- **A strong hit without synteny")),
    }
    for where, text in statements.items():
        text = _norm(text)
        assert "retrocopy_untested" in text, where
        assert "canonical splice sites" in text, where
        assert "one exon" in text, where
    removed = next(line for line in read("docs", "PARAMETERS.md").splitlines()
                   if line.startswith("| `goi_genome_wide_without_blocks` |"))
    assert "always on" not in removed
