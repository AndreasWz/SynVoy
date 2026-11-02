"""Fix-phase item 12: SynVoy's fragment ids have one reader (sequence_utils.fragment_of).

generate_fragments writes <id>|frag_<s>_<e>_<total> and, for a protein above 100 aa, 50-aa sliding
windows <id>|slide_<s>_<e>_<total>. The search's hit mapping gave |frag_ and |exon_ an offset but
not |slide_: a hit's qstart and qend stayed window-relative and fed union coverage and the chains.
And the full-length test read a |slide_ id as a full-length GOI query.
"""
import ast
import os
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import build_home_paralog_panel as panel  # noqa: E402
import fragment_query as fq  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
from sequence_utils import fragment_of  # noqa: E402


@pytest.mark.parametrize("qid,want", [
    ("GOI_x|frag_1_100_200", ("frag", 1, 100, 200)),
    ("GOI_x|frag_101_200_200", ("frag", 101, 200, 200)),
    ("GOI_x|slide_26_75_200", ("slide", 26, 75, 200)),
    ("GOI_x~rel|slide_151_200_200", ("slide", 151, 200, 200)),
    ("GOI_x~rev|slide_1_50_200", ("slide", 1, 50, 200)),
    ("GOI_lcl|ORF1|slide_1_50_200", ("slide", 1, 50, 200)),
    # An id that merely holds the tag is a fragment too, as |frag_ always read; it has no positions.
    ("GOI_X|frag_2", ("frag", None, None, None)),
    ("GOI_lcl|ORF1|frag_1", ("frag", None, None, None)),
    ("GOI_x|slide_7", ("slide", None, None, None)),
    ("GOI_x", None),
    ("GOI_x|exon_2", None),
    ("GOI_x|Genome_fna_b0_l1_exon_ann", None),
    ("GOI_x~rel", None),
    ("GOI_slide_1_50_200", None),
    ("GOI_frag_1_50_200", None),
    ("", None),
    (None, None),
])
def test_fragment_of(qid, want):
    assert fragment_of(qid) == want


def test_every_generated_id_reads_back_its_own_residues():
    seq = "".join("ACDEFGHIKLMNPQRSTVWY"[(i * 7 + i // 3) % 20] for i in range(437))
    frags = fq.generate_fragments(seq, "GOI_x", min_size=15)
    assert {fragment_of(fid)[0] for fid, _s, _d in frags} == {"frag", "slide"}
    for fid, fseq, _desc in frags:
        _kind, start, end, total = fragment_of(fid)
        assert total == len(seq) and seq[start - 1:end] == fseq, fid


def test_a_short_protein_has_no_sliding_windows():
    frags = fq.generate_fragments("A" * 100, "GOI_x", min_size=15)
    assert all(fragment_of(fid)[0] == "frag" for fid, _s, _d in frags)


@pytest.mark.parametrize("qid,kind,start", [
    ("g|frag_1_100_100", "full", 1),
    ("g|frag_1_60_100", "half", 1),
    ("g|frag_34_66_100", "third", 34),
    ("g|frag_76_100_100", "quarter", 76),
    ("g|slide_26_75_100", "sliding", 26),
    ("g|exon_3", "exon", 1),
    ("g|frag_2", "full", 1),
    ("g|slide_2", "sliding", 1),
    ("g", "full", 1),
])
def test_parse_fragment_id_reads_through_the_helper(qid, kind, start):
    parsed = fq.parse_fragment_id(qid)
    assert (parsed["gene_id"], parsed["fragment_type"], parsed["start"]) == ("g", kind, start)


@pytest.mark.parametrize("qid,full", [
    ("GOI_x", True),
    ("GOI_x|Genome_fna_b0_l1_exon_ann", True),
    ("GOI_x|exon_1", False),
    ("GOI_x|frag_1_50_100", False),
    ("GOI_x|slide_1_50_120", False),
    ("gene-a", False),
])
def test_a_sliding_window_is_no_full_length_query(qid, full):
    assert isr.is_full_length_goi_query_id(qid) is full


def test_a_sliding_window_is_never_the_home_query():
    """PREPARE_INITIAL_DB's fallback writes fragments and windows into the query set when the GOI
    has no exon records: the query of a window is the protein, never a window."""
    db = {"GOI_x": {"id": "GOI_x", "seq": "M" * 120},
          "GOI_x|slide_1_50_120": {"id": "GOI_x|slide_1_50_120", "seq": "A" * 50},
          "GOI_x|frag_1_60_120": {"id": "GOI_x|frag_1_60_120", "seq": "C" * 60}}
    queries, _by_parent, _alias = isr.select_goi_queries(db, set(db))
    assert [q["id"] for q in queries] == ["GOI_x"]


def test_the_panels_home_record_is_no_sliding_window(tmp_path):
    """Same rule, build_home_paralog_panel.home_goi_record: a window with no full record present."""
    f = tmp_path / "goi.faa"
    f.write_text(">GOI_a|slide_1_50_120\n" + "A" * 50 + "\n>GOI_a|frag_1_60_120\n" + "C" * 60 + "\n")
    assert panel.home_goi_record(str(f)) is None
    f.write_text(">GOI_a|slide_1_50_120\n" + "A" * 50 + "\n>GOI_a\n" + "M" * 40 + "\n")
    assert panel.home_goi_record(str(f))[0] == "GOI_a"


def test_no_reader_of_the_old_forms_is_left():
    """By the source: every test for a fragment id goes through the helper."""
    for name in ("iterative_search_runner.py", "build_home_paralog_panel.py", "fragment_query.py",
                 "annotate_goi_exons.py", "generate_report.py", "rescue_goi_hull.py"):
        with open(os.path.join(ROOT, "bin", name)) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Compare) and isinstance(node.left, ast.Constant):
                assert node.left.value not in ("|frag_", "|slide_", "frag_", "slide_"), f"{name}:{node.lineno}"


def test_a_tree_has_no_fragment_leaves(tmp_path):
    """COMPUTE_TREE's own filter, cut from the module: the pieces of PREPARE_INITIAL_DB's fallback reach
    goi_for_tree.faa, and they used to become leaves drawn as species ('slide 26 75 340')."""
    with open(os.path.join(ROOT, "modules", "compute_tree.nf")) as fh:
        text = fh.read()
    script = text[text.index('python3 -c "') + len('python3 -c "'):text.index('" > goi_only.faa')]
    script = script.replace("${projectDir}", ROOT)
    (tmp_path / "all_sequences.faa").write_text(
        ">GOI_x\nMKV\n>GOI_x|frag_1_60_120\nMK\n>GOI_x|slide_26_75_120\nKV\n>GOI_x|exon_1\nM\n"
        ">GOI_x|Genome_fna_b0_l1_exon_ann extra words\nMKL\n>gene-a\nAAA\n")
    done = subprocess.run([sys.executable, "-c", script], cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert [ln[1:] for ln in done.stdout.splitlines() if ln.startswith(">")] == [
        "GOI_x", "GOI_x|Genome_fna_b0_l1_exon_ann extra words"]


class _Sentinel(BaseException):
    """Passes the window's own except Exception."""


def _window_args(tmp_path):
    """The attributes a window reads up to its first coverage: everything else is None."""
    class Args(types.SimpleNamespace):
        def __getattr__(self, name):
            return None
    return Args(output_dir=str(tmp_path), region_padding=100, adaptive_padding=False, max_intron=1000,
                goi_fallback_max_gap=0, goi_margin=5000, gap_search_window=500, goi_home_span=0, expansion=1.0)


def _mapped_hits(monkeypatch, tmp_path, window_hits, seq_len=120):
    """What the window hands to union coverage, for hits as run_augmented_search returns them."""
    seen = []

    def spy(hits, query_len):
        seen.append((sorted((h["qstart"], h["qend"]) for h in hits), query_len))
        raise _Sentinel()

    monkeypatch.setattr(isr, "_query_union_coverage", spy)
    monkeypatch.setattr(isr, "run_augmented_search", lambda *a, **k: list(window_hits))
    monkeypatch.setattr(isr, "record_block_row", lambda *a, **k: None, raising=False)
    q = "M" * seq_len
    db = {"GOI_x": {"id": "GOI_x", "seq": q, "header": "GOI_x"}}
    block = {"chrom": "c1", "start": 1000, "end": 2000, "genes_count": 2, "hits": []}
    with pytest.raises(_Sentinel):
        isr.process_region_block(0, block, [], {"c1": "ACGT" * 2000}, db, "G.fna", _window_args(tmp_path), "u", 1)
    return seen[0]


def _hit(qid, qstart, qend, start, end):
    return {"query": qid, "qstart": qstart, "qend": qend, "start": start, "end": end, "evalue": 1e-5,
            "pident": 60.0, "alnlen": qend - qstart + 1, "bits": 50.0, "strand": "+", "chrom": "c1"}


def test_a_windows_hit_positions_are_positions_in_the_protein(monkeypatch, tmp_path):
    """The defect, through the window's own mapping: a hit of residues 5-45 of the window at 71-120
    is query 75-115. The old mapping gave it 5-45, so it overlapped the first hit in the query, and
    the union read 45 of 120 residues instead of 86."""
    hits = [_hit("GOI_x|frag_1_60_120", 5, 45, 100, 223), _hit("GOI_x|slide_71_120_120", 5, 45, 400, 523)]
    positions, qlen = _mapped_hits(monkeypatch, tmp_path, hits)
    assert qlen == 120 and positions == [(5, 45), (75, 115)]


@pytest.mark.parametrize("qid,want", [
    ("GOI_x", (5, 45)),
    ("GOI_x|frag_1_120_120", (5, 45)),
    ("GOI_x|frag_61_120_120", (65, 105)),
    ("GOI_x|slide_1_50_120", (5, 45)),
    ("GOI_x|slide_26_75_120", (30, 70)),
    ("GOI_x|frag_2", (5, 45)),
])
def test_each_kind_of_query_gets_its_offset(monkeypatch, tmp_path, qid, want):
    positions, _ = _mapped_hits(monkeypatch, tmp_path, [_hit(qid, 5, 45, 100, 223)])
    assert positions == [want]
