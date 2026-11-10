"""Fix item 46 (9 Nov 2026): the tree does not read the report.

COMPUTE_TREE infers the tree from every HIGH or MEDIUM model of the search, before the report adjudicates them, so a
model the report merged into another, demoted, relabelled or overruled stayed a leaf, drawn as a copy of its species
(all four pieces of c7's luciferase case were leaves on Round 2 and r3). The tree stays the inferred tree; the
report lists the leaves it does not count (trees.<locus>.not_counted, from the leaves COMPUTE_TREE records), and the
tree page draws them hollow with the reason and numbers copies over the counted leaves only.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import generate_report as gr  # noqa: E402


def _row(mid, start, end, conf="HIGH", cls="confident_goi", source="locus_1", genome="G", chrom="c1"):
    return {"genome": genome, "source": source, "chrom": chrom, "start": start, "end": end, "mrna_id": mid,
            "confidence": conf, "goi_class": cls}


def _rec(mid, start, end, conf="HIGH", cls="confident_goi", prov=("locus_1",), genome="G", chrom="c1",
         members=None, overruled=(), **kw):
    """A record as dedupe_goi_annotations builds it, with its members (locus, model id, start, end), which the
    report step drops before it writes the report; by default its representative alone."""
    return dict({"genome": genome, "chrom": chrom, "start": start, "end": end, "mrna_id": mid, "confidence": conf,
                 "goi_class": cls, "provenance": list(prov), "_rep": (prov[0], mid),
                 "_members": list(members or [(prov[0], mid, start, end)]), "_overruled": list(overruled)}, **kw)


ROWS = [_row("GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000),
        _row("GOI_X|G_fna_b0_l1_extra1", 1100, 2900, conf="MEDIUM", cls="probable_goi"),   # merged into A
        _row("GOI_X|G_fna_b1_l1_exon_ann", 9000, 11000),                                    # relabelled
        _row("GOI_X|G_fna_b2_l1_extra1", 20500, 21000, conf="MEDIUM", cls="probable_goi"),  # overruled piece
        _row("GOI_X|G_fna_b2_l1_exon_ann", 20000, 24000, conf="AMBIGUOUS", cls="syntenic_candidate_unconfirmed"),
        _row("GOI_X|G_fna_b0_l1_exon_ann", 50000, 52000, source="locus_2")]                 # same id, other locus
RECORDS = [_rec("GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000,
                members=[("locus_1", "GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000),
                         ("locus_1", "GOI_X|G_fna_b0_l1_extra1", 1100, 2900)]),
           _rec("GOI_X|G_fna_b1_l1_exon_ann", 9000, 11000, cls="paralog_not_goi"),
           _rec("GOI_X|G_fna_b2_l1_exon_ann", 20000, 24000, conf="AMBIGUOUS", cls="syntenic_candidate_unconfirmed",
                overruled_counted_piece=True,
                members=[("locus_1", "GOI_X|G_fna_b2_l1_exon_ann", 20000, 24000),
                         ("locus_1", "GOI_X|G_fna_b2_l1_extra1", 20500, 21000)],
                overruled=[("locus_1", "GOI_X|G_fna_b2_l1_extra1")]),
           _rec("GOI_X|G_fna_b0_l1_exon_ann", 50000, 52000, prov=("locus_2",))]
LEAVES = ["GOI_X", "GOI_X|G_fna_b0_l1_exon_ann", "GOI_X|G_fna_b0_l1_extra1", "GOI_X|G_fna_b1_l1_exon_ann",
          "GOI_X|G_fna_b2_l1_extra1"]


def test_the_leaves_the_report_does_not_count():
    out = gr.tree_leaves_not_counted(ROWS, RECORDS, "locus_1", LEAVES)
    assert [(r["leaf"], r["reason"]) for r in out] == [
        ("GOI_X|G_fna_b0_l1_extra1", "merged"),
        ("GOI_X|G_fna_b1_l1_exon_ann", "paralog_not_goi"),
        ("GOI_X|G_fna_b2_l1_extra1", "overruled_piece")]
    assert out[0]["record"] == "G c1:1000-3000"


def test_a_leaf_is_judged_in_its_own_locus():
    """Model ids repeat across loci: locus_2's model of the same id is its own counted record."""
    assert gr.tree_leaves_not_counted(ROWS, RECORDS, "locus_2", ["GOI_X|G_fna_b0_l1_exon_ann"]) == []


def test_a_gene_counted_from_another_locus_is_counted_in_this_tree():
    """A record keeps one model id, its representative's; the same model found from locus_3 has another id
    there (the block index). In locus_3's tree that leaf is the counted gene, not 'merged' (r3: cow decorin and
    the luciferase cross-locus duplicates were drawn hollow)."""
    rows = [_row("GOI_X|G_fna_b1_l1_exon_ann", 1000, 3000), _row("GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000, source="locus_3")]
    rec = _rec("GOI_X|G_fna_b1_l1_exon_ann", 1000, 3000, cls="cross_locus_duplicate", prov=("locus_1", "locus_3"),
               members=[("locus_1", "GOI_X|G_fna_b1_l1_exon_ann", 1000, 3000),
                        ("locus_3", "GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000)])
    assert gr.tree_leaves_not_counted(rows, [rec], "locus_3", ["GOI_X", "GOI_X|G_fna_b0_l1_exon_ann"]) == []
    assert gr.tree_leaves_not_counted(rows, [rec], "locus_1", ["GOI_X", "GOI_X|G_fna_b1_l1_exon_ann"]) == []


def test_a_merged_model_says_merged_whatever_it_overlaps_most():
    """A wide model merged into a compact counted record (superseded_wide_model) overlaps a paralog's record
    more than its own; its reason is its record's, not the paralog's (r3 luciferase, Pyrocoelia)."""
    rows = [_row("GOI_X|G_fna_b0_l1_exon_ann", 185000, 187000), _row("GOI_X|G_fna_b0_l1m7_exon_ann", 154000, 255000),
            _row("GOI_X|G_fna_b3_l1_exon_ann", 239000, 290000)]
    counted = _rec("GOI_X|G_fna_b0_l1_exon_ann", 185000, 187000, superseded_wide_model=True,
                   members=[("locus_1", "GOI_X|G_fna_b0_l1_exon_ann", 185000, 187000),
                            ("locus_1", "GOI_X|G_fna_b0_l1m7_exon_ann", 154000, 255000)])
    paralog = _rec("GOI_X|G_fna_b3_l1_exon_ann", 239000, 290000, cls="paralog_not_goi")
    out = gr.tree_leaves_not_counted(rows, [counted, paralog], "locus_1",
                                     [r["mrna_id"] for r in rows])
    assert [(r["leaf"], r["reason"], r["record"]) for r in out] == [
        ("GOI_X|G_fna_b0_l1m7_exon_ann", "merged", "G c1:185000-187000"),
        ("GOI_X|G_fna_b3_l1_exon_ann", "paralog_not_goi", "G c1:239000-290000")]


def test_a_cross_locus_duplicate_through_the_report(tmp_path):
    """The same model from three home loci, one id each (the real dedup): every locus's tree counts its leaf,
    and the report keeps no member list."""
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    attrs = "SynVoyRole=goi;EvidenceType=exon_annotation;Identity=69.0;QueryCoverage=0.98;Bits=760.8"
    for locus, b in (("locus_1", 1), ("locus_3", 0), ("locus_4", 1)):
        (tmp_path / "regions" / f"{locus}__G.fna.gff").write_text(
            "##gff-version 3\n" f"c1\tx\tmRNA\t1000\t3000\t69\t+\t.\tID=GOI_X|G_fna_b{b}_l1_exon_ann;"
            f"SynVoy_Parent=GOI_X;Confidence=HIGH;GOIClass=confident_goi;{attrs}\n")
    infos = {loc: {"locus": loc, "leaves": ["GOI_X", f"GOI_X|G_fna_b{b}_l1_exon_ann"]}
             for loc, b in (("locus_1", 1), ("locus_3", 0), ("locus_4", 1))}
    report = gr.build_report(str(tmp_path), tree_infos=infos)
    (rec,) = report["goi_dedup"]["records"]
    assert rec["cross_locus_duplicate"] and rec["provenance"] == ["locus_1", "locus_3", "locus_4"]
    assert all(info["not_counted"] == [] for info in report["trees"].values())
    assert not any(k.startswith("_") for k in rec)


def test_the_report_lists_them_only_with_a_tree(tmp_path):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    attrs = "SynVoyRole=goi;EvidenceType=exon_annotation;Identity=80;QueryCoverage=0.95"
    (tmp_path / "regions" / "locus_1__G.fna.gff").write_text("##gff-version 3\n" + "".join(
        f"c1\tx\tmRNA\t{s}\t{e}\t80\t+\t.\tID={mid};SynVoy_Parent=GOI_X;Confidence={c};GOIClass={k};{attrs}\n"
        for mid, s, e, c, k in (("GOI_X|G_fna_b0_l1_exon_ann", 1000, 3000, "HIGH", "confident_goi"),
                                ("GOI_X|G_fna_b0_l1_extra1", 1100, 2900, "MEDIUM", "probable_goi"))))
    report = gr.build_report(str(tmp_path))
    assert "trees" not in report
    info = {"locus": "locus_1", "leaves": ["GOI_X", "GOI_X|G_fna_b0_l1_exon_ann", "GOI_X|G_fna_b0_l1_extra1"]}
    report = gr.build_report(str(tmp_path), tree_infos={"locus_1": info})
    assert report["trees"]["locus_1"]["not_counted"] == [
        {"leaf": "GOI_X|G_fna_b0_l1_extra1", "reason": "merged", "record": "G c1:1000-3000"}]
    assert any("trees.locus_1.not_counted" in w for w in report["warnings"])
    assert any(w.startswith("1 of the 2 model leaves of the locus_1 tree is a model the report does not count")
               and "draws it hollow" in w for w in report["warnings"])


def test_compute_tree_records_its_leaves():
    src = open(os.path.join(ROOT, "bin", "compute_tree.py")).read()
    assert "leaves=leaves" in src and 'leaves = [name.split()[0] for name, _seq in entries' in src


def test_the_tree_page_draws_them_hollow_and_numbers_the_counted(tmp_path):
    import plot_synteny as ps
    nwk = tmp_path / "locus_1_tree.nwk"
    nwk.write_text("(GOI_X:0.1,(GOI_X|Gen_a_fna_b0_l1_exon_ann:0.1,GOI_X|Gen_a_fna_b0_l1_extra1:0.2):0.1,"
                   "(GOI_X|Gen_b_fna_b0_l1_exon_ann:0.1,GOI_X|Gen_b_fna_b1_l1_exon_ann:0.3):0.2);\n")
    out = tmp_path / "synteny_block_locus_1_tree.html"
    ps._render_tree_svg(str(nwk), {}, str(out), species_map={}, gene_legend=False,
                        not_counted={"GOI_X|Gen_a_fna_b0_l1_extra1": "merged"})
    html = out.read_text()
    assert "Gen a</tspan> (not counted: merged)" in html and html.count('class="tree-leaf-not-counted"') == 1
    assert "Gen a</tspan> · copy" not in html             # one counted leaf left: no copy number
    assert "Gen b</tspan> · copy 1" in html and "Gen b</tspan> · copy 2" in html
    assert "the tree is inferred before the report" in html
    report = tmp_path / "synvoy_report.json"
    report.write_text(json.dumps({"trees": {"locus_1": {"not_counted": [
        {"leaf": "GOI_X|Gen_a_fna_b0_l1_extra1", "reason": "merged", "record": "x"}]}}}))
    assert ps._tree_not_counted(str(report), str(out)) == {"GOI_X|Gen_a_fna_b0_l1_extra1": "merged"}
    assert ps._tree_not_counted(str(report), "synteny_block_locus_2_tree.html") == {}


def test_the_tree_page_finds_a_leaf_iqtree_renamed(tmp_path):
    """IQ-TREE writes every byte outside [A-Za-z0-9_./|-] as '_' (IQ-TREE 3.1.2: '+', '(', '~', 'ä' as '__'),
    and the report's leaves are the FASTA's ids: the page looks a leaf up under both names."""
    import plot_synteny as ps
    report = tmp_path / "synvoy_report.json"
    report.write_text(json.dumps({"trees": {"locus_1": {"not_counted": [
        {"leaf": "GOI_Na+K-ATP(ä)|G_fna_b0_l1_extra1", "reason": "merged", "record": "x"}]}}}))
    marks = ps._tree_not_counted(str(report), str(tmp_path / "synteny_block_locus_1_tree.html"))
    assert marks["GOI_Na+K-ATP(ä)|G_fna_b0_l1_extra1"] == "merged"
    assert marks["GOI_Na_K-ATP____|G_fna_b0_l1_extra1"] == "merged"
    assert ps.iqtree_leaf_name("GOI_k/l.m-n|s_1") == "GOI_k/l.m-n|s_1"


def test_the_rule_is_iqtrees(tmp_path):
    import random
    import shutil
    import subprocess
    import plot_synteny as ps
    import pytest
    binary = shutil.which("iqtree3") or shutil.which("iqtree2") or shutil.which("iqtree")
    if binary is None:
        pytest.skip("IQ-TREE not on PATH")
    names = ["GOI_a+b|x_fna_b0_l1", "GOI_c(x)|y", "GOI_d:e~f|z", "GOI_g#1[2]|v", "GOI_k/l.m-n|s", "GOI_ä'|l"]
    rng = random.Random(3)
    base = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))
    with open(tmp_path / "aln.faa", "w") as fh:
        for i, n in enumerate(names):
            s = list(base)
            for _ in range(5 + 3 * i):
                s[rng.randrange(80)] = rng.choice("ACDEFGHIKLMNPQRSTVWY")
            fh.write(f">{n}\n{''.join(s)}\n")
    subprocess.run([binary, "-s", str(tmp_path / "aln.faa"), "-m", "LG", "-seed", "1", "-nt", "1", "-quiet", "-redo"],
                   check=True, capture_output=True, cwd=tmp_path)
    import re
    tree = (tmp_path / "aln.faa.treefile").read_text()
    written = set(re.findall(r"[(,]([^(),:;]+):", tree))
    assert written == {ps.iqtree_leaf_name(n) for n in names}
