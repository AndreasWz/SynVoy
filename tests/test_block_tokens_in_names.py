"""Fix item 26 (9 Nov 2026): a GOI name that holds a block-like token (`_b<d>_l<d>`, e.g. Ly6e_b3_l1).

The model ids are `<parent>|<genome>_b<block>_<l|r><locus>_<kind>`, and a seed-form tandem copy's id ends in
`_by_<parent>`. Readers that looked for the token in the whole id read the GOI's own name: models found with a
relative's query lost their parent (and stayed LOW and unjudged), the region scorer filed every model under the
name's block, the plot's richest-block step took every GOI model and no flanking gene, and the tree made a
phantom species of a seed-form copy. Each reader now reads the first token behind the first '|' (as
block_index_of does) or cuts the copy's parent first. Ordinary names read exactly as before.
"""
import os
import random
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import cluster_grs  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import sequence_utils as su  # noqa: E402
import synvoy_tree  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


@pytest.mark.parametrize("goi,tag,given", [("GOI_Ly6e_b3_l1", "Dsim_fna", False), ("GOI_Ly6e", "Dsim_fna", False),
                                           ("GOI_Ly6e_b3_l1", "Dsim_fna", True),
                                           ("GOI_Ly6e", "Sample_b2_l3_fna", True)])   # a genome named Sample_b2_l3.fna
def test_a_relatives_model_keeps_its_home_parent(goi, tag, given):
    pytest.importorskip("parasail")
    rng = random.Random(1)
    home = "M" + "".join(rng.choice(AA) for _ in range(120))
    rel = home[:60] + "".join(rng.choice(AA) for _ in range(61))
    parent = goi + "~rel"
    mid = f"{parent}|{tag}_b0_l1_exon_ann"
    attrs = (f"ID={mid};Name={parent};SynVoy_Parent={parent};EvidenceType=exon_annotation;Identity=80.0;Exons=1;"
             "BlockFlankingSupport=4")
    cand = {"start": 1000, "end": 1363,
            "gff": [f"chr1\texon_annotation\tmRNA\t1000\t1363\t80.0\t+\t.\t{attrs}",
                    f"chr1\texon_annotation\tCDS\t1000\t1363\t.\t+\t0\tID={mid}_CDS1;Parent={mid}"],
            "record": {"id": mid, "seq": home, "description": f"parent:{parent}"}}
    by_parent = {goi: {"id": goi, "seq": home},
                 parent: {"id": parent, "seq": rel, "model_id": goi + "|Dyak_fna_b2_l1_exon_ann"}}
    out = isr.rebase_relative_candidates([cand], [parent], by_parent, genome_tag=tag if given else None)
    a = isr._parse_gff_attr_string(out[0]["gff"][0].split("\t")[8])
    assert a["SynVoy_Parent"] == goi and a["ID"] == f"{goi}|{tag}_b0_r1_exon_ann"
    cds = isr._parse_gff_attr_string(out[0]["gff"][1].split("\t")[8])
    assert cds["Parent"] == f"{goi}|{tag}_b0_r1_exon_ann"
    assert out[0]["record"]["id"] == f"{goi}|{tag}_b0_r1_exon_ann"
    assert isr._model_genome_tag(a["ID"]) == tag and su.block_index_of(a["ID"], tag) == 0


@pytest.mark.parametrize("goi", ["GOI_Ly6e_b3_l1", "GOI_Ly6e"])
def test_the_region_scorer_reads_each_models_own_block(tmp_path, goi):
    rows = []
    for b, s in ((0, 10000), (5, 900000)):
        rows.append(f"chr1\texon_annotation\tmRNA\t{s}\t{s + 900}\t80\t+\t.\tID={goi}|Dsim_fna_b{b}_l1_exon_ann;"
                    f"Name={goi};SynVoy_Parent={goi};SynVoyRole=goi;Confidence=HIGH;Identity=80")
        for k in range(3):
            rows.append(f"chr1\tflanking_miniprot\tmRNA\t{s + 2000 + k * 3000}\t{s + 2900 + k * 3000}\t90\t+\t.\t"
                        f"ID=gene-f{b}{k}|Dsim_fna_b{b}_fl1_flank_ann;Name=gene-f{b}{k};SynVoy_Parent=gene-f{b}{k};"
                        "SynVoyRole=flanking;Confidence=HIGH;Identity=90")
    gff = tmp_path / "Dsim.fna.gff"
    gff.write_text("\n".join(rows) + "\n")
    for tag in ("Dsim_fna", None):
        anchors = cluster_grs.load_goi_block_anchors_from_gff(str(gff), padding_bp=0, genome_tag=tag)
        assert sorted((a["start"], a["end"], a["flanking_support"]) for a in anchors) == \
            [(10000, 18900, 3), (900000, 908900, 3)]
    assert cluster_grs.genome_tag_of("Dsim.fna") == "Dsim_fna"


@pytest.mark.parametrize("leaf,species,tag", [
    ("GOI_copy_3|Dsim_fna_b0_l1_by_GOI_Ly6e_b3_l1", "Dsim", "Dsim_fna"),
    ("GOI_copy_3|Dsim_fna_b0_l1_by_GOI_Ly6e", "Dsim", "Dsim_fna"),
    ("GOI_copy_3|Dsim_fna_b0_l1m2_by_GOI_x_b1_r2", "Dsim", "Dsim_fna"),
    ("GOI_Ly6e_b3_l1|Dsim_fna_b0_l1_exon_ann", "Dsim", "Dsim_fna"),
    ("GOI_x|Yeast_r1_fna_b2_l1_exon_ann", "Yeast_r1", "Yeast_r1_fna"),   # a genome tag that holds a token-like part
])
def test_a_seed_form_copy_names_its_own_genome(leaf, species, tag):
    assert synvoy_tree.species_from_leaf(leaf) == species
    assert isr._model_genome_tag(leaf) == tag


def test_strip_copy_parent():
    assert su.strip_copy_parent("Dsim_fna_b0_l1_by_GOI_Ly6e_b3_l1") == "Dsim_fna_b0_l1"
    assert su.strip_copy_parent("Dsim_fna_b0_l1_exon_ann") == "Dsim_fna_b0_l1_exon_ann"
    assert su.strip_copy_parent("") == ""


def test_the_plots_richest_block_step_is_gone():
    """It read the block from the Name attribute, which holds none for ordinary names, so it never
    acted; with a GOI named like Ly6e_b3_l1 it took every GOI model and no flanking gene."""
    src = open(os.path.join(ROOT, "bin", "plot_synteny.py")).read()
    assert "using richest block" not in src


def test_a_copy_parent_with_a_token_is_no_longer_warned_about(caplog):
    import logging
    isr._WARNED_COPY_PARENTS.clear()
    with caplog.at_level(logging.WARNING, logger=isr.logger.name):
        rid = isr.tandem_copy_record_id("GOI_copy_1", "GOI_Ly6e_b3_l1", "Dsim_fna", 0, 1)
    assert rid == "GOI_copy_1|Dsim_fna_b0_l1_by_GOI_Ly6e_b3_l1" and "tandem copy parent" not in caplog.text
    with caplog.at_level(logging.WARNING, logger=isr.logger.name):
        isr.tandem_copy_record_id("GOI_copy_1", "GOI_a|b", "Dsim_fna", 0, 1)
    assert "holds whitespace or '|'" in caplog.text
