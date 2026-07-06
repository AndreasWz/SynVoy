"""Tests for the §1.1 gene-family detect-and-warn advisory in generate_report.

SynVoy is validated for divergent LOW-COPY genes via conserved flanking synteny; it is
not a general ortholog finder for large paralog families. The advisory says so honestly
when the query looks family-like, instead of asserting the GOI name on ambiguous calls.
"""
import importlib.util
import os

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")
spec = importlib.util.spec_from_file_location(
    "generate_report", os.path.join(BIN, "generate_report.py"))
gr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gr)


def test_single_copy_query_does_not_fire():
    adv = gr.build_gene_family_advisory(
        {"confident_goi": 1, "probable_goi": 3}, n_home_loci=1)
    assert adv["is_gene_family_query"] is False
    assert adv["reasons"] == []
    assert adv["advice"] == ""


def test_low_copy_toxin_with_some_ambiguous_does_not_fire():
    # Melittin-like: 1 home locus, ambiguous share below the ratio floor → NOT a family.
    counts = {"probable_goi": 11, "tandem_goi_copy": 47,
              "ambiguous_goi_family_member": 29, "confident_goi": 1}
    adv = gr.build_gene_family_advisory(counts, n_home_loci=1)
    assert adv["is_gene_family_query"] is False


def test_many_home_loci_fires():
    # Decorin/SLRP-like: query fanned out to a paralog cluster of home loci.
    adv = gr.build_gene_family_advisory(
        {"confident_goi": 1, "probable_goi": 4}, n_home_loci=5)
    assert adv["is_gene_family_query"] is True
    assert any("home loci" in r for r in adv["reasons"])
    assert "candidates" in adv["advice"].lower()


def test_high_ambiguous_share_fires():
    counts = {"ambiguous_goi_family_member": 29, "probable_goi": 5, "confident_goi": 1}
    adv = gr.build_gene_family_advisory(counts, n_home_loci=1)
    assert adv["is_gene_family_query"] is True
    assert any("ambiguous" in r for r in adv["reasons"])


def test_ratio_needs_both_floor_and_share():
    # 4 ambiguous of 5 total = 80% share but below the absolute floor (5) → no fire,
    # so a tiny run with one stray ambiguous call isn't mislabeled a family query.
    adv = gr.build_gene_family_advisory(
        {"ambiguous_goi_family_member": 4, "confident_goi": 1}, n_home_loci=1)
    assert adv["is_gene_family_query"] is False


def test_empty_counts_safe():
    adv = gr.build_gene_family_advisory({}, n_home_loci=0)
    assert adv["is_gene_family_query"] is False
