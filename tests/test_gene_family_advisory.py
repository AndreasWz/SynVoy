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


def _staged_region(path, genome, chroms, goi=True):
    rows = ["##gff-version 3"]
    for k, chrom in enumerate(chroms):
        rows.append("\t".join([chrom, "x", "mRNA", str(1000 + 10 * k), str(2000 + 10 * k), "90", "+", ".",
                               f"ID=gene-f{k}|{genome}_b0_l{k};SynVoy_Parent=gene-f{k};SynVoyRole=flanking;"
                               "Confidence=HIGH;EvidenceType=exon_annotation;Identity=90"]))
    if goi:
        rows.append("\t".join([chroms[0], "x", "mRNA", "5000", "6000", "90", "+", ".",
                               f"ID=GOI_q|{genome}_b0_l1_exon_ann;SynVoy_Parent=GOI_q;SynVoyRole=goi;"
                               "Confidence=HIGH;GOIClass=confident_goi;EvidenceType=exon_annotation;"
                               "Identity=90;QueryCoverage=0.95"]))
    path.write_text("\n".join(rows) + "\n")


def test_the_run_counts_home_loci_not_cells(tmp_path):
    """Fix item 38: the run passed len(flanking_per_locus), one per (genome, staging source,
    chromosome) cell with a HIGH flanking gene: four genomes with flanking genes on three
    chromosomes each read as '12 home loci' for a gene with one."""
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    for g in ("G1", "G2", "G3", "G4"):
        _staged_region(tmp_path / "regions" / f"locus_1__{g}.fna.gff", g, ("c1", "c2", "c3"))
    adv = gr.build_report(str(tmp_path))["gene_family_advisory"]
    assert not any("home loci" in r for r in adv["reasons"]), adv["reasons"]
    for g in ("G1", "G2"):
        for locus in ("locus_2", "locus_3"):
            _staged_region(tmp_path / "regions" / f"{locus}__{g}.fna.gff", g, ("c9",))
    adv = gr.build_report(str(tmp_path))["gene_family_advisory"]
    assert "the query resolved to 3 home loci (paralog cluster)" in adv["reasons"]


def test_count_home_loci_reads_only_locus_tags():
    files = ["/x/locus_1__G1.fna.gff", "locus_1__G2.fna.gff", "locus_2__G1.hull_rescue.gff",
             "G3.fna.gff", "7__G4.fna.gff"]
    assert gr.count_home_loci(files) == 2
    assert gr.count_home_loci([]) == 0
