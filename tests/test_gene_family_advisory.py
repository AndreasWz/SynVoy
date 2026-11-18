"""Tests for the gene-family detect-and-warn advisory in generate_report.

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


def _rec(confidence="HIGH", goi_class="confident_goi", basis="significant_specific", **kw):
    return {"genome": kw.pop("genome", "G1"), "chrom": "c1", "start": 1, "end": 2,
            "confidence": confidence, "goi_class": goi_class, "basis": basis, **kw}


def test_single_copy_query_does_not_fire():
    adv = gr.build_gene_family_advisory([_rec(), _rec("MEDIUM", "probable_goi")] * 2, n_home_loci=1)
    assert adv["is_gene_family_query"] is False
    assert adv["reasons"] == []
    assert adv["advice"] == ""


def test_many_home_loci_fires():
    # Decorin/SLRP-like: query fanned out to a paralog cluster of home loci.
    adv = gr.build_gene_family_advisory([_rec(), _rec("MEDIUM", "probable_goi")], n_home_loci=5)
    assert adv["is_gene_family_query"] is True
    assert any("home loci" in r for r in adv["reasons"])
    assert "candidates" in adv["advice"].lower()


def test_records_another_home_gene_explains_fire():
    """5 of 10 records given to another home gene (paralog check): the family trigger fires,
    with the reason in records."""
    recs = [_rec(goi_class="paralog_not_goi") for _ in range(5)] + [_rec() for _ in range(5)]
    adv = gr.build_gene_family_advisory(recs, n_home_loci=1)
    assert adv["is_gene_family_query"] is True
    assert adv["reasons"] == ["5 of 10 GOI records (50%) are explained by another home gene at least as "
                              "well (paralog check or specificity)"]


def test_the_floor_of_five_records():
    # 4 of 5 = 80 % but below the absolute floor (5): a small run is no family query.
    recs = [_rec(goi_class="paralog_not_goi") for _ in range(4)] + [_rec()]
    assert gr.build_gene_family_advisory(recs, n_home_loci=1)["is_gene_family_query"] is False


def test_the_share_of_forty_percent():
    recs = [_rec(goi_class="paralog_not_goi") for _ in range(5)] + [_rec() for _ in range(8)]
    assert gr.build_gene_family_advisory(recs, n_home_loci=1)["is_gene_family_query"] is False


def test_not_specific_counts_gene_order_only_does_not():
    not_specific = [_rec("AMBIGUOUS", "syntenic_candidate_unconfirmed", "not_specific") for _ in range(5)]
    adv = gr.build_gene_family_advisory(not_specific, n_home_loci=1)
    assert adv["is_gene_family_query"] is True
    gene_order = [_rec("AMBIGUOUS", "syntenic_candidate_unconfirmed", "gene_order_only") for _ in range(5)]
    assert gr.build_gene_family_advisory(gene_order, n_home_loci=1)["is_gene_family_query"] is False


def test_answers_of_the_genome_wide_path_are_left_out():
    recs = [_rec(goi_class="paralog_not_goi") for _ in range(5)]
    recs += [_rec("MEDIUM", "genome_wide_strong_hit", genome_wide_strong_hit=True) for _ in range(20)]
    adv = gr.build_gene_family_advisory(recs, n_home_loci=1)
    assert adv["reasons"][0].startswith("5 of 5 GOI records")


def test_empty_records_safe():
    adv = gr.build_gene_family_advisory([], n_home_loci=0)
    assert adv["is_gene_family_query"] is False
    assert gr.build_gene_family_advisory(None)["is_gene_family_query"] is False


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


def test_low_rows_of_a_single_copy_query_do_not_fire(tmp_path):
    """Six LOW rows and one HIGH call of a one-locus query: the old trigger counted the LOW class
    ambiguous_goi_family_member and fired ("6/7 GOI calls ... are ambiguous family members")."""
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    rows = ["##gff-version 3"]
    for k in range(6):
        rows.append("\t".join(["c1", "x", "mRNA", str(10000 + 1000 * k), str(10500 + 1000 * k), "30", "+", ".",
                               f"ID=GOI_q|G1_b0_l1_m{k};SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=LOW;"
                               "GOIClass=ambiguous_goi_family_member;EvidenceType=fallback_hit_span;Identity=30;"
                               "Basis=not_significant"]))
    rows.append("\t".join(["c1", "x", "mRNA", "5000", "6000", "90", "+", ".",
                           "ID=GOI_q|G1_b0_l1_exon_ann;SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=HIGH;"
                           "GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=90;QueryCoverage=0.95"]))
    (tmp_path / "regions" / "locus_1__G1.fna.gff").write_text("\n".join(rows) + "\n")
    report = gr.build_report(str(tmp_path))
    assert report["gene_family_advisory"]["is_gene_family_query"] is False
    assert report["summary"]["gene_family_query"] is False
    assert report["summary"]["low_confidence_goi"] == 6
