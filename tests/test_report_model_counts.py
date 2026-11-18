"""The report's model counts.

(a) the class-based counts keep their values for this release; goi_low_only_genomes is the same
    reading by label, and a genome's confidence counts hold no zero-valued labels of its flanking
    rows;
(b) a rescue model (a gene row and an mRNA row, the CDS rows on the mRNA) is one model, not two;
(c) total_raw_search_hits, genomes_with_raw_search_hits and synteny_hits_count are gone
    (always 0: the raw search hits do not reach the report);
(d) synvoy_label.BASIS_TEXT, read by nothing, is gone.
"""
import importlib.util
import os
import sys

import pytest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")
sys.path.insert(0, BIN)


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gr = _load("generate_report", "generate_report.py")
hull = _load("rescue_goi_hull", "rescue_goi_hull.py")
strong = _load("rescue_strong_synteny", "rescue_strong_synteny.py")

MRNA = ["w", "miniprot", "mRNA", "1", "900", "100", "+", ".", "ID=MP1;Identity=0.9100"]
CDS = [["w", "miniprot", "CDS", "1", "400", "100", "+", "0", "Parent=MP1"],
       ["w", "miniprot", "CDS", "501", "900", "100", "+", "2", "Parent=MP1"]]


def _row(chrom, ftype, start, end, attrs, strand="+"):
    return "\t".join([chrom, "SynVoy", ftype, str(start), str(end), ".", strand, ".", attrs])


@pytest.fixture
def run(tmp_path):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    return tmp_path


def _write(run, name, rows):
    (run / "regions" / name).write_text("##gff-version 3\n" + "\n".join(rows) + "\n")


def test_a_rescue_model_is_one_model(run):
    hull_rows = hull._build_hull_gff_rows(MRNA, CDS, "c1", 10_000, "GOI_q", "GOI_q|G_c1_hull_rescue", "HIGH",
                                          query_coverage=0.95, bits=310.0)
    strong_rows = strong._build_rescue_gff_rows(MRNA, CDS, "c2", 50_000, "GOI_q", "GOI_q|G_c2_rescue")
    assert [r.split("\t")[2] for r in hull_rows] == ["gene", "mRNA", "CDS", "CDS"]
    _write(run, "locus_1__G.hull_rescue.gff", hull_rows)
    _write(run, "locus_1__G.rescue.gff", strong_rows)
    report = gr.build_report(str(run))
    s = report["summary"]
    assert s["total_goi_annotations"] == 2
    assert s["high_confidence_goi_pre_dedup"] == 1
    assert s["low_confidence_goi"] == 1
    assert s["goi_hits_collapsed_by_dedup"] == 0
    assert report["annotations"]["goi_class_counts"] == {"synteny_hull_rescue": 1, "relaxed_rescue": 1}
    rec = report["goi_dedup"]["records"][0]
    assert rec["n_merged_hits"] == 1
    assert rec["mrna_id"] == "GOI_q|G_c1_hull_rescue.mRNA"
    low = report["goi_low_records"]["records"]
    assert [r["mrna_id"] for r in low] == ["GOI_q|G_c2_rescue.mRNA"]
    assert low[0]["cds"] == [[50_000, 50_399], [50_500, 50_899]]


def test_a_lone_tandem_gene_row_counts_once(run):
    _write(run, "locus_1__G.fna.gff", [
        _row("c1", "gene", 1001, 1300, "ID=GOI_copy_1|G_fna_b0_l1_by_GOI_q;SynVoy_Parent=GOI_q;SynVoyRole=goi;"
                                       "Confidence=MEDIUM;GOIClass=tandem_goi_copy;EvidenceType=tandem_copy;"
                                       "Identity=80;QueryCoverage=0.9"),
    ])
    s = gr.build_report(str(run))["summary"]
    assert s["total_goi_annotations"] == 1 and s["medium_confidence_goi"] == 1


def _goi(chrom, start, conf, goi_class, mid):
    return _row(chrom, "mRNA", start, start + 500,
                f"ID={mid};SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence={conf};GOIClass={goi_class};"
                "EvidenceType=fallback_hit_span;Identity=30;QueryCoverage=0.5")


def test_class_counts_keep_their_values_and_the_low_reading_by_label(run):
    _write(run, "locus_1__A.fna.gff", [
        _goi("c1", 1001, "LOW", "ambiguous_goi_family_member", "GOI_q|A_fna_b0_l1_m1"),
        _goi("c1", 3001, "LOW", "ambiguous_goi_family_member", "GOI_q|A_fna_b0_l1_m2"),
        _row("c1", "mRNA", 9001, 9500, "ID=gene-F1|A_fna_b0_l1;SynVoy_Parent=gene-F1;SynVoyRole=flanking;"
                                       "Confidence=HIGH;EvidenceType=exon_annotation;Identity=90"),
    ])
    _write(run, "locus_1__B.fna.gff", [
        _goi("c1", 1001, "LOW", "tandem_goi_copy", "GOI_copy_1|B_fna_b0_l1_by_GOI_q"),
        _goi("c1", 5001, "AMBIGUOUS", "syntenic_candidate_unconfirmed", "GOI_q|B_fna_b0_l1_m1"),
    ])
    report = gr.build_report(str(run))
    s = report["summary"]
    # The deprecated class-based keys, as before.
    assert s["ambiguous_goi_annotations"] == 2
    assert s["goi_ambiguous_only_genomes"] == ["A"]
    assert s["resolved_goi_annotations"] == 2
    # The same reading by label: B also holds an AMBIGUOUS record, so it is not LOW-only.
    assert s["goi_low_only_genomes"] == ["A"]
    assert report["annotations"]["genomes_with_only_low_goi"] == ["A"]
    per = {row["genome"]: row for row in report["annotations"]["per_genome"]}
    assert per["A"]["goi_confidence_counts"] == {"LOW": 2}       # no "HIGH": 0 of the flanking row
    assert per["B"]["goi_confidence_counts"] == {"LOW": 1, "AMBIGUOUS": 1}


def test_the_raw_hit_keys_are_gone(run):
    (run / "hits" / "x.m8").write_text("q\tc\t90\t10\t0\t0\t1\t10\t1\t10\t1e-9\t40\n")
    _write(run, "locus_1__A.fna.gff", [_goi("c1", 1001, "MEDIUM", "probable_goi", "GOI_q|A_fna_b0_l1_m1")])
    report = gr.build_report(str(run))
    assert "total_raw_search_hits" not in report["summary"]
    assert "genomes_with_raw_search_hits" not in report["summary"]
    assert report["synteny_results"] == {"genes_discovered": {}}
    assert "hit_files" not in report["staging_diagnostics"]["match_counts"]
    assert not hasattr(gr, "summarize_hits")


def test_basis_text_is_gone():
    import synvoy_label
    assert not hasattr(synvoy_label, "BASIS_TEXT")
