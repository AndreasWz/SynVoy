"""Bug list item 17: the report exported no LOW models (only counts of them). It lists them apart
now (goi_low_records): one entry per distinct model, over all home loci and rescue passes, with
what overlaps it; no count, label or flag reads them."""
import importlib.util
import os

import pytest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")
spec = importlib.util.spec_from_file_location("generate_report", os.path.join(BIN, "generate_report.py"))
gr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gr)


def _row(chrom, ftype, start, end, attrs, strand="+"):
    return "\t".join([chrom, "SynVoy", ftype, str(start), str(end), ".", strand, ".", attrs]) + "\n"


def _goi(start, end, conf, goi_class="ambiguous_goi_family_member", mid="GOI_q|G_fna_b0_l1_m1",
         chrom="c1", extra="", ftype="mRNA"):
    return _row(chrom, ftype, start, end,
                f"ID={mid};SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence={conf};GOIClass={goi_class};"
                f"EvidenceType=fallback_hit_span;Identity=31.5;QueryCoverage=0.42;Bits=33.1;Eeff=4.2;"
                f"Place=window;BlockEsyn=0.01;ModelStatus=partial;Basis=not_significant{extra}")


@pytest.fixture
def run(tmp_path):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    return tmp_path


def _write(run, name, rows):
    (run / "regions" / name).write_text("##gff-version 3\n" + "".join(rows))


def test_a_low_row_is_exported_with_its_fields(run):
    _write(run, "locus_1__G.fna.gff", [
        _goi(1001, 1600, "LOW"),
        _row("c1", "CDS", 1001, 1600, "ID=GOI_q|G_fna_b0_l1_m1.cds;Parent=GOI_q|G_fna_b0_l1_m1"),
        _goi(50001, 50600, "HIGH", "confident_goi", mid="GOI_q|G_fna_b1_l1_m1"),
        _goi(70001, 70600, "MEDIUM", "probable_goi", mid="GOI_q|G_fna_b2_l1_m1"),
        _goi(90001, 90600, "AMBIGUOUS", "syntenic_candidate_unconfirmed", mid="GOI_q|G_fna_b3_l1_m1"),
    ])
    block = gr.build_report(str(run))["goi_low_records"]
    assert block["count"] == 1 and len(block["records"]) == 1
    assert "counted nowhere" in block["note"]
    r = block["records"][0]
    assert (r["genome"], r["chrom"], r["start"], r["end"], r["strand"], r["confidence"]) == \
        ("G", "c1", 1001, 1600, "+", "LOW")
    assert (r["basis"], r["eeff"], r["place"], r["block_esyn"], r["bits"]) == \
        ("not_significant", "4.2", "window", "0.01", "33.1")
    assert (r["identity"], r["query_cov"], r["evidence"], r["model_status"]) == \
        ("31.5", "0.42", "fallback_hit_span", "partial")
    assert r["cds"] == [[1001, 1600]]
    assert r["provenance"] == ["locus_1"] and r["n_source_loci"] == 1
    assert r["overlapping_record"] is None
    assert r["genome_has_counted_call"] is True


def test_a_model_from_two_loci_is_one_record(run):
    _write(run, "locus_1__G.fna.gff", [_goi(1001, 1600, "LOW")])
    _write(run, "locus_2__G.fna.gff", [_goi(1001, 1600, "LOW", mid="GOI_q|G_fna_b0_l2_m1")])
    block = gr.build_report(str(run))["goi_low_records"]
    assert block["count"] == 1
    assert block["records"][0]["provenance"] == ["locus_1", "locus_2"]
    assert block["records"][0]["n_source_loci"] == 2


def test_a_hull_rescue_model_is_one_record_on_its_mrna(run):
    gid = "GOI_q|G_c1_hull_rescue"
    _write(run, "locus_1__G.hull_rescue.gff", [
        _row("c1", "gene", 1001, 3000, f"ID={gid};SynVoyRole=goi;Confidence=LOW;GOIClass=synteny_hull_rescue;"
                                       "EvidenceType=synteny_hull_rescue;Identity=41"),
        _row("c1", "mRNA", 1001, 3000, f"ID={gid}.mRNA;Parent={gid};SynVoyRole=goi;Confidence=LOW;"
                                       "GOIClass=synteny_hull_rescue;EvidenceType=synteny_hull_rescue;Identity=41"),
        _row("c1", "CDS", 1001, 1500, f"Parent={gid}.mRNA"),
        _row("c1", "CDS", 2001, 3000, f"Parent={gid}.mRNA"),
    ])
    report = gr.build_report(str(run))
    recs = report["goi_low_records"]["records"]
    assert len(recs) == 1
    assert recs[0]["mrna_id"] == f"{gid}.mRNA"
    assert recs[0]["cds"] == [[1001, 1500], [2001, 3000]]
    assert recs[0]["basis"] == ""                   # a rescue model is not judged
    assert report["summary"]["low_confidence_goi"] == 1


def test_overlapping_record_and_counted_genome(run):
    _write(run, "locus_1__G.fna.gff", [
        _goi(1001, 1600, "LOW"),
        _goi(1201, 1900, "AMBIGUOUS", "syntenic_candidate_unconfirmed", mid="GOI_q|G_fna_b0_l1_m2"),
    ])
    _write(run, "locus_1__H.fna.gff", [_goi(1001, 1600, "LOW", mid="GOI_q|H_fna_b0_l1_m1")])
    recs = {r["genome"]: r for r in gr.build_report(str(run))["goi_low_records"]["records"]}
    assert recs["G"]["overlapping_record"] == {"confidence": "AMBIGUOUS", "goi_class": "syntenic_candidate_unconfirmed",
                                               "start": 1201, "end": 1900, "counted": False}
    assert recs["G"]["genome_has_counted_call"] is False
    assert recs["H"]["overlapping_record"] is None


def test_the_export_changes_nothing_else(run, monkeypatch):
    _write(run, "locus_1__G.fna.gff", [_goi(1001, 1600, "LOW"),
                                       _goi(50001, 50600, "MEDIUM", "probable_goi", mid="GOI_q|G_fna_b1_l1_m1")])
    with_export = gr.build_report(str(run))
    monkeypatch.setattr(gr, "goi_low_records", lambda *a, **k: {})
    without = gr.build_report(str(run))
    assert with_export["summary"] == without["summary"]
    assert with_export["goi_dedup"] == without["goi_dedup"]
    assert [r["confidence"] for r in with_export["goi_dedup"]["records"]] == ["MEDIUM"]


def test_ownership_rows_are_matched_by_overlap():
    rows = [{"genome": "locus_1__G.fna", "chrom": "c1", "start_i": 1001, "end_i": 1600, "best_bit_f": 120.0,
             "gap_f": 40.0, "best_paralog": "home_gene-B|gene-B"},
            {"genome": "G.fna", "chrom": "c1", "start_i": 1001, "end_i": 1100, "best_bit_f": 300.0,
             "gap_f": 1.0, "best_paralog": "locus_1|gene-Q"}]
    ann = [{"genome": "G", "chrom": "c1", "start": 1001, "end": 1600, "strand": "+", "confidence": "LOW",
            "source": "locus_1", "mrna_id": "m1", "cds": []}]
    rec = gr.goi_low_records(ann, [], {"paralog_not_goi"}, rows)["records"][0]
    assert (rec["ownership_best_paralog"], rec["ownership_bit"], rec["ownership_gap"]) == \
        ("home_gene-B|gene-B", 120.0, 40.0)
    assert "ownership_bit" not in gr.goi_low_records(ann, [], set(), [])["records"][0]
