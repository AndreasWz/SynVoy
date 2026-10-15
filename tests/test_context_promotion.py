"""Context promotion: a weak GOI model on a block with conserved gene order becomes MEDIUM
when the genome has no confident call for the query (iterative_search_runner
_promote_context_supported_goi). Off by default; the report names such calls apart
(generate_report: record field context_promoted, summary gene_order_only_goi)."""
import os
import re
import sys

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402


def _goi(mid, conf="LOW", ev="fallback_hit_span", ident=25.0, qcov=0.8, coll=6, parent="GOI_q1"):
    attrs = (f"ID={mid};Name={parent};SynVoy_Parent={parent};Identity={ident};SynVoyRole=goi;"
             f"EvidenceType={ev};Confidence={conf};GOIClass=ambiguous_goi_family_member;"
             f"BlockFlankingSupport={coll + 1};BlockCollinearSupport={coll};"
             f"InferenceReason=fallback_span_only;QueryCoverage={qcov}")
    return f"chr1\tfallback_hits\tmRNA\t100\t400\t{ident}\t+\t.\t{attrs}"


def _cds(mid):
    return f"chr1\tfallback_hits\tCDS\t100\t400\t.\t+\t0\tID={mid}_CDS1;Parent={mid}"


def _tandem(mid, conf="MEDIUM", parent="GOI_q1"):
    """A tandem copy as the search writes it: a `gene` row, not an mRNA."""
    attrs = (f"ID={mid};Name=GOI_copy_1;SynVoy_Parent={parent};Type=tandem_copy;Identity=55.0;"
             f"SynVoyRole=goi;EvidenceType=tandem_copy;Confidence={conf};GOIClass=tandem_goi_copy;"
             f"BlockFlankingSupport=9;BlockCollinearSupport=9;QueryCoverage=0.6")
    return f"chr1\ttandem\tgene\t5000\t5200\t55.0\t+\t.\t{attrs}"


def _attrs(line):
    return isr._parse_gff_attributes(line.split("\t")[8])


@pytest.fixture(autouse=True)
def _thresholds(monkeypatch):
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "context_min_collinear", 4)
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "context_min_identity", 15.0)
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "context_min_qcov", 0.5)
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "fragment_max_qcov", 0.4)


def test_the_promotion_is_off_as_shipped():
    """Config, script option and threshold table agree on 0: on for one day from 2026-10-14,
    switched off again because decoys meet the rule."""
    config = open(os.path.join(REPO, "nextflow.config")).read()
    assert re.search(r"^\s*classify_context_min_collinear\s*=\s*0\b", config, re.M)
    script = open(os.path.join(REPO, "bin", "iterative_search_runner.py")).read()
    assert re.search(r'"--classify_context_min_collinear", type=lambda v: int\(float\(v\)\), default=0\b', script)
    assert re.search(r'^\s*"context_min_collinear": 0,', script, re.M)


def test_best_ordered_model_is_promoted_once():
    lines = [_goi("m1", coll=5), _cds("m1"), _goi("m2", coll=9, ident=22.0), _goi("m3", coll=9, ident=30.0)]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 1
    a = _attrs(out[3])
    assert a["Confidence"] == "MEDIUM"
    assert a["GOIClass"] == "probable_goi"
    assert a["InferenceReason"] == "low_identity_in_conserved_gene_order"
    assert a["ContextPromoted"] == "LOW"
    assert _attrs(out[0])["Confidence"] == "LOW" and _attrs(out[2])["Confidence"] == "LOW"
    assert out[1] == lines[1]


def test_no_promotion_when_genome_has_a_confident_call():
    lines = [_goi("m1", coll=12), _goi("m2", conf="MEDIUM", ident=45.0, coll=2)]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 0 and out == lines


def test_each_query_is_judged_on_its_own():
    lines = [_goi("a1", coll=8, parent="GOI_a"), _goi("b1", conf="HIGH", parent="GOI_b"),
             _goi("b2", coll=8, parent="GOI_b")]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 1
    assert _attrs(out[0])["Confidence"] == "MEDIUM" and _attrs(out[2])["Confidence"] == "LOW"


@pytest.mark.parametrize("kw", [
    {"coll": 3},                    # too few anchors in home order
    {"qcov": 0.3},                  # fragment
    {"qcov": 0.45},                 # not a fragment, below the coverage floor
    {"ident": 12.0},                # below the identity floor
    {"ev": "rescued_exon"},         # isolated exon
    {"ev": "tandem_copy"},
])
def test_not_promoted(kw):
    out, n = isr._promote_context_supported_goi([_goi("m1", **kw)])
    assert n == 0


def test_fragment_bound_above_the_floor_wins(monkeypatch):
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "fragment_max_qcov", 0.6)
    assert isr._promote_context_supported_goi([_goi("m1", qcov=0.55)])[1] == 0
    assert isr._promote_context_supported_goi([_goi("m1", qcov=0.65)])[1] == 1


def test_a_confident_tandem_copy_counts_as_a_call():
    """Tandem copies are `gene` rows. A genome whose confident call is one has a call."""
    lines = [_goi("m1", coll=12), _tandem("t1", conf="MEDIUM")]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 0 and out == lines
    # a LOW tandem copy is no call, and it is never the model that is promoted
    lines = [_tandem("t1", conf="LOW"), _goi("m1", coll=12)]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 1
    assert _attrs(out[0])["Confidence"] == "LOW" and _attrs(out[1])["Confidence"] == "MEDIUM"
    # a tandem copy of another query does not stand in for this one
    lines = [_goi("m1", coll=12), _tandem("t1", conf="MEDIUM", parent="GOI_other")]
    assert isr._promote_context_supported_goi(lines)[1] == 1


def test_ambiguous_is_promoted_and_switch_off(monkeypatch):
    lines = [_goi("m1", conf="AMBIGUOUS", ident=33.0, coll=10)]
    out, n = isr._promote_context_supported_goi(lines)
    assert n == 1 and _attrs(out[0])["ContextPromoted"] == "AMBIGUOUS"
    monkeypatch.setitem(isr.CLASSIFY_THRESHOLDS, "context_min_collinear", 0)
    assert isr._promote_context_supported_goi(lines) == (lines, 0)


def test_ml_rescoring_keeps_a_promoted_call():
    line, _ = isr._promote_context_supported_goi([_goi("m1", ident=22.0, coll=10)])
    new, _boosted = isr._apply_ml_scores_to_mrna_line(line[0], None, None)
    assert _attrs(new)["Confidence"] == "MEDIUM"


def test_promoted_call_feeds_the_tree_but_never_seeds():
    genes = [{"id": "GOI_q1|g_b0_l1", "seq": "M" * 50}, {"id": "GOI_q2|g_b1_l1", "seq": "M" * 50}]
    meta = {
        "GOI_q1|g_b0_l1": {"role": "goi", "confidence": "MEDIUM", "goi_class": "probable_goi",
                            "block_flanking_support": "9", "query_coverage": "0.8",
                            "context_promoted": "LOW"},
        "GOI_q2|g_b1_l1": {"role": "goi", "confidence": "MEDIUM", "goi_class": "probable_goi",
                            "block_flanking_support": "9", "query_coverage": "0.8",
                            "context_promoted": ""},
    }
    seed, tree_extra, suppressed = isr._classify_goi_for_seed_and_tree(
        genes, meta, lambda gid: gid.startswith("GOI_"))
    assert [g["id"] for g in seed] == ["GOI_q2|g_b1_l1"]
    assert [g["id"] for g in tree_extra] == ["GOI_q1|g_b0_l1"]
    assert suppressed == 1


# --- the report -----------------------------------------------------------------------

def _record(conf="MEDIUM", promoted="", source="locus_1", start=100, end=400, ident="25.0"):
    return {"genome": "g1", "source": source, "chrom": "chr1", "start": start, "end": end,
            "confidence": conf, "goi_class": "probable_goi", "identity": ident, "query_cov": "0.8",
            "target_gene": "", "model_status": "partial", "mrna_id": f"m_{source}_{start}",
            "context_promoted": promoted}


def test_report_record_carries_the_old_label():
    recs = gr.dedupe_goi_annotations([_record(promoted="LOW"),
                                      _record(start=9000, end=9400, ident="45.0")])["records"]
    by_start = {r["start"]: r for r in recs}
    assert by_start[100]["context_promoted"] == "LOW"
    assert by_start[9000]["context_promoted"] == ""


def test_same_gene_called_on_sequence_from_another_locus_is_not_order_only():
    recs = gr.dedupe_goi_annotations([
        _record(promoted="AMBIGUOUS", source="locus_1", ident="53.0"),
        _record(promoted="", source="locus_2", ident="45.0"),
    ])["records"]
    assert len(recs) == 1 and recs[0]["context_promoted"] == ""


def test_headline_names_the_calls_by_gene_order():
    plain = gr._format_goi_headline(3, 2, 4, 5, ambiguous=1)
    named = gr._format_goi_headline(3, 2, 4, 5, ambiguous=1, order_only=1)
    assert "gene order" not in plain
    assert "3 high-confidence + 2 medium-confidence" in named
    assert "1 of the medium-confidence by gene order only" in named
    assert named.replace(", 1 of the medium-confidence by gene order only", "") == plain


def test_report_counts_calls_by_gene_order(tmp_path):
    """End to end: a promoted row in a result GFF is one of the MEDIUM calls of the summary
    and is counted apart in gene_order_only_goi."""
    promoted, _ = isr._promote_context_supported_goi([_goi("GOI_q1|g1_b0_l1", ident=25.0, coll=9)])
    regular = _goi("GOI_q1|g2_b0_l1", conf="MEDIUM", ident=45.0).replace("chr1\t", "chr2\t", 1)
    regions = tmp_path / "regions"
    regions.mkdir()
    (regions / "locus_1__g1.fna.gff").write_text("##gff-version 3\n" + promoted[0] + "\n")
    (regions / "locus_1__g2.fna.gff").write_text("##gff-version 3\n" + regular + "\n")
    report = gr.build_report(str(tmp_path))
    summary = report["summary"]
    assert summary["medium_confidence_goi"] == 2
    assert summary["gene_order_only_goi"] == 1
    assert "1 of the medium-confidence by gene order only" in summary["headline"]
    flagged = [r for r in report["goi_dedup"]["records"] if r["context_promoted"]]
    assert [r["genome"] for r in flagged] == ["g1"] and flagged[0]["context_promoted"] == "LOW"
