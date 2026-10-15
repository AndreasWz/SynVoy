"""The new GOI label in the search and the report (label rework, step 4, 2026-10-14).

Every GOI candidate is judged by synvoy_label (E-value in the expected place, specificity
against the home proteome, completeness): provisionally in its window, on its E in the
place, to choose among overlapping candidates; finally at the end of its genome, on Eeff.
Confidence keeps its four codes and Basis says why. The context promotion is gone: a weak
match in the slot of a significant block is AMBIGUOUS with Basis=gene_order_only (needs
curation, never counted). The report counts needs_curation_goi by basis and offers a
gene-order-only call only in a genome without a counted call.
"""
import os
import re
import sys

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

CODON = {"M": "ATG", "A": "GCT", "K": "AAA", "L": "CTT", "G": "GGT"}


def _goi_line(mid="l1", conf="LOW", ev="exon_annotation", extra="", start=101, end=400, strand="+"):
    attrs = (f"ID=GOI_q|g_b0_{mid};Name=GOI_q;SynVoy_Parent=GOI_q;Identity=31.0;SynVoyRole=goi;"
             f"EvidenceType={ev};Confidence={conf};GOIClass=ambiguous_goi_family_member;"
             f"InferenceReason=modeled_goi_but_family_context_is_weak;QueryCoverage=0.92" + extra)
    return f"chr1\texon_annotation\tmRNA\t{start}\t{end}\t31.0\t{strand}\t.\t{attrs}"


def _attrs(line):
    return isr._parse_gff_attr_string(line.split("\t")[8])


# ---------------------------------------------------------------- frame guard
def _model(cds_pieces, w_start=0, strand="+"):
    """A miniprot-like candidate whose CDS rows cover these (start, end, phase) pieces."""
    rows = [f"chr1\texon_annotation\tCDS\t{s}\t{e}\t.\t{strand}\t{ph}\tID=x_CDS{i};Parent=x"
            for i, (s, e, ph) in enumerate(cds_pieces, 1)]
    return {"start": cds_pieces[0][0], "end": cds_pieces[-1][1], "record": {"id": "x", "seq": "MKL"},
            "gff": [_goi_line(start=cds_pieces[0][0], end=cds_pieces[-1][1], strand=strand)] + rows}


def test_frame_guard_on_a_clean_model():
    window = "".join(CODON[a] for a in "MAKLGAKL") + "TAA" + "GGG" * 10
    cand = _model([(1, 27, 0)])
    isr.prepare_label_proteins([cand], window, 0)
    assert cand["label_protein"] == "MAKLGAKL"
    assert "Frame" not in cand["gff"][0]


def test_frame_guard_fixes_a_shifted_frame_and_marks_a_broken_one():
    window = "C" + "".join(CODON[a] for a in "MAKLGAKL") + "GGG"
    cand = _model([(1, 28, 0)])                    # one base too many in front
    isr.prepare_label_proteins([cand], window, 0)
    assert cand["label_protein"] == "MAKLGAKLG" and "FrameFixed=1" in cand["gff"][0]
    broken = "TAAATAAATAAAGGGCCC"
    cand = _model([(1, 18, 0)])
    isr.prepare_label_proteins([cand], broken, 0)
    assert cand["label_protein"] == "" and "FrameBroken=1" in cand["gff"][0]


def test_hit_chains_keep_their_translated_hits():
    cand = _model([(1, 27, 0)])
    cand["gff"][0] = cand["gff"][0].replace("EvidenceType=exon_annotation", "EvidenceType=fallback_hit_span")
    isr.prepare_label_proteins([cand], "A" * 60, 0)
    assert cand["label_protein"] == "MKL"


# ---------------------------------------------------------------- the window's choice
def _judged(mid, e, start, end, margin="40.0", real="1", qcov="0.9", place="slot", bits=100.0):
    extra = (f";Place={place};BlockEsyn=1e-20;EPlace={e};SpecMargin={margin};RealCompetitor={real};"
             f"BestOtherInNeighbourhood=0;QueryCoverage={qcov}")
    line = _goi_line(mid=mid, extra=extra, start=start, end=end).replace(";QueryCoverage=0.92", "")
    return {"start": start, "end": end, "bits": bits, "record": {"id": f"GOI_q|g_b0_{mid}", "seq": "M"},
            "gff": [line]}


def test_window_choice_goes_by_label_then_e():
    specific_weak = _judged("a", "1e-5", 100, 900)                    # MEDIUM
    relative_strong = _judged("b", "1e-40", 120, 880, margin="-50")   # AMBIGUOUS (not specific)
    isr.judge_candidates([specific_weak, relative_strong])
    assert [c["record"]["id"] for c in isr.greedy_select([relative_strong, specific_weak])] == ["GOI_q|g_b0_a"]
    strong = _judged("c", "1e-30", 100, 900)
    weak = _judged("d", "1e-8", 150, 850, bits=500.0)
    isr.judge_candidates([strong, weak])
    assert isr.greedy_select([weak, strong])[0] is strong


# ---------------------------------------------------------------- the final label
def test_relabel_on_eeff():
    lines = [
        _goi_line("h", extra=";Place=slot;BlockEsyn=1e-20;Eeff=1e-40;SpecMargin=200;RealCompetitor=1;BestOtherInNeighbourhood=0"),
        _goi_line("m", extra=";Place=core;BlockEsyn=1e-20;Eeff=2e-4;SpecMargin=12;RealCompetitor=1;BestOtherInNeighbourhood=0"),
        _goi_line("n", extra=";Place=slot;BlockEsyn=1e-20;Eeff=1e-30;SpecMargin=-40;RealCompetitor=1;BestOtherInNeighbourhood=0"),
        _goi_line("g", extra=";Place=slot;BlockEsyn=1e-20;Eeff=0.2;SpecMargin=8;RealCompetitor=1;BestOtherInNeighbourhood=0"),
        _goi_line("l", conf="MEDIUM", extra=";Place=genome;BlockEsyn=3;Eeff=0.2;SpecMargin=8;RealCompetitor=0"),
        _goi_line("t", ev="tandem_copy", extra=";Place=slot;BlockEsyn=1e-20;Eeff=1e-12;SpecMargin=60;RealCompetitor=1;BestOtherInNeighbourhood=0"),
        _goi_line("old"),                                     # no Eeff: left as it is
    ]
    out, counts = isr.relabel_goi_lines(lines)
    got = [(_attrs(x)["Confidence"], _attrs(x).get("Basis"), _attrs(x)["GOIClass"]) for x in out]
    assert got == [("HIGH", "significant_specific", "confident_goi"),
                   ("MEDIUM", "significant_specific", "probable_goi"),
                   ("AMBIGUOUS", "not_specific", "syntenic_candidate_unconfirmed"),
                   ("AMBIGUOUS", "gene_order_only", "syntenic_candidate_unconfirmed"),
                   ("LOW", "not_significant", "ambiguous_goi_family_member"),
                   ("HIGH", "significant_specific", "tandem_goi_copy"),
                   ("LOW", None, "ambiguous_goi_family_member")]
    assert _attrs(out[1])["InferenceReason"] == "significant_specific_in_core"
    assert counts["AMBIGUOUS/gene_order_only"] == 1
    again, _ = isr.relabel_goi_lines(out)                 # idempotent: one Basis per row
    assert again == out


def test_ml_rescoring_keeps_a_judged_label():
    line = _goi_line(extra=";Basis=not_significant;Eeff=3")
    new, boosted = isr._apply_ml_scores_to_mrna_line(line, 0.95, None)
    assert boosted == 0 and _attrs(new)["Confidence"] == "LOW" and "EmbeddingSimilarity=0.950" in new


def test_the_context_promotion_is_gone():
    config = open(os.path.join(REPO, "nextflow.config")).read()
    assert "classify_context_min" not in config
    assert not hasattr(isr, "_promote_context_supported_goi")
    main_nf = open(os.path.join(REPO, "main.nf")).read()
    assert re.search(r"classify_context_min_collinear\s*:", main_nf)          # replaced: the run stops
    assert "'classify_context_min_identity', 'classify_context_min_qcov'" in main_nf


# ---------------------------------------------------------------- the report
def _rec(genome, conf, basis="significant_specific", eeff="1e-20", start=100, end=900, source="locus_1",
         bits="200", esyn="1e-20"):
    return {"genome": genome, "source": source, "chrom": "chr1", "start": start, "end": end,
            "confidence": conf, "goi_class": isr.GOI_CLASS_OF_LABEL[conf], "identity": "40.0",
            "query_cov": "0.9", "bits": bits, "basis": basis, "eeff": eeff, "place": "slot",
            "block_esyn": esyn, "target_gene": "", "model_status": "complete", "mrna_id": f"m{start}"}


def test_report_keeps_the_record_the_judge_ranks_first():
    weak = _rec("g1", "HIGH", eeff="1e-8", bits="900")
    strong = _rec("g1", "HIGH", eeff="1e-30", bits="100", start=110, end=890)
    rep = gr.dedupe_goi_annotations([weak, strong])["records"][0]
    assert rep["eeff"] == "1e-30" and rep["basis"] == "significant_specific"


def test_gene_order_only_is_offered_only_where_nothing_was_counted():
    recs = [_rec("g1", "MEDIUM"), _rec("g1", "AMBIGUOUS", basis="gene_order_only", eeff="0.3", start=5000, end=5400),
            _rec("g2", "AMBIGUOUS", basis="gene_order_only", eeff="0.3"),
            _rec("g3", "AMBIGUOUS", basis="not_specific", eeff="1e-9")]
    out = gr.dedupe_goi_annotations(recs)
    kept = sorted((r["genome"], r["basis"]) for r in out["records"])
    assert kept == [("g1", "significant_specific"), ("g2", "gene_order_only"), ("g3", "not_specific")]
    assert out["gene_order_only_dropped"] == 1


def test_headline_names_what_needs_curation():
    h = gr._format_goi_headline(9, 3, 40, 33, ambiguous=4,
                                curation_by_basis={"gene_order_only": 3, "not_specific": 1})
    assert h == ("12 ortholog(s) (9 high, 3 medium confidence); 4 need curation (3 gene order only, "
                 "1 another home gene fits as well); 40 no significant match; across 33 genome(s).")
    assert gr._format_goi_headline(0, 0, 5, 3, ambiguous=2).startswith("NO ortholog; 2 need curation")
