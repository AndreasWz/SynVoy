"""Bug list item 9: the region tables said "dispersed macrosynteny, MEDIUM" where no gene was
modelled. The report lists every gene-bearing region of MEDIUM or better (and every
goi_hit_no_call region) that holds no counted call, with what it holds and why
(regions.goi_regions_without_call, summary.goi_regions_without_call)."""
import importlib.util
import json
import os

import pytest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gr = _load("generate_report", "generate_report.py")
isr = _load("iterative_search_runner", "iterative_search_runner.py")

SCORE_COLS = ["region_rank", "region_name", "species", "chrom", "start", "end", "strand", "score",
              "goi_overlap", "is_goi_anchor", "region_class", "confidence", "selection_reason",
              "goi_model_confidence", "goi_model_class", "goi_hit_best_identity", "goi_hit_best_evalue"]


def _region(name, chrom, start, end, conf="MEDIUM", cls="goi_overlap", overlap=True, score=0.45,
            model="", hit_id="", hit_e=""):
    return {"region_rank": 1, "region_name": name, "species": "", "chrom": chrom, "start": start,
            "end": end, "strand": "+", "score": score, "goi_overlap": overlap, "is_goi_anchor": False,
            "region_class": cls, "confidence": conf, "selection_reason": "score_floor",
            "goi_model_confidence": model, "goi_model_class": "", "goi_hit_best_identity": hit_id,
            "goi_hit_best_evalue": hit_e}


def _scores(tmp, name, rows, cols=SCORE_COLS):
    with open(tmp / "scores" / name, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _gff(tmp, name, rows):
    with open(tmp / "regions" / name, "w") as fh:
        fh.write("##gff-version 3\n" + "".join(rows))


def _goi(chrom, start, end, conf, goi_class, identity=70, qcov=0.9, mid="GOI_q|G_fna_b0_l1_m1"):
    return "\t".join([chrom, "x", "mRNA", str(start), str(end), ".", "+", ".",
                      f"ID={mid};SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence={conf};GOIClass={goi_class};"
                      f"EvidenceType=exon_annotation;Identity={identity};QueryCoverage={qcov}"]) + "\n"


def _refused(tmp, name, rows):
    cols = isr.NONSYNTENIC_COLUMNS
    with open(tmp / "hits" / f"{name}{gr.NONSYNTENIC_SUFFIX}", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _hit(genome, chrom, start, identity=62.0, outcome="refused", rejected_by="chromosome_lacks_min_anchors",
         locus="locus_1"):
    return {"genome": genome, "query": "GOI_q", "chrom": chrom, "start": start, "end": start + 600,
            "strand": "+", "qstart": 1, "qend": 200, "bits": 150.0, "identity": identity, "evalue": 1e-40,
            "alnlen": 600, "rejected_by": rejected_by, "outcome": outcome, "locus": locus}


@pytest.fixture
def run(tmp_path):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    return tmp_path


def _listed(run):
    report = gr.build_report(str(run))
    json.loads(json.dumps(report))
    return report, report["regions"]["goi_regions_without_call"]


def test_a_medium_region_with_a_refused_hit_and_no_model(run):
    _scores(run, "locus_1__G.fna.scores.tsv",
            [_region("Reg1_NOCALL_CMEDIUM_S0.45", "c1", 9_000, 12_000, cls="goi_hit_no_call",
                     hit_id=62.0, hit_e=1e-40)])
    _refused(run, "locus_1_G.fna", [_hit("G.fna", "c1", 10_000)])
    _gff(run, "locus_1__G.fna.gff", [_goi("c9", 1, 600, "MEDIUM", "probable_goi")])
    report, listed = _listed(run)
    assert len(listed) == 1
    e = listed[0]
    assert (e["genome"], e["locus"], e["region_class"], e["reason"]) == \
        ("G", "locus_1", "goi_hit_no_call", "refused_by_synteny_gate")
    assert e["model"] is None
    assert e["refused_hit"]["rejected_by"] == "chromosome_lacks_min_anchors"
    assert e["refused_hit"]["identity"] == 62.0
    assert e["goi_hit_best_identity"] == 62.0
    assert report["summary"]["goi_regions_without_call"] == 1


def test_a_counted_record_in_the_region_is_not_listed(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000, conf="HIGH")])
    # The counted call comes from another home locus.
    _gff(run, "locus_2__G.fna.gff", [_goi("c1", 10_001, 10_600, "MEDIUM", "probable_goi")])
    _, listed = _listed(run)
    assert listed == []


def test_only_a_low_model(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)])
    _gff(run, "locus_1__G.fna.gff", [_goi("c1", 10_001, 10_600, "LOW", "ambiguous_goi_family_member", 30)])
    _, listed = _listed(run)
    assert [(e["reason"], e["model"]["confidence"]) for e in listed] == [("model_not_counted", "LOW")]


def test_a_demoted_call(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000, conf="HIGH")])
    _gff(run, "locus_1__G.fna.gff", [_goi("c1", 10_001, 10_600, "HIGH", "confident_goi", 60, 0.2)])
    _, listed = _listed(run)
    assert [(e["reason"], e["model"]["goi_class"]) for e in listed] == \
        [("call_not_counted", "identity_coverage_decoupled")]


def test_a_low_region_of_another_class_is_not_listed(run):
    _scores(run, "locus_1__G.fna.scores.tsv",
            [_region("Reg1", "c1", 9_000, 12_000, conf="LOW", score=0.1),
             _region("Reg2", "c2", 9_000, 12_000, conf="HIGH", overlap=False, cls="synteny_block")])
    _, listed = _listed(run)
    assert listed == []


def test_a_low_no_call_region_is_listed(run):
    _scores(run, "locus_1__G.fna.scores.tsv",
            [_region("Reg1_NOCALL_CLOW_S0.45", "c1", 9_000, 12_000, conf="LOW", cls="goi_hit_no_call")])
    _refused(run, "locus_1_G.fna", [_hit("G.fna", "c1", 10_000, outcome="dropped", rejected_by="dispersed_identity",
                                         identity=35.0)])
    _, listed = _listed(run)
    assert [(e["region_class"], e["reason"], e["refused_hit"]["rejected_by"]) for e in listed] == \
        [("goi_hit_no_call", "dropped_by_rescue", "dispersed_identity")]


def test_no_record_at_all(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)])
    _, listed = _listed(run)
    assert [(e["reason"], e["model"], e["refused_hit"]) for e in listed] == [("no_record", None, None)]


def test_a_refused_hit_of_another_locus_is_not_this_region_s(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)])
    _refused(run, "locus_2_G.fna", [_hit("G.fna", "c1", 10_000, locus="locus_2")])
    _, listed = _listed(run)
    assert listed[0]["reason"] == "no_record"


def test_a_scores_table_without_the_new_columns(run):
    old = [c for c in SCORE_COLS if not c.startswith("goi_model") and not c.startswith("goi_hit")]
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)], cols=old)
    _, listed = _listed(run)
    assert len(listed) == 1
    assert listed[0]["goi_model_confidence"] == "" and listed[0]["goi_hit_best_identity"] is None


def test_the_order_is_fixed_and_the_count_matches(run):
    _scores(run, "locus_2__H.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)])
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg2", "c2", 50_000, 52_000),
                                               _region("Reg1", "c1", 9_000, 12_000)])
    report, listed = _listed(run)
    assert [(e["genome"], e["chrom"]) for e in listed] == [("G", "c1"), ("G", "c2"), ("H", "c1")]
    assert report["summary"]["goi_regions_without_call"] == len(listed) == 3
    assert gr.build_report(str(run))["regions"]["goi_regions_without_call"] == listed


# ---------------------------------------------------------------- a record the report dropped (review round 2)

def test_a_region_whose_only_model_the_report_dropped_names_it(run):
    """The genome has a counted call on c9; the AMBIGUOUS gene-order-only model inside the MEDIUM
    region on c1 is dropped by the report (drop_gene_order_only). The region still names that model,
    and the reason is model_not_counted, not no_record."""
    _scores(run, "locus_1__G.fna.scores.tsv",
            [_region("Reg1_G5_CMEDIUM_S0.42", "c1", 90_000, 120_000, model="AMBIGUOUS")])
    go_row = _goi("c1", 100_001, 101_000, "AMBIGUOUS", "syntenic_candidate_unconfirmed", 35,
                  mid="GOI_q|G_fna_b1_l1_m1").rstrip("\n") + ";Basis=gene_order_only\n"
    _gff(run, "locus_1__G.fna.gff", [_goi("c9", 5_001, 6_000, "HIGH", "confident_goi"), go_row])
    report, listed = _listed(run)
    assert [r["mrna_id"] for r in report["goi_dedup"]["gene_order_only_dropped_records"]] == ["GOI_q|G_fna_b1_l1_m1"]
    assert len(listed) == 1
    e = listed[0]
    assert e["reason"] == "model_not_counted"
    assert e["model"]["mrna_id"] == "GOI_q|G_fna_b1_l1_m1"
    assert (e["model"]["confidence"], e["model"]["report_drop"]) == ("AMBIGUOUS", "gene_order_only")


def test_a_region_whose_answer_the_report_dropped_is_a_call_not_counted(run):
    """A MEDIUM answer of the genome-wide path in a genome with a counted call elsewhere is dropped
    (genome_has_counted_call): the region near it names it, call_not_counted."""
    _scores(run, "locus_1__G.fna.scores.tsv",
            [_region("Reg1_G4_CMEDIUM_S0.40", "c1", 90_000, 120_000, model="MEDIUM")])
    answer = _goi("c1", 100_001, 101_000, "MEDIUM", "genome_wide_strong_hit", 55, mid="GOI_q|G_fna_gw1")
    _gff(run, "locus_1__G.fna.gff", [_goi("c9", 5_001, 6_000, "HIGH", "confident_goi"), answer])
    report, listed = _listed(run)
    assert [d["reason"] for d in report["goi_dedup"]["genome_wide_strong_hits_dropped_records"]] == \
        ["genome_has_counted_call"]
    assert [(e["reason"], e["model"]["mrna_id"], e["model"]["report_drop"]) for e in listed] == \
        [("call_not_counted", "GOI_q|G_fna_gw1", "genome_has_counted_call")]


def test_a_model_that_stands_has_no_report_drop(run):
    _scores(run, "locus_1__G.fna.scores.tsv", [_region("Reg1", "c1", 9_000, 12_000)])
    _gff(run, "locus_1__G.fna.gff", [_goi("c1", 10_001, 10_600, "AMBIGUOUS", "ambiguous_goi_family_member", 30)])
    _, listed = _listed(run)
    assert [(e["reason"], e["model"]["report_drop"]) for e in listed] == [("model_not_counted", "")]
