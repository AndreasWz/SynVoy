"""The copy-lineage field: with --tandem_copy_parent seed, a model found with a
seeded tandem copy (SynVoy_Parent=GOI_copy_<n>) has Identity, QueryCoverage and Bits against that
copy. HomeIdentity, HomeQueryCoverage and HomeBits give the same against the home protein; no
label reads them, and the report carries them (copy_lineage, home_*)."""
import importlib.util
import inspect
import os

import pytest

BIN = os.path.join(os.path.dirname(__file__), os.pardir, "bin")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(BIN, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


isr = _load("iterative_search_runner", "iterative_search_runner.py")
gr = _load("generate_report", "generate_report.py")
pytest.importorskip("parasail")

HOME = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACRESAPAEFRDALSR"
COPY = HOME[:60] + "W" * 10 + HOME[70:200]
MODEL = HOME[5:150].replace("L", "I")


def _cand(parent, protein=MODEL, attrs=""):
    row = "\t".join(["c1", "SynVoy", "mRNA", "1001", "1450", "80.0", "+", ".",
                     f"ID={parent}|G_fna_b1_l1_m1;SynVoy_Parent={parent};SynVoyRole=goi;Confidence=LOW;"
                     f"EvidenceType=exon_annotation;Identity=88.0;QueryCoverage=0.950;Bits=250.0{attrs}"])
    return {"gff": [row], "record": {"id": "m1", "seq": protein}, "label_subjects": [protein],
            "start": 1000, "end": 1450}


def _attrs(cand):
    return isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8])


BY_PARENT = {"GOI_q": {"seq": HOME}, "GOI_copy_1": {"seq": COPY}}


def test_a_copy_lineage_model_gets_the_home_terms():
    cand = _cand("GOI_copy_1")
    isr.annotate_model_extras([cand], BY_PARENT, None, "", home_parents=["GOI_q"])
    a = _attrs(cand)
    aln = isr.sw_align_span(HOME, MODEL)
    assert a["HomeIdentity"] == f"{aln['identity']:.1f}"
    assert a["HomeQueryCoverage"] == f"{aln['query_cov']:.3f}"
    assert a["HomeBits"] == f"{isr.sw_bits(HOME, MODEL):.1f}"
    assert "HomeQuery" not in a
    # The copy-measured terms stay as they are.
    assert (a["Identity"], a["QueryCoverage"], a["Bits"]) == ("88.0", "0.950", "250.0")


def test_with_two_home_gois_the_better_one_is_named():
    by_parent = dict(BY_PARENT, GOI_a={"seq": "MSTNPKPQRKTKRNTNRRPQDVKFPGG" * 3})
    cand = _cand("GOI_copy_1")
    isr.annotate_model_extras([cand], by_parent, None, "", home_parents=["GOI_a", "GOI_q"])
    assert _attrs(cand)["HomeQuery"] == "GOI_q"


def test_home_and_decoy_candidates_get_nothing():
    home, decoy = _cand("GOI_q"), _cand("GOI_copy_1~shuf1")
    isr.annotate_model_extras([home, decoy], BY_PARENT, None, "", home_parents=["GOI_q"])
    for cand in (home, decoy):
        assert not any(k.startswith("Home") for k in _attrs(cand))


def test_without_home_parents_nothing_is_written():
    cand = _cand("GOI_copy_1")
    before = cand["gff"][0]
    isr.annotate_model_extras([cand], BY_PARENT, None, "")
    assert cand["gff"][0] == before


def test_without_an_aligner_nothing_is_written(monkeypatch):
    def boom(*_a, **_k):
        raise RuntimeError("parasail-python is required")
    monkeypatch.setattr(isr, "first_best_subject", boom)
    cand = _cand("GOI_copy_1")
    before = cand["gff"][0]
    isr.annotate_model_extras([cand], BY_PARENT, None, "", home_parents=["GOI_q"])
    assert cand["gff"][0] == before


def test_the_judge_reads_none_of_the_terms():
    attrs = {"Eeff": "1e-30", "Place": "slot", "BlockEsyn": "0.001", "SpecMargin": "20",
             "RealCompetitor": "1", "BestOtherInNeighbourhood": "0", "QueryCoverage": "0.9"}
    with_terms = dict(attrs, HomeIdentity="40.0", HomeQueryCoverage="0.300", HomeBits="50.0")
    assert isr._judge_attrs(attrs, "Eeff") == isr._judge_attrs(with_terms, "Eeff")


def test_the_window_passes_the_home_gois():
    src = inspect.getsource(isr.process_region_block)
    assert "home_parents=home_goi_parents(args)" in src


def test_the_report_carries_the_lineage(tmp_path):
    for d in ("regions", "hits", "scores"):
        (tmp_path / d).mkdir()
    rows = [
        "\t".join(["c1", "x", "mRNA", "1001", "1450", ".", "+", ".",
                   "ID=GOI_copy_1|G_fna_b1_l1_m1;SynVoy_Parent=GOI_copy_1;SynVoyRole=goi;Confidence=HIGH;"
                   "GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=93.3;QueryCoverage=1.000;"
                   "HomeIdentity=53.0;HomeQueryCoverage=0.430;HomeBits=120.5"]),
        "\t".join(["c2", "x", "mRNA", "1001", "1450", ".", "+", ".",
                   "ID=GOI_q|G_fna_b2_l1_m1;SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=HIGH;"
                   "GOIClass=confident_goi;EvidenceType=exon_annotation;Identity=80;QueryCoverage=0.9"]),
    ]
    (tmp_path / "regions" / "locus_1__G.fna.gff").write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    report = gr.build_report(str(tmp_path))
    recs = {r["chrom"]: r for r in report["goi_dedup"]["records"]}
    assert (recs["c1"]["copy_lineage"], recs["c1"]["home_identity"], recs["c1"]["home_query_cov"],
            recs["c1"]["home_bits"]) == ("GOI_copy_1", "53.0", "0.430", "120.5")
    assert "copy_lineage" not in recs["c2"]
    assert report["summary"]["counted_copy_lineage_calls"] == 1
