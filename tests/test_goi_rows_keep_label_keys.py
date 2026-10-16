"""PR-0 review, amendment 3 (16 Oct): every final GOI row has Confidence, GOIClass and Basis.

relabel_goi_lines rewrites only the label keys a row already carries and appends only Basis,
so a GOI row written without Confidence or GOIClass would lose them in silence, and at least
nine readers compare Confidence as an exact string. Every writer of the search's GOI rows goes
through _goi_feature_attrs, which writes the provisional label (_provisional_goi_label); this
checks that for each evidence type, and that the genome-end judge leaves all three keys on the
row."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402

EVIDENCE = ("exon_annotation", "tandem_copy", "fallback_hit_span", "rescued_exon", "raw_hit")
LABELS = {"HIGH", "MEDIUM", "AMBIGUOUS", "LOW"}
# Terms for the judge: from significant and specific (HIGH) to no E at all (LOW).
JUDGED = {
    "high": "Place=slot;BlockEsyn=1e-08;EPlace=1e-30;SigBlocks=1;Eeff=1e-30;SpecMargin=50.0;RealCompetitor=0",
    "gene_order": "Place=slot;BlockEsyn=1e-08;EPlace=0.05;SigBlocks=1;Eeff=0.05;SpecMargin=50.0;RealCompetitor=0",
    "none": "Place=genome;BlockEsyn=3;EPlace=inf;SigBlocks=1;Eeff=inf",
}


def _row(evidence, extra=""):
    attrs = isr._goi_feature_attrs({"ID": f"GOI_q|g_b0_l1_{evidence}", "Name": "GOI_q",
                                    "SynVoy_Parent": "GOI_q"},
                                   evidence_type=evidence, identity=60.0, exon_count=2,
                                   query_cov=0.9, flanking_support=6)
    ftype = "gene" if evidence == "tandem_copy" else "mRNA"
    text = isr._compose_gff_attrs(attrs)
    return f"c1\tsynvoy\t{ftype}\t100\t900\t60.0\t+\t.\t{text}" + (f";{extra}" if extra else "")


@pytest.mark.parametrize("evidence", EVIDENCE)
def test_the_provisional_row_carries_the_label_keys(evidence):
    a = isr._parse_gff_attr_string(_row(evidence).split("\t")[8])
    assert a["SynVoyRole"] == "goi"
    assert a["Confidence"] == "LOW" and a["GOIClass"] and a["InferenceReason"]


@pytest.mark.parametrize("terms", sorted(JUDGED))
@pytest.mark.parametrize("evidence", EVIDENCE)
def test_every_judged_goi_row_ends_with_confidence_class_and_basis(evidence, terms):
    out, _counts = isr.relabel_goi_lines([_row(evidence, JUDGED[terms])])
    a = isr._parse_gff_attr_string(out[0].split("\t")[8])
    assert a.get("Confidence") in LABELS, out[0]
    assert a.get("GOIClass"), out[0]
    assert a.get("Basis"), out[0]
    assert a["InferenceReason"].startswith(a["Basis"])
    # Each key once: rewritten in place, never appended a second time.
    keys = [kv.split("=", 1)[0] for kv in out[0].split("\t")[8].split(";")]
    assert keys.count("Confidence") == keys.count("GOIClass") == keys.count("Basis") == 1
