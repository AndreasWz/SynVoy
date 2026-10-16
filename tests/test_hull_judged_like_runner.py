"""Review fix RF-1d (batch of 2026-10-15): the hull rescue is judged like the runner.

* the label protein goes through synvoy_label.frame_guard (FrameFixed / FrameBroken; a
  broken frame is judged on the translation's pieces between stops);
* completeness is the query span of the miniprot alignment (Target=), not CDS aa / query
  length, which reaches 1.0 for any long model;
* the specificity gets the GOI's home neighbourhood (it got None);
* the 40 % identity gate stays as a cost filter in front of the judge (decision, 15 Oct).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import rescue_goi_hull as hull  # noqa: E402

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'}


def _dna(p):
    return "".join(CODON[a] for a in p)


def _rows(spans, strand="+"):
    return [["w", "miniprot", "CDS", str(s), str(e), ".", strand, "0", "Parent=m"] for s, e in spans]


def test_coverage_is_the_query_span_of_the_alignment():
    mrna = ["w", "miniprot", "mRNA", "1", "3000", "100", "+", ".", "ID=m;Identity=0.9;Target=q 11 59"]
    rows = _rows([(1, 1500), (2001, 3000)])                   # 833 aa of CDS
    assert hull._query_span_coverage(mrna, rows, 100) == 0.49
    assert hull._cds_coverage(rows, 100) == 1.0              # the old measure
    no_target = mrna[:8] + ["ID=m;Identity=0.9"]
    assert hull._query_span_coverage(no_target, rows, 100) == 1.0


def test_the_label_protein_goes_through_the_frame_guard():
    good = "AAA" + _dna("MKTAYIAKQRQISF") + "TAA"
    assert hull.hull_label_protein(good, _rows([(4, 3 + 42)])) == ("MKTAYIAKQRQISF", "")
    cds = "ATGTAGTTGAATTACTCTAAATGA"             # a stop in every frame (table 1)
    protein, flag = hull.hull_label_protein("AAA" + cds + "GGG", _rows([(4, 3 + len(cds))]))
    assert flag == "FrameBroken=1"
    assert protein == hull.translate(cds).rstrip("*") and "*" in protein   # judged in pieces


def test_the_flanking_anchors_of_the_genome_give_the_neighbourhood(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n"
                   "c1\tx\tmRNA\t100\t900\t.\t+\t.\tID=a|g_b0;SynVoy_Parent=gene-A;SynVoyRole=flanking;Confidence=HIGH\n"
                   "c1\tx\tmRNA\t2000\t2900\t.\t+\t.\tID=b|g_b0;SynVoy_Parent=gene-C;SynVoyRole=flanking;Confidence=LOW\n"
                   "c1\tx\tmRNA\t1000\t1900\t.\t+\t.\tID=GOI_q|g_b0;SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=HIGH\n")
    assert hull.flanking_parent_ids(str(gff)) == ["gene-A", "gene-C"]


def test_the_judge_gets_the_neighbourhood_and_the_span_coverage(monkeypatch):
    seen = {}
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", lambda q, s, sp: (1e-30, 200.0))
    monkeypatch.setattr(hull, "sw_bits", lambda q, x: 200.0)

    def fake_spec(proteome, qid, q, x, bits, hood, n_flanking=0):
        seen["hood"], seen["n_flanking"] = hood, n_flanking
        return {"margin": 50.0, "real_competitor": False, "best_other": "", "best_other_in_neighbourhood": None}
    monkeypatch.setattr(hull.synvoy_specificity, "specificity", fake_spec)
    monkeypatch.setattr(hull.synvoy_specificity, "gff_attrs", lambda spec: "SpecMargin=50")
    rows = _rows([(1, 300)])
    conf, _key, attrs, _bits = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, 10 ** 7, object(),
                                                     label_protein="M" * 100, coverage=0.49,
                                                     neighbourhood={"gene-A", "gene-B"}, n_flanking=10)
    assert seen["hood"] == {"gene-A", "gene-B"} and seen["n_flanking"] == 10   # RF-2
    assert conf == "MEDIUM"                                  # significant and specific, 0.49 < 0.7
    conf, *_ = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, 10 ** 7, object(),
                                     label_protein="M" * 100, coverage=0.9, neighbourhood=None)
    assert conf == "HIGH"


def test_the_hull_process_hands_over_the_home_database_folder():
    root = os.path.join(os.path.dirname(__file__), "..")
    module = open(os.path.join(root, "modules", "rescue_goi_hull.nf")).read()
    assert "--home_db_dir ${home_db}" in module
    sub = open(os.path.join(root, "subworkflows", "adjudicate_and_report.nf")).read()
    assert "RESCUE_GOI_HULL(paralog_inputs_ch, genomes_dir, rescue_query_ch, home_db)" in sub
    assert "home_db.map { d -> d.resolve('home_proteome.faa') }" in sub
