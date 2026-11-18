"""What the docs say about the region step's selection (cluster_grs.py), checked against it.

* PARAMETERS.md max_regions: a GOI-anchor region is put in only when no cluster lies near a GOI
  model labelled HIGH, MEDIUM or AMBIGUOUS, as LOW rows and answers of the genome-wide path are no
  GOI intervals (load_goi_intervals_from_gff). So a counted model 25 kb past its neighbourhood gets
  a region of its own when the only other GOI model, the one inside a cluster, is LOW or an answer,
  and none when that model is AMBIGUOUS. The text said "near any GOI model".
* ALGORITHM.md section 14 said that the selection stops at the first region below both floors;
  it adds every region that passes a floor, in its order, until the cap (section 8).
All synthetic."""
import contextlib
import csv
import io
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "bin"))

import cluster_grs as cg  # noqa: E402


def _hit(q, chrom, s, e, strand="+", ident=80.0):
    ts, te = (s, e) if strand == "+" else (e, s)
    return f"{q}\t{chrom}\t{ident}\t300\t0\t0\t1\t100\t{ts}\t{te}\t1e-30\t200\n"


def _model(chrom, start, end, confidence, cls, mid):
    return "\t".join([chrom, "SynVoy", "mRNA", str(start), str(end), ".", "+", ".",
                      f"ID=GOI_X|{mid};SynVoy_Parent=GOI_X;SynVoyRole=goi;Confidence={confidence};"
                      f"GOIClass={cls};EvidenceType=exon_annotation;BlockFlankingSupport=5"])


def _regions(tmp, other_label, other_class="x"):
    """chrA: five neighbours in home order and a HIGH model 25 kb past the last of them; chrB: a
    GOI hit beside two neighbours, holding a model labelled ``other_label`` of class
    ``other_class``."""
    with open(tmp / "home.bed", "w") as fh:
        for i in range(1, 11):
            fh.write(f"chrH\t{i * 10000}\t{i * 10000 + 2000}\tF{i}\t0\t{'+' if i % 2 else '-'}\n")
    rows = [_hit(f"F{i}", "chrA", 100000 + i * 10000, 100000 + i * 10000 + 600, "+" if i % 2 else "-")
            for i in range(1, 6)]
    rows += [_hit("GOI_X", "chrB", 500000, 500600, ident=40.0), _hit("F9", "chrB", 510000, 510600),
             _hit("F10", "chrB", 520000, 520600, "-")]
    (tmp / "g.m8").write_text("".join(rows))
    (tmp / "t.gff").write_text("\n".join([
        "##gff-version 3", _model("chrA", 175601, 176600, "HIGH", "confident_goi", "g_b1_l1m1"),
        _model("chrB", 500001, 500600, other_label, other_class, "g_b2_l1m1")]) + "\n")
    argv = sys.argv
    sys.argv = ["cluster_grs.py", "--hits", str(tmp / "g.m8"), "--synteny_bed", str(tmp / "home.bed"),
                "--target_gff", str(tmp / "t.gff"), "--output", str(tmp / "r.bed"),
                "--scores_output", str(tmp / "s.tsv"), "--flanking_count", "10",
                "--cluster_distance", "50000"]
    try:
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            cg.main()
    finally:
        sys.argv = argv
    with open(tmp / "s.tsv", newline="") as fh:
        return [(r["chrom"], r["region_class"]) for r in csv.DictReader(fh, delimiter="\t")]


def test_a_low_model_in_a_cluster_does_not_stop_the_goi_anchor_region(tmp_path):
    (tmp_path / "low").mkdir()
    (tmp_path / "ambiguous").mkdir()
    low = _regions(tmp_path / "low", "LOW")
    assert ("chrA", "goi_anchor") in low and not any(chrom == "chrB" for chrom, _ in low)
    ambiguous = _regions(tmp_path / "ambiguous", "AMBIGUOUS")
    assert ("chrB", "goi_overlap") in ambiguous and ("chrA", "goi_anchor") not in ambiguous


def test_an_answer_in_a_cluster_does_not_stop_the_goi_anchor_region(tmp_path):
    """An answer of the genome-wide path (MEDIUM, no synteny) is no GOI interval either: the
    regions are those of the LOW case, as in a run without the path."""
    (tmp_path / "answer").mkdir()
    (tmp_path / "low").mkdir()
    answer = _regions(tmp_path / "answer", "MEDIUM", cg.synvoy_label.GENOME_WIDE_CLASS)
    assert ("chrA", "goi_anchor") in answer and not any(chrom == "chrB" for chrom, _ in answer)
    assert answer == _regions(tmp_path / "low", "LOW")


def test_the_max_regions_text_names_the_labels_that_count():
    text = open(os.path.join(ROOT, "docs", "PARAMETERS.md"), encoding="utf-8").read()
    section = " ".join(text[text.index("### `max_regions`"):].split("\n### ", 1)[0].split())
    rule = ("near any GOI model labelled HIGH, MEDIUM or AMBIGUOUS; LOW models and answers of the "
            "genome-wide path take no part in this rule")
    assert rule in section
    assert "lies near another GOI model labelled HIGH, MEDIUM or AMBIGUOUS" in section
    assert not re.search(r"near (?:any|another) GOI model\)", section)


def test_the_known_limits_describe_the_selection_as_it_is():
    text = open(os.path.join(ROOT, "docs", "ALGORITHM.md"), encoding="utf-8").read()
    limits = " ".join(text[text.index("## 14. Properties and known limits"):].split())
    assert "stops at the first region below both floors" not in " ".join(text.split())
    assert ("Region selection adds every region that passes a floor, in the order of §8, until "
            "the cap.") in limits
    section8 = " ".join(text[text.index("## 8. Region scores"):text.index("## 9.")].split())
    assert "is added, in this order, until 6 regions are written" in section8
