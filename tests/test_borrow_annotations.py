"""Gene models transferred from an annotated target when the home genome has no annotation
(batch of 2026-10-15, HG-5b; bin/borrow_annotations.py, lib/HomeLoci.groovy).

Before, the borrowed genes were gene rows without CDS (the walk read them as non-coding: drawn,
never anchors), placed by the first row of the hit BED, borrowed from the best-annotated target
rather than the closest one, from a target region of min..max of all hits. Now: the primary
locus's best interval, the closest annotated target in the phylogenetic order, the hit cluster
with the most home proteins, miniprot models with CDS rows, a reciprocal check, one model per
place, one merged GFF and a home proteome built from it without the genes over the query's own
model. All synthetic (no MMseqs2 / miniprot run)."""
import json
import os
import random
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import borrow_annotations as ba  # noqa: E402
from sequence_utils import cds_models_by_gene  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def test_the_home_locus_is_the_best_interval_not_the_first_row(tmp_path):
    bed = tmp_path / "locus_1.bed"
    # First row: a weak hit on another scaffold (the alphabetically first one).
    bed.write_text("chrA\t946896\t947487\tgene_loc\t8e-08\t-\t54.7\n"
                   "chrB\t14616338\t14617241\tgene_loc\t2e-134\t+\t271\n"
                   "chrB\t9426937\t9427759\tgene_loc\t4e-33\t+\t129\n")
    assert ba.goi_position(str(bed)) == ("chrB", 14616338, 14617241)
    bed.write_text("chrB\t9426937\t9427759\tbest_interval\t4e-33\t+\t129\t9426937\t9427750\n"
                   "chrB\t14616338\t14617241\tgene_loc\t2e-134\t+\t271\n")
    assert ba.goi_position(str(bed)) == ("chrB", 9426937, 9427759)    # split_loci's choice
    bed.write_text("")
    assert ba.goi_position(str(bed)) is None


def test_the_target_is_the_closest_annotated_one(tmp_path):
    targets = [("/g/A.fna", "/g/A.gff", 900), ("/g/B.fna", "/g/B.gff", 50), ("/g/C.fna", "/g/C.gff", 70)]
    fna, gff, why = ba.choose_target(targets, {"X.fna": 1, "B.fna": 2, "A.fna": 5, "C.fna": 7})
    assert fna == "/g/B.fna" and "closest annotated target (distance 2)" in why   # not the best annotated
    ident = {"/g/A.fna": 61.0, "/g/B.fna": 55.0, "/g/C.fna": 74.5}.get
    fna, _g, why = ba.choose_target(targets, {"A.fna": 2, "B.fna": 2, "C.fna": 5}, ident)
    assert fna == "/g/A.fna" and "2 at that distance" in why                      # a tie: identity
    fna, _g, why = ba.choose_target(targets, {}, ident)
    assert fna == "/g/C.fna" and "74.5 %" in why                                  # no distances
    assert ba.choose_target(targets, {})[0] == "/g/A.fna"
    assert ba.choose_target([], {"A.fna": 1}) is None
    (tmp_path / "sorted.txt").write_text("B.fna\t2\nA.fna\t5\n")
    assert ba.read_sorted_genomes(str(tmp_path / "sorted.txt")) == {"B.fna": 2.0, "A.fna": 5.0}
    (tmp_path / "alpha.txt").write_text("A.fna\t0\nB.fna\t0\n")                  # no order at all
    assert ba.read_sorted_genomes(str(tmp_path / "alpha.txt")) == {}


def test_annotated_targets_need_ten_cds_rows(tmp_path):
    (tmp_path / "A.fna").write_text(">c\nACGT\n")
    (tmp_path / "A.gff").write_text("".join(f"c\tx\tCDS\t{i}\t{i + 2}\t.\t+\t0\tParent=m\n" for i in range(1, 31, 3)))
    (tmp_path / "B.fna").write_text(">c\nACGT\n")
    (tmp_path / "B.gff").write_text("c\tx\tCDS\t1\t3\t.\t+\t0\tParent=m\n")
    (tmp_path / "C.fna").write_text(">c\nACGT\n")
    assert [os.path.basename(t[0]) for t in ba.annotated_targets(str(tmp_path))] == ["A.fna"]


def _hit(q, chrom, start, end):
    return {"query": q, "chrom": chrom, "start": start, "end": end, "pident": 70.0,
            "evalue": 1e-20, "bits": 100.0}


def test_the_target_region_is_the_cluster_with_most_home_proteins():
    hits = [_hit("p1", "c1", 1000, 2000), _hit("p2", "c1", 30000, 31000), _hit("p3", "c1", 60000, 61000),
            _hit("p1", "c2", 5000, 6000), _hit("p1", "c2", 7000, 8000), _hit("p4", "c2", 9000, 9500),
            _hit("p5", "c1", 200000, 201000)]                     # beyond the 50 kb gap
    assert ba.target_region(hits, 100000) == ("c1", 0, 161000, 3)
    assert ba.target_region([], 100000) is None


def _model(source, start, end, identity, chrom="h"):
    return {"source": source, "name": source, "chrom": chrom, "strand": "+", "identity": identity,
            "coverage": 0.9, "protein": "M", "exons": [(start, start + 100, 0), (end - 100, end, 0)]}


def test_one_model_per_place():
    borrowed = [_model("t1", 1000, 3000, 0.62), _model("t2", 1200, 3100, 0.81), _model("t3", 8000, 9000, 0.7)]
    predicted = {"p1": {"chrom": "h", "strand": "+", "exons": [(1101, 2900, 0)]},
                 "p2": {"chrom": "h", "strand": "+", "exons": [(5001, 6000, 0)]},
                 "p3": {"chrom": "h", "strand": "+", "exons": [(8901, 9500, 0)]}}
    kept, dropped = ba.one_model_per_place(borrowed, predicted)
    assert [m["source"] for m in kept] == ["t2", "t3"]            # t2 beats t1 at the same place
    assert dropped == {"p1"}                                        # p3 overlaps t3 by less than half


def test_the_reciprocal_check(tmp_path):
    rng = random.Random(3)
    p1 = "M" + "".join(rng.choice(AA) for _ in range(150))
    p2 = "M" + "".join(rng.choice(AA) for _ in range(150))
    model = "".join(rng.choice(AA) if rng.random() < 0.2 else c for c in p1)
    proteins = [("t1", "T1", p1), ("t2", "T2", p2)]
    good = dict(_model("t1", 0, 500, 0.8), protein=model)
    wrong = dict(_model("t2", 0, 500, 0.8), protein=model)        # its best match is t1
    assert [m["source"] for m in ba.reciprocal([good, wrong], proteins)] == ["t1"]


def _gff_gene(gid, chrom, start, end):
    return [f"{chrom}\tp\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}",
            f"{chrom}\tp\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.t1;Parent={gid}",
            f"{chrom}\tp\tCDS\t{start}\t{end}\t.\t+\t0\tID={gid}.cds;Parent={gid}.t1"]


def test_the_merged_gff_and_the_home_proteome(tmp_path):
    rng = random.Random(4)
    seq = "".join(rng.choice("ACG") for _ in range(12000))    # no T: no stop codon anywhere
    genome = {"h": seq}
    pred = tmp_path / "pred.gff"
    pred.write_text("\n".join(["##gff-version 3"] + _gff_gene("p1", "h", 1001, 1900)
                              + _gff_gene("p2", "h", 5001, 5900) + _gff_gene("p3", "h", 9001, 9900)) + "\n")
    m = _model("gene-T7", 990, 1950, 0.75)
    m["exons"] = [(990, 1200, 0), (1500, 1950, 2)]
    out = tmp_path / "merged.gff"
    ba.write_merged_gff(str(out), str(pred), {"p1"}, [m], "Drosophila_simulans")
    text = out.read_text()
    assert "ID=p1" not in text and "ID=p2" in text and "ID=p3" in text
    assert "ID=borrowed_gene-T7;Name=gene-T7;gene_biotype=protein_coding;Note=transferred from Drosophila_simulans" in text
    models = cds_models_by_gene(str(out))
    assert sorted(models) == ["borrowed_gene-T7", "p2", "p3"]
    assert models["borrowed_gene-T7"]["exons"] == [(991, 1200, 0), (1501, 1950, 2)]   # 1-based, phases
    record = {"source": "predicted+borrowed", "partial": True}
    n = ba.write_home_proteome(str(out), genome, [("h", 9050, 9100)], str(tmp_path / "db"), record)
    assert n == 2                                                      # p3 lies over the GOI model
    faa = (tmp_path / "db" / "home_proteome.faa").read_text()
    assert ">borrowed_gene-T7" in faa and ">p2" in faa and ">p3" not in faa
    info = json.loads((tmp_path / "db" / "home_annotation.json").read_text())
    assert info["partial"] is True and info["left_out_over_goi"] == ["p3"] and info["proteins"] == 2
    assert (tmp_path / "db" / "home_genes.tsv").read_text().splitlines()[0] == "gene_id\tchrom\tstart\tend\tstrand"


def test_the_window_grows_until_it_holds_enough_predictions():
    preds = {f"p{i}": {"chrom": "h", "strand": "+", "exons": [(i * 20000 + 1, i * 20000 + 900, 0)]}
             for i in range(1, 40)}
    lo, hi = ba.home_window(("h", 400000, 401000), preds, 2_000_000, 50000, 10)
    inside = sum(1 for m in preds.values() if ba.model_span(m)[0] >= lo and ba.model_span(m)[1] <= hi)
    assert inside >= 20 and hi - lo <= ba.WINDOW_MAX_BP + 1000
    lo2, hi2 = ba.home_window(("h", 400000, 401000), {}, 2_000_000, 50000, 10)
    assert hi2 - lo2 >= ba.WINDOW_MAX_BP                              # nothing to find: the cap


def test_the_report_names_a_partial_home_proteome(tmp_path):
    import generate_report as gr
    assert gr.load_home_annotation(str(tmp_path)) is None              # a home GFF: no file
    (tmp_path / "home_annotation.json").write_text(json.dumps(
        {"source": "predicted+borrowed", "partial": True, "target": "Drosophila_simulans",
         "borrowed": 18, "predicted": 40, "predictions_dropped": 12, "proteins": 58,
         "left_out_over_goi": ["p3"]}))
    out = gr.load_home_annotation(str(tmp_path))
    assert out["borrowed"] == 18 and out["genes_left_out_over_goi"] == 1
    assert "proteome is partial" in out["warning"] and "transferred from Drosophila_simulans" in out["warning"]


NEXTFLOW = shutil.which("nextflow")


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_the_primary_locus_bed_is_taken_by_its_index(tmp_path):
    script = tmp_path / "t.nf"
    script.write_text(
        "workflow {\n"
        "  println 'OLD ' + HomeLoci.primary([file('/x/locus_2.bed'), file('/x/locus_10.bed'), file('/x/locus_1.bed')]).getFileName()\n"
        "  println 'NEW ' + HomeLoci.primary([file('/x/home_locus_3.bed'), file('/x/home_locus_2.bed')]).getFileName()\n"
        "  println 'ONE ' + HomeLoci.primary(file('/x/locus_1.bed')).getFileName()\n"
        "}\n")
    p = subprocess.run([NEXTFLOW, "-q", "run", str(script), "-lib", os.path.join(ROOT, "lib")],
                       capture_output=True, text=True, cwd=tmp_path, timeout=300)
    assert p.returncode == 0, p.stdout + p.stderr
    assert "OLD locus_1.bed" in p.stdout and "NEW home_locus_2.bed" in p.stdout and "ONE locus_1.bed" in p.stdout
