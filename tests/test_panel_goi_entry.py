"""The GOI entry of the ownership panel (fix-phase item 17, 2026-10-28).

The panel's GOI entry used to be the panel gene the query aligned to best, with no minimum.
Without a home GFF the GOI's own gene is never in the home proteome
(borrow_annotations.write_home_proteome leaves the genes over the query's home model out), so a
neighbour became the "GOI", another neighbour owned the true orthologs, and the report relabelled
them paralog_not_goi (stage 5c's Candida run without its GFF: 19 HIGH or MEDIUM calls). Now the GOI entry is q's own home gene,
found as the search's specificity check finds it (synvoy_specificity: the gene the home step took,
from goi_info.json, else the home protein >= 95 % identical over >= 90 % of q), and q itself, at the
primary home locus, when no home gene is q's own. q is the search's q, the GOI's home protein (the
full-length GOI record of goi_exons.faa), with the query file as the fallback: with the query, a
query of another species found no own gene where the search found one (review of 28 Oct). How the
entry was found is in the log, the meta.tsv and the report's locus_ownership_summary.
"""
from __future__ import annotations

import csv
import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import build_home_paralog_panel as bhp  # noqa: E402
import generate_report as gr  # noqa: E402

pytest.importorskip("parasail")

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(seed: int, n: int = 160) -> str:
    rng = random.Random(seed)
    return "M" + "".join(rng.choice(_AA) for _ in range(n - 1))


def _mutate(seq: str, share: float, seed: int) -> str:
    """``seq`` with about ``share`` of its residues replaced (a homolog at 1 - share identity)."""
    rng = random.Random(seed)
    return "".join(rng.choice(_AA) if i and rng.random() < share else c for i, c in enumerate(seq))


Q = _protein(1)
NEIGHBOUR_PARALOG = _mutate(Q, 0.55, 2)      # a relative of q next to it (~45 % identity)
NEIGHBOUR_OTHER = _protein(3)                # an unrelated neighbour
OTHER_COPY = Q                               # a protein identical to q elsewhere in the genome


def _write(path: Path, text: str) -> str:
    path.write_text(text)
    return str(path)


def _fasta(path: Path, records) -> str:
    return _write(path, "".join(f">{rid}\n{seq}\n" for rid, seq in records))


def _gff(path: Path, genes) -> str:
    """genes: (chrom, start, end, gene id)."""
    return _write(path, "".join(f"{c}\tsrc\tgene\t{s}\t{e}\t.\t+\t.\tID={g};Name={g}\n"
                                for c, s, e, g in genes))


def _goi_info(path: Path, gene_id=None, found_by="miniprot", ambiguous=False) -> str:
    return _write(path, json.dumps({
        "query_id": "q1", "found_by": found_by, "gene_id": gene_id,
        "home_gene_ambiguous": ambiguous, "runner_up_gene_id": "gene-R" if ambiguous else None,
        "runner_up_score_fraction": 0.97 if ambiguous else None}))


def _run_panel(tmp: Path, *, genes, proteins, beds, query=(("q1", Q),), goi_info=None, home_goi=None,
               extra=()):
    gff = _gff(tmp / "home.gff", genes)
    prot = _fasta(tmp / "home_proteome.faa", proteins)
    bed_paths = []
    for name, rows in beds:
        bed_paths.append(_write(tmp / name, "".join(f"{c}\t{s}\t{e}\tgene_loc\t1e-30\t+\n"
                                                     for c, s, e in rows)))
    qfa = _fasta(tmp / "query.faa", query)
    out = tmp / "panel.faa"
    cmd = [sys.executable, str(BIN / "build_home_paralog_panel.py"), "--home_gff", gff,
           "--home_proteome", prot, "--locus_beds", *bed_paths, "--query", qfa,
           "--output", str(out), *extra]
    if goi_info:
        cmd += ["--goi_info", goi_info]
    if home_goi:
        cmd += ["--home_goi", _fasta(tmp / "goi_exons.faa", home_goi)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(str(out) + ".meta.tsv") as fh:
        meta = {row["panel_id"]: row for row in csv.DictReader(fh, delimiter="\t")}
    labels = [ln[1:].split()[0] for ln in open(out) if ln.startswith(">")]
    return meta, labels, r.stderr


def _goi(meta):
    return {pid: row["goi_found_by"] for pid, row in meta.items() if row["is_goi_gene"] == "1"}


def test_the_home_steps_gene_is_the_goi_entry(tmp_path):
    """goi_info.json names q's own gene (found by name); it is the GOI entry although a copy
    elsewhere is identical to q and the named gene's protein is a partial model."""
    partial = Q[:80]
    meta, labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-Q"), ("c1", 2500, 3500, "gene-N"), ("c9", 50000, 51000, "gene-C")],
        proteins=[("gene-Q", partial), ("gene-N", NEIGHBOUR_OTHER), ("gene-C", OTHER_COPY)],
        beds=[("locus_1.bed", [("c1", 1100, 1900)])],
        goi_info=_goi_info(tmp_path / "goi_info.json", gene_id="gene-Q", found_by="name"))
    assert _goi(meta) == {"locus_1|gene-Q": "goi_info"}
    assert meta["home_gene-C|gene-C"]["is_goi_gene"] == "0"
    assert "GOI entry: gene-Q, q's own home gene, the one the home step took" in err


def test_the_identity_rule_without_a_gene_from_the_home_step(tmp_path):
    """No gene in goi_info.json (a miniprot model): the home protein near-identical to q is q's own."""
    meta, _labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-Q"), ("c1", 2500, 3500, "gene-N")],
        proteins=[("gene-Q", Q), ("gene-N", NEIGHBOUR_PARALOG)],
        beds=[("locus_1.bed", [("c1", 1100, 3400)])],
        goi_info=_goi_info(tmp_path / "goi_info.json"))
    assert _goi(meta) == {"locus_1|gene-Q": "identity_rule"}
    assert meta["locus_1|gene-N"]["is_goi_gene"] == "0"
    assert "by the identity rule: >= 95 % identity over >= 90 % of q" in err


def test_an_ambiguous_choice_of_the_home_step_falls_to_the_identity_rule(tmp_path):
    """goi_info.json flags its gene as ambiguous (a relative may have been taken): as in the search,
    the gene is not used; the identity rule names the near-identical protein instead."""
    meta, _labels, _err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-R"), ("c1", 2500, 3500, "gene-Q")],
        proteins=[("gene-R", NEIGHBOUR_PARALOG), ("gene-Q", Q)],
        beds=[("locus_1.bed", [("c1", 1100, 3400)])],
        goi_info=_goi_info(tmp_path / "goi_info.json", gene_id="gene-R", found_by="name",
                           ambiguous=True))
    assert _goi(meta) == {"locus_1|gene-Q": "identity_rule"}
    assert meta["locus_1|gene-R"]["is_goi_gene"] == "0"


def test_an_ambiguous_choice_and_no_near_identical_protein_gives_q(tmp_path):
    meta, _labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-R")],
        proteins=[("gene-R", NEIGHBOUR_PARALOG)],
        beds=[("locus_1.bed", [("c1", 900, 2100)])],
        goi_info=_goi_info(tmp_path / "goi_info.json", gene_id="gene-R", found_by="position",
                           ambiguous=True))
    assert _goi(meta) == {"locus_1|q1": "query"}
    # the log says why goi_info.json gives no gene (it named one, set aside as ambiguous)
    assert "goi_info.json names gene-R, which the home step flagged as ambiguous" in err


def test_a_query_of_another_species_finds_the_searchs_own_gene(tmp_path):
    """The review's case: a home GFF, a query of another species (~75 % identical to its home
    gene R), and the home step's choice of R flagged ambiguous. The search runs the identity rule
    on the GOI's home protein (R's annotated protein, the GOI record of goi_exons.faa) and names R
    its own gene; with the query the panel named none, kept R as a paralog entry, and the report
    would relabel R's orthologs paralog_not_goi. With --home_goi R is the GOI entry."""
    R = _protein(21, 300)
    query = _mutate(R, 0.25, 22)                      # the query of another species
    paralog = _mutate(R, 0.45, 23)
    case = dict(genes=[("c1", 1000, 2000, "gene-R"), ("c1", 2500, 3500, "gene-P")],
                proteins=[("gene-R", R), ("gene-P", paralog)],
                beds=[("locus_1.bed", [("c1", 1100, 3400)])], query=(("P99999", query),))
    info = _goi_info(tmp_path / "goi_info.json", gene_id="gene-R", found_by="position", ambiguous=True)
    with_home, _l, err = _run_panel(tmp_path, goi_info=info, home_goi=(("GOI_P99999", R),
                                    ("GOI_P99999|exon_1", R[:100])), **case)
    assert _goi(with_home) == {"locus_1|gene-R": "identity_rule"}
    assert "q for the GOI entry: GOI_P99999 (300 aa, the GOI's home protein from goi_exons.faa)" in err
    query_only, _l, _err = _run_panel(tmp_path, goi_info=info, **case)
    assert _goi(query_only) == {"locus_1|P99999": "query"}         # what the query alone gave
    assert query_only["locus_1|gene-R"]["is_goi_gene"] == "0"


def test_a_query_id_with_a_pipe_keeps_the_home_protein_as_q(tmp_path):
    """A query id that keeps a '|' (NCBI's lcl|ORF1): goi_exons.faa holds GOI_lcl|ORF1 and its
    exons, and the search takes GOI_lcl|ORF1 as its home protein; so does the panel (it dropped
    every id with a '|' and fell back to the query). With the home step's choice refused (a query
    of another species), the identity rule on the home protein names R, as in the search."""
    R = _protein(21, 300)
    meta, _labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-R"), ("c1", 2500, 3500, "gene-P")],
        proteins=[("gene-R", R), ("gene-P", _mutate(R, 0.45, 23))],
        beds=[("locus_1.bed", [("c1", 1100, 3400)])], query=(("lcl|ORF1", _mutate(R, 0.45, 22)),),
        goi_info=_goi_info(tmp_path / "goi_info.json"),
        home_goi=(("GOI_lcl|ORF1", R), ("GOI_lcl|ORF1|exon_1", R[:100])))
    assert _goi(meta) == {"locus_1|gene-R": "identity_rule"}
    assert "q for the GOI entry: GOI_lcl|ORF1 (300 aa, the GOI's home protein from goi_exons.faa)" in err


def test_without_an_own_gene_the_goi_entry_is_the_home_protein(tmp_path):
    """No home gene is the GOI's own (a home without a GFF): the entry is the GOI's home protein,
    named by its record id, at the primary locus."""
    model = _mutate(Q, 0.01, 31)                      # the home model, one or two residues off
    meta, labels, _err = _run_panel(
        tmp_path,
        genes=[("c1", 400, 900, "gene-N1")],
        proteins=[("gene-N1", NEIGHBOUR_PARALOG)],
        beds=[("locus_1.bed", [("c1", 500, 800)])],
        goi_info=_goi_info(tmp_path / "goi_info.json"),
        home_goi=(("GOI_q1|exon_1", model[:60]), ("GOI_q1", model), ("GOI_copy_q1_2", model + "W")))
    assert _goi(meta) == {"locus_1|GOI_q1": "query"}
    faa = (tmp_path / "panel.faa").read_text().split(">locus_1|GOI_q1\n", 1)[1].split(">", 1)[0]
    assert faa.replace("\n", "") == model


def test_without_qs_gene_a_neighbour_no_longer_becomes_the_goi(tmp_path):
    """Today's failure (stage 5c, Candida without its GFF): q's own gene (gene-G) is in the merged
    GFF but not in the proteome; the neighbour that aligns best (a relative of q) was the "GOI".
    Now q itself is the GOI entry at the locus, and both neighbours are paralog entries. The log
    names the genes kept and the one without a protein (it printed gids[:kept])."""
    meta, labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-G"), ("c1", 400, 900, "gene-N1"), ("c1", 2100, 2600, "gene-N2")],
        proteins=[("gene-N1", NEIGHBOUR_PARALOG), ("gene-N2", NEIGHBOUR_OTHER)],
        beds=[("locus_1.bed", [("c1", 1100, 1900)])],
        goi_info=_goi_info(tmp_path / "goi_info.json"),
        extra=["--pad", "1000"])
    assert _goi(meta) == {"locus_1|q1": "query"}
    assert meta["locus_1|gene-N1"]["is_goi_gene"] == "0"
    assert meta["locus_1|gene-N2"]["is_goi_gene"] == "0"
    q_row = meta["locus_1|q1"]
    assert (q_row["locus_id"], q_row["gene"], q_row["chrom"], q_row["start"], q_row["end"]) == \
        ("locus_1", "q1", "c1", "1100", "1900")
    kept = [ln for ln in err.splitlines() if "[locus_1]" in ln][0]
    assert "-> 2 panel gene(s):" in kept
    assert "without a protein in the home proteome: ['gene-G']" in kept
    assert "gene-G" not in kept.split("without")[0]
    assert "GOI entry: q itself (q1), as locus_1|q1; no home gene is q's own" in err
    assert "q for the GOI entry: q1 (160 aa, the query; no home GOI protein given)" in err


def test_the_query_goes_to_the_primary_locus_with_one_pipe(tmp_path):
    """The primary home locus is the lowest index (lib/HomeLoci.groovy: locus_2 before locus_10),
    and a query id holding '|' keeps exactly one '|' in the panel id."""
    meta, labels, _err = _run_panel(
        tmp_path,
        genes=[("c1", 400, 900, "gene-N1"), ("c2", 400, 900, "gene-N2")],
        proteins=[("gene-N1", NEIGHBOUR_PARALOG), ("gene-N2", NEIGHBOUR_OTHER)],
        beds=[("locus_10.bed", [("c1", 500, 800)]), ("locus_2.bed", [("c2", 500, 800)])],
        query=(("myq|frag 1", Q),))
    goi = _goi(meta)
    assert len(goi) == 1
    pid = next(iter(goi))
    assert pid.count("|") == 1 and pid.startswith("locus_2|")
    assert meta[pid]["locus_id"] == "locus_2" and meta[pid]["chrom"] == "c2"
    assert gr._panel_locus(pid) == "locus_2"
    assert pid in labels


@pytest.mark.parametrize("family_cap", [None, "0"])
def test_qs_own_gene_outside_every_locus_span_is_at_its_pseudo_locus(tmp_path, family_cap):
    """q's own gene from the home step, in no locus span: it is the GOI entry at the family
    augmentation's pseudo-locus home_<gene>, as the parent had it, whether the family members
    hold it (the default cap) or not (cap 0); no search started from there, so the report marks a
    call it owns as found from another home locus."""
    meta, _labels, err = _run_panel(
        tmp_path,
        genes=[("c1", 400, 900, "gene-N1"), ("c7", 90000, 91000, "gene-Q")],
        proteins=[("gene-N1", NEIGHBOUR_PARALOG), ("gene-Q", Q)],
        beds=[("locus_1.bed", [("c1", 500, 800)])],
        goi_info=_goi_info(tmp_path / "goi_info.json", gene_id="gene-Q", found_by="position"),
        extra=["--max_family_paralogs", family_cap] if family_cap else [])
    assert _goi(meta) == {"home_gene-Q|gene-Q": "goi_info"}
    row = meta["home_gene-Q|gene-Q"]
    assert (row["locus_id"], row["chrom"], row["start"]) == ("home_gene-Q", "c7", "90000")
    assert ("GOI entry: gene-Q, q's own home gene, the one the home step took (goi_info.json); "
            "in no locus span, at home_gene-Q.") in err


def test_qs_own_gene_left_out_by_the_gene_cap_joins_its_locus(tmp_path):
    """Five longer genes fill locus_1's gene cap (5) and leave q's own gene out, although the
    locus span holds it; the family augmentation would enter it at the pseudo-locus home_gene-Q.
    It is the GOI entry at locus_1 instead, the locus it lies in."""
    genes = [("c1", 1000 + 1500 * i, 2000 + 1500 * i, f"gene-L{i}") for i in range(5)]
    meta, _labels, err = _run_panel(
        tmp_path,
        genes=genes + [("c1", 4000, 4300, "gene-Q")],
        proteins=[(g, _protein(40 + i)) for i, (_c, _s, _e, g) in enumerate(genes)] + [("gene-Q", Q)],
        beds=[("locus_1.bed", [("c1", 1000, 8000)])],
        goi_info=_goi_info(tmp_path / "goi_info.json", gene_id="gene-Q", found_by="name"))
    assert _goi(meta) == {"locus_1|gene-Q": "goi_info"}
    assert "home_gene-Q|gene-Q" not in meta
    assert sum(row["locus_id"] == "locus_1" for row in meta.values()) == 6
    assert "left out by the gene cap per locus, added at locus_1, whose span holds it" in err
    assert "Panel: 5 locus + 0 family + q's own gene = 6 entries." in err


def test_the_report_keeps_the_true_ortholog_and_records_the_entry(tmp_path):
    """End to end through the ownership check: with q as the GOI entry, an ortholog of q is owned
    by the GOI (kept, not re-attributed), an ortholog of the neighbour is relabelled
    paralog_not_goi, and locus_ownership_summary.goi_entry says how the entry was found."""
    _meta, _labels, _err = _run_panel(
        tmp_path,
        genes=[("c1", 1000, 2000, "gene-G"), ("c1", 400, 900, "gene-N1"), ("c1", 2100, 2600, "gene-N2")],
        proteins=[("gene-N1", NEIGHBOUR_PARALOG), ("gene-N2", NEIGHBOUR_OTHER)],
        beds=[("locus_1.bed", [("c1", 1100, 1900)])],
        goi_info=_goi_info(tmp_path / "goi_info.json"),
        extra=["--pad", "1000"])
    ortholog_q, ortholog_n1 = _mutate(Q, 0.1, 11), _mutate(NEIGHBOUR_PARALOG, 0.1, 12)
    target_faa = _fasta(tmp_path / "t.faa", [("GOI_a", ortholog_q), ("GOI_b", ortholog_n1)])
    target_gff = _write(tmp_path / "t.gff", "".join(
        f"t1\tsynvoy\tmRNA\t{s}\t{e}\t.\t+\t.\tID={i};SynVoyRole=goi;Confidence=HIGH\n"
        for i, s, e in (("GOI_a", 100, 600), ("GOI_b", 5000, 5500))))
    out = tmp_path / "own.tsv"
    r = subprocess.run([sys.executable, str(BIN / "reciprocal_best_paralog_check.py"),
                        "--home_query", str(tmp_path / "panel.faa"), "--target_faa", target_faa,
                        "--target_gff", target_gff, "--locus_id", "locus_1", "--genome_name", "tg",
                        "--output", str(out)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    for row in rows:
        row.update(best_bit_f=float(row["best_bit"]), second_bit_f=float(row["second_bit"] or 0),
                   gap_f=float(row["bitscore_gap"]), start_i=int(row["start"]), end_i=int(row["end"]))
    assert {row["mrna_id"]: row["best_paralog"] for row in rows} == \
        {"GOI_a": "locus_1|q1", "GOI_b": "locus_1|gene-N1"}
    records = [{"genome": "tg", "chrom": "t1", "start": s, "end": e, "confidence": "HIGH",
                "identity": "90", "goi_class": "confident_goi", "provenance": ["locus_1"]}
               for s, e in ((100, 600), (5000, 5500))]
    panel_meta = gr.load_panel_meta(str(tmp_path / "panel.faa") + ".meta.tsv")
    _flags, summary = gr.build_locus_ownership({"records": records}, rows, panel_meta, {})
    assert records[0]["goi_class"] == "confident_goi" and records[0]["owning_is_goi_gene"]
    assert records[0]["locus_reattributed"] is False
    assert records[1]["goi_class"] == "paralog_not_goi" and records[1]["owning_gene"] == "gene-N1"
    assert summary["goi_entry"] == {"gene": "q1", "panel_ids": ["locus_1|q1"], "found_by": "query"}


def test_the_summary_records_the_entry_without_rows_and_for_an_old_panel(tmp_path):
    meta_path = tmp_path / "old.meta.tsv"
    meta_path.write_text("panel_id\tlocus_id\tgene\tchrom\tstart\tend\tis_goi_gene\n"
                         "locus_1|gene-Q\tlocus_1\tgene-Q\tc1\t1\t9\t1\n"
                         "locus_1|gene-N\tlocus_1\tgene-N\tc1\t20\t30\t0\n")
    old = gr.load_panel_meta(str(meta_path))
    _flags, summary = gr.build_locus_ownership({"records": []}, [], old, {})
    assert summary["goi_entry"] == {"gene": "gene-Q", "panel_ids": ["locus_1|gene-Q"],
                                    "found_by": "unrecorded"}
    assert gr.build_locus_ownership({"records": []}, [], {}, {})[1]["goi_entry"] is None


def test_primary_locus_panel_name_and_home_record_helpers(tmp_path):
    exons = tmp_path / "goi_exons.faa"
    exons.write_text(">GOI_x|exon_1\nMKV\n>GOI_x\nMKVLLA\n>GOI_copy_x_2\nMKVLLAAAA\n>gene-y\nMKVLLAAAAAA\n")
    assert bhp.home_goi_record(str(exons)) == ("GOI_x", "MKVLLA")
    exons.write_text(">GOI_lcl|ORF1\nMKVLLA\n>GOI_lcl|ORF1|exon_1\nMKVLLAW\n>GOI_lcl|ORF1|frag_1\nMKVLLAWW\n")
    assert bhp.home_goi_record(str(exons)) == ("GOI_lcl|ORF1", "MKVLLA")
    exons.write_text(">GOI_x|exon_1\nMKV\n")
    assert bhp.home_goi_record(str(exons)) is None
    beds = [_write(tmp_path / "locus_2.bed", "c1\t5000\t6000\tg\t1e-30\t+\n"),
            _write(tmp_path / "locus_1.bed", "c1\t5500\t9000\tg\t1e-30\t+\n")]
    assert bhp.locus_holding(beds, "c1", 5800, 5900, 0) == "locus_1"      # the lowest index first
    assert bhp.locus_holding(beds, "c1", 4500, 4700, 0) is None
    assert bhp.locus_holding(beds, "c1", 4500, 4700, 400) == "locus_2"    # padded as genes_for_locus pads
    assert bhp.locus_holding(beds, "c2", 5800, 5900, 0) is None
    assert bhp.locus_holding(beds, "", 0, 0, 0) is None
    assert bhp.home_goi_record(str(tmp_path / "missing.faa")) is None
    assert bhp.home_goi_record(None) is None
    assert bhp.primary_locus_bed(["a/locus_3.bed", "a/locus_1.bed", "a/locus_12.bed"]) == "a/locus_1.bed"
    assert bhp.primary_locus_bed(["x/odd.bed", "x/locus_7.bed"]) == "x/locus_7.bed"
    assert bhp.primary_locus_bed([]) is None
    assert bhp.panel_name("sp|P07585|PGS2_HUMAN") == "sp_P07585_PGS2_HUMAN"
    assert bhp.panel_name("orf19.3940.1") == "orf19.3940.1"


def test_the_panel_gets_goi_info_and_the_home_protein():
    module = (ROOT / "modules" / "assign_locus_ownership.nf").read_text()
    panel = module.split("process BUILD_HOME_PARALOG_PANEL", 1)[1].split("process ", 1)[0]
    assert "path goi_info" in panel and "path home_goi" in panel
    assert "--goi_info ${goi_info}" in panel and "--home_goi ${home_goi}" in panel
    sub = (ROOT / "subworkflows" / "adjudicate_and_report.nf").read_text()
    call = sub.split("BUILD_HOME_PARALOG_PANEL(", 1)[1].split("\n    )", 1)[0]
    args = [ln.strip().rstrip(",") for ln in call.splitlines() if ln.strip()]
    assert args[-3:] == ["effective_params", "goi_info", "goi_exons"]
    take = sub.split("take:", 1)[1].split("main:", 1)[0]
    names = [ln.split("//")[0].strip() for ln in take.splitlines() if ln.split("//")[0].strip()]
    assert names[-2:] == ["goi_info", "goi_exons"]
    main = (ROOT / "main.nf").read_text()
    call = main.split("ADJUDICATE_AND_REPORT(", 1)[1].split("\n        )", 1)[0]
    args = [ln.strip().rstrip(",") for ln in call.splitlines() if ln.strip()]
    assert args[-2:] == ["ANNOTATE_GOI.out.info.first()", "ANNOTATE_GOI.out.exons.first()"]
    assert len(args) == len(names)
