"""The figures choose what a target row shows by the search's blocks (plot selection, step 3,
2026-10-21; switch plot_block_selection).

With block tables (blocks/<locus>/<locus>_<genome>.blocks.tsv, copied into plot_inputs_* when
the switch is on) a row shows the block of every HIGH, MEDIUM and AMBIGUOUS GOI call, then the
significant blocks (E_syn <= 0.1) without one, up to 6 blocks: rule "Calls+significant". A LOW
call draws no block of its own and a non-significant block without a call is not drawn.
Fragment GOI models take part in the choice and stay hidden in the default figures, so the
'_with_fragments' figures draw the same blocks. A genome without its table, and every genome of
a run without tables, keeps the candidate regions, which only the console says.

The fixture's first target genome holds four blocks: 0, the neighbourhood with a HIGH call;
1, a non-significant block with a LOW call inside a candidate region; 2, a significant block
without a call outside every region; 3, a non-significant block whose only call is a MEDIUM
fragment. The candidate regions show block 0 alone; the rule shows blocks 0, 2 and 3.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import plot_synteny as ps  # noqa: E402

REPLOT = ROOT / "scripts" / "replot.sh"
LOCUS = "plot_inputs_synteny_block_locus_1"
G1, G2, G3 = "GCA_900000001.1", "GCA_900000002.1", "GCA_900000003.1"
HOME = [("FLA", 1000, 2000, "+"), ("FLB", 4000, 5000, "-"), ("QRY1", 7000, 8200, "+"),
        ("FLC", 10000, 11000, "+"), ("FLD", 13000, 14000, "-")]
COLUMNS = ["genome", "block", "chrom", "start", "end", "hood_start", "hood_end", "esyn",
           "significant"]


def _tag(genome):
    return genome.replace(".", "_") + "_fna"


def _flank(genome, gene, chrom, s, e, st, block):
    mid = f"gene-{gene}|{_tag(genome)}_b{block}_fl1_flank_ann" if block is not None \
        else f"gene-{gene}|{_tag(genome)}_rearranged"
    return [f"{chrom}\tflanking_annotation\tmRNA\t{s}\t{e}\t80.0\t{st}\t.\tID={mid};Name=gene-{gene};"
            f"SynVoy_Parent=gene-{gene};Type=flanking_miniprot;Identity=80.0;SynVoyRole=flanking;"
            "EvidenceType=flanking_miniprot;Confidence=HIGH;QueryCoverage=0.950;Exons=1",
            f"{chrom}\tflanking_annotation\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={mid}_CDS1;Parent={mid}"]


def _goi(genome, chrom, s, e, block, conf, status, esyn, suffix="l1_exon_ann"):
    mid = f"GOI_QRY1|{_tag(genome)}_b{block}_{suffix}"
    return [f"{chrom}\texon_annotation\tmRNA\t{s}\t{e}\t70.0\t+\t.\tID={mid};Name=GOI_QRY1;"
            f"SynVoy_Parent=GOI_QRY1;Identity=70.0;SynVoyRole=goi;EvidenceType=exon_annotation;"
            f"Confidence={conf};Basis=significant_specific;GOIClass=confident_goi;ModelStatus={status};"
            f"QueryCoverage=0.950;Exons=1;Bits=200;Eeff=1e-30;BlockEsyn={esyn}",
            f"{chrom}\texon_annotation\tCDS\t{s}\t{e}\t.\t+\t0\tID={mid}_CDS1;Parent={mid}"]


def _table(path, genome, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        for b, chrom, s, e, esyn in rows:
            w.writerow({"genome": f"{genome}.fna", "block": b, "chrom": chrom, "start": s, "end": e,
                        "hood_start": s - 5000, "hood_end": e + 5000, "esyn": repr(esyn),
                        "significant": int(esyn <= 0.1)})


def _folder(d: Path, tables=True) -> Path:
    d.mkdir(parents=True)
    (d / "synteny_block_locus_1.bed").write_text("".join(
        f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t{st}\t{g}\n" for g, s, e, st in HOME))
    (d / "locus_1.bed").write_text("hchr1\t7100\t8100\tgene_loc\t1e-40\t+\n")
    home = ["##gff-version 3"]
    for g, s, e, st in HOME:
        a = f"gene={g};product={g} protein"
        home += [f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{g};Name={g};{a}",
                 f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID=rna-{g};Parent=gene-{g};{a}",
                 f"hchr1\tRefSeq\tCDS\t{s}\t{e}\t.\t{st}\t0\tID=cds-{g};Parent=rna-{g};{a}"]
    (d / "home_genome.gff").write_text("\n".join(home) + "\n")
    g1 = ["##gff-version 3"]
    for g, s, e, st in HOME:                                    # block 0: the neighbourhood
        if g == "QRY1":
            g1 += _goi(G1, "tchr5", s + 500_000, e + 500_000, 0, "HIGH", "complete", "1e-20")
        else:
            g1 += _flank(G1, g, "tchr5", s + 500_000, e + 500_000, st, 0)
    g1 += _flank(G1, "FLA", "tchr9", 901_000, 902_000, "+", 1)     # block 1: E_syn 2, a LOW call
    g1 += _flank(G1, "FLB", "tchr9", 904_000, 905_000, "-", 1)
    g1 += _goi(G1, "tchr9", 907_000, 908_000, 1, "LOW", "partial", "2", "l1_fallback")
    g1 += _flank(G1, "FLC", "tchr7", 310_000, 311_000, "+", 2)     # block 2: E_syn 0.01
    g1 += _flank(G1, "FLD", "tchr7", 313_000, 314_000, "-", 2)
    g1 += _flank(G1, "FLB", "tchr8", 104_000, 105_000, "-", 3)     # block 3: a MEDIUM fragment
    g1 += _goi(G1, "tchr8", 107_000, 107_400, 3, "MEDIUM", "fragment", "5", "r1_fallback")
    (d / f"{G1}.fna.gff").write_text("\n".join(g1) + "\n")
    (d / f"{G1}.fna.regions.bed").write_text(
        "tchr5\t500000\t515000\tReg1_G5_CHIGH_S0.90\t0.90\t+\n"
        "tchr9\t900000\t909000\tReg2_G2_CLOW_S0.20\t0.20\t+\n")
    # G2: no block (header-only table); one cross-chromosome flanking model inside a region
    (d / f"{G2}.fna.gff").write_text("\n".join(
        ["##gff-version 3"] + _flank(G2, "FLA", "uchr1", 51_000, 52_000, "+", None)) + "\n")
    (d / f"{G2}.fna.regions.bed").write_text("uchr1\t50000\t53000\tReg1_G1_CLOW_S0.10\t0.10\t+\n")
    # G3: no table at all
    (d / f"{G3}.fna.gff").write_text("\n".join(
        ["##gff-version 3"] + _goi(G3, "vchr1", 20_000, 21_000, 0, "HIGH", "complete", "1e-9")
        + _flank(G3, "FLC", "vchr1", 23_000, 24_000, "+", 0)) + "\n")
    (d / f"{G3}.fna.regions.bed").write_text("vchr1\t19000\t25000\tReg1_G2_CHIGH_S0.80\t0.80\t+\n")
    (d / "species_mapping.tsv").write_text("".join(
        f"{g}\tNonexistus fixturae{k}\ttarget\n" for k, g in enumerate((G1, G2, G3))))
    (d / "locus_1_tree.nwk").write_text("(GOI_placeholder:0.0);\n")
    (d / "gene_names.tsv").write_text("gene_id\tlabel\n" + "".join(
        f"gene-{g}\t{g}\n" for g, _s, _e, _st in HOME))
    if tables:
        _table(d / f"locus_1_{G1}.fna.blocks.tsv", G1,
               [(0, "tchr5", 501_000, 514_000, 1e-20), (1, "tchr9", 901_000, 908_000, 2.0),
                (2, "tchr7", 310_000, 314_000, 0.01), (3, "tchr8", 104_000, 107_400, 5.0)])
        _table(d / f"locus_1_{G2}.fna.blocks.tsv", G2, [])
    return d


def _redraw(folder, out, *extra):
    env = dict(os.environ, PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""))
    proc = subprocess.run(["bash", str(REPLOT), str(folder), str(out), *map(str, extra)],
                          capture_output=True, text=True, env=env)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout


def _blocks_lines(stdout):
    return {m.group(1): m.group(2) for m in re.finditer(r"\[plot\] (\S+): blocks drawn: (.*)", stdout)}


def test_the_table_of_a_gff_is_named_as_the_search_names_it(tmp_path):
    assert ps.block_table_name("x/synteny_block_locus_2.bed", "y/X.fna.gff") == "locus_2_X.fna.blocks.tsv"
    assert ps.genome_tag("y/GCA_1.2.fna.gff") == "GCA_1_2_fna"
    _table(tmp_path / "empty.blocks.tsv", "G", [])
    assert ps.read_block_table(str(tmp_path / "empty.blocks.tsv")) == {}   # a table all the same
    assert ps.read_block_table(str(tmp_path / "missing.blocks.tsv")) is None


def test_calls_then_significant_blocks(tmp_path):
    d = _folder(tmp_path / LOCUS)
    gff = str(d / f"{G1}.fna.gff")
    table = ps.read_block_table(str(d / f"locus_1_{G1}.fna.blocks.tsv"))
    choice = ps.choose_blocks(ps.choice_models(gff), ps.genome_tag(gff), table)
    assert choice.drawn == {0: "call_block", 3: "call_block", 2: "predicted_no_call"}
    assert choice.predicted == {2}
    # the cap counts the call blocks
    assert set(ps.choose_blocks(ps.choice_models(gff), ps.genome_tag(gff), table, max_blocks=2).drawn) == {0, 3}
    # without the fragment the block of the fragment is not drawn (E_syn 5, no other call)
    hidden = [g for g in ps.choice_models(gff) if not ps._is_fragment_goi(g)]
    assert ps.choose_blocks(hidden, ps.genome_tag(gff), table).drawn == {0: "call_block", 2: "predicted_no_call"}


def test_a_run_with_tables_shows_the_blocks_of_the_rule(tmp_path):
    out = _redraw(_folder(tmp_path / LOCUS), tmp_path / "fig")
    lines = _blocks_lines(out)
    # default figures: block 0 (5 models), block 2 (2), block 3 (1: the fragment is hidden);
    # the '_with_fragments' figures draw the same blocks and the fragment
    assert re.findall(rf"\[plot\] {re.escape(G1)}: blocks drawn: (.*)", out) == [
        "b0 call_block, b2 predicted_no_call, b3 call_block; 8 models",
        "b0 call_block, b2 predicted_no_call, b3 call_block; 9 models"]
    assert lines[G2] == "none; 0 models"         # header-only table: an empty row
    assert G3 not in lines
    assert f"[plot] {G3}: no block table locus_1_{G3}.fna.blocks.tsv: its row shows the candidate regions" in out
    with open(next((tmp_path / "fig").glob("*_anchor_grid_print.source.tsv")), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["genome_id"] for r in rows} >= {G1, G2, G3}            # every genome keeps its row
    g1_scaffolds = {r["scaffold"] for r in rows if r["genome_id"] == G1 and r["scaffold"]}
    assert "tchr9" not in g1_scaffolds                               # block 1 is not drawn


def test_a_run_without_tables_keeps_the_candidate_regions(tmp_path):
    out = _redraw(_folder(tmp_path / LOCUS, tables=False), tmp_path / "fig")
    assert "[plot] no block tables (--block_tables): rows show the candidate regions" in out
    assert not _blocks_lines(out)


def test_tables_of_another_locus_change_nothing(tmp_path):
    plain = _folder(tmp_path / "a" / LOCUS, tables=False)
    _redraw(plain, tmp_path / "fig_a")
    other = tmp_path / "b" / LOCUS
    shutil.copytree(plain, other)
    _table(other / f"locus_2_{G1}.fna.blocks.tsv", G1, [(0, "tchr5", 501_000, 514_000, 1e-20)])
    out = _redraw(other, tmp_path / "fig_b")
    assert f"[plot] {G1}: no block table locus_1_{G1}.fna.blocks.tsv" in out
    a = {p.name: p.read_bytes() for p in (tmp_path / "fig_a").iterdir()}
    b = {p.name: p.read_bytes() for p in (tmp_path / "fig_b").iterdir()}
    assert a.keys() == b.keys() and len(a) > 5
    assert [n for n in a if a[n] != b[n]] == []


def test_the_module_copies_the_tables_only_with_the_switch():
    src = (ROOT / "modules" / "plot_synteny.nf").read_text()
    assert "path block_tables" in src
    assert ("def tables_str = (params.plot_block_selection.toString().toBoolean() && block_tables) "
            "? block_tables.join(' ') : ''") in src
    # copied without '|| true': a table that does not arrive must not quietly give the old figure
    assert '        cp ${tables_str} "\\$inputs_dir/"\n' in src
    assert "plot_block_selection = false" in (ROOT / "nextflow.config").read_text()
    main = (ROOT / "main.nf").read_text()
    assert ".join(blocks_by_locus_ch, remainder: true)" in main
    assert "plot_inputs_split.block_tables" in main
    replot = REPLOT.read_text()
    assert 'block_tables=("$IN"/*.blocks.tsv)' in replot
    assert "<run>/blocks/" not in replot.split("set -euo pipefail", 1)[1]   # never the run's blocks/ folder
