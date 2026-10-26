"""scripts/replot_run.sh: every locus of a finished run re-drawn in one command (2026-10-20).

Runs made with older SynVoy code get today's figures without a new search. Two formats of
those runs break a direct re-draw with replot.sh, each in a different figure, so the
wrapper draws such a locus from a renamed copy:

* GFF keys from before the rename (SynTerra_Parent=, SynTerraRole=): the plotter reads only
  SynVoy_Parent and SynVoyRole, no target gene is tied to its home gene, and the ribbon plot
  loses its ribbons (a March LY6E run: 22 of 46 drawn);
* the query BED under its old name (home_gene_location.bed, not locus_<N>.bed): replot.sh
  passes no --query_bed, the plotter finds no home GOI, and the anchor grid has no GOI
  column (the same run: 30 of 54 cells placed, a GOI in 2 of its 9 genomes).

The fixture is a run of that format with one target genome. Drawn as it is, it has neither
ribbons nor a GOI column; that control is what makes the first test bite.

Added 2026-10-26, from a test as a student would run it: the home species is read from the
run when it is not given (a forgotten --home_species left "Home genome" in the home row, a
copied one put a firefly on a human run), the missing taxonomy file is named once before
drawing, a read-only run folder (as unpacked from a student's zip) gets a message instead of
a bare mkdir error, the console ends with the slide figure of every locus, and a copy outside
its own git history names no commit.
"""
from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPLOT = ROOT / "scripts" / "replot.sh"
REPLOT_RUN = ROOT / "scripts" / "replot_run.sh"

LOCUS = "plot_inputs_synteny_block_locus_1"
TARGET = "GCA_900000001.1"
HOME = [("FLA", 1000, 2000, "+"), ("FLB", 4000, 5000, "-"), ("QRY1", 7000, 8200, "+"),
        ("FLC", 10000, 11000, "+"), ("FLD", 13000, 14000, "-")]
SHIFT = 500_000         # the target copy of the neighbourhood lies here on tchr5


def _march_run(run: Path) -> Path:
    """A run folder in the format of March 2026: SynTerra keys, home_gene_location.bed."""
    d = run / LOCUS
    d.mkdir(parents=True)
    (d / "synteny_block_locus_1.bed").write_text("".join(
        f"hchr1\t{s - 1}\t{e}\tgene-{g}\t.\t{st}\t{g}\n" for g, s, e, st in HOME))
    (d / "home_gene_location.bed").write_text("hchr1\t7100\t8100\tgene_loc\t1e-40\t+\n")
    home = ["##gff-version 3"]
    for g, s, e, st in HOME:
        a = f"gene={g};product={g} protein"
        home.append(f"hchr1\tRefSeq\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{g};Name={g};{a}")
        home.append(f"hchr1\tRefSeq\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID=rna-{g};Parent=gene-{g};{a}")
        exons = [(s, s + 299), (e - 299, e)] if g == "QRY1" else [(s, e)]
        home += [f"hchr1\tRefSeq\tCDS\t{a0}\t{b0}\t.\t{st}\t0\tID=cds-{g};Parent=rna-{g};{a}"
                 for a0, b0 in exons]
    (d / "home_genome.gff").write_text("\n".join(home) + "\n")
    tag = TARGET.replace(".", "_") + "_fna"
    target = ["##gff-version 3"]
    for g, s, e, st in HOME:
        s, e = s + SHIFT, e + SHIFT
        if g == "QRY1":
            mid = f"GOI_{g}|{tag}_b0_l1_exon_ann"
            target.append(f"tchr5\texon_annotation\tmRNA\t{s}\t{e}\t70.0\t{st}\t.\tID={mid};"
                          f"Name=GOI_{g};SynTerra_Parent=GOI_{g};Identity=70.0;SynTerraRole=goi;"
                          "EvidenceType=exon_annotation;Confidence=HIGH;GOIClass=confident_goi;"
                          "BlockFlankingSupport=4;QueryCoverage=0.950;Exons=2")
            target += [f"tchr5\texon_annotation\tCDS\t{a0}\t{b0}\t.\t{st}\t0\tID={mid}_CDS{k};Parent={mid}"
                       for k, (a0, b0) in enumerate([(s, s + 299), (e - 299, e)], 1)]
        else:
            mid = f"gene-{g}|{tag}_b0_fl1_flank_ann"
            target.append(f"tchr5\tflanking_annotation\tmRNA\t{s}\t{e}\t80.0\t{st}\t.\tID={mid};"
                          f"Name=gene-{g};SynTerra_Parent=gene-{g};Type=flanking_miniprot;"
                          "Identity=80.0;SynTerraRole=flanking;EvidenceType=flanking_miniprot;"
                          "Confidence=HIGH;QueryCoverage=0.950;Exons=1")
            target.append(f"tchr5\tflanking_annotation\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={mid}_CDS1;Parent={mid}")
    (d / f"{TARGET}.fna.gff").write_text("\n".join(target) + "\n")
    (d / f"{TARGET}.fna.regions.bed").write_text(
        f"tchr5\t{SHIFT}\t{SHIFT + 15000}\tReg1_G5_CHIGH_S0.90\t0.90\t+\n")
    # Not in the NCBI taxonomy: the rows follow the GOI gene tree with or without ete3.
    (d / "species_mapping.tsv").write_text(f"{TARGET}\tNonexistus fixturae\ttarget\n")
    (d / "locus_1_tree.nwk").write_text("(GOI_placeholder:0.0);\n")
    return run


def _current_run(run: Path) -> Path:
    """The same run in today's format: SynVoy keys, locus_1.bed, labels with Basis=."""
    d = _march_run(run) / LOCUS
    for gff in d.glob("*.gff"):
        gff.write_text(gff.read_text().replace("SynTerra_Parent=", "SynVoy_Parent=")
                       .replace("SynTerraRole=", "SynVoyRole=")
                       .replace("Confidence=HIGH;GOIClass", "Confidence=HIGH;Basis=significant_specific;GOIClass"))
    (d / "home_gene_location.bed").rename(d / "locus_1.bed")
    return run


def _redraw(*args, **env_extra) -> subprocess.CompletedProcess:
    """Run a script with this interpreter first on PATH (replot.sh calls `python3`)."""
    env = dict(os.environ, PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", ""),
               **env_extra)
    return subprocess.run(["bash", *map(str, args)], capture_output=True, text=True, env=env)


def _files(run: Path) -> dict:
    """Every file of the run folder except the re-drawn figures: path -> (bytes, mode)."""
    return {str(p.relative_to(run)): (p.read_bytes(), p.stat().st_mode)
            for p in sorted(run.rglob("*"))
            if p.is_file() and not p.relative_to(run).parts[0].startswith("figures_redrawn_")}


def _figures(out: Path) -> tuple[list[dict], int]:
    """The anchor grid's source rows and the number of ribbons in the ribbon plot."""
    (src,) = out.glob("*_anchor_grid_print.source.tsv")
    with open(src, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    (view,) = out.glob("*_synteny_plot_view.svg")
    return rows, view.read_text().count('class="ribbon"')


def _goi_states(rows):
    return [r["state"] for r in rows if r["genome_id"] == TARGET and r["column"] == "QRY1"]


def test_an_older_run_is_drawn_from_a_renamed_copy(tmp_path):
    run = _march_run(tmp_path / "run")
    before = _files(run)
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert f"[{LOCUS}] older run format: drawing from a renamed copy" in proc.stdout
    (out,) = run.glob("figures_redrawn_*")
    rows, ribbons = _figures(out)
    assert _goi_states(rows) == ["placed"]          # needs the query BED under its new name
    assert ribbons == 5                              # 4 flanking genes + the GOI: needs the keys
    redrawn = (out / "REDRAWN.txt").read_text()
    assert "only the drawing is new" in redrawn
    assert f"Labels from code before the label rework of 2026-10-15, not today's labels: {LOCUS}\n" in redrawn
    assert f"Older run format, drawn from a renamed copy: {LOCUS}\n" in redrawn
    assert f"Rows ordered by the GOI gene tree, not the species tree: {LOCUS}\n" in redrawn
    assert "Loci: 1, failed: 0" in redrawn
    assert (out / f"{LOCUS}.log").is_file()
    assert proc.stdout.endswith(
        f"  {out}/synteny_block_locus_1_anchor_grid_print.svg\n"
        "The labels of every locus come from code before 15 Oct 2026: check numbered copies before counting them.\n")
    assert _files(run) == before                     # the run itself is unchanged


def test_drawn_as_it_is_the_fixture_has_neither_ribbons_nor_a_goi(tmp_path):
    run = _march_run(tmp_path / "run")
    proc = _redraw(REPLOT, run / LOCUS, tmp_path / "raw")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    rows, ribbons = _figures(tmp_path / "raw")
    assert _goi_states(rows) == [] and ribbons == 0, (
        "the plotter now reads this format itself: the conversion in replot_run.sh "
        "and this test can go")


def test_a_run_in_todays_format_is_drawn_from_the_run_itself(tmp_path):
    run = _current_run(tmp_path / "run")
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "older run format" not in proc.stdout and "label rework" not in proc.stdout
    (out,) = run.glob("figures_redrawn_*")
    rows, ribbons = _figures(out)
    assert _goi_states(rows) == ["placed"] and ribbons == 5
    redrawn = (out / "REDRAWN.txt").read_text()
    assert "Older run format" not in redrawn and "label rework" not in redrawn


def test_a_failed_locus_fails_the_run_and_the_others_are_drawn(tmp_path):
    run = _march_run(tmp_path / "run")
    (run / "plot_inputs_synteny_block_locus_2").mkdir()    # no synteny_block_*.bed
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode != 0
    assert "[plot_inputs_synteny_block_locus_2] FAILED" in proc.stderr
    (out,) = run.glob("figures_redrawn_*")
    assert "Loci: 2, failed: 1" in (out / "REDRAWN.txt").read_text()
    assert "no synteny_block_*.bed" in (out / "plot_inputs_synteny_block_locus_2.log").read_text()
    assert _goi_states(_figures(out)[0]) == ["placed"]
    assert (f"  [plot_inputs_synteny_block_locus_2] no print grid was drawn, see {out}/plot_inputs_synteny_block_locus_2.log"
            in proc.stdout)


# --- 2026-10-26: what a student runs into ----------------------------------------------------

def _logs(run: Path, plot=(), phylo=()) -> Path:
    """Task logs as main.nf (collectTaskLogs) keeps them: one folder per task, a counter for a
    repeat. A work folder shared by several runs leaves the tasks of all of them there."""
    for kind, values, line in (("PLOT_SYNTENY__1_____synteny_plot.html", plot, '        --home_species "{}" \\'),
                               ("PHYLO_SORT__sorting___sorted_genomes.txt", phylo, 'HOME_QUERY="{}"')):
        for i, value in enumerate(values, 1):
            d = run / "logs" / (kind + (f"__{i}" if i > 1 else ""))
            d.mkdir(parents=True)
            (d / "script.sh").write_text("#!/bin/bash -ue\n" + line.format(value) + "\n")
    return run


def _home_row(out: Path) -> set:
    """The species the print grid gives its home row."""
    return {r["species"] for r in _figures(out)[0] if r["genome_id"] == "home"}


def test_the_home_species_is_read_from_the_runs_logs(tmp_path):
    # The value the pipeline used, here with the file path that leaked into it in May.
    run = _logs(_current_run(tmp_path / "run"), phylo=["data/genomes/Apis_mellifera.fna"])
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Home species: Apis mellifera (read from the run: its task logs;" in proc.stdout
    (out,) = run.glob("figures_redrawn_*")
    assert _home_row(out) == {"Apis mellifera"}
    assert "Home species: Apis mellifera (its task logs)\n" in (out / "REDRAWN.txt").read_text()


def test_several_logged_species_are_told_apart_by_the_home_gff_name(tmp_path):
    run = _logs(_current_run(tmp_path / "run"), plot=["Homo sapiens"], phylo=["Apis_mellifera", "human"])
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    # 'human' is no species name; the two that are cannot be told apart by a home_genome.gff
    assert "the run's task logs name several (Apis mellifera, Homo sapiens)" in proc.stdout
    (out,) = run.glob("figures_redrawn_*")
    assert _home_row(out) == {"Home genome"}
    shutil.rmtree(out)
    (run / LOCUS / "home_genome.gff").rename(run / LOCUS / "Apis_mellifera.gff")
    proc = _redraw(REPLOT_RUN, run)
    assert "Home species: Apis mellifera (read from the run: its task logs and the file name of its home GFF;" in proc.stdout
    (out,) = run.glob("figures_redrawn_*")
    assert _home_row(out) == {"Apis mellifera"}
    shutil.rmtree(out)
    shutil.rmtree(run / "logs")                     # no logs: the file name alone
    proc = _redraw(REPLOT_RUN, run)
    assert "Home species: Apis mellifera (read from the run: the file name of its home GFF;" in proc.stdout


def test_a_given_home_species_is_used_and_a_different_one_in_the_run_is_flagged(tmp_path):
    plain = _current_run(tmp_path / "plain")
    logged = _logs(_current_run(tmp_path / "logged"), plot=["Homo sapiens"])
    for run in (plain, logged):
        proc = _redraw(REPLOT_RUN, run, "--home_species", "Apis mellifera")
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "Home species: Apis mellifera (as given)" in proc.stdout
    assert 'WARNING: the run names its home species Homo sapiens (its task logs), not "Apis mellifera"' in proc.stdout
    (a,) = plain.glob("figures_redrawn_*")
    (b,) = logged.glob("figures_redrawn_*")
    assert _home_row(b) == {"Apis mellifera"}

    def figures(d):
        return {p.name: p.read_bytes() for p in d.iterdir() if p.suffix != ".log" and p.name != "REDRAWN.txt"}
    assert figures(a) == figures(b)                  # given, the run's own species changes nothing


def test_a_missing_taxonomy_file_is_named_once_with_the_command_that_gets_it(tmp_path):
    run = _logs(_current_run(tmp_path / "run"), phylo=["Apis mellifera"])
    home = tmp_path / "home"
    home.mkdir()
    proc = _redraw(REPLOT_RUN, run, HOME=str(home))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    if importlib.util.find_spec("ete3") is None:     # CI: no ete3
        assert "cannot import ete3" in proc.stdout and "conda activate synvoy_env" in proc.stdout
        cause = "no ete3"
    else:
        assert "the NCBI taxonomy file ~/.etetoolkit/taxa.sqlite is missing" in proc.stdout
        assert "(cd ~ && python3 -c 'from ete3 import NCBITaxa; NCBITaxa()' && rm -f taxdump.tar.gz)" in proc.stdout
        cause = "no taxonomy file"
    assert proc.stdout.count("NOTE:") == 1
    assert f"[{LOCUS}] rows ordered by the GOI gene tree, not the species tree ({cause}, see the NOTE above)" in proc.stdout


def test_the_console_ends_with_the_slide_figure_of_every_locus(tmp_path):
    run = _current_run(tmp_path / "run")
    for n in (2, 10):                                # in folder order 10 comes first
        d = run / f"plot_inputs_synteny_block_locus_{n}"
        shutil.copytree(run / LOCUS, d)
        (d / "synteny_block_locus_1.bed").rename(d / f"synteny_block_locus_{n}.bed")
        (d / "locus_1.bed").rename(d / f"locus_{n}.bed")
    shutil.rmtree(run / LOCUS)
    gff = run / "plot_inputs_synteny_block_locus_10" / f"{TARGET}.fna.gff"
    gff.write_text(gff.read_text().replace("Basis=significant_specific;", ""))      # labels of before 15 Oct
    proc = _redraw(REPLOT_RUN, run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (out,) = run.glob("figures_redrawn_*")
    listing = proc.stdout.split("For a slide, the anchor grid at print size of each locus:\n", 1)[1]
    assert listing.splitlines() == [
        f"  {out}/synteny_block_locus_2_anchor_grid_print.svg",
        f"  {out}/synteny_block_locus_10_anchor_grid_print.svg"
        "  (labels from before 15 Oct 2026: check numbered copies before counting them)"]


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root writes into a read-only folder")
def test_a_read_only_run_folder_gets_a_message_with_the_fix(tmp_path):
    run = _march_run(tmp_path / "run")
    before = _files(run)
    folders = [run, *(p for p in run.rglob("*") if p.is_dir())]
    for d in folders:
        d.chmod(0o555)                               # dr-xr-xr-x, as a student's zip unpacks
    try:
        proc = _redraw(REPLOT_RUN, run)
    finally:
        for d in folders:
            d.chmod(0o755)
    assert proc.returncode == 1
    assert "the folder is read-only" in proc.stderr and f'chmod -R u+w "{run}"' in proc.stderr
    assert not list(run.glob("figures_redrawn_*")) and _files(run) == before


@pytest.mark.skipif(shutil.which("git") is None, reason="needs git")
def test_a_copy_outside_its_own_git_history_names_no_commit(tmp_path):
    # An archive of the code unpacked inside another repository: that repository's commit
    # is not the figure code's.
    code = tmp_path / "outer" / "SynVoy-archive"
    (code / "scripts").mkdir(parents=True)
    for script in ("replot.sh", "replot_run.sh"):
        shutil.copy2(ROOT / "scripts" / script, code / "scripts" / script)
    shutil.copy2(ROOT / "VERSION", code / "VERSION")
    (code / "bin").symlink_to(ROOT / "bin")
    subprocess.run(["git", "init", "-q", str(tmp_path / "outer")], check=True)
    run = _march_run(tmp_path / "run")
    proc = _redraw(code / "scripts" / "replot_run.sh", run)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (out,) = run.glob("figures_redrawn_*")
    version = (ROOT / "VERSION").read_text().strip()
    assert (f"Figure code: SynVoy {version}, commit unknown (not a git checkout: SynVoy-archive)\n"
            in (out / "REDRAWN.txt").read_text())
