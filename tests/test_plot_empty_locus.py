"""A home locus without flanking genes gets a page that says so, not a failed run.

On a short contig a home locus can have no annotated gene beside the GOI, and its
synteny-block BED is then empty. plot_synteny.py stopped there with exit 2, after the
report was written, and PLOT_SYNTENY's failure ended the whole run (three cluster runs
that started from short contigs, 1-2 Nov 2026). Now it warns, writes the plot page with
the reason, draws nothing else and exits 0. A home BED that does not exist is still an
error.
"""
import os
import re
import subprocess
import sys

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(REPO, "bin", "plot_synteny.py")


def _plot(tmp_path, bed_text=None):
    bed = tmp_path / "synteny_block_locus_5.bed"
    if bed_text is not None:
        bed.write_text(bed_text)
    out = tmp_path / "synteny_block_locus_5_synteny_plot.html"
    proc = subprocess.run([sys.executable, SCRIPT, "--home_bed", str(bed), "--output", str(out),
                           "--no_network"], capture_output=True, text=True, cwd=tmp_path)
    return proc, out


@pytest.mark.parametrize("bed_text", ["", "# synteny block of locus_5\n\n"])
def test_a_locus_without_flanking_genes_gets_a_page_and_the_run_goes_on(tmp_path, bed_text):
    proc, out = _plot(tmp_path, bed_text)
    assert proc.returncode == 0, proc.stderr
    assert "WARNING: no flanking genes in the home BED" in proc.stderr
    assert "its figures are skipped" in proc.stderr
    page = out.read_text()
    assert "No synteny plot for this locus" in page and "nothing to draw" in page
    # Nothing else is drawn (no anchor grid, names table or _with_fragments variant).
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "synteny_block_locus_5.bed", "synteny_block_locus_5_synteny_plot.html"]


def test_a_home_bed_that_does_not_exist_still_stops(tmp_path):
    proc, out = _plot(tmp_path)
    assert proc.returncode == 2
    assert "ERROR: home BED not found" in proc.stderr
    assert not out.exists()


def test_the_module_needs_only_the_plot_page():
    """The skip writes the plot page alone, so every other output of PLOT_SYNTENY must stay
    optional; a required one would end the run again."""
    module = open(os.path.join(REPO, "modules", "plot_synteny.nf")).read()
    block = module.split("output:", 1)[1].split("script:", 1)[0]
    outputs = [ln.strip() for ln in block.splitlines() if ln.strip().startswith("path ")]
    assert outputs, "no outputs found"
    required = [ln for ln in outputs if not re.search(r"optional:\s*true", ln)]
    assert required == ['path "*_synteny_plot.html", emit: plot']
