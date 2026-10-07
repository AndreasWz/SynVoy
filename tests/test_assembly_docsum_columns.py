"""Assembly quality columns must stay aligned in every NCBI docsum query.

An NCBI assembly docsum has no top-level <ScaffoldCount>/<ContigCount>. ``xtract
-element`` drops fields that are absent, so without ``-def NA`` the N50 values
slide into the count columns. tests/test_fetch_home_quality.py pinned this for
get_assembly_quality(); two other call sites had the same defect until 2026-10-05:

* fetch_related_genomes.get_related_species (every easy-mode target search), and
* fetch_home_genome.find_any_genome (home species without a reference genome).

Seen on a real run (yeast, easy mode): *S. paradoxus* GCF_002079055.1 was logged
as ``scaf=903028, contigs=835812, N50=NA/NA`` -- those are its scaffold and contig
N50. The effects: a scaffold-level assembly with an N50 above 500 kb failed the
quality gate ("scaffold_count > 500000") and was dropped, and within one species
the assembly with the SMALLEST contig N50 ranked first.
"""
from __future__ import annotations

import ast
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import fetch_home_genome as fh  # noqa: E402
import fetch_related_genomes as fr  # noqa: E402

FIELDS = ["AssemblyAccession", "SpeciesName", "RefSeq_category", "AssemblyStatus",
          "ScaffoldCount", "ContigCount", "ScaffoldN50", "ContigN50",
          "ScaffoldN80", "ContigN80"]

# Real `xtract -def NA` rows (2026-10-05).
PARADOXUS = ("GCF_002079055.1\tSaccharomyces paradoxus\treference genome\tChromosome"
             "\tNA\tNA\t903028\t835812\tNA\tNA")
PHOTINUS = ("GCF_008802855.1\tPhotinus pyralis\treference genome\tScaffold"
            "\tNA\tNA\t47017841\t170308\tNA\tNA")
# What the old call printed for the same record: six columns, N50s in the count slots.
PHOTINUS_SHIFTED = "GCF_008802855.1\tPhotinus pyralis\treference genome\tScaffold\t47017841\t170308"


def _gate_args():
    return types.SimpleNamespace(bad_max_contigs=500000, bad_max_scaffolds=500000,
                                 bad_min_n50=5000)


def _docsum_xtract_lists(path: Path):
    """Every list literal in *path* that is an xtract call on assembly docsums."""
    found = []
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.List):
            items = [e.value for e in node.elts
                     if isinstance(e, ast.Constant) and isinstance(e.value, str)]
            if "DocumentSummary" in items and "ScaffoldCount" in items:
                found.append(items)
    return found


@pytest.mark.parametrize("script", ["fetch_related_genomes.py", "fetch_home_genome.py"])
def test_every_docsum_query_pads_absent_fields(script):
    calls = _docsum_xtract_lists(ROOT / "bin" / script)
    assert calls, f"{script}: no assembly-docsum xtract call found"
    for items in calls:
        assert "-def" in items and items[items.index("-def") + 1] == "NA", items
        assert items.index("-def") < items.index("-element"), items
        assert items[items.index("-element") + 1:] == FIELDS, items


def test_related_row_keeps_n50_out_of_the_count_columns():
    e = fr.parse_assembly_row(PARADOXUS, "genus (Saccharomyces)")
    assert e["scaffold_count"] is None and e["contig_count"] is None
    assert e["scaffold_n50"] == 903028 and e["contig_n50"] == 835812
    assert e["assembly_status"] == "Chromosome"
    assert e["category"] == "reference genome"
    assert e["tax_level"] == "genus (Saccharomyces)"


def test_related_row_treats_placeholders_as_absent():
    e = fr.parse_assembly_row("GCA_1.1\tGenus species\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")
    assert e["category"] == "" and e["assembly_status"] == ""
    assert e["scaffold_n50"] is None and e["contig_n50"] is None
    assert fr.parse_assembly_row("NA\tGenus species") is None
    assert fr.parse_assembly_row("GCA_1.1\tNA") is None


def test_contiguous_scaffold_level_assembly_passes_the_gate():
    good = fr.parse_assembly_row(PHOTINUS)
    assert fr.is_bad_quality(good, _gate_args()) == (False, [])
    # The shifted row is what the gate used to see: a 47 Mb N50 read as a count.
    shifted = fr.parse_assembly_row(PHOTINUS_SHIFTED)
    bad, reasons = fr.is_bad_quality(shifted, _gate_args())
    assert bad and "scaffold_count=47017841" in reasons[0]


@pytest.mark.parametrize("mode", ["hybrid", "nstats", "counts"])
def test_more_contiguous_assembly_ranks_first(mode):
    rows = ["GCA_000000001.1\tGenus species\tna\tScaffold\tNA\tNA\t400000\t100000\tNA\tNA",
            "GCA_000000002.1\tGenus species\tna\tScaffold\tNA\tNA\t9000000\t5000000\tNA\tNA"]
    for module, parse in ((fr, fr.parse_assembly_row), (fh, fh.parse_quality_line)):
        entries = [parse(r) for r in rows]
        best = min(entries, key=lambda x: module.assembly_rank_tuple(x, mode))
        assert best["accession"] == "GCA_000000002.1", (module.__name__, mode)


@pytest.mark.skipif(shutil.which("xtract") is None, reason="entrez-direct not installed")
def test_real_xtract_emits_ten_columns_for_a_docsum_without_counts():
    xml = ("<eSummaryResult><DocumentSummarySet><DocumentSummary>"
           "<AssemblyAccession>GCF_002079055.1</AssemblyAccession>"
           "<SpeciesName>Saccharomyces paradoxus</SpeciesName>"
           "<AssemblyStatus>Chromosome</AssemblyStatus>"
           "<ContigN50>835812</ContigN50><ScaffoldN50>903028</ScaffoldN50>"
           "<RefSeq_category>reference genome</RefSeq_category>"
           "</DocumentSummary></DocumentSummarySet></eSummaryResult>")
    out = subprocess.run(fr.ASSEMBLY_DOCSUM_XTRACT, input=xml, capture_output=True,
                         text=True, check=True).stdout.strip()
    assert out.split("\t") == ["GCF_002079055.1", "Saccharomyces paradoxus",
                               "reference genome", "Chromosome", "NA", "NA",
                               "903028", "835812", "NA", "NA"]
    assert fr.parse_assembly_row(out)["scaffold_n50"] == 903028
