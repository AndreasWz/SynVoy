"""Tests for bin/rescue_strong_synteny.py (docs/TODO.md §1e follow-up).

The rescue script is invoked per (locus, genome) by RESCUE_STRONG_SYNTENY when
cluster_grs.py has classified at least one block as ``goi_missing_but_strong_synteny``.
We don't exercise the actual ``miniprot`` subprocess here (it needs a real binary +
genome); we pin the deterministic Python pieces:

1. ``scores.tsv`` with NO ``goi_missing_but_strong_synteny`` rows → empty (header-
   only) GFF (the channel must still produce a valid file).
2. Miniprot --gff parsing rebuilds the mRNA + CDS structure correctly.
3. ``_build_rescue_gff_rows`` maps window-relative coords back to genome coords
   and tags rows with the right SynVoy attributes so generate_report.py picks
   them up as LOW-confidence GOI rescues.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import rescue_strong_synteny as rss  # noqa: E402


def _write_scores(path: Path, rows: list[dict]) -> None:
    fields = [
        "region_rank", "region_name", "species", "chrom", "start", "end",
        "strand", "score", "quality_score", "coverage_score", "unique_genes",
        "total_genes_expected", "consistency", "strand_consistency", "p_value",
        "goi_overlap", "is_goi_anchor", "high_flanking_count", "region_class",
        "goi_missing", "confidence", "selection_reason",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def test_empty_scores_produces_header_only_gff(tmp_path: Path, monkeypatch) -> None:
    """If no block carries goi_missing_but_strong_synteny, emit a stub GFF and
    NEVER touch samtools / miniprot. This is the common case in healthy runs."""
    scores = tmp_path / "GenomeA.scores.tsv"
    _write_scores(scores, [
        {"region_rank": "1", "region_name": "r1", "species": "X", "chrom": "chr1",
         "start": "100", "end": "500", "region_class": "goi_anchor", "goi_missing": "False",
         "confidence": "HIGH"},
        {"region_rank": "2", "region_name": "r2", "species": "X", "chrom": "chr2",
         "start": "1000", "end": "1500", "region_class": "goi_overlap", "goi_missing": "False",
         "confidence": "MEDIUM"},
    ])
    query = tmp_path / "q.faa"
    query.write_text(">myGOI\nMKRLILILV\n")
    out_gff = tmp_path / "out.gff"

    def _fail_if_called(*a, **kw):
        raise AssertionError("samtools / miniprot must not be invoked when no block needs rescue")

    monkeypatch.setattr(rss, "_read_fasta_subseq", _fail_if_called)
    monkeypatch.setattr(rss, "_run_miniprot_relaxed", _fail_if_called)

    monkeypatch.setattr(sys, "argv", [
        "rescue_strong_synteny.py",
        "--scores", str(scores),
        "--target_genome", str(tmp_path / "no_such.fna"),
        "--query", str(query),
        "--genome_name", "GenomeA",
        "--output", str(out_gff),
    ])
    rss.main()

    text = out_gff.read_text()
    assert text.startswith("##gff-version 3")
    assert "0 blocks to rescue" in text
    # No data rows.
    assert not any(line for line in text.splitlines()
                   if line and not line.startswith(("#", "##")))


def test_parse_miniprot_gff_yields_mrna_plus_cds_blocks() -> None:
    """The parser groups CDS rows under the preceding mRNA. Two mRNAs → two tuples."""
    fixture = """##gff-version 3
##PAF\tdummy
window\tminiprot\tmRNA\t10\t90\t100\t+\t.\tID=MP1;Rank=1;Identity=0.812
window\tminiprot\tCDS\t10\t40\t.\t+\t0\tParent=MP1
window\tminiprot\tCDS\t60\t90\t.\t+\t0\tParent=MP1
window\tminiprot\tmRNA\t200\t300\t50\t-\t.\tID=MP2;Rank=2;Identity=0.555
window\tminiprot\tCDS\t200\t300\t.\t-\t0\tParent=MP2
"""
    pairs = list(rss._parse_miniprot_gff(fixture))
    assert len(pairs) == 2
    mrna1, cds1 = pairs[0]
    assert mrna1[2] == "mRNA"
    assert [c[2] for c in cds1] == ["CDS", "CDS"]
    mrna2, cds2 = pairs[1]
    assert mrna2[6] == "-"
    assert len(cds2) == 1


def test_build_rescue_gff_rows_maps_window_coords_to_genome() -> None:
    """Window-relative coords (mRNA at win 10..90) at window_start=5000 → genome
    coords 5009..5089. Attributes tagged for LOW-confidence rescue GOI."""
    mrna = [
        "win", "miniprot", "mRNA", "10", "90", "100", "+", ".",
        "ID=MP1;Rank=1;Identity=0.812",
    ]
    cds_rows = [
        ["win", "miniprot", "CDS", "10", "40", ".", "+", "0", "Parent=MP1"],
        ["win", "miniprot", "CDS", "60", "90", ".", "+", "0", "Parent=MP1"],
    ]
    rows = rss._build_rescue_gff_rows(
        mrna, cds_rows, chrom="chr3", window_start=5000,
        query_id="MyGOI", parent_id="GOI_rescue_GenomeA_b1_1",
    )
    # gene + mRNA + 2 CDS = 4 rows.
    assert len(rows) == 4
    gene_row = rows[0].split("\t")
    assert gene_row[0] == "chr3"
    assert gene_row[2] == "gene"
    assert gene_row[3:5] == ["5009", "5089"]
    assert gene_row[6] == "+"
    # Critical attribute tags so generate_report.py classifies the row correctly:
    assert "SynVoyRole=goi" in gene_row[8]
    assert "Confidence=LOW" in gene_row[8]
    assert "EvidenceType=relaxed_miniprot_rescue" in gene_row[8]
    assert "GOIClass=relaxed_rescue" in gene_row[8]

    # CDS coords are window-relative too — 10..40 → 5009..5039.
    cds_split = rows[2].split("\t")
    assert cds_split[2] == "CDS"
    assert cds_split[3:5] == ["5009", "5039"]


def test_identity_parse_handles_bad_input() -> None:
    """Bad Identity attribute (or missing) shouldn't crash the row builder."""
    assert rss._identity_from_mrna_attrs("ID=foo;Identity=") == 0.0
    assert rss._identity_from_mrna_attrs("ID=foo;Identity=not_a_number") == 0.0
    assert rss._identity_from_mrna_attrs("ID=foo;Identity=0.91") == 0.91
    assert rss._identity_from_mrna_attrs("ID=foo") == 0.0
