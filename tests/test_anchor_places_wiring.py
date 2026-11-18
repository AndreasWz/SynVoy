"""The genome-wide places of each anchor are counted with the gene's own gap limit.

process_single_genome hands synvoy_evalue.places_per_gene every flanking hit and the gap limit of
one occurrence of the gene in this target (_gene_limit: 1.5 x the gene's home span x the expansion,
at least max_intron). A translated search reports one hit per exon, so with a constant limit each
exon of a gene with long introns was a place of its own, and its E_syn term grew with it. Here one
flanking gene spans 300 kb at home and has 8 exon hits 30 kb apart in the target: one place under
its own limit (450 kb), eight under a constant 5 kb or the 20 kb floor (max_intron).
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import iterative_search_runner as isr  # noqa: E402
import synvoy_scale  # noqa: E402


class _Reached(Exception):
    pass


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n")
    return str(path)


def _hits():
    # gene-A: 8 exons 30 kb apart (gaps above the 20 kb floor, far below its own 450 kb limit)
    exons = [{"query": f"gene-A|exon_{i + 1}", "chrom": "contig1", "start": 1_000 + i * 30_000,
              "end": 1_300 + i * 30_000, "identity": 80.0, "evalue": 1e-20, "bits": 90.0}
             for i in range(8)]
    other = [{"query": "gene-B|exon_1", "chrom": "contig1", "start": 300_000, "end": 300_300,
              "identity": 80.0, "evalue": 1e-20, "bits": 90.0}]
    return exons + other


def test_the_places_of_an_anchor_follow_its_own_gap_limit(tmp_path, monkeypatch):
    geometry = synvoy_scale.home_geometry_from_headers([
        "gene-A home=100000-400000",
        "gene-B home=500000-502000",
    ])
    assert synvoy_scale.gene_gap_limit(geometry, "gene-A", 1.0, 20_000) == 450_000
    genome = _write_fasta(tmp_path / "Target.fna", [("contig1", "A" * 320_000)])
    db = _write_fasta(tmp_path / "db.faa", [("gene-A|exon_1", "MKV"), ("gene-B|exon_1", "MKT"),
                                            ("GOI_x", "MKL")])
    (tmp_path / "hits").mkdir()
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (True, "", False))
    monkeypatch.setattr(isr, "parse_hits", lambda *a, **k: [dict(h) for h in _hits()])
    seen = {}
    real_count_limit = isr.count_limit

    def count_limit(args, key, n=1):
        # called right after the genome's context is set up: stop there
        if key == "mmseqs_maxseqs_genome":
            seen["places"] = dict(args.genome_ctx["places"])
            seen["expansion"] = args.expansion
            raise _Reached()
        return real_count_limit(args, key, n)
    monkeypatch.setattr(isr, "count_limit", count_limit)
    args = SimpleNamespace(output_dir=str(tmp_path), min_identity=20.0, min_length=10, evalue=1.0,
                           max_intron=20_000, home_geometry=geometry)
    isr.process_single_genome(genome, db, args, None, "", 1)
    assert seen["expansion"] == pytest.approx(1.0)
    # with places_per_gene(..., lambda g: 5000) or lambda g: max_intron, gene-A has 8 places
    assert seen["places"] == {"gene-A": 1, "gene-B": 1}
