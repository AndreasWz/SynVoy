"""Easy mode's choice of related genomes (item 57 of the bug list of 16 Nov 2026).

The ranking began with NCBI's status (reference, representative, RefSeq, GenBank) and only then the
assembly level, so within a species a scaffold-level RefSeq assembly beat a chromosome-level GenBank
one, and across species a contig-level RefSeq assembly took a level's place before a chromosome-level
GenBank one. The search reads the target's DNA, and contiguity decides the synteny it can find: the
level comes first now, the status breaks ties. And the quality policy ran after the whole budget was
filled, so a dropped low-quality pick left its place empty; it runs as each candidate is added now,
and the next candidate of the level takes the place. Rows are synthetic (made-up accessions and
species), in the column order of ASSEMBLY_DOCSUM_XTRACT.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import fetch_related_genomes as fr  # noqa: E402


def _row(acc, species, category, level, scaffolds, contigs, scaf_n50, contig_n50):
    return "\t".join([acc, species, category, level, scaffolds, contigs, scaf_n50, contig_n50, "NA", "NA"])


ONE_SPECIES = "\n".join([
    _row("GCF_900000001.1", "Genus alpha", "na", "Scaffold", "900", "4000", "180000", "40000"),
    _row("GCA_900000002.1", "Genus alpha", "na", "Chromosome", "30", "300", "20000000", "900000"),
])

TWO_SPECIES = "\n".join([
    _row("GCF_900000011.1", "Genus beta", "na", "Contig", "NA", "9000", "NA", "60000"),
    _row("GCA_900000012.1", "Genus gamma", "na", "Chromosome", "25", "200", "15000000", "700000"),
    _row("GCA_900000013.1", "Genus delta", "na", "Scaffold", "800", "3000", "400000", "80000"),
])


def _related(monkeypatch, rows, n, mode="hybrid"):
    monkeypatch.setattr(fr, "run_piped_command", lambda cmds, timeout=None: rows)
    monkeypatch.setattr(fr, "enrich_quality_metadata", lambda entry, cache: entry)
    return fr.get_related_species("1", n, exclude_species=None, tax_level="genus", ranking_mode=mode)


@pytest.mark.parametrize("mode", ["hybrid", "counts", "nstats"])
def test_a_chromosome_level_assembly_wins_within_its_species(monkeypatch, mode):
    found = _related(monkeypatch, ONE_SPECIES, 5, mode)
    assert [a["accession"] for a in found] == ["GCA_900000002.1"]


@pytest.mark.parametrize("mode", ["hybrid", "counts", "nstats"])
def test_the_level_orders_the_species_of_a_level(monkeypatch, mode):
    found = _related(monkeypatch, TWO_SPECIES, 5, mode)
    assert [a["species"] for a in found] == ["Genus gamma", "Genus delta", "Genus beta"]
    assert [a["species"] for a in _related(monkeypatch, TWO_SPECIES, 2, mode)] == ["Genus gamma", "Genus delta"]


def _entry(acc, category, level, **counts):
    e = {"accession": acc, "species": "Genus alpha", "category": category, "assembly_status": level}
    e.update(counts)
    return e


def test_the_status_breaks_a_tie_of_level():
    ranked = sorted([_entry("GCA_900000021.1", "na", "Chromosome"),
                     _entry("GCF_900000022.1", "na", "Chromosome"),
                     _entry("GCA_900000023.1", "reference genome", "Chromosome"),
                     _entry("GCF_900000024.1", "representative genome", "Chromosome")],
                    key=lambda e: fr.assembly_rank_tuple(e, "hybrid"))
    assert [e["accession"] for e in ranked] == ["GCA_900000023.1", "GCF_900000024.1", "GCF_900000022.1",
                                                "GCA_900000021.1"]


def test_the_contiguity_tails_of_the_three_modes_are_kept():
    a = _entry("GCF_900000031.1", "na", "Chromosome", chromosome_count=10, scaffold_count=50, contig_count=500,
               scaffold_n50=8e6, contig_n50=1e5)
    b = _entry("GCF_900000032.1", "na", "Chromosome", chromosome_count=12, scaffold_count=40, contig_count=900,
               scaffold_n50=9e6, contig_n50=2e5)
    assert fr.assembly_rank_tuple(a, "hybrid")[2:] < fr.assembly_rank_tuple(b, "hybrid")[2:]   # fewer contigs
    assert fr.assembly_rank_tuple(a, "counts")[2:] < fr.assembly_rank_tuple(b, "counts")[2:]   # fewer chromosomes
    assert fr.assembly_rank_tuple(b, "nstats")[2:] < fr.assembly_rank_tuple(a, "nstats")[2:]   # larger N80/N50
    for mode in ("hybrid", "counts", "nstats"):
        assert fr.assembly_rank_tuple(a, mode)[:2] == (0, 2)            # (level, status) first


def _args(policy="auto"):
    return types.SimpleNamespace(bad_quality_policy=policy, bad_max_contigs=500000, bad_max_scaffolds=500000,
                                 bad_min_n50=5000, bad_quality_timeout=0)


def _candidate(species, level, n50):
    return {"accession": f"GCA_9000001{len(species):02d}.1", "species": species, "category": "na",
            "assembly_status": level, "scaffold_count": 100, "contig_count": 1000,
            "scaffold_n50": n50, "contig_n50": n50}


def test_a_dropped_low_quality_pick_is_refilled_from_its_level(capsys):
    candidates = [_candidate("Genus epsilon", "Scaffold", 800),         # below bad_min_n50: dropped
                  _candidate("Genus zeta", "Scaffold", 90000),
                  _candidate("Genus eta", "Contig", 70000),
                  _candidate("Genus theta", "Contig", 60000)]
    seen = set()
    picked, dropped = fr.fill_level(candidates, 2, seen, _args())
    assert [c["species"] for c in picked] == ["Genus zeta", "Genus eta"] and dropped == 1
    assert "[drop] Excluding low-quality assembly" in capsys.readouterr().out
    assert fr.normalize_species("Genus epsilon") in seen                 # not offered again at the next level
    assert picked[0]["bad_quality"] is False and candidates[0]["bad_quality"] is True


def test_keep_keeps_the_low_quality_pick_with_its_flag():
    candidates = [_candidate("Genus epsilon", "Scaffold", 800), _candidate("Genus zeta", "Scaffold", 90000)]
    picked, dropped = fr.fill_level(candidates, 2, set(), _args("keep"))
    assert [c["species"] for c in picked] == ["Genus epsilon", "Genus zeta"] and dropped == 0
    assert picked[0]["bad_quality"] is True and "best_N50=800 < 5000" in picked[0]["bad_quality_reasons"]


def test_a_species_seen_at_a_closer_level_is_skipped():
    seen = {fr.normalize_species("Genus zeta")}
    picked, _ = fr.fill_level([_candidate("Genus zeta", "Scaffold", 90000),
                               _candidate("Genus eta", "Contig", 70000)], 2, seen, _args())
    assert [c["species"] for c in picked] == ["Genus eta"]
