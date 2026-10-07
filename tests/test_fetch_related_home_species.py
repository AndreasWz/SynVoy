"""The home species must never be offered as a related genome.

Easy mode takes the home species from the query's database record. UniProt appends
the strain ("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"), an NCBI protein
record names the strain taxon ("Saccharomyces cerevisiae S288C"). The candidate filter
compared that string for equality with the species name of each assembly record
("Saccharomyces cerevisiae"), so it never matched: on the real genus records
(2026-10-06) the home assembly GCF_000146045.2 itself was the fourth candidate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import fetch_related_genomes as fr  # noqa: E402

# Real `xtract -def NA` rows, genus Saccharomyces, 2026-10-06.
ROWS = "\n".join([
    "GCF_000146045.2\tSaccharomyces cerevisiae\treference genome\tComplete Genome\tNA\tNA\t924431\t924431\tNA\tNA",
    "GCF_002079055.1\tSaccharomyces paradoxus\treference genome\tChromosome\tNA\tNA\t903028\t835812\tNA\tNA",
    "GCF_001298625.1\tSaccharomyces eubayanus\treference genome\tChromosome\tNA\tNA\t896107\t197643\tNA\tNA",
    "GCA_947243795.1\tSaccharomyces uvarum\treference genome\tComplete Genome\tNA\tNA\t920918\t920918\tNA\tNA",
])

HOME_NAMES = [
    "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",   # UniProt
    "Saccharomyces cerevisiae S288C",                          # NCBI strain taxon
    "Saccharomyces cerevisiae",
    "saccharomyces cerevisiae",
]


@pytest.mark.parametrize("home", HOME_NAMES)
def test_home_species_is_recognised_in_every_spelling(home):
    assert fr.is_home_species("Saccharomyces cerevisiae", home)
    assert not fr.is_home_species("Saccharomyces paradoxus", home)


def test_a_genus_or_an_empty_name_is_not_the_home_species():
    assert not fr.is_home_species("Saccharomyces", "Saccharomyces cerevisiae S288C")
    assert not fr.is_home_species("", "Saccharomyces cerevisiae")
    assert not fr.is_home_species("Saccharomyces cerevisiae", "")
    # Two unnamed species of one genus stay different species.
    assert not fr.is_home_species("Drosophila sp. 14028", "Drosophila sp. 14030")


@pytest.mark.parametrize("home", HOME_NAMES)
def test_related_search_never_returns_the_home_species(monkeypatch, home):
    monkeypatch.setattr(fr, "run_piped_command", lambda cmds, timeout=None: ROWS)
    monkeypatch.setattr(fr, "enrich_quality_metadata", lambda entry, cache: entry)
    found = fr.get_related_species("4930", 10, exclude_species=home, tax_level="genus")
    species = [a["species"] for a in found]
    assert "Saccharomyces cerevisiae" not in species
    assert species == ["Saccharomyces paradoxus", "Saccharomyces eubayanus", "Saccharomyces uvarum"]


def test_without_an_exclusion_every_species_is_a_candidate(monkeypatch):
    monkeypatch.setattr(fr, "run_piped_command", lambda cmds, timeout=None: ROWS)
    monkeypatch.setattr(fr, "enrich_quality_metadata", lambda entry, cache: entry)
    found = fr.get_related_species("4930", 10, exclude_species=None, tax_level="genus")
    assert len(found) == 4
