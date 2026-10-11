"""Every gene model of a locus, not only the best (bin/annotate_goi_exons.py, 2026-10-09).

miniprot reports several alignments of the query in a region; until now the first was
used and the rest thrown away, so a tandem array gave one model. Now every model that
does not lie on a better one is kept, and where hits of the window search remain
unexplained the region is searched again with the models found so far masked.
"""
import random
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bin"))
sys.path.insert(0, str(REPO / "tests"))

import annotate_goi_exons as ag  # noqa: E402
from test_miniprot_model_translation import _gene, _protein, _random_dna  # noqa: E402

needs_miniprot = pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
AA = "ARNDCQEGHILKFPSTWYV"


def _mutate(seq, share, seed):
    rnd = random.Random(seed)
    out = list(seq)
    for i in rnd.sample(range(1, len(seq)), int(len(seq) * share)):
        out[i] = rnd.choice(AA.replace(out[i], ""))
    return "".join(out)


def _array(proteins, spacer=1500, lead=0):
    """Genes one after another on the plus strand; returns (dna, [(start, end), ...])."""
    genome, spans = _random_dna(lead, 900), []
    for i, prot in enumerate(proteins):
        g, coords, _pieces = _gene(prot, [100, 302], seed=31 + 7 * i)
        spans.append((len(genome) + coords[0][0], len(genome) + coords[-1][1]))
        genome += g + _random_dna(spacer, 500 + i)
    return genome, spans


def _span(model):
    exons, _protein_seq = model
    return min(e["gstart"] for e in exons), max(e["gend"] for e in exons)


@pytest.fixture
def count_runs(monkeypatch):
    calls = []
    real = ag._run_miniprot

    def counting(query_seq, target, *a, **k):
        calls.append(target)
        return real(query_seq, target, *a, **k)
    monkeypatch.setattr(ag, "_run_miniprot", counting)
    return calls


P = _protein(180, seed=21)


# ----------------------------------------------------------------------- real miniprot

@needs_miniprot
def test_three_close_copies_give_three_models(count_runs):
    copies = [P, _mutate(P, 0.05, 1), _mutate(P, 0.10, 2)]
    genome, spans = _array(copies)
    models = ag.miniprot_gene_models(P, genome, "chr", sensitive=True)
    assert [_span(m) for m in models] == spans                 # best first = copy order here
    assert [protein for _ex, protein in models] == copies
    assert all(len(ex) == 3 for ex, _p in models)
    assert len(count_runs) == 1


@needs_miniprot
def test_the_number_of_models_can_be_capped():
    genome, spans = _array([P, _mutate(P, 0.05, 1), _mutate(P, 0.10, 2)])
    assert [_span(m) for m in ag.miniprot_gene_models(P, genome, "chr", sensitive=True, max_models=1)] == spans[:1]
    assert [_span(m) for m in ag.miniprot_gene_models(P, genome, "chr", sensitive=True, max_models=2)] == spans[:2]


@needs_miniprot
def test_best_model_function_is_unchanged():
    genome, spans = _array([_mutate(P, 0.10, 2), P])
    exons, protein = ag.annotate_using_miniprot(P, genome, "chr", sensitive=True)
    assert (exons[0]["gstart"], exons[-1]["gend"]) == spans[1] and protein == P


@needs_miniprot
def test_one_gene_gives_one_model_in_one_run(count_runs):
    genome, spans = _array([P])
    models = ag.miniprot_gene_models(P, genome, "chr", sensitive=True,
                                     evidence=[(spans[0][0] + 5, spans[0][0] + 95)])
    assert [_span(m) for m in models] == spans and len(count_runs) == 1


@needs_miniprot
def test_diverged_neighbour_copy_needs_the_masked_round(count_runs):
    """80 % copy 2 kb downstream, maximum intron = the region: one run reports only the
    best copy. With a hit on the second copy as evidence it is found."""
    copies = [P, _mutate(P, 0.20, 1)]
    genome, spans = _array(copies)
    one_run = ag.miniprot_gene_models(P, genome, "chr", max_intron=len(genome), sensitive=True)
    assert [_span(m) for m in one_run] == spans[:1]
    assert len(count_runs) == 1

    count_runs.clear()
    evidence = [(s + 10, s + 90) for s, _e in spans]
    models = ag.miniprot_gene_models(P, genome, "chr", max_intron=len(genome), sensitive=True,
                                     evidence=evidence)
    assert [_span(m) for m in models] == spans
    assert [protein for _ex, protein in models] == copies
    assert len(count_runs) == 2
    assert count_runs[0] == genome
    first = spans[0]
    assert set(count_runs[1][first[0]:first[0] + 100]) == {"N"}          # exon 1 of copy 1 masked
    assert count_runs[1][first[0] + 100:first[0] + 110] == genome[first[0] + 100:first[0] + 110]  # its intron not


@needs_miniprot
def test_an_unexplained_hit_on_empty_dna_adds_no_model(count_runs):
    genome, spans = _array([P], spacer=6000)
    models = ag.miniprot_gene_models(P, genome, "chr", sensitive=True,
                                     evidence=[(spans[0][0] + 5, spans[0][0] + 95), (5000, 5090)])
    assert [_span(m) for m in models] == spans
    assert len(count_runs) == 2


@needs_miniprot
def test_models_from_hits_are_in_chromosome_coordinates(count_runs):
    copies = [P, _mutate(P, 0.20, 1)]
    genome, spans = _array(copies, lead=40_000)
    hits = [{"query": "GOI_X", "gstart": s + 10, "gend": s + 90} for s, _e in spans]
    models = ag.annotate_models_from_hit_list(hits, P, genome, "chr", gap_search_window=3000)
    assert [_span(m) for m in models] == spans
    assert len(count_runs[0]) == (spans[1][0] + 90 + 3000) - (spans[0][0] + 10 - 3000)   # the region only

    best_exons, best_protein = ag.annotate_exons_from_hit_list(hits, P, genome, "chr", gap_search_window=3000)
    assert (best_exons[0]["gstart"], best_exons[-1]["gend"]) == spans[0] and best_protein == P


@needs_miniprot
def test_best_only_never_runs_twice(count_runs):
    genome, spans = _array([P, _mutate(P, 0.20, 1)])
    hits = [{"query": "GOI_X", "gstart": s + 10, "gend": s + 90} for s, _e in spans]
    models = ag.annotate_models_from_hit_list(hits, P, genome, "chr", gap_search_window=3000, max_models=1)
    assert [_span(m) for m in models] == spans[:1] and len(count_runs) == 1


# ------------------------------------------------------------- the rules, without miniprot

def _model(score, cds, rank=1):
    return {"id": f"MP{rank}", "rank": rank, "score": score, "strand": "+", "identity": 0.9,
            "cds_list": [{"gstart": b, "gend": e, "strand": "+", "phase": 0, "qstart": 1, "qend": 9}
                         for b, e in cds], "protein": "MK"}


@pytest.fixture
def fake_runs(monkeypatch):
    """Feed prepared miniprot results, one list per run; record the sequence of each run."""
    def install(*runs):
        queue = list(runs)
        seen = []

        def fake(_query, target, *_a, **_k):
            seen.append(target)
            return queue.pop(0) if queue else []
        monkeypatch.setattr(ag, "_run_miniprot", fake)
        monkeypatch.setattr(ag, "_build_model", lambda model, q, seq, name, offset: (
            [{"gstart": c["gstart"] + offset, "gend": c["gend"] + offset} for c in model["cds_list"]],
            model["id"]))
        return seen
    return install


DNA = "ACGT" * 500


def test_second_alignment_to_the_same_gene_is_dropped(fake_runs):
    fake_runs([_model(900, [(100, 200), (300, 400)]),
               _model(850, [(100, 200), (320, 400), (600, 620)], rank=2),     # 180 of 200 bp shared
               _model(600, [(1000, 1100), (1200, 1300)], rank=3)])
    models = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True)
    assert [protein for _ex, protein in models] == ["MP1", "MP3"]


def test_a_model_sharing_less_than_half_of_its_coding_sequence_is_kept(fake_runs):
    fake_runs([_model(900, [(100, 200), (300, 400)]),
               _model(700, [(350, 400), (500, 700)], rank=2)])               # 50 of 250 bp shared
    assert len(ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True)) == 2


def test_overlap_is_counted_on_either_strand(fake_runs):
    minus = _model(800, [(100, 200), (300, 400)], rank=2)
    minus["strand"] = "-"
    fake_runs([_model(900, [(100, 200), (300, 400)]), minus])
    assert len(ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True)) == 1


def test_masked_round_applies_the_score_floor_of_the_first_run(fake_runs):
    """A leftover alignment below 30 % of the best score is not a copy of the gene."""
    seen = fake_runs([_model(1000, [(100, 200)])],
                     [_model(290, [(1000, 1100)]), _model(310, [(1500, 1600)], rank=2)],
                     [_model(280, [(1000, 1100)], rank=3)])
    models = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True,
                                     evidence=[(120, 180), (1020, 1080), (1520, 1580)])
    assert [protein for _ex, protein in models] == ["MP1", "MP2"]
    assert len(seen) == 3                    # round 3: the hit at 1020 is still unexplained
    assert seen[1][100:200] == "N" * 100 and seen[1][:100] == DNA[:100]
    assert seen[2][1500:1600] == "N" * 100 and seen[2][100:200] == "N" * 100


def test_outside_the_sensitive_mode_only_near_ties_are_further_models(fake_runs):
    fake_runs([_model(1000, [(100, 200)])], [_model(700, [(1000, 1100)])])
    models = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=False, evidence=[(1020, 1080)])
    assert len(models) == 1


def test_rounds_stop_when_a_round_adds_nothing(fake_runs):
    seen = fake_runs([_model(1000, [(100, 200)])], [], [_model(900, [(1000, 1100)])])
    models = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True, evidence=[(1020, 1080)])
    assert len(models) == 1 and len(seen) == 2


def test_rounds_stop_at_the_cap(fake_runs):
    seen = fake_runs(*[[_model(1000 - i, [(100 + 300 * i, 200 + 300 * i)])] for i in range(6)])
    models = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True, max_models=3,
                                     evidence=[(1900, 1950)])
    assert len(models) == 3 and len(seen) == 3


def test_no_model_in_the_first_run_means_no_further_runs(fake_runs):
    seen = fake_runs([], [_model(900, [(1000, 1100)])])
    assert ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True, evidence=[(1020, 1080)]) == []
    assert len(seen) == 1


def test_region_offset_is_passed_to_the_model(fake_runs):
    fake_runs([_model(1000, [(100, 200)])])
    (exons, _protein_seq), = ag.miniprot_gene_models("MK", DNA, "chr", sensitive=True, region_offset=50_000)
    assert (exons[0]["gstart"], exons[0]["gend"]) == (50_100, 50_200)


def test_hits_outside_the_region_are_not_evidence(fake_runs, monkeypatch):
    got = {}

    def spy(query, seq, name, **kw):
        got.update(kw, length=len(seq))
        return []
    monkeypatch.setattr(ag, "miniprot_gene_models", spy)
    hits = [{"query": "GOI_X", "gstart": 5000, "gend": 5300}, {"query": "GOI_X", "start": 5400, "end": 5600}]
    ag.annotate_models_from_hit_list(hits, "MK", "A" * 20_000, "chr", gap_search_window=1000)
    assert got["region_offset"] == 4000 and got["length"] == 2600
    assert got["evidence"] == [(1000, 1300), (1400, 1600)]
    assert got["sensitive"] is True and got["max_models"] == 10
    ag.annotate_models_from_hit_list(hits, "MK", "A" * 20_000, "chr", max_models=1)
    assert got["evidence"] is None
