"""A genome with any hit of the genome-wide search has a search block, with the genome-wide path
on or off.

process_single_genome ends the search of a genome whose block build gives None (no block, the path
off) before it writes the genome-end row of target_policy.tsv; with the path such a genome is
searched on. That branch is reached by no genome that has a hit: the seeds are the flanking hits or,
without any, every hit; identify_synteny_blocks makes at least one block of any seed; the
min_block_genes filter falls back to the unfiltered blocks; and the dispersed-GOI rescue only adds
seeds. A genome without a hit ends before the blocks, at the same place with and without the path.
So the limit counters and rescue_dropped_outside_zones of target_policy.tsv do not depend on the
path through that branch. All synthetic."""
import argparse
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402

PITCH = 3000


def _args(**kw):
    a = argparse.Namespace(max_intron=20_000, cluster_distance=7 * PITCH, region_padding=7 * PITCH,
                           min_block_genes=2, max_blocks_per_genome=80, genome_axis=False,
                           anchor_max_places=None, genome_ctx={}, goi_margin=2 * PITCH,
                           gap_search_window=50_000, dispersed_envelope_margin=9 * PITCH)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def _opts(gw_any_genome):
    return dict(home_rank=None, bridge={}, gene_gap=None, proxy_parents=set(),
                goi_seeded=lambda b: any(isr.is_goi_query_id(g) for g in b.get("genes", [])),
                gw_any_genome=gw_any_genome)


def _hit(query, chrom, start, pident=80.0, evalue=1e-50, bits=500.0):
    return {"query": query, "chrom": chrom, "start": start, "end": start + 900, "strand": "+",
            "pident": pident, "evalue": evalue, "alnlen": 900, "bits": bits, "qstart": 1, "qend": 300}


def _random_hits(rng):
    hits = [_hit(f"gene-{rng.choice('ABCDEFG')}", rng.choice(["c1", "c2", "c3", "c4"]),
                 rng.randint(0, 5_000_000)) for _ in range(rng.randint(0, 6))]
    for _ in range(rng.randint(0, 4)):
        weak = rng.random() < 0.6
        hits.append(_hit("GOI_q", rng.choice(["c1", "c2", "c5"]), rng.randint(0, 5_000_000),
                         pident=20.0 if weak else 70.0, evalue=1e-3 if weak else 1e-100,
                         bits=40.0 if weak else 400.0))
    hits.sort(key=lambda h: (h["query"], h["chrom"], h["start"], h["end"], h["strand"], -h["bits"]))
    return isr.resolve_anchor_hits(hits)                # one anchor per place, as the search does


def test_no_block_only_without_a_seed():
    assert isr.build_genome_blocks([], _args(), "g", **_opts(False)) == (None, 0)
    assert isr.build_genome_blocks([], _args(), "g", **_opts(True)) == ([], 0)


@pytest.mark.parametrize("min_block_genes, axis", [(2, False), (5, False), (2, True)])
def test_every_genome_with_a_hit_has_a_block_without_the_path(min_block_genes, axis):
    rng = random.Random(7)
    seen = 0
    for _ in range(400):
        hits = _random_hits(rng)
        if not hits:
            continue
        seen += 1
        args = _args(min_block_genes=min_block_genes, genome_axis=axis,
                     max_blocks_per_genome=rng.choice([1, 80]))
        flanking = [h for h in hits if not isr.is_goi_query_id(h["query"])]
        seeds = flanking if flanking else hits             # as process_single_genome seeds them
        select = isr.rescue_selector(hits, flanking, args, "g") if flanking else None
        blocks, _values, _rescue = isr.genome_axis_blocks(seeds, flanking, args, "g",
                                                         select_rescue=select, **_opts(False))
        assert blocks, hits
    assert seen > 300
