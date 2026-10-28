"""A limit the user sets on anchor places keeps the GOI proxies (fix-phase item 15, 2026-10-28).

With --anchor_max_places N each flanking gene seeds blocks only at its N best places. A GOI
proxy, the GOI's own home gene among the anchors, is a flanking query too, so the user's limit
cut its places as well: a copy of the GOI at the proxy's (N+1)-th place seeded no block. The
genome axis's rebuild (stage 5a) already exempts the proxies (limit_anchor_places, keep=); the
user path in process_single_genome now passes the same set. Proxies occur in 0 of 876 archived
search logs, and the changed call runs only with a number the user sets (the default 'auto' limits
places only in the genome axis's rebuild, which already keeps them), so --modules legacy and auto
runs cannot change.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def _hit(query, chrom, start, bits):
    return {"query": query, "chrom": chrom, "start": start, "end": start + 900, "strand": "+",
            "pident": 60.0, "alnlen": 300, "evalue": 1e-20, "bits": float(bits),
            "qstart": 1, "qend": 300}


def _genome():
    """chr1: the neighbourhood (gene-A, gene-B, and the proxy gene-P next to them); chr9-chr12:
    four further places of gene-P, then its 6th place, the weakest, on chr14. The flanking gene
    gene-A has seven further places (chr20-chr26), so a limit of 5 cuts it to five."""
    hits = [_hit("gene-A", "chr1", 100_000, 500), _hit("gene-B", "chr1", 110_000, 500),
            _hit("gene-P", "chr1", 120_000, 600)]
    hits += [_hit("gene-A", f"chr{c}", 50_000, 120 - c) for c in range(20, 27)]
    hits += [_hit("gene-P", f"chr{c}", 50_000, 80 - c) for c in range(9, 13)]
    hits.append(_hit("gene-P", "chr14", 50_000, 20))
    return hits


def _args(**kw):
    a = argparse.Namespace(max_intron=20_000, cluster_distance=50_000, region_padding=30_000,
                           min_block_genes=2, max_blocks_per_genome=50, genome_axis=False,
                           anchor_max_places=5, genome_ctx={})
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def _opts(proxies):
    proxies = set(proxies)

    def seeded(b):
        return bool(set(b.get("genes", [])) & proxies)

    return dict(home_rank=None, bridge={}, gene_gap=None, goi_seeded=seeded,
                proxy_parents=proxies, gw_any_genome=False)


def _goi_seeded_chroms(seeds, hits):
    opts = _opts({"gene-P"})
    blocks, _values = isr.genome_axis_blocks(seeds, hits, [], _args(), "g", **opts)
    return sorted({b["chrom"] for b in blocks if opts["goi_seeded"](b)}, key=lambda c: int(c[3:]))


def _user_path_call() -> str:
    """The call by which process_single_genome limits the seeds to a user-set number of places."""
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    body = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    m = re.search(r"seed_hits = (limit_anchor_places\(flanking_seed_hits,[^)]*\))", body)
    assert m, "no user-path call of limit_anchor_places"
    return m.group(1)


def test_a_proxy_at_its_sixth_place_seeds_a_block_under_a_limit_of_five():
    """The user path's own call, evaluated on a genome where the proxy's 6th place holds a copy
    of the GOI: under --anchor_max_places 5 that place seeds a GOI-seeded block (the parent's
    call, without keep=, gives the first five places only)."""
    hits = _genome()
    five = ["chr1"] + [f"chr{c}" for c in range(9, 13)]
    before = isr.limit_anchor_places(hits, 5)                    # the call before item 15
    after = eval(_user_path_call(), {"limit_anchor_places": isr.limit_anchor_places},
                 {"flanking_seed_hits": hits, "_max_places": 5, "_gene_gap": None,
                  "goi_proxy_flanking_parents": {"gene-P"}})
    assert _goi_seeded_chroms(before, hits) == five              # the 6th place, chr14, cut
    assert _goi_seeded_chroms(after, hits) == five + ["chr14"]   # kept now
    # the flanking genes are limited as before: gene-A keeps its 5 best places of 8
    assert [h for h in after if h["query"] != "gene-P"] == [h for h in before if h["query"] != "gene-P"]
    assert sorted({h["chrom"] for h in after if h["query"] == "gene-A"}, key=lambda c: int(c[3:])) == \
        ["chr1", "chr20", "chr21", "chr22", "chr23"]


def test_the_user_path_passes_the_proxies():
    """process_single_genome limits the seeds of a user-set value with keep= the GOI proxies,
    the set it builds before that call (and passes to the blocks as proxy_parents)."""
    src = (ROOT / "bin" / "iterative_search_runner.py").read_text()
    body = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    call = re.search(r"seed_hits = limit_anchor_places\(flanking_seed_hits, _max_places, _gene_gap,"
                     r"\s*keep=goi_proxy_flanking_parents\)", body)
    assert call, "the user path does not pass keep=goi_proxy_flanking_parents"
    assert body.index("goi_proxy_flanking_parents: set = set()") < call.start()
    assert len(re.findall(r"limit_anchor_places\(", body)) == 1
