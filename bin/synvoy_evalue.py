#!/usr/bin/env python3
"""E-values of the label rework (import-only helper of iterative_search_runner.py).

Two quantities of tmp/formula_check_2026-10-14/FINAL_for_vscode.md, section 3.2:

E_syn of a search block: the expected number of windows of this length in the genome that
hold, by chance, a run of >= k anchors in home order (either direction), given how many
places each home neighbour has in the genome. A neighbour with m places lies in a random
window of length W with probability 1 - exp(-m W / G); the number present is a Poisson
binomial; for n present anchors in random order the chance of a monotone run >= k is
exact (Robinson-Schensted-Knuth: the longest increasing run is the first row of the shape,
the longest decreasing run its number of rows). E_syn = (G / W) x sum_n P(n) P(run >= k | n).
A block with E_syn <= 0.1 predicted the gene's place; elsewhere the genome is the search
space.

E of a candidate in its place: blastp of the home protein against the candidate's protein
and against each of its search-hit regions (cut at the aligned span, translated in the
hit's frame, split at stop codons), composition-based statistics, an explicit search space
S = |q| x floor(place_bp / 3) (-searchsp; -dbsize gave 3.33 % false calls at a nominal 1 %
on 16,380 window decoys, -searchsp 0.35 %), read at full precision from blastp's JSON. The
place is the smallest pre-specified one that holds the candidate: the slot between the
nearest found home neighbours on either side (widened by 1 gene pitch), else the core
(block +- GOI margin), else the window; each only in a block with E_syn <= 0.1.
"""
import functools
import json
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

SIGNIFICANT_BLOCK_ESYN = 0.1    # a block with E_syn at or below this predicted the place
PLACE_MERGE_BP = 5000           # hits of one anchor closer than this are one place

PLACES = ("slot", "core", "window", "genome")


# ------------------------------------------------------------------ E_syn
def _partitions(n: int, maxpart: Optional[int] = None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in _partitions(n - first, first):
            yield (first,) + rest


def _hook_product(shape: Tuple[int, ...]) -> int:
    conj = [sum(1 for r in shape if r > j) for j in range(shape[0])] if shape else []
    prod = 1
    for i, r in enumerate(shape):
        for j in range(r):
            prod *= (r - j - 1) + (conj[j] - i - 1) + 1
    return prod


@functools.lru_cache(maxsize=None)
def p_monotone_run(k: int, n: int) -> float:
    """P(longest increasing or longest decreasing subsequence >= k) for a uniformly random
    permutation of n (exact; a shape lambda occurs (f^lambda)^2 times among the n!)."""
    if k <= 1:
        return 1.0 if n >= 1 else 0.0
    if n < k:
        return 0.0
    fact = math.factorial(n)
    below = 0
    for lam in _partitions(n):
        if lam[0] < k and len(lam) < k:
            f = fact // _hook_product(lam)
            below += f * f
    return 1.0 - below / fact


def poisson_binomial(ps: Sequence[float]) -> List[float]:
    """Distribution of the number of successes of independent trials with these chances."""
    dist = [1.0] + [0.0] * len(ps)
    for p in ps:
        for n in range(len(dist) - 1, 0, -1):
            dist[n] = dist[n] * (1 - p) + dist[n - 1] * p
        dist[0] *= (1 - p)
    return dist


def e_syn(k: int, places: Sequence[int], window_len: int, genome_len: int) -> float:
    """Expected number of windows of ``window_len`` in a genome of ``genome_len`` with a
    run of >= k anchors in home order (either direction) by chance; ``places`` = the
    number of genome-wide places of each home neighbour. k <= 0: genome_len / window_len."""
    window_len = max(1, int(window_len))
    genome_len = max(window_len, int(genome_len))
    if k <= 0:
        return genome_len / window_len
    ps = [1.0 - math.exp(-max(0, m) * window_len / genome_len) for m in places]
    dist = poisson_binomial(ps)
    tail = sum(dist[n] * p_monotone_run(k, n) for n in range(len(dist)))
    return (genome_len / window_len) * tail


def places_per_gene(hits: Iterable[Dict], merge_bp: int = PLACE_MERGE_BP) -> Dict[str, int]:
    """{anchor gene: number of places in the genome}: its hits (``query`` base id, ``chrom``,
    ``start``, ``end``) merged when closer than ``merge_bp``."""
    by_gene: Dict[str, List[Tuple[str, int, int]]] = {}
    for h in hits:
        gene = str(h.get("query", "")).split("|")[0]
        if gene:
            by_gene.setdefault(gene, []).append((h["chrom"], int(h["start"]), int(h["end"])))
    out = {}
    for gene, spans in by_gene.items():
        spans.sort()
        n, last_chrom, last_end = 0, None, None
        for chrom, s, e in spans:
            if chrom != last_chrom or s - last_end > merge_bp:
                n += 1
                last_chrom, last_end = chrom, e
            else:
                last_end = max(last_end, e)
        out[gene] = n
    return out


# ------------------------------------------------------------------ places
def ranks_around_goi(genes: Dict[str, Dict[str, int]],
                     goi: Optional[Tuple[int, int]]) -> Dict[str, int]:
    """Home rank of each anchor relative to the gene of interest: -1 for the nearest
    neighbour on its left, -2 for the next, +1 for the nearest on its right. ``genes`` is
    synvoy_scale's {gene: {'start', 'end', ...}} on the home scaffold, ``goi`` the home
    locus span. A gene that overlaps the locus (one that holds the GOI in an intron) has
    no side and no rank. Empty without the locus span."""
    if not goi:
        return {}
    lo, hi = goi
    left = sorted((g for g, c in genes.items() if c['end'] <= lo), key=lambda g: -genes[g]['end'])
    right = sorted((g for g, c in genes.items() if c['start'] >= hi), key=lambda g: genes[g]['start'])
    ranks = {g: -i for i, g in enumerate(left, 1)}
    ranks.update({g: i for i, g in enumerate(right, 1)})
    return ranks


def slot_span(placed: Dict[int, Tuple[int, int]], pitch_bp: float) -> Optional[Tuple[int, int]]:
    """The slot of a block: between the target positions of the nearest found home
    neighbour on the left (rank -1, -2, ...) and on the right (+1, +2, ...), facing edges,
    widened by one gene pitch. None unless both sides have a placed neighbour. ``placed``
    maps home rank (relative to the GOI) to the neighbour's span in the target."""
    lefts = [r for r in placed if r < 0]
    rights = [r for r in placed if r > 0]
    if not lefts or not rights:
        return None
    (ls, le), (rs, re_) = placed[max(lefts)], placed[min(rights)]
    if le <= rs:
        lo, hi = le, rs
    elif re_ <= ls:                 # the neighbourhood is inverted in this genome
        lo, hi = re_, ls
    else:                           # the two neighbours overlap: the slot is their union
        lo, hi = min(ls, rs), max(le, re_)
    pad = int(round(max(0.0, pitch_bp)))
    return max(0, lo - pad), hi + pad


def place_of(span: Tuple[int, int], significant: bool, slot: Optional[Tuple[int, int]],
             core: Tuple[int, int], window: Tuple[int, int], genome_len: int) -> Tuple[str, int]:
    """(place, length in bp): the smallest pre-specified place that holds ``span``; the
    genome when the block is not significant or nothing holds it."""
    s, e = span
    if significant:
        for name, bounds in (("slot", slot), ("core", core), ("window", window)):
            if not bounds:
                continue
            lo, hi = bounds
            if lo <= s and e <= hi:
                return name, max(1, hi - lo)
    return "genome", max(1, int(genome_len))


def search_space(query_len: int, place_bp: int) -> int:
    """S = |q| x floor(place_bp / 3), the residues of one frame of the place."""
    return max(1, int(query_len)) * max(1, int(place_bp) // 3)


# ------------------------------------------------------------------ candidate E-value
def split_at_stops(protein: str, min_len: int = 5) -> List[str]:
    """Pieces of a translation between stop codons ('*'), as separate subjects: joining
    them across stops would let a broken locus look like a gene."""
    return [p for p in protein.split("*") if len(p) >= min_len]


logger = logging.getLogger(__name__)
# Failed blastp calls of this process (review fix RF-1c): each one is logged as a warning,
# because its candidate is judged with E = inf; the runner stops at the start without blastp.
BLASTP_FAILURES = {"count": 0}


def _blastp_failed(reason: str) -> Tuple[float, float]:
    BLASTP_FAILURES["count"] += 1
    logger.warning("blastp for a place E-value failed (%s); the candidate is judged with E = inf "
                   "(failure %d of this process)", reason, BLASTP_FAILURES["count"])
    return float("inf"), 0.0


def _blastp() -> Optional[str]:
    found = shutil.which("blastp")
    if found:
        return found
    sibling = os.path.join(os.path.dirname(sys.executable), "blastp")
    return sibling if os.access(sibling, os.X_OK) else None


def blastp_evalue(query: str, subjects: Sequence[str], searchsp: int) -> Tuple[float, float]:
    """(smallest E, its bit score) of the query against each subject, with
    -comp_based_stats 2 -seg no -word_size 2 -threshold 11 and an explicit -searchsp,
    read from blastp's JSON output (tabular output prints E with 3 decimals between 0.001
    and 0.1: 0.00146 reads as 0.001). (inf, 0.0) when nothing aligns; a missing blastp, a
    nonzero exit, an empty or unparsable output also give (inf, 0.0), each with a warning."""
    blastp = _blastp()
    subjects = [s for s in subjects if s]
    if not query or not subjects:
        return float("inf"), 0.0
    if not blastp:
        return _blastp_failed("blastp not found")
    best = (float("inf"), 0.0)
    with tempfile.TemporaryDirectory(prefix="synvoy_eeff_") as tmp:
        qf = os.path.join(tmp, "q.faa")
        with open(qf, "w") as fh:
            fh.write(f">q\n{query}\n")
        sf = os.path.join(tmp, "s.faa")
        with open(sf, "w") as fh:
            for i, s in enumerate(subjects):
                fh.write(f">s{i}\n{s}\n")
        proc = subprocess.run([blastp, "-query", qf, "-subject", sf, "-comp_based_stats", "2",
                               "-seg", "no", "-word_size", "2", "-threshold", "11",
                               "-searchsp", str(int(searchsp)), "-outfmt", "15"],
                              capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            return _blastp_failed(f"exit {proc.returncode}: {proc.stderr.strip()[-300:]}")
        if not proc.stdout.strip():
            return _blastp_failed("no output")
        try:
            data = json.loads(proc.stdout)
        except ValueError:
            return _blastp_failed("unparsable JSON output")
        for rep in data.get("BlastOutput2", []):
            results = rep.get("report", {}).get("results", {})
            for search in results.get("bl2seq", []) or ([results["search"]] if "search" in results else []):
                for hit in search.get("hits", []):
                    for hsp in hit.get("hsps", []):
                        e, b = float(hsp["evalue"]), float(hsp["bit_score"])
                        if (e, -b) < (best[0], -best[1]):
                            best = (e, b)
    return best
