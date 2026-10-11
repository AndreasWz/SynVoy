#!/usr/bin/env python3
"""Distances of the search in natural units (import-only helper).

A distance that describes a neighbourhood -- how far apart two anchors may lie, how
much DNA around a block is searched -- means something else in every genome. Ten
neighbouring genes span about 20 kb in yeast, 150 kb in the honeybee and 1 Mb in human,
so one number of base pairs is right for one kind of genome at most. The search
therefore states these distances in *gene pitches* and converts them to base pairs for
each home locus and each target genome:

    pitch (home)    the mean step from one anchor of the home neighbourhood to the next:
                    span of the anchors / (number of anchors - 1)
    expansion       how much larger or smaller the same neighbourhood is in a target
                    genome: the length of the longest run of anchors that lie in home
                    order on one target sequence, divided by the length of the same
                    run at home
    pitch (target)  pitch (home) x expansion

The anchors and their home coordinates come from the query set of the search: the
flanking proteins written by extract_flanking_genes.py carry ``home=start-end`` (the
gene) and ``coords=start-end`` (each coding exon) in their FASTA headers, and
PREPARE_INITIAL_DB keeps both. No further input is needed, and a query set without
these fields simply has no geometry: every caller then falls back to fixed base pairs.

Used by iterative_search_runner.py (blocks, windows, bridges), by cluster_grs.py and
rescue_goi_hull.py (which read the values the search wrote into the result GFF), and by
split_loci.py and profile_hits.py (the home step, where the unit is counted in
annotated genes directly).
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple

# Fixed distances of the search before 2026-10-08, in base pairs. They are what a run
# uses when no geometry is available (a query set without home coordinates), and what
# `-profile legacy_distances` sets to reproduce an older run.
LEGACY_BP = {
    'cluster_distance': 150000,
    'region_padding': 150000,
    'padding_min': 50000,
    'padding_max': 200000,
    'gap_search_window': 50000,
    'synteny_bridge_max_gap': 6000000,
    'dispersed_envelope_margin': 200000,
    'max_intron': 20000,
}

# Distances in gene pitches. Calibrated where the pipeline was validated: at the median
# honeybee locus of the family benchmark the anchor pitch is 23 kb, so 7 pitches are the
# former 150 kb and 2 pitches the former 50 kb. The bridge limit was set for one case,
# cow decorin, where 6 Mb are 18 pitches of that neighbourhood; 20 pitches are the length
# of a whole home neighbourhood of 10 genes on each side.
DEFAULT_GENES = {
    'block_gap_genes': 7.0,
    'window_margin_genes': 7.0,
    'goi_margin_genes': 2.0,
    'bridge_max_gap_genes': 20.0,
    'envelope_margin_genes': 9.0,
}

EXPANSION_BOUNDS = (0.2, 8.0)
# The measurement counts only where the gene order is conserved: this many anchors in a row
# along one target sequence, in one direction of home order. Measured on 23 genome-locus
# cases (2026-10-09): 6-19 where the neighbourhood is conserved (yeasts, flies, bees),
# at most 3 where it is not (fly -> mosquito, firefly -> fly / mosquito).
EXPANSION_MIN_CHAIN = 5
EXPANSION_MAX_RANK_STEP = 3
GFF_HEADER = '##synvoy-distances'

_HOME_RE = re.compile(r'(?:^|\s)home=(\d+)-(\d+)')
_GENOME_RE = re.compile(r'(?:^|\s)genome=(\d+)')
_COORDS_RE = re.compile(r'(?:^|\s)coords=(\d+)-(\d+)')


def _base_id(record_id: str) -> str:
    return record_id.split('|', 1)[0]


def _is_goi(record_id: str) -> bool:
    return _base_id(record_id).startswith('GOI_') or bool(re.fullmatch(r'exon_\d+', record_id))


def home_geometry_from_headers(headers: Iterable[str]) -> Optional[Dict[str, Any]]:
    """Geometry of the home neighbourhood from FASTA headers (without the '>').

    Returns None when fewer than two anchors carry coordinates. Otherwise::

        {'order':   [gene id, ...]            anchors in home order
         'genes':   {gene id: {'start', 'end', 'mid', 'span', 'max_intron'}}
         'records': {record id: midpoint}     of every query record with coordinates
                                              (an exon record: that exon; a gene record: the gene)
         'n', 'span', 'pitch', 'max_gap'      of the anchors, in bp
         'home_genome_bp'                     length of the home genome, or None}
    """
    spans: Dict[str, List[int]] = {}
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    records: Dict[str, int] = {}
    genome_bp = None
    for header in headers:
        header = header.strip()
        if not header:
            continue
        record_id = header.split()[0]
        if _is_goi(record_id):
            continue
        if genome_bp is None:
            m = _GENOME_RE.search(header)
            if m:
                genome_bp = int(m.group(1))
        gene = _base_id(record_id)
        m = _HOME_RE.search(header)
        if m:
            spans[gene] = [int(m.group(1)), int(m.group(2))]
            records.setdefault(record_id, (int(m.group(1)) + int(m.group(2))) // 2)
        m = _COORDS_RE.search(header)
        if m and '|exon_' in record_id:
            exons[gene].append((int(m.group(1)), int(m.group(2))))
            records[record_id] = (int(m.group(1)) + int(m.group(2))) // 2

    genes: Dict[str, Dict[str, int]] = {}
    for gene in set(spans) | set(exons):
        parts = sorted(exons.get(gene, []))
        if gene in spans:
            start, end = spans[gene]
        else:
            start, end = parts[0][0], max(e for _s, e in parts)
        if end <= start:
            continue
        introns = [b[0] - a[1] for a, b in zip(parts, parts[1:])]
        genes[gene] = {
            'start': start, 'end': end, 'mid': (start + end) // 2, 'span': end - start,
            'max_intron': max([0] + introns),
        }
    if len(genes) < 2:
        return None
    order = sorted(genes, key=lambda g: (genes[g]['start'], genes[g]['end'], g))
    gaps = [max(0, genes[b]['start'] - genes[a]['end']) for a, b in zip(order, order[1:])]
    span = max(genes[g]['end'] for g in order) - genes[order[0]]['start']
    return {
        'order': order,
        'genes': genes,
        'records': {r: mid for r, mid in records.items() if _base_id(r) in genes},
        'n': len(order),
        'span': span,
        'pitch': span / (len(order) - 1),
        'max_gap': max(gaps),
        'home_genome_bp': genome_bp,
    }


def home_geometry(fasta_path: str) -> Optional[Dict[str, Any]]:
    """Geometry of the home neighbourhood read from the query set of a search."""
    try:
        with open(fasta_path) as fh:
            return home_geometry_from_headers(line[1:] for line in fh if line.startswith('>'))
    except OSError:
        return None


def anchor_run(points: Dict[Any, List[Tuple[int, int, int]]],
               bounds: Tuple[float, float] = EXPANSION_BOUNDS,
               max_home_step: Optional[float] = None) -> Tuple[int, int, int]:
    """The longest run of anchors that lie in home order on one target sequence.

    Returns ``(anchors, length in the target, length at home)``, the two lengths in bp.

    ``points``: target sequence -> [(position in the target, position at home, home rank)],
    one entry per hit of an anchor. A run is a series of hits with rising target position
    whose home ranks all rise or all fall, each step at most EXPANSION_MAX_RANK_STEP ranks
    and with a target / home distance ratio inside ``bounds``. Hits that chance has put
    near each other form runs of two or three anchors; a conserved neighbourhood forms a
    long one. The lengths are sums over the steps of the run, so their quotient lies
    inside ``bounds`` as well. Among runs with the same number of anchors the one that is
    longest at home is taken (then the shorter one in the target), whatever the order of
    the hits. ``max_home_step`` (the largest home distance a step can have) only limits how
    far back the search looks; without it every earlier hit is tried.
    """
    lo, hi = bounds
    reach = hi * max_home_step if max_home_step else None
    best = (0, 0, 0)                       # anchors, length at home, -length in the target
    for key in sorted(points, key=str):
        pts = sorted(points[key])
        if not pts:
            continue
        best = max(best, (1, 0, 0))
        for sign in (1, -1):
            run = [(1, 0, 0)] * len(pts)
            first = 0
            for j, (tj, hj, rj) in enumerate(pts):
                if reach is not None:
                    while tj - pts[first][0] > reach:
                        first += 1
                for i in range(first, j):
                    ti, hi_home, ri = pts[i]
                    step = (rj - ri) * sign
                    if step < 1 or step > EXPANSION_MAX_RANK_STEP or run[i][0] + 1 < run[j][0]:
                        continue
                    d_home = abs(hj - hi_home)
                    if d_home <= 0 or not lo <= (tj - ti) / d_home <= hi:
                        continue
                    cand = (run[i][0] + 1, run[i][1] + d_home, run[i][2] - (tj - ti))
                    if cand > run[j]:
                        run[j] = cand
            best = max(best, max(run))
    return best[0], -best[2], best[1]


def longest_anchor_chain(points: Dict[Any, List[Tuple[int, int, int]]],
                         bounds: Tuple[float, float] = EXPANSION_BOUNDS,
                         max_home_step: Optional[float] = None) -> int:
    """Largest number of anchors that lie in a row, in home order, on one target sequence
    (the first value of ``anchor_run``)."""
    return anchor_run(points, bounds, max_home_step)[0]


def target_expansion(geometry: Optional[Dict[str, Any]], hits: Iterable[Dict[str, Any]],
                     size_ratio: Optional[float] = None,
                     bounds: Tuple[float, float] = EXPANSION_BOUNDS,
                     min_chain: int = EXPANSION_MIN_CHAIN) -> Tuple[float, int, str]:
    """How much larger the home neighbourhood is in one target genome.

    Returns ``(factor, anchors, source)``; ``anchors`` is the number of anchors in the
    longest run found (``anchor_run``), whether or not it was long enough to be used.

    Measured (``source`` = 'anchors'): the length of the longest run of anchors that lie in
    home order on one target sequence, divided by the length of the same run at home. Each
    step of the run goes from one hit to the next and is measured at home between the two
    query records those hits belong to: with exon-level queries between two exons,
    otherwise between two genes, so both lengths are taken between the same things.

    The measurement is used only where the gene order is conserved: the run must hold at
    least ``min_chain`` anchors. In a genome where the order is lost, hits of two or three
    anchors still lie near each other by chance and would give any factor up to the bound.

    Not measured: ``size_ratio``, the length of the target genome divided by the length of
    the home genome, kept inside ``bounds`` (``source`` = 'genome size'); without it the
    factor is 1 (``source`` = 'none').

    What this replaced, and why:
    - Until 2026-10-09 the nearest hits of two anchors were divided by the distance of the
      two GENE midpoints. For genes with several exons the nearest exons are closer than
      the gene midpoints, so the factor came out too small: 0.74 for *D. simulans* at a
      locus whose anchors span the same 100 kb as in *D. melanogaster*.
    - Until 2026-10-09 the measurement was accepted when a quarter of the anchor pairs on
      one sequence lay near each other. Firefly against fruit fly passed with 11 chance
      pairs and a factor of 3.7, where the genome sizes say 0.3; one search window grew to
      megabases and used 10 GB of memory.
    - Until 2026-10-11 the factor was the MEDIAN of the distance ratios of anchor pairs up
      to three steps apart. A neighbourhood rarely grows evenly: most neighbouring genes
      keep their distance and a few large insertions make up the difference, which a
      median does not see. *A. lyrata* at the FT locus of *A. thaliana*: 1.5, where the
      same 17 genes span 387 kb instead of 78 kb (4.9); the bee *Colletes gigas* at the
      melittin locus: 2.0 instead of 4.5. The quotient of two lengths is the quantity the
      distances need: how long the neighbourhood is there.
    """
    def fallback(anchors):
        if size_ratio and size_ratio > 0:
            return max(bounds[0], min(bounds[1], float(size_ratio))), anchors, 'genome size'
        return 1.0, anchors, 'none'

    if not geometry:
        return fallback(0)
    genes = geometry['genes']
    record_mid = geometry.get('records') or {}
    rank = {g: i for i, g in enumerate(geometry['order'])}
    # target sequence -> [(midpoint of the hit, midpoint at home of its query record, rank)]
    points: Dict[Any, List[Tuple[int, int, int]]] = defaultdict(list)
    for h in hits:
        record = str(h.get('query', ''))
        gene = _base_id(record)
        if gene in genes:
            target_mid = (int(h.get('start', 0)) + int(h.get('end', 0))) // 2
            points[h.get('chrom')].append((target_mid, record_mid.get(record, genes[gene]['mid']), rank[gene]))
    full = geometry['order']
    max_home_step = max(
        (max(genes[g]['end'] for g in full[i:i + EXPANSION_MAX_RANK_STEP + 1]) - genes[full[i]]['start']
         for i in range(len(full) - 1)), default=0)
    anchors, target_length, home_length = anchor_run(points, bounds, max_home_step or None)
    if anchors < min_chain or home_length <= 0:
        return fallback(anchors)
    return target_length / home_length, anchors, 'anchors'


def derive_distances(geometry: Optional[Dict[str, Any]], expansion: float = 1.0,
                     genes: Optional[Dict[str, float]] = None,
                     fixed_bp: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """Base-pair distances of one home locus in one target genome.

    ``genes`` are the distances in gene pitches (missing keys: DEFAULT_GENES).
    ``fixed_bp`` holds base-pair values set by the user; a value > 0 is used as it is
    and the pitch is ignored for that distance. Without ``geometry`` every distance not
    fixed by the user takes its former fixed value (LEGACY_BP).

    The block gap is never smaller than 1.25 times the largest gap between two
    neighbouring anchors at home, so that a neighbourhood that is unchanged in the
    target comes out as one block.
    """
    units = dict(DEFAULT_GENES)
    units.update({k: float(v) for k, v in (genes or {}).items() if v is not None})
    fixed = {k: float(v) for k, v in (fixed_bp or {}).items() if v is not None and float(v) > 0}

    def pick(bp_key: str, derived: Optional[float]) -> Tuple[int, str]:
        if bp_key in fixed:
            return int(fixed[bp_key]), 'set in bp'
        if derived is None:
            return int(LEGACY_BP[bp_key]), 'fixed (no home coordinates)'
        return max(1, int(round(derived))), 'derived'

    pitch_t = geometry['pitch'] * expansion if geometry else None
    out: Dict[str, Any] = {
        'expansion': round(float(expansion), 3),
        'home_pitch_bp': int(round(geometry['pitch'])) if geometry else None,
        'home_span_bp': int(geometry['span']) if geometry else None,
        'home_anchors': geometry['n'] if geometry else 0,
        'source': {},
    }
    plan = (
        ('cluster_distance', None if pitch_t is None else
         max(units['block_gap_genes'] * pitch_t, 1.25 * geometry['max_gap'] * expansion)),
        ('region_padding', None if pitch_t is None else units['window_margin_genes'] * pitch_t),
        ('goi_margin', None if pitch_t is None else units['goi_margin_genes'] * pitch_t),
        ('synteny_bridge_max_gap', None if pitch_t is None else units['bridge_max_gap_genes'] * pitch_t),
        ('dispersed_envelope_margin', None if pitch_t is None else units['envelope_margin_genes'] * pitch_t),
    )
    for key, derived in plan:
        bp_key = 'gap_search_window' if key == 'goi_margin' else key
        out[key], out['source'][key] = pick(bp_key, derived)
    return out


def gene_gap_limit(geometry: Optional[Dict[str, Any]], gene: str, expansion: float,
                   floor_bp: int, margin: float = 1.5) -> int:
    """Largest gap between two hits of one gene that still belong to one occurrence.

    A gene's exons lie within the gene, so its own span at home (times the expansion of
    the target, times a margin for longer introns) is the natural limit; ``floor_bp``
    is used for compact genes and for genes without home coordinates.
    """
    if not geometry or gene not in geometry['genes']:
        return int(floor_bp)
    return int(max(floor_bp, margin * geometry['genes'][gene]['span'] * expansion))


def format_gff_header(distances: Dict[str, Any]) -> str:
    """One comment line for the result GFF of a genome, read back by parse_gff_header."""
    return f"{GFF_HEADER} {json.dumps(distances, sort_keys=True, separators=(',', ':'))}"


def parse_gff_header(gff_path: str) -> Optional[Dict[str, Any]]:
    """The distances the search used for this genome, or None (older GFF, no file)."""
    try:
        with open(gff_path) as fh:
            for line in fh:
                if not line.startswith('#'):
                    return None
                if line.startswith(GFF_HEADER):
                    return json.loads(line[len(GFF_HEADER):].strip())
    except (OSError, ValueError):
        return None
    return None


# Two places where the query aligns at home are one locus up to this gap when the home
# genome has no annotation to count genes in (split_loci.py, profile_hits.py).
FALLBACK_LOCUS_GAP_BP = 50000


def load_gene_spans(gff_path: Optional[str], chroms) -> Dict[str, List[Tuple[int, int]]]:
    """Sorted ``(start, end)`` of the annotated genes on the given sequences (0-based
    start). Protein-coding genes only where the annotation states a biotype; every
    ``gene`` row otherwise. Returns {} when the file cannot be used ('NO_GFF', missing)."""
    import os
    if not gff_path or os.path.basename(str(gff_path)) == 'NO_GFF' or not os.path.exists(gff_path):
        return {}
    coding: Dict[str, List[Tuple[int, int]]] = {}
    any_gene: Dict[str, List[Tuple[int, int]]] = {}
    try:
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.split('\t', 9)
                if len(parts) < 9 or parts[2] != 'gene' or parts[0] not in chroms:
                    continue
                span = (int(parts[3]) - 1, int(parts[4]))
                any_gene.setdefault(parts[0], []).append(span)
                if 'gene_biotype=protein_coding' in parts[8] or 'biotype=protein_coding' in parts[8]:
                    coding.setdefault(parts[0], []).append(span)
    except (OSError, ValueError):
        return {}
    genes = coding if coding else any_gene
    return {chrom: sorted(spans) for chrom, spans in genes.items()}


def genes_between(sorted_genes: List[Tuple[int, int]], left_end: int, right_start: int) -> int:
    """Number of annotated genes that lie entirely between two positions.

    ``sorted_genes`` are ``(start, end)`` of the genes of one sequence, sorted."""
    return sum(1 for s, e in sorted_genes if s >= left_end and e <= right_start)
