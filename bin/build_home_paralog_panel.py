#!/usr/bin/env python3
"""Build the home-paralog panel for locus ownership (docs/TODO.md §1m).

SynVoy splits the home GOI hits into several home loci (split_loci.py): for a
gene family like the SLRPs, that means one locus per paralog cluster — e.g. the
decorin query (P07585) lands on chr12 (DCN, locus_1) AND on the chr9 OMD/OGN/
ASPN cluster (locus_4). Each home locus runs its own independent iterative
search in other genomes, and because flanking genes can cross-match between
sister paralogs, the *wrong* locus can end up modelling a gene (the chr9 OMD
locus modelled real decorin on cow BTA5, while the chr12 decorin locus produced
nothing there — see docs/TODO.md §1m diagnosis).

This script builds a small reference panel: the home GOI-paralog protein at each
home locus, labelled ``>locus_<id>|<gene-name>``. assign_locus_ownership (which
reuses reciprocal_best_paralog_check.py) then Smith-Waterman aligns each
recovered GOI target against this panel; the best-matching panel entry tells us
which home locus a target *truly* belongs to, independent of which locus's search
happened to recover it. generate_report.py uses that to re-attribute mis-filed
orthologs and relabel paralogs that were carried under the GOI name.

The panel gene at a locus is the annotated home gene(s) overlapping that locus's
GOI-hit span (the narrow ``locus_<id>.bed`` from split_loci, NOT the wider
flanking window — so we pick the central paralog, not its neighbours). A
tandem cluster locus (chr9 OGN+OMD+ASPN within one split locus) contributes one
panel entry per overlapping gene, all sharing the same locus label.

A sidecar ``<output>.meta.tsv`` records panel_id -> (locus, gene, coords,
is_goi_gene) so the report knows which panel gene is the GOI itself (the one the
query aligns to best) and which are paralogs.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import (  # noqa: E402
    parse_gff,
    parse_bed,
    parse_fasta,
    sw_align,
    setup_logging,
)

logger = setup_logging(name="build_home_paralog_panel")

_LOCUS_ID_RE = re.compile(r"locus_([A-Za-z0-9]+)")


def locus_id_from_path(path: str) -> str:
    """Recover ``locus_3`` from ``.../locus_3.bed`` / ``intermediate/.../locus_3.bed``."""
    base = os.path.basename(path)
    m = _LOCUS_ID_RE.search(base)
    return f"locus_{m.group(1)}" if m else os.path.splitext(base)[0]


def load_gene_index(gff_path: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """chrom -> sorted [(start, end, gene_id)] for every ``gene`` feature."""
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for feat in parse_gff(gff_path, feature_types=["gene"]):
        attrs = feat.get("attributes", {})
        gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
        if not gene_id:
            continue
        s, e = feat["start"], feat["end"]
        if e < s:
            s, e = e, s
        by_chrom[feat["seqid"]].append((s, e, gene_id))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def locus_span(bed_path: str) -> Tuple[str, int, int]:
    """Dominant (chrom, min_start, max_end) of a per-locus GOI-hit BED.

    split_loci emits single-chromosome loci; if a BED somehow mixes chroms we
    take the one with the most rows so a stray off-locus hit can't blow up the
    span.
    """
    rows = parse_bed(bed_path)
    if not rows:
        return "", 0, 0
    by_chrom: Dict[str, List[Dict]] = defaultdict(list)
    for r in rows:
        by_chrom[r["chrom"]].append(r)
    chrom = max(by_chrom, key=lambda c: len(by_chrom[c]))
    starts = [r["start"] for r in by_chrom[chrom]]
    ends = [r["end"] for r in by_chrom[chrom]]
    return chrom, min(starts), max(ends)


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def find_family_paralogs(
    query_seq: str,
    proteome: Dict[str, str],
    already: set,
    min_identity: float,
    max_n: int,
    min_len: int = 40,
    min_score_frac: float = 0.25,
) -> List[Tuple[str, str, float]]:
    """The GOI's home paralog FAMILY by homology: align the query against the whole
    home proteome and return the top genes that are both (a) >= ``min_identity`` and
    (b) cover a meaningful fraction of the query (SW score >= ``min_score_frac`` of the
    query self-score). Genes already in the panel are excluded.

    This is the fix for the decorin/SLRP false-positive: ``split_loci`` spans capture
    whatever the noisy self-alignment found (LRR cross-hits like FLRT/LRRC4), NOT the
    real paralogs (biglycan/asporin/osteomodulin). Sourcing the panel from sequence
    homology puts BGN/ASPN/OMD in it, so the RBH check can flag a recovered ``GOI_DCN``
    call whose best home match is actually biglycan.

    The score-fraction guard is essential: local SW %identity over a SHORT spurious
    alignment can exceed ``min_identity`` for unrelated proteins (a few matching
    residues), so identity alone would pull junk LRR proteins into the panel.
    """
    self_score, _ = sw_align(query_seq, query_seq)
    min_score = max(1.0, min_score_frac * self_score)
    scored: List[Tuple[float, float, str]] = []  # (identity, score, gid)
    for gid, seq in proteome.items():
        if gid in already or not seq or len(seq) < min_len:
            continue
        score, ident = sw_align(query_seq, seq)
        if ident >= min_identity and score >= min_score:
            scored.append((ident, score, gid))
    # Rank by identity (closest paralogs first), then SW score.
    scored.sort(reverse=True)
    return [(gid, proteome[gid], ident) for ident, _score, gid in scored[:max_n]]


def genes_for_locus(
    chrom: str, start: int, end: int,
    gene_index: Dict[str, List[Tuple[int, int, str]]],
    pad: int, fallback_window: int, max_genes: int,
) -> List[str]:
    """Pick the home gene(s) representing this locus.

    Primary: genes whose body overlaps the (padded) locus span, ranked by overlap.
    Fallback: if nothing overlaps (sparse/mispadded annotation), the single nearest
    gene within ``fallback_window`` bp. Returns gene ids (possibly empty).
    """
    genes = gene_index.get(chrom, [])
    if not genes:
        return []
    lo, hi = start - pad, end + pad
    overlapping = []
    for gs, ge, gid in genes:
        ov = _overlap(lo, hi, gs, ge)
        if ov > 0:
            overlapping.append((ov, gid))
    if overlapping:
        overlapping.sort(reverse=True)
        return [gid for _, gid in overlapping[:max_genes]]
    # Fallback: nearest gene within window
    best_gid, best_dist = None, fallback_window + 1
    mid = (start + end) // 2
    for gs, ge, gid in genes:
        if ge < lo:
            dist = lo - ge
        elif gs > hi:
            dist = gs - hi
        else:
            dist = 0
        if dist < best_dist:
            best_dist, best_gid = dist, gid
    return [best_gid] if best_gid else []


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the home-paralog panel for locus ownership")
    ap.add_argument("--home_gff", required=True, help="Home genome GFF (gene coords + names)")
    ap.add_argument("--home_proteome", required=True,
                    help="home_proteome.faa from PREPARE_HOME_PROTEOME (>gene-NAME per gene)")
    ap.add_argument("--locus_beds", required=True, nargs="+",
                    help="Per-locus GOI-hit BEDs from split_loci (locus_*.bed)")
    ap.add_argument("--query", default=None,
                    help="Normalised GOI query FASTA — used to flag which panel gene IS the GOI "
                         "(the panel entry the query aligns to best). Optional.")
    ap.add_argument("--output", required=True, help="Panel FASTA (>locus_<id>|<gene>)")
    ap.add_argument("--pad", type=int, default=2000,
                    help="bp padding added to each locus span when intersecting genes")
    ap.add_argument("--fallback_window", type=int, default=20000,
                    help="bp window for the nearest-gene fallback when nothing overlaps")
    ap.add_argument("--max_genes_per_locus", type=int, default=5,
                    help="Cap panel entries per locus (tandem clusters)")
    ap.add_argument("--paralog_min_identity", type=float, default=28.0,
                    help="Min %%identity (query vs a home protein) to include a gene in the "
                         "homology-based paralog FAMILY panel (default 28). Captures the SLRP "
                         "family (BGN/ASPN/OMD ~38-55%%) that split_loci spans miss.")
    ap.add_argument("--max_family_paralogs", type=int, default=15,
                    help="Cap homology-based family entries added beyond the per-locus genes "
                         "(default 15; ranked by identity → the closest paralogs).")
    args = ap.parse_args()

    gene_index = load_gene_index(args.home_gff)
    proteome = {clean_id: seq for _h, clean_id, seq in parse_fasta(args.home_proteome)}
    logger.info(
        f"Loaded {sum(len(v) for v in gene_index.values())} home genes across "
        f"{len(gene_index)} contigs; {len(proteome)} home proteins."
    )

    # Load the query once (longest record = the full protein, not an exon fragment).
    query_seq = ""
    if args.query:
        q_records = [(cid, seq) for _h, cid, seq in parse_fasta(args.query) if seq]
        if q_records:
            q_records.sort(key=lambda kv: len(kv[1]), reverse=True)
            query_seq = q_records[0][1]

    # Fail loud on a missing parasail *before* doing the panel work: with a query present we
    # will Smith-Waterman-align (find_family_paralogs + GOI-gene identification), and a silent
    # crash here (swallowed by the module's errorStrategy) leaves the §1m ownership safety net
    # dead with no signal. The classic cause is a `.venv` shadowing the conda env on PATH.
    if query_seq:
        try:
            import parasail  # type: ignore  # noqa: F401
        except ImportError:
            venv = os.environ.get("VIRTUAL_ENV")
            logger.error(
                "parasail is not importable, but the §1m home-paralog panel needs it for "
                "Smith-Waterman alignment. The locus-ownership / paralog safety net cannot run."
                + (f"\n  NOTE: a virtualenv is active (VIRTUAL_ENV={venv}); it is likely shadowing "
                   f"the conda env that has parasail — deactivate it / strip it from PATH." if venv else "")
                + "\n  Fix: conda install -c bioconda parasail-python, "
                  "or --disable_locus_ownership to skip §1m deliberately."
            )
            sys.exit(1)

    # panel_id -> (gene_id, seq, locus_id, chrom, gstart, gend)
    panel: "Dict[str, Tuple[str, str, str, str, int, int]]" = {}
    gene_meta: Dict[str, Tuple[str, int, int]] = {}
    for gid_list in gene_index.values():
        for gs, ge, gid in gid_list:
            gene_meta.setdefault(gid, (None, gs, ge))  # chrom filled below
    # chrom for each gene
    gid_chrom: Dict[str, str] = {}
    for chrom, gid_list in gene_index.items():
        for gs, ge, gid in gid_list:
            gid_chrom.setdefault(gid, chrom)

    for bed in args.locus_beds:
        lid = locus_id_from_path(bed)
        chrom, s, e = locus_span(bed)
        if not chrom:
            logger.warning(f"[{lid}] empty/unreadable BED — skipped.")
            continue
        gids = genes_for_locus(chrom, s, e, gene_index,
                               args.pad, args.fallback_window, args.max_genes_per_locus)
        kept = 0
        for gid in gids:
            seq = proteome.get(gid)
            if not seq:
                continue
            panel_id = f"{lid}|{gid}"
            gchrom = gid_chrom.get(gid, chrom)
            _, gstart, gend = gene_meta.get(gid, (None, s, e))
            panel[panel_id] = (gid, seq, lid, gchrom, gstart, gend)
            kept += 1
        logger.info(f"[{lid}] {chrom}:{s}-{e} -> {kept} panel gene(s): {gids[:kept]}")

    n_locus_entries = len(panel)

    # Homology-based FAMILY augmentation (the decorin/biglycan fix): split_loci spans
    # only capture what the noisy self-alignment found (LRR cross-hits), missing the real
    # paralogs. Add the GOI's closest home homologs from the whole proteome so the RBH
    # check has BGN/ASPN/OMD to flag a mislabeled GOI call against.
    if query_seq:
        already = {gid for _pid, (gid, *_rest) in panel.items()}
        family = find_family_paralogs(
            query_seq, proteome, already,
            args.paralog_min_identity, args.max_family_paralogs,
        )
        for gid, seq, ident in family:
            gchrom = gid_chrom.get(gid, "")
            _, gstart, gend = gene_meta.get(gid, (None, 0, 0))
            panel[f"home_{gid}|{gid}"] = (gid, seq, f"home_{gid}", gchrom, gstart, gend)
        if family:
            logger.info(
                f"Family augmentation: +{len(family)} home paralog(s) by homology "
                f"(>= {args.paralog_min_identity}% id): "
                f"{', '.join(f'{g}({i:.0f}%)' for g, _s, i in family)}"
            )

    if not panel:
        logger.warning("No panel entries built (no annotated home genes overlap any locus "
                       "and no homology-based paralogs found). Writing empty panel — locus "
                       "ownership will be a no-op.")

    # Identify the GOI's own home gene: the panel gene the query aligns to best
    # (now across BOTH locus and family entries — DCN itself wins at ~100%).
    goi_genes: set = set()
    if query_seq and panel:
        best_gene, best_score = None, -1.0
        for _pid, (gid, seq, _lid, _c, _gs, _ge) in panel.items():
            score, _ident = sw_align(query_seq, seq)
            if score > best_score:
                best_score, best_gene = score, gid
        if best_gene is not None:
            goi_genes.add(best_gene)
            logger.info(f"GOI home gene identified as {best_gene} (SW score {best_score:.0f}).")
    logger.info(f"Panel: {n_locus_entries} locus + {len(panel) - n_locus_entries} family = "
                f"{len(panel)} entries.")

    # Write panel FASTA + meta sidecar.
    with open(args.output, "w") as fh:
        for panel_id, (gid, seq, _lid, _c, _gs, _ge) in sorted(panel.items()):
            fh.write(f">{panel_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    meta_path = args.output + ".meta.tsv"
    with open(meta_path, "w") as fh:
        fh.write("panel_id\tlocus_id\tgene\tchrom\tstart\tend\tis_goi_gene\n")
        for panel_id, (gid, _seq, lid, gchrom, gstart, gend) in sorted(panel.items()):
            is_goi = "1" if gid in goi_genes else "0"
            fh.write(f"{panel_id}\t{lid}\t{gid}\t{gchrom}\t{gstart}\t{gend}\t{is_goi}\n")

    logger.info(f"Wrote {len(panel)} panel entries -> {args.output} (+ {os.path.basename(meta_path)}).")


if __name__ == "__main__":
    main()
