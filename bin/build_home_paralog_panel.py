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
is_goi_gene, goi_found_by) so the report knows which panel entry is the GOI itself and
which are paralogs. The GOI entry is q's own home gene, found as the search's
specificity check finds it (synvoy_specificity: the gene the home step took, from
goi_info.json, else the home protein >= 95 % identical over >= 90 % of q); when no home
gene is q's own, q itself is added as the GOI entry at the primary home locus (fix-phase
item 17, 2026-10-28). q is the search's q: the GOI's home protein, the full-length GOI
record of goi_exons.faa (--home_goi), with the query as the fallback. Until then the GOI
was the panel gene the query aligned to best, with no minimum: without a home GFF the
GOI's own gene is never in the home proteome (borrow_annotations.write_home_proteome
leaves it out on purpose), so a neighbour became the "GOI" and the true orthologs were
relabelled paralog_not_goi.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import (  # noqa: E402
    parse_gff,
    parse_bed,
    parse_fasta,
    sw_align,
    setup_logging,
)
import synvoy_settings  # noqa: E402
import synvoy_specificity  # noqa: E402

logger = setup_logging(name="build_home_paralog_panel")

_LOCUS_ID_RE = re.compile(r"locus_([A-Za-z0-9]+)")


def locus_id_from_path(path: str) -> str:
    """Recover ``locus_3`` from ``.../locus_3.bed`` / ``intermediate/.../locus_3.bed``."""
    base = os.path.basename(path)
    m = _LOCUS_ID_RE.search(base)
    return f"locus_{m.group(1)}" if m else os.path.splitext(base)[0]


def _bed_index(path: str) -> int:
    """A locus BED's index, as lib/HomeLoci.groovy reads it: the number at the end of the file
    name without its extension; a name without one comes last."""
    m = re.search(r"(\d+)$", os.path.splitext(os.path.basename(path))[0])
    return int(m.group(1)) if m else sys.maxsize


def primary_locus_bed(beds: List[str]) -> Optional[str]:
    """The primary home locus BED, as lib/HomeLoci.groovy picks it: the lowest index, the
    first of equal ones."""
    return min(beds, key=_bed_index) if beds else None


def locus_holding(beds: List[str], chrom: str, start: int, end: int, pad: int) -> Optional[str]:
    """The id of the home locus whose span, padded as genes_for_locus pads it, overlaps a gene
    at chrom:start-end; the lowest index first. None when no locus span holds the gene."""
    for bed in sorted(beds, key=_bed_index):
        b_chrom, b_s, b_e = locus_span(bed)
        if chrom and b_chrom == chrom and _overlap(b_s - pad, b_e + pad, start, end) > 0:
            return locus_id_from_path(bed)
    return None


def panel_name(query_id: str) -> str:
    """q's id as the gene part of a panel id: no '|' (the report splits a panel id at its
    first '|', and a UniProt-style id holds two) and no whitespace (a FASTA header ends at
    the first blank)."""
    return re.sub(r"[|\s]+", "_", query_id.strip()) or "query"


def home_goi_record(path: Optional[str]) -> Optional[Tuple[str, str]]:
    """The GOI's home protein in goi_exons.faa of ANNOTATE_GOI, picked as the search picks its
    home GOI query (iterative_search_runner.select_goi_queries): a full-length GOI record, a
    GOI_ id that is no exon (|exon_), fragment (|frag_) or GOI_copy_ record (a query id can hold
    '|': GOI_lcl|ORF1), an id without '|' first, then the longest. ANNOTATE_GOI writes one, the
    annotated gene's protein or the home model, before its exons. It is the q of the search's
    specificity check. None without a file or such a record."""
    if not path or not os.path.isfile(path):
        return None
    records = [(cid, seq) for _h, cid, seq in parse_fasta(path)
               if seq and cid.startswith("GOI_") and not cid.startswith("GOI_copy_")
               and "|exon_" not in cid and "|frag_" not in cid]
    if not records:
        return None
    records.sort(key=lambda kv: ("|" in kv[0], -len(kv[1]), kv[0]))
    return records[0]


def _goi_info_gives_no_gene(goi_info: Optional[str], info_gene: Optional[str]) -> str:
    """Why goi_info.json names no gene of q's own for the search (for the log)."""
    if info_gene:
        return f"goi_info.json names {info_gene}, which the home proteome does not hold"
    try:
        with open(goi_info) as fh:
            raw = json.load(fh)
    except (OSError, TypeError, ValueError):
        raw = {}
    gene = raw.get("gene_id") if isinstance(raw, dict) else None
    if not gene:
        return "goi_info.json names no home gene"
    if raw.get("home_gene_ambiguous"):
        return f"goi_info.json names {gene}, which the home step flagged as ambiguous"
    return f"goi_info.json names {gene}, found by {raw.get('found_by')}, which the search does not take"


def goi_entry(proteome: Dict[str, str], query_id: str, query_seq: str,
              goi_info: Optional[str]) -> Tuple[Optional[str], str, str]:
    """q's own home gene, found as the search's specificity check finds it
    (synvoy_specificity.HomeProteome.own_gene): the gene the home step took as q's
    (own_gene_id_from_goi_info: found by name or by position, not flagged ambiguous) when
    the proteome holds it, else the home protein >= 95 % identical over >= 90 % of q.
    Returns (gene id or None, how: 'goi_info' / 'identity_rule' / 'query', a note for the
    log); None = no home gene is q's own, so q itself is the GOI entry."""
    info_gene = (synvoy_specificity.own_gene_id_from_goi_info(goi_info)
                 if goi_info and os.path.exists(goi_info) else None)
    home = synvoy_specificity.HomeProteome(
        [(gid, seq.upper().rstrip("*")) for gid, seq in proteome.items()], {})
    own = home.own_gene(query_id, query_seq, info_gene)
    rule = (f">= {synvoy_specificity.OWN_GENE_MIN_IDENTITY:.0f} % identity over >= "
            f"{synvoy_specificity.OWN_GENE_MIN_COVERAGE * 100:.0f} % of q")
    if own and own == info_gene:
        return own, "goi_info", "q's own home gene, the one the home step took (goi_info.json)"
    if own:
        return own, "identity_rule", f"q's own home gene by the identity rule: {rule}"
    why = _goi_info_gives_no_gene(goi_info, info_gene)
    return None, "query", f"no home gene is q's own: {why}, and no home protein has {rule}"


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
                    help="Normalised GOI query FASTA: the homology family, and q of the GOI entry "
                         "when --home_goi gives no home GOI protein. Optional.")
    ap.add_argument("--goi_info", default=None,
                    help="goi_info.json of ANNOTATE_GOI: the home gene the home step took as q's "
                         "is the panel's GOI entry (as for the search's specificity check); "
                         "without one, the home protein near-identical to q, else q itself")
    ap.add_argument("--home_goi", default=None,
                    help="goi_exons.faa of ANNOTATE_GOI: its full-length GOI record, the GOI's home "
                         "protein, is q for the GOI entry, as in the search's specificity check; "
                         "without it, the query")
    ap.add_argument("--output", required=True, help="Panel FASTA (>locus_<id>|<gene>)")
    ap.add_argument("--pad", type=int, default=2000,
                    help="bp padding added to each locus span when intersecting genes")
    ap.add_argument("--fallback_window", type=int, default=20000,
                    help="bp window for the nearest-gene fallback when nothing overlaps")
    ap.add_argument("--max_genes_per_locus", type=int, default=None,
                    help="Cap panel entries per locus (tandem clusters). Default: the run's "
                         "locus_ownership_max_genes_per_locus (--effective_params), else 5")
    ap.add_argument("--paralog_min_identity", type=float, default=28.0,
                    help="Min %%identity (query vs a home protein) to include a gene in the "
                         "homology-based paralog FAMILY panel (default 28). Captures the SLRP "
                         "family (BGN/ASPN/OMD ~38-55%%) that split_loci spans miss.")
    ap.add_argument("--max_family_paralogs", type=int, default=None,
                    help="Cap homology-based family entries added beyond the per-locus genes "
                         "(ranked by identity → the closest paralogs). Default: the run's "
                         "paralog_panel_max_family (--effective_params), else 15")
    ap.add_argument("--effective_params", default=None,
                    help="effective_params.json of the run (RESOLVE_EFFECTIVE_PARAMS): the two "
                         "caps above come from its settings unless given here; a sentinel or "
                         "a missing file = the defaults")
    args = ap.parse_args()
    settings = synvoy_settings.load_settings(args.effective_params)
    max_genes_per_locus, genes_from = synvoy_settings.pick(
        args.max_genes_per_locus, settings, "locus_ownership_max_genes_per_locus", 5,
        synvoy_settings.whole)
    max_family_paralogs, family_from = synvoy_settings.pick(
        args.max_family_paralogs, settings, "paralog_panel_max_family", 15, synvoy_settings.whole)
    logger.info(f"Panel caps: {max_genes_per_locus} genes per locus ({genes_from}), "
                f"{max_family_paralogs} family paralogs ({family_from}).")

    gene_index = load_gene_index(args.home_gff)
    proteome = {clean_id: seq for _h, clean_id, seq in parse_fasta(args.home_proteome)}
    logger.info(
        f"Loaded {sum(len(v) for v in gene_index.values())} home genes across "
        f"{len(gene_index)} contigs; {len(proteome)} home proteins."
    )

    # Load the query once (longest record = the full protein, not an exon fragment).
    query_id, query_seq = "", ""
    if args.query:
        q_records = [(cid, seq) for _h, cid, seq in parse_fasta(args.query) if seq]
        if q_records:
            q_records.sort(key=lambda kv: len(kv[1]), reverse=True)
            query_id, query_seq = q_records[0]

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
                               args.pad, args.fallback_window, max_genes_per_locus)
        kept, no_protein = [], []
        for gid in gids:
            seq = proteome.get(gid)
            if not seq:
                no_protein.append(gid)
                continue
            panel_id = f"{lid}|{gid}"
            gchrom = gid_chrom.get(gid, chrom)
            _, gstart, gend = gene_meta.get(gid, (None, s, e))
            panel[panel_id] = (gid, seq, lid, gchrom, gstart, gend)
            kept.append(gid)
        logger.info(f"[{lid}] {chrom}:{s}-{e} -> {len(kept)} panel gene(s): {kept}"
                    + (f"; without a protein in the home proteome: {no_protein}" if no_protein else ""))

    n_locus_entries = len(panel)

    # Homology-based FAMILY augmentation (the decorin/biglycan fix): split_loci spans
    # only capture what the noisy self-alignment found (LRR cross-hits), missing the real
    # paralogs. Add the GOI's closest home homologs from the whole proteome so the RBH
    # check has BGN/ASPN/OMD to flag a mislabeled GOI call against.
    if query_seq:
        already = {gid for _pid, (gid, *_rest) in panel.items()}
        family = find_family_paralogs(
            query_seq, proteome, already,
            args.paralog_min_identity, max_family_paralogs,
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

    n_family_entries = len(panel) - n_locus_entries

    # The GOI entry (fix-phase item 17): q's own home gene as the search's specificity check
    # finds it (goi_entry), not the panel gene q aligns to best. Without a home GFF the GOI's
    # own gene is never in this proteome, so the best-aligning gene was a neighbour, and the
    # neighbours then "owned" the true orthologs. q's own gene that the gene cap per locus left
    # out joins the locus whose span holds it; one outside every locus span is at the family
    # augmentation's pseudo-locus home_<gene>, as before, so a call it owns is marked as found
    # from another home locus. q itself, when no home gene is q's own, is added at the primary
    # home locus.
    goi_genes: set = set()
    goi_found_by = ""
    added = ""
    if query_seq:
        # q is the search's q: the GOI's home protein (goi_exons.faa), not the query file. The two
        # are the same protein for a query of the home species; for a query of another species
        # the identity rule on the query names no gene where the search names the home gene.
        home = home_goi_record(args.home_goi)
        q_id, q_seq = home if home else (query_id, query_seq)
        if home:
            q_from = "the GOI's home protein from goi_exons.faa"
        elif args.home_goi and os.path.isfile(args.home_goi):
            q_from = "the query; goi_exons.faa holds no full-length GOI record"
        else:
            q_from = "the query; no home GOI protein given"
        logger.info(f"q for the GOI entry: {q_id} ({len(q_seq)} aa, {q_from})")
        own, goi_found_by, how = goi_entry(proteome, q_id, q_seq, args.goi_info)
        primary = primary_locus_bed(args.locus_beds)
        p_lid = locus_id_from_path(primary) if primary else "locus_1"
        p_chrom, p_s, p_e = locus_span(primary) if primary else ("", 0, 0)
        if own:
            goi_genes.add(own)
            where = ""
            family_id = f"home_{own}|{own}"
            if not any(gid == own and pid != family_id for pid, (gid, *_rest) in panel.items()):
                o_chrom = gid_chrom.get(own, "")
                _, o_s, o_e = gene_meta.get(own, (None, 0, 0))
                o_lid = locus_holding(args.locus_beds, o_chrom, o_s, o_e, args.pad)
                if o_lid:
                    if panel.pop(family_id, None) is not None:
                        n_family_entries -= 1
                    panel[f"{o_lid}|{own}"] = (own, proteome[own], o_lid, o_chrom, o_s, o_e)
                    added = " + q's own gene"
                    where = f"; left out by the gene cap per locus, added at {o_lid}, whose span holds it"
                else:
                    if family_id not in panel:
                        panel[family_id] = (own, proteome[own], f"home_{own}", o_chrom, o_s, o_e)
                        added = " + q's own gene"
                    where = f"; in no locus span, at home_{own}"
            logger.info(f"GOI entry: {own}, {how}{where}.")
        else:
            name = panel_name(q_id)
            while f"{p_lid}|{name}" in panel:
                name += "_query"
            panel[f"{p_lid}|{name}"] = (name, q_seq, p_lid, p_chrom, p_s, p_e)
            goi_genes.add(name)
            added = " + q"
            logger.info(f"GOI entry: q itself ({q_id}), as {p_lid}|{name}; {how}.")
    logger.info(f"Panel: {n_locus_entries} locus + {n_family_entries} family{added} = "
                f"{len(panel)} entries.")
    if len(panel) < 2:
        logger.warning(f"The panel holds {len(panel)} entr{'y' if len(panel) == 1 else 'ies'}: "
                       "locus ownership compares a call with two or more, so it changes nothing "
                       "in this run (no other home gene in a locus span, and no home protein "
                       "close enough to the query).")

    # Write panel FASTA + meta sidecar.
    with open(args.output, "w") as fh:
        for panel_id, (gid, seq, _lid, _c, _gs, _ge) in sorted(panel.items()):
            fh.write(f">{panel_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    meta_path = args.output + ".meta.tsv"
    with open(meta_path, "w") as fh:
        fh.write("panel_id\tlocus_id\tgene\tchrom\tstart\tend\tis_goi_gene\tgoi_found_by\n")
        for panel_id, (gid, _seq, lid, gchrom, gstart, gend) in sorted(panel.items()):
            is_goi = gid in goi_genes
            found_by = goi_found_by if is_goi else ""
            fh.write(f"{panel_id}\t{lid}\t{gid}\t{gchrom}\t{gstart}\t{gend}\t"
                     f"{'1' if is_goi else '0'}\t{found_by}\n")

    logger.info(f"Wrote {len(panel)} panel entries -> {args.output} (+ {os.path.basename(meta_path)}).")


if __name__ == "__main__":
    main()
