#!/usr/bin/env python3
"""
§1n — Two-path home-locus establishment (name-lookup path).

When the GOI is a *named, annotated* gene in the home genome, re-discovering its
location by noisy self-alignment (`locate_gene` + `split_loci`) is solving a problem
you don't have — and for gene families it actively misfires (decorin's LRR cross-hits
outrank the true SLRP locus; the OMD/OGN/ASPN cluster ranks 6th and is dropped at
`max_loci=5`; the target search then anchors on a *paralog's* neighbourhood).

This front-end instead:
  1. Looks the GOI symbol up in the home GFF → a SINGLE, unambiguous annotated locus.
  2. Runs a query↔gene **consistency check** (Smith-Waterman of the query protein vs the
     looked-up gene's translated CDS) so a name typo / wrong species / wrong symbol
     fails loudly instead of silently anchoring the wrong gene.
  3. Emits one locus BED (drop-in for `split_loci` output → `EXTRACT_FLANKING`).

On any miss (no symbol given, symbol absent from the GFF, or consistency below the
threshold) it writes an EMPTY BED and a status so the caller falls back to the existing
alignment-locate path — SynVoy's real domain for genuinely novel/divergent genes.

Output BED columns match `split_loci.py` / `locate_gene.py`:
    chrom  start(0-based)  end  name  evalue  strand  bitscore
"""
import argparse
import sys
from typing import Dict, List, Optional, Tuple

try:
    from sequence_utils import (
        parse_gff,
        parse_fasta,
        translate,
        reverse_complement,
        sw_align,
        setup_logging,
        write_json,
    )
except ImportError:
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from sequence_utils import (  # type: ignore
        parse_gff,
        parse_fasta,
        translate,
        reverse_complement,
        sw_align,
        setup_logging,
        write_json,
    )

logger = setup_logging(name="resolve_home_locus")

# GFF attribute keys that may carry a human gene symbol, in priority order.
SYMBOL_ATTR_KEYS = ("Name", "gene", "gene_name", "gene_id", "locus_tag", "Alias")


def _attr_symbols(attrs: Dict[str, str]) -> List[str]:
    """All candidate symbol strings on a feature (Alias may be comma-separated)."""
    out: List[str] = []
    for key in SYMBOL_ATTR_KEYS:
        val = attrs.get(key)
        if not val:
            continue
        # Alias / Name can be a comma list (e.g. "DCN,PG40,SLRR1B").
        out.extend(p.strip() for p in val.split(",") if p.strip())
    return out


def find_gene_features(gff_path: str, symbol: str) -> List[Dict]:
    """Return all `gene` features whose symbol attributes match `symbol` (case-insensitive)."""
    target = symbol.strip().lower()
    matches = []
    for feat in parse_gff(gff_path, feature_types=["gene"]):
        syms = {s.lower() for s in _attr_symbols(feat["attributes"])}
        if target in syms:
            matches.append(feat)
    return matches


def _gene_id(feat: Dict) -> Optional[str]:
    return feat["attributes"].get("ID") or feat["attributes"].get("gene_id")


def collect_cds_by_gene(gff_path: str, gene_ids: set) -> Dict[str, List[Dict]]:
    """
    Single forward pass: map each matched gene_id -> its CDS features (best transcript).

    Handles the standard gene→mRNA→CDS hierarchy (CDS.Parent = transcript,
    transcript.Parent = gene) and the flatter gene→CDS case (CDS.Parent = gene).
    Picks, per gene, the transcript with the largest total CDS length.
    """
    tid_to_gene: Dict[str, str] = {}
    # gene_id -> transcript_id -> list of CDS records
    cds: Dict[str, Dict[str, List[Dict]]] = {gid: {} for gid in gene_ids}

    for feat in parse_gff(gff_path, feature_types=["mRNA", "transcript", "mrna", "CDS"]):
        attrs = feat["attributes"]
        ftype = feat["type"].lower()
        if ftype in ("mrna", "transcript"):
            parent = attrs.get("Parent")
            tid = attrs.get("ID")
            if parent in gene_ids and tid:
                tid_to_gene[tid] = parent
        elif feat["type"] == "CDS":
            parent = attrs.get("Parent")
            if not parent:
                continue
            # CDS.Parent may be a transcript (normal) or the gene itself (flat GFF).
            for p in parent.split(","):
                gid = tid_to_gene.get(p) or (p if p in gene_ids else None)
                if gid is None:
                    continue
                key = p
                cds[gid].setdefault(key, []).append(feat)

    # Keep, per gene, only the longest transcript's CDS list.
    best: Dict[str, List[Dict]] = {}
    for gid, by_tx in cds.items():
        if not by_tx:
            continue
        best_tx = max(
            by_tx.values(),
            key=lambda recs: sum(r["end"] - r["start"] + 1 for r in recs),
        )
        best[gid] = sorted(best_tx, key=lambda r: r["start"])
    return best


def translate_cds(cds_recs: List[Dict], genome: Dict[str, str], strand: str) -> str:
    """Concatenate CDS exons (genome-ordered), reverse-complement for '-' strand, translate."""
    if not cds_recs:
        return ""
    chrom = cds_recs[0]["seqid"]
    if chrom not in genome:
        return ""
    seq_parts = []
    # GFF is 1-based inclusive; python slice is 0-based half-open.
    for rec in sorted(cds_recs, key=lambda r: r["start"]):
        seq_parts.append(genome[chrom][rec["start"] - 1 : rec["end"]])
    dna = "".join(seq_parts)
    if strand == "-":
        dna = reverse_complement(dna)
    return translate(dna).rstrip("*")


def load_target_chroms(genome_path: str, chroms: set) -> Dict[str, str]:
    """Stream the genome FASTA, keeping only the sequences we need (avoids a 3 GB load)."""
    want = set(chroms)
    out: Dict[str, str] = {}
    if not want:
        return out
    for _header, clean_id, seq in parse_fasta(genome_path):
        if clean_id in want:
            out[clean_id] = seq
            if len(out) == len(want):
                break
    return out


def read_first_protein(faa_path: str) -> str:
    """Read the first sequence from a FASTA (the query GOI protein)."""
    seq = []
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    return "".join(seq)


def pick_best_gene(
    matches: List[Dict], query_seq: str, gff_path: str, genome: Dict[str, str]
) -> Tuple[Optional[Dict], Optional[List[Dict]], float]:
    """
    Among same-symbol gene features, pick the one whose translated protein best aligns
    to the query. Returns (gene_feature, its_cds_records, percent_identity).
    """
    gene_ids = {gid for f in matches if (gid := _gene_id(f))}
    cds_by_gene = collect_cds_by_gene(gff_path, gene_ids)

    best_feat = None
    best_cds = None
    best_ident = -1.0
    for feat in matches:
        gid = _gene_id(feat)
        cds_recs = cds_by_gene.get(gid, [])
        protein = translate_cds(cds_recs, genome, feat["strand"]) if cds_recs else ""
        if not protein:
            # No CDS / unreadable chrom: can't consistency-check this candidate.
            ident = 0.0
        else:
            _score, ident = sw_align(query_seq, protein)
        if ident > best_ident:
            best_feat, best_cds, best_ident = feat, cds_recs, ident
    return best_feat, best_cds, best_ident


def write_empty(out_bed: str, out_status: str, status: str, **extra) -> None:
    open(out_bed, "w").close()
    payload = {"status": status, "matched": False}
    payload.update(extra)
    write_json(out_status, payload)


def main() -> int:
    ap = argparse.ArgumentParser(description="§1n name-lookup home-locus resolver")
    ap.add_argument("--query", required=True, help="GOI protein FASTA (for consistency check)")
    ap.add_argument("--gene_symbol", default="", help="GOI gene symbol (e.g. DCN). Empty => fall back.")
    ap.add_argument("--home_gff", required=True, help="Home genome GFF (annotated)")
    ap.add_argument("--home_genome", required=True, help="Home genome FASTA (for CDS translation)")
    ap.add_argument("--out_bed", required=True, help="Output single-locus BED")
    ap.add_argument("--out_status", required=True, help="Output status JSON")
    ap.add_argument("--min_consistency_identity", type=float, default=60.0,
                    help="Min %% identity (query vs looked-up gene protein) to accept the "
                         "name match (default: 60). Below this => loud fallback to alignment-locate.")
    ap.add_argument("--pad", type=int, default=0, help="Pad (bp) added either side of the gene span in the BED.")
    args = ap.parse_args()

    symbol = (args.gene_symbol or "").strip()
    if not symbol or symbol.upper() in ("NO_SYMBOL", "NA", "NONE", "."):
        logger.info("§1n: no GOI gene symbol provided — falling back to alignment-locate.")
        write_empty(args.out_bed, args.out_status, "no_symbol")
        return 0

    if args.home_gff in ("NO_GFF", "") or not _exists(args.home_gff):
        logger.info("§1n: no usable home GFF — falling back to alignment-locate.")
        write_empty(args.out_bed, args.out_status, "no_gff", symbol=symbol)
        return 0

    matches = find_gene_features(args.home_gff, symbol)
    if not matches:
        logger.info("§1n: symbol '%s' not found in home GFF — falling back to alignment-locate.", symbol)
        write_empty(args.out_bed, args.out_status, "not_in_gff", symbol=symbol)
        return 0

    query_seq = read_first_protein(args.query)
    if not query_seq:
        logger.warning("§1n: empty query protein in %s — falling back.", args.query)
        write_empty(args.out_bed, args.out_status, "empty_query", symbol=symbol)
        return 0

    # Load only the chromosome(s) carrying a same-symbol gene (not the whole 3 GB genome).
    cand_chroms = {f["seqid"] for f in matches}
    genome = load_target_chroms(args.home_genome, cand_chroms)
    feat, _cds, ident = pick_best_gene(matches, query_seq, args.home_gff, genome)

    if feat is None:
        write_empty(args.out_bed, args.out_status, "no_cds", symbol=symbol)
        return 0

    chrom = feat["seqid"]
    g_start, g_end, strand = feat["start"], feat["end"], feat["strand"]
    gene_id = _gene_id(feat) or symbol

    if ident < args.min_consistency_identity:
        logger.warning(
            "§1n: consistency check FAILED for symbol '%s' (%s %s:%d-%d): query vs gene "
            "protein identity %.1f%% < %.1f%% threshold. Likely a name typo, wrong species, "
            "or a paralog symbol — NOT anchoring this locus. Falling back to alignment-locate.",
            symbol, gene_id, chrom, g_start, g_end, ident, args.min_consistency_identity,
        )
        write_empty(args.out_bed, args.out_status, "failed_consistency",
                    symbol=symbol, gene_id=gene_id,
                    chrom=chrom, start=g_start, end=g_end,
                    consistency_identity=round(ident, 2),
                    min_consistency_identity=args.min_consistency_identity,
                    n_symbol_matches=len(matches))
        return 0

    # 0-based half-open BED, padded, clamped to chromosome bounds.
    chrom_len = len(genome.get(chrom, ""))
    bed_start = max(0, (g_start - 1) - args.pad)
    bed_end = min(chrom_len, g_end + args.pad) if chrom_len else g_end + args.pad
    with open(args.out_bed, "w") as fh:
        # evalue=0.0 (definitive annotation); bitscore column carries the SW score proxy.
        fh.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\t0.0\t{strand}\t9999\n")

    logger.info(
        "§1n: anchored single annotated locus for '%s' → %s %s:%d-%d (consistency %.1f%%, "
        "%d symbol match%s). Bypassing split_loci.",
        symbol, gene_id, chrom, bed_start, bed_end, ident,
        len(matches), "" if len(matches) == 1 else "es",
    )
    write_json(args.out_status, {
        "status": "matched", "matched": True, "symbol": symbol, "gene_id": gene_id,
        "chrom": chrom, "start": bed_start, "end": bed_end, "strand": strand,
        "consistency_identity": round(ident, 2),
        "min_consistency_identity": args.min_consistency_identity,
        "n_symbol_matches": len(matches),
        "padded_bp": args.pad,
    })
    return 0


def _exists(path: str) -> bool:
    import os
    return os.path.exists(path)


if __name__ == "__main__":
    sys.exit(main())
