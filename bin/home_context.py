#!/usr/bin/env python3
"""The home gene's own neighbourhood and scaffold, for the preset modules (select_modules.py).

  relatives_nearby      protein-coding home genes within +-n_flanking_genes gene ranks of
                        the gene of interest on its scaffold that are FAMILY MEMBERS of the
                        query: the label's own family test (synvoy_specificity.pairwise_evalue,
                        blastp -comp_based_stats 2) at E <= 1e-3. Each entry: gene_id,
                        signed rank (-1 = nearest on the left), distance_bp, evalue_seg_no,
                        evalue_seg_yes. The list follows the SEG setting the label ships
                        with (FAMILY_TEST_SEG). null without a home annotation, and null
                        without blastp: the test's fallback (plain Karlin-Altschul) is the
                        statistic that calls low-complexity look-alikes a family.
  relative_adjacent     a relative at rank +-1 (MCScanX: tandem)
  relative_proximal     relatives present, none adjacent (MCScanX: proximal)
  scaffold_len, dist_to_scaffold_end, home_genome_bp
  telomeric_repeat_at_end  the scaffold end nearest to the gene carries a run of a common
                        telomeric repeat (TTAGGG, TTAGG, TTTAGGG; either strand)

Reads only the scaffold of the gene of interest from the genome and the annotation, and
streams the sequence lengths. Not a label input: it decides which preset module a run
uses (proposal of 2026-10-15, A1 local_family and C2 home_context).

    home_context.py --query q.faa --genome home.fna --gff home.gff --goi_bed goi.bed \\
                    --n_flanking_genes 10 --output home_context.json
"""
import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

try:
    import synvoy_specificity
    from sequence_utils import (cds_model_protein, cds_models_by_gene, load_genome_subset,
                                parse_bed, stream_fasta_lengths)
    from query_features import first_sequence
except ImportError:  # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import synvoy_specificity
    from sequence_utils import (cds_model_protein, cds_models_by_gene, load_genome_subset,
                                parse_bed, stream_fasta_lengths)
    from query_features import first_sequence

TELOMERE_MOTIFS = ("TTAGGG", "TTAGG", "TTTAGGG")
TELOMERE_END_BP = 5000         # how much of the scaffold end is looked at
TELOMERE_MIN_COPIES = 10       # copies of the motif (either strand) that make a telomeric run


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def goi_span(bed_path: str) -> Optional[Tuple[str, int, int]]:
    """(chrom, start, end) of the gene of interest (1-based, inclusive): the span of the
    BED rows on the chromosome that holds most of them."""
    rows = parse_bed(bed_path) if bed_path and os.path.exists(bed_path) else []
    if not rows:
        return None
    by_chrom: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        by_chrom.setdefault(r["chrom"], []).append(r)
    chrom = sorted(by_chrom.items(), key=lambda kv: (-len(kv[1]), kv[0]))[0][0]
    return (chrom, min(int(r["start"]) for r in by_chrom[chrom]) + 1,
            max(int(r["end"]) for r in by_chrom[chrom]))


def telomere_motif(seq_end: str) -> Optional[str]:
    """The telomeric repeat found at least TELOMERE_MIN_COPIES times in ``seq_end`` (either
    strand), or None."""
    up = seq_end.upper()
    for motif in TELOMERE_MOTIFS:
        if up.count(motif) >= TELOMERE_MIN_COPIES or up.count(_revcomp(motif)) >= TELOMERE_MIN_COPIES:
            return motif
    return None


def family_evalues(query: str, protein: str) -> Tuple[float, float]:
    """(E with SEG off, E with SEG on) of the label's family test. synvoy_specificity reads
    FAMILY_TEST_SEG when it runs; the shipped value is restored after each call."""
    shipped = synvoy_specificity.FAMILY_TEST_SEG
    out = []
    try:
        for seg in ("no", "yes"):
            synvoy_specificity.FAMILY_TEST_SEG = seg
            out.append(synvoy_specificity.pairwise_evalue(query, protein))
    finally:
        synvoy_specificity.FAMILY_TEST_SEG = shipped
    return out[0], out[1]


def neighbours(models: Dict[str, Dict[str, Any]], chrom: str, start: int, end: int,
               n: int) -> Tuple[Optional[str], List[Tuple[str, int, int, Dict[str, Any]]]]:
    """(gene of interest, [(gene_id, signed rank, distance_bp, model)]) of the coding genes
    on ``chrom`` within +-n ranks of [start, end]. The gene of interest is the gene whose
    coding span overlaps the span most (None when none does); it is not a neighbour."""
    genes = []
    for gid, info in models.items():
        if info["chrom"] != chrom:
            continue
        g0 = min(e[0] for e in info["exons"])
        g1 = max(e[1] for e in info["exons"])
        genes.append((g0, g1, gid, info))
    genes.sort(key=lambda t: (t[0], t[1], t[2]))
    best, best_ov = None, 0
    for g0, g1, gid, _info in genes:
        ov = min(g1, end) - max(g0, start) + 1
        if ov > best_ov:
            best, best_ov = gid, ov
    left = [g for g in genes if g[2] != best and g[1] < start]
    right = [g for g in genes if g[2] != best and g[0] > end]
    out = []
    for i, (g0, g1, gid, info) in enumerate(reversed(left[-n:] if n else []), 1):
        out.append((gid, -i, start - g1, info))
    for i, (g0, g1, gid, info) in enumerate(right[:n], 1):
        out.append((gid, i, g0 - end, info))
    return best, out


def home_context(query: str, genome: str, gff: Optional[str], goi_bed: str,
                 n_flanking: int) -> Dict[str, Any]:
    lengths = stream_fasta_lengths(genome)
    ctx: Dict[str, Any] = {
        "goi_gene": None, "relatives_nearby": None, "n_neighbours_tested": 0,
        "relative_adjacent": None, "relative_proximal": None,
        "family_test_seg": synvoy_specificity.FAMILY_TEST_SEG,
        "statistic": "blastp" if synvoy_specificity._blastp() else "fallback",
        "scaffold": None, "scaffold_len": None, "dist_to_scaffold_end": None,
        "home_genome_bp": sum(lengths.values()), "telomeric_repeat_at_end": None,
    }
    span = goi_span(goi_bed)
    if span is None:
        return ctx
    chrom, start, end = span
    ctx["scaffold"] = chrom
    if chrom not in lengths:
        return ctx
    scaffold_len = lengths[chrom]
    ctx["scaffold_len"] = scaffold_len
    ctx["dist_to_scaffold_end"] = max(0, min(start - 1, scaffold_len - end))
    seq = load_genome_subset(genome, [chrom]).get(chrom, "")
    near_left = (start - 1) <= (scaffold_len - end)
    piece = seq[:TELOMERE_END_BP] if near_left else seq[-TELOMERE_END_BP:]
    ctx["telomeric_repeat_at_end"] = telomere_motif(piece) is not None
    if not gff or os.path.basename(gff) == "NO_GFF" or not os.path.exists(gff):
        return ctx
    models = cds_models_by_gene(gff, chroms={chrom})
    goi_gene, near = neighbours(models, chrom, start, end, n_flanking)
    ctx["goi_gene"] = goi_gene
    ctx["n_neighbours_tested"] = len(near)
    if ctx["statistic"] != "blastp":
        return ctx
    shipped_seg = synvoy_specificity.FAMILY_TEST_SEG
    relatives = []
    for gid, rank, dist, info in near:
        protein = cds_model_protein(info, seq)
        if not protein:
            continue
        e_no, e_yes = family_evalues(query, protein)
        e_used = e_yes if shipped_seg == "yes" else e_no
        if e_used <= synvoy_specificity.FAMILY_MEMBER_EVALUE:
            relatives.append({"gene_id": gid, "rank": rank, "distance_bp": max(0, dist),
                              "evalue_seg_no": float(f"{e_no:.3g}"),
                              "evalue_seg_yes": float(f"{e_yes:.3g}")})
    relatives.sort(key=lambda r: (abs(r["rank"]), r["rank"]))
    ctx["relatives_nearby"] = relatives
    ctx["relative_adjacent"] = any(abs(r["rank"]) == 1 for r in relatives)
    ctx["relative_proximal"] = bool(relatives) and not ctx["relative_adjacent"]
    return ctx


def main() -> None:
    ap = argparse.ArgumentParser(description="Home neighbourhood and scaffold of the gene of interest")
    ap.add_argument("--query", required=True, help="Query protein FASTA (first record)")
    ap.add_argument("--genome", required=True, help="Home genome FASTA")
    ap.add_argument("--gff", default="NO_GFF", help="Home annotation (GFF3), or NO_GFF")
    ap.add_argument("--goi_bed", required=True, help="goi_annotation.bed of ANNOTATE_GOI")
    ap.add_argument("--n_flanking_genes", type=int, default=10,
                    help="Neighbours tested on each side of the gene of interest (default 10)")
    ap.add_argument("--output", required=True, help="home_context.json")
    args = ap.parse_args()
    ctx = home_context(first_sequence(args.query), args.genome, args.gff, args.goi_bed,
                       args.n_flanking_genes)
    with open(args.output, "w") as fh:
        json.dump(ctx, fh, indent=2, sort_keys=True)
        fh.write("\n")
    rel = ctx["relatives_nearby"]
    print(f"home_context: gene {ctx['goi_gene']} on {ctx['scaffold']} "
          f"({ctx['dist_to_scaffold_end']} bp from its end); relatives within "
          f"+-{args.n_flanking_genes} genes: "
          + ("not tested" if rel is None else ", ".join(f"{r['gene_id']} ({r['rank']:+d})" for r in rel) or "none"))


if __name__ == "__main__":
    main()
