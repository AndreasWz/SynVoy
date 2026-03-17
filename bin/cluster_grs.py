#!/usr/bin/env python3

import argparse
import sys
import csv
import os
import re
import random
from collections import defaultdict

# No BioPython needed - we parse FASTA manually for genome length

def parse_args():
    parser = argparse.ArgumentParser(description="Cluster and Score Synteny Regions (Prioritizer Mode)")
    parser.add_argument("--hits", required=True, help="Input MMseqs hits (m8)")
    parser.add_argument("--synteny_bed", required=True, help="BED file defining the expected synteny block (genes in order)")
    parser.add_argument("--genome", required=False, default=None,
                        help="Target Genome FASTA (optional; used for approximate p-value context)")
    parser.add_argument("--output", required=True, help="Output Region BED")
    parser.add_argument(
        "--scores_output",
        required=False,
        default=None,
        help="Optional TSV with structured scores for emitted regions",
    )
    parser.add_argument(
        "--target_gff",
        required=False,
        default=None,
        help="Optional target GFF from iterative search (used to prioritize GOI-overlapping regions)",
    )
    parser.add_argument(
        "--goi_padding",
        type=int,
        default=20000,
        help="Padding (bp) around GOI intervals when injecting fallback candidate regions",
    )
    parser.add_argument("--flanking_count", type=int, default=10, help="Expected number of flanking genes (fallback)")
    parser.add_argument("--cluster_distance", type=int, default=50000, help="Max distance to cluster hits (bp)")
    parser.add_argument("--min_score", type=float, default=0.5, help="Score threshold for High Confidence")
    parser.add_argument("--weight_base", type=float, default=0.4, help="Base weight for coverage")
    parser.add_argument("--weight_consistency", type=float, default=0.3, help="Weight for order consistency")
    parser.add_argument("--weight_strand", type=float, default=0.3, help="Weight for strand consistency")
    parser.add_argument(
        "--goi_overlap_bonus",
        type=float,
        default=0.4,
        help="Additive score bonus for clusters overlapping GOI intervals",
    )
    parser.add_argument(
        "--max_regions",
        type=int,
        default=0,
        help="Max regions to output (0 = adaptive: emit all above relative threshold, capped at 6)",
    )
    return parser.parse_args()

def load_synteny_map(bed_file):
    """
    Load the expected order of genes from the Home/Query BED file.
    Returns: dict { 'gene_name': rank_index }
    """
    gene_map = {}
    if not os.path.exists(bed_file):
        return gene_map
        
    try:
        rank = 0
        with open(bed_file) as f:
            for line in f:
                if not line.strip(): continue
                parts = line.strip().split('\t')
                # BED: chrom start end name ...
                if len(parts) >= 4:
                    gene_name = parts[3]
                    # Handle pipe format if present in inputs (clean ID)
                    clean_name = gene_name.split('|')[0]
                    gene_map[clean_name] = rank
                    # Also map full name if distinct
                    if gene_name != clean_name:
                        gene_map[gene_name] = rank
                    rank += 1
    except Exception as e:
        print(f"Error loading synteny map: {e}", file=sys.stderr)
    return gene_map


def parse_gff_attrs(attr_field):
    attrs = {}
    for kv in (attr_field or "").split(";"):
        if "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def merge_intervals(intervals):
    if not intervals:
        return []

    by_chrom = defaultdict(list)
    for iv in intervals:
        by_chrom[iv["chrom"]].append((iv["start"], iv["end"]))

    merged = []
    for chrom, spans in by_chrom.items():
        spans.sort()
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                merged.append({"chrom": chrom, "start": cur_s, "end": cur_e})
                cur_s, cur_e = s, e
        merged.append({"chrom": chrom, "start": cur_s, "end": cur_e})
    return merged


def load_goi_intervals_from_gff(gff_file, padding_bp=20000):
    """
    Parse GOI intervals from SynTerra iterative target GFF.
    GOI records are detected via Name/ID/SynTerra_Parent containing GOI_.
    """
    if not gff_file or not os.path.exists(gff_file):
        return []

    goi_intervals = []
    accepted_types = {"mRNA", "gene", "tandem_copy", "transcript", "mrna"}

    try:
        with open(gff_file) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] not in accepted_types:
                    continue

                attrs = parse_gff_attrs(parts[8])
                tokens = [
                    attrs.get("Name", ""),
                    attrs.get("ID", ""),
                    attrs.get("SynTerra_Parent", ""),
                    attrs.get("Parent", ""),
                ]
                marker = "|".join(tokens)
                if "GOI_" not in marker:
                    continue

                try:
                    s = int(parts[3])
                    e = int(parts[4])
                except ValueError:
                    continue
                if e < s:
                    s, e = e, s

                s = max(0, s - max(0, int(padding_bp)))
                e = e + max(0, int(padding_bp))
                goi_intervals.append({"chrom": parts[0], "start": s, "end": e})
    except Exception as exc:
        print(f"WARNING: Could not parse GOI intervals from target GFF: {exc}", file=sys.stderr)
        return []

    return merge_intervals(goi_intervals)

def _is_goi_query_name(query_name):
    """
    Identify GOI-derived query names in m8 hits.
    """
    if not query_name:
        return False
    q = str(query_name).strip()
    if q.startswith("GOI_") or q.startswith("GOI|"):
        return True
    base = q.split("|")[0]
    if base.startswith("GOI_"):
        return True
    # Backward compatibility for legacy bare exon IDs from GOI expansion.
    if re.fullmatch(r"exon_\d+", base):
        return True
    return False

def load_goi_intervals_from_hits(hits, padding_bp=20000):
    """
    Build GOI intervals from m8 hits when target GFF GOI models are missing.
    """
    if not hits:
        return []

    intervals = []
    pad = max(0, int(padding_bp))
    for h in hits:
        if not _is_goi_query_name(h.get("query", "")):
            continue
        chrom = h.get("chrom")
        if not chrom:
            continue
        start = int(h.get("start", 0))
        end = int(h.get("end", 0))
        if end < start:
            start, end = end, start
        intervals.append(
            {
                "chrom": chrom,
                "start": max(0, start - pad),
                "end": end + pad,
            }
        )

    return merge_intervals(intervals)


def overlaps_interval(chrom, start, end, interval):
    if chrom != interval["chrom"]:
        return False
    return not (end < interval["start"] or start > interval["end"])


def cluster_overlaps_goi(cluster, goi_intervals):
    if not goi_intervals:
        return False
    for iv in goi_intervals:
        if overlaps_interval(cluster["chrom"], cluster["start"], cluster["end"], iv):
            return True
    return False


def build_goi_anchor_clusters(goi_intervals, existing_clusters):
    """
    Inject GOI anchor regions when score-ranked clusters miss GOI loci entirely.
    """
    anchors = []
    for iv in goi_intervals:
        covered = False
        for cl in existing_clusters:
            if overlaps_interval(cl["chrom"], cl["start"], cl["end"], iv):
                covered = True
                break
        if covered:
            continue
        anchors.append(
            {
                "cluster": [],
                "unique": 1,
                "consistency": 1.0,
                "strand_cons": 1.0,
                "score": 1.0,
                "p_value": 0.0,
                "start": iv["start"],
                "end": iv["end"],
                "chrom": iv["chrom"],
                "is_goi_anchor": True,
                "goi_overlap": True,
            }
        )
    return anchors

def get_genome_length(genome_file):
    """Get total length of genome (sum of sequence lengths)."""
    total_len = 0
    try:
        fai = genome_file + ".fai"
        if os.path.exists(fai):
            with open(fai) as f:
                for line in f:
                    parts = line.split('\t')
                    total_len += int(parts[1])
        else:
            with open(genome_file) as f:
                for line in f:
                    if not line.startswith('>'):
                        total_len += len(line.strip())
    except Exception as e:
        print(f"WARNING: Could not determine genome length from {genome_file}: {e}", file=sys.stderr)
    return max(1, total_len)

def cluster_hits_proximity(hits, gene_map, max_dist):
    """
    Cluster hits based on genomic proximity in Target.
    hits: list of dicts {'query', 'chrom', 'start', 'end', ...}
    """
    if not hits: return []
    
    # Sort by Chrom, Start
    hits.sort(key=lambda x: (x['chrom'], x['start']))
    
    clusters = []
    current_cluster = [hits[0]]
    
    for i in range(1, len(hits)):
        h = hits[i]
        prev = current_cluster[-1]
        
        # Check if same chrom and close enough
        dist = h['start'] - prev['end']
        
        if h['chrom'] == prev['chrom'] and dist < max_dist:
            current_cluster.append(h)
        else:
            clusters.append(current_cluster)
            current_cluster = [h]
            
    clusters.append(current_cluster)
    return clusters

def score_flexible_synteny(cluster, gene_map):
    """
    Score a cluster based on synteny preservation.
    Returns: (unique_genes_count, consistency_score, strand_score)
    """
    # 1. Unique Genes Coverage
    # Map hits to ranks
    hit_ranks = []
    seen_genes = set()
    
    strand_matches = 0
    total_strand_ops = 0
    
    for h in cluster:
        q = h['query']
        # Try to find rank
        rank = -1
        # Try full match or split
        if q in gene_map:
            rank = gene_map[q]
        else:
            q_clean = q.split('|')[0]
            if q_clean in gene_map:
                 rank = gene_map[q_clean]
        
        if rank != -1:
            hit_ranks.append({'rank': rank, 't_pos': h['start'], 'h_strand': h['strand']})
            seen_genes.add(rank)
            
    unique_genes = len(seen_genes)
    if not hit_ranks:
        return unique_genes, 0.0, 0.0
        
    # 2. Consistency (Order preservation)
    # Check if ranks are increasing (or decreasing if inverted)
    # LIS (Longest Increasing Subsequence) could be used
    # But simple pair consistency is faster/robust
    
    # Sort by target position (already filtered/clustered this way)
    # hit_ranks is naturally sorted by t_pos
    
    if len(hit_ranks) < 2:
        return unique_genes, 1.0, 1.0 # Single gene is consistent
        
    # Determine dominant direction
    increasing = 0
    decreasing = 0
    for i in range(len(hit_ranks)-1):
        r1 = hit_ranks[i]['rank']
        r2 = hit_ranks[i+1]['rank']
        if r2 > r1: increasing += 1
        elif r2 < r1: decreasing += 1
        
    total_pairs = len(hit_ranks) - 1
    if total_pairs > 0:
        consistency = max(increasing, decreasing) / total_pairs
    else:
        consistency = 1.0
        
    # 3. Strand Consistency
    # If Increasing (Home + -> Target +), Strands should match?
    # Not necessarily depending on annotation.
    # Usually: if Inversion (ranks decreasing), Strands should be opposite of Home?
    # Simple check: Are all hits on same strand? Or consistent with direction?
    # Let's use simple "majority strand" logic for the block relative to query.
    # Assuming Query genes are all + in the BED (usually).
    # If Target Block is +, then query strands should be Match.
    # If Target Block is -, query strands should be Inverted?
    # Simplify: Fraction of hits on the Majority Strand of the cluster.
    
    plus_cnt = sum(1 for h in cluster if h['strand'] == '+')
    minus_cnt = len(cluster) - plus_cnt
    majority = max(plus_cnt, minus_cnt)
    strand_cons = majority / len(cluster)
    
    return unique_genes, consistency, strand_cons

def estimate_pvalue(observed_score, cluster_hits, all_hits, genome_len, cluster_distance, score_func, gene_map, total_genes_expected, n=200):
    """
    Estimate P-value via label-shuffling permutation test.

    Shuffles gene-label assignments among the cluster's hits and re-scores,
    counting how often a random labelling achieves the observed score or better.
    This tests whether the *identity* of genes in the cluster is non-random
    (i.e., they match the expected synteny map better than chance).
    """
    if not cluster_hits or not gene_map or total_genes_expected <= 0:
        return 1.0

    # Pool of all query labels seen across the genome
    all_queries = [h['query'].split('|')[0] for h in all_hits]
    if not all_queries:
        return 1.0

    cluster_size = len(cluster_hits)
    better_or_equal = 0

    for _ in range(n):
        # Create a synthetic cluster by shuffling query labels
        shuffled = list(cluster_hits)  # shallow copy
        sampled_labels = random.choices(all_queries, k=cluster_size)
        fake_cluster = []
        for hit, label in zip(shuffled, sampled_labels):
            fake_hit = dict(hit)
            fake_hit['query'] = label
            fake_cluster.append(fake_hit)

        rand_unique, rand_consistency, rand_strand = score_func(fake_cluster, gene_map)
        rand_coverage = rand_unique / total_genes_expected if total_genes_expected > 0 else 0
        # Use same weighted formula as the real scoring
        rand_score = rand_coverage * (rand_consistency * 0.5 + rand_strand * 0.5)
        if rand_score >= observed_score:
            better_or_equal += 1

    return (better_or_equal + 1) / (n + 1)  # Laplace-corrected

def main():
    args = parse_args()
    
    gene_map = load_synteny_map(args.synteny_bed)
    # gene_map contains alias keys (raw + cleaned IDs) that may point to the
    # same rank; coverage must use unique ranks, not raw key count.
    # Exclude GOI-derived entries from the denominator — those are the targets
    # of the search, not flanking synteny anchors. Without this, max coverage
    # can never reach 1.0.
    flanking_ranks = set()
    for key, rank in gene_map.items():
        if not key.startswith('GOI_'):
            flanking_ranks.add(rank)
    total_genes_expected = len(flanking_ranks) if flanking_ranks else len(set(gene_map.values()))
    if total_genes_expected == 0: total_genes_expected = args.flanking_count
    
    hits = []
    if os.path.exists(args.hits):
        with open(args.hits) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('query'):
                    continue
                try:
                    # MMseqs FMT: query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
                    #              0     1      2      3       4        5       6      7      8      9     10     11
                    t_start = int(row[8])
                    t_end = int(row[9])

                    # Detect strand and normalize coordinates.
                    strand = "+"
                    if t_start > t_end:
                        strand = "-"
                        t_start, t_end = t_end, t_start

                    # Convert from 1-based (m8 format) to 0-based half-open (BED/Python)
                    t_start = t_start - 1

                    h = {
                        'query': row[0],
                        'chrom': row[1],
                        'start': t_start,
                        'end': t_end,
                        'strand': strand,
                        'evalue': float(row[10])
                    }
                    hits.append(h)
                except ValueError:
                    continue
    else:
        print(f"INFO: Hits file {args.hits} not found. Continuing with GOI-anchor fallback if possible.", file=sys.stderr)

    goi_intervals_gff = load_goi_intervals_from_gff(args.target_gff, padding_bp=args.goi_padding)
    goi_intervals_hits = load_goi_intervals_from_hits(hits, padding_bp=args.goi_padding)
    goi_intervals = merge_intervals(goi_intervals_gff + goi_intervals_hits)
    if goi_intervals_hits and not goi_intervals_gff:
        print(
            f"INFO: Derived {len(goi_intervals_hits)} GOI interval(s) from hit file (no GOI GFF intervals).",
            file=sys.stderr,
        )
    elif goi_intervals_hits:
        print(
            f"INFO: Added GOI support from hits ({len(goi_intervals_hits)} interval(s)) on top of GFF GOI intervals.",
            file=sys.stderr,
        )

    # Cluster
    clusters = cluster_hits_proximity(hits, gene_map, args.cluster_distance)
    
    genome_len = get_genome_length(args.genome) if args.genome else 1
    
    scored_clusters = []
    for cl in clusters:
        unique_genes, consistency, strand_cons = score_flexible_synteny(cl, gene_map)
        
        # Composite Score
        if total_genes_expected > 0:
            coverage_score = unique_genes / total_genes_expected
        else:
            coverage_score = 0
        
        # Weighted quality: each weight is a true proportion of the quality score.
        # Uses coverage as base metric, modulated by order consistency and strand coherence.
        quality_score = (
            args.weight_base * coverage_score +
            args.weight_consistency * consistency +
            args.weight_strand * strand_cons
        )
        
        final_score = quality_score * coverage_score
        cluster_chrom = cl[0]['chrom']
        cluster_start = min(h['start'] for h in cl)
        cluster_end = max(h['end'] for h in cl)
        goi_overlap = cluster_overlaps_goi(
            {
                "chrom": cluster_chrom,
                "start": cluster_start,
                "end": cluster_end,
            },
            goi_intervals,
        )
        if goi_overlap:
            # Additive bonus capped at 0.15 to avoid GOI signal dominating synteny evidence
            final_score += min(0.15, max(0.0, float(args.goi_overlap_bonus)))
        
        p_val = estimate_pvalue(final_score, cl, hits, genome_len, args.cluster_distance, score_flexible_synteny, gene_map, total_genes_expected, n=200)
        
        scored_clusters.append({
            'cluster': cl,
            'unique': unique_genes,
            'consistency': consistency,
            'strand_cons': strand_cons,
            'coverage_score': coverage_score,
            'quality_score': quality_score,
            'score': final_score,
            'p_value': p_val,
            'start': cluster_start,
            'end': cluster_end,
            'chrom': cluster_chrom,
            'goi_overlap': goi_overlap,
        })

    # Sort: Score desc, P-value asc
    scored_clusters.sort(key=lambda x: (-x['score'], x['p_value']))

    # GOI-aware prioritization:
    # 1) Prefer clusters overlapping GOI intervals from iterative target GFF.
    # 2) If none overlap, inject GOI-anchor region(s) so the true locus is never dropped.
    ordered_clusters = list(scored_clusters)
    if goi_intervals:
        goi_clusters = [cl for cl in scored_clusters if cluster_overlaps_goi(cl, goi_intervals)]
        non_goi_clusters = [cl for cl in scored_clusters if not cluster_overlaps_goi(cl, goi_intervals)]
        if goi_clusters:
            for cl in goi_clusters:
                cl["goi_overlap"] = True
            for cl in non_goi_clusters:
                cl["goi_overlap"] = False
            ordered_clusters = goi_clusters + non_goi_clusters
        else:
            anchors = build_goi_anchor_clusters(goi_intervals, scored_clusters)
            if anchors:
                print(
                    f"INFO: No score-ranked cluster overlapped GOI; injected {len(anchors)} GOI-anchor region(s).",
                    file=sys.stderr,
                )
                ordered_clusters = anchors + scored_clusters

    selected_rows = []

    with open(args.output, 'w') as f_out:
        if not ordered_clusters:
            # Case 1: No clusters formed and no GOI fallback available.
            print(f"INFO: No synteny clusters could be formed for this genome.", file=sys.stderr)
        else:
            # Case 2: Clusters found — adaptive region selection
            if args.max_regions > 0:
                # User-specified hard cap
                num_to_output = min(args.max_regions, len(ordered_clusters))
            else:
                # Adaptive: emit all regions whose score is at least 30% of
                # the best region's score OR that have >=3 unique genes,
                # guaranteeing at least 1 and capping at 6.
                best_score = ordered_clusters[0]['score'] if ordered_clusters else 0
                score_floor = max(0.03, best_score * 0.3)
                num_to_output = 0
                for cl in ordered_clusters:
                    if cl['score'] >= score_floor or cl.get('unique', 0) >= 3:
                        num_to_output += 1
                    else:
                        break  # clusters are sorted by score desc
                num_to_output = max(1, min(num_to_output, 6))

            print(
                f"INFO: Emitting {num_to_output}/{len(ordered_clusters)} regions "
                f"(adaptive selection, best_score={ordered_clusters[0]['score']:.2f}).",
                file=sys.stderr,
            )
            
            for i in range(num_to_output):
                best = ordered_clusters[i]
                
                # Determine Confidence
                if best['score'] >= args.min_score:
                    confidence = "HIGH"
                elif best['score'] >= (args.min_score * 0.5):
                    confidence = "MEDIUM"
                else:
                    confidence = "LOW"
                
                # Determine Region Strand
                plus_cnt = sum(1 for h in best['cluster'] if h['strand'] == '+')
                minus_cnt = len(best['cluster']) - plus_cnt
                region_strand = "-" if minus_cnt > plus_cnt else "+"
                
                if best.get("is_goi_anchor"):
                    name = f"Reg{i+1}_GOI_anchor_C{confidence}_S{best['score']:.2f}"
                else:
                    name = f"Reg{i+1}_G{best['unique']}_C{confidence}_S{best['score']:.2f}"

                if args.max_regions > 0:
                    selected_reason = "user_cap"
                elif best['score'] >= score_floor:
                    selected_reason = "score_floor"
                elif best.get('unique', 0) >= 3:
                    selected_reason = "unique_gene_floor"
                else:
                    selected_reason = "adaptive_backstop"
                
                f_out.write(f"{best['chrom']}\t{best['start']}\t{best['end']}\t{name}\t{best['score']:.2f}\t{region_strand}\n")

                selected_rows.append({
                    "region_rank": i + 1,
                    "region_name": name,
                    "chrom": best["chrom"],
                    "start": best["start"],
                    "end": best["end"],
                    "strand": region_strand,
                    "score": best["score"],
                    "quality_score": best["quality_score"],
                    "coverage_score": best["coverage_score"],
                    "unique_genes": best["unique"],
                    "total_genes_expected": total_genes_expected,
                    "consistency": best["consistency"],
                    "strand_consistency": best["strand_cons"],
                    "p_value": best["p_value"],
                    "goi_overlap": bool(best.get("goi_overlap")),
                    "is_goi_anchor": bool(best.get("is_goi_anchor")),
                    "confidence": confidence,
                    "selection_reason": selected_reason,
                })

                # Emit compact summary to stdout for testability and quick diagnostics.
                print(
                    f"Score: {best['score']:.2f} | "
                    f"Consistency: {best['consistency']:.2f} | "
                    f"Strand: {best['strand_cons']:.2f}"
                )
                
                # Log Low Confidence
                goi_tag = " [GOI]" if best.get("goi_overlap") or best.get("is_goi_anchor") else ""
                if confidence == "LOW":
                    print(f"WARNING: Region {i+1} has LOW confidence (score={best['score']:.2f}), "
                          f"Genes: {best['unique']}/{total_genes_expected}, "
                          f"Consistency: {best['consistency']:.2f}{goi_tag}", file=sys.stderr)
                else:
                    print(f"Region {i+1}: {best['chrom']}:{best['start']}-{best['end']} "
                          f"(Score: {best['score']:.2f}, Conf: {confidence}){goi_tag}", file=sys.stderr)

    if args.scores_output:
        fieldnames = [
            "region_rank",
            "region_name",
            "chrom",
            "start",
            "end",
            "strand",
            "score",
            "quality_score",
            "coverage_score",
            "unique_genes",
            "total_genes_expected",
            "consistency",
            "strand_consistency",
            "p_value",
            "goi_overlap",
            "is_goi_anchor",
            "confidence",
            "selection_reason",
        ]
        with open(args.scores_output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in selected_rows:
                writer.writerow(row)

if __name__ == "__main__":
    main()
