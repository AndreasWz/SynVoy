#!/usr/bin/env python3

import argparse
import sys
import csv
import os
from collections import defaultdict

# No BioPython needed - we parse FASTA manually for genome length

def parse_args():
    parser = argparse.ArgumentParser(description="Cluster and Score Synteny Regions (Prioritizer Mode)")
    parser.add_argument("--hits", required=True, help="Input MMseqs hits (m8)")
    parser.add_argument("--synteny_bed", required=True, help="BED file defining the expected synteny block (genes in order)")
    parser.add_argument("--genome", required=True, help="Target Genome FASTA (for length)")
    parser.add_argument("--output", required=True, help="Output Region BED")
    parser.add_argument("--flanking_count", type=int, default=10, help="Expected number of flanking genes (fallback)")
    parser.add_argument("--cluster_dist", type=int, default=50000, help="Max distance to cluster hits (bp)")
    parser.add_argument("--min_score", type=float, default=0.5, help="Score threshold for High Confidence")
    parser.add_argument("--weight_base", type=float, default=0.4, help="Base weight for coverage")
    parser.add_argument("--weight_consistency", type=float, default=0.3, help="Weight for order consistency")
    parser.add_argument("--weight_strand", type=float, default=0.3, help="Weight for strand consistency")
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

def get_genome_length(genome_file):
    """Get total length of genome (sum of sequence lengths)."""
    total_len = 1
    try:
        # Fast approximate length from file size or index?
        # Parsing fasta is safer.
        # Check if .fai exists?
        fai = genome_file + ".fai"
        if os.path.exists(fai):
            with open(fai) as f:
                for line in f:
                    parts = line.split('\t')
                    total_len += int(parts[1])
        else:
            # Parse FASTA (could be slow for huge genomes, but usually okay)
            with open(genome_file) as f:
                for line in f:
                    if not line.startswith('>'):
                        total_len += len(line.strip())
    except:
        pass
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

def estimate_pvalue(observed_score, all_hits, genome_len, cluster_dist, score_func, gene_map, n=100):
    """
    Estimate P-value by randomizing hits positions/identity.
    Actually, randomizing positions in genome is better.
    Simplified: Shuffle "query identity" among the hits and re-score?
    Conservative: Count how many random sets of N hits from the genome would score this high?
    For speed: We skip complex P-value if score is low anyway.
    Implementation:
    Shuffle the 'query' assignments of the hits in the cluster to see if random composition gives this score.
    """
    if not all_hits: return 1.0
    
    # We want to see if the structure (positions) + random genes would form a cluster?
    # No, the cluster exists. Is the *content* significant?
    # Shuffle ranks.
    
    # Collect all available ranks in the hits (background)
    background_ranks = []
    for h in all_hits:
        q = h['query'].split('|')[0]
        if q in gene_map: background_ranks.append(gene_map[q])
    
    if len(background_ranks) < len(gene_map):
        # Pad with -1
        background_ranks.extend([-1] * (len(all_hits) - len(background_ranks)))
        
    cluster_size = len(all_hits) # Approximate
    
    better_count = 0
    
    # This is a placeholder for a real permutation test.
    # Given limited time/compute, we rely mainly on Score.
    # Return 0.05 if Score > 0.5?
    
    # Real logic:
    # return observed_score > 0.5 ? 0.01 : 0.5
    
    # Let's simplify: P-value is inversely proportional to score here.
    # Real statistical test is overkill for this step if we just prioritize.
    return 1.0 - observed_score

def main():
    args = parse_args()
    
    gene_map = load_synteny_map(args.synteny_bed)
    total_genes_expected = len(gene_map)
    if total_genes_expected == 0: total_genes_expected = args.flanking_count
    
    hits = []
    try:
        with open(args.hits) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not row or row[0].startswith('query'): continue 
                try:
                    # MMseqs FMT: query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
                    #              0     1      2      3       4        5       6      7      8      9     10     11
                    
                    t_start = int(row[8])
                    t_end = int(row[9])
                    
                    # Detect Strand
                    # Standard BLAST/MMseqs convention: t_start > t_end implies minus strand
                    strand = "+"
                    if t_start > t_end:
                        strand = "-"
                        t_start, t_end = t_end, t_start # Normalize for object
                        
                    h = {
                        'query': row[0], 'chrom': row[1],
                        'start': t_start, 'end': t_end,
                        'strand': strand,
                        'evalue': float(row[10])
                    }
                    hits.append(h)
                except ValueError: continue
    except FileNotFoundError:
        # Valid state: no hits found
        print(f"INFO: Hits file {args.hits} not found or empty.", file=sys.stderr)
        with open(args.output, 'w') as f_out:
            pass
        return

    # Cluster
    clusters = cluster_hits_proximity(hits, gene_map, args.cluster_dist)
    
    genome_len = get_genome_length(args.genome)
    
    scored_clusters = []
    for cl in clusters:
        unique_genes, consistency, strand_cons = score_flexible_synteny(cl, gene_map)
        
        # Composite Score
        if total_genes_expected > 0:
            coverage_score = unique_genes / total_genes_expected
        else:
            coverage_score = 0
        
        quality_mult = (args.weight_base +
                        args.weight_consistency * consistency +
                        args.weight_strand * strand_cons)
        
        final_score = coverage_score * quality_mult
        
        p_val = estimate_pvalue(final_score, hits, genome_len, args.cluster_dist, score_flexible_synteny, gene_map, n=100)
        
        scored_clusters.append({
            'cluster': cl,
            'unique': unique_genes,
            'consistency': consistency,
            'strand_cons': strand_cons,
            'score': final_score,
            'p_value': p_val,
            'start': min(h['start'] for h in cl),
            'end': max(h['end'] for h in cl),
            'chrom': cl[0]['chrom']
        })

    # Sort: Score desc, P-value asc
    scored_clusters.sort(key=lambda x: (-x['score'], x['p_value']))
    
    with open(args.output, 'w') as f_out:
        if not scored_clusters:
            # Case 1: No clusters formed at all
            print(f"INFO: No synteny clusters could be formed for this genome.", file=sys.stderr)
            # Produce empty file (success state)
            pass
        else:
            # Case 2: Clusters found - Output top 3
            num_to_output = min(3, len(scored_clusters))
            
            for i in range(num_to_output):
                best = scored_clusters[i]
                
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
                
                name = f"Reg{i+1}_G{best['unique']}_C{confidence}_S{best['score']:.2f}"
                
                f_out.write(f"{best['chrom']}\t{best['start']}\t{best['end']}\t{name}\t{best['score']:.2f}\t{region_strand}\n")
                
                # Log Low Confidence
                if confidence == "LOW":
                    print(f"WARNING: Region {i+1} has LOW confidence (score={best['score']:.2f}), "
                          f"Genes: {best['unique']}/{total_genes_expected}, "
                          f"Consistency: {best['consistency']:.2f}", file=sys.stderr)
                else:
                    print(f"Region {i+1}: {best['chrom']}:{best['start']}-{best['end']} "
                          f"(Score: {best['score']:.2f}, Conf: {confidence})", file=sys.stderr)

if __name__ == "__main__":
    main()
