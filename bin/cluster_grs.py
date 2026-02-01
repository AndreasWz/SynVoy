#!/usr/bin/env python3
import argparse
import csv
import random
import os
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

def parse_args():
    parser = argparse.ArgumentParser(description="Cluster hits into Genomic Regions (GRs) with collinearity checks")
    parser.add_argument("--hits", required=True, help="Input MMseqs2 hits (m8 format)")
    parser.add_argument("--synteny_bed", required=True, help="Original synteny block BED (provides gene order)")
    parser.add_argument("--flanking_count", type=int, required=True, help="Total number of flanking genes in query")
    parser.add_argument("--genome", help="Target genome FASTA (for length estimation/p-value)")
    parser.add_argument("--cluster_dist", type=int, default=50000, help="Initial clustering distance (bp)")
    parser.add_argument("--min_score", type=float, default=0.6, help="Minimum synteny score (0.0-1.0)")
    parser.add_argument("--output", required=True, help="Output BED file")
    return parser.parse_args()

def load_synteny_map(bed_file):
    """
    Parse synteny BED to map query gene IDs to their index (0 to N-1) AND strand.
    """
    gene_map = {}
    try:
        with open(bed_file) as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    strand = parts[5] if len(parts) > 5 else "+" 
                    gene_map[name] = {'index': i, 'strand': strand}
    except Exception as e:
        print(f"Error reading synteny BED: {e}")
    return gene_map

def get_genome_length(fasta_path):
    """Get total length of genome."""
    total_len = 0
    if not fasta_path or not os.path.exists(fasta_path):
        return 0
    try:
        if os.path.exists(fasta_path + ".fai"):
             with open(fasta_path + ".fai") as f:
                 for line in f:
                     total_len += int(line.split('\t')[1])
        elif SeqIO:
            for record in SeqIO.parse(fasta_path, "fasta"):
                total_len += len(record.seq)
        else:
            with open(fasta_path) as f:
                for line in f:
                    if not line.startswith('>'):
                        total_len += len(line.strip())
    except Exception as e:
        print(f"Error reading genome: {e}")
    return total_len

def cluster_hits_proximity(hits, gene_map, dist_threshold):
    """
    Cluster hits purely based on genomic proximity.
    """
    hits.sort(key=lambda x: (x['chrom'], x['start']))
    clusters = []
    if not hits:
        return clusters
        
    current_cluster = []
    
    for hit in hits:
        q_id = hit['query']
        if '|' in q_id: q_id = q_id.split('|')[0]
        
        if q_id not in gene_map: continue
        
        hit['order_index'] = gene_map[q_id]['index']
        # We don't store strand in 'hit' here yet? We should.
        # It's already in hit['strand'] from parsing.
        
        if not current_cluster:
            current_cluster.append(hit)
            continue
            
        last_hit = current_cluster[-1]
        
        # Chromosome check
        if hit['chrom'] != last_hit['chrom']:
            clusters.append(current_cluster)
            current_cluster = [hit]
            continue
            
        # Distance check
        gap = max(0, hit['start'] - last_hit['end'])
        if gap > dist_threshold:
            clusters.append(current_cluster)
            current_cluster = [hit]
            continue
            
        current_cluster.append(hit)
            
    if current_cluster:
        clusters.append(current_cluster)
        
    return clusters

def score_flexible_synteny(cluster, gene_map):
    """
    Score based on pairwise adjacency compatibility (allowing inversions/jumps).
    Also checks REF-STRAND consistency.
    """
    if len(cluster) < 2:
        return 1, 0.0, 1.0 # Min score
        
    # Sort by genomic position
    sorted_cluster = sorted(cluster, key=lambda x: x['start'])
    indices = [h['order_index'] for h in sorted_cluster]
    
    # 1. Adjacency Consistency
    pairwise_score = 0
    for i in range(len(indices) - 1):
        diff = abs(indices[i] - indices[i+1])
        if diff <= 3: # Allow up to 3 skipped genes or inversion
            pairwise_score += 1
            
    unique_genes = len(set(indices))
    consistency = pairwise_score / (len(indices) - 1) if len(indices) > 1 else 0
    
    # 2. Strand Consistency
    # We expect the relationship (QueryStrand vs HitStrand) to be constant.
    # Case A: Synteny (All Same): Q(+) -> H(+), Q(-) -> H(-) => Relation (+)
    # Case B: Inversion (All Inverted): Q(+) -> H(-), Q(-) -> H(+) => Relation (-)
    
    same_count = 0
    total_valid = 0
    
    for h in cluster:
        q_id = h['query'].split('|')[0]
        if q_id in gene_map:
            q_strand = gene_map[q_id]['strand'] # + or -
            h_strand = h['strand'] # + or -
            
            if q_strand == h_strand:
                same_count += 1
            total_valid += 1
            
    if total_valid > 0:
        # We take the MAX of (Same vs Diff) ratio.
        # If 90% are same -> 0.9. If 90% are diff -> 0.9 (Inversion).
        # If 50/50 -> 0.5 (Bad).
        diff_count = total_valid - same_count
        strand_consistency = max(same_count, diff_count) / total_valid
    else:
        strand_consistency = 1.0
    
    return unique_genes, consistency, strand_consistency

def estimate_pvalue(observed_score, all_hits, genome_len, window_size, score_func, gene_map, n=100):
    """
    Permutation test: Pick random windows in genome, check hits inside, score them.
    """
    if genome_len <= 0:
        max_coord = max(h['end'] for h in all_hits) if all_hits else 1000000
        # Ensure genome len is at least window size + buffer if we fall back
        genome_len = max(max_coord * 1.5, window_size * 2)
        
    random_scores = []
    
    for _ in range(n):
        search_space = int(genome_len - window_size)
        if search_space <= 0:
            rand_start = 0
        else:
            rand_start = random.randint(0, search_space)
        rand_end = rand_start + window_size
        
        hits_in_window = [h for h in all_hits if h['start'] >= rand_start and h['end'] <= rand_end]
        
        if not hits_in_window:
            random_scores.append(0)
            continue
            
        u, c, S = score_func(hits_in_window, gene_map)
        
        # Calculate random score using same formula
        # coverage = u / len(gene_map)
        # random_final = coverage * (0.5 + 0.3*c + 0.2*S)
        # But for p-value usually we just care about "Did we find as many unique genes?" (simplest proxy)
        # Or do we use the full score? 
        # Using full score is better.
        
        coverage = u / (len(gene_map) if len(gene_map) > 0 else 1)
        r_score = coverage # Simplify for p-val to just coverage? Or use the real metric
        
        random_scores.append(r_score) 
        
    count_better = sum(1 for s in random_scores if s >= observed_score)
    p_value = (count_better + 1) / (n + 1)
    
    return p_value

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
        open(args.output, 'w').close()
        return

    # Cluster
    clusters = cluster_hits_proximity(hits, gene_map, args.cluster_dist)
    
    genome_len = get_genome_length(args.genome)
    
    scored_clusters = []
    for cl in clusters:
        unique_genes, consistency, strand_cons = score_flexible_synteny(cl, gene_map)
        
        # Composite Score
        coverage_score = unique_genes / total_genes_expected
        
        # New Scoring Formula: Coverage * (WeightA + WeightB*Consistency + WeightC*Strand)
        # We punish bad strand/order but want Coverage to be king.
        # If strand is random (0.5), we penalize.
        
        quality_mult = (0.4 + 0.3 * consistency + 0.3 * strand_cons) 
        # range: 0.4 (worst) to 1.0 (best)
        
        final_score = coverage_score * quality_mult
        
        # Use coverage score for p-value to keep it simple?
        # Or use final score. Let's use coverage for p-value to represent "Cluster Density Significance"
        p_val = estimate_pvalue(coverage_score, hits, genome_len, args.cluster_dist, score_flexible_synteny, gene_map, n=100)
        
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
    
    if not scored_clusters:
        open(args.output, 'w').close()
        return

    best = scored_clusters[0]
    
    is_significant = best['p_value'] < 0.1 
    passes_score = best['score'] >= (args.min_score * 0.5) # Lower threshold since we multiply by quality
    
    print(f"Best Region: {best['chrom']}:{best['start']}-{best['end']}")
    print(f"  Score: {best['score']:.2f} (Unique: {best['unique']}/{total_genes_expected})")
    print(f"  Consistency: {best['consistency']:.2f}, Strand: {best['strand_cons']:.2f}")
    print(f"  P-value: {best['p_value']:.4f}")

    with open(args.output, 'w') as f_out:
        # We output if it's decent. Let user filter strictness in nextflow params?
        # The args.min_score is passed from main.nf (default 0.6).
        # We should respect it.
        
        if passes_score: # Relaxed p-value check for now?
             name = f"Reg_G{best['unique']}_S{best['strand_cons']:.1f}_P{best['p_value']:.3f}"
             # Determine overall strand of region based on majority?
             # For BED output, strand is useful.
             # If consistency is High on Inverted, set region strand to - ?
             
             region_strand = "."
             # Count +'s and -'s in hits
             plus_cnt = sum(1 for h in best['cluster'] if h['strand'] == '+')
             minus_cnt = len(best['cluster']) - plus_cnt
             if minus_cnt > plus_cnt: region_strand = "-"
             else: region_strand = "+"
             
             f_out.write(f"{best['chrom']}\t{best['start']}\t{best['end']}\t{name}\t{best['score']:.2f}\t{region_strand}\n")

if __name__ == "__main__":
    main()
