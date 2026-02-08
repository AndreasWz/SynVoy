#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
import concurrent.futures
import uuid
import sys
import json
import logging
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Use our own sequence utilities (no BioPython dependency)
try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )
    from annotate_goi_exons import annotate_exons_from_hit_list
except ImportError:
    # Fallback if not in path - add bin directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )
    from annotate_goi_exons import annotate_exons_from_hit_list

# Import fragment utilities if available
try:
    from fragment_query import generate_fragments, parse_fragment_id, merge_fragment_hits
    FRAGMENT_SUPPORT = True
except ImportError:
    FRAGMENT_SUPPORT = False

def run_command(cmd):
    subprocess.check_call(cmd)

def normalize_coordinates(start: int, end: int) -> Tuple[int, int]:
    return min(start, end), max(start, end)

def parse_hits(hits_file: str, min_identity: float, min_length: int, evalue_thresh: float) -> List[Dict[str, Any]]:
    """
    Parse MMseqs2 hits and return a list of hit dictionaries.
    Filters by basic quality metrics.
    Preserves qstart/qend (query protein positions) and strand for exon annotation.
    """
    hits = []
    if not os.path.exists(hits_file):
        return hits
        
    with open(hits_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            try:
                # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
                # 0      1       2       3       4         5        6       7     8       9     10      11
                if len(parts) < 11: continue
                
                pident = float(parts[2])
                alnlen = int(parts[3])
                evalue = float(parts[10])
                
                if (evalue <= evalue_thresh and 
                    pident >= min_identity and 
                    alnlen >= min_length):
                    
                    t_start = int(parts[8])
                    t_end = int(parts[9])
                    start, end = normalize_coordinates(t_start, t_end)
                    # Convert 1-based mmseqs/BLAST coordinates to 0-based half-open
                    # for Python slicing: start-1 becomes 0-based, end stays (exclusive)
                    start -= 1
                    strand = '+' if t_start <= t_end else '-'
                    
                    q_start = int(parts[6])
                    q_end = int(parts[7])
                    
                    hits.append({
                        'query': parts[0],
                        'target': parts[1], # Chromosome/Scaffold
                        'chrom': parts[1],
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'qstart': min(q_start, q_end),
                        'qend': max(q_start, q_end),
                        'evalue': evalue,
                        'pident': pident,
                        'alnlen': alnlen
                    })
            except Exception as e:
                continue
    return hits

def create_locus_object(query_id, hits):
    chrom = hits[0]['chrom']
    start = min(h['start'] for h in hits)
    end = max(h['end'] for h in hits)
    return {
        'query': query_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'hits': hits
    }

def extract_base_gene_id(query_id: str) -> str:
    """
    Extract the base gene ID from a query ID that may contain exon info.
    
    Handles formats:
    - "gene-LOC726866" -> "gene-LOC726866"
    - "gene-LOC726866|exon_1" -> "gene-LOC726866"
    - "gene-LOC726866|var1" -> "gene-LOC726866"
    - "gene-LOC726866|GCA_xxx_MP000001" -> "gene-LOC726866"
    """
    # Split on | and take the first part
    parts = query_id.split('|')
    base_id = parts[0]
    
    # Also handle cases where gene ID itself contains underscores but not exon info
    # The exon suffix is always "|exon_N"
    
    return base_id

def identify_best_synteny_block(hits, max_intron=20000, cluster_dist=50000):
    """
    Identify the best synteny block from hits.
    
    Key Logic:
    - Group hits by query gene (handling multi-exon genes AND exon-level queries)
    - Cluster loci that are close together (likely same gene/region)
    - Score blocks by number of unique flanking genes found
    
    Updated to handle exon-level queries like "gene-LOC726866|exon_1"
    
    Args:
        hits: List of hit dictionaries
        max_intron: Maximum distance between exons of same gene (bp)
        cluster_dist: Maximum distance to cluster genes into synteny block (bp)
    
    Returns:
        Dictionary with best synteny block info or None
    """
    if not hits:
        return None
        
    # --- Step 1: Group hits by Query Gene (Base ID) ---
    # Multiple hits from same query = different exons or duplicates
    # With exon_mode, query may be "gene-LOC726866|exon_1" -> base is "gene-LOC726866"
    hits_by_query = defaultdict(list)
    for h in hits:
        # Extract base query ID (handles both |var and |exon_ suffixes)
        base_query = extract_base_gene_id(h['query'])
        hits_by_query[base_query].append(h)
        
    # --- Step 2: Define Gene Loci per Query ---
    # Each query gene may have multiple loci (paralogs/duplications)
    # Hits close together (<max_intron) = same locus (multi-exon gene)
    all_loci = []
    for query_id, q_hits in hits_by_query.items():
        # Sort hits by genomic position
        q_hits.sort(key=lambda x: (x['chrom'], x['start']))
        
        current_locus_hits = []
        for h in q_hits:
            if not current_locus_hits:
                current_locus_hits.append(h)
                continue
            
            last_hit = current_locus_hits[-1]
            
            # Same chromosome and close enough = same locus (exons)
            if (h['chrom'] == last_hit['chrom'] and 
                h['start'] - last_hit['end'] < max_intron):
                current_locus_hits.append(h)
            else:
                # Start new locus
                all_loci.append(create_locus_object(query_id, current_locus_hits))
                current_locus_hits = [h]
        
        # Don't forget last locus
        if current_locus_hits:
            all_loci.append(create_locus_object(query_id, current_locus_hits))

    # --- Step 3: Cluster Loci into Synteny Blocks ---
    # Loci from different genes that are close = synteny block
    all_loci.sort(key=lambda x: (x['chrom'], x['start']))
    if not all_loci: 
        return None
        
    synteny_blocks = []
    current_block = [all_loci[0]]
    
    for locus in all_loci[1:]:
        last_locus = current_block[-1]
        
        # Same chromosome and within clustering distance
        if (locus['chrom'] == last_locus['chrom'] and 
            locus['start'] - last_locus['end'] < cluster_dist):
            current_block.append(locus)
        else:
            synteny_blocks.append(current_block)
            current_block = [locus]
    synteny_blocks.append(current_block)
    
    # --- Step 4: Score and Select Best Block ---
    best_block = None
    best_gene_count = -1
    
    for block in synteny_blocks:
        # Score = Count of Unique Query Genes (base IDs only)
        unique_genes = set(extract_base_gene_id(l['query']) for l in block)
        gene_count = len(unique_genes)
        
        # Prefer block with most unique genes
        if gene_count > best_gene_count:
            best_gene_count = gene_count
            best_block = block
        # Tie-breaker: more loci = better (more complete gene models)
        elif gene_count == best_gene_count:
            if len(block) > len(best_block):
                best_block = block

    if not best_block: 
        return None

    # Compile block metadata
    chrom = best_block[0]['chrom']
    start = min(l['start'] for l in best_block)
    end = max(l['end'] for l in best_block)
    genes_list = list(set(extract_base_gene_id(l['query']) for l in best_block))
    
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'genes_count': best_gene_count,
        'loci_count': len(best_block),
        'genes': genes_list
    }

def calculate_adaptive_padding(hits: List[Dict[str, Any]], best_region: Dict[str, Any], default: int = 100000) -> int:
    """
    Calculate region padding based on gene spacing in hits.
    Returns padding distance in base pairs.
    """
    # Filter hits to the best region's chromosome
    region_hits = [h for h in hits if h['chrom'] == best_region['chrom']]
    
    if len(region_hits) < 2:
        return default
    
    # Sort by position
    sorted_hits = sorted(region_hits, key=lambda h: h['start'])
    
    # Calculate inter-gene gaps
    gaps = []
    for i in range(len(sorted_hits) - 1):
        gap = sorted_hits[i+1]['start'] - sorted_hits[i]['end']
        if gap > 0:  # Only positive gaps
            gaps.append(gap)
    
    if not gaps:
        return default
    
    # Average gap * 2 (to cover one gene on each side)
    avg_gap = sum(gaps) / len(gaps)
    adaptive_padding = int(avg_gap * 2)
    
    # Clamp to reasonable range
    final_padding = max(50000, min(200000, adaptive_padding))
    
    return final_padding


def estimate_cluster_dist(genome_file: str, gff_file: Optional[str] = None, default_dist: int = 50000) -> int:
    """
    Estimate gene density to adjust cluster_dist intelligently.
    
    Strategy:
    1. If GFF provided: Calculate actual inter-gene distances
    2. Else: Use genome size heuristic (improved)
    3. Return 2-3x median inter-gene distance as cluster threshold
    """
    
    # Method 1: Use GFF if available (most accurate)
    if gff_file and os.path.exists(gff_file) and gff_file != "NO_GFF":
        try:
            genes = parse_gff(gff_file)
            if len(genes) > 10:  # Need reasonable sample size
                # Sort by chromosome and position
                by_chrom = defaultdict(list)
                for gene in genes:
                    by_chrom[gene['chrom']].append(gene['start'])
                
                # Calculate inter-gene distances per chromosome
                all_distances = []
                for chrom, positions in by_chrom.items():
                    sorted_pos = sorted(positions)
                    for i in range(len(sorted_pos) - 1):
                        dist = sorted_pos[i+1] - sorted_pos[i]
                        if dist > 0:  # Skip overlapping genes
                            all_distances.append(dist)
                
                if all_distances:
                    # Use median distance * 2.5 as clustering threshold
                    all_distances.sort()
                    median_dist = all_distances[len(all_distances) // 2]
                    cluster_dist = int(median_dist * 2.5)
                    # Clamp to reasonable range
                    cluster_dist = max(10000, min(200000, cluster_dist))
                    print(f"Estimated cluster distance from GFF: {cluster_dist} bp "
                          f"(median inter-gene: {median_dist} bp)", file=sys.stderr)
                    return cluster_dist
        except Exception as e:
            print(f"Warning: Could not parse GFF for gene density: {e}", file=sys.stderr)
    
    # Method 2: Improved genome size heuristic
    try:
        size = os.path.getsize(genome_file)
        
        # More refined heuristics based on typical genomes
        if size < 5_000_000:  # < 5MB: Bacteria/Archaea
            return 15000  # Dense gene packing
        elif size < 20_000_000:  # 5-20MB: Large bacteria, fungi
            return 25000
        elif size < 100_000_000:  # 20-100MB: Small eukaryotes
            return 40000
        elif size < 500_000_000:  # 100-500MB: Insects, small vertebrates
            return 70000
        elif size < 2_000_000_000:  # 0.5-2GB: Mammals, birds
            return 100000
        else:  # > 2GB: Plants, large genomes
            return 150000
    except:
        pass
    
    return default_dist

def run_augmented_search(region_fasta: str, goi_queries: List[Dict[str, str]], 
                        genome_name: str, args, unique_id: str, threads: int) -> List[Dict[str, Any]]:
    """
    Run augmented search (MMseqs2 + Smith-Waterman) for GOI queries.
    
    Uses both methods for maximum sensitivity:
    1. MMseqs2 with query fragments (fast, good for similar sequences)
    2. Smith-Waterman via parasail/ssearch36 (slower, better for divergent sequences)
    
    Args:
        region_fasta: Path to extracted region FASTA
        goi_queries: List of GOI query dicts with 'id' and 'seq'
        genome_name: Name of current genome
        args: Command line arguments
        unique_id: Unique ID for temp files
        threads: Number of threads to use
        
    Returns:
        List of hit dictionaries combining MMseqs2 and Smith-Waterman results
    """
    all_hits = []
    
    try:
        # Generate variants for each GOI query
        variants_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_variants.faa"
        query_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_goi_query.faa"
        all_variants = []
        
        # Write full query sequences for Smith-Waterman
        write_fasta([(q['id'], q['seq']) for q in goi_queries], query_fasta)
        
        if not FRAGMENT_SUPPORT:
            print(f"[{genome_name}] Warning: fragment_query module not available, using full sequences only", flush=True)
            # Just use the original sequences
            all_variants = [(q['id'], q['seq']) for q in goi_queries]
        else:
            for query in goi_queries:
                # Generate fragments (halves, thirds, quarters)
                fragments = generate_fragments(query['seq'], query['id'], min_size=15)
                all_variants.extend([(f[0], f[1]) for f in fragments])
        
        write_fasta(all_variants, variants_fasta)
        
        # ========== 1. MMseqs2 Search ==========
        aug_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_hits.m8"
        aug_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_mmseqs"
        
        os.makedirs(aug_tmp_dir, exist_ok=True)
        
        # CRITICAL: Use VERY relaxed e-value for augmented search (100x more permissive)
        # This ensures divergent hits are not filtered out by MMseqs2 before we can parse them
        relaxed_evalue = min(10.0, args.evalue * 1000)  # Much more permissive, cap at 10
        
        subprocess.run([
            "mmseqs", "easy-search",
            variants_fasta, region_fasta, aug_hits_file, aug_tmp_dir,
            "--search-type", "2",  # Protein search
            "--threads", str(threads),
            "-s", str(args.mmseqs_sens),  # Use same sensitivity as main search
            "-e", str(relaxed_evalue),  # RELAXED e-value to capture divergent hits
            "--min-seq-id", "0.0",  # NO identity filtering at search time - we filter in parse_hits
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        ], check=True, stderr=subprocess.DEVNULL)
        
        # Parse hits - use relaxed thresholds for augmented search
        relaxed_identity = max(25.0, args.min_identity * 0.6)  # 60% of normal threshold
        relaxed_length = max(15, args.min_length // 2)  # Half of normal length
        
        mmseqs_hits = parse_hits(aug_hits_file, relaxed_identity, relaxed_length, args.evalue * 10)
        if mmseqs_hits:
            print(f"[{genome_name}] MMseqs2 augmented search found {len(mmseqs_hits)} hits.", flush=True)
            all_hits.extend(mmseqs_hits)
        
        # ========== 2. Smith-Waterman Search ==========
        # Use Smith-Waterman for very divergent sequences (more sensitive than MMseqs2)
        sw_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_sw_hits.m8"
        
        try:
            # Use smith_waterman_search.py script
            sw_cmd = [
                "python3", os.path.join(os.path.dirname(__file__), "smith_waterman_search.py"),
                "--query", query_fasta,
                "--target", region_fasta,
                "--output", sw_hits_file,
                "--min_score", "30",
                "--min_identity", "15.0",  # Very relaxed for divergent sequences
                "--threads", str(threads)
            ]
            
            result = subprocess.run(sw_cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0 and os.path.exists(sw_hits_file):
                # Parse Smith-Waterman hits (BLAST m8 format)
                sw_hits = parse_hits(sw_hits_file, 15.0, 10, 100.0)  # Very relaxed thresholds
                if sw_hits:
                    print(f"[{genome_name}] Smith-Waterman found {len(sw_hits)} additional hits.", flush=True)
                    # Mark hits as from Smith-Waterman
                    for hit in sw_hits:
                        hit['method'] = 'smith_waterman'
                    all_hits.extend(sw_hits)
            elif result.stderr:
                print(f"[{genome_name}] Smith-Waterman warning: {result.stderr[:200]}", flush=True)
                
        except subprocess.TimeoutExpired:
            print(f"[{genome_name}] Smith-Waterman timed out after 5 minutes, using MMseqs2 only.", flush=True)
        except FileNotFoundError:
            print(f"[{genome_name}] Smith-Waterman script not found, using MMseqs2 only.", flush=True)
        except Exception as sw_err:
            print(f"[{genome_name}] Smith-Waterman failed: {sw_err}, using MMseqs2 only.", flush=True)
        
        # Clean up temp files
        for f in [variants_fasta, query_fasta, aug_hits_file, sw_hits_file]:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(aug_tmp_dir):
            shutil.rmtree(aug_tmp_dir, ignore_errors=True)
        
        # Deduplicate hits by position (prefer higher identity)
        if all_hits:
            deduped = {}
            for hit in all_hits:
                key = (hit.get('query', ''), hit.get('start', 0) // 100, hit.get('end', 0) // 100)
                if key not in deduped or hit.get('pident', 0) > deduped[key].get('pident', 0):
                    deduped[key] = hit
            all_hits = list(deduped.values())
            print(f"[{genome_name}] Combined augmented search: {len(all_hits)} unique hits.", flush=True)
        
        return all_hits
        
    except Exception as e:
        print(f"[{genome_name}] Augmented search failed: {e}", flush=True)
        return []

def batch_rbh_check(candidates, home_db, cand_map, threads=1, evalue=1e-5, min_coverage=0.5):
    """
    Perform Reciprocal Best Hit check for multiple candidates at once.
    
    Enhanced validation:
    1. RBH to home genome (traditional)
    2. Coverage check: alignment must cover >50% of both query and target
    3. Identity must be reasonable (>25%)
    
    candidates: list of dicts with 'id' and 'seq' keys
    """
    if not candidates: return []
    
    unique_id = uuid.uuid4().hex
    query_fasta = f"batch_candidates_{unique_id}.fasta"
    rbh_out = f"batch_rbh_{unique_id}.m8"
    tmp_subdir = f"tmp_rbh_batch_{unique_id}"
    
    valid_ids = set()

    try:
        # Write FASTA using our utility
        records = [(c['id'], c['seq']) for c in candidates]
        write_fasta(records, query_fasta)
            
        db_path = home_db
        if os.path.isdir(home_db):
            db_path = os.path.join(home_db, "db")
            
        cmd = [
            "mmseqs", "easy-search",
            query_fasta, db_path, rbh_out, tmp_subdir,
            "-e", str(evalue),
            "--format-output", "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen",
            "--max-seqs", "1", # Top hit only
            "--threads", str(threads)
        ]
        
        subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
        
        if os.path.exists(rbh_out):
            with open(rbh_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 9: continue
                    
                    cand_id = parts[0]
                    target_id = parts[1]
                    pident = float(parts[2])
                    qcov = float(parts[3]) if len(parts) > 3 else 100
                    tcov = float(parts[4]) if len(parts) > 4 else 100
                    
                    # Enhanced validation
                    if cand_id not in cand_map:
                        continue
                    
                    parent = cand_map[cand_id]
                    parent_base = extract_base_gene_id(parent).strip()
                    target_base = extract_base_gene_id(target_id).strip()
                    
                    # Check 1: ID matching (exact or close)
                    ids_match = (parent_base == target_base or 
                                parent == target_id or 
                                target_id == parent or
                                parent_base in target_base or
                                target_base in parent_base)
                    
                    # Check 2: Coverage (both query and target must be well-covered)
                    coverage_ok = (qcov >= min_coverage * 100 and 
                                  tcov >= min_coverage * 100)
                    
                    # Check 3: Identity must be reasonable
                    identity_ok = pident >= 25.0
                    
                    if ids_match and coverage_ok and identity_ok:
                        valid_ids.add(cand_id)
                    elif ids_match and not coverage_ok:
                        print(f"RBH: {cand_id} matches {target_id} but low coverage "
                              f"(qcov={qcov:.0f}%, tcov={tcov:.0f}%). Likely fragment/paralog.",
                              file=sys.stderr)
                    elif ids_match and not identity_ok:
                        logger.debug(f"RBH: {cand_id} matches {target_id} but very low identity "
                              f"({pident:.1f}%). Possible pseudogene.")
                              
    except Exception as e:
        logger.error(f"RBH check failed: {e}")
        return set()
    finally:
        # Cleanup
        if os.path.exists(query_fasta): os.remove(query_fasta)
        if os.path.exists(rbh_out): os.remove(rbh_out)
        if os.path.exists(tmp_subdir):
             shutil.rmtree(tmp_subdir, ignore_errors=True)
                          
    return valid_ids

def process_single_genome(genome_path, db_path, args, home_db_dir, prefix, threads_per_job):
    """
    Worker function to search a single genome.
    Returns: (genome_name, list_of_new_genes)
    """
    genome_name = os.path.basename(genome_path)
    if not os.path.exists(genome_path):
        logger.warning(f"[{genome_name}] Genome file not found. Skipping.")
        return genome_name, []
    
    # Create unique temp space
    unique_id = uuid.uuid4().hex
    hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
    tmp_dir = f"{args.output_dir}/tmp_mmseqs_{unique_id}_{genome_name}"
    
    new_genes = []
    
    try:
        # Auto-param
        c_dist = args.cluster_dist
        if c_dist <= 0:
            c_dist = estimate_cluster_dist(genome_path)
            
        # 1. Search (MMseqs)
        subprocess.run([
            "mmseqs", "easy-search",
            db_path, genome_path, hits_file, tmp_dir,
            "--search-type", "2", 
            "--threads", str(threads_per_job),
            "-s", str(args.mmseqs_sens),  # Use configurable sensitivity
            "-e", str(args.evalue),
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        ], check=True, stderr=subprocess.DEVNULL)

        hits = parse_hits(hits_file, args.min_identity, args.min_length, args.evalue)
        if not hits:
            logger.info(f"[{genome_name}] No hits found in MMseqs output.")
            return genome_name, []
            
        logger.info(f"[{genome_name}] Parsed {len(hits)} hits.")

        # 2. Identify Synteny
        best_region = identify_best_synteny_block(hits, cluster_dist=c_dist)
        
        if not best_region:
            logger.info(f"[{genome_name}] No valid syntenic region found.")
            return genome_name, []
            
        logger.info(f"[{genome_name}] Found syntenic region: {best_region['chrom']}:{best_region['start']}-{best_region['end']} with {len(best_region['genes'])} genes.")
            # print(f"[{genome_name}] Found Region: {best_region['chrom']}:{best_region['start']}-{best_region['end']}")
            
        # 3. Extract Region - using our FASTA parser
        genome_seqs = load_genome(genome_path)
        chrom = best_region['chrom']
        
        if chrom in genome_seqs:
            slen = len(genome_seqs[chrom])
            
            # ADAPTIVE PADDING - CRITICAL FIX: Increased from 20kb to 150kb default
            # This ensures query gene is captured even if distant from flanking genes
            padding = calculate_adaptive_padding(hits, best_region, default=150000)
            
            w_start = max(0, best_region['start'] - padding)
            w_end = min(slen, best_region['end'] + padding)
            subseq = genome_seqs[chrom][w_start:w_end]
            
            temp_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_reg.fasta"
            
            # Write region FASTA
            write_fasta([("region_seq", subseq)], temp_fa)
            
            # 4. Prepare Hits & Query Sequences
            relevant_hits = [h for h in hits if h['chrom'] == chrom]
            unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)
            
            # CRITICAL FIX: Always include GOI queries (marked with GOI_ prefix)
            # This ensures the query gene of interest is searched iteratively
            print(f"[{genome_name}] Scanning database for GOI queries...", flush=True)
            
            # Extract sequences from DB - using our parser
            found_queries = []  # List of {'id': ..., 'seq': ...}
            db_sequences = {}
            goi_queries = set()  # Track GOI queries
            
            try:
                # Parse DB once into dict
                for header, clean_id, seq in parse_fasta(db_path):
                    db_sequences[clean_id] = {'id': clean_id, 'seq': seq, 'header': header}
                    # CRITICAL FIX: Also index by base gene ID (stripping |exon_N suffix)
                    # This matches how unique_queries is built via extract_base_gene_id
                    base = extract_base_gene_id(clean_id)
                    if base != clean_id and base not in db_sequences:
                        db_sequences[base] = {'id': clean_id, 'seq': seq, 'header': header}
                    
                    # Track GOI queries (these MUST always be searched)
                    if 'GOI_' in clean_id or clean_id.startswith('GOI_'):
                        goi_queries.add(clean_id)
                
                # CRITICAL: Force include all GOI queries
                # But prefer full-length GOI over fragments for miniprot
                full_length_goi = set()
                fragment_goi = set()
                for goi_id in goi_queries:
                    if '|frag_' in goi_id:
                        fragment_goi.add(goi_id)
                    else:
                        full_length_goi.add(goi_id)
                
                # Use full-length if available, otherwise use all
                if full_length_goi:
                    unique_queries.update(full_length_goi)
                    print(f"[{genome_name}] Using {len(full_length_goi)} full-length GOI queries (excluding {len(fragment_goi)} fragments).", flush=True)
                else:
                    unique_queries.update(goi_queries)
                    print(f"[{genome_name}] Using all {len(goi_queries)} GOI queries (fragments only).", flush=True)
                
                for query_id in unique_queries:
                    if query_id in db_sequences:
                        found_queries.append(db_sequences[query_id])
                    else:
                        # Fallback: try base gene ID (consistent with database indexing)
                        base = extract_base_gene_id(query_id)
                        if base in db_sequences:
                            found_queries.append(db_sequences[base])
            except Exception as dex:
                print(f"[{genome_name}] Warning: DB parsing failed: {dex}")

            if found_queries:
                print(f"[{genome_name}] Annotating {len(found_queries)} queries using MMseqs2/Smith-Waterman...", flush=True)
                query_mini_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_query.faa"
                # Write query FASTA
                write_fasta([(q['id'], q['seq']) for q in found_queries], query_mini_fa)
                
                try:
                    # 5. EXON-AWARE ANNOTATION
                    # Uses same logic as annotate_goi_exons.py: each hit is a
                    # candidate exon, with splice-site detection, start/stop
                    # codon checks, coverage-gap search, and deduplication.
                    
                    annotated_records_raw = []  # List of {'id': ..., 'seq': ...}
                    valid_gff_lines = []
                    clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')
                    
                    # 5a. Run augmented search (MMseqs2 + Smith-Waterman)
                    print(f"[{genome_name}] Running augmented search (MMseqs2 + SW) for {len(found_queries)} queries...", flush=True)
                    
                    augmented_hits_mmseqs = run_augmented_search(
                        temp_fa, found_queries, genome_name, 
                        args, unique_id, threads_per_job
                    )
                    
                    # Collect search hits; fall back to original MMseqs2 hits
                    all_search_hits = list(augmented_hits_mmseqs) if augmented_hits_mmseqs else []
                    if not all_search_hits:
                        print(f"[{genome_name}] No augmented hits. Using original MMseqs2 hits...", flush=True)
                        for hit in relevant_hits:
                            region_hit = dict(hit)
                            region_hit['start'] = max(0, hit['start'] - w_start)
                            region_hit['end'] = max(0, hit['end'] - w_start)
                            all_search_hits.append(region_hit)
                    
                    if all_search_hits:
                        print(f"[{genome_name}] Total search hits: {len(all_search_hits)}", flush=True)
                        
                        # Group hits by parent gene, converting to exon format
                        hits_by_gene = defaultdict(list)
                        for hit in all_search_hits:
                            query_id = hit['query']
                            parent = query_id.split('|frag_')[0] if '|frag_' in query_id else query_id
                            parent = extract_base_gene_id(parent)
                            
                            # Map fragment qstart/qend to full protein coordinates
                            frag_offset = 0
                            if '|frag_' in query_id and FRAGMENT_SUPPORT:
                                try:
                                    frag_info = parse_fragment_id(query_id)
                                    frag_offset = frag_info['start'] - 1
                                except Exception:
                                    frag_offset = 0
                            
                            hits_by_gene[parent].append({
                                'qstart': hit.get('qstart', 1) + frag_offset,
                                'qend': hit.get('qend', 100) + frag_offset,
                                'gstart': hit['start'],
                                'gend': hit['end'],
                                'evalue': hit.get('evalue', 1),
                                'pident': hit.get('pident', 0),
                                'alnlen': hit.get('alnlen', 0),
                                'strand': hit.get('strand', '+'),
                                'chrom': chrom
                            })
                        
                        # 5b. Exon-aware annotation for each gene
                        for parent_id, gene_hits in hits_by_gene.items():
                            # Find query protein for this gene
                            parent_query_seq = None
                            for q in found_queries:
                                if extract_base_gene_id(q['id']) == parent_id:
                                    parent_query_seq = q['seq']
                                    break
                            if not parent_query_seq:
                                for q in found_queries:
                                    if parent_id in q['id'] or q['id'].startswith(parent_id):
                                        parent_query_seq = q['seq']
                                        break
                            if not parent_query_seq:
                                continue
                            
                            # Run exon-aware annotation (same logic as annotate_goi_exons)
                            try:
                                # Check for tandem duplications first (GOI genes)
                                is_goi = parent_id.startswith('GOI_') or 'GOI_' in parent_id
                                is_tandem = False
                                if is_goi:
                                    from annotate_goi_exons import detect_tandem_duplications
                                    is_tandem, tandem_copies = detect_tandem_duplications(
                                        gene_hits, parent_query_seq, subseq, chrom
                                    )
                                
                                if is_tandem and tandem_copies:
                                    exons = tandem_copies
                                    print(f"[{genome_name}] TANDEM: {len(exons)} copies of {parent_id}", flush=True)
                                else:
                                    exons, _ = annotate_exons_from_hit_list(
                                        gene_hits, parent_query_seq, subseq, chrom,
                                        search_missing=True
                                    )
                            except Exception as ann_err:
                                print(f"[{genome_name}] Exon annotation failed for {parent_id}: {ann_err}", flush=True)
                                exons = []
                            
                            if exons:
                                # Check if this is tandem (copies, not exons)
                                is_tandem_result = any(e.get('id', '').startswith('GOI_copy_') for e in exons)
                                
                                if is_tandem_result:
                                    # TANDEM: output each copy as separate record
                                    for copy in exons:
                                        global_start = w_start + copy['gstart'] + 1
                                        global_end = w_start + copy['gend']
                                        strand = copy.get('strand', '+')
                                        
                                        copy_id = f"{copy['id']}|{clean_gname}"
                                        annotated_records_raw.append({
                                            'id': copy_id,
                                            'seq': copy['seq'],
                                            'description': (f"coords:{global_start}-{global_end} "
                                                           f"parent:{parent_id} tandem_copy "
                                                           f"identity:{copy.get('pident', 0):.1f}")
                                        })
                                        
                                        valid_gff_lines.append(
                                            f"{chrom}\ttandem_copy\tgene\t{global_start}\t{global_end}\t"
                                            f"{copy.get('pident', 0):.1f}\t{strand}\t.\t"
                                            f"ID={copy_id};Name={copy['id']};"
                                            f"SynTerra_Parent={parent_id};SynTerra_ID={copy_id};"
                                            f"Identity={copy.get('pident', 0):.1f};Type=tandem_copy"
                                        )
                                        
                                        print(f"[{genome_name}]   {copy['id']}: {len(copy['seq'])} aa, "
                                              f"{copy.get('pident', 0):.1f}% id (tandem copy)", flush=True)
                                else:
                                    # EXONS: concatenate into single protein
                                    exons.sort(key=lambda e: e.get('qstart', 0))
                                    exon_protein = ''.join(e['seq'] for e in exons)
                                    strand = exons[0].get('strand', '+')
                                    avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                                
                                    # Convert 0-based coords to 1-based for GFF output
                                    global_start = w_start + min(e['gstart'] for e in exons) + 1
                                    global_end = w_start + max(e['gend'] for e in exons)
                                    
                                    new_id = f"{parent_id}|{clean_gname}_exon_ann"
                                    annotated_records_raw.append({
                                        'id': new_id,
                                        'seq': exon_protein,
                                        'description': (f"coords:{global_start}-{global_end} "
                                                       f"parent:{parent_id} exons:{len(exons)} "
                                                       f"identity:{avg_pident:.1f}")
                                    })
                                    
                                    # GFF: mRNA line
                                    valid_gff_lines.append(
                                        f"{chrom}\texon_annotation\tmRNA\t{global_start}\t{global_end}\t"
                                        f"{avg_pident:.1f}\t{strand}\t.\t"
                                        f"ID={new_id};Name={parent_id};"
                                        f"SynTerra_Parent={parent_id};SynTerra_ID={new_id};"
                                        f"Identity={avg_pident:.1f};Exons={len(exons)}"
                                    )
                                    
                                    # GFF: CDS lines per exon (with splice site metadata)
                                    for eidx, exon in enumerate(exons, 1):
                                        exon_gs = w_start + exon['gstart'] + 1  # 1-based for GFF
                                        exon_ge = w_start + exon['gend']        # 0-based excl = 1-based incl
                                        attrs = f"ID={new_id}_CDS{eidx};Parent={new_id}"
                                        if exon.get('splice_acceptor'):
                                            attrs += f";SpliceAcceptor={exon['splice_acceptor']}"
                                        if exon.get('splice_donor'):
                                            attrs += f";SpliceDonor={exon['splice_donor']}"
                                        if exon.get('has_start_codon'):
                                            attrs += ";StartCodon=ATG"
                                        if exon.get('has_stop_codon'):
                                            attrs += ";StopCodon=yes"
                                        valid_gff_lines.append(
                                            f"{chrom}\texon_annotation\tCDS\t{exon_gs}\t{exon_ge}\t"
                                            f".\t{strand}\t0\t{attrs}"
                                        )
                                    
                                    print(f"[{genome_name}]   {parent_id}: {len(exons)} exon(s), "
                                          f"{len(exon_protein)} aa, {avg_pident:.1f}% id", flush=True)
                            else:
                                # Fallback: translate best raw hit directly
                                best_hit = min(gene_hits, key=lambda h: h.get('evalue', 1))
                                g_s, g_e = best_hit['gstart'], best_hit['gend']
                                strand = best_hit.get('strand', '+')
                                
                                region_dna = subseq[g_s:g_e]
                                if strand == '-':
                                    region_dna = reverse_complement(region_dna)
                                region_dna = region_dna[:len(region_dna) - len(region_dna) % 3]
                                
                                if len(region_dna) >= 9:
                                    hit_protein = translate(region_dna).replace('*', '')
                                    if hit_protein:
                                        nt_s = w_start + g_s + 1  # 1-based for GFF
                                        nt_e = w_start + g_e      # 0-based excl = 1-based incl
                                        new_id = f"{parent_id}|{clean_gname}_raw"
                                        annotated_records_raw.append({
                                            'id': new_id,
                                            'seq': hit_protein,
                                            'description': f"coords:{nt_s}-{nt_e} parent:{parent_id} identity:{best_hit.get('pident', 0):.1f}"
                                        })
                                        valid_gff_lines.append(
                                            f"{chrom}\traw_hit\tmRNA\t{nt_s}\t{nt_e}\t"
                                            f"{best_hit.get('pident', 0):.1f}\t{strand}\t.\t"
                                            f"ID={new_id};Name={parent_id};"
                                            f"SynTerra_Parent={parent_id};SynTerra_ID={new_id}"
                                        )
                                        valid_gff_lines.append(
                                            f"{chrom}\traw_hit\tCDS\t{nt_s}\t{nt_e}\t"
                                            f".\t{strand}\t0\t"
                                            f"ID={new_id}_CDS1;Parent={new_id}"
                                        )
                        
                        print(f"[{genome_name}] Exon-aware annotation: {len(annotated_records_raw)} genes.", flush=True)
                    
                    # 6. Deduplicate: remove entries with identical coordinates
                    # (e.g., GOI_P01501 and GOI_Melt often annotate the same locus)
                    seen_coords = {}
                    deduped_records = []
                    deduped_gff = []
                    for rec in annotated_records_raw:
                        desc = rec.get('description', '')
                        coords_key = None
                        if 'coords:' in desc:
                            coords_key = desc.split('coords:')[1].split(' ')[0]
                        if coords_key and coords_key in seen_coords:
                            print(f"[{genome_name}] Removing duplicate: {rec['id']} (same coords as {seen_coords[coords_key]})", flush=True)
                            # Remove corresponding GFF lines
                            rid = rec['id']
                            valid_gff_lines = [g for g in valid_gff_lines if f"ID={rid}" not in g and f"Parent={rid}" not in g]
                            continue
                        if coords_key:
                            seen_coords[coords_key] = rec['id']
                        deduped_records.append(rec)
                    annotated_records_raw = deduped_records
                    
                    new_genes = annotated_records_raw
                    print(f"[{genome_name}] Keeping {len(new_genes)} candidates.", flush=True)
                    
                    # Write GFF, FASTA, and Homology TSV
                    if valid_gff_lines:
                        gff_out = f"{args.output_dir}/regions/{genome_name}.gff"
                        faa_out = f"{args.output_dir}/regions/{genome_name}.faa"
                        tsv_out = f"{args.output_dir}/regions/{genome_name}.homology.tsv"
                        
                        with open(gff_out, 'w') as gf:
                            gf.write("##gff-version 3\n")
                            for gl in valid_gff_lines:
                                gf.write(gl + "\n")

                        write_fasta([(g['id'], g['seq']) for g in new_genes], faa_out)
                            
                        with open(tsv_out, 'w') as tf:
                            for rec in new_genes:
                                parent = extract_base_gene_id(rec['id'])
                                tf.write(f"{rec['id']}\t{parent}\n")

                except Exception as ann_err:
                    print(f"[{genome_name}] Error during annotation: {ann_err}")
                    new_genes = [] # Fail safe
                    
                # Cleanup temp
                if os.path.exists(temp_fa): os.remove(temp_fa)
                try:
                    if query_mini_fa and os.path.exists(query_mini_fa): os.remove(query_mini_fa)
                except NameError:
                    pass
            else:
                print(f"[{genome_name}] Warning: Could not find query sequences for relevant hits.")
        else:
            print(f"[{genome_name}] Warning: Chromosome {chrom} not found in genome file.")

    except Exception as e:
        print(f"[{genome_name}] Error processing: {e}")
    finally:
        # Cleanup mmseqs tmp dir
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)
            
    return genome_name, new_genes

def main():
    parser = argparse.ArgumentParser(description="Iterative Genome Search Runner (Wavefront Parallel)")
    parser.add_argument("--initial_db", required=True)
    parser.add_argument("--sorted_genomes", required=True, 
                        help="Tab-separated file: genome_path\\tdistance. Genomes sorted by distance.")
    parser.add_argument("--genomes_dir", help="Directory containing genome files (if paths in sorted_genomes are relative)")
    parser.add_argument("--home_db_dir", help="Home Proteome MMseqs DB for RBH")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min_identity", type=float, default=40.0)
    parser.add_argument("--min_length", type=int, default=50)
    parser.add_argument("--threads", type=int, default=4, help="Total threads available for parallel processing")
    parser.add_argument("--cluster_dist", type=int, default=-1, help="Auto-detect if -1")
    parser.add_argument("--mmseqs_sens", type=float, default=7.5, help="MMseqs2 sensitivity (higher = more sensitive but slower)")
    parser.add_argument("--prefix", default="", help="Prefix for output files (e.g. locus ID)")
    parser.add_argument("--resume", action="store_true", help="Resume from previous checkpoint if available")
    
    args = parser.parse_args()
    
    # INPUT VALIDATION
    # 1. Validate required files exist
    if not os.path.exists(args.initial_db):
        logger.error(f"Initial database file not found: {args.initial_db}")
        sys.exit(1)
    
    if not os.path.exists(args.sorted_genomes):
        logger.error(f"Sorted genomes file not found: {args.sorted_genomes}")
        sys.exit(1)
    
    # 2. Validate initial_db is not empty
    if os.path.getsize(args.initial_db) == 0:
        logger.error("Initial database file is empty")
        sys.exit(1)
    
    # 3. Validate parameters are in valid ranges
    if args.min_identity < 0 or args.min_identity > 100:
        logger.error(f"Invalid min_identity: {args.min_identity}. Must be between 0 and 100")
        sys.exit(1)
    
    if args.min_length < 1:
        logger.error(f"Invalid min_length: {args.min_length}. Must be >= 1")
        sys.exit(1)
    
    if args.evalue <= 0:
        logger.error(f"Invalid evalue: {args.evalue}. Must be > 0")
        sys.exit(1)
    
    if args.threads < 1:
        logger.error(f"Invalid threads: {args.threads}. Must be >= 1")
        sys.exit(1)
    
    if args.mmseqs_sens < 1 or args.mmseqs_sens > 9:
        logger.warning(f"MMseqs sensitivity {args.mmseqs_sens} outside typical range (1-9)")
    
    logger.info(f"Starting iterative search with {args.threads} threads")
    logger.info(f"Parameters: identity>={args.min_identity}%, length>={args.min_length}, evalue<={args.evalue}")
    
    prefix = f"{args.prefix}_" if args.prefix else ""
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/hits", exist_ok=True)
    os.makedirs(f"{args.output_dir}/regions", exist_ok=True)
    
    # CHECKPOINTING: Check for resume
    checkpoint_file = f"{args.output_dir}/.checkpoint"
    start_wave = 0
    current_db = f"{args.output_dir}/current_db.faa"
    
    if args.resume and os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as cf:
            checkpoint_data = json.loads(cf.read())
            start_wave = checkpoint_data.get('completed_waves', 0)
            last_db = checkpoint_data.get('last_db', None)
            
            if last_db and os.path.exists(last_db):
                logger.info(f"Resuming from wave {start_wave + 1}, using DB: {last_db}")
                current_db = last_db
            else:
                logger.warning("Checkpoint found but DB missing, starting from beginning")
                start_wave = 0
                shutil.copyfile(args.initial_db, current_db)
    else:
        shutil.copyfile(args.initial_db, current_db)
    
    # Parse Genomes and Distances
    genome_entries = []
    with open(args.sorted_genomes, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            gname = parts[0]
            dist = float(parts[1]) if len(parts) > 1 else 0.0
            
            gpath = gname
            if args.genomes_dir:
                if not os.path.isabs(gname): # If gname is not an absolute path, assume it's relative to genomes_dir
                    gpath = os.path.join(args.genomes_dir, os.path.basename(gname))
            
            genome_entries.append({'name': gname, 'path': gpath, 'dist': dist})
    
    # Validate we loaded genomes
    if not genome_entries:
        logger.error("No genomes found in sorted_genomes file")
        sys.exit(1)
            
    logger.info(f"Loaded {len(genome_entries)} genomes.")
    
    # Define Waves
    waves = []

    # IMPROVED WAVEFRONT STRATEGY:
    # - Closest genomes (dist < 0.05): Process strictly serially for maximum sensitivity
    # - Medium distance (0.05 - 0.15): Small waves (2-3 genomes)
    # - Distant genomes (> 0.15): Larger waves (can parallelize more)
    
    i = 0
    while i < len(genome_entries):
        curr = genome_entries[i]
        
        if curr['dist'] < 0.05:
            # Very close: Serial processing (wave of 1)
            waves.append([curr])
            i += 1
        elif curr['dist'] < 0.15:
            # Medium distance: Small waves of 2-3 genomes with similar distance
            wave = [curr]
            i += 1
            while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.01:
                wave.append(genome_entries[i])
                i += 1
                if len(wave) >= 3:  # Max 3 per wave for medium distance
                    break
            waves.append(wave)
        else:
            # Distant: Can parallelize more (waves of up to 5)
            wave = [curr]
            i += 1
            while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.02:
                wave.append(genome_entries[i])
                i += 1
                if len(wave) >= 5:  # Max 5 per wave for distant genomes
                    break
            waves.append(wave)
    
    logger.info(f"Defined {len(waves)} waves of execution.")
    
    latest_db = current_db
    
    for i, wave in enumerate(waves):
        # Skip already completed waves
        if i < start_wave:
            logger.info(f"Skipping wave {i+1}/{len(waves)} (already completed)")
            continue
            
        logger.info(f"=== Starting Wave {i+1}/{len(waves)} ({len(wave)} genomes, dist={wave[0]['dist']:.3f}) ===")
        
        # Parallel Execution
        max_workers = min(len(wave), args.threads)
        threads_per_job = max(1, args.threads // max_workers)
        
        logger.info(f"  Running {len(wave)} jobs in parallel with {max_workers} workers, each using {threads_per_job} threads.")
        
        wave_results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for entry in wave:
                futures.append(
                    executor.submit(process_single_genome, 
                                    entry['path'], latest_db, args, args.home_db_dir, prefix, threads_per_job)
                )
            
            for future in concurrent.futures.as_completed(futures):
                try:
                    gname, new_genes = future.result()
                    if new_genes:
                        wave_results.extend(new_genes)
                except Exception as exc:
                    logger.error(f"Wave execution generated an exception for one genome: {exc}")

        # Update DB after Wave
        if wave_results:
            logger.info(f"Wave {i+1} completed. Found {len(wave_results)} new genes. Updating DB.")
            
            new_genes_fasta = f"{args.output_dir}/iter_{i+1}_new_genes.faa"
            # wave_results is list of {'id': ..., 'seq': ...}
            write_fasta([(g['id'], g['seq']) for g in wave_results], new_genes_fasta)
            
            next_db = f"{args.output_dir}/db_iter_{i+1}.faa"
            with open(next_db, 'w') as ndb:
                with open(latest_db, 'r') as old_db:
                    shutil.copyfileobj(old_db, ndb)
                with open(new_genes_fasta, 'r') as new_g:
                    shutil.copyfileobj(new_g, ndb)
            
            # Clean up previous DB if it's not the initial one
            if i > 0 and latest_db != current_db:
                try:
                    os.remove(latest_db)
                except OSError as e:
                    logger.warning(f"Could not remove old DB file {latest_db}: {e}")
            
            latest_db = next_db
        else:
            logger.info(f"Wave {i+1} completed. No new genes found.")
        
        # CHECKPOINT: Save progress after each wave
        with open(checkpoint_file, 'w') as cf:
            json.dump({
                'completed_waves': i + 1,
                'last_db': latest_db,
                'total_waves': len(waves)
            }, cf)
            
    expanded_db = f"{args.output_dir}/expanded_db.faa"
    if os.path.exists(latest_db):
        shutil.move(latest_db, expanded_db)
        
    logger.info(f"Iterative wavefront search complete. Final DB: {expanded_db}")

if __name__ == "__main__":
    main()
