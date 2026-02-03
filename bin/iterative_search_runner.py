#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
import concurrent.futures
import math
import uuid
import sys
from collections import defaultdict

# Use our own sequence utilities (no BioPython dependency)
try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )
except ImportError:
    # Fallback if not in path - add bin directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )

# Import fragment utilities if available
try:
    from fragment_query import generate_fragments, parse_fragment_id, merge_fragment_hits
    FRAGMENT_SUPPORT = True
except ImportError:
    FRAGMENT_SUPPORT = False

def run_command(cmd):
    subprocess.check_call(cmd)

def normalize_coordinates(start, end):
    return min(start, end), max(start, end)

def parse_hits(hits_file, min_identity, min_length, evalue_thresh):
    """
    Parse MMseqs2 hits and return a list of hit dictionaries.
    Filters by basic quality metrics.
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
                    
                    hits.append({
                        'query': parts[0],
                        'target': parts[1], # Chromosome/Scaffold
                        'chrom': parts[1],
                        'start': start,
                        'end': end,
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

def extract_base_gene_id(query_id):
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

def calculate_adaptive_padding(hits, best_region, default=100000):
    """
    Calculate region padding based on gene spacing in hits.
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

def run_miniprot(target_fasta, query_protein, output_paf):
    """
    Run miniprot to align protein query to target DNA.
    Returns list of hit objects (parsed from PAF/GFF).
    """
    cmd = [
        "miniprot", "-I", "--gff", 
        target_fasta, query_protein
    ]
    
    hits = []
    
    try:
        # Capture stdout
        with open(output_paf, "w") as outfile:
            subprocess.run(cmd, stdout=outfile, check=True, stderr=subprocess.DEVNULL)
            
        # Parse GFF
        current_hit = None
        with open(output_paf, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9: continue
                
                feat_type = parts[2]
                
                # mRNA line starts a new hit
                if feat_type == "mRNA":
                    try:
                        info = {}
                        for item in parts[8].split(";"):
                            if "=" in item:
                                k, v = item.split("=", 1)
                                info[k] = v
                        
                        target_name = "Unknown"
                        if "Target" in info:
                            target_name = info["Target"].split()[0]
                        
                        hit = {
                            'id': info.get('ID', 'unknown'),
                            'parent_query': target_name,
                            'chrom': parts[0],
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'strand': parts[6],
                            'identity': float(info.get('Identity', 0)) * 100,
                            'score': float(parts[5]) if parts[5] != '.' else 0,
                            'cds_parts': [],
                            'gff_lines': [line.strip()] # Store raw lines
                        }
                        hits.append(hit)
                        current_hit = hit
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Failed to parse mRNA line: {e}", file=sys.stderr)
                        current_hit = None
                        continue
                    
                elif current_hit:
                     # Associate following lines (CDS, etc.) with current hit
                     current_hit['gff_lines'].append(line.strip())
                     if feat_type == "CDS":
                         try:
                             current_hit['cds_parts'].append((int(parts[3]), int(parts[4])))
                         except (ValueError, IndexError) as e:
                             print(f"Warning: Failed to parse CDS coordinates: {e}", file=sys.stderr)
                     
    except Exception as e:
        print(f"Miniprot failed: {e}")
        raise subprocess.CalledProcessError(1, cmd)
        
    return hits

def extract_cds_sequence(genome_seq, hit):
    """
    Extracts and concatenates CDS sequences based on Miniprot alignment.
    Handles strand orientation.
    
    Returns dict with 'id' and 'seq' keys.
    """
    seq_str = ""
    # Sort CDS by position (GFF usually is sorted but be safe)
    sorted_cds = sorted(hit['cds_parts'], key=lambda x: x[0])
    
    for start, end in sorted_cds:
        # GFF is 1-based, inclusive. Python slice is 0-based.
        # start-1 to end
        exon = genome_seq[start-1:end]
        seq_str += str(exon)
        
    # Reverse complement if negative strand
    if hit['strand'] == '-':
        seq_str = reverse_complement(seq_str)
    
    # Translate to protein
    protein_seq = translate(seq_str)
    # Remove stop codon if present at end
    if protein_seq.endswith('*'):
        protein_seq = protein_seq[:-1]
        
    return {'id': hit['id'], 'seq': protein_seq}


def estimate_cluster_dist(genome_file, default_dist=50000):
    """
    Estimate gene density to adjust cluster_dist.
    Approximation: If genome is small/dense (Bacteria), use 20kb.
    If large/sparse (Euk), use 100kb+?
    Simple heuristic: If size < 10Mb -> assume Bacteria -> 20kb.
    If size > 100Mb -> assume Euk -> 100kb.
    """
    try:
        size = os.path.getsize(genome_file)
        if size < 10_000_000: # < 10MB
            return 20000
        elif size > 100_000_000: # > 100MB
            return 100000
    except: pass
    return default_dist

def batch_rbh_check(candidates, home_db, cand_map, threads=1, evalue=1e-5):
    """
    Perform Reciprocal Best Hit check for multiple candidates at once.
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
            "--format-output", "query,target,pident,evalue,bits",
            "--max-seqs", "1", # Top hit only
            "--threads", str(threads)
        ]
        
        subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
        
        if os.path.exists(rbh_out):
            with open(rbh_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 2: continue
                    cand_id = parts[0]
                    target_id = parts[1]
                    
                    # Debug print
                    # print(f"DEBUG RBH: cand={cand_id} target={target_id}")
                    
                    # Check match - Use exact matching to avoid false positives
                    if cand_id in cand_map:
                        parent = cand_map[cand_id]
                        # Extract base IDs without suffixes/variants
                        parent_base = extract_base_gene_id(parent).strip()
                        target_base = extract_base_gene_id(target_id).strip()
                        
                        # Exact match or one is contained as full word
                        if parent_base == target_base or parent == target_id or target_id == parent:
                             valid_ids.add(cand_id)
    except Exception as e:
        print(f"RBH check failed: {e}")
        return set() # Fail safe? Or return empty
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
        print(f"[{genome_name}] Warning: Genome file not found. Skipping.")
        return genome_name, []
    
    # Create unique temp space
    unique_id = uuid.uuid4().hex
    hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
    tmp_dir = f"{args.output_dir}/tmp_mmseqs_{unique_id}_{genome_name}"
    
    new_genes = []
    
    # Database Index Cache
    db_index = None
    
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
            "-s", "7.5",
            "-e", str(args.evalue),
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        ], check=True, stderr=subprocess.DEVNULL)

        hits = parse_hits(hits_file, args.min_identity, args.min_length, args.evalue)
        if not hits:
            print(f"[{genome_name}] No hits found in MMseqs output.", flush=True)
            return genome_name, []
            
        print(f"[{genome_name}] Parsed {len(hits)} hits.", flush=True)

        # 2. Identify Synteny
        best_region = identify_best_synteny_block(hits, cluster_dist=c_dist)
        
        if not best_region:
            print(f"[{genome_name}] No valid syntenic region found.", flush=True)
            return genome_name, []
            
        print(f"[{genome_name}] Found syntenic region: {best_region['chrom']}:{best_region['start']}-{best_region['end']} with {len(best_region['genes'])} genes.", flush=True)
            # print(f"[{genome_name}] Found Region: {best_region['chrom']}:{best_region['start']}-{best_region['end']}")
            
        # 3. Extract Region - using our FASTA parser
        genome_seqs = load_genome(genome_path)
        chrom = best_region['chrom']
        
        if chrom in genome_seqs:
            slen = len(genome_seqs[chrom])
            
            # ADAPTIVE PADDING
            padding = calculate_adaptive_padding(hits, best_region, default=100000)
            
            w_start = max(0, best_region['start'] - padding)
            w_end = min(slen, best_region['end'] + padding)
            subseq = genome_seqs[chrom][w_start:w_end]
            
            temp_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_reg.fasta"
            miniprot_paf = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_reg.paf"
            
            # Write region FASTA
            write_fasta([("region_seq", subseq)], temp_fa)
            
            # 4. Prepare Hits & Query Sequences
            relevant_hits = [h for h in hits if h['chrom'] == chrom]
            unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)
            
            # Extract sequences from DB - using our parser
            found_queries = []  # List of {'id': ..., 'seq': ...}
            db_sequences = {}
            try:
                # Parse DB once into dict
                for header, clean_id, seq in parse_fasta(db_path):
                    db_sequences[clean_id] = {'id': clean_id, 'seq': seq, 'header': header}
                    # Also index by base ID
                    base = extract_base_id(clean_id)
                    if base != clean_id and base not in db_sequences:
                        db_sequences[base] = {'id': clean_id, 'seq': seq, 'header': header}
                
                for query_id in unique_queries:
                    if query_id in db_sequences:
                        found_queries.append(db_sequences[query_id])
                    else:
                        base = extract_base_id(query_id)
                        if base in db_sequences:
                            found_queries.append(db_sequences[base])
            except Exception as dex:
                print(f"[{genome_name}] Warning: DB parsing failed: {dex}")

            if found_queries:
                print(f"[{genome_name}] Running Miniprot with {len(found_queries)} queries...", flush=True)
                query_mini_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_query.faa"
                # Write query FASTA
                write_fasta([(q['id'], q['seq']) for q in found_queries], query_mini_fa)
                
                try:
                    # 5. Run Miniprot
                    print(f"[{genome_name}] CMD: miniprot -I --gff {temp_fa} {query_mini_fa}", flush=True)
                    miniprot_hits = run_miniprot(temp_fa, query_mini_fa, miniprot_paf)
                    print(f"[{genome_name}] Miniprot found {len(miniprot_hits)} raw hits.", flush=True)
                    
                    # 6. Process Miniprot Hits
                    # Group hits by parent_query to handle multi-exon genes
                    hits_by_parent = defaultdict(list)
                    for hit in miniprot_hits:
                        hits_by_parent[hit['parent_query']].append(hit)
                    
                    annotated_records_raw = []  # List of {'id': ..., 'seq': ...}
                    valid_gff_lines = []
                    
                    # Clean genome name for use in IDs (remove special chars)
                    clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')
                    
                    # Process each gene (may have multiple exons/hits)
                    for parent_id, gene_hits in hits_by_parent.items():
                        try:
                            # Sort hits by position
                            gene_hits.sort(key=lambda h: h['start'])
                            
                            # Check if hits are close enough to be same gene
                            # If they span >max_intron, they might be paralogs
                            first_hit = gene_hits[0]
                            last_hit = gene_hits[-1]
                            span = last_hit['end'] - first_hit['start']
                            
                            # If span is huge (>500kb), treat as separate genes
                            MAX_GENE_SPAN = 500000
                            if span > MAX_GENE_SPAN and len(gene_hits) > 1:
                                print(f"[{genome_name}] Warning: Hits for {parent_id} span {span}bp, "
                                      f"treating as {len(gene_hits)} separate loci.", flush=True)
                                # Process each hit individually as potential paralogs
                                for idx, hit in enumerate(gene_hits):
                                    try:
                                        cds_seq_record = extract_cds_sequence(subseq, hit)
                                        global_start = w_start + hit['start']
                                        global_end = w_start + hit['end']
                                        
                                        new_id = f"{parent_id}|{clean_gname}_{hit['id']}_paralog{idx+1}"
                                        cds_seq_record.id = new_id
                                        cds_seq_record.description = f"coords:{global_start}-{global_end} parent:{parent_id} score:{hit['score']:.1f} identity:{hit['identity']:.1f}"
                                        annotated_records_raw.append(cds_seq_record)
                                        
                                        # Shift GFF
                                        for gline in hit['gff_lines']:
                                            gp = gline.split('\t')
                                            if len(gp) < 9: continue
                                            gp[0] = chrom
                                            gp[3] = str(int(gp[3]) + w_start)
                                            gp[4] = str(int(gp[4]) + w_start)
                                            if gp[2] == "mRNA":
                                                gp[8] += f";SynTerra_Parent={parent_id};SynTerra_ID={new_id}"
                                            valid_gff_lines.append('\t'.join(gp))
                                    except Exception as ex:
                                        print(f"[{genome_name}] Failed to process paralog {idx+1}: {ex}", flush=True)
                            else:
                                # Consolidate into single gene annotation
                                # Use first hit's ID but combine all CDS
                                all_cds = []
                                all_gff = []
                                for hit in gene_hits:
                                    all_cds.extend(hit['cds_parts'])
                                    all_gff.extend(hit['gff_lines'])
                                
                                # Create consolidated hit
                                consolidated_hit = {
                                    'id': f"{first_hit['id']}_consolidated",
                                    'parent_query': parent_id,
                                    'start': first_hit['start'],
                                    'end': last_hit['end'],
                                    'strand': first_hit['strand'],
                                    'score': sum(h['score'] for h in gene_hits),
                                    'identity': sum(h['identity'] for h in gene_hits) / len(gene_hits),
                                    'cds_parts': all_cds,
                                    'gff_lines': all_gff
                                }
                                
                                # Extract CDS sequence
                                cds_result = extract_cds_sequence(subseq, consolidated_hit)
                                global_start = w_start + consolidated_hit['start']
                                global_end = w_start + consolidated_hit['end']
                                
                                # Global unique ID for the sequence
                                new_id = f"{parent_id}|{clean_gname}_{consolidated_hit['id']}"
                                description = f"coords:{global_start}-{global_end} parent:{parent_id} score:{consolidated_hit['score']:.1f} identity:{consolidated_hit['identity']:.1f}"
                                
                                # Store record as dict (no BioPython)
                                cds_seq_record = {
                                    'id': new_id,
                                    'seq': cds_result['seq'],
                                    'description': description
                                }
                                annotated_records_raw.append(cds_seq_record)
                                
                                # Store and Shift GFF lines
                                shifted_lines = []
                                for gline in consolidated_hit['gff_lines']:
                                    gp = gline.split('\t')
                                    if len(gp) < 9: continue
                                    # Update seqid to real chrom
                                    gp[0] = chrom
                                    # Update coords
                                    gp[3] = str(int(gp[3]) + w_start)
                                    gp[4] = str(int(gp[4]) + w_start)
                                    # Update attributes
                                    if gp[2] == "mRNA":
                                        gp[8] += f";SynTerra_Parent={parent_id};SynTerra_ID={new_id}"
                                    
                                    shifted_lines.append('\t'.join(gp))
                                
                                valid_gff_lines.extend(shifted_lines)
                                
                        except Exception as ex:
                            print(f"[{genome_name}] Failed to process gene {parent_id}: {ex}", flush=True)
                    
                    # 7. RBH Validation
                    if home_db_dir and annotated_records_raw:
                        print(f"[{genome_name}] Running RBH for {len(annotated_records_raw)} candidates...", flush=True)
                        cand_map = {rec['id']: extract_base_gene_id(rec['id']) for rec in annotated_records_raw}
                        valid_ids = batch_rbh_check(annotated_records_raw, home_db_dir, cand_map, threads=threads_per_job)
                        new_genes = [rec for rec in annotated_records_raw if rec['id'] in valid_ids]
                        
                        # Filter GFF lines to only RBH-validated genes
                        valid_parent_ids = set(extract_base_gene_id(rec['id']) for rec in new_genes)
                        final_gff_lines = []
                        for gline in valid_gff_lines:
                            # Check if this GFF line belongs to a validated gene
                            # Look for SynTerra_Parent in attributes
                            if 'SynTerra_Parent=' in gline:
                                for parent_id in valid_parent_ids:
                                    if f'SynTerra_Parent={parent_id}' in gline:
                                        final_gff_lines.append(gline)
                                        break
                            else:
                                # Fallback: include all if no parent annotation
                                final_gff_lines.append(gline)
                        
                        valid_gff_lines = final_gff_lines
                        print(f"[{genome_name}] RBH Kept {len(new_genes)} genes.", flush=True)
                    else:
                        print(f"[{genome_name}] Skipping RBH (records={len(annotated_records_raw)})", flush=True)
                        new_genes = annotated_records_raw
                    
                    # Write GFF and Homology TSV for the genome
                    if valid_gff_lines:
                        gff_out = f"{args.output_dir}/regions/{genome_name}.gff"
                        faa_out = f"{args.output_dir}/regions/{genome_name}.faa"
                        tsv_out = f"{args.output_dir}/regions/{genome_name}.homology.tsv"
                        
                        with open(gff_out, 'w') as gf:
                            gf.write("##gff-version 3\n")
                            for gl in valid_gff_lines:
                                gf.write(gl + "\n")

                        # Write FASTA - new_genes is list of {'id': ..., 'seq': ...}
                        write_fasta([(g['id'], g['seq']) for g in new_genes], faa_out)
                            
                        with open(tsv_out, 'w') as tf:
                            for rec in new_genes:
                                # rec['id'] is the full unique ID, parent is the Home ID
                                parent = extract_base_gene_id(rec['id'])
                                tf.write(f"{rec['id']}\t{parent}\n")

                except subprocess.CalledProcessError as miniprot_err:
                    print(f"[{genome_name}] Error: Miniprot failed. Skipping genome. {miniprot_err}")
                    new_genes = [] # Fail safe
                    
                # Cleanup temp
                if os.path.exists(temp_fa): os.remove(temp_fa)
                if os.path.exists(miniprot_paf): os.remove(miniprot_paf)
                if os.path.exists(query_mini_fa): os.remove(query_mini_fa)
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
        if db_index:
            db_index.close()
            
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
    parser.add_argument("--prefix", default="", help="Prefix for output files (e.g. locus ID)")
    
    args = parser.parse_args()
    prefix = f"{args.prefix}_" if args.prefix else ""
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/hits", exist_ok=True)
    os.makedirs(f"{args.output_dir}/regions", exist_ok=True)
    
    current_db = f"{args.output_dir}/current_db.faa"
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
            
    print(f"Loaded {len(genome_entries)} genomes.")
    
    # Define Waves
    waves = []
    if not genome_entries:
        print("No genomes to process.")
        return

    current_wave = [genome_entries[0]]
    for i in range(1, len(genome_entries)):
        prev = current_wave[-1]
        curr = genome_entries[i]
        
        # Wavefront heuristic: group genomes with the same or very similar distance
        if abs(curr['dist'] - prev['dist']) < 0.001: # Using a small epsilon for float comparison
            current_wave.append(curr)
        else:
            waves.append(current_wave)
            current_wave = [curr]
    waves.append(current_wave) # Add the last wave
    
    print(f"Defined {len(waves)} waves of execution.")
    
    latest_db = current_db
    
    for i, wave in enumerate(waves):
        print(f"=== Starting Wave {i+1}/{len(waves)} ({len(wave)} genomes, dist={wave[0]['dist']:.3f}) ===")
        
        # Parallel Execution
        max_workers = min(len(wave), args.threads)
        threads_per_job = max(1, args.threads // max_workers)
        
        print(f"  Running {len(wave)} jobs in parallel with {max_workers} workers, each using {threads_per_job} threads.")
        
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
                    print(f"Wave execution generated an exception for one genome: {exc}")

        # Update DB after Wave
        if wave_results:
            print(f"Wave {i+1} completed. Found {len(wave_results)} new genes. Updating DB.")
            
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
                    print(f"Warning: Could not remove old DB file {latest_db}: {e}")
            
            latest_db = next_db
        else:
            print(f"Wave {i+1} completed. No new genes found.")
            
    expanded_db = f"{args.output_dir}/expanded_db.faa"
    if os.path.exists(latest_db):
        shutil.move(latest_db, expanded_db)
        
    print(f"Iterative wavefront search complete. Final DB: {expanded_db}")

if __name__ == "__main__":
    main()
