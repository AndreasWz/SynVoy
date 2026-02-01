#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def run_command(cmd):
    # print(f"CMD: {' '.join(cmd)}")
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

def identify_best_synteny_block(hits, max_intron=20000, cluster_dist=50000):
    if not hits:
        return None
        
    # --- Step 1: Group hits by Query ---
    hits_by_query = defaultdict(list)
    for h in hits:
        hits_by_query[h['query']].append(h)
        
    # --- Step 2: Define Loci per Query ---
    all_loci = []
    for query_id, q_hits in hits_by_query.items():
        q_hits.sort(key=lambda x: (x['chrom'], x['start']))
        current_locus_hits = []
        for h in q_hits:
            if not current_locus_hits:
                current_locus_hits.append(h)
                continue
            last_hit = current_locus_hits[-1]
            if (h['chrom'] == last_hit['chrom'] and 
                h['start'] - last_hit['end'] < max_intron):
                current_locus_hits.append(h)
            else:
                all_loci.append(create_locus_object(query_id, current_locus_hits))
                current_locus_hits = [h]
        if current_locus_hits:
            all_loci.append(create_locus_object(query_id, current_locus_hits))

    # --- Step 3: Cluster Loci into Synteny Blocks ---
    all_loci.sort(key=lambda x: (x['chrom'], x['start']))
    if not all_loci: return None
        
    synteny_blocks = []
    current_block = [all_loci[0]]
    
    for locus in all_loci[1:]:
        last_locus = current_block[-1]
        if (locus['chrom'] == last_locus['chrom'] and 
            locus['start'] - last_locus['end'] < cluster_dist):
            current_block.append(locus)
        else:
            synteny_blocks.append(current_block)
            current_block = [locus]
    synteny_blocks.append(current_block)
    
    # --- Step 4: Score Blocks ---
    best_block = None
    best_gene_count = -1
    
    for block in synteny_blocks:
        # Score = Count of Unique Query Genes (stripping variants)
        # Assuming original queries don't have pipes, or we want the base ID
        unique_genes = set(l['query'].split('|')[0] for l in block)
        gene_count = len(unique_genes)
        
        if gene_count > best_gene_count:
            best_gene_count = gene_count
            best_block = block
        elif gene_count == best_gene_count:
            if len(block) > len(best_block):
                best_block = block

    if not best_block: return None

    chrom = best_block[0]['chrom']
    start = min(l['start'] for l in best_block)
    end = max(l['end'] for l in best_block)
    
    genes_list = list(set(l['query'].split('|')[0] for l in best_block))
    
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'genes_count': best_gene_count,
        'loci_count': len(best_block),
        'genes_list': genes_list
    }

def annotate_and_filter_predicted_genes(fasta_file, hits, region_offset, genome_name):
    """
    Reads Prodigal predicted proteins.
    Renames them if they overlap with a Search Hit (inheriting the Query ID).
    Returns list of SeqRecords to save.
    """
    annotated_records = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Parse Prodigal Header: >id # start # end # ...
        desc_parts = record.description.split(' # ')
        if len(desc_parts) >= 3:
            local_start = int(desc_parts[1])
            local_end = int(desc_parts[2])
            
            # Prodigal Output is relative to the extracted sequence
            global_start = region_offset + local_start
            global_end = region_offset + local_end
            
            # Find overlapping hits
            best_hit_query = None
            max_overlap = 0
            
            for h in hits:
                # Basic interval intersection
                overlap_start = max(global_start, h['start'])
                overlap_end = min(global_end, h['end'])
                overlap_len = max(0, overlap_end - overlap_start)
                
                if overlap_len > 0:
                    if overlap_len > max_overlap:
                        max_overlap = overlap_len
                        best_hit_query = h['query'].split('|')[0] # Parent ID
            
            # New ID
            # Clean genome name
            clean_gname = genome_name.replace(".fna", "").replace(".fasta", "")
            
            if best_hit_query:
                # Inherit Parent ID
                new_id = f"{best_hit_query}|{clean_gname}_{record.id}"
                record.id = new_id
                record.description = f"coords:{global_start}-{global_end} parent:{best_hit_query}"
                annotated_records.append(record)
            else:
                # Orphan gene - do we want to add non-syntenic neighbors?
                # The user said: "locate the gene". Usually we only want the orthologs.
                # Adding neighbors is risky for iteration as it creates drift.
                # Let's SKIP orphans for the DB update, to keep the profile clean.
                pass
                
    return annotated_records


# ... (Previous parts) ...

def estimate_cluster_dist(genome_file, default_dist=50000):
    """
    Estimate gene density to adjust cluster_dist.
    Approximation: If genome is small/dense (Bacteria), use 20kb.
    If large/sparse (Euk), use 100kb+?
    Simple heuristic: Size / 1000 genes?
    We don't know gene count.
    Just use file size as proxy for now?
    Bacteria ~5Mb. Mammal ~3Gb.
    If size < 10Mb -> assume Bacteria -> 20kb.
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

def batch_rbh_check(candidates, home_db, unique_id_map, threads=1, evalue=1e-5):
    """
    Perform Reciprocal Best Hit check for multiple candidates at once.
    candidates: list of SeqRecords
    home_db: path to MMseqs DB
    unique_id_map: dict mapping candidate.id -> parent_query_id
    """
    if not candidates: return []
    
    # Write all to one file
    query_fasta = "batch_candidates.fasta"
    with open(query_fasta, 'w') as f:
        SeqIO.write(candidates, f, "fasta")
        
    # Run MMseqs
    rbh_out = "batch_rbh.m8"
    
    db_path = home_db
    if os.path.isdir(home_db):
        db_path = os.path.join(home_db, "db")
        
    cmd = [
        "mmseqs", "easy-search",
        query_fasta, db_path, rbh_out, "tmp_rbh_batch",
        "-e", str(evalue),
        "--format-output", "query,target,pident,evalue,bits",
        "--max-seqs", "1", # Top hit only
        "--threads", str(threads)
    ]
    
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
    
    valid_ids = set()
    if os.path.exists(rbh_out):
        with open(rbh_out) as f:
            for line in f:
                parts = line.strip().split('\t')
                cand_id = parts[0]
                target_id = parts[1]
                
                # Check match
                if cand_id in unique_id_map:
                    parent = unique_id_map[cand_id]
                    # Target should match parent
                    # Parent might be "GeneA". Target might be "GeneA|..." or just "GeneA"
                    if parent in target_id or target_id in parent:
                         valid_ids.add(cand_id)
                         
    return valid_ids

def main():
    parser = argparse.ArgumentParser(description="Iterative Genome Search Runner")
    parser.add_argument("--initial_db", required=True)
    parser.add_argument("--sorted_genomes", required=True)
    parser.add_argument("--genomes_dir", help="Directory containing genome files")
    parser.add_argument("--home_db_dir", help="Home Proteome MMseqs DB for RBH")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--diamond_sensitivity", default="very-sensitive")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--expand_threshold", type=float, default=1e-10)
    parser.add_argument("--min_identity", type=float, default=40.0)
    parser.add_argument("--min_length", type=int, default=50)
    parser.add_argument("--threads", default="4")
    parser.add_argument("--cluster_dist", type=int, default=-1, help="Auto-detect if -1")
    parser.add_argument("--prefix", default="", help="Prefix for output files (e.g. locus ID)")
    
    args = parser.parse_args()
    
    prefix = f"{args.prefix}_" if args.prefix else ""
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/hits", exist_ok=True)
    os.makedirs(f"{args.output_dir}/regions", exist_ok=True)
    
    current_db = f"{args.output_dir}/current_db.faa"
    shutil.copyfile(args.initial_db, current_db)
    
    genomes = []
    with open(args.sorted_genomes, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            gname = parts[0]
            if args.genomes_dir:
                gname = os.path.basename(gname)
                gpath = os.path.join(args.genomes_dir, gname)
                genomes.append(gpath)
            else:
                genomes.append(gname)

    print(f"Loaded {len(genomes)} genomes for iteration.")
    
    latest_db = current_db
    
    for i, genome_file in enumerate(genomes):
        genome_name = os.path.basename(genome_file)
        print(f"=== Processing Genome {i+1}/{len(genomes)}: {genome_name} ===")
        
        if not os.path.exists(genome_file):
            print(f"Warning: Genome file {genome_file} not found. Skipping.")
            continue
            
        hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
        
        # Auto-param
        c_dist = args.cluster_dist
        if c_dist <= 0:
            c_dist = estimate_cluster_dist(genome_file)
            
        # Run MMseqs (Search)
        subprocess.run([
            "mmseqs", "easy-search",
            latest_db, genome_file, hits_file, f"{args.output_dir}/tmp_mmseqs",
            "--search-type", "2", 
            "--threads", str(args.threads),
            "-s", "7.5",
            "-e", str(args.evalue),
            "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        ], check=True, stderr=subprocess.DEVNULL)

        hits = parse_hits(hits_file, args.min_identity, args.min_length, args.evalue)
        
        # Identify Synteny
        best_region = identify_best_synteny_block(hits, cluster_dist=c_dist)
        
        if not best_region:
            print(f"No valid syntenic region found in {genome_name}.")
            continue
            
        print(f"Found Best Syntenic Region: {best_region['chrom']}:{best_region['start']}-{best_region['end']} "
              f"(Matched Genes: {best_region['genes_count']})")
        
        # Extract Region & Predict
        genome_seqs = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
        chrom = best_region['chrom']
        if chrom not in genome_seqs: continue
        
        slen = len(genome_seqs[chrom])
        w_start = max(0, best_region['start'] - 20000)
        w_end = min(slen, best_region['end'] + 20000)
        
        subseq = genome_seqs[chrom].seq[w_start:w_end]
        
        temp_fa = f"{args.output_dir}/tmp_reg.fasta"
        temp_faa = f"{args.output_dir}/tmp_reg.faa"
        
        with open(temp_fa, 'w') as tf:
            SeqIO.write(SeqRecord(subseq, id="region_seq", description=""), tf, "fasta")
        
        try:
            subprocess.run(["prodigal", "-i", temp_fa, "-a", temp_faa, "-p", "meta", "-q"], check=True)
            
            # Annotate Candidates
            relevant_hits = [h for h in hits if h['chrom'] == chrom]
            annotated_records_raw = annotate_and_filter_predicted_genes(temp_faa, relevant_hits, w_start, genome_name)
            
            # --- BATCH ORTHOLOGY VALIDATION (RBH) ---
            annotated_genes = []
            if args.home_db_dir and annotated_records_raw:
                # Prepare map
                cand_map = {}
                for rec in annotated_records_raw:
                    parent_id = rec.id.split('|')[0]
                    cand_map[rec.id] = parent_id
                
                # Run Batch
                valid_ids = batch_rbh_check(annotated_records_raw, args.home_db_dir, cand_map, threads=args.threads)
                
                annotated_genes = [rec for rec in annotated_records_raw if rec.id in valid_ids]
                print(f"RBH Filter: kept {len(annotated_genes)}/{len(annotated_records_raw)} candidates.")
            else:
                annotated_genes = annotated_records_raw
            
            # Update DB
            if annotated_genes:
                new_genes_fasta = f"{args.output_dir}/regions/{prefix}{genome_name}_new_genes.faa"
                with open(new_genes_fasta, 'w') as out_faa:
                    SeqIO.write(annotated_genes, out_faa, "fasta")
                
                print(f"Adding {len(annotated_genes)} sequences to stats...")
                
                # Append to DB
                next_db = f"{args.output_dir}/db_iter_{i+1}.faa"
                with open(next_db, 'w') as ndb:
                    with open(latest_db, 'r') as old_db:
                        shutil.copyfileobj(old_db, ndb)
                    with open(new_genes_fasta, 'r') as new_g:
                        shutil.copyfileobj(new_g, ndb)
                
                if i > 0 and latest_db != current_db:
                    try:
                        os.remove(latest_db)
                    except: pass
                
                latest_db = next_db
            else:
                print("No validated genes found. Skipping DB update.")
                
        except Exception as e:
            print(f"Processing failed for {genome_name}: {e}")

    expanded_db = f"{args.output_dir}/expanded_db.faa"
    if os.path.exists(latest_db):
        shutil.move(latest_db, expanded_db)
        
    print(f"Iterative search complete. Final DB: {expanded_db}")

if __name__ == "__main__":
    main()

