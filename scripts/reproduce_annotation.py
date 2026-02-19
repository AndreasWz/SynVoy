#!/usr/bin/env python3
import sys
import os
import shutil
import argparse
import subprocess
import json
import re
from collections import defaultdict

# Add bin to path
sys.path.insert(0, os.path.abspath("bin"))

from iterative_search_runner import (
    parse_hits,
    identify_synteny_blocks,
    process_region_block,
    extract_base_gene_id,
    deduplicate_flanking_models,
    collapse_flanking_cds_to_gene_span
)
from sequence_utils import load_genome, parse_gff, parse_fasta, write_fasta
from flanking_query_utils import collapse_flanking_query_records

def write_bed(hits, output_file, chrom_offset=0, chrom_override=None):
    """Write hits to BED format."""
    with open(output_file, 'w') as f:
        for hit in hits:
            # Input hits are 1-based inclusive. BED needs 0-based half-open.
            start = max(0, hit['start'] + chrom_offset - 1)
            end = hit['end'] + chrom_offset
            chrom = chrom_override or hit['chrom']
            name = hit['query']
            score = hit['bits'] if 'bits' in hit else 0
            strand = hit['strand']
            # BED: chrom, start, end, name, score, strand
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

def first_truth_chrom(truth_gff):
    """Return first non-comment seqid/chrom from ground-truth GFF."""
    if not truth_gff or not os.path.exists(truth_gff):
        return None
    with open(truth_gff) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 1 and parts[0]:
                return parts[0]
    return None

def build_query_fasta(args):
    """
    Resolve query FASTA path:
    1) Use --queries if provided and exists.
    2) Otherwise combine --goi_queries + --flanking_queries into one FASTA.
    GOI records are prefixed with GOI_ so iterative logic can treat them specially.
    """
    if args.queries and os.path.exists(args.queries):
        return args.queries

    if not os.path.exists(args.goi_queries):
        raise FileNotFoundError(f"GOI query file not found: {args.goi_queries}")
    if not os.path.exists(args.flanking_queries):
        raise FileNotFoundError(f"Flanking query file not found: {args.flanking_queries}")

    combined_queries = os.path.join(args.outdir, "queries_combined.fasta")
    records = []
    seen = set()
    seen_goi_seq = set()

    for _, clean_id, seq in parse_fasta(args.goi_queries):
        goi_id = clean_id if clean_id.startswith("GOI_") else f"GOI_{clean_id}"
        if seq in seen_goi_seq:
            continue
        seen_goi_seq.add(seq)
        key = (goi_id, seq)
        if key in seen:
            continue
        seen.add(key)
        records.append((goi_id, seq))

    flanking_raw = list(parse_fasta(args.flanking_queries))
    flanking_records, flanking_stats = collapse_flanking_query_records(flanking_raw)
    print(
        "Normalized flanking queries: {count} genes from {inp} input records "
        "(exon_reconstructed={recon}, fragment_collapsed={frag}, dropped_empty={drop})".format(
            count=len(flanking_records),
            inp=flanking_stats.get("input_records", 0),
            recon=flanking_stats.get("exon_reconstructed", 0),
            frag=flanking_stats.get("fragment_collapsed", 0),
            drop=flanking_stats.get("dropped_empty", 0),
        )
    )

    for clean_id, seq in flanking_records:
        key = (clean_id, seq)
        if key in seen:
            continue
        seen.add(key)
        records.append((clean_id, seq))

    write_fasta(records, combined_queries)
    print(f"Built combined queries FASTA: {combined_queries} ({len(records)} sequences)")
    return combined_queries

def compare_annotations(predicted_hits, truth_gff, offset=0, tolerance=10):
    """
    Compare predicted hits (list of dicts) against ground truth GFF.
    """
    print(f"\n[Comparison] Comparing {len(predicted_hits)} predictions vs {truth_gff}...")
    
    # Load Truth
    truth_genes = parse_gff(truth_gff)
    truth_exons = [f for f in truth_genes if f['type'] == 'CDS']
    
    print(f"Ground Truth: {len(truth_exons)} CDS features")
    
    tp = 0
    fp = 0
    fn = 0
    
    matched_truth = set()
    matched_pred = set()
    
    def _norm_chrom(chrom):
        if chrom is None:
            return ""
        return str(chrom).replace("region_seq", "region").split('.')[0]

    for i, t in enumerate(truth_exons):
        match_found = False
        t_chrom = t.get('chrom', t.get('seqid')) # Handle both
        
        for j, p in enumerate(predicted_hits):
            # Adjust prediction to genomic coords
            p_start = p['start'] + offset
            p_end = p['end'] + offset
            p_chrom = p.get('chrom', p.get('target', 'unknown'))
            
            # Check overlap or boundary match
            # Relaxed chrom check: if p_chrom is 'region' (from extracted FASTA), assume it matches
            chrom_match = (_norm_chrom(t_chrom) == _norm_chrom(p_chrom)) or (p_chrom == 'region') or (p_chrom == 'region_seq')
            
            if chrom_match and t['strand'] == p['strand']:
                # Exact boundary match (with tolerance)
                start_diff = abs(t['start'] - p_start)
                end_diff = abs(t['end'] - p_end)
                
                # Or check substantial overlap
                overlap_start = max(t['start'], p_start)
                overlap_end = min(t['end'], p_end)
                overlap_len = max(0, overlap_end - overlap_start)
                t_len = t['end'] - t['start']
                p_len = p_end - p_start
                
                if (start_diff <= tolerance and end_diff <= tolerance) or \
                   (overlap_len > 0.8 * t_len and overlap_len > 0.8 * p_len):
                    match_found = True
                    matched_pred.add(j)
                    # print(f"  Match: Truth {t['start']}-{t['end']} vs Pred {p_start}-{p_end}")
                    break
        
        if match_found:
            tp += 1
            matched_truth.add(i)
        else:
            fn += 1
            # Report missed item
            t_len = t['end'] - t['start']
            t_name = t.get('attributes', {}).get('name', 'unknown')
            t_note = t.get('attributes', {}).get('note', '')
            t_chrom = t.get('chrom', t.get('seqid'))
            print(f"  [MISS] {t_name} ({t_note}) at {t_chrom}:{t['start']}-{t['end']} ({t['strand']}), len={t_len}bp")
            
    fp = len(predicted_hits) - len(matched_pred)
    
    # Report false positives
    for j, p in enumerate(predicted_hits):
        if j not in matched_pred:
            p_name = p['query']
            print(f"  [EXTRA] {p_name} at {p['chrom']}:{p['start']+offset}-{p['end']+offset} ({p['strand']})")
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"Precision: {precision:.2f}")
    print(f"Recall:    {recall:.2f}")
    print(f"F1 Score:  {f1:.2f}")
    
    return {
        'precision': precision,
        'recall': recall,
        'f1': f1
    }

def main():
    parser = argparse.ArgumentParser(description="Reproduce Annotation Logic")
    parser.add_argument("--genome", default="tests/ground_truth_test/input/genome.fasta")
    parser.add_argument("--queries", default=None,
                        help="Combined query FASTA. If omitted, GOI/flanking inputs are merged.")
    parser.add_argument("--goi_queries", default="tests/ground_truth_test/input/GOI_query.fasta")
    parser.add_argument("--flanking_queries", default="tests/ground_truth_test/input/flanking_genes.fasta")
    parser.add_argument("--ground_truth", default="tests/ground_truth_test/ground_truth/ground_truth.gff")
    parser.add_argument("--outdir", default="tests/ground_truth_test/output")
    parser.add_argument("--offset", type=int, default=0, help="Genomic offset for coordinates")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    query_fasta = build_query_fasta(args)
    truth_chrom = first_truth_chrom(args.ground_truth)
    
    # Mock args object for iterative_search_runner functions
    class MockArgs:
        def __init__(self):
            self.output_dir = args.outdir
            self.mmseqs_sens = 9.5 # Max sensitivity
            self.evalue = 10.0
            self.min_identity = 30
            self.min_length = 5
            self.max_intron = 20000
            self.cluster_dist = 50000
            self.region_padding = 20000
            self.padding_min = 5000
            self.padding_max = 50000
            # Keep harness portable: many environments here lack parasail/ssearch36.
            self.enable_smith_waterman = False
            self.sw_method = "auto"
            self.sw_min_score = 15
            self.sw_min_identity = 20
            self.sw_timeout_seconds = 300
            self.aug_relaxed_evalue_mult = 2000000000.0
            self.aug_relaxed_evalue_cap = 20000.0
            self.aug_relaxed_parse_evalue_mult = 10
            self.aug_relaxed_identity_factor = 0.6
            self.aug_relaxed_identity_min = 25.0
            self.aug_relaxed_length_div = 2
            self.aug_relaxed_length_min = 2
            self.aug_dedup_bin_bp = 100
            self.gap_search_window = 50000
            self.gap_min_size = 10
            self.gap_evalue = 10.0
            self.gap_min_identity = 25.0
            self.gap_min_alnlen = 10
            self.gap_max_hits = 5
            self.min_gene_identity = 30
            self.min_exon_query_cov = 0.4
            self.min_exon_alnlen = 5

    runner_args = MockArgs()
    
    # 1. MMseqs Search (Genome vs Queries)
    print(f"Searching {args.genome} with {query_fasta}...")
    hits_file = os.path.join(args.outdir, "initial_hits.m8")
    tmp_dir = os.path.join(args.outdir, "tmp_mmseqs")
    
    # Clean up previous run
    if os.path.exists(tmp_dir): shutil.rmtree(tmp_dir)
    
    cmd = [
        "mmseqs", "easy-search",
        query_fasta, args.genome, hits_file, tmp_dir,
        "--search-type", "2",
        "--threads", "4",
        "-s", str(runner_args.mmseqs_sens),
        "-e", str(runner_args.evalue),
        "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    ]
    
    subprocess.run(cmd, check=True)
    
    # 2. Parse Hits
    hits = parse_hits(hits_file, runner_args.min_identity, runner_args.min_length, runner_args.evalue)
    print(f"Found {len(hits)} initial hits")
    
    if not hits:
        print("No hits found! Exiting.")
        return

    # 3. Identify Synteny Block
    print(f"\n[Analysis] Analyzing {len(hits)} hits for synteny...")
    
    # --- GLOBAL EVALUATION ---
    # Check if the raw hits (after basic parsing/filtering) contain the ground truth
    # This checks if MMseqs FOUND them, regardless of whether the pipeline SELECTS them.
    from copy import deepcopy
    global_hits = deepcopy(hits)
    
    # Filter global hits by toxin names (for fairness)
    toxin_names = set()
    if os.path.exists(args.ground_truth):
        with open(args.ground_truth) as f:
            for line in f:
                if "name=" in line:
                    match = re.search(r'name=([^;]+)', line)
                    if match:
                        toxin_names.add(match.group(1))
    
    if toxin_names:
        global_hits = [h for h in global_hits if any(n in h['query'] for n in toxin_names)]
    
    print(f"[Global Check] {len(global_hits)} toxin hits found by MMseqs (genome-wide).")
    
    if os.path.exists(args.ground_truth):
        print("\n[Global Check] Comparison (Ignoring Region Selection):")
        compare_annotations(global_hits, args.ground_truth, offset=args.offset)
    # -------------------------

    # 3. Identify Synteny Blocks
    # Match pipeline intent: use flanking genes as anchors, then search GOI inside regions.
    flanking_seed_hits = [h for h in hits if not h.get('query', '').startswith('GOI_')]
    goi_seed_hits = [h for h in hits if h.get('query', '').startswith('GOI_')]
    if flanking_seed_hits:
        seed_hits = flanking_seed_hits
        print(f"[Analysis] Using {len(flanking_seed_hits)} flanking-anchor hits for block seeding.")
    elif goi_seed_hits:
        seed_hits = goi_seed_hits
        print(f"[Analysis] No flanking anchors found; using {len(goi_seed_hits)} GOI hits for block seeding.")
    else:
        seed_hits = hits
    synteny_blocks = identify_synteny_blocks(seed_hits, max_intron=runner_args.max_intron, cluster_dist=runner_args.cluster_dist)
    
    if not synteny_blocks:
        print("No valid synteny blocks found.")
        return

    print(f"[Analysis] Found {len(synteny_blocks)} synteny blocks.")

    # Prepare DB sequences
    db_sequences = {}
    for h, i, s in parse_fasta(query_fasta):
        db_sequences[i] = {'id': i, 'seq': s, 'header': h}
        base = extract_base_gene_id(i)
        if base != i:
            db_sequences[base] = {'id': i, 'seq': s, 'header': h}

    genome_seqs = load_genome(args.genome)
    
    # 4. Process Blocks
    all_gff_lines = []
    
    # Filter global hits by toxin names FIRST to define the set of interesting toxins
    toxin_names = set()
    if os.path.exists(args.ground_truth):
        with open(args.ground_truth) as f:
            for line in f:
                if "name=" in line:
                    match = re.search(r'name=([^;]+)', line)
                    if match:
                        toxin_names.add(match.group(1))
    print(f"Loaded {len(toxin_names)} toxin names from ground truth for filtering.")

    for i, block in enumerate(synteny_blocks):
        print(f"  Block {i+1}: {block['chrom']}:{block['start']}-{block['end']} ({block['genes_count']} genes)")
        try:
             # Call helper
            _, valid_gff = process_region_block(i, block, hits, genome_seqs, db_sequences, "test_genome", runner_args, "test_run", 4)
            all_gff_lines.extend(valid_gff)
        except Exception as e:
            print(f"Error processing block {i}: {e}")

    # Match pipeline behavior: flanking duplicates (same parent protein ID) are collapsed.
    _, all_gff_lines = deduplicate_flanking_models(
        [],
        all_gff_lines,
        genome_name="test_genome",
        locus_gap_bp=max(5000, int(runner_args.cluster_dist)),
    )
    all_gff_lines = collapse_flanking_cds_to_gene_span(all_gff_lines)
        
    # 5. Parse GFF lines to hits for comparison
    final_hits = []
    all_cds_hits = []
    all_gene_hits = []
    
    for line in all_gff_lines:
        parts = line.split('\t')
        if len(parts) >= 9 and parts[2] in {'mRNA', 'CDS'}:
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attr = parts[8]
            
            # Extract ID or Name
            name = "unknown"
            if "Name=" in attr:
                m = re.search(r'Name=([^;]+)', attr)
                if m: name = m.group(1)
            elif "ID=" in attr:
                 m = re.search(r'ID=([^;]+)', attr)
                 if m: name = m.group(1)
            
            if parts[2] == 'mRNA':
                all_gene_hits.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'query': name
                })
                continue

            all_cds_hits.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'query': name
            })

            is_toxin = False
            if toxin_names:
                for t in toxin_names:
                    if t in name:
                        is_toxin = True
                        break
            else:
                is_toxin = True 
            
            # Keep benchmark metrics focused on GOI predictions.
            if is_toxin and name.startswith("GOI_"):
                final_hits.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'query': name
                })

    aug_hits = final_hits
    print(f"Collected {len(aug_hits)} final predicted toxin hits across all blocks.")
    print(f"Collected {len(all_cds_hits)} total CDS annotations (toxin + flanking).")
    print(f"Collected {len(all_gene_hits)} total gene-span annotations (toxin + flanking).")
    # 6. Output BED
    bed_file = os.path.join(args.outdir, "predicted_annotations.bed")
    write_bed(aug_hits, bed_file, chrom_offset=args.offset, chrom_override=truth_chrom)
    print(f"Written predictions to {bed_file}")

    all_bed_file = os.path.join(args.outdir, "predicted_annotations_all.bed")
    write_bed(all_gene_hits, all_bed_file, chrom_offset=args.offset, chrom_override=truth_chrom)
    print(f"Written all gene-span annotations to {all_bed_file}")

    all_cds_bed_file = os.path.join(args.outdir, "predicted_annotations_all_cds.bed")
    write_bed(all_cds_hits, all_cds_bed_file, chrom_offset=args.offset, chrom_override=truth_chrom)
    print(f"Written all CDS annotations to {all_cds_bed_file}")
    
    # 7. Compare
    if os.path.exists(args.ground_truth):
        results = compare_annotations(aug_hits, args.ground_truth, offset=args.offset)
        
        # Save metrics
        metrics_file = os.path.join(args.outdir, "metrics.json")
        with open(metrics_file, 'w') as f:
            json.dump(results, f, indent=2)
    else:
        print("Ground truth file not found, skipping comparison.")

if __name__ == "__main__":
    main()
