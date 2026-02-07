#!/usr/bin/env python3
"""
Integrate augmented search results as GOI annotations.

This script takes candidates found by augmented search (using variants + MMseqs2)
and creates proper gene annotations for them. These become the new GOI sequences
for the next iteration.

The key insight: augmented search finds hits that miniprot might miss due to
high divergence or short length. We need to annotate these and add them to
the expanding database.
"""

import argparse
import sys
import os
from typing import List, Dict, Tuple

# Import sequence utilities
try:
    from sequence_utils import parse_fasta, write_fasta, load_genome, translate, reverse_complement
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta, load_genome, translate, reverse_complement


def annotate_candidates_simple(candidate_bed: str, genome_fasta: str, query_fasta: str) -> List[Dict]:
    """
    Create simple gene annotations from augmented search candidates.
    
    Strategy:
    1. Read BED coordinates of candidates
    2. Extract sequences from genome
    3. Translate in best frame (frame with longest ORF)
    4. Create gene annotation
    
    Args:
        candidate_bed: BED file with candidate regions
        genome_fasta: Target genome FASTA
        query_fasta: Original query for reference
        
    Returns:
        List of gene records with 'id', 'seq', 'chrom', 'start', 'end', 'strand'
    """
    genes = []
    
    # Load genome
    genome_seqs = load_genome(genome_fasta)
    
    # Get genome name from filename
    genome_name = os.path.basename(genome_fasta).replace('.fna', '').replace('.fasta', '').replace('.fa', '')
    
    # Read candidates
    candidates = []
    with open(candidate_bed) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            candidates.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3],
                'score': float(parts[4]),
                'strand': parts[5]
            })
    
    # Process each candidate
    for cand in candidates:
        chrom = cand['chrom']
        if chrom not in genome_seqs:
            continue
        
        # Extract sequence
        seq = genome_seqs[chrom][cand['start']:cand['end']]
        
        if cand['strand'] == '-':
            seq = reverse_complement(seq)
        
        # Find best frame (longest ORF without stops)
        best_frame = 0
        best_protein = ""
        best_length = 0
        
        for frame in range(3):
            frame_seq = seq[frame:]
            # Make it divisible by 3
            frame_seq = frame_seq[:len(frame_seq) - len(frame_seq) % 3]
            if not frame_seq:
                continue
                
            protein = translate(frame_seq)
            # Count length until first stop or end
            length = 0
            for aa in protein:
                if aa == '*':
                    break
                length += 1
            
            if length > best_length:
                best_length = length
                best_frame = frame
                best_protein = protein[:length]  # Exclude stop
        
        if best_protein and len(best_protein) >= 20:  # Min 20 amino acids
            gene_id = f"GOI_{genome_name}_{cand['name']}_f{best_frame}"
            
            genes.append({
                'id': gene_id,
                'seq': best_protein,
                'chrom': chrom,
                'start': cand['start'],
                'end': cand['end'],
                'strand': cand['strand'],
                'score': cand['score'],
                'method': 'augmented_simple'
            })
    
    return genes


def annotate_candidates_with_prodigal(candidate_bed: str, genome_fasta: str, output_dir: str) -> List[Dict]:
    """
    Use Prodigal for more sophisticated ORF prediction on candidate regions.
    
    This is more accurate than simple 6-frame translation but requires Prodigal.
    """
    import subprocess
    import tempfile
    
    genes = []
    
    # Load genome
    genome_seqs = load_genome(genome_fasta)
    genome_name = os.path.basename(genome_fasta).replace('.fna', '').replace('.fasta', '').replace('.fa', '')
    
    # Read candidates and extract regions
    candidates = []
    with open(candidate_bed) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            candidates.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3],
                'score': float(parts[4]),
                'strand': parts[5]
            })
    
    # Create temp FASTA with candidate regions
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fna', delete=False) as temp_fna:
        temp_fna_path = temp_fna.name
        for cand in candidates:
            if cand['chrom'] in genome_seqs:
                seq = genome_seqs[cand['chrom']][cand['start']:cand['end']]
                temp_fna.write(f">{cand['name']}\n{seq}\n")
    
    # Run Prodigal
    temp_faa = temp_fna_path.replace('.fna', '.faa')
    temp_gff = temp_fna_path.replace('.fna', '.gff')
    
    try:
        subprocess.run([
            'prodigal',
            '-i', temp_fna_path,
            '-a', temp_faa,
            '-f', 'gff',
            '-o', temp_gff,
            '-p', 'meta',  # Metagenomic mode for short sequences
            '-q'  # Quiet
        ], check=True, capture_output=True)
        
        # Parse Prodigal output
        for header, clean_id, seq in parse_fasta(temp_faa):
            # Map back to candidate
            for cand in candidates:
                if cand['name'] in clean_id:
                    gene_id = f"GOI_{genome_name}_{clean_id}"
                    genes.append({
                        'id': gene_id,
                        'seq': seq,
                        'chrom': cand['chrom'],
                        'start': cand['start'],
                        'end': cand['end'],
                        'strand': cand['strand'],
                        'score': cand['score'],
                        'method': 'augmented_prodigal'
                    })
                    break
        
        # Cleanup
        os.unlink(temp_fna_path)
        os.unlink(temp_faa)
        os.unlink(temp_gff)
        
    except subprocess.CalledProcessError:
        # Prodigal failed, return empty
        print("Warning: Prodigal failed, returning empty results", file=sys.stderr)
        os.unlink(temp_fna_path)
    
    return genes


def main():
    parser = argparse.ArgumentParser(
        description='Integrate augmented search results as GOI annotations'
    )
    parser.add_argument('--candidate_bed', required=True,
                       help='BED file with augmented search candidates')
    parser.add_argument('--genome', required=True,
                       help='Target genome FASTA')
    parser.add_argument('--query', required=True,
                       help='Original query FASTA (for reference)')
    parser.add_argument('--output', required=True,
                       help='Output FASTA with annotated GOI sequences')
    parser.add_argument('--method', choices=['simple', 'prodigal'], default='simple',
                       help='Annotation method')
    parser.add_argument('--output_dir', default='.',
                       help='Output directory for temp files')
    
    args = parser.parse_args()
    
    # Check if candidate BED is empty
    if not os.path.exists(args.candidate_bed) or os.path.getsize(args.candidate_bed) == 0:
        print("Candidate BED is empty. No annotations to create.")
        # Create empty output
        open(args.output, 'w').close()
        return
    
    # Annotate candidates
    if args.method == 'simple':
        genes = annotate_candidates_simple(args.candidate_bed, args.genome, args.query)
    else:
        genes = annotate_candidates_with_prodigal(args.candidate_bed, args.genome, args.output_dir)
    
    # Write output
    if genes:
        write_fasta([(g['id'], g['seq']) for g in genes], args.output)
        print(f"Annotated {len(genes)} GOI candidates")
    else:
        print("No candidates passed annotation filters")
        open(args.output, 'w').close()


if __name__ == '__main__':
    main()
