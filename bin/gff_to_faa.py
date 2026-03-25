#!/usr/bin/env python3

import argparse
import os
import sys

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import parse_fasta, write_fasta, load_genome, reverse_complement, translate
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta, load_genome, reverse_complement, translate

def parse_gff(gff_file):
    relationships = {} # child -> parent
    mRNAs = {} # mRNA_id -> info
    
    if not os.path.exists(gff_file):
        return {}
        
    # First pass: collect relationships
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            attr_str = parts[8].strip()
            attrs = {}
            for item in attr_str.split(';'):
                if '=' in item:
                    k, v = item.strip().split('=', 1)
                    attrs[k] = v
            
            fid = attrs.get('ID')
            parent = attrs.get('Parent')
            if fid and parent:
                relationships[fid] = parent
    
    # Second pass: collect CDS for each mRNA
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            feat_type = parts[2]
            if feat_type != 'CDS':
                continue
                
            attr_str = parts[8].strip()
            attrs = {}
            for item in attr_str.split(';'):
                if '=' in item:
                    k, v = item.strip().split('=', 1)
                    attrs[k] = v
            
            mrna_id = attrs.get('Parent')
            if not mrna_id:
                continue
                
            if mrna_id not in mRNAs:
                mRNAs[mrna_id] = {
                    'chrom': parts[0],
                    'strand': parts[6],
                    'exons': []
                }
            mRNAs[mrna_id]['exons'].append((int(parts[3]), int(parts[4]),
                                              int(parts[7]) if parts[7] != '.' else 0))

    # Third pass: Group mRNAs by gene and pick longest
    gene_to_longest_mrna = {}
    
    for mrna_id, info in mRNAs.items():
        # Find the 'gene' ancestor
        ancestor = mrna_id
        visited = set()
        while ancestor in relationships and ancestor not in visited:
            visited.add(ancestor)
            if ancestor.startswith('gene-') or ancestor.startswith('gene:'):
                break
            ancestor = relationships[ancestor]
        
        gene_id = ancestor
        length = sum(e[1] - e[0] + 1 for e in info['exons'])
        
        if gene_id not in gene_to_longest_mrna or length > gene_to_longest_mrna[gene_id]['length']:
            gene_to_longest_mrna[gene_id] = {
                'id': mrna_id,
                'length': length,
                'info': info
            }
                
    return {gid: d['info'] for gid, d in gene_to_longest_mrna.items()}

def main():
    parser = argparse.ArgumentParser(description="Extract proteins from GFF and Genome")
    parser.add_argument("--gff", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    
    genome = load_genome(args.genome)
    genes = parse_gff(args.gff)
    
    records = []  # List of (id, seq) tuples
    for gene_id, info in genes.items():
        chrom = info['chrom']
        if chrom not in genome:
            continue
            
        exons = sorted(info['exons'], key=lambda x: x[0])
        dna_seq = ""
        # GFF phase = number of bases to skip at the 5' coding end.
        # For + strand: first CDS genomically is first in coding order → trim from its left.
        # For - strand: last CDS genomically is first in coding order → trim after revcomp.
        # Collect phase from the first CDS in coding direction; apply after concat+revcomp.
        coding_first_phase = exons[0][2] if info['strand'] == '+' else exons[-1][2]
        for exon in exons:
            start, end = exon[0], exon[1]
            exon_seq = genome[chrom][start-1:end]
            dna_seq += exon_seq

        if info['strand'] == '-':
            dna_seq = reverse_complement(dna_seq)

        # Apply phase: trim leading bases at coding 5' end
        if coding_first_phase > 0:
            dna_seq = dna_seq[coding_first_phase:]
            
        try:
            prot_seq = translate(dna_seq)
            # Stop at first stop codon
            if '*' in prot_seq:
                prot_seq = prot_seq.split('*')[0]
            if len(prot_seq) > 0:
                records.append((gene_id, prot_seq))
        except Exception as e:
            print(f"Warning: Failed to translate {gene_id}: {e}", file=sys.stderr)
        
    write_fasta(records, args.output)

if __name__ == "__main__":
    main()
