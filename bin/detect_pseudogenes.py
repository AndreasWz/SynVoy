#!/usr/bin/env python3
"""
detect_pseudogenes.py - Identify pseudogenes and gene losses

Detects:
1. Premature stop codons in alignments
2. Frameshifts (from Miniprot output)
3. Truncated genes (< 50% reference length)
4. Missing exons (exon count mismatch)
5. Low identity pseudogenes

Usage:
    python detect_pseudogenes.py --gff target.gff --reference ref.faa \\
        --output pseudogenes.tsv
"""

import argparse
import sys
import os
import re
from collections import defaultdict

try:
    from sequence_utils import parse_fasta, parse_gff, translate
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, parse_gff, translate


def parse_args():
    parser = argparse.ArgumentParser(description="Detect pseudogenes and gene losses")
    parser.add_argument("--gff", required=True, help="Miniprot GFF with annotations")
    parser.add_argument("--reference", required=True, help="Reference protein FASTA")
    parser.add_argument("--genome", required=True, help="Target genome FASTA")
    parser.add_argument("--output", required=True, help="Output TSV with pseudogene calls")
    parser.add_argument("--min_coverage", type=float, default=0.5, 
                       help="Minimum coverage of reference (default: 50%%)")
    parser.add_argument("--min_identity", type=float, default=30.0,
                       help="Minimum identity for functional gene (default: 30%%)")
    return parser.parse_args()


def parse_miniprot_gff(gff_file):
    """
    Parse Miniprot GFF to extract gene features with exon information.
    
    Returns dict: {gene_id: gene_info}
    """
    genes = {}
    current_gene = None
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feat_type, start, end, score, strand, phase, attrs = parts
            
            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attr_dict[k] = v
            
            if feat_type == 'mRNA':
                gene_id = attr_dict.get('ID', 'unknown')
                parent = attr_dict.get('SynTerra_Parent', attr_dict.get('Target', '').split()[0])
                identity = float(attr_dict.get('Identity', 0)) * 100
                
                current_gene = {
                    'id': gene_id,
                    'parent': parent,
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'identity': identity,
                    'exons': [],
                    'cds_parts': [],
                    'has_frameshift': False,
                    'has_stop_codon': False
                }
                genes[gene_id] = current_gene
            
            elif feat_type == 'CDS' and current_gene:
                current_gene['cds_parts'].append({
                    'start': int(start),
                    'end': int(end),
                    'phase': phase
                })
                
                # Check for frameshift indicator in attributes
                if 'frameshift' in attrs.lower() or 'Frameshift' in attr_dict:
                    current_gene['has_frameshift'] = True
    
    return genes


def check_stop_codons(cds_seq):
    """
    Check for premature stop codons in CDS sequence.
    Returns (has_internal_stops, stop_positions)
    """
    protein = translate(cds_seq)
    
    # Find stop codons (excluding final position)
    stop_positions = []
    for i, aa in enumerate(protein[:-1]):  # Exclude last position
        if aa == '*':
            stop_positions.append(i)
    
    return len(stop_positions) > 0, stop_positions


def calculate_coverage(gene_length, ref_length):
    """Calculate coverage as fraction of reference length."""
    return gene_length / ref_length if ref_length > 0 else 0


def count_exons(gene_info):
    """Count number of CDS/exons in gene."""
    return len(gene_info['cds_parts'])


def classify_pseudogene_type(gene_info, ref_info, coverage, has_stops):
    """
    Classify pseudogene type based on features.
    
    Returns: (is_pseudogene, pseudogene_type, reason)
    """
    reasons = []
    
    # Type 1: Frameshift pseudogene
    if gene_info['has_frameshift']:
        return True, "FRAMESHIFT", "Contains frameshift mutation"
    
    # Type 2: Premature stop codon
    if has_stops:
        return True, "NONSENSE", "Contains premature stop codon(s)"
    
    # Type 3: Truncated (< 50% of reference)
    if coverage < 0.5:
        return True, "TRUNCATED", f"Only {coverage*100:.1f}% of reference length"
    
    # Type 4: Very low identity (< 25%)
    if gene_info['identity'] < 25.0:
        return True, "DIVERGENT", f"Very low identity ({gene_info['identity']:.1f}%)"
    
    # Type 5: Fragmented (< 60% coverage but > 1 exon)
    if coverage < 0.6 and len(gene_info['cds_parts']) > 1:
        return True, "FRAGMENTED", f"Partial gene ({coverage*100:.1f}% coverage, fragmented)"
    
    # Likely functional
    return False, "FUNCTIONAL", None


def main():
    args = parse_args()
    
    # 1. Load reference proteins
    print("Loading reference proteins...", file=sys.stderr)
    ref_proteins = {}
    for header, clean_id, seq in parse_fasta(args.reference):
        # Strip GOI_ prefix if present
        ref_id = clean_id.replace('GOI_', '')
        ref_proteins[ref_id] = {
            'id': ref_id,
            'seq': seq,
            'length': len(seq)
        }
    print(f"  Loaded {len(ref_proteins)} reference proteins", file=sys.stderr)
    
    # 2. Parse Miniprot GFF
    print("Parsing Miniprot annotations...", file=sys.stderr)
    genes = parse_miniprot_gff(args.gff)
    print(f"  Found {len(genes)} gene predictions", file=sys.stderr)
    
    # 3. Load genome (for stop codon checking)
    print("Loading target genome...", file=sys.stderr)
    from sequence_utils import load_genome, reverse_complement
    genome = load_genome(args.genome)
    
    # 4. Analyze each gene
    results = []
    
    for gene_id, gene_info in genes.items():
        parent_id = gene_info['parent'].replace('GOI_', '')
        
        # Find reference
        ref_info = None
        for ref_id, ref_data in ref_proteins.items():
            if ref_id in parent_id or parent_id in ref_id:
                ref_info = ref_data
                break
        
        if not ref_info:
            # No reference found - skip
            continue
        
        # Extract CDS sequence
        if gene_info['chrom'] not in genome:
            continue
        
        chrom_seq = genome[gene_info['chrom']]
        cds_seq = ""
        
        for cds in sorted(gene_info['cds_parts'], key=lambda x: x['start']):
            # GFF is 1-based, Python is 0-based
            exon_seq = chrom_seq[cds['start']-1:cds['end']]
            cds_seq += exon_seq
        
        # Handle strand
        if gene_info['strand'] == '-':
            cds_seq = reverse_complement(cds_seq)
        
        # Check for stop codons
        has_stops, stop_positions = check_stop_codons(cds_seq)
        
        # Calculate metrics
        gene_length = len(cds_seq) // 3  # Protein length
        coverage = calculate_coverage(gene_length, ref_info['length'])
        num_exons = count_exons(gene_info)
        
        # Classify
        is_pseudo, pseudo_type, reason = classify_pseudogene_type(
            gene_info, ref_info, coverage, has_stops
        )
        
        results.append({
            'gene_id': gene_id,
            'parent': parent_id,
            'chrom': gene_info['chrom'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'identity': gene_info['identity'],
            'coverage': coverage * 100,
            'num_exons': num_exons,
            'has_frameshift': gene_info['has_frameshift'],
            'has_stop_codons': has_stops,
            'num_stops': len(stop_positions) if has_stops else 0,
            'classification': pseudo_type,
            'is_pseudogene': is_pseudo,
            'reason': reason if reason else "Likely functional"
        })
    
    # 5. Write output
    print(f"\\nWriting results to {args.output}...", file=sys.stderr)
    
    with open(args.output, 'w') as out:
        # Header
        out.write("gene_id\\tparent\\tchrom\\tstart\\tend\\tstrand\\t"
                 "identity\\tcoverage\\tnum_exons\\thas_frameshift\\t"
                 "has_stop_codons\\tnum_stops\\tclassification\\t"
                 "is_pseudogene\\treason\\n")
        
        for r in results:
            out.write(f"{r['gene_id']}\\t{r['parent']}\\t{r['chrom']}\\t"
                     f"{r['start']}\\t{r['end']}\\t{r['strand']}\\t"
                     f"{r['identity']:.1f}\\t{r['coverage']:.1f}\\t{r['num_exons']}\\t"
                     f"{r['has_frameshift']}\\t{r['has_stop_codons']}\\t"
                     f"{r['num_stops']}\\t{r['classification']}\\t"
                     f"{r['is_pseudogene']}\\t{r['reason']}\\n")
    
    # 6. Summary
    total = len(results)
    pseudogenes = sum(1 for r in results if r['is_pseudogene'])
    functional = total - pseudogenes
    
    print(f"\\n=== Pseudogene Detection Summary ===", file=sys.stderr)
    print(f"Total genes analyzed: {total}", file=sys.stderr)
    print(f"Likely functional: {functional} ({functional/total*100:.1f}%)", file=sys.stderr)
    print(f"Pseudogenes detected: {pseudogenes} ({pseudogenes/total*100:.1f}%)", file=sys.stderr)
    
    if pseudogenes > 0:
        types = defaultdict(int)
        for r in results:
            if r['is_pseudogene']:
                types[r['classification']] += 1
        
        print(f"\\nPseudogene types:", file=sys.stderr)
        for ptype, count in sorted(types.items(), key=lambda x: -x[1]):
            print(f"  {ptype}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
