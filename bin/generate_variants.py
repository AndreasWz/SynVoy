#!/usr/bin/env python3
"""
generate_variants.py - Generate protein variants and fragments for augmented search

This script generates:
1. Mutated variants using BLOSUM-like substitution groups
2. Systematic fragments (halves, thirds, quarters) for detecting partial matches
3. Overlapping fragments for edge-case detection

Usage:
    python generate_variants.py --query input.faa --output variants.faa \
        --num_variants 10 --mutation_rate 0.05
"""

import argparse
import random
import os
import sys

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import parse_fasta, write_fasta
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta

# Simplified substitution groups (BLOSUM-like)
# AA within same group are likely substitutions
SUB_GROUPS = [
    "LIVM", # Hydrophobic
    "FYW",  # Aromatic
    "DE",   # Acidic
    "KR",   # Basic
    "ST",   # Hydroxyl
    "NQ",   # Amide
    "G",    # Small
    "A",    # Small
    "P",    # Proline
    "C",    # Cysteine
    "H"     # Histidine
]

def parse_args():
    parser = argparse.ArgumentParser(description="Generate protein variants using substitution groups")
    parser.add_argument("--query", required=True, help="Input query protein FASTA")
    parser.add_argument("--output", required=True, help="Output variants FASTA")
    parser.add_argument("--num_variants", type=int, default=10, help="Number of mutant variants to generate")
    parser.add_argument("--mutation_rate", type=float, default=0.05, help="Mutation rate (fractions of residues)")
    parser.add_argument("--min_fragment_size", type=int, default=20, help="Minimum fragment size in AA")
    parser.add_argument("--no_fragments", action="store_true", help="Skip fragment generation")
    parser.add_argument("--no_mutations", action="store_true", help="Skip mutation variants")
    return parser.parse_args()

def get_substitution(aa):
    """Return a likely substitution for the given amino acid."""
    # Find group
    for group in SUB_GROUPS:
        if aa in group:
            if len(group) > 1:
                # Pick other AA from same group
                choices = [c for c in group if c != aa]
                return random.choice(choices)
            else:
                # Unique AA (e.g. Proline), maybe keep or very rare random?
                # For now, 90% keep, 10% completely random
                if random.random() < 0.9:
                   return aa
                else: 
                   return random.choice("ACDEFGHIKLMNPQRSTVWY")
    
    # Fallback random
    return random.choice("ACDEFGHIKLMNPQRSTVWY")

def mutate_sequence(seq, rate):
    mutated = list(seq)
    length = len(seq)
    num_mutations = int(length * rate)
    if num_mutations == 0: num_mutations = 1
    
    positions = random.sample(range(length), num_mutations)
    
    for pos in positions:
        curr = mutated[pos]
        new_aa = get_substitution(curr)
        mutated[pos] = new_aa
        
    return "".join(mutated)


def generate_systematic_fragments(seq, seq_id, min_size=20):
    """
    Generate systematic sequence fragments at different granularities.
    
    Args:
        seq: Sequence string
        seq_id: Original sequence ID
        min_size: Minimum fragment size in amino acids
    
    Returns:
        List of (id, seq) tuples for each fragment
    """
    fragments = []
    length = len(seq)
    
    # Halves
    if length // 2 >= min_size:
        half = length // 2
        fragments.append((
            f"{seq_id}|frag_half_1",
            seq[:half]
        ))
        fragments.append((
            f"{seq_id}|frag_half_2",
            seq[half:]
        ))
    
    # Thirds
    if length // 3 >= min_size:
        third = length // 3
        fragments.append((
            f"{seq_id}|frag_third_1",
            seq[:third]
        ))
        fragments.append((
            f"{seq_id}|frag_third_2",
            seq[third:2*third]
        ))
        fragments.append((
            f"{seq_id}|frag_third_3",
            seq[2*third:]
        ))
    
    # Quarters
    if length // 4 >= min_size:
        quarter = length // 4
        for i in range(4):
            start = i * quarter
            end = (i + 1) * quarter if i < 3 else length
            fragments.append((
                f"{seq_id}|frag_quarter_{i+1}",
                seq[start:end]
            ))
    
    return fragments
    
    return fragments


def main():
    args = parse_args()
    
    # parse_fasta yields (raw_header, clean_id, seq) - we want (clean_id, seq)
    records = [(clean_id, seq) for _, clean_id, seq in parse_fasta(args.query)]
    if not records:
        print("Error: No sequences found.")
        return
        
    output_records = []
    
    for rec_id, rec_seq in records:
        # Always include original
        output_records.append((rec_id, rec_seq))
        
        # Generate mutation variants
        if not args.no_mutations:
            for i in range(args.num_variants):
                mut_seq = mutate_sequence(rec_seq, args.mutation_rate)
                mut_id = f"{rec_id}|var_{i+1}"
                output_records.append((mut_id, mut_seq))
        
        # Generate systematic fragments
        if not args.no_fragments:
            fragments = generate_systematic_fragments(
                rec_seq, 
                rec_id, 
                min_size=args.min_fragment_size
            )
            output_records.extend(fragments)
            
            # Also add overlapping legacy fragments for backwards compatibility
            slen = len(rec_seq)
            if slen > 50:
                # Split in 3 overlapping parts (legacy behavior)
                part1 = rec_seq[:slen//2 + 10]
                part2 = rec_seq[slen//3 : 2*slen//3 + 10]
                part3 = rec_seq[slen//2 - 10:]
                
                output_records.append((f"{rec_id}|overlap_1", part1))
                output_records.append((f"{rec_id}|overlap_2", part2))
                output_records.append((f"{rec_id}|overlap_3", part3))

    write_fasta(output_records, args.output)
    print(f"Generated {len(output_records)} variants/fragments from {len(records)} input sequence(s).")


if __name__ == "__main__":
    main()
