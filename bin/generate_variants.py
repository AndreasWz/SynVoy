#!/usr/bin/env python3
import argparse
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

def main():
    args = parse_args()
    
    records = list(SeqIO.parse(args.query, "fasta"))
    if not records:
        print("Error: No sequences found.")
        return
        
    output_records = []
    
    for rec in records:
        output_records.append(rec)
        
        # Variants
        for i in range(args.num_variants):
            mut_seq = mutate_sequence(str(rec.seq), args.mutation_rate)
            mut_id = f"{rec.id}_var_{i+1}"
            output_records.append(SeqRecord(Seq(mut_seq), id=mut_id, description="simulated variant"))
            
        # Fragments (Domain shuffling/split simulation)
        slen = len(rec.seq)
        if slen > 50:
             # Split in 3 overlapping parts
             part1 = rec.seq[:slen//2 + 10]
             part2 = rec.seq[slen//3 : 2*slen//3 + 10]
             part3 = rec.seq[slen//2 - 10:]
             
             output_records.append(SeqRecord(part1, id=f"{rec.id}_part1", description="fragment"))
             output_records.append(SeqRecord(part2, id=f"{rec.id}_part2", description="fragment"))
             output_records.append(SeqRecord(part3, id=f"{rec.id}_part3", description="fragment"))

    SeqIO.write(output_records, args.output, "fasta")
    print(f"Generated {len(output_records)} variants/fragments.")

if __name__ == "__main__":
    main()
