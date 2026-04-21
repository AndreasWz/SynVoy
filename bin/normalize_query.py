#!/usr/bin/env python3
"""
normalize_query.py - Ensure query is protein FASTA.

If the input appears nucleotide-only, translate in 6 frames and
select the longest ORF (longest segment between stop codons).
Otherwise, pass through the original protein sequence.
"""

import argparse
import os
import sys

try:
    from sequence_utils import parse_fasta, write_fasta, translate, reverse_complement
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta, translate, reverse_complement


NUC_ALPHABET = set("ACGTNRYWSKMBDHV-")


def is_nucleotide(seq: str) -> bool:
    seq = seq.upper()
    if not seq:
        return False
    return all(ch in NUC_ALPHABET for ch in seq)


def best_orf_from_protein(prot: str) -> str:
    """Return the longest segment between stop codons."""
    if not prot:
        return ""
    segments = prot.split("*")
    if not segments:
        return ""
    return max(segments, key=len)


def translate_best_orf(dna: str) -> str:
    """Translate DNA in 6 frames and return the longest ORF protein."""
    dna = dna.upper().replace("U", "T")
    best = ""

    # Forward frames
    for frame in range(3):
        prot = translate(dna[frame:])
        orf = best_orf_from_protein(prot)
        if len(orf) > len(best):
            best = orf

    # Reverse frames
    rev = reverse_complement(dna)
    for frame in range(3):
        prot = translate(rev[frame:])
        orf = best_orf_from_protein(prot)
        if len(orf) > len(best):
            best = orf

    return best


def main():
    parser = argparse.ArgumentParser(description="Normalize query FASTA to protein")
    parser.add_argument("--input", required=True, help="Input FASTA (DNA or protein)")
    parser.add_argument("--output", required=True, help="Output protein FASTA")
    parser.add_argument("--min_length", type=int, default=30,
                        help="Minimum protein length (aa). Queries below this "
                             "are rejected — searches on very short queries "
                             "return noise. Set to 0 to disable.")
    args = parser.parse_args()

    records = list(parse_fasta(args.input))
    if not records:
        print("ERROR: No sequences found in query FASTA", file=sys.stderr)
        sys.exit(1)

    if len(records) > 1:
        print("WARNING: Multiple sequences found; using the first record only.", file=sys.stderr)

    header, clean_id, seq = records[0]

    if is_nucleotide(seq):
        print(f"[normalize_query] Detected nucleotide query: {clean_id}", file=sys.stderr)
        prot = translate_best_orf(seq)
        if not prot:
            print("ERROR: Could not translate nucleotide query into a protein ORF", file=sys.stderr)
            sys.exit(1)
        final_seq = prot
        print(f"[normalize_query] Translated query length: {len(prot)} aa", file=sys.stderr)
    else:
        print(f"[normalize_query] Detected protein query: {clean_id}", file=sys.stderr)
        final_seq = seq

    if args.min_length > 0 and len(final_seq) < args.min_length:
        print(
            f"ERROR: Query '{clean_id}' is {len(final_seq)} aa, below the "
            f"minimum of {args.min_length} aa. Short queries produce noisy "
            f"MMseqs2 / tblastn hits. Override with --min_length 0 if you "
            f"really want to search with a fragment this short.",
            file=sys.stderr,
        )
        sys.exit(2)

    write_fasta([(clean_id, final_seq)], args.output)


if __name__ == "__main__":
    main()
