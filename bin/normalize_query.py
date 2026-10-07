#!/usr/bin/env python3
"""
normalize_query.py - Ensure query is a clean protein FASTA.

Formatting that is not sequence is removed first: case is folded to upper,
and whitespace, digits (numbered GenBank-style lines), alignment gaps
(``-`` ``.``) and a trailing ``*`` are dropped. Anything else that is not an
IUPAC letter is an error rather than something to search with.

If the cleaned input is nucleotide, translate in 6 frames and select the
longest ORF (longest segment between stop codons). Otherwise the cleaned
protein sequence is written out.
"""

import argparse
import os
import sys

try:
    from sequence_utils import parse_fasta, write_fasta, translate, reverse_complement
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, write_fasta, translate, reverse_complement


NUC_ALPHABET = set("ACGTUNRYWSKMBDHV")
# IUPAC amino acids incl. ambiguity codes (B, J, Z, X) and the rare U (Sec) / O (Pyl).
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")
_FORMATTING = set("-.")


def clean_sequence(seq: str):
    """Fold case and drop characters that are formatting rather than sequence.

    Returns ``(cleaned, n_removed)``. A query used to pass through verbatim, so a
    lowercase FASTA, an aligned sequence with gaps, or numbered GenBank-style lines
    ("1 mkflv nvalv ...") reached the search with those characters still in it --
    counted in the query length and compared residue-by-residue downstream.
    """
    out = []
    removed = 0
    for ch in seq.upper():
        if ch.isspace() or ch.isdigit() or ch in _FORMATTING:
            removed += 1
            continue
        out.append(ch)
    return "".join(out), removed


def is_nucleotide(seq: str) -> bool:
    """True for DNA/RNA. Every character must be an IUPAC nucleotide code AND at
    least 90 % must be plain A/C/G/T/U/N: the ambiguity codes alone overlap 14 amino
    acids, so a Gly/Ser-rich peptide such as "GSGSGSGS" is not nucleotide."""
    seq = seq.upper()
    if not seq:
        return False
    if not all(ch in NUC_ALPHABET for ch in seq):
        return False
    plain = sum(1 for ch in seq if ch in "ACGTUN")
    return plain >= 0.9 * len(seq)


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
        print(
            f"ERROR: No sequences found in query FASTA '{args.input}'. "
            f"The file is either empty or does not contain a valid '>header' "
            f"line. Try: verify the file with `head -2 {args.input}` and "
            f"ensure the first non-empty line starts with '>'.",
            file=sys.stderr,
        )
        sys.exit(1)

    header, clean_id, seq = records[0]

    if len(records) > 1:
        print(
            f"WARNING: {len(records)} sequences found in '{args.input}'; using only the "
            f"first ('{clean_id}'). SynVoy searches one gene of interest per run -- the "
            f"other {len(records) - 1} are ignored.",
            file=sys.stderr,
        )

    seq, n_removed = clean_sequence(seq)
    if n_removed:
        print(f"[normalize_query] Removed {n_removed} formatting character(s) "
              f"(whitespace, digits, alignment gaps) from '{clean_id}'", file=sys.stderr)
    if not seq:
        print(
            f"ERROR: Query '{clean_id}' in '{args.input}' has a header but no sequence.",
            file=sys.stderr,
        )
        sys.exit(1)

    if is_nucleotide(seq):
        print(f"[normalize_query] Detected nucleotide query: {clean_id}", file=sys.stderr)
        prot = translate_best_orf(seq)
        if not prot:
            print(
                f"ERROR: Could not translate nucleotide query '{clean_id}' "
                f"(length {len(seq)} nt) into a protein ORF. Every 6-frame "
                f"translation was empty, which usually means the sequence is "
                f"shorter than 3 nt or consists only of stop codons. "
                f"Try: supply the coding sequence (CDS) of your gene instead "
                f"of a partial genomic region, or translate it manually and "
                f"pass the protein FASTA directly.",
                file=sys.stderr,
            )
            sys.exit(1)
        final_seq = prot
        print(f"[normalize_query] Translated query length: {len(prot)} aa", file=sys.stderr)
    else:
        print(f"[normalize_query] Detected protein query: {clean_id}", file=sys.stderr)
        final_seq = seq.rstrip("*")          # a trailing stop symbol is not a residue
        if "*" in final_seq:
            n_stop = final_seq.count("*")
            print(
                f"WARNING: Query '{clean_id}' contains {n_stop} internal stop symbol(s) "
                f"('*'); they were removed. If this is a pseudogene or several ORFs "
                f"joined together, the search results will reflect that.",
                file=sys.stderr,
            )
            final_seq = final_seq.replace("*", "")
        bad = sorted(set(final_seq) - PROTEIN_ALPHABET)
        if bad:
            print(
                f"ERROR: Query '{clean_id}' contains character(s) that are not amino-acid "
                f"codes: {' '.join(repr(c) for c in bad)}. Check that '{args.input}' is "
                f"a plain FASTA file (one '>header' line, then sequence only).",
                file=sys.stderr,
            )
            sys.exit(1)

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
