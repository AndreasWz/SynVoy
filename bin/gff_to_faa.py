#!/usr/bin/env python3

import argparse
import os
import sys

# Use our own sequence utilities (no BioPython)
try:
    from sequence_utils import (cds_model_protein, cds_models_by_gene, load_genome,
                                set_genetic_code, write_fasta)
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (cds_model_protein, cds_models_by_gene, load_genome,
                                set_genetic_code, write_fasta)


def main():
    parser = argparse.ArgumentParser(description="Extract proteins from GFF and Genome")
    parser.add_argument("--gff", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--coords_tsv",
                        help="Also write gene_id, chrom, start, end, strand (1-based, coding span) "
                             "of every protein written, for the search's specificity check")
    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the genome (default 1, the standard code)")
    args = parser.parse_args()
    set_genetic_code(args.genetic_code)
    
    genome = load_genome(args.genome)
    genes = cds_models_by_gene(args.gff)

    records = []  # List of (id, seq) tuples
    coords = []   # (id, chrom, start, end, strand) of the written proteins
    for gene_id, info in genes.items():
        chrom = info['chrom']
        if chrom not in genome:
            continue
        exons = sorted(info['exons'], key=lambda x: x[0])
        try:
            prot_seq = cds_model_protein(info, genome[chrom])
            if len(prot_seq) > 0:
                records.append((gene_id, prot_seq))
                coords.append((gene_id, chrom, exons[0][0], max(e[1] for e in exons), info['strand']))
        except Exception as e:
            print(f"Warning: Failed to translate {gene_id}: {e}", file=sys.stderr)

    write_fasta(records, args.output)
    if args.coords_tsv:
        with open(args.coords_tsv, "w") as fh:
            fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
            for row in coords:
                fh.write("\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()
