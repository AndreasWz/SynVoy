#!/usr/bin/env python3
"""Rewrite an MMseqs2 table written with the two header columns (sequence_utils.MMSEQS_FORMAT_OUTPUT)
into the 12-column BLAST-tabular table the pipeline's readers parse, the query and the target named by
the first word of their FASTA headers: MMseqs2's own name columns rewrite names that start with a
database tag or with uc or consensus_ (lcl|chr1 as chr1, gnl|a|1 and gnl|b|1 both as 1; fix item 29).

    mmseqs_first_words.py --input hits_mmseqs.raw --output hits_mmseqs.m8
"""
import argparse
import os
import sys

try:
    from sequence_utils import exit_step_failed, mmseqs_table_first_words
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import exit_step_failed, mmseqs_table_first_words


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--input", required=True, help="MMseqs2 table with --format-output ...,bits,qheader,theader")
    ap.add_argument("--output", required=True, help="the 12-column table, names = first words of the headers")
    args = ap.parse_args()
    n = mmseqs_table_first_words(args.input, args.output)
    print(f"[mmseqs_first_words] {n} row(s) written to {args.output}")
    return 0


if __name__ == "__main__":
    exit_step_failed(main)
