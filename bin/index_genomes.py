#!/usr/bin/env python3
"""Index every target genome of a folder once, before the steps that read it run side by side.

    index_genomes.py --dir staged_genomes [--samtools samtools]

The rescue passes call ``samtools faidx`` on the genomes of the shared target folder, and
samtools writes ``<genome>.fai`` beside the file it is given when there is none: several rescue
tasks then wrote the same index at once while others read it (fix-phase d5, item 21). Each
FASTA (.fa .fna .fasta, also .gz) gets its index here, in genome order. A bgzip file is
indexed too (.fai and .gzi); a plain-gzip file cannot be (samtools refuses it), and the rescues
read it without samtools (sequence_utils.genome_window): it is named on stderr. Without
samtools nothing is indexed, with one line. Exit 0 unless the folder is missing."""
import argparse
import os
import shutil
import subprocess
import sys

FASTA_EXT = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True, help="the folder of target genomes")
    ap.add_argument("--samtools", default="samtools")
    args = ap.parse_args()
    if not os.path.isdir(args.dir):
        sys.exit(f"index_genomes.py: no folder {args.dir}")
    if not shutil.which(args.samtools):
        print(f"[index_genomes] {args.samtools} not found: the genomes are not indexed, and the rescues "
              "read them without it.", file=sys.stderr)
        return 0
    for name in sorted(os.listdir(args.dir)):
        path = os.path.join(args.dir, name)
        if not name.endswith(FASTA_EXT) or not os.path.isfile(path):
            continue
        proc = subprocess.run([args.samtools, "faidx", path], capture_output=True, text=True)
        if proc.returncode != 0:
            for stale in (path + ".fai", path + ".gzi"):
                if os.path.exists(stale):
                    os.remove(stale)
            why = (proc.stderr.strip().splitlines() or [""])[-1][:200]
            print(f"[index_genomes] {name}: no index ({'a plain-gzip file' if name.endswith('.gz') else why}); "
                  "the rescues read it without samtools.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
