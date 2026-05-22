#!/usr/bin/env python3
"""Sanity-check that fetched benchmark GFFs parse cleanly with the same logic
SynVoy uses for flanking-gene extraction.

For each (species_dir) under local_data/benchmark/genomes/, count:
- gene features
- mRNA features
- CDS features
- which attribute scheme is used (ID=/Parent= chains, RefSeq locus_tag, etc.)

Fails loud if any GFF has zero gene features (which would break SynVoy).
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


def scan_gff(path: Path) -> dict:
    counts = Counter()
    has_id = has_parent = has_locus_tag = 0
    total = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            counts[fields[2]] += 1
            attrs = fields[8]
            if "ID=" in attrs:
                has_id += 1
            if "Parent=" in attrs:
                has_parent += 1
            if "locus_tag=" in attrs:
                has_locus_tag += 1
            total += 1
    return {
        "path": str(path),
        "total_features": total,
        "gene": counts.get("gene", 0),
        "mRNA": counts.get("mRNA", 0),
        "CDS": counts.get("CDS", 0),
        "exon": counts.get("exon", 0),
        "with_ID": has_id,
        "with_Parent": has_parent,
        "with_locus_tag": has_locus_tag,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes-dir", type=Path,
                    default=Path("local_data/benchmark/genomes"))
    args = ap.parse_args()

    if not args.genomes_dir.exists():
        sys.exit(f"genomes dir not found: {args.genomes_dir}\n"
                 f"Run scripts/benchmark/fetch_benchmark_genomes.sh first.")

    cols = ["species", "tier", "gene", "mRNA", "CDS", "exon",
            "with_ID", "with_Parent", "with_locus_tag", "path"]
    print("\t".join(cols))

    failures = []
    for tier_dir in sorted(args.genomes_dir.iterdir()):
        if not tier_dir.is_dir():
            continue
        for species_dir in sorted(tier_dir.iterdir()):
            if not species_dir.is_dir():
                continue
            gffs = list(species_dir.glob("*.gff*"))
            if not gffs:
                failures.append(f"NO GFF: {species_dir}")
                continue
            stats = scan_gff(gffs[0])
            row = [species_dir.name, tier_dir.name,
                   str(stats["gene"]), str(stats["mRNA"]),
                   str(stats["CDS"]), str(stats["exon"]),
                   str(stats["with_ID"]), str(stats["with_Parent"]),
                   str(stats["with_locus_tag"]), stats["path"]]
            print("\t".join(row))
            if stats["gene"] == 0:
                failures.append(f"ZERO GENES: {gffs[0]}")

    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    print("\nAll GFFs OK.", file=sys.stderr)


if __name__ == "__main__":
    main()
