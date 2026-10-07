#!/usr/bin/env python3
"""Build tests/benchmark_truth/melittin_loci.tsv — the §1y COORDINATE truth set.

One row per curated MODEL (not per species, as `melittin_orthologs.tsv` is), because
"did the tool find the gene" is a question about a locus, not a species. See
docs/TODO.md §1y and docs/STATE_OF_THE_PROJECT.md D.1b.

Sources, and why each is trusted:

  * Bee rows come from `melittin_orthologs.tsv`'s `syntenic_locus` column — already
    curated, reused verbatim so the two truth files cannot drift apart.
  * Ant rows are collapsed from the curated UGENE/manual GFFs in the Ant_Venoms data
    tree. Every scaffold's assembly membership was verified against the downloaded
    FASTAs (OV788322.1 in GCA_928718305.1, CM0208xx in GCA_009859135.1,
    CM0797xx in GCA_019399895.2) rather than inferred from the file name.

Regenerate after editing either source:

    python3 scripts/benchmark/build_melittin_loci.py [--ant-dir DIR]

NOTE tests/benchmark_truth/ is gitignored (unpublished annotation data), so this script
is the committed record of how that file is produced.
"""
import argparse
import collections
import csv
import os
import re
import sys

# scripts/benchmark/<this file>  ->  repo root is three levels up
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_ANT_DIR = "/home/faw/dev/projects/Ant_Venoms/data/annotations"
BEE_SRC = os.path.join(ROOT, "tests", "benchmark_truth", "melittin_orthologs.tsv")
OUT = os.path.join(ROOT, "tests", "benchmark_truth", "melittin_loci.tsv")

# GenBank <-> RefSeq scaffold mirrors, VERIFIED BY SEQUENCE LENGTH, not assumed:
#   CM079762.1 == NC_091864.1   13,819,712 bp   Cobs3.1 LG01
#   CM079767.1 == NC_091869.1    9,247,653 bp   Cobs3.1 LG06
# The Cardiocondyla models are curated on the GenBank accessions, while a run fetches
# the RefSeq mirror -- which left 13 of 28 models unscoreable as "assembly_mismatch".
# Same assembly, same sequence, so coordinates transfer 1:1.
CHROM_ALIASES = {
    "CM079762.1": ["NC_091864.1"],
    "CM079767.1": ["NC_091869.1"],
}

# A bee row's syntenic_locus is a REGION span (Colletes' is 110 kb), not a curated gene
# model, so a hit there means "landed in the right neighbourhood" and an overlap
# FRACTION against it would be meaningless. Only Apis has a curated gene.
KIND = {"PRESENT": "locus_gene_present", "LOST": "locus_gene_lost"}

# (species, assembly, gff basename, default kind, confidence, note)
ANT_SOURCES = [
    ("Tetramorium_bicarinatum", "GCA_928718305.1",
     "toxins_in_Tbic_OV788322_annotations.gff", "gene", "gold",
     "melittin-family MYRTX (U11-MYRTX-Tb1a clade)"),
    ("Cardiocondyla_obscurior", "GCA_019399895.2",
     "toxins_in_Cobs_CM079767_annotations.gff", "gene", "gold",
     "GR2 Clade_H tandem array"),
    ("Cardiocondyla_obscurior", "GCA_019399895.2",
     "toxins_in_Cobs_CM079762_annotations.gff", "gene", "putative",
     "GR2-related tandem array, single-exon models"),
    ("Formica_selysi", "GCA_009859135.1",
     "toxins_in_Fsel_CM020806_annotations.gff", "gene", "gold", "GR2 Clade_M/N"),
    ("Formica_selysi", "GCA_009859135.1",
     "toxins_in_Fsel_CM020811_annotations.gff", "fragment", "putative",
     "GR2-related fragments"),
]


def bee_rows():
    rows = []
    with open(BEE_SRC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            loc = (r.get("syntenic_locus") or "").strip()
            if loc in ("", "—"):
                continue
            m = re.match(r"^(.+?):(\d+)-(\d+)\(([+-])\)$", loc)
            if not m:
                sys.exit(f"unparsed syntenic_locus for {r['species']}: {loc}")
            is_apis = r["species"] == "Apis_mellifera"
            rows.append(dict(
                species=r["species"], assembly=r["accession"], chrom=m.group(1),
                start=int(m.group(2)), end=int(m.group(3)), strand=m.group(4),
                model_id="gene-Melt" if is_apis else "syntenic_locus",
                truth_kind="gene" if is_apis else KIND[r["status"]],
                confidence=r["confidence"],
                source="melittin_orthologs.tsv:syntenic_locus",
                notes=("reference melittin, 2 exons" if is_apis else
                       "GR2 syntenic locus; " + ("gene present" if r["status"] == "PRESENT"
                       else "gene LOST - a call here is a FALSE POSITIVE")),
            ))
    return rows


def ant_rows(ant_dir):
    rows = []
    for species, asm, fname, kind, conf, note in ANT_SOURCES:
        path = os.path.join(ant_dir, fname)
        if not os.path.isfile(path):
            print(f"  WARNING: missing {path} — skipping", file=sys.stderr)
            continue
        spans = collections.OrderedDict()
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
                m = re.search(r"name=([^;]+)", attrs)
                name = m.group(1) if m else "unnamed"
                # A "fragment" row is a standalone single-exon model, and several can
                # share one name (CM020811 has five, 4-12 kb apart). Key those on their
                # own coordinates so they stay distinct instead of collapsing into a
                # single 21 kb span; named multi-exon models still merge across exons.
                key = (name, start) if "fragment" in name else (name,)
                if key in spans:
                    spans[key]["start"] = min(spans[key]["start"], start)
                    spans[key]["end"] = max(spans[key]["end"], end)
                else:
                    spans[key] = dict(chrom=chrom, start=start, end=end,
                                      strand=strand, name=name)
        n_frag = sum(1 for v in spans.values() if "fragment" in v["name"])
        frag_i = 0
        for v in spans.values():
            is_frag = "fragment" in v["name"]
            if is_frag:
                frag_i += 1
            rows.append(dict(
                species=species, assembly=asm, chrom=v["chrom"],
                start=v["start"], end=v["end"],
                strand=v["strand"] if v["strand"] in "+-" else ".",
                model_id=f"{v['name']}_{frag_i}" if (is_frag and n_frag > 1) else v["name"],
                # truth_kind is a property of the MODEL, not the file: a fragment row in
                # an otherwise-gold gene file is still a fragment.
                truth_kind="fragment" if is_frag else kind,
                confidence=conf, source=f"Ant_Venoms/data/annotations/{fname}", notes=note,
            ))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ant-dir", default=DEFAULT_ANT_DIR)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    rows = bee_rows() + ant_rows(args.ant_dir)
    for r in rows:
        r["chrom_aliases"] = ";".join(CHROM_ALIASES.get(r["chrom"], []))
    rows.sort(key=lambda r: (r["species"], r["chrom"], r["start"]))

    cols = ["species", "assembly", "chrom", "chrom_aliases", "start", "end", "strand",
            "model_id", "truth_kind", "confidence", "source", "notes"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kinds = collections.Counter(r["truth_kind"] for r in rows)
    print(f"wrote {len(rows)} truth models -> {args.out}")
    print("  kinds:", dict(kinds))
    print("  aliased chroms:", sorted({r["chrom"] for r in rows if r["chrom_aliases"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
