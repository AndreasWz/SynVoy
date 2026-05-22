#!/usr/bin/env python3
"""Build tests/benchmark_truth/tp53_orthologs.tsv from NCBI annotations.

Queries the NCBI `datasets` CLI for TP53 (and TP63, TP73 where applicable)
per species in the benchmark list. Emits the same column schema as
melittin_orthologs_gr2.tsv so score_benchmark.py can consume it directly.

Run:
    mamba run -n synvoy_env python scripts/benchmark/fetch_tp53_truth.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

# (snake_case, scientific_name, assembly_accession, tier, notes)
SPECIES = [
    ("Homo_sapiens",            "Homo sapiens",             "GCF_000001405.40", "home",  "Reference. TP53 / TP63 / TP73 all present."),
    ("Pan_troglodytes",         "Pan troglodytes",          "GCF_028858775.2",  "tier1", "Closest primate outgroup. Expect HIGH identity."),
    ("Macaca_mulatta",          "Macaca mulatta",           "GCF_049350105.2",  "tier1", "Old World monkey. Using T2T-MMU8v2.0 because that's where NCBI annotates TP53."),
    ("Mus_musculus",            "Mus musculus",             "GCF_000001635.27", "tier1", "Trp53 on chr11. Single MEDIUM hit in prior run."),
    ("Rattus_norvegicus",       "Rattus norvegicus",        "GCF_036323735.1",  "tier1", "Standard rodent comparator."),
    ("Heterocephalus_glaber",   "Heterocephalus glaber",    "GCF_000247695.1",  "tier1", "Naked mole rat. HetGla_female_1.0 (only annotated assembly). Single TP53, stabilised protein (longevity)."),
    ("Canis_lupus_familiaris",  "Canis lupus familiaris",   "GCF_011100685.1",  "tier1", "Dog. Laurasiatheria."),
    ("Bos_taurus",              "Bos taurus",               "GCF_002263795.3",  "tier1", "Cow. Laurasiatheria."),
    ("Loxodonta_africana",      "Loxodonta africana",       "GCF_030014295.1",  "tier1", "African elephant. KNOWN ~20 TP53 copies (Peto's paradox, Sulak 2016)."),
    ("Monodelphis_domestica",   "Monodelphis domestica",    "GCF_027887165.1",  "tier1", "Opossum. Outside placentals."),
    ("Gallus_gallus",           "Gallus gallus",            "GCF_016699485.2",  "tier2", "Chicken. In GCF_016699485.2 TP53 is on UNSCAFFOLDED contig NW_024096016.1."),
    ("Anolis_carolinensis",     "Anolis carolinensis",      "GCF_035594765.1",  "tier2", "Squamate reptile."),
    ("Crocodylus_porosus",      "Crocodylus porosus",       "GCF_001723895.1",  "tier2", "Archosaur sister to birds. TP53 NOT annotated by symbol — homology-blind test for SynVoy."),
    ("Xenopus_tropicalis",      "Xenopus tropicalis",       "GCF_000004195.4",  "tier2", "Frog. tp53 on chr3 NC_030679.2."),
    ("Latimeria_chalumnae",     "Latimeria chalumnae",      "GCF_037176945.1",  "tier3", "Coelacanth. fLatCha1.pri — newest annotated assembly. Sarcopterygian fish-tetrapod transition."),
    ("Protopterus_annectens",   "Protopterus annectens",    "GCF_019279795.1",  "tier3", "African lungfish (40 Gb genome). May fail if assembly not RefSeq-annotated."),
    ("Danio_rerio",             "Danio rerio",              "GCF_049306965.1",  "tier3", "Zebrafish. tp53 on chr5 NC_133180.1 24.74-24.75 Mb."),
    ("Oryzias_latipes",         "Oryzias latipes",          "GCF_053564925.1",  "tier3", "Medaka. tp53 on chr18 NC_136254.1 20.92-20.94 Mb."),
    ("Takifugu_rubripes",       "Takifugu rubripes",        "GCF_901000725.3",  "tier3", "Fugu. tp53 on chr8 NC_042292.1 17.32 Mb."),
    ("Callorhinchus_milii",     "Callorhinchus milii",      "GCF_018977255.1",  "tier4", "Elephant shark. TP53 NOT annotated by symbol (Pan 2011 confirmed all 3 paralogs present). Homology-blind test."),
    ("Petromyzon_marinus",      "Petromyzon marinus",       "GCF_010993605.1",  "tier4", "Sea lamprey. TP53 NOT annotated. Deep root; pre-jaw vertebrate. Expect low recall — stress test."),
]

PARALOGS = ["TP53", "TP63", "TP73"]


def datasets_summary(symbol: str, taxon: str) -> dict | None:
    """Call `datasets summary gene symbol`. Strip the deprecation banner if present."""
    try:
        out = subprocess.check_output(
            ["datasets", "summary", "gene", "symbol", symbol, "--taxon", taxon],
            stderr=subprocess.DEVNULL,
            timeout=30,
        ).decode()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    # Banner sometimes prepended to stdout. Find the first '{'.
    idx = out.find("{")
    if idx < 0:
        return None
    try:
        return json.loads(out[idx:])
    except json.JSONDecodeError:
        return None


def pick_location(payload: dict, want_assembly: str) -> dict | None:
    """Return the genomic_locations entry whose assembly_accession matches
    (ignoring patch suffix). If none match, return the first available."""
    reports = (payload or {}).get("reports") or []
    for report in reports:
        gene = report.get("gene") or {}
        for ann in gene.get("annotations") or []:
            if ann.get("assembly_accession", "").split(".")[0] == want_assembly.split(".")[0]:
                locs = ann.get("genomic_locations") or []
                if locs:
                    return {
                        "assembly": ann.get("assembly_accession"),
                        "loc": locs[0],
                        "gene_id": gene.get("gene_id"),
                        "symbol": gene.get("symbol"),
                    }
    # Fallback: first annotation with locations
    for report in reports:
        gene = report.get("gene") or {}
        for ann in gene.get("annotations") or []:
            locs = ann.get("genomic_locations") or []
            if locs:
                return {
                    "assembly": ann.get("assembly_accession"),
                    "loc": locs[0],
                    "gene_id": gene.get("gene_id"),
                    "symbol": gene.get("symbol"),
                    "fallback": True,
                }
    return None


def syntenic_locus(entry: dict) -> str:
    loc = entry["loc"]
    rng = loc.get("genomic_range") or {}
    chrom = loc.get("genomic_accession_version", "?")
    start = rng.get("begin", "?")
    end = rng.get("end", "?")
    strand_raw = rng.get("orientation", "")
    strand = "+" if strand_raw == "plus" else "-" if strand_raw == "minus" else "."
    suffix = " [diff-assembly]" if entry.get("fallback") else ""
    return f"{chrom}:{start}-{end}({strand}){suffix}"


def main() -> int:
    out_path = Path("tests/benchmark_truth/tp53_orthologs.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = ["\t".join(["species", "accession", "tier", "status", "confidence",
                       "source", "syntenic_locus", "notes"])]

    for snake, sci, assembly, tier, base_note in SPECIES:
        # Only score against TP53. TP63/TP73 looked up purely for the notes
        # field so we know whether the paralogs were annotated in this genome.
        tp53 = datasets_summary("TP53", sci)
        tp53_loc = pick_location(tp53, assembly) if tp53 else None

        paralog_notes = []
        for sym in ("TP63", "TP73"):
            data = datasets_summary(sym, sci)
            picked = pick_location(data, assembly) if data else None
            if picked:
                lc = syntenic_locus(picked)
                paralog_notes.append(f"{sym}={lc}")
            else:
                paralog_notes.append(f"{sym}=not_annotated")

        if tp53_loc:
            status = "PRESENT"
            confidence = "gold" if not tp53_loc.get("fallback") else "medium"
            locus = syntenic_locus(tp53_loc)
            source = f"NCBI_Gene_{tp53_loc.get('gene_id','?')}"
        else:
            status = "UNCERTAIN"
            confidence = "low"
            locus = "-"
            source = "NCBI_no_annotation"

        notes = base_note + " | " + "; ".join(paralog_notes)
        rows.append("\t".join([snake, assembly, tier, status, confidence,
                               source, locus, notes]))
        print(f"  {snake}: TP53={locus}", flush=True)

    out_path.write_text("\n".join(rows) + "\n")
    print(f"\nwrote {out_path} ({len(rows)-1} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
