#!/usr/bin/env python3
"""Family-recovery benchmark on the curated ant venom families.

Design: docs/NEXT_SESSION_FAMILY_BENCHMARK.md. This is its Phase 0 harness.

Why this exists
---------------
Every SynVoy claim so far rests on one gene at a time (melittin, oskar, STE2, decorin),
each with a hand-built truth set. The Ant_Venoms `final_dataset` is an independent,
labelled truth set: 939 curated gene models, each tagged `family=` and `locus=`, across
34 species. This script turns it into a benchmark that says not only how much was
found but **which pipeline layer lost what**, and which misses were never reachable by
synteny in the first place.

Subcommands
-----------
  plan    Pre-registration. For every family and every home-genome seed candidate,
          the tier of each target gene relative to that seed (below). Run it and commit
          the output BEFORE looking at any result: the tiers are the ceiling, fixed in
          advance, so a later "fix" cannot quietly redefine what counted as reachable.
  stage   Write the seed's query FASTA and a species-named target directory of
          symlinks, so SynVoy's genome names ARE the truth table's species keys.
  seedcheck  Is each curated home seed the best home-proteome homolog of the target
          genes it is meant to find? `seed_not_rbh` is a lead, not a verdict: it
          names the home gene the targets pick, and the address audit (below)
          decides which of the two sits at the label's neighbourhood.
  score   Score one SynVoy run for one family: locus-, address- and family-level
          recall, per-call precision, and an A-D layer attribution for every miss.
  callcheck  For every confident call, the home gene that best explains the called
          window (translated search): an off-truth call on the QUERY's own gene is
          an ortholog the curation lacks, one on another gene is a real false call.

Tiers -- relative to the home loci the run actually searched
------------------------------------------------------------
  home_locus       same `locus=` label as a searched home locus: the syntenic ortholog.
                   Locus-level recall is measured on these.
  home_address     a different label, but LOCUS_ANCHORS.tsv gives it the same flanking
                   address (anchor_up + anchor_dn) as a searched home locus -- e.g. the
                   copies of a tandem array. Synteny has no signal to separate
                   these, so a hit is scored `unresolvable_by_synteny`, never as an error,
                   and changing the code to "fix" one is overfitting.
  unsearched_home  the locus has a member in the home genome but the run never searched
                   that home locus (split_loci dropped it or never found it).
  no_home_anchor   no home-genome member at this locus or at its address: a lineage-
                   specific locus with no home neighbourhood to anchor on. Unreachable by
                   synteny from this home BY CONSTRUCTION; excluded from locus recall.
  unassigned_label the curation assigns the gene to NO orthologous locus (per-species
                   survey labels, per-scaffold labels): counted at family level only,
                   excluded from locus recall AND locus precision. (Amendment 1.)

Truth genes whose protein is not detectably homologous to the searched query (E-value >
MAX_SEED_EVALUE) are flagged `seed_unrelated`; locus recall is reported both raw and
excluding them. (Amendment 1 introduced the flag after a curated Vespa model scored 38
against its own family's seed; Amendment 2 made it length-aware.) Likewise genes on a
contig shorter than MIN_CONTIG_BP are flagged `assembly_limited` (Amendment 2): no search
can anchor them, so a miss there is not a layer-A search failure.

Seeds -- the home locus must BE the label's address (Amendment 2)
------------------------------------------------------------------
A seed stands for its curated `locus=` label, but the outgroup (Apis) labels were
assigned by sequence, not by position: in four families the curated Apis gene does not
sit at its label's curated address (LOCUS_ANCHORS anchor_up/anchor_dn), so its run
searches a neighbourhood the ant members do not share. Such a family is also run from
the address seed -- the Apis gene or ORF that DOES sit at the address -- via
`stage --seed <home gene | name> --family F [--label L] [--seed-faa F --seed-coords C]`,
and a curated seed that sits at ANOTHER label's address is scored under that label
(`stage --seed GM --label L`). The audit is by hand, before the runs; results in the
gitignored tests/benchmark_truth/family_benchmark_private.md.

The `no_home_anchor` tier is the big one: with Apis as home, most ant genes of the
multi-copy families (SP, PLA2, APH, VA) sit at loci Apis does not have.

Layers -- every miss in the home_locus / home_address tiers gets exactly one
----------------------------------------------------------------------------
  A_neighbourhood   no search block from a responsible home locus covers the gene
                    (block hull of that home locus's models, +/- --block-pad)
  B_modelling       a block covers it, but no GOI model overlaps it
  C_classification  a GOI model overlaps it, but below --min-confidence
  D_adjudication    a GOI model at/above the floor overlaps it, but the report demoted
                    or dropped it (paralog_not_goi / identity_coverage_decoupled / absent)

A-miss is a search bug; C is where a parameter argument is legitimate; D is the
adjudication layer. Do not tune anything that has not been attributed to a layer.

What counts as a call
---------------------
GOI models are read from the run's per-home-locus GFFs (`SynVoyRole=goi`, BOTH `mRNA`
and `gene` features -- tandem copies are emitted as `gene` only). Rescue models
(hull / strong-synteny) never reach the published GFFs, so the report's
`goi_dedup.records` are merged in. A call is a FINAL confident call when the report
lists it HIGH/MEDIUM and did not relabel it `paralog_not_goi` or
`identity_coverage_decoupled` -- the same rule generate_report.py uses for its headline.

Usage
-----
    family_recovery.py plan  [--data DIR] [--out plan.tsv]
    family_recovery.py seedcheck [--out check.tsv]
    family_recovery.py stage --seed <HOME_GM_ID | HOME_GENE --family F> [--label L]
                             [--seed-faa ORF.faa --seed-coords C] --out RUN_INPUTS/
    family_recovery.py score --run RESULTS/ --family FAM_APYR [--out-prefix P]
    family_recovery.py callcheck --run RESULTS/ --targets RUN_INPUTS/targets

`--data` is the Ant_Venoms `data/` directory holding `final_dataset/` and `genomes/`
(default: $SYNVOY_FAMILY_DATA, else ../Ant_Venoms_git/data next to this repo).
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from score_coordinates import _attrs, distance_bp, overlap_bp  # noqa: E402

DEFAULT_DATA = (os.environ.get("SYNVOY_FAMILY_DATA")
                or os.path.normpath(os.path.join(ROOT, "..", "Ant_Venoms_git", "data")))
DEFAULT_HOME = "Apis_mellifera"

# GR1/GR2/GR3/GR7 are the paper's showcase regions, benchmarked separately; VA2 (7 genes,
# 2 species) and VA4 (2 genes, 1 species) are too small to say anything.
DEFAULT_EXCLUDE = ("FAM_GR1", "FAM_GR2", "FAM_GR3", "FAM_GR7", "FAM_VA2", "FAM_VA4")

# Pre-registered difficulty ladder and tune / hold-out split (design doc §2.2, §4).
# Changing a family's split AFTER seeing results defeats the purpose of having one.
LADDER = {
    "FAM_APYR": 1, "FAM_SPIN": 1, "FAM_DPP4": 1,
    "FAM_HYAL": 2, "FAM_KAZA": 2,
    "FAM_SECP": 3,
    "FAM_VA": 4, "FAM_SP": 4, "FAM_PLA2": 4,
    "FAM_APH": 5,
}
SPLIT = {
    "FAM_SECP": "tune", "FAM_PLA2": "tune",
    "FAM_VA": "holdout", "FAM_SP": "holdout", "FAM_APH": "holdout",
}

CONF_RANK = {"LOW": 0, "AMBIGUOUS": 1, "MEDIUM": 2, "HIGH": 3}
DEMOTED_CLASSES = {"paralog_not_goi", "identity_coverage_decoupled"}
REACHABLE_TIERS = ("home_locus", "home_address")
# Labels that assign a gene to NO orthologous locus: per-species survey labels (added
# with one late species) and per-scaffold labels. Neither credited nor penalised at locus
# level -- the truth set does not say which locus they are. Added 2026-09-10 after the
# first run, as a documented deviation from the hashed pre-registration (design doc §0).
UNASSIGNED_LABEL_RE = re.compile(r"^(?:TETLOC_|scaf:)|unresolved")
# A truth protein whose best local alignment to the searched query has an E-value above
# this is not detectably homologous to it: no sequence- or synteny-anchored search from
# this seed can be expected to call it, and it may be a curation error. Length-aware on
# purpose (Amendment 2): the Amendment-1 raw-score floor of 100 flagged EVERY member of a
# 77-aa family whose true orthologs score 53-99. BLASTP defaults -- BLOSUM62, gap 11/1,
# which is parasail open 12 / extend 1 (parasail charges `open` for the first gap
# residue) -- with the matching Karlin-Altschul lambda and K.
MAX_SEED_EVALUE = 1e-3
SEED_SW_GAP = (12, 1)
KA_LAMBDA, KA_K = 0.267, 0.041
# A truth gene on a contig shorter than this cannot be reached by synteny whatever the
# search does: SynVoy's block filter needs >= 2 flanking-gene hits on the same sequence
# (min_block_genes), and at insect gene density a contig this short holds the target and
# at most one neighbour. Flagged `assembly_limited` and reported both ways, like
# seed_unrelated (Amendment 2: a KAZA gene on a 3.5-kb contig read as a search miss).
MIN_CONTIG_BP = 20000
LAYERS = ("A_neighbourhood", "B_modelling", "C_classification", "D_adjudication")
LOCUS_RE = re.compile(r"(locus_\d+)")
BLOCK_RE = re.compile(r"_b(\d+)_")


# ---------------------------------------------------------------- truth ---
def species_key(name):
    """'Probolomyrmex sp. GAGA-0580' -> 'Probolomyrmex_sp_GAGA-0580'.

    The manifest spells species with spaces and punctuation; ACCESSIONS.tsv, the
    per-species protein files and SynVoy's genome names use this underscore form.
    """
    return re.sub(r"[^A-Za-z0-9_-]", "", name.strip().replace(" ", "_"))


def load_manifest(path, exclude_families=DEFAULT_EXCLUDE):
    """One dict per curated gene model from EXPORT_MANIFEST.tsv.

    The manifest was verified identical to gff/ALL_SPECIES.gff3 (scaffold, start, end,
    strand, family, locus) for all 939 genes on 2026-09-10; it is used because it is
    one row per gene with no parsing ambiguity.
    """
    genes = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"] in exclude_families:
                continue
            genes.append({
                "gm": r["gm"], "species": species_key(r["species"]),
                "family": r["family"], "locus": r["locus"], "chrom": r["scaffold"],
                "start": int(r["start"]), "end": int(r["end"]), "strand": r["strand"],
                "status": r.get("structural_status", ""),
                "prot_len": int(r.get("prot_len") or 0),
            })
    if not genes:
        sys.exit(f"no truth genes in {path}")
    return genes


def load_truth_cds(data_dir):
    """{gm: [(start, end), ...]} of the curated CDS intervals, from the per-species GFF3s.

    A call is matched to a truth gene by overlap. Span overlap alone is too weak where a
    curated gene is intron-rich: in the 2026-09-15 boundary audit one VA call sat on the
    OPPOSITE strand inside a 10.2-kb first intron, sharing no coding sequence with the gene
    it was credited to (1 of 228 pairings). Returns {} when the GFF3s are unavailable, and
    `_ov` then falls back to span overlap.
    """
    out = {}
    pattern = os.path.join(data_dir, "final_dataset", "browser_gff", "*.gff3")
    for path in sorted(glob.glob(pattern)):
        if os.path.basename(path).startswith("ALL_SPECIES"):
            continue
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                m = re.search(r"GM_\d+", f[8])
                if m:
                    out.setdefault(m.group(0), []).append((int(f[3]), int(f[4])))
    return out


def load_addresses(path):
    """locus -> (family, anchor_up, anchor_dn), or None when either side is uncurated.

    Two loci are synteny-equivalent only when BOTH flanks match: a shared one-sided
    anchor is not an address. Keyed within family, so an SP locus that happens to sit
    in the HYAL neighbourhood is not made equivalent to HYAL.
    """
    addr = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            up, dn = (r.get("anchor_up") or "").strip(), (r.get("anchor_dn") or "").strip()
            addr[r["stable_id"]] = (r["family_stable_id"], up, dn) if up and dn else None
    return addr


def assign_tiers(genes, family, home_species, searched_labels, addresses):
    """Tier every non-home gene of `family` against the home loci that were searched."""
    fam = [g for g in genes if g["family"] == family]
    home_labels = {g["locus"] for g in fam if g["species"] == home_species}
    searched_addr = {addresses.get(lb) for lb in searched_labels} - {None}
    home_addr = {addresses.get(lb) for lb in home_labels} - {None}
    tiered = []
    for g in fam:
        if g["species"] == home_species:
            continue
        a = addresses.get(g["locus"])
        if UNASSIGNED_LABEL_RE.search(g["locus"]):
            tier = "unassigned_label"
        elif g["locus"] in searched_labels:
            tier = "home_locus"
        elif a is not None and a in searched_addr:
            tier = "home_address"
        elif g["locus"] in home_labels or (a is not None and a in home_addr):
            tier = "unsearched_home"
        else:
            tier = "no_home_anchor"
        tiered.append(dict(g, tier=tier))
    return tiered


def label_is_spread(genes, species, locus, max_span=1_000_000):
    """True when one species' genes carrying `locus` sit on >1 scaffold or span >1 Mb.

    Then the label is an orthogroup, not a genomic address, and locus-level scoring
    against it is weaker. In the curated set, one family's Apis label covers three
    genes on two chromosomes.
    """
    rows = [g for g in genes if g["species"] == species and g["locus"] == locus]
    if len({g["chrom"] for g in rows}) > 1:
        return True
    return bool(rows) and (max(g["end"] for g in rows) - min(g["start"] for g in rows)
                           > max_span)


# ----------------------------------------------------------------- plan ---
def build_plan(genes, addresses, home_species=DEFAULT_HOME):
    """One row per (family, home-genome seed candidate): the pre-registered ceiling."""
    rows = []
    for family in sorted({g["family"] for g in genes}):
        fam = [g for g in genes if g["family"] == family]
        target_species = {g["species"] for g in fam if g["species"] != home_species}
        seeds = sorted((g for g in fam if g["species"] == home_species),
                       key=lambda g: (g["locus"], g["gm"]))
        if not seeds:
            rows.append({"family": family, "seed_gm": "", "warnings": "no_home_member"})
            continue
        for s in seeds:
            tiers = assign_tiers(genes, family, home_species, {s["locus"]}, addresses)
            c = Counter(t["tier"] for t in tiers)
            reach_sp = {t["species"] for t in tiers if t["tier"] in REACHABLE_TIERS}
            warn = []
            if label_is_spread(genes, home_species, s["locus"]):
                warn.append("home_label_not_an_address")
            if c["home_address"]:
                warn.append("shared_address_unresolvable")
            rows.append({
                "family": family, "ladder": LADDER.get(family, ""),
                "split": SPLIT.get(family, "ladder"),
                "seed_gm": s["gm"], "seed_locus": s["locus"],
                "seed_coords": f"{s['chrom']}:{s['start']}-{s['end']}{s['strand']}",
                "seed_prot_len": s["prot_len"], "target_genes": len(tiers),
                "home_locus": c["home_locus"], "home_address": c["home_address"],
                "unsearched_home": c["unsearched_home"],
                "no_home_anchor": c["no_home_anchor"],
                "unassigned_label": c["unassigned_label"],
                "reachable_species": len(reach_sp),
                "target_species_with_family": len(target_species),
                "warnings": ",".join(warn),
            })
    return rows


# ---------------------------------------------------------------- stage ---
def find_genome_files(genomes_root, accession):
    """(fasta, gff or None) for one assembly directory under data/genomes/."""
    d = os.path.join(genomes_root, accession)
    fastas = sorted(p for p in glob.glob(os.path.join(d, "*"))
                    if p.endswith((".fna", ".fa", ".fasta")))
    gffs = sorted(glob.glob(os.path.join(d, "*.gff")) + glob.glob(os.path.join(d, "*.gff3")))
    return (fastas[0] if fastas else None), (gffs[0] if gffs else None)


def read_fasta_record(path, record_id):
    seq, name, grab = [], None, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if grab:
                    break
                name = line[1:].split()[0]
                grab = name == record_id
                header = line.rstrip("\n")
                continue
            if grab:
                seq.append(line.strip())
    if not seq:
        sys.exit(f"{record_id} not found in {path}")
    return header, "".join(seq)


def home_gene_record(gff_path, gene):
    """(chrom, start, end, strand, [protein ids]) of one home-genome gene, by its
    `gene=` / `Name=` attribute in a RefSeq GFF; None when the GFF has no such gene."""
    span, pids = None, []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "CDS"):
                continue
            a = _attrs(f[8])
            if gene not in (a.get("gene"), a.get("Name")):
                continue
            if f[2] == "gene" and span is None:
                span = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == "CDS" and a.get("protein_id") and a["protein_id"] not in pids:
                pids.append(a["protein_id"])
    return (*span, pids) if span else None


def seed_label(genes, family, label, home_species=DEFAULT_HOME, what="seed"):
    """The curated label a non-manifest seed stands for: `label` (any curated label of
    the family -- the address audit decides which), else the family's only home label."""
    fam_labels = {g["locus"] for g in genes if g["family"] == family}
    if label is None:
        home = sorted({g["locus"] for g in genes
                       if g["family"] == family and g["species"] == home_species})
        if len(home) != 1:
            sys.exit(f"{family} has {len(home)} home labels: pass --label for {what}")
        return home[0]
    if label not in fam_labels:
        sys.exit(f"--label {label} is not a curated label of {family}")
    return label


def parse_coords(text):
    """'CHROM:START-END[:STRAND]' -> (chrom, start, end, strand), 1-based inclusive."""
    m = re.match(r"^(.+):(\d+)-(\d+)(?::([+-]))?$", text or "")
    if not m:
        sys.exit(f"--seed-coords {text!r}: expected CHROM:START-END[:STRAND]")
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4) or "."


def home_gene_seed(data, gene, family, genes, home_species=DEFAULT_HOME, label=None):
    """The address seed (Amendment 2): a home-genome gene that is NOT a curated member.

    The query is that gene's longest RefSeq isoform. The run's home locus has no curated
    label, so the seed takes one (see seed_label). Returns (seed_info, header, seq).
    """
    accessions = {}
    with open(os.path.join(data, "genomes", "ACCESSIONS.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accessions[species_key(r["species"])] = r["accession"]
    hdir = os.path.join(data, "genomes", accessions[home_species])
    home_fa, home_gff = find_genome_files(os.path.join(data, "genomes"),
                                          accessions[home_species])
    rec = home_gene_record(home_gff, gene)
    if rec is None:
        sys.exit(f"seed {gene}: neither a manifest gene model nor a gene in {home_gff}")
    chrom, start, end, strand, pids = rec
    proteome = read_fasta(os.path.join(hdir, "protein.faa"))
    isoforms = sorted((p for p in pids if p in proteome),
                      key=lambda p: (-len(proteome[p]), p))
    if not isoforms:
        sys.exit(f"seed {gene}: none of its proteins {pids} is in {hdir}/protein.faa")
    label = seed_label(genes, family, label, home_species, what=gene)
    pid = isoforms[0]
    info = {"kind": "home_gene", "seed": gene, "protein_id": pid, "family": family,
            "locus": label, "chrom": chrom, "start": start, "end": end, "strand": strand,
            "prot_len": len(proteome[pid])}
    return info, f">{pid} {gene} home-gene seed for {family}", proteome[pid]


def stage(data, seed_gm, out, home_species=DEFAULT_HOME, exclude=DEFAULT_EXCLUDE,
          family=None, label=None, seed_faa=None, seed_coords=None):
    """Write query.faa, seed.json, home/, targets/<Species>.fna symlinks.

    `seed_gm` is a curated home gene model -- with `label`, scored under that label
    instead of its manifest one -- or, with `family`, a home-genome gene id (see
    home_gene_seed), or, with `seed_faa` + `seed_coords`, the name of an unannotated
    home ORF (Amendment 2).
    """
    fd = os.path.join(data, "final_dataset")
    genes = load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"), exclude_families=())
    seed = next((g for g in genes if g["gm"] == seed_gm), None)
    if seed_faa:
        if seed is not None or family is None or not seed_coords:
            sys.exit("an ORF seed needs a new name, --family and --seed-coords")
        recs = read_fasta(seed_faa)
        if len(recs) != 1:
            sys.exit(f"{seed_faa}: expected exactly one protein, found {len(recs)}")
        (pid, seq), = recs.items()
        chrom, start, end, strand = parse_coords(seed_coords)
        info = {"kind": "orf", "seed": seed_gm, "family": family,
                "locus": seed_label(genes, family, label, home_species, what=seed_gm),
                "chrom": chrom, "start": start, "end": end, "strand": strand,
                "prot_len": len(seq), "source_faa": os.path.basename(seed_faa)}
        header = f">{seed_gm} home ORF seed for {family}"
        seed = {"gm": seed_gm, "species": home_species, **info}
    elif seed is None:
        if family is None:
            sys.exit(f"seed {seed_gm} not in the manifest (a home-gene seed needs --family)")
        info, header, seq = home_gene_seed(data, seed_gm, family, genes, home_species,
                                           label)
        seed = {"gm": seed_gm, "species": home_species, **info}
    else:
        if seed["species"] != home_species:
            sys.exit(f"seed {seed_gm} is from {seed['species']}, not the home genome "
                     f"{home_species}: SynVoy locates the query in the home genome")
        if family is not None and seed["family"] != family:
            sys.exit(f"seed {seed_gm} is a {seed['family']} gene, not {family}")
        info = {"kind": "curated", "seed": seed_gm, "family": seed["family"],
                "locus": seed["locus"], "chrom": seed["chrom"], "start": seed["start"],
                "end": seed["end"], "strand": seed["strand"], "prot_len": seed["prot_len"]}
        if label is not None and label != seed["locus"]:
            # The address audit put this gene at ANOTHER label's neighbourhood.
            info.update(locus=seed_label(genes, seed["family"], label, home_species,
                                         what=seed_gm), manifest_locus=seed["locus"])
            seed = dict(seed, locus=info["locus"])
        header, seq = read_fasta_record(
            os.path.join(fd, "proteins", f"{home_species}.faa"), seed_gm)

    accessions = {}
    with open(os.path.join(data, "genomes", "ACCESSIONS.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accessions[species_key(r["species"])] = r["accession"]

    os.makedirs(os.path.join(out, "targets"), exist_ok=True)
    os.makedirs(os.path.join(out, "home"), exist_ok=True)
    with open(os.path.join(out, "query.faa"), "w") as fh:
        fh.write(header + "\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    with open(os.path.join(out, "seed.json"), "w") as fh:
        json.dump(info, fh, indent=2, sort_keys=True)

    home_fa, home_gff = find_genome_files(os.path.join(data, "genomes"),
                                          accessions[home_species])
    if not home_fa or not home_gff:
        sys.exit(f"home genome {accessions[home_species]} needs both a FASTA and a GFF")
    _relink(home_fa, os.path.join(out, "home", f"{home_species}.fna"))
    _relink(home_gff, os.path.join(out, "home", f"{home_species}.gff"))

    # Every curated species is a target for EVERY family -- same targets, only the query
    # changes. A species lacking the family is the absence test, not a skipped cell.
    staged, missing = [], []
    for sp in sorted({g["species"] for g in genes} - {home_species}):
        fa, _gff = find_genome_files(os.path.join(data, "genomes"), accessions.get(sp, ""))
        if not fa:
            missing.append(sp)
            continue
        _relink(fa, os.path.join(out, "targets", f"{sp}.fna"))
        staged.append(sp)
    return {"seed": seed, "staged": staged, "missing": missing,
            "home_fasta": home_fa, "home_gff": home_gff}


def _relink(src, dst):
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.abspath(src), dst)


# ------------------------------------------------------------ run loading ---
def genome_of(path):
    """Species key from a SynVoy per-genome file name (drops locus prefix + suffixes)."""
    name = re.sub(r"^locus_\d+__", "", os.path.basename(path))
    prev = None
    while prev != name:
        prev = name
        name = re.sub(r"\.(?:gff3?|fna|fa|fasta|hull_rescue|rescue)$", "", name)
    return name


def home_locus_of(path):
    """`locus_N` of the home locus a GFF came from, or None when the layout lacks one."""
    m = re.match(r"^(locus_\d+)__", os.path.basename(path))
    if m:
        return m.group(1)
    found = LOCUS_RE.findall(os.path.dirname(path))
    return found[-1] if found else None


def load_home_loci(run_dir):
    """{locus_id: (chrom, start, end)} of every home locus, 1-based inclusive.

    Read from the query's home-genome hits (`locus_N.bed`), published under
    plot_inputs_*/ and, with --keep_intermediate, intermediate/split_loci/.
    """
    loci = {}
    for pat in (os.path.join(run_dir, "plot_inputs_*", "locus_*.bed"),
                os.path.join(run_dir, "intermediate", "split_loci", "locus_*.bed")):
        for p in sorted(glob.glob(pat)):
            m = re.match(r"^(locus_\d+)\.bed$", os.path.basename(p))
            if not m or m.group(1) in loci:
                continue
            rows = []
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 3 and not line.startswith(("#", "track")):
                        rows.append((f[0], int(f[1]) + 1, int(f[2])))
            if rows:
                chrom = rows[0][0]
                same = [r for r in rows if r[0] == chrom]
                loci[m.group(1)] = (chrom, min(r[1] for r in same), max(r[2] for r in same))
    return loci


def label_home_loci(home_loci, genes, family, home_species, slack=5000):
    """{locus_id: curated locus label, or None} for each home locus of the run.

    None means SynVoy searched a home region that is not a curated member of the
    family -- a home-locus selection error in its own right, reported as such.
    """
    home = [g for g in genes if g["family"] == family and g["species"] == home_species]
    labels = {}
    for lid, (chrom, s, e) in sorted(home_loci.items()):
        best = None
        for g in home:
            if g["chrom"] != chrom:
                continue
            ov = overlap_bp(s - slack, e + slack, g["start"], g["end"])
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        labels[lid] = best[1]["locus"] if best else None
    return labels


def gather_run_gffs(run_dir):
    """Every region GFF under a run. Files only: `logs/` holds per-task DIRECTORIES
    named after their output, e.g. `..._RESCUE_GOI_HULL__..._X.hull_rescue.gff/`."""
    if os.path.isfile(run_dir):
        return [run_dir]
    return sorted(p for p in glob.glob(os.path.join(run_dir, "**", "*.gff*"), recursive=True)
                  if os.path.isfile(p) and os.path.basename(p) != "home_genome.gff")


def load_run_models(paths, home_species=DEFAULT_HOME):
    """GOI calls and per-block model spans from SynVoy region GFFs.

    Returns (goi_calls, block_spans). `block_spans` maps
    (home_locus, species, chrom, block_idx) -> [start, end] over EVERY model carrying
    that block's `_b<N>_` id -- the footprint of one search block, used for layer A.
    """
    calls, seen, spans = [], set(), {}
    for p in paths:
        sp = genome_of(p)
        if sp in (home_species, "home_genome"):
            continue
        hl = home_locus_of(p)
        with open(p) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh
                    if not line.startswith("#") and line.strip()]
        # Tandem copies are `gene`-only GOI models; a gene that parents an mRNA is not.
        parents = {_attrs(f[8]).get("Parent") for f in rows
                   if len(f) >= 9 and f[2] == "mRNA"}
        for f in rows:
            if len(f) < 9 or f[2] not in ("mRNA", "gene"):
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            a = _attrs(f[8])
            mid = a.get("ID", "")
            m = BLOCK_RE.search(mid.split("|", 1)[-1])
            if m:
                span = spans.setdefault((hl, sp, f[0], int(m.group(1))), [start, end])
                span[0], span[1] = min(span[0], start), max(span[1], end)
            if a.get("SynVoyRole") != "goi" or (f[2] == "gene" and mid in parents):
                continue
            key = (hl, sp, f[0], start, end, mid)
            if key in seen:
                continue
            seen.add(key)
            calls.append({
                "species": sp, "chrom": f[0], "start": start, "end": end,
                "strand": f[6], "confidence": a.get("Confidence", "").upper(),
                "goi_class": a.get("GOIClass", ""), "identity": a.get("Identity", f[5]),
                "query_cov": a.get("QueryCoverage", ""), "id": mid,
                "evidence": a.get("EvidenceType", f[1]), "home_loci": {hl} - {None},
                "source": "gff",
            })
    return calls, spans


def load_report_records(run_dir):
    """goi_dedup.records from synvoy_report.json -- the ADJUDICATED calls.

    These carry the post-ownership / post-coverage-demotion labels and include rescue
    models that never reach the published GFFs. Empty when the report is absent.
    """
    path = os.path.join(run_dir, "synvoy_report.json") if os.path.isdir(run_dir) else ""
    if not path or not os.path.isfile(path):
        return []
    with open(path) as fh:
        rep = json.load(fh)
    out = []
    for r in (rep.get("goi_dedup") or {}).get("records", []) or []:
        try:
            start, end = int(r["start"]), int(r["end"])
        except (KeyError, TypeError, ValueError):
            continue
        conf = str(r.get("confidence", "")).upper()
        cls = r.get("goi_class", "")
        out.append({
            "species": species_key(str(r.get("genome", ""))), "chrom": r.get("chrom", ""),
            "start": start, "end": end, "confidence": conf, "goi_class": cls,
            "identity": r.get("identity", ""), "query_cov": r.get("query_cov", ""),
            "id": r.get("mrna_id", ""),
            # The report's FINAL attribution: the §1m owning locus when ownership ran
            # (it may reassign a call to the paralogous home locus it really belongs
            # to), else the home loci whose search produced it.
            "home_loci": ({r["owning_locus"]} if r.get("owning_locus")
                          else set(r.get("provenance") or [])),
            "final": conf in ("HIGH", "MEDIUM") and cls not in DEMOTED_CLASSES,
            "phylo_verdict": r.get("phylo_verdict", ""),
            "source": "report",
        })
    return out


def load_seed_info(run_dir):
    """seed.json written by `stage` and copied into the run by run_family.sh, or None."""
    p = os.path.join(run_dir, "seed.json") if os.path.isdir(run_dir) else ""
    if not p or not os.path.isfile(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def seed_home_member(info, home_species=DEFAULT_HOME):
    """A home-gene or ORF seed as a pseudo truth gene, so the run's home locus -- which
    no curated gene overlaps -- is labelled with the curated label it was assigned. None
    for a curated seed, which is in the manifest already."""
    if not info or info.get("kind") not in ("home_gene", "orf"):
        return None
    return {"gm": f"SEED:{info['seed']}", "species": home_species,
            "family": info["family"], "locus": info["locus"], "chrom": info["chrom"],
            "start": int(info["start"]), "end": int(info["end"]),
            "strand": info.get("strand", "."), "status": f"{info['kind']}_seed",
            "prot_len": int(info.get("prot_len") or 0)}


def apply_seed(genes, info, home_species=DEFAULT_HOME):
    """The truth genes as this run's seed defines them (Amendment 2).

    A home-gene / ORF seed joins as a pseudo home member under its label; a curated seed
    relabelled by the address audit (`manifest_locus` recorded) carries its new label.
    Returns a new list; `genes` is untouched."""
    member = seed_home_member(info, home_species)
    if member is not None:
        return genes + [member]
    if info and info.get("kind") == "curated" and info.get("manifest_locus"):
        return [dict(g, locus=info["locus"]) if g["gm"] == info["seed"] else g
                for g in genes]
    return genes


def load_searched_species(run_dir):
    """Target species of a run, from regions/<species>.<ext>.regions.bed (every target
    gets one whether or not anything was found). None when the layout is absent."""
    d = os.path.join(run_dir, "regions") if os.path.isdir(run_dir) else ""
    if not d or not os.path.isdir(d):
        return None
    sp = {genome_of(re.sub(r"\.(?:regions\.bed|scores\.tsv)$", "", n))
          for n in os.listdir(d) if n.endswith((".regions.bed", ".scores.tsv"))}
    return sp or None


def load_contig_lengths(data, species, home_species=DEFAULT_HOME):
    """{(species, contig): length} from the `.fai` next to each target genome FASTA.
    Species without an index are simply absent (nothing is flagged for them)."""
    accessions = {}
    with open(os.path.join(data, "genomes", "ACCESSIONS.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accessions[species_key(r["species"])] = r["accession"]
    out = {}
    for sp in sorted(set(species) - {home_species}):
        fa, _gff = find_genome_files(os.path.join(data, "genomes"), accessions.get(sp, ""))
        if not fa or not os.path.isfile(fa + ".fai"):
            continue
        with open(fa + ".fai") as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) >= 2:
                    out[(sp, f[0])] = int(f[1])
    return out


def read_fasta(path):
    """{first header token: sequence}."""
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name is not None:
                out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def load_run_query(run_dir):
    """The protein the run actually searched with (normalized query), or None."""
    p = os.path.join(run_dir, "intermediate", "query", "normalized_query.faa")
    if not os.path.isfile(p):
        return None
    seqs = read_fasta(p)
    return next(iter(seqs.values()), None)


def sw_evalue(score, m, n):
    """Karlin-Altschul E-value of a raw local-alignment score for an m x n comparison."""
    bits = (KA_LAMBDA * score - math.log(KA_K)) / math.log(2)
    return m * n * 2.0 ** -bits


def seed_similarity(seed_seq, proteins, gms):
    """{gm: E-value of the truth protein's best local alignment to the seed}.

    Returns {} when parasail is unavailable -- then nothing is flagged, rather than
    everything. See MAX_SEED_EVALUE.
    """
    try:
        import parasail
    except ImportError:
        print("  ** parasail not importable: seed-homology flags skipped")
        return {}
    return {gm: sw_evalue(parasail.sw_striped_16(seed_seq, proteins[gm], *SEED_SW_GAP,
                                                 parasail.blosum62).score,
                          len(seed_seq), len(proteins[gm]))
            for gm in gms if proteins.get(gm)}


# ---------------------------------------------------------------- score ---
def _ov(a, b):
    """Does call `b` land on truth gene `a`?

    Coding overlap when the gene's curated CDS intervals are known (`a["cds"]`), span
    overlap otherwise — an intron-rich gene's span can be overlapped by a model that shares
    no coding sequence with it (see load_truth_cds).
    """
    if a["chrom"] != b["chrom"]:
        return False
    if not overlap_bp(a["start"], a["end"], b["start"], b["end"]) > 0:
        return False
    cds = a.get("cds")
    if not cds:
        return True
    return any(overlap_bp(cs, ce, b["start"], b["end"]) > 0 for cs, ce in cds)


def score_run(genes, family, home_species, home_labels, addresses, gff_calls,
              block_spans, report_calls, searched_species=None, min_confidence="MEDIUM",
              block_pad=50000, seed_evalues=None, max_seed_evalue=MAX_SEED_EVALUE,
              contig_lengths=None, min_contig_bp=MIN_CONTIG_BP):
    """Score one family run. Returns (gene_rows, call_rows, summary).

    `seed_evalues` maps gm -> E-value of the truth protein against the searched query
    (see seed_similarity); genes above `max_seed_evalue` are flagged `seed_unrelated`.
    Without it nothing is flagged. `contig_lengths` ({(species, contig): bp}, see
    load_contig_lengths) flags genes on contigs below `min_contig_bp` `assembly_limited`.
    """
    floor = CONF_RANK[min_confidence]
    seed_evalues = seed_evalues or {}
    contig_lengths = contig_lengths or {}
    searched_labels = {lb for lb in home_labels.values() if lb}
    tiered = assign_tiers(genes, family, home_species, searched_labels, addresses)
    fam_by_species = defaultdict(list)
    for g in tiered:
        fam_by_species[g["species"]].append(g)

    def label_address(lb):
        return addresses.get(lb) if lb else None

    def responsible(gene):
        """Home loci whose label is this gene's locus or shares its address."""
        a = addresses.get(gene["locus"])
        return {lid for lid, lb in home_labels.items()
                if lb and (lb == gene["locus"] or (a is not None and label_address(lb) == a))}

    # A record without provenance can only have come from the one home locus of a
    # single-locus run; with several, leave it unattributed rather than guess.
    if len(home_labels) == 1:
        only = set(home_labels)
        gff_calls, report_calls = (
            [dict(c, home_loci=c["home_loci"] or set(only)) for c in cs]
            for cs in (gff_calls, report_calls))
    final_calls = [c for c in report_calls if c["final"]]
    gene_rows = []
    for g in tiered:
        row = {k: g[k] for k in ("gm", "species", "locus", "chrom", "start", "end", "tier")}
        ev = seed_evalues.get(g["gm"])
        row["seed_evalue"] = "" if ev is None else f"{ev:.2g}"
        row["seed_unrelated"] = "yes" if (ev is not None and ev > max_seed_evalue) else ""
        clen = contig_lengths.get((g["species"], g["chrom"]))
        row["contig_bp"] = "" if clen is None else clen
        row["assembly_limited"] = "yes" if (clen is not None and clen < min_contig_bp) else ""
        if searched_species is not None and g["species"] not in searched_species:
            gene_rows.append(dict(row, verdict="not_searched", layer="", best_call="",
                                  nearest_block_bp=""))
            continue
        resp = responsible(g)
        hits = [c for c in final_calls if c["species"] == g["species"] and _ov(g, c)]
        locus_hit = any(lb == g["locus"] for c in hits for lid in c["home_loci"]
                        for lb in [home_labels.get(lid)])
        addr_hit = any(lid in resp for c in hits for lid in c["home_loci"])
        if g["tier"] == "home_locus":
            if locus_hit:
                verdict = "recovered"
            elif addr_hit:
                verdict = "unresolvable_by_synteny"
            elif hits:
                verdict = "family_only"
            else:
                verdict = "missed"
        elif g["tier"] == "home_address":
            verdict = ("unresolvable_by_synteny" if (addr_hit or locus_hit)
                       else "family_only" if hits else "missed")
        elif g["tier"] == "unassigned_label":
            verdict = "unassigned_found" if hits else "unassigned_not_found"
        else:
            verdict = "found_without_anchor" if hits else "unreachable"

        layer, nearest = "", ""
        if verdict == "missed":
            layer, nearest = attribute_layer(g, gff_calls, report_calls, block_spans,
                                             floor, block_pad)
        best = max(hits, key=lambda c: CONF_RANK.get(c["confidence"], -1), default=None)
        gene_rows.append(dict(
            row, verdict=verdict, layer=layer, nearest_block_bp=nearest,
            best_call=(f"{best['chrom']}:{best['start']}-{best['end']} {best['confidence']}"
                       if best else "")))

    def classify(c):
        on = [g for g in fam_by_species.get(c["species"], []) if _ov(g, c)]
        lbs = {home_labels.get(lid) for lid in c["home_loci"]} - {None}
        addrs = {label_address(lb) for lb in lbs} - {None}
        other = sorted({x["family"] for x in genes
                        if x["family"] != family and x["species"] == c["species"]
                        and _ov(x, c)})
        if any(g["locus"] in lbs for g in on):
            cls = "on_locus"
        elif any(addresses.get(g["locus"]) in addrs for g in on):
            cls = "on_address"
        elif on and all(g["tier"] == "unassigned_label" for g in on):
            cls = "on_unassigned_label"      # family-correct; no locus to judge it by
        elif on:
            cls = "on_other_family_locus"
        elif other:
            cls = "other_family"
        else:
            cls = "off_truth"
        return cls, on, other

    # Every adjudicated record, whatever its tier: what the AMBIGUOUS tier and the F9 phylo
    # verdict are worth is measured against the truth here (design doc §5, deliverable 3).
    census = {"by_confidence": defaultdict(Counter), "by_phylo_verdict": defaultdict(Counter)}
    for c in report_calls:
        if searched_species is not None and c["species"] not in searched_species:
            continue
        cls = classify(c)[0]
        tier = c["confidence"] + ("" if c["final"] or c["confidence"] not in ("HIGH", "MEDIUM")
                                  else f"->{c['goi_class']}")
        census["by_confidence"][tier][cls] += 1
        census["by_phylo_verdict"][c.get("phylo_verdict") or "none"][cls] += 1

    call_rows = []
    for c in final_calls:
        if searched_species is not None and c["species"] not in searched_species:
            continue
        cls, on, other = classify(c)
        call_rows.append({
            "species": c["species"], "chrom": c["chrom"], "start": c["start"],
            "end": c["end"], "confidence": c["confidence"], "goi_class": c["goi_class"],
            "identity": c["identity"], "home_loci": ",".join(sorted(c["home_loci"])),
            "call_class": cls, "truth_loci": ",".join(sorted({g["locus"] for g in on})),
            "other_families": ",".join(other),
        })

    summary = summarize(gene_rows, call_rows, home_labels, fam_by_species, final_calls,
                        searched_species)
    summary["record_census"] = {k: {t: dict(c) for t, c in sorted(v.items())}
                                for k, v in census.items()}
    return gene_rows, call_rows, summary


def attribute_layer(gene, gff_calls, report_calls, block_spans, floor, block_pad):
    """Which layer lost a reachable gene. Returns (layer, nearest_block_distance_bp).

    Layer A asks whether ANY search block reached the gene, whichever home locus seeded
    it: a neighbourhood searched from a paralogous home locus that still produced no
    model is a modelling failure, not a search one. A block only leaves a trace when it
    produced at least one model, so A strictly means "no model-bearing block".
    """
    spans = [(v[0], v[1]) for k, v in block_spans.items()
             if k[1] == gene["species"] and k[2] == gene["chrom"]]
    nearest = min((distance_bp(gene["start"], gene["end"], s, e) for s, e in spans),
                  default="")
    if nearest == "" or nearest > block_pad:
        return "A_neighbourhood", nearest
    models = [c for c in gff_calls + report_calls
              if c["species"] == gene["species"] and _ov(gene, c)]
    if not models:
        return "B_modelling", nearest
    if max(CONF_RANK.get(c["confidence"], -1) for c in models) < floor:
        return "C_classification", nearest
    return "D_adjudication", nearest


def summarize(gene_rows, call_rows, home_labels, fam_by_species, final_calls,
              searched_species):
    scored = [r for r in gene_rows if r["verdict"] != "not_searched"]
    locus_tier = [r for r in scored if r["tier"] == "home_locus"]
    reach = [r for r in scored if r["tier"] in REACHABLE_TIERS]
    ok_addr = ("recovered", "unresolvable_by_synteny")

    species = [s for s in fam_by_species
               if searched_species is None or s in searched_species]
    reach_sp = {r["species"] for r in reach}
    family_ok = ("on_locus", "on_address", "on_unassigned_label", "on_other_family_locus")
    found_sp = {c["species"] for c in call_rows if c["call_class"] in family_ok}

    def frac(n, d):
        return round(n / d, 4) if d else None

    cc = Counter(c["call_class"] for c in call_rows)
    absent_sp = (set(searched_species) - set(fam_by_species)) if searched_species else set()
    related = [r for r in locus_tier if not r.get("seed_unrelated")]
    assembled = [r for r in locus_tier if not r.get("assembly_limited")]
    # A call on an unassigned-label gene cannot be judged at locus level: it leaves the
    # locus-precision denominator instead of counting either way.
    locus_judged = len(call_rows) - cc["on_unassigned_label"]
    return {
        "home_loci": {lid: (lb or "UNCURATED_HOME_REGION") for lid, lb in home_labels.items()},
        "tier_counts": dict(Counter(r["tier"] for r in scored)),
        "verdict_counts": dict(Counter(r["verdict"] for r in scored)),
        "layer_census": {ly: sum(1 for r in reach if r["layer"] == ly) for ly in LAYERS},
        "locus_recall": frac(sum(r["verdict"] == "recovered" for r in locus_tier),
                             len(locus_tier)),
        "locus_recall_excl_seed_unrelated": frac(
            sum(r["verdict"] == "recovered" for r in related), len(related)),
        "locus_recall_excl_assembly_limited": frac(
            sum(r["verdict"] == "recovered" for r in assembled), len(assembled)),
        "assembly_limited_truth": sorted(r["gm"] for r in scored if r.get("assembly_limited")),
        "address_recall": frac(sum(r["verdict"] in ok_addr for r in reach), len(reach)),
        "reachable_genes": len(reach),
        "seed_unrelated_truth": sorted(r["gm"] for r in scored if r.get("seed_unrelated")),
        "family_species_recall": frac(len(found_sp & set(species)), len(species)),
        "reachable_species_recall": frac(len(found_sp & reach_sp), len(reach_sp)),
        "final_calls": len(call_rows),
        "call_classes": dict(cc),
        "locus_precision": frac(cc["on_locus"] + cc["on_address"], locus_judged),
        "family_precision": frac(sum(cc[k] for k in family_ok), len(call_rows)),
        # The absence test: the curated set records these species as lacking the
        # family, so any confident call there is a false positive by construction.
        "species_lacking_family": sorted(absent_sp),
        "final_calls_in_species_lacking_family": sum(1 for c in call_rows
                                                     if c["species"] in absent_sp),
    }


# ------------------------------------------------------------------ cli ---
def _write_tsv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def cmd_plan(args):
    fd = os.path.join(args.data, "final_dataset")
    genes = load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"))
    rows = build_plan(genes, load_addresses(os.path.join(fd, "LOCUS_ANCHORS.tsv")),
                      args.home_species)
    print(f"{'family':<10}{'rung':>5} {'split':<8}{'seed':<12}{'seed locus':<40}"
          f"{'locus':>6}{'addr':>6}{'unsrch':>7}{'noanch':>7}{'unasg':>6}{'reach_sp':>9}  warnings")
    for r in rows:
        if not r.get("seed_gm"):
            print(f"{r['family']:<10}  -- {r['warnings']}")
            continue
        print(f"{r['family']:<10}{r['ladder']!s:>5} {r['split']:<8}{r['seed_gm']:<12}"
              f"{r['seed_locus'][:39]:<40}{r['home_locus']:>6}{r['home_address']:>6}"
              f"{r['unsearched_home']:>7}{r['no_home_anchor']:>7}{r['unassigned_label']:>6}"
              f"{r['reachable_species']:>4}/{r['target_species_with_family']:<4}  "
              f"{r['warnings']}")
    if args.out:
        _write_tsv(args.out, rows)
        print(f"\nwrote {args.out}")
    return 0


# ------------------------------------------------------------- seedcheck ---
def protein_gene_map(gff_path):
    """{protein_id: (gene, chrom, start, product)} from a RefSeq GFF's CDS lines."""
    out = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = _attrs(f[8])
            pid = a.get("protein_id")
            if pid and pid not in out:
                prod = re.sub(r"%2C.*", "", a.get("product", ""))
                out[pid] = (a.get("gene", ""), f[0], int(f[3]), prod)
    return out


def rbh_summary(targets, top_gene, seed_gene):
    """Is the seed the best home-proteome homolog of the genes it should find?

    targets   : list of gene-model ids (the seed's reachable target genes)
    top_gene  : {gm: home gene id of that protein's best hit in the home proteome}
    seed_gene : the home gene id the seed protein maps to

    Returns (fraction of targets whose best home hit IS the seed gene, the most common
    other best hit or None, how many targets chose it). A fraction below 0.5 means the
    seed is not the reciprocal-best home homolog of its own locus -- the benchmark
    would be searching the wrong home neighbourhood by construction.
    """
    hits = [top_gene.get(t) for t in targets if top_gene.get(t)]
    if not hits:
        return None, None, 0
    same = sum(h == seed_gene for h in hits)
    others = Counter(h for h in hits if h != seed_gene)
    alt, n_alt = (others.most_common(1)[0] if others else (None, 0))
    return round(same / len(hits), 3), alt, n_alt


def run_mmseqs_hits(query_fa, target_faa, workdir, threads=4):
    """{query: [(target, bits, fident), ...]} via mmseqs easy-search. A nucleotide
    query is searched translated (blastx-like) against the protein target."""
    out = os.path.join(workdir, "hits.m8")
    subprocess.run(["mmseqs", "easy-search", query_fa, target_faa, out,
                    os.path.join(workdir, "tmp"), "-s", "7.5", "--max-seqs", "300",
                    "--threads", str(threads), "-v", "1",
                    "--format-output", "query,target,bits,fident"],
                   check=True, stdout=subprocess.DEVNULL)
    hits = defaultdict(list)
    with open(out) as fh:
        for line in fh:
            q, t, bits, fid = line.rstrip("\n").split("\t")
            hits[q].append((t, float(bits), float(fid)))
    return hits


def run_mmseqs_top_hits(query_faa, target_faa, workdir, threads=4):
    """{query: (target, bits, fident)} best hit per query, via mmseqs easy-search."""
    return {q: max(hs, key=lambda h: h[1])
            for q, hs in run_mmseqs_hits(query_faa, target_faa, workdir, threads).items()}


def cmd_seedcheck(args):
    fd = os.path.join(args.data, "final_dataset")
    genes = load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"))
    addresses = load_addresses(os.path.join(fd, "LOCUS_ANCHORS.tsv"))
    accessions = {}
    with open(os.path.join(args.data, "genomes", "ACCESSIONS.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accessions[species_key(r["species"])] = r["accession"]
    hdir = os.path.join(args.data, "genomes", accessions[args.home_species])
    proteome = args.home_proteome or os.path.join(hdir, "protein.faa")
    pmap = protein_gene_map(args.home_gff or glob.glob(os.path.join(hdir, "*.gff"))[0])

    proteins = {}
    for p in sorted(glob.glob(os.path.join(fd, "proteins", "*.faa"))):
        proteins.update(read_fasta(p))
    with tempfile.TemporaryDirectory() as tmp:
        q = os.path.join(tmp, "q.faa")
        with open(q, "w") as fh:
            for g in genes:
                if proteins.get(g["gm"]):
                    fh.write(f">{g['gm']}\n{proteins[g['gm']]}\n")
        best = run_mmseqs_top_hits(q, proteome, tmp, threads=args.threads)

    def home_gene(gm):
        hit = best.get(gm)
        return pmap.get(hit[0], ("",))[0] if hit else None

    rows = []
    for fam in sorted({g["family"] for g in genes}):
        for s in sorted((g for g in genes if g["family"] == fam
                         and g["species"] == args.home_species), key=lambda g: g["gm"]):
            tiers = assign_tiers(genes, fam, args.home_species, {s["locus"]}, addresses)
            reach = [t["gm"] for t in tiers if t["tier"] in REACHABLE_TIERS]
            sg = home_gene(s["gm"])
            frac, alt, n_alt = rbh_summary(reach, {t: home_gene(t) for t in reach}, sg)
            alt_info = ""
            if alt:
                pid = next((p for p, v in pmap.items() if v[0] == alt), None)
                if pid:
                    g_, c_, st_, prod = pmap[pid]
                    alt_info = f"{alt} {c_}:{st_} {prod}"
            rows.append({
                "family": fam, "seed_gm": s["gm"], "seed_home_gene": sg or "",
                "seed_self_identity": best.get(s["gm"], ("", 0, ""))[2],
                "reachable": len(reach), "seed_is_best_home_hit_frac": frac,
                "best_alternative": alt_info, "alternative_n": n_alt,
                "flag": "seed_not_rbh" if (frac is not None and frac < 0.5) else "",
            })
    for r in rows:
        print(f"{r['family']:<10}{r['seed_gm']:<12}{r['seed_home_gene']:<14}"
              f"reach={r['reachable']:<4}seed_best={r['seed_is_best_home_hit_frac']!s:<6}"
              f"{r['flag']:<13}{r['best_alternative'][:70]} (n={r['alternative_n']})")
    if args.out:
        _write_tsv(args.out, rows)
        print(f"\nwrote {args.out}")
    return 0


# ------------------------------------------------------------- callcheck ---
# Which home gene does a confident call really belong to? A call off the truth set is
# either an ortholog of the QUERY the curation did not record (or a seed problem), or
# a hit on a different home gene's ortholog: a genuine SynVoy precision error. The test
# is the half of RBH that needs no model: the home protein that best explains the
# called genomic window, by translated search.
MIN_CALL_BITS = 50


def extract_windows(fasta, wanted):
    """{key: sequence} of 1-based inclusive windows, reading the FASTA once.

    wanted: {chrom: [(start, end, key), ...]}
    """
    out, state = {}, {"name": None, "buf": []}

    def flush():
        if state["name"] in wanted:
            seq = "".join(state["buf"])
            for s, e, key in wanted[state["name"]]:
                out[key] = seq[s - 1:e]

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                state["name"], state["buf"] = line[1:].split()[0], []
            elif state["name"] in wanted:
                state["buf"].append(line.strip())
    flush()
    return out


def best_home_genes(hits, pmap, n=2):
    """[(home gene, best bits), ...] for one query's hits: top `n` DISTINCT genes, so
    two isoforms of one gene never read as a runner-up."""
    by_gene = {}
    for t, bits, _fid in hits:
        g = pmap.get(t, (t,))[0] or t
        by_gene[g] = max(bits, by_gene.get(g, 0.0))
    return sorted(by_gene.items(), key=lambda kv: (-kv[1], kv[0]))[:n]


def call_origin(best, seed_gene, min_bits=MIN_CALL_BITS):
    """query_gene: the call's window is best explained by the seed's own home gene.
    other_home_gene: by a different home gene (a paralog or an unrelated gene).
    no_home_hit: nothing in the home proteome explains it."""
    if not best or best[0][1] < min_bits:
        return "no_home_hit"
    return "query_gene" if best[0][0] == seed_gene else "other_home_gene"


def cmd_callcheck(args):
    accessions = {}
    with open(os.path.join(args.data, "genomes", "ACCESSIONS.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            accessions[species_key(r["species"])] = r["accession"]
    hdir = os.path.join(args.data, "genomes", accessions[args.home_species])
    proteome = args.home_proteome or os.path.join(hdir, "protein.faa")
    pmap = protein_gene_map(args.home_gff or glob.glob(os.path.join(hdir, "*.gff"))[0])
    calls_tsv = args.calls or os.path.join(args.run, "family_score.calls.tsv")
    with open(calls_tsv) as fh:
        calls = list(csv.DictReader(fh, delimiter="\t"))
    query = os.path.join(args.run, "intermediate", "query", "normalized_query.faa")

    wanted = defaultdict(lambda: defaultdict(list))
    for i, c in enumerate(calls):
        wanted[c["species"]][c["chrom"]].append((int(c["start"]), int(c["end"]), f"c{i}"))
    windows, missing = {}, set()
    for sp, per_chrom in sorted(wanted.items()):
        fa = os.path.join(args.targets, f"{sp}.fna")
        if not os.path.exists(fa):
            missing.add(sp)
            continue
        windows.update(extract_windows(fa, per_chrom))

    with tempfile.TemporaryDirectory() as tmp:
        qd, wd = os.path.join(tmp, "q"), os.path.join(tmp, "w")
        os.makedirs(qd)
        os.makedirs(wd)
        seed_hit = run_mmseqs_top_hits(query, proteome, qd, threads=args.threads)
        fna = os.path.join(tmp, "windows.fna")
        with open(fna, "w") as fh:
            for key, seq in windows.items():
                fh.write(f">{key}\n{seq}\n")
        hits = (run_mmseqs_hits(fna, proteome, wd, threads=args.threads)
                if windows else {})
    seed_gene = next((pmap.get(h[0], ("",))[0] for h in seed_hit.values()), "")

    rows = []
    for i, c in enumerate(calls):
        best = best_home_genes(hits.get(f"c{i}", []), pmap)
        origin = ("no_genome" if c["species"] in missing
                  else call_origin(best, seed_gene))
        rows.append(dict(
            c, seed_home_gene=seed_gene,
            best_home_gene=best[0][0] if best else "",
            best_bits=best[0][1] if best else "",
            second_home_gene=best[1][0] if len(best) > 1 else "",
            second_bits=best[1][1] if len(best) > 1 else "",
            origin=origin))

    census = Counter((r["call_class"], r["origin"]) for r in rows)
    print(f"callcheck {args.run}   seed home gene: {seed_gene or '?'}   calls: {len(rows)}")
    for (cls, origin), n in sorted(census.items()):
        print(f"  {cls:<24}{origin:<18}{n}")
    out = args.out or os.path.join(args.run, "family_score.callcheck.tsv")
    _write_tsv(out, rows)
    with open(re.sub(r"\.tsv$", "", out) + ".json", "w") as fh:
        json.dump({"seed_home_gene": seed_gene,
                   "census": {f"{k[0]}|{k[1]}": n for k, n in sorted(census.items())}},
                  fh, indent=2, sort_keys=True)
    print(f"wrote {out}")
    return 0


def cmd_stage(args):
    res = stage(args.data, args.seed, args.out, args.home_species,
                family=args.family, label=args.label, seed_faa=args.seed_faa,
                seed_coords=args.seed_coords)
    s = res["seed"]
    print(f"seed      {s['gm']}  {s['family']}  {s['locus']}  "
          f"{s['chrom']}:{s['start']}-{s['end']}  {s['prot_len']} aa")
    print(f"targets   {len(res['staged'])} staged in {args.out}/targets/")
    if res["missing"]:
        print(f"  ** no genome FASTA for: {', '.join(res['missing'])}")
    out = os.path.abspath(args.out)
    print("\nrun (pro mode; the home species must be set, the targets are species-named):")
    print(f"  nextflow run main.nf --mode pro --query {out}/query.faa \\\n"
          f"    --home_genome {out}/home/{args.home_species}.fna "
          f"--home_gff {out}/home/{args.home_species}.gff \\\n"
          f"    --home_species '{args.home_species.replace('_', ' ')}' "
          f"--target_genomes '{out}/targets/*' \\\n"
          f"    --keep_intermediate true --outdir <RESULTS>")
    return 0


def cmd_score(args):
    fd = os.path.join(args.data, "final_dataset")
    genes = load_manifest(os.path.join(fd, "EXPORT_MANIFEST.tsv"), exclude_families=())
    addresses = load_addresses(os.path.join(fd, "LOCUS_ANCHORS.tsv"))
    home_loci = load_home_loci(args.run)
    if not home_loci:
        sys.exit(f"no locus_N.bed under {args.run}: cannot tell which home loci were "
                 f"searched (was the run published with plot inputs?)")
    if args.seed_json:
        with open(args.seed_json) as fh:
            seed_info = json.load(fh)
    else:
        seed_info = load_seed_info(args.run)
    if seed_info and seed_info.get("family") not in (None, args.family):
        sys.exit(f"run was seeded for {seed_info['family']}, not {args.family}")
    genes = apply_seed(genes, seed_info, args.home_species)
    truth_cds = load_truth_cds(args.data)
    if truth_cds:
        for g in genes:
            g["cds"] = truth_cds.get(g["gm"], [])
    else:
        print("  ** no final_dataset/browser_gff/*.gff3: calls are matched to truth genes "
              "by span overlap, not coding overlap")
    home_labels = label_home_loci(home_loci, genes, args.family, args.home_species,
                                  slack=args.home_slack)
    gff_calls, spans = load_run_models(gather_run_gffs(args.run), args.home_species)
    report = load_report_records(args.run)
    if not report:
        print("  ** no synvoy_report.json records: every call will read as unadjudicated")
    seed = load_run_query(args.run)
    seed_evalues = {}
    if seed:
        proteins = {}
        for p in sorted(glob.glob(os.path.join(fd, "proteins", "*.faa"))):
            proteins.update(read_fasta(p))
        seed_evalues = seed_similarity(
            seed, proteins, [g["gm"] for g in genes if g["family"] == args.family])
    else:
        print("  ** no intermediate/query/normalized_query.faa: seed-homology flags skipped")
    gene_rows, call_rows, summary = score_run(
        genes, args.family, args.home_species, home_labels, addresses, gff_calls, spans,
        report, searched_species=load_searched_species(args.run),
        min_confidence=args.min_confidence, block_pad=args.block_pad,
        seed_evalues=seed_evalues,
        contig_lengths=load_contig_lengths(
            args.data, {g["species"] for g in genes if g["family"] == args.family},
            args.home_species))
    summary["seed_arm"] = (seed_info or {}).get("kind", "unrecorded")
    summary["seed"] = (seed_info or {}).get("seed", "")

    print(f"family {args.family}   seed {summary['seed'] or '?'} ({summary['seed_arm']})"
          f"   home loci: "
          + ", ".join(f"{k}={v}" for k, v in summary["home_loci"].items()))
    print(f"GOI models in GFFs: {len(gff_calls)}   adjudicated records: {len(report)}   "
          f"final confident calls: {summary['final_calls']}")
    print("\n=== MISSES BY LAYER (reachable genes only)")
    for ly, n in summary["layer_census"].items():
        print(f"  {ly:<18}{n}")
    print("\n=== SUMMARY")
    for k in ("tier_counts", "verdict_counts", "call_classes", "locus_recall",
              "locus_recall_excl_seed_unrelated", "locus_recall_excl_assembly_limited",
              "address_recall", "reachable_genes", "seed_unrelated_truth",
              "assembly_limited_truth", "family_species_recall", "reachable_species_recall",
              "locus_precision", "family_precision",
              "final_calls_in_species_lacking_family"):
        print(f"  {k:<34}{summary[k]}")
    print("\n=== ALL ADJUDICATED RECORDS: confidence (-> demoted class) x call class")
    for group, rows in summary["record_census"].items():
        print(f"  -- {group}")
        for tier, cc in rows.items():
            print(f"  {tier:<40}" + "  ".join(f"{k}={v}" for k, v in sorted(cc.items())))
    if args.out_prefix:
        _write_tsv(args.out_prefix + ".genes.tsv", gene_rows)
        _write_tsv(args.out_prefix + ".calls.tsv", call_rows)
        with open(args.out_prefix + ".summary.json", "w") as fh:
            json.dump(summary, fh, indent=2, sort_keys=True)
        print(f"\nwrote {args.out_prefix}.{{genes.tsv,calls.tsv,summary.json}}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DEFAULT_DATA,
                    help="Ant_Venoms data/ dir with final_dataset/ and genomes/")
    ap.add_argument("--home-species", default=DEFAULT_HOME)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plan", help="pre-registered tiers per family and seed")
    p.add_argument("--out", help="write the plan TSV here")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("stage", help="write query.faa + species-named targets/")
    p.add_argument("--seed", required=True,
                   help="curated home gene-model id, or a home-genome gene id for the "
                        "RBH arm (seeds: family_preregistration*.tsv)")
    p.add_argument("--family", help="required for a home-gene / ORF seed; checked otherwise")
    p.add_argument("--label", help="curated label the seed stands for, from the address "
                                   "audit (default: the family's only home label / the "
                                   "curated seed's own label)")
    p.add_argument("--seed-faa", help="one-protein FASTA of an unannotated home ORF seed")
    p.add_argument("--seed-coords", help="CHROM:START-END[:STRAND] of that ORF")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_stage)

    p = sub.add_parser("callcheck", help="which home gene best explains each confident call?")
    p.add_argument("--run", required=True, help="a scored run (family_score.calls.tsv)")
    p.add_argument("--targets", required=True, help="the staged targets/ dir of the run")
    p.add_argument("--calls", help="default: <run>/family_score.calls.tsv")
    p.add_argument("--home-proteome", help="default: genomes/<home acc>/protein.faa")
    p.add_argument("--home-gff", help="default: genomes/<home acc>/*.gff")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", help="default: <run>/family_score.callcheck.tsv")
    p.set_defaults(func=cmd_callcheck)

    p = sub.add_parser("seedcheck", help="is each seed the best home homolog of its locus?")
    p.add_argument("--home-proteome", help="default: genomes/<home acc>/protein.faa")
    p.add_argument("--home-gff", help="default: genomes/<home acc>/*.gff")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", help="write the check TSV here")
    p.set_defaults(func=cmd_seedcheck)

    p = sub.add_parser("score", help="score one run for one family")
    p.add_argument("--run", required=True, help="SynVoy --outdir of the run")
    p.add_argument("--family", required=True, help="e.g. FAM_APYR")
    p.add_argument("--min-confidence", choices=list(CONF_RANK), default="MEDIUM",
                   help="layer-C floor for GOI models (default MEDIUM)")
    p.add_argument("--block-pad", type=int, default=50000,
                   help="a block covers a gene within this distance of its model hull "
                        "(default = gap_search_window, the GOI search reach)")
    p.add_argument("--home-slack", type=int, default=5000,
                   help="bp slack when matching a home locus to a curated home gene")
    p.add_argument("--out-prefix", help="write <prefix>.genes.tsv/.calls.tsv/.summary.json")
    p.add_argument("--seed-json", help="score as if seeded by this seed.json (default: "
                                       "<run>/seed.json) -- e.g. one run under both its "
                                       "manifest label and its address label")
    p.set_defaults(func=cmd_score)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
