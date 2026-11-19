#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import parse_gff_attributes as _parse_gff_attrs  # noqa: E402,F401
from sequence_utils import clopper_pearson  # noqa: E402
from query_features import completeness_coverage  # noqa: E402
import synvoy_label  # noqa: E402
import synvoy_settings  # noqa: E402
import synvoy_specificity  # noqa: E402


KNOWN_SUFFIXES = [
    ".homology.tsv",
    ".scores.tsv",
    ".regions.bed",
    ".gff3",
    ".gff",
    ".faa",
    ".fna",
    # ".fasta"/".fa" must follow ".faa"/".fna" so the longer FASTA suffixes match
    # first. Without ".fa", Pro-mode "<genome>.fa" inputs canonicalize to
    # "<genome>.fa" while the (bare-stem) hull-rescue GFF "<genome>.hull_rescue.gff"
    # canonicalizes to "<genome>" — splitting one genome into a real record + an
    # empty phantom that falsely lands in goi_absent_genomes (docs/TODO_JUN.md QW1).
    ".fasta",
    ".fa",
    ".m8",
    # Rescue-pass GFFs (§1e strong-synteny, §1m/§17 hull) are named
    # "<genome>.rescue.gff" / "<genome>.hull_rescue.gff". Strip the infix so they
    # collapse to the real genome and are deduped/owned alongside the main calls —
    # otherwise a redundant rescue model leaks as a phantom genome "<genome>.hull_rescue".
    ".hull_rescue",
    ".rescue",
]


def _clean_json_text(text):
    cleaned = (text or "").strip()
    if cleaned.endswith(",]"):
        cleaned = cleaned[:-2] + "]"
    elif cleaned.endswith(","):
        cleaned = cleaned[:-1]
    return cleaned


def _safe_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def _safe_int(value, default=0):
    try:
        return int(float(value))
    except Exception:
        return default


def _is_true(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


# Per-locus region files are staged for the report with a "<locus>__" prefix
# (modules/stage_for_report.nf) so that same-named files from different home loci
# (e.g. GCA_036250285.1.fna.gff produced once per locus) no longer overwrite each
# other during the GENERATE_REPORT cp. The prefix is stripped for genome-id purposes
# and recovered as provenance for cross-locus GOI dedup (docs/TODO.md §1h).
_STAGING_PREFIX_RE = re.compile(r"^(locus_[A-Za-z0-9.]+|\d+)__")


def _strip_staging_prefix(name):
    return _STAGING_PREFIX_RE.sub("", name, count=1)


def staging_source_label(path_or_name):
    """Recover the per-locus staging tag (e.g. ``locus_3``) a region file came from.

    Falls back to the bare filename when the file carries no staging prefix
    (older runs / direct invocation), so provenance is always non-empty.
    """
    base = os.path.basename(path_or_name or "")
    m = _STAGING_PREFIX_RE.match(base)
    if m:
        return m.group(1)
    return base


def rescue_pass_of(path_or_name):
    """``hull_rescue`` or ``rescue`` for a GFF of a rescue pass (``<genome>.hull_rescue.gff``,
    ``<genome>.rescue.gff``), '' for one of the search. The figures are drawn from the search's
    GFFs alone (figure_marks)."""
    base = _strip_staging_prefix(os.path.basename(path_or_name or ""))
    for tag in ("hull_rescue", "rescue"):
        if base.endswith((f".{tag}.gff", f".{tag}.gff3")):
            return tag
    return ""


def canonical_genome_id(path_or_name):
    name = _strip_staging_prefix(os.path.basename(path_or_name or ""))
    m = re.search(r"(GC[AF]_\d+\.\d+)", name)
    if m:
        return m.group(1)

    changed = True
    while changed:
        changed = False
        for suffix in KNOWN_SUFFIXES:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                changed = True
    name = re.sub(r"(_new_genes|\.regions|\.candidates)$", "", name)
    return name


def load_qc_records(path):
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        content = _clean_json_text(fh.read())
    if not content:
        return []
    data = json.loads(content)
    return data if isinstance(data, list) else []


def summarize_qc(records):
    counts = Counter()
    failed_genomes = []
    thresholds = {}
    for rec in records:
        status = str(rec.get("status", "UNKNOWN")).upper()
        counts[status] += 1
        if status == "FAIL":
            failed_genomes.append({
                "genome": canonical_genome_id(rec.get("genome", "")),
                "raw_genome": rec.get("genome", ""),
                "msg": rec.get("msg", ""),
            })
        if not thresholds and isinstance(rec.get("thresholds"), dict):
            thresholds = rec["thresholds"]

    return {
        "total_genomes": sum(counts.values()),
        "pass": counts.get("PASS", 0),
        "fail": counts.get("FAIL", 0),
        "unknown": counts.get("UNKNOWN", 0),
        "failed_genomes": failed_genomes,
        "thresholds": thresholds,
    }


def count_fasta_records(path):
    try:
        with open(path) as fh:
            return sum(1 for line in fh if line.startswith(">"))
    except Exception:
        return 0


def summarize_fasta_outputs(fasta_files):
    genes_added_per_genome = {}
    for fasta_path in fasta_files:
        genome = canonical_genome_id(fasta_path)
        genes_added_per_genome[genome] = genes_added_per_genome.get(genome, 0) + count_fasta_records(fasta_path)
    return genes_added_per_genome


def _process_gffs_unified(gff_files):
    """Single read-pass over every per-locus region GFF that simultaneously produces
    the three outputs build_report needs: the annotation summary, the per-row GOI
    record list (used for dedup), and the (genome, locus, chrom)→HIGH-flanking-set map.

    Previously these three were computed in three separate functions, each iterating
    the entire ``gff_files`` list and re-opening every file (docs/TODO.md §1c). For the
    SP-family runs (~50K GFFs) that meant 3× the disk I/O and explained the timeouts
    Ivan's pipeline was hitting. The three public wrappers below remain for tests and
    external callers (each just discards the two outputs it doesn't need).
    """
    # ── summarize_annotations accumulators ────────────────────────────────
    per_genome_ann = {}
    role_counts = Counter()
    goi_confidence_counts = Counter()
    goi_class_counts = Counter()
    evidence_type_counts = Counter()
    model_status_counts = Counter()
    goi_evidence_counts = Counter()
    fallback_goi_annotations = 0
    total_annotations = 0
    # The rows of the genome-wide path's answers (GOIClass=genome_wide_strong_hit): a reading of
    # their own, counted here only, so that the counts above describe the other rows alone.
    genome_wide_answer_rows = 0
    # ── collect_goi_annotations accumulator ───────────────────────────────
    goi_records = []
    # ── collect_high_flanking_per_locus accumulator ───────────────────────
    flanking_seen = defaultdict(set)  # (genome, locus, chrom) -> set of gene names

    for gff_path in gff_files:
        genome = canonical_genome_id(gff_path)
        source = staging_source_label(gff_path)
        rescue_pass = rescue_pass_of(gff_path)
        ann_stats = per_genome_ann.setdefault(
            genome,
            {
                "genome": genome,
                "total_annotations": 0,
                "role_counts": Counter(),
                "goi_annotations": 0,
                "resolved_goi_annotations": 0,
                "ambiguous_goi_annotations": 0,
                "goi_confidence_counts": Counter(),
                "goi_class_counts": Counter(),
                "evidence_type_counts": Counter(),
                "model_status_counts": Counter(),
                "genome_wide_answer_annotations": 0,
            },
        )

        try:
            with open(gff_path) as fh:
                # NB: rstrip("\n") matches the collect_goi_annotations behaviour
                # which trimmed only the newline (not surrounding whitespace).
                # summarize_annotations used .strip() but parts of the GFF aren't
                # whitespace-sensitive past the 9-tab field split, so the looser
                # rstrip preserves identical semantics for the line.split("\t")
                # downstream.
                rows = [line.rstrip("\n").split("\t") for line in fh
                        if line.strip() and not line.startswith("#")]
                # One model is one row: an mRNA row, or a gene row that is no mRNA's parent.
                # The rescue passes write a gene row and an mRNA row per model (the CDS rows
                # hang on the mRNA); the search writes one row per model (a tandem copy is a
                # lone gene row).
                mrna_parents = {p for parts in rows if len(parts) >= 9 and parts[2] == "mRNA"
                                for p in (_parse_gff_attrs(parts[8]).get("Parent", "") or "").split(",") if p}
                cds_of_model = defaultdict(list)
                frames_of_model = defaultdict(list)   # the genomic frame of each CDS row (_same_orf)
                for parts in rows:
                    if len(parts) >= 9 and parts[2] == "CDS":
                        # A model's coding rows, for the dedup's test of a piece (_is_piece_of).
                        parent = _parse_gff_attrs(parts[8]).get("Parent", "")
                        try:
                            cs, ce = int(parts[3]), int(parts[4])
                            cds_of_model[parent].append((cs, ce))
                            phase = int(parts[7]) if parts[7].isdigit() else 0
                            frames_of_model[parent].append((ce - phase) % 3 if parts[6] == "-" else (cs + phase) % 3)
                        except ValueError:
                            pass
                        continue
                    if len(parts) < 9 or parts[2] not in {"mRNA", "gene"}:
                        continue
                    attrs = _parse_gff_attrs(parts[8])
                    model_id = attrs.get("ID", "")
                    if parts[2] == "gene" and model_id in mrna_parents:
                        continue
                    role = (attrs.get("SynVoyRole") or "").strip().lower()
                    if not role:
                        role = "goi" if model_id.startswith("GOI_") or attrs.get("SynVoy_Parent", "").startswith("GOI_") else "flanking"
                    confidence = (attrs.get("Confidence", "") or "UNKNOWN").upper()
                    goi_class = attrs.get("GOIClass", "")
                    evidence_type = attrs.get("EvidenceType", attrs.get("Type", "")) or "unknown"
                    model_status = attrs.get("ModelStatus", "")
                    genome_wide_answer = role == "goi" and goi_class == synvoy_label.GENOME_WIDE_CLASS

                    # — summarize_annotations bookkeeping (not for an answer of the genome-wide
                    # path, which is counted apart) —
                    if genome_wide_answer:
                        ann_stats["genome_wide_answer_annotations"] += 1
                        genome_wide_answer_rows += 1
                    else:
                        ann_stats["total_annotations"] += 1
                        ann_stats["role_counts"][role] += 1
                        ann_stats["evidence_type_counts"][evidence_type] += 1
                        if model_status:
                            ann_stats["model_status_counts"][model_status] += 1
                            model_status_counts[model_status] += 1
                        if role == "goi":
                            ann_stats["goi_confidence_counts"][confidence] += 1
                            if goi_class:
                                ann_stats["goi_class_counts"][goi_class] += 1

                        total_annotations += 1
                        role_counts[role] += 1
                        evidence_type_counts[evidence_type] += 1

                    if role == "goi":
                        if not genome_wide_answer:
                            ann_stats["goi_annotations"] += 1
                            goi_confidence_counts[confidence] += 1
                            if goi_class:
                                goi_class_counts[goi_class] += 1
                            goi_evidence_counts[evidence_type] += 1
                            if goi_class == "ambiguous_goi_family_member":
                                ann_stats["ambiguous_goi_annotations"] += 1
                            else:
                                ann_stats["resolved_goi_annotations"] += 1
                            if evidence_type in {"fallback_hit_span", "raw_hit", "rescued_exon"}:
                                fallback_goi_annotations += 1

                        # — collect_goi_annotations bookkeeping —
                        try:
                            start, end = int(parts[3]), int(parts[4])
                        except ValueError:
                            continue
                        if end < start:
                            start, end = end, start
                        goi_records.append({
                            "genome": genome,
                            "source": source,
                            # A rescue pass's row ('hull_rescue', 'rescue'), '' for the search's.
                            "rescue_pass": rescue_pass,
                            "chrom": parts[0],
                            "start": start,
                            "end": end,
                            "strand": parts[6],
                            "evidence": evidence_type,
                            "cds": cds_of_model[model_id],
                            "cds_frames": frames_of_model[model_id],
                            "confidence": confidence,
                            "goi_class": goi_class,
                            "identity": (attrs.get("Identity", "") or "").strip(),
                            "query_cov": (attrs.get("QueryCoverage", "") or "").strip(),
                            "core_cov": (attrs.get("CoreCoverage", "") or "").strip(),
                            "bits": (attrs.get("Bits", "") or "").strip(),
                            # An own copy of q left out of p* (synvoy_specificity, fix-phase d2).
                            "own_copy_set_aside": (attrs.get("OwnCopySetAside", "") or "").strip(),
                            "own_copy_bits": (attrs.get("OwnCopyBits", "") or "").strip(),
                            "target_gene": attrs.get("TargetGene", ""),
                            "model_status": model_status,
                            # A hit chain that the fallback gate refused, written because its
                            # block is significant: the tests it failed ("" otherwise).
                            "refusal_lifted": (attrs.get("RefusalLifted", "") or "").strip(),
                            "mrna_id": model_id,
                            # The judge's terms (label rework 2026-10-14, synvoy_label):
                            # why the label (Basis), the E-value in the expected place, the
                            # place and its block's E-value of gene order.
                            "basis": (attrs.get("Basis", "") or "").strip(),
                            "eeff": (attrs.get("Eeff", "") or "").strip(),
                            "place": (attrs.get("Place", "") or "").strip(),
                            "block_esyn": (attrs.get("BlockEsyn", "") or "").strip(),
                            # --report_copies_in_slot (preset plan PR-6): "" otherwise.
                            "copies_in_slot": (attrs.get("CopiesInSlot", "") or "").strip(),
                            "slot_rank": (attrs.get("SlotRank", "") or "").strip(),
                            # A model found only because the dispersed-GOI rescue made the
                            # search look there (DispersedRescue=1), with the rescue's rule, the
                            # anchors on the hit's sequence and the nearest one's distance.
                            "dispersed_rescue": (attrs.get("DispersedRescue", "") or "").strip() == "1",
                            "rescue_rule": (attrs.get("RescueRule", "") or "").strip(),
                            "rescue_anchors": (attrs.get("RescueAnchors", "") or "").strip(),
                            "rescue_anchor_bp": (attrs.get("RescueAnchorBp", "") or "").strip(),
                            # A model judged on its hit regions alone (no clean frame).
                            "frame_broken": (attrs.get("FrameBroken", "") or "").strip() == "1",
                            # A call of the genome-wide path (GenomeWide=1): the place's union
                            # coverage, its retrocopy test ("" when untested), the genome's median
                            # neighbour identity (recorded only) and the block whose window holds it.
                            "genome_wide": (attrs.get("GenomeWide", "") or "").strip() == "1",
                            "genome_wide_coverage": (attrs.get("GenomeWideCoverage", "") or "").strip(),
                            "retrocopy_candidate": (attrs.get("RetrocopyCandidate", "") or "").strip(),
                            "neighbour_median_identity": (attrs.get("NeighbourMedianIdentity", "") or "").strip(),
                            "near_block": (attrs.get("NearBlock", "") or "").strip(),
                            # A model of a tandem-copy lineage (SynVoy_Parent=GOI_copy_<n>): its
                            # Identity, QueryCoverage and Bits refer to the copy; these to the
                            # home protein ("" when not written).
                            "copy_lineage": (attrs.get("SynVoy_Parent", "")
                                             if attrs.get("SynVoy_Parent", "").startswith("GOI_copy_") else ""),
                            "home_identity": (attrs.get("HomeIdentity", "") or "").strip(),
                            "home_query_cov": (attrs.get("HomeQueryCoverage", "") or "").strip(),
                            "home_bits": (attrs.get("HomeBits", "") or "").strip(),
                        })

                    elif role == "flanking" and confidence == "HIGH":
                        # — collect_high_flanking_per_locus bookkeeping —
                        gene = (attrs.get("Name") or attrs.get("SynVoy_Parent")
                                or attrs.get("ID") or f"{parts[0]}:{parts[3]}")
                        flanking_seen[(genome, source, parts[0])].add(gene)
        except Exception as exc:
            print(f"Warning: Could not parse GFF {gff_path}: {exc}")

    # Build the summarize_annotations return shape
    per_genome_list = []
    genomes_without_goi = []
    genomes_with_only_ambiguous_goi = []
    genomes_with_only_low_goi = []
    for genome in sorted(per_genome_ann):
        stats = per_genome_ann[genome]
        row = {
            "genome": genome,
            "total_annotations": stats["total_annotations"],
            "role_counts": dict(stats["role_counts"]),
            "goi_annotations": stats["goi_annotations"],
            "resolved_goi_annotations": stats["resolved_goi_annotations"],
            "ambiguous_goi_annotations": stats["ambiguous_goi_annotations"],
            "goi_confidence_counts": dict(stats["goi_confidence_counts"]),
            "goi_class_counts": dict(stats["goi_class_counts"]),
            "evidence_type_counts": dict(stats["evidence_type_counts"]),
            "model_status_counts": dict(stats["model_status_counts"]),
            **({"genome_wide_answer_annotations": stats["genome_wide_answer_annotations"]}
               if stats["genome_wide_answer_annotations"] else {}),
        }
        per_genome_list.append(row)
        if stats["goi_annotations"] == 0:
            genomes_without_goi.append(genome)
        elif stats["resolved_goi_annotations"] == 0:
            genomes_with_only_ambiguous_goi.append(genome)
        # Every GOI model of the genome LOW by its label (an answer of the genome-wide path is
        # a MEDIUM model).
        if (stats["goi_annotations"] and set(stats["goi_confidence_counts"]) == {"LOW"}
                and not stats["genome_wide_answer_annotations"]):
            genomes_with_only_low_goi.append(genome)

    annotation_summary = {
        "per_genome": per_genome_list,
        "total_annotations": total_annotations,
        "role_counts": dict(role_counts),
        "goi_confidence_counts": dict(goi_confidence_counts),
        "goi_class_counts": dict(goi_class_counts),
        "goi_evidence_counts": dict(goi_evidence_counts),
        "evidence_type_counts": dict(evidence_type_counts),
        "model_status_counts": dict(model_status_counts),
        "fallback_goi_annotations": fallback_goi_annotations,
        "genomes_without_goi": genomes_without_goi,
        # Class-based (ambiguous_goi_family_member, the LOW class): deprecated, use
        # genomes_with_only_low_goi.
        "genomes_with_only_ambiguous_goi": genomes_with_only_ambiguous_goi,
        "genomes_with_only_low_goi": genomes_with_only_low_goi,
        **({"genome_wide_answer_rows": genome_wide_answer_rows} if genome_wide_answer_rows else {}),
    }
    flanking_per_locus = {key: len(genes) for key, genes in flanking_seen.items()}
    return annotation_summary, goi_records, flanking_per_locus


def summarize_annotations(gff_files):
    """Backwards-compat wrapper. Prefer ``_process_gffs_unified`` for build_report
    so the same GFF isn't read three times (docs/TODO.md §1c)."""
    return _process_gffs_unified(gff_files)[0]


def summarize_region_scores(score_files):
    per_genome = {}
    confidence_counts = Counter()
    selection_reason_counts = Counter()
    goi_anchor_regions = 0
    total_regions = 0

    for score_path in score_files:
        genome = canonical_genome_id(score_path)
        stats = per_genome.setdefault(
            genome,
            {
                "genome": genome,
                "total_regions": 0,
                "confidence_counts": Counter(),
                "selection_reason_counts": Counter(),
                "goi_anchor_regions": 0,
                "best_score": None,
            },
        )
        try:
            with open(score_path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    if not row:
                        continue
                    confidence = (row.get("confidence") or "UNKNOWN").upper()
                    selection_reason = row.get("selection_reason", "") or "unknown"
                    score = _safe_float(row.get("score"), default=0.0)

                    total_regions += 1
                    stats["total_regions"] += 1
                    stats["confidence_counts"][confidence] += 1
                    stats["selection_reason_counts"][selection_reason] += 1
                    confidence_counts[confidence] += 1
                    selection_reason_counts[selection_reason] += 1

                    if _is_true(row.get("is_goi_anchor")) or _is_true(row.get("goi_overlap")):
                        stats["goi_anchor_regions"] += 1
                        goi_anchor_regions += 1

                    if stats["best_score"] is None or score > stats["best_score"]:
                        stats["best_score"] = score
        except Exception as exc:
            print(f"Warning: Could not parse scores TSV {score_path}: {exc}")

    per_genome_list = []
    for genome in sorted(per_genome):
        stats = per_genome[genome]
        per_genome_list.append({
            "genome": genome,
            "total_regions": stats["total_regions"],
            "confidence_counts": dict(stats["confidence_counts"]),
            "selection_reason_counts": dict(stats["selection_reason_counts"]),
            "goi_anchor_regions": stats["goi_anchor_regions"],
            "best_score": stats["best_score"],
        })

    return {
        "per_genome": per_genome_list,
        "total_regions": total_regions,
        "confidence_counts": dict(confidence_counts),
        "selection_reason_counts": dict(selection_reason_counts),
        "goi_anchor_regions": goi_anchor_regions,
    }


def _staging_locus(path_or_name):
    """The ``locus_N`` staging tag of a staged file, '' without one."""
    m = _STAGING_PREFIX_RE.match(os.path.basename(path_or_name or ""))
    return m.group(1) if m and m.group(1).startswith("locus_") else ""


def collect_region_rows(score_files):
    """Every row of every region scores TSV (cluster_grs.py), as dicts: genome, locus (the
    staging tag, '' without one), region_name, region_class, chrom, start, end (the BED's 0-based,
    end-exclusive span), confidence (upper case), score, goi_overlap, and the columns that name
    what the region holds (goi_model_confidence, goi_model_class, goi_hit_best_identity,
    goi_hit_best_evalue; '' or None for a table written without them). Unreadable files are
    skipped."""
    rows = []
    for path in sorted(score_files):
        genome, locus = canonical_genome_id(path), _staging_locus(path)
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if not row or not row.get("chrom"):
                        continue
                    rows.append({
                        "genome": genome, "locus": locus,
                        "region_name": row.get("region_name", "") or "",
                        "region_class": row.get("region_class", "") or "",
                        "chrom": row["chrom"],
                        "start": _safe_int(row.get("start")), "end": _safe_int(row.get("end")),
                        "confidence": (row.get("confidence") or "").upper(),
                        "score": _optional_float(row.get("score")),
                        "goi_overlap": _is_true(row.get("goi_overlap")) or _is_true(row.get("is_goi_anchor")),
                        "goi_model_confidence": (row.get("goi_model_confidence") or "").upper(),
                        "goi_model_class": row.get("goi_model_class", "") or "",
                        "goi_hit_best_identity": _optional_float(row.get("goi_hit_best_identity")),
                        "goi_hit_best_evalue": _optional_float(row.get("goi_hit_best_evalue")),
                    })
        except (OSError, csv.Error):
            continue
    return rows


# Why a gene-bearing region of MEDIUM or better holds no counted call (goi_regions_without_call).
REGION_NO_CALL_REASONS = (
    "call_not_counted",          # a HIGH or MEDIUM model the report does not count (ownership gave
                                 # it to another home gene, the coverage rule, or an answer of the
                                 # genome-wide path)
    "model_not_counted",         # only AMBIGUOUS or LOW models
    "refused_by_synteny_gate",   # no model; the search refused the hit (rejected_by names the gate)
    "dropped_by_rescue",         # no model; the rescue dropped the hit (rejected_by names the test)
    "no_record",                 # no model, no record of the hit
)


def goi_regions_without_call(region_rows, goi_annotations, records, refused, dropped,
                             demoted_classes, pad=20000, report_dropped=()):
    """The gene-bearing regions (goi_overlap: within ``pad`` of a GOI model of the search, or of a
    raw GOI hit where the search wrote no GOI model) labelled HIGH or MEDIUM, or of the class
    goi_hit_no_call, that hold no counted call: no counted record (_is_counted, after ownership,
    the coverage rule and the drops; of any home locus) of the same genome and sequence within
    ``pad`` of the region. ``pad`` is cluster_grs's default --goi_padding. Each entry names the
    best uncounted model there (``model``: a record first, else a record the report dropped, of
    ``report_dropped`` ((record, why) pairs of drop_genome_wide_strong_hits and
    drop_gene_order_only; ``model.report_drop`` names why, '' for a record that stands), else a
    LOW row of ``goi_annotations``), the strongest record of the hit the search kept
    (``refused_hit``: a refused hit of ``refused``, else a dropped one of ``dropped``, overlapping
    the region, of the same genome and, when both name one, the same home locus) and the first of
    REGION_NO_CALL_REASONS that applies. Sorted by genome, locus, sequence, start, end, name."""
    def near(rec, row):
        return int(rec["start"]) - 1 < row["end"] + pad and row["start"] - pad < int(rec["end"])

    recs_by_gc = defaultdict(list)
    for r in records or ():
        recs_by_gc[(r.get("genome"), r.get("chrom"))].append(r)
    drop_why = {}
    report_dropped_by_gc = defaultdict(list)
    for r, why in report_dropped or ():
        drop_why[id(r)] = why
        report_dropped_by_gc[(r.get("genome"), r.get("chrom"))].append(r)
    low_by_gc = defaultdict(list)
    for a in goi_annotations or ():
        if a.get("confidence") == "LOW":
            low_by_gc[(a.get("genome"), a.get("chrom"))].append(a)
    hits_by_gc = defaultdict(list)
    for outcome, rows in (("refused", refused), ("dropped", dropped)):
        for h in rows or ():
            hits_by_gc[(h.get("genome"), h.get("chrom"))].append((outcome, h))

    out = []
    for row in region_rows or ():
        if not row.get("goi_overlap"):
            continue
        if row["confidence"] not in ("HIGH", "MEDIUM") and row["region_class"] != "goi_hit_no_call":
            continue
        key = (row["genome"], row["chrom"])
        here = [r for r in recs_by_gc.get(key, ()) if near(r, row)]
        if any(_is_counted(r, demoted_classes) for r in here):
            continue
        order = lambda m: (_record_rank(m), int(m["start"]), m.get("mrna_id", ""))  # noqa: E731
        model = min(here, key=order) if here else None
        if model is None:
            gone = [r for r in report_dropped_by_gc.get(key, ()) if near(r, row)]
            model = min(gone, key=order) if gone else None
        if model is None:
            lows = [a for a in low_by_gc.get(key, ()) if near(a, row)]
            model = min(lows, key=order) if lows else None
        hits = [(o, h) for o, h in hits_by_gc.get(key, ())
                if h["start"] < row["end"] and row["start"] < h["end"]
                and not (row["locus"] and h.get("locus") and h["locus"] != row["locus"])]
        hit = (min(hits, key=lambda t: (t[0] != "refused", -(t[1].get("identity") or 0.0),
                                        -(t[1].get("bits") or 0.0), t[1]["start"], t[1].get("query", "")))
               if hits else None)
        if model is not None and model.get("confidence") in ("HIGH", "MEDIUM"):
            reason = "call_not_counted"
        elif model is not None:
            reason = "model_not_counted"
        elif hit is not None:
            reason = "refused_by_synteny_gate" if hit[0] == "refused" else "dropped_by_rescue"
        else:
            reason = "no_record"
        out.append({
            "genome": row["genome"], "locus": row["locus"], "region_name": row["region_name"],
            "region_class": row["region_class"], "chrom": row["chrom"],
            "start": row["start"], "end": row["end"],
            "region_confidence": row["confidence"], "score": row["score"],
            "goi_model_confidence": row["goi_model_confidence"],
            "goi_hit_best_identity": row["goi_hit_best_identity"],
            "goi_hit_best_evalue": row["goi_hit_best_evalue"],
            "model": ({"confidence": model.get("confidence", ""), "goi_class": model.get("goi_class", ""),
                       "mrna_id": model.get("mrna_id", ""), "start": model["start"], "end": model["end"],
                       "report_drop": drop_why.get(id(model), "")}
                      if model is not None else None),
            "refused_hit": ({"outcome": hit[0], "identity": hit[1].get("identity"),
                             "evalue": hit[1].get("evalue", ""), "bits": hit[1].get("bits"),
                             "rejected_by": hit[1].get("rejected_by", ""), "query": hit[1].get("query", ""),
                             "start": hit[1]["start"], "end": hit[1]["end"]}
                            if hit is not None else None),
            "reason": reason,
        })
    out.sort(key=lambda e: (e["genome"], e["locus"], e["chrom"], e["start"], e["end"], e["region_name"]))
    return out


def scan_dir_by_suffix(dir_path, suffix_buckets):
    """Single ``os.scandir()`` pass that classifies entries by filename suffix.

    Replaces the previous N×``glob.glob(dir + '*.ext')`` pattern, which made one
    directory walk per extension (4 walks of ``staged_results/regions/`` per report
    at SP-family scale — ~50K files each). At that scale the redundant walks were a
    real fraction of the report's wall time (docs/TODO.md §1c).

    Args:
        dir_path: Directory to scan. Missing directories return empty buckets.
        suffix_buckets: ``{bucket_name: (suffix1, suffix2, ...)}``. First matching
            suffix wins; pass the more-specific suffix first if there's any overlap
            (e.g. ".scores.tsv" before ".tsv").

    Returns:
        ``{bucket_name: [str]}`` of full paths.
    """
    result = {name: [] for name in suffix_buckets}
    if not os.path.isdir(dir_path):
        return result
    try:
        with os.scandir(dir_path) as it:
            for entry in it:
                try:
                    if not entry.is_file(follow_symlinks=True):
                        continue
                except OSError:
                    continue
                name = entry.name
                for bucket, suffixes in suffix_buckets.items():
                    if any(name.endswith(s) for s in suffixes):
                        result[bucket].append(entry.path)
                        break
    except OSError:
        pass
    # scandir order is file-system dependent (ext4 hash order vs. e.g. GPFS); the
    # report must not depend on which machine staged the files.
    for bucket in result:
        result[bucket].sort()
    return result


def _dir_diagnostics(dir_path, patterns, sample_size=5):
    exists = os.path.isdir(dir_path)
    entries = []
    matches = []
    if exists:
        try:
            entries = sorted(os.listdir(dir_path))
        except OSError as exc:
            entries = []
            return {
                "path": dir_path,
                "exists": True,
                "readable": False,
                "error": str(exc),
                "entry_count": 0,
                "pattern_match_count": 0,
                "patterns": list(patterns),
                "sample_entries": [],
                "sample_matches": [],
            }
        for pat in patterns:
            # sorted: glob follows the file system's listing order, so two runs on LRZ
            # scratch never gave byte-identical reports (fix-phase item 14)
            matches.extend(sorted(glob.glob(os.path.join(dir_path, pat))))
    return {
        "path": dir_path,
        "exists": exists,
        "readable": exists,
        "entry_count": len(entries),
        "pattern_match_count": len(matches),
        "patterns": list(patterns),
        "sample_entries": entries[:sample_size],
        "sample_matches": [os.path.basename(m) for m in matches[:sample_size]],
    }


NONSYNTENIC_SUFFIX = ".nonsyntenic.tsv"
GENOME_WIDE_SUFFIX = ".genome_wide.tsv"
NOT_SEARCHED_SUFFIX = ".not_searched.tsv"
TARGET_POLICY_SUFFIX = ".target_policy.tsv"
_CONTEXT_ORDER = {"complete": 0, "left": 1, "right": 1, "both": 1, "NA": 2}


def collect_target_policy(paths):
    """{genome: row} of the per-genome tables the search writes
    (iterative_search_runner._write_target_policy): median flanking identity, blocks,
    significant blocks, the best block's E_syn and window context. With several home loci a
    genome has one row per locus; the row kept is the one with the most informative
    context (complete over truncated over no block), then the most significant block. The
    limit_* columns are left out here; collect_limits sums them."""
    out = {}
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    genome = canonical_genome_id(
                        row.get("genome") or os.path.basename(path).split(TARGET_POLICY_SUFFIX)[0])
                    row = {k: v for k, v in row.items()
                           if k and not k.startswith(LIMIT_PREFIX) and k not in SW_BACKEND_COLUMNS}
                    row["genome"] = genome
                    key = (_CONTEXT_ORDER.get(row.get("best_block_context") or "NA", 2),
                           _safe_float(row.get("best_block_esyn")) if row.get("best_block_esyn") not in (None, "", "NA") else float("inf"))
                    old = out.get(genome)
                    if old is None or key < old[0]:
                        out[genome] = (key, row)
        except (OSError, csv.Error):
            continue
    return {g: row for g, (_key, row) in sorted(out.items())}


SW_BACKEND_COLUMNS = ("sw_backend", "sw_parasail_version", "sw_simd")


def collect_sw_backend(paths):
    """What ran the Smith-Waterman search (fix-phase item 18): the sw_backend, sw_parasail_version
    and sw_simd columns of the target-policy rows, as {backend, parasail_version, simd} when
    every row agrees, else {'mixed': [...]} (one entry per combination, e.g. the tasks of a
    cluster run on different nodes). None when no row carries them (a run before 2026-10-29)."""
    seen = []
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if "sw_backend" not in row:
                        continue
                    entry = {"backend": row.get("sw_backend") or "",
                             "parasail_version": row.get("sw_parasail_version") or "",
                             "simd": row.get("sw_simd") or ""}
                    if entry not in seen:
                        seen.append(entry)
        except (OSError, csv.Error):
            continue
    if not seen:
        return None
    return seen[0] if len(seen) == 1 else {"mixed": seen}


DECOY_COUNT_COLUMNS = ("decoy_cases", "decoy_high", "decoy_medium_plus", "decoy_ambiguous_plus")
DECOY_GW_COUNT_COLUMNS = ("decoy_gw_cases", "decoy_gw_answers")
LIMIT_PREFIX = "limit_"


def collect_limits(paths):
    """What each limit of the search dropped (profile handoff, stage 2; counts only): the
    limit_* columns of every target-policy row (one row per home locus and searched genome),
    summed, with the genomes where a count is not 0. iterative_search_runner.LIMIT_KEYS and
    docs/OUTPUT.md say what each counts. None when no row carries them."""
    totals, genomes, order = {}, {}, []
    rows = 0
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    keys = [k for k in row if k and k.startswith(LIMIT_PREFIX)]
                    if not keys:
                        continue
                    rows += 1
                    genome = canonical_genome_id(
                        row.get("genome") or os.path.basename(path).split(TARGET_POLICY_SUFFIX)[0])
                    for k in keys:
                        name = k[len(LIMIT_PREFIX):]
                        if name not in totals:
                            order.append(name)
                            totals[name], genomes[name] = 0, set()
                        n = int(float(row.get(k) or 0))
                        totals[name] += n
                        if n:
                            genomes[name].add(genome)
        except (OSError, csv.Error, ValueError):
            continue
    if not rows:
        return None
    return {"genome_searches": rows,
            "counts": {k: {"total": totals[k], "genomes": sorted(genomes[k])} for k in order}}


def collect_decoy_rate(paths):
    """The run's false-call rate from decoy mode (preset plan PR-8): the decoy cases of every
    target-policy row (one row per home locus and searched genome; a case = home GOI x decoy
    kind, as iterative_search_runner.decoy_case_counts counts it), summed. Per label: the
    calls and the upper end of the exact (Clopper-Pearson) 95 % interval; no bare rate, so
    that 0 calls is never read as "0 %". None when no row carries decoy counts.

    With the genome-wide path, ``genome_wide``: its decoy pass apart (decoy_gw_cases: a case =
    home GOI x decoy kind searched in the genome place; decoy_gw_answers: cases with a decoy
    answer and no window decoy call of MEDIUM or better), the answers with their interval."""
    totals = dict.fromkeys(DECOY_COUNT_COLUMNS + DECOY_GW_COUNT_COLUMNS, 0)
    rows = gw_rows = 0
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if not str(row.get("decoy_cases") or "").strip():
                        continue
                    rows += 1
                    for k in DECOY_COUNT_COLUMNS:
                        totals[k] += int(float(row.get(k) or 0))
                    if str(row.get("decoy_gw_cases") or "").strip():
                        gw_rows += 1
                        for k in DECOY_GW_COUNT_COLUMNS:
                            totals[k] += int(float(row.get(k) or 0))
        except (OSError, csv.Error, ValueError):
            continue
    if not rows:
        return None
    n = totals["decoy_cases"]
    out = {"cases": n, "genome_searches": rows}
    for key, label in (("decoy_high", "high"), ("decoy_medium_plus", "medium_or_better"),
                       ("decoy_ambiguous_plus", "ambiguous_or_better")):
        k = totals[key]
        out[label] = {"calls": k,
                      "ci95_upper": round(clopper_pearson(k, n)[1], 4) if n else None}
    if gw_rows:
        n_gw, k_gw = totals["decoy_gw_cases"], totals["decoy_gw_answers"]
        out["genome_wide"] = {"cases": n_gw, "genome_searches": gw_rows,
                              "answers": {"calls": k_gw,
                                          "ci95_upper": round(clopper_pearson(k_gw, n_gw)[1], 4) if n_gw else None}}
    return out


def genome_absence_categories(records, target_policy, demoted_classes):
    """What a genome without a counted call means (preset plan PR-2):

      nearest_paralog_only        its best GOI record (the judge's order) is AMBIGUOUS
                                  with Basis=not_specific: only a model that another
                                  home gene explains as well
      absent_context_complete     the best block (smallest E_syn) was searched with its
                                  whole window inside the scaffold
      unknown_context_truncated   that window ran past a scaffold end: the place where
                                  the gene should sit is partly missing from the assembly
      absent_no_block             the search found no block in this genome
      strong_hit_without_synteny  it holds an answer of the genome-wide path (a strong hit
                                  outside every searched window; not counted)
      strong_hit_near_block       it holds such an answer inside a block's searched window,
                                  outside its GOI zone (NearBlock)
    The first is a property of the calls, the others of the search; a genome can be
    nearest_paralog_only (read on the records other than the answers) and one of the others. A
    genome with an answer is not absent and is in none of the three absent lists."""
    counted = counted_genomes(records, demoted_classes)
    answered = {}
    for r in genome_wide_answer_records(records, demoted_classes):
        answered.setdefault(r.get("genome"), r)
    by_genome = {}
    for r in records:
        if not r.get("genome_wide_strong_hit"):
            by_genome.setdefault(r.get("genome"), []).append(r)
    genomes = sorted((set(target_policy) | set(by_genome) | set(answered)) - counted - {None})
    out = {"nearest_paralog_only": [], "absent_context_complete": [],
           "unknown_context_truncated": [], "absent_no_block": [],
           "strong_hit_without_synteny": [], "strong_hit_near_block": []}
    for g in genomes:
        recs = by_genome.get(g, [])
        if recs:
            best = min(recs, key=_record_rank)
            if best.get("confidence") == "AMBIGUOUS" and best.get("basis") == "not_specific":
                out["nearest_paralog_only"].append(g)
        if g in answered:
            out["strong_hit_near_block" if str(answered[g].get("near_block") or "").strip()
                else "strong_hit_without_synteny"].append(g)
            continue
        ctx = (target_policy.get(g) or {}).get("best_block_context") or "NA"
        if ctx == "complete":
            out["absent_context_complete"].append(g)
        elif ctx in ("left", "right", "both"):
            out["unknown_context_truncated"].append(g)
        else:
            out["absent_no_block"].append(g)
    return out


def collect_not_searched(paths):
    """Genomes the search could not search (iterative_search_runner._write_not_searched),
    as {genome: [reason, ...]}. Unreadable files are skipped: diagnostic, not fatal."""
    out = {}
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    genome = canonical_genome_id(
                        row.get("genome") or os.path.basename(path).split(NOT_SEARCHED_SUFFIX)[0])
                    reason = row.get("reason") or "unknown"
                    if reason not in out.setdefault(genome, []):
                        out[genome].append(reason)
        except (OSError, csv.Error):
            continue
    return out


NONSYNTENIC_REFUSAL_NOTE = (
    "GOI hits that passed the identity/e-value/alignment-length bar but were refused by a "
    "synteny gate. Not orthology calls: at long divergences a translocated ortholog and a "
    "paralog are indistinguishable by gene order. Listed so a genome with a strong "
    "unplaceable hit is not read as a genome with no hit. `places` merges the rows of one "
    "home locus that lie within 5 kb of each other; a place that a counted call of the run "
    "overlaps (`modelled_in_run`) is left out of `goi_found_but_not_syntenic`.")
NONSYNTENIC_DROPPED_NOTE = (
    "GOI hits of the genome-wide search outside every searched GOI zone that the dispersed-GOI "
    "rescue dropped, with the test each failed (`rejected_by`: dispersed_identity, "
    "dispersed_evalue or dispersed_alnlen; dead_zone: the hit lay in the GOI zone of a block "
    "that was not searched). The strongest 200 per genome search; `hits` counts them all. "
    "Diagnostics, not candidates.")


def _optional_int(value):
    """An int, or None for a missing or empty field (a table written before the column)."""
    if value is None or str(value).strip() == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _optional_float(value):
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _nonsyntenic_locus(path, raw_genome):
    """The home locus of a nonsyntenic table without a locus column: the file-name prefix
    before ``<genome>.nonsyntenic.tsv`` (iterative_search_runner names it
    ``<locus>_<genome>.nonsyntenic.tsv``); '' when it cannot be read."""
    base = os.path.basename(path)
    stem = base[:-len(NONSYNTENIC_SUFFIX)] if base.endswith(NONSYNTENIC_SUFFIX) else base
    if raw_genome and stem.endswith(raw_genome):
        return stem[:-len(raw_genome)].rstrip("_")
    return ""


def collect_nonsyntenic_candidates(paths, dropped_out=None):
    """§1x: load the GOI hits that cleared the quality bar and were then refused by a
    synteny gate.

    These are the run's *false negatives that it knows about*. Refusing a hit whose
    neighbourhood cannot be confirmed is the correct call — at 250 My a translocated
    ortholog is genuinely indistinguishable from a paralog by synteny alone — but the
    old behaviour reported such a genome as "no ortholog" while holding a 54-72 %
    identity hit sitting inside the real gene. `goi_absent` and
    `goi_found_but_not_syntenic` are scientifically opposite claims and must not share
    a field.

    The same tables list the GOI hits the rescue dropped (column ``outcome`` = dropped:
    below its quality bar, or in the GOI zone of a block that was not searched, outside every
    searched GOI zone); they are appended to ``dropped_out`` when it is a list, and are never
    refused candidates. A row with ``outcome`` = searched (a refused hit that the GOI zone of a
    block seeded by the rescue then held) is left out altogether. A table without the column
    (before 2026-11-16) holds refused hits only; one without a ``locus`` column takes the locus
    from its file name.

    Returns (records, per_genome_map) of the refused hits. Malformed or unreadable files are
    skipped rather than failing the report — this is diagnostic output, never a reason to lose
    a run.
    """
    records, per_genome = [], {}
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if not row.get("chrom"):
                        continue
                    raw_genome = row.get("genome") or ""
                    genome = canonical_genome_id(
                        raw_genome or os.path.basename(path).split(NONSYNTENIC_SUFFIX)[0])
                    rec = {
                        "genome": genome,
                        "locus": row.get("locus") or _nonsyntenic_locus(path, raw_genome),
                        "query": row.get("query", ""),
                        "chrom": row.get("chrom", ""),
                        "start": _safe_int(row.get("start")),
                        "end": _safe_int(row.get("end")),
                        "strand": row.get("strand", "") or "",
                        "qstart": _optional_int(row.get("qstart")),
                        "qend": _optional_int(row.get("qend")),
                        "bits": _optional_float(row.get("bits")),
                        "identity": _safe_float(row.get("identity")),
                        "evalue": row.get("evalue", ""),
                        "alnlen": _safe_int(row.get("alnlen")),
                        "rejected_by": row.get("rejected_by", ""),
                        "nearest_anchor_gene": row.get("nearest_anchor_gene", ""),
                        # Empty when the hit's chromosome carries no anchor: null, not 0
                        # ("overlaps an anchor"; fix item 6).
                        "nearest_anchor_bp": (_safe_int(row.get("nearest_anchor_bp"))
                                              if str(row.get("nearest_anchor_bp", "")).strip() else None),
                        "nearest_anchor_same_chrom":
                            str(row.get("nearest_anchor_same_chrom", "")).lower() == "true",
                        "envelope": row.get("envelope", "") or "",
                        "anchors_on_sequence": _optional_int(row.get("anchors_on_sequence")),
                        "genes_on_sequence": _optional_int(row.get("genes_on_sequence")),
                        "modelled_here": row.get("modelled_here", "") or "",
                        "modelled_by": row.get("modelled_by", "") or "",
                    }
                    outcome = row.get("outcome") or "refused"
                    if outcome == "searched":
                        # A refused hit that a block seeded by the rescue then searched: no
                        # candidate, no count, no place.
                        continue
                    if outcome != "refused":
                        if dropped_out is not None:
                            dropped_out.append(rec)
                        continue
                    records.append(rec)
                    per_genome.setdefault(genome, []).append(rec)
        except (OSError, csv.Error):
            continue
    records.sort(key=lambda r: (-r["identity"], r["genome"], r["chrom"], r["start"]))
    if dropped_out is not None:
        dropped_out.sort(key=lambda r: (r["genome"], r["locus"], -(r["bits"] or 0.0),
                                        _safe_float(r["evalue"], float("inf")), r["chrom"],
                                        r["start"], r["end"], r["query"]))
    return records, per_genome


# Refused GOI rows of one home locus and sequence closer than this are one place
# (nonsyntenic_places).
PLACE_MERGE_BP = 5000


def _row_strength(r):
    """Best first: most bits, then the smallest E, the highest identity, the query, the start."""
    return (-(r.get("bits") or 0.0), _safe_float(r.get("evalue"), float("inf")),
            -(r.get("identity") or 0.0), r.get("query", ""), r.get("start", 0))


def nonsyntenic_places(records, merge_bp=PLACE_MERGE_BP):
    """The refused hits of one home locus, genome and sequence that overlap or lie within
    ``merge_bp`` of each other, as one place (several GOI queries, such as the home protein,
    its exons and earlier genomes' models, hit one place, one row each). A place: its span
    (0-based, end exclusive), the strand, identity, E and bits of its best row (most bits,
    then smallest E), the nearest flanking anchor of that row, its rows (``n_rows``), queries
    and tests (``rejected_by``), and ``modelled_here``, the best label of a model of the same
    search that overlaps one of its rows (with ``modelled_by``)."""
    groups = defaultdict(list)
    for r in records:
        groups[(r["genome"], r.get("locus", ""), r["chrom"])].append(r)
    places = []
    for (genome, locus, chrom), rows in sorted(groups.items()):
        rows = sorted(rows, key=lambda r: (r["start"], r["end"], r.get("query", "")))
        runs, current, run_end = [], [], None
        for r in rows:
            if current and r["start"] - run_end <= merge_bp:
                current.append(r)
                run_end = max(run_end, r["end"])
            else:
                if current:
                    runs.append(current)
                current, run_end = [r], r["end"]
        if current:
            runs.append(current)
        for run in runs:
            best = min(run, key=_row_strength)
            modelled = [r for r in run if r.get("modelled_here")]
            top_model = (min(modelled, key=lambda r: (-synvoy_label.LABEL_RANK.get(r["modelled_here"], 0),
                                                      r.get("modelled_by", "")))
                         if modelled else None)
            places.append({
                "genome": genome, "locus": locus, "chrom": chrom,
                "start": min(r["start"] for r in run), "end": max(r["end"] for r in run),
                "strand": best.get("strand", ""), "identity": best.get("identity"),
                "evalue": best.get("evalue", ""), "bits": best.get("bits"),
                "nearest_anchor_gene": best.get("nearest_anchor_gene", ""),
                "nearest_anchor_bp": best.get("nearest_anchor_bp"),
                "nearest_anchor_same_chrom": best.get("nearest_anchor_same_chrom", False),
                "n_rows": len(run),
                "queries": sorted({r.get("query", "") for r in run}),
                "rejected_by": sorted({r.get("rejected_by", "") for r in run}),
                "modelled_here": top_model["modelled_here"] if top_model else "",
                "modelled_by": top_model.get("modelled_by", "") if top_model else "",
                "_rows": run,
            })
    return places


def annotate_place_calls(places, counted_records):
    """``modelled_in_run`` on each place: the first counted record (``counted_records`` in the
    report's order; HIGH or MEDIUM, not demoted) of the same genome and sequence that overlaps
    it, with its model id, label, home loci (``provenance``) and ``other_locus`` (the place's
    home locus is not among them; null when the place has no locus), else null. Record spans
    are 1-based and inclusive, place spans 0-based and half-open."""
    by_key = defaultdict(list)
    for rec in counted_records or ():
        by_key[(rec.get("genome"), rec.get("chrom"))].append(rec)
    for p in places:
        hit = next((rec for rec in by_key.get((p["genome"], p["chrom"]), ())
                    if int(rec.get("start", 0)) - 1 < p["end"] and p["start"] < int(rec.get("end", 0))), None)
        if hit is None:
            p["modelled_in_run"] = None
            continue
        provenance = list(hit.get("provenance") or [])
        p["modelled_in_run"] = {"mrna_id": hit.get("mrna_id", ""), "confidence": hit.get("confidence", ""),
                                "provenance": provenance,
                                "other_locus": (p["locus"] not in provenance) if p.get("locus") else None}
    return places


def summarize_nonsyntenic_candidates(records, per_genome, counted_records=None, dropped=None,
                                     dropped_hits=None):
    """Report block for §1x. `best_per_genome` is what a reader should look at: the
    strongest hit we found and declined to call, per genome, among the places no counted call
    of the run overlaps (``counted_records``: the HIGH and MEDIUM records that are not
    demoted). ``dropped``: the rows the rescue dropped (collect_nonsyntenic_candidates'
    dropped_out), listed under ``dropped`` with ``dropped_hits``, how many there were before the
    search's per-genome cap (None when no table says)."""
    by_gate = {}
    for r in records:
        by_gate[r["rejected_by"]] = by_gate.get(r["rejected_by"], 0) + 1
    places = annotate_place_calls(nonsyntenic_places(records), counted_records)
    unmodelled = [p for p in places if not p.get("modelled_in_run")]
    best = {}
    for p in unmodelled:
        top = max(p["_rows"], key=lambda r: r["identity"])
        old = best.get(p["genome"])
        if old is None or top["identity"] > old["identity"]:
            best[p["genome"]] = {
                "chrom": top["chrom"], "start": top["start"], "end": top["end"],
                "identity": top["identity"], "rejected_by": top["rejected_by"],
                "nearest_anchor_bp": top["nearest_anchor_bp"],
                "nearest_anchor_same_chrom": top["nearest_anchor_same_chrom"],
            }
    for p in places:
        p.pop("_rows", None)
    block = {
        "total": len(records),
        "genomes": sorted(per_genome),
        "by_gate": by_gate,
        "best_per_genome": {g: best[g] for g in sorted(best)},
        "note": NONSYNTENIC_REFUSAL_NOTE,
        "candidates": records,
        "places": places,
        "place_count": len(places),
        "unmodelled_place_count": len(unmodelled),
        "unmodelled_genomes": sorted({p["genome"] for p in unmodelled}),
    }
    if dropped is not None:
        by_test = {}
        for r in dropped:
            by_test[r["rejected_by"]] = by_test.get(r["rejected_by"], 0) + 1
        block["dropped"] = {
            "total": len(dropped),
            "hits": dropped_hits,
            "by_test": {k: by_test[k] for k in sorted(by_test)},
            "genomes": sorted({r["genome"] for r in dropped}),
            "note": NONSYNTENIC_DROPPED_NOTE,
            "candidates": dropped,
        }
    return block


def sum_target_policy_column(paths, column):
    """The sum of one column over every target-policy row (one per home locus and searched
    genome); None when no row carries it (a run before the column)."""
    total, seen = 0, False
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if column not in row:
                        continue
                    seen = True
                    total += int(float(row.get(column) or 0))
        except (OSError, csv.Error, ValueError):
            continue
    return total if seen else None


def genomes_without_flanking_hits(paths):
    """The genomes where a search found no hit of any flanking gene and still searched blocks
    (target-policy rows with n_flanking_genes 0 and blocks_searched > 0): their blocks were
    seeded by the GOI's own hits, with no neighbour test (RescueRule=no_flanking_hits)."""
    out = set()
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    try:
                        no_flanking = int(float(row.get("n_flanking_genes") or 0)) == 0
                        searched = int(float(row.get("blocks_searched") or 0)) > 0
                    except ValueError:
                        continue
                    if no_flanking and searched:
                        out.add(canonical_genome_id(
                            row.get("genome") or os.path.basename(path).split(TARGET_POLICY_SUFFIX)[0]))
        except (OSError, csv.Error):
            continue
    return sorted(out)


GENOME_WIDE_NOTE = ("Places of GOI hits outside every searched GOI zone and every written row, from "
                    "the genome-wide path (--goi_genome_wide_path). A place is tried when the home "
                    "protein's blastp E against its hits over the whole genome (-searchsp = |q| x "
                    "genome/3) is <= 1e-6, its HSPs at E <= 1 cover >= 0.7 of the query and no "
                    "other home protein is a real competitor within 10 bits; per home GOI the places "
                    "are tried best first until one gives an answer: a model the judge labels HIGH "
                    "in the genome place, with no exon across a home exon junction (no retrocopy; an "
                    "exon >= 30 nt per spanned junction longer than 3 x its aligned query span holds "
                    "an intron read as coding and does not count, also when the excess is one "
                    "read-through intron of >= 60 nt next to a lost one), not frame-broken without an "
                    "intron, overlapping no row. An answer is written as "
                    "MEDIUM, GOIClass=genome_wide_strong_hit, GenomeWide=1; it never seeds and is "
                    "never counted; the report keeps one per genome, only in a genome without a "
                    "counted call (summary.genome_wide_strong_hits; one inside a block's window is "
                    "near that block, summary.near_block_strong_hits). The genome place is calibrated "
                    "by decoys only in a --decoy_mode run (decoy_rate.genome_wide; decoy_places here).")


def collect_genome_wide_places(paths):
    """The places of the genome-wide path (hits/<genome>.genome_wide.tsv,
    iterative_search_runner.genome_wide_goi_pass): one record per place with its terms over the
    genome and what became of it (``decision``, and for a tried place its window's
    ``model_decision``); a decoy place (decoy mode) has its decoy kind in ``decoy``. Tables written
    before 2026-11-16 load with the columns they have. Unreadable files are skipped."""
    records = []
    for path in sorted(paths):
        try:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if not row.get("chrom"):
                        continue
                    records.append({
                        "genome": canonical_genome_id(
                            row.get("genome") or os.path.basename(path).split(GENOME_WIDE_SUFFIX)[0]),
                        "parent": row.get("parent", ""), "decoy": row.get("decoy", "") or "",
                        "chrom": row.get("chrom", ""),
                        "start": _safe_int(row.get("start")), "end": _safe_int(row.get("end")),
                        "hits": _safe_int(row.get("hits")),
                        "best_hit_identity": _safe_float(row.get("best_hit_identity")),
                        "evalue_genome": row.get("evalue_genome", ""),
                        "bits": row.get("bits", "") or "",
                        "coverage": row.get("coverage", "") or "",
                        "best_other": row.get("best_other", ""), "margin": row.get("margin", ""),
                        "real_competitor": row.get("real_competitor", ""),
                        "rank": row.get("rank", ""),
                        "decision": row.get("decision", ""),
                        "block": row.get("block", "") or "",
                        "model_decision": row.get("model_decision", "") or "",
                        "near_block": row.get("near_block", "") or "",
                        "significant_blocks": row.get("significant_blocks", "") or "",
                        "median_flank_identity": row.get("median_flank_identity", "") or "",
                        **({"cap": row["cap"]} if row.get("cap") else {}),
                        **({"own_copy_set_aside": row["own_copy_set_aside"],
                            "own_copy_bits": row.get("own_copy_bits", "")}
                           if row.get("own_copy_set_aside") else {}),
                    })
        except (OSError, csv.Error):
            continue
    records.sort(key=lambda r: (r["genome"], r["decoy"], r["parent"], r["chrom"], r["start"]))
    return records


def summarize_genome_wide_places(records):
    """Report block of the genome-wide path: the real places per decision, per window outcome
    (model_decision) and per genome, every place, and the decoy places (decoy mode) apart, with
    the note on the rule and its calibration."""
    def tally(rows):
        by_decision, by_model, per_genome = {}, {}, {}
        for r in rows:
            by_decision[r["decision"]] = by_decision.get(r["decision"], 0) + 1
            if r.get("model_decision"):
                by_model[r["model_decision"]] = by_model.get(r["model_decision"], 0) + 1
            g = per_genome.setdefault(r["genome"], {})
            g[r["decision"]] = g.get(r["decision"], 0) + 1
        return by_decision, by_model, per_genome
    real = [r for r in records if not r.get("decoy")]
    decoys = [r for r in records if r.get("decoy")]
    by_decision, by_model, per_genome = tally(real)
    block = {"total": len(real), "by_decision": by_decision, "by_model_decision": by_model,
             "per_genome": per_genome, "note": GENOME_WIDE_NOTE, "places": real}
    if decoys:
        d_decision, d_model, _per = tally(decoys)
        block["decoy_places"] = {"total": len(decoys), "by_decision": d_decision,
                                 "by_model_decision": d_model, "places": decoys}
    return block


def build_staging_diagnostics(results_dir, dir_patterns, match_counts):
    diagnostics = {
        "results_dir": results_dir,
        "results_dir_exists": os.path.isdir(results_dir),
        "dirs": {name: _dir_diagnostics(path, patterns) for name, (path, patterns) in dir_patterns.items()},
        "match_counts": dict(match_counts),
    }
    diagnostics["total_matches"] = sum(match_counts.values())
    diagnostics["empty"] = diagnostics["total_matches"] == 0
    return diagnostics


LOW_QUALITY_CAVEAT = ("A gene that is not found in such a genome may be missing from the "
                      "assembly, not from the species.")


def low_quality_warnings(failed_genomes, searched):
    """One sentence per target genome that failed the assembly quality check and was
    searched anyway (a genome the user chose is kept; see ``qc_fail_policy``).

    ``failed_genomes`` are the ``failed_genomes`` of ``summarize_qc``; ``searched`` the
    genomes with results downstream. Returns ``(genome_ids, sentences)``."""
    searched = set(searched)
    ids, sentences = [], []
    for rec in sorted(failed_genomes, key=lambda r: r.get("genome", "")):
        genome = rec.get("genome", "")
        if not genome or genome not in searched or genome in ids:
            continue
        ids.append(genome)
        reason = (rec.get("msg") or "").strip()
        sentences.append(
            f"LOW-QUALITY ASSEMBLY: {genome} failed the quality check"
            + (f" ({reason})" if reason else "")
            + f" and was searched anyway. {LOW_QUALITY_CAVEAT}")
    return ids, sentences


def _format_goi_headline(high, medium, low, n_genomes, ambiguous=0, curation_by_basis=None,
                         other_home_gene=0, genome_wide=0, near_block=0):
    """One-line summary of the GOI yield, for the summary block (docs/TODO.md §1k; names of
    the label rework, 2026-10-14):

        12 orthologs (9 high, 3 medium confidence); 4 need curation (3 gene order only,
        1 another home gene fits as well); 2 of another home gene; 1 strong hit(s) without
        synteny (not counted); 40 no significant match; across 33 genome(s).

    Only HIGH and MEDIUM are orthologs. AMBIGUOUS calls are named apart and never counted:
    ``gene_order_only`` is a weak match in the slot where the gene was expected,
    ``not_specific`` a significant match that another home gene explains as well.
    ``other_home_gene``: AMBIGUOUS calls that ownership gave to another home gene.
    ``genome_wide`` / ``near_block``: the answers of the genome-wide path that the report kept,
    outside every searched window / inside a block's window outside its GOI zone; never counted.
    ``n_genomes`` counts the genomes with a GOI row other than such an answer."""
    if not (high or medium or low or ambiguous or other_home_gene or genome_wide or near_block):
        return "No GOI ortholog annotations were produced."
    parts = []
    if high or medium:
        conf = []
        if high:
            conf.append(f"{high} high")
        if medium:
            conf.append(f"{medium} medium")
        parts.append(f"{high + medium} ortholog(s) ({', '.join(conf)} confidence)")
    else:
        parts.append("NO ortholog")
    if ambiguous:
        by = curation_by_basis or {}
        detail = []
        if by.get("gene_order_only"):
            detail.append(f"{by['gene_order_only']} gene order only")
        if by.get("not_specific"):
            detail.append(f"{by['not_specific']} another home gene fits as well")
        rest = ambiguous - sum(v for v in by.values() if v)
        if rest > 0 and detail:
            detail.append(f"{rest} by sequence and synteny alone")
        parts.append(f"{ambiguous} need curation" + (f" ({', '.join(detail)})" if detail else ""))
    if other_home_gene:
        parts.append(f"{other_home_gene} of another home gene")
    if genome_wide:
        parts.append(f"{genome_wide} strong hit(s) without synteny (not counted)")
    if near_block:
        parts.append(f"{near_block} strong hit(s) near a block, outside its GOI zone (not counted)")
    if low:
        parts.append(f"{low} no significant match")
    return "; ".join(parts) + f"; across {n_genomes} genome(s)."


def _confidence_rank(confidence):
    """HIGH > MEDIUM > AMBIGUOUS > LOW/unknown — used to keep the best member in dedup."""
    return {"HIGH": 3, "MEDIUM": 2, "AMBIGUOUS": 1}.get((confidence or "").upper(), 0)


def _record_rank(rec):
    """Order of records of one gene, best first (smaller is better): label, then the
    judge's order (synvoy_label.rank_key: E-value in the place with ties below 1e-10, the
    block's E_syn, bits, shorter span) for records that carry Eeff, else bits or identity
    (_model_strength) for records written before 2026-10-14 or by a rescue pass."""
    def num(key, default):
        try:
            return float(str(rec.get(key) or "").strip())
        except ValueError:
            return default
    eeff = num("eeff", None)
    if eeff is not None:
        judged = synvoy_label.rank_key(rec.get("confidence", ""), eeff, num("block_esyn", None),
                                       num("bits", 0.0), _span(rec), int(rec.get("start", 0) or 0))
        return (-_confidence_rank(rec.get("confidence")), 0) + judged[1:]
    strength = _model_strength(rec)
    return (-_confidence_rank(rec.get("confidence")), 1, -strength[0], -strength[1], _span(rec),
            int(rec.get("start", 0) or 0))


def _model_strength(rec):
    """How strongly a model matches the query, for choosing among models of one gene: its
    bit score against the home protein (``Bits=``, written since 2026-10-14) when it has
    one, else its identity (older runs). A model with a bit score ranks above one without,
    so the two units are never compared. On the Candida pilot, label then bits picked the
    model on the curated gene in 28 of 31 genome pairs with several models, label then
    identity in 16."""
    bits = str(rec.get("bits") or "").strip()
    if bits:
        return (1, _safe_float(bits))
    return (0, _safe_float(rec.get("identity")))


def collect_goi_annotations(gff_files):
    """Flatten every GOI mRNA across all per-locus region GFFs into structured
    records (genome, target chrom/start/end, confidence, identity, provenance).

    Backwards-compat wrapper around ``_process_gffs_unified`` — ``build_report`` calls
    the unified pass directly so the same GFFs aren't read three times at SP-family
    scale (docs/TODO.md §1c).
    """
    return _process_gffs_unified(gff_files)[1]


def collect_high_flanking_per_locus(gff_files):
    """Count distinct HIGH-confidence flanking genes per (genome, locus, chrom) from
    the staged per-locus GFFs. Used by build_self_consistency to spot blocks with
    strong conserved synteny but no GOI model (docs/TODO.md §1j identity-decay sanity).

    Backwards-compat wrapper around ``_process_gffs_unified`` (docs/TODO.md §1c).
    """
    return _process_gffs_unified(gff_files)[2]


def _reciprocal_overlap(a, b):
    """min(overlap/len_a, overlap/len_b) — 1.0 for identical spans, 0.0 if disjoint."""
    inter = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    if inter <= 0:
        return 0.0
    len_a = a["end"] - a["start"] + 1
    len_b = b["end"] - b["start"] + 1
    return min(inter / len_a, inter / len_b)


def _containment(a, b):
    """overlap / length of the SHORTER span — 1.0 when one span sits wholly inside the other."""
    inter = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    if inter <= 0:
        return 0.0
    shorter = min(a["end"] - a["start"] + 1, b["end"] - b["start"] + 1)
    return inter / shorter if shorter > 0 else 0.0


def _span(r):
    return max(1, r["end"] - r["start"] + 1)


def _prefers_contained_model(wide, compact, span_ratio, qcov_margin):
    """True when ``compact`` should represent the pair instead of the wider ``wide``.

    A model that spans several times more genome than a model of the same gene at the same
    place, while aligning barely more of the query, has chained in sequence that does not
    belong to the gene — a fused model. Family benchmark 2026-09-15: one serine-protease
    call was reported as an 11.4-kb hull rescue (query coverage 0.988) while the 2.6-kb
    in-block model of the same gene, matching the curated annotation exon for exon, covered
    0.904 and lost the tie on confidence; the extra 8.8 kb bought 8 points of coverage.

    Refuses to decide when either coverage is unrecorded: then the wider model may be the
    complete gene and the compact one a fragment of it.
    """
    if _span(wide) < span_ratio * _span(compact):
        return False
    try:
        wide_cov = float(str(wide.get("query_cov") or "").strip())
        compact_cov = float(str(compact.get("query_cov") or "").strip())
    except ValueError:
        return False
    return (wide_cov - compact_cov) < qcov_margin


def _not_specific_whole_model(rec):
    """AMBIGUOUS because another home protein explains the model as well (synvoy_label: Basis)."""
    return rec.get("confidence") == "AMBIGUOUS" and rec.get("basis") == "not_specific"


# What can overrule a piece as the whole gene: a gene model, or a hit the search itself left as one
# (the luciferase runs hold both). A hit chain or a tandem row can lie over several genes, and a hull
# rescue's not_specific can come from its place alone (it is judged in the genome place, with no block).
WHOLE_MODEL_EVIDENCE = {"exon_annotation", "rescued_exon"}


def _is_piece_of(piece, wide, min_overlap=0.8):
    """A real piece of ``wide``: the pipeline itself calls it a fragment (ModelStatus), it lies on the
    wide model's strand and shares coding bases with it, it is no near-identical twin (the dedup's
    reciprocal overlap), and the wide model covers at least as much of the query and is a gene model
    or a rescued hit. A complete gene inside another gene's intron, a hit chain over
    several copies, or a model and its rescue twin at the same coordinates is no case of this."""
    if piece.get("model_status") != "fragment":
        return False
    if wide.get("evidence") and wide["evidence"] not in WHOLE_MODEL_EVIDENCE:
        return False
    if _reciprocal_overlap(piece, wide) > min_overlap:
        return False
    strands = {piece.get("strand"), wide.get("strand")}
    if strands <= {"+", "-"} and len(strands) == 2:
        return False
    if piece.get("cds") and wide.get("cds") and not any(
            min(e1, e2) >= max(s1, s2) for s1, e1 in piece["cds"] for s2, e2 in wide["cds"]):
        return False
    return _safe_float(wide.get("query_cov"), 1.0) >= _safe_float(piece.get("query_cov"), 0.0)


def _same_orf(a, b, min_bp=3):
    """Two gene models that are the two halves of one reading frame, miniprot's models of one gene
    (fix of 9 Nov, item D2: C. lusitaniae orf19.968, a HIGH and a MEDIUM model sharing 30 bp of one
    reading frame of an expanded repeat, counted twice in both parent arms). Both exon_annotation,
    on one strand, and meeting only where the upstream model's last CDS row and the downstream
    model's first CDS row overlap, partly (neither row inside the other), by at least ``min_bp``
    bases read in one frame by both: the shared bases encode the same residues in both models. A
    model with a spurious terminal exon inside its neighbour's exon (a tandem copy) stays apart; a
    stop between the two halves is not looked for (the report has no genome)."""
    if a.get("evidence") != "exon_annotation" or b.get("evidence") != "exon_annotation":
        return False
    if a.get("strand") != b.get("strand") or a.get("strand") not in ("+", "-"):
        return False
    fa, fb = a.get("cds_frames") or [], b.get("cds_frames") or []
    if not fa or len(fa) != len(a.get("cds") or []) or not fb or len(fb) != len(b.get("cds") or []):
        return False
    up, down = sorted((sorted(zip(a["cds"], fa)), sorted(zip(b["cds"], fb))), key=lambda rows: rows[0][0])
    (us, ue), uf = max(up, key=lambda r: r[0][1])
    (ds, de), df = down[0]
    if not (us < ds and ue < de) or ue - ds + 1 < min_bp or uf != df:
        return False
    return (all(e < ds for (s, e), _f in up if (s, e) != (us, ue))
            and all(s > ue for (s, e), _f in down if (s, e) != (ds, de)))


def _whole_model_of(wide, members, piece, min_overlap=0.8):
    """The whole gene a piece is tested against: ``wide`` itself when it is a gene model or a
    rescued hit, else (a hull rescue, a hit chain, a tandem row) its not-specific twin of that kind
    among the cluster's members (the same gene, each covering more than ``min_overlap`` of the
    other): a hull rescue that outranked its exon twin left c7's piece counted (fix item 44).
    ``wide`` when there is none."""
    if wide.get("evidence") in WHOLE_MODEL_EVIDENCE or not wide.get("evidence"):
        return wide
    return next((m for m in members if m is not piece and m is not wide
                 and m.get("evidence") in WHOLE_MODEL_EVIDENCE and _not_specific_whole_model(m)
                 and _reciprocal_overlap(m, wide) > min_overlap), wide)


def dedupe_goi_annotations(annotations, min_overlap=0.8, contain_min=0.9,
                           fused_span_ratio=2.0, fused_qcov_margin=0.15):
    """Collapse HIGH/MEDIUM/AMBIGUOUS GOI hits that are the same target gene found from
    multiple home loci (docs/TODO.md §1h).

    Key: ``(genome, target_chrom)`` plus reciprocal coordinate overlap > ``min_overlap``.
    The kept record is the highest-confidence member, then the one with the highest bit
    score against the home protein (identity for models without one), with a
    ``provenance`` list of the source home loci it was recovered from. When the
    identity is bit-for-bit identical across >=2 loci (a near-certain "one gene,
    many seeds" signal) the record is reclassified ``goi_class=cross_locus_duplicate``.

    LOW-confidence hits are left untouched (they are noisy fallback spans, not
    ortholog calls), matching the headline metric which only counts HIGH/MEDIUM.
    AMBIGUOUS calls ARE deduped: they are reported (separately from the ortholog
    count), so the same candidate must not appear once per home locus.
    """
    considered = [a for a in annotations
                  if a["confidence"] in {"HIGH", "MEDIUM", "AMBIGUOUS"}]
    groups = defaultdict(list)
    for a in considered:
        groups[(a["genome"], a["chrom"])].append(a)

    records = []
    for (genome, chrom), items in groups.items():
        items.sort(key=lambda r: (r["start"], r["end"]))
        clusters = []  # each: {"members": [...], "rep": <best member>}
        for a in items:
            placed = False
            for cl in clusters:
                rep = cl["rep"]
                # Same target gene: near-identical spans (one gene found from several home
                # loci), or one span wholly inside the other (a fused or rescue model that
                # swallowed the compact model of the same gene).
                contained = _containment(a, rep) >= contain_min
                same_orf = _same_orf(a, rep)
                if _reciprocal_overlap(a, rep) > min_overlap or contained or same_orf:
                    if same_orf and not contained and _reciprocal_overlap(a, rep) <= min_overlap:
                        cl["same_orf"] = True
                    cl["members"].append(a)
                    wide, compact = (a, rep) if _span(a) >= _span(rep) else (rep, a)
                    # The compact model must itself be a counted call (MEDIUM or better):
                    # an AMBIGUOUS hit chain inside a HIGH model of the same gene took the
                    # gene out of the headline (bug hunt 2026-10-13, A2; never seen in runs).
                    counted_piece = (contained and wide is not compact
                                     and _confidence_rank(compact["confidence"]) >= _confidence_rank("MEDIUM"))
                    if counted_piece and _prefers_contained_model(
                            wide, compact, fused_span_ratio, fused_qcov_margin):
                        if cl["rep"] is not compact:
                            compact["superseded_wide_model"] = True
                            cl["rep"] = compact
                    elif counted_piece and _not_specific_whole_model(wide) and _is_piece_of(
                            compact, _whole_model_of(wide, cl["members"], compact, min_overlap), min_overlap):
                        # A piece of a gene another home protein explains as well is no ortholog
                        # call: the whole gene's verdict stands. Without this rule a paralog's
                        # exon, HIGH or MEDIUM on its own, won on confidence (the luciferase runs).
                        # Counted models merged under the piece earlier are judged against the
                        # gene themselves.
                        if cl["rep"] is compact:
                            for m in cl["members"]:
                                if m is not compact and m is not wide and _confidence_rank(m["confidence"]) >= _confidence_rank("MEDIUM"):
                                    items.append(m)
                            cl["members"] = [m for m in cl["members"] if m is compact or m is wide
                                             or _confidence_rank(m["confidence"]) < _confidence_rank("MEDIUM")]
                        cl.setdefault("overruled", []).append(compact)
                        cl["rep"] = wide
                    else:
                        if _record_rank(a) < _record_rank(cl["rep"]):
                            cl["rep"] = a
                    # A piece was overruled by the record's whole gene. When that gene stops being
                    # the record (the fused-model rule, or the judge's order chose another model),
                    # the piece goes back into the pool: it must not vanish into a record that says
                    # nothing about it.
                    for piece in list(cl.get("overruled", ())):
                        rep = cl["rep"]
                        if not (_not_specific_whole_model(rep)
                                and _is_piece_of(piece, _whole_model_of(rep, cl["members"], piece, min_overlap),
                                                 min_overlap)
                                and _containment(piece, rep) >= contain_min):
                            cl["overruled"] = [p for p in cl["overruled"] if p is not piece]
                            cl["members"] = [m for m in cl["members"] if m is not piece]
                            items.append(piece)
                    placed = True
                    break
            if not placed:
                clusters.append({"members": [a], "rep": a})

        for cl in clusters:
            members, rep = cl["members"], cl["rep"]
            sources = sorted({m["source"] for m in members})
            identity_strs = {m["identity"] for m in members}
            is_cross_dup = (
                len(sources) >= 2
                and len(identity_strs) == 1
                and "" not in identity_strs
            )
            records.append({
                "genome": genome,
                "chrom": chrom,
                "start": rep["start"],
                "end": rep["end"],
                "confidence": rep["confidence"],
                "identity": rep["identity"],
                # Bit score against the home protein ("" for a model written before
                # 2026-10-14); the member that ranks first (_record_rank) is kept.
                "bits": rep.get("bits", ""),
                "basis": rep.get("basis", ""),
                "eeff": rep.get("eeff", ""),
                "place": rep.get("place", ""),
                "block_esyn": rep.get("block_esyn", ""),
                "goi_class": "cross_locus_duplicate" if is_cross_dup else rep["goi_class"],
                "original_goi_class": rep["goi_class"],
                "target_gene": rep["target_gene"],
                "model_status": rep["model_status"],
                # Only on a hit chain the fallback gate refused and the search wrote because its
                # block is significant (RefusalLifted): the tests it failed.
                **({"refusal_lifted": rep["refusal_lifted"]} if rep.get("refusal_lifted") else {}),
                "provenance": sources,
                "n_source_loci": len(sources),
                "n_merged_hits": len(members),
                "cross_locus_duplicate": is_cross_dup,
                # Carried so the coverage demotion can distinguish "coverage is low"
                # from "coverage was never recorded" (see apply_coverage_demotion).
                "query_cov": rep.get("query_cov", ""),
                # A model's core coverage (core_completeness), for the same demotion.
                **({"core_cov": rep["core_cov"]} if rep.get("core_cov") else {}),
                # The own copy of q that the specificity left out of p*, and by how many bits
                # the model favours it over q (fix-phase d2); absent when there is none.
                **({"own_copy_set_aside": rep["own_copy_set_aside"],
                    "own_copy_bits": rep.get("own_copy_bits", "")}
                   if rep.get("own_copy_set_aside") else {}),
                # Set when this compact model replaced a wider model of the same gene that
                # aligned barely more of the query (see _prefers_contained_model).
                "superseded_wide_model": bool(rep.get("superseded_wide_model")),
                # Set when this not-specific whole gene overruled a counted piece of itself.
                **({"overruled_counted_piece": True} if cl.get("overruled") else {}),
                # Set when two models of one open reading frame were merged (_same_orf).
                **({"same_orf_merged": True} if cl.get("same_orf") else {}),
                # Join key for the F9 phylo-placement verdicts, which are computed from
                # goi_for_tree.faa and so are keyed by model id, not coordinates.
                "mrna_id": rep.get("mrna_id", ""),
                # What the record is made of, for the readers of one home locus's model ids (F9's
                # verdicts, the tree's leaves; fix item 46). Model ids restart per home locus, so a
                # member is its (locus, model id, start, end). Private: build_report drops them.
                "_rep": (rep["source"], rep.get("mrna_id", "")),
                "_members": [(m["source"], m.get("mrna_id", ""), m["start"], m["end"]) for m in members],
                "_overruled": [(m["source"], m.get("mrna_id", "")) for m in cl.get("overruled", ())],
                # Only in a run with --report_copies_in_slot (preset plan PR-6).
                **({"copies_in_slot": rep["copies_in_slot"], "slot_rank": rep.get("slot_rank", "")}
                   if rep.get("copies_in_slot") else {}),
                # Only on a model of a tandem-copy lineage: the copy it was found with, and its
                # identity, coverage and bits against the home protein.
                **({"copy_lineage": rep["copy_lineage"], "home_identity": rep.get("home_identity", ""),
                    "home_query_cov": rep.get("home_query_cov", ""), "home_bits": rep.get("home_bits", "")}
                   if rep.get("copy_lineage") else {}),
                # Only on a call the dispersed-GOI rescue alone made the search find: every
                # model merged into the record carries DispersedRescue=1. One home locus's rescue
                # and another's own block finding the same gene is no rescue-born call (which
                # member was the representative depended on the order of the loci).
                **({"dispersed_rescue": True, "rescue_rule": rep.get("rescue_rule", ""),
                    "rescue_anchors": rep.get("rescue_anchors", ""),
                    "rescue_anchor_bp": rep.get("rescue_anchor_bp", "")}
                   if all(m.get("dispersed_rescue") for m in members) else {}),
                # Only on a call of the genome-wide path (GenomeWide=1): frame_broken says that its
                # model has no clean frame (an answer then has an intron).
                **({"genome_wide": True,
                    "genome_wide_coverage": rep.get("genome_wide_coverage", ""),
                    "retrocopy_candidate": rep.get("retrocopy_candidate", ""),
                    "neighbour_median_identity": rep.get("neighbour_median_identity", ""),
                    "near_block": rep.get("near_block", ""),
                    "frame_broken": bool(rep.get("frame_broken"))}
                   if rep.get("genome_wide") else {}),
            })

    records.sort(key=lambda r: (_record_rank(r), r["genome"], r["chrom"]))
    return {
        "records": records,
        "high_confidence_goi_deduped": sum(1 for r in records if r["confidence"] == "HIGH"),
        "medium_confidence_goi_deduped": sum(1 for r in records if r["confidence"] == "MEDIUM"),
        "ambiguous_goi_deduped": sum(1 for r in records if r["confidence"] == "AMBIGUOUS"),
        "cross_locus_duplicates": sum(1 for r in records if r["cross_locus_duplicate"]),
        "hits_collapsed_by_dedup": len(considered) - len(records),
        "wide_models_superseded": sum(1 for r in records if r["superseded_wide_model"]),
        "pre_dedup_high_medium": len(considered),
        "post_dedup_records": len(records),
    }


def _is_counted(r, demoted_classes):
    """A counted call: HIGH or MEDIUM, of no demoted class (ownership, the coverage rule), and no
    answer of the genome-wide path (``genome_wide_strong_hit``, a reading of its own; the field
    says so, not the class, which the cross-locus relabel can overwrite)."""
    return (r.get("confidence") in ("HIGH", "MEDIUM") and r.get("goi_class") not in demoted_classes
            and not r.get("genome_wide_strong_hit"))


def counted_genomes(records, demoted_classes):
    """The genomes that hold a counted call (_is_counted) among ``records``."""
    return {r.get("genome") for r in records or () if _is_counted(r, demoted_classes)}


def dedupe_goi_records(annotations):
    """dedupe_goi_annotations of the GOI rows, the answers of the genome-wide path
    (GOIClass=genome_wide_strong_hit) apart: they are deduped among themselves, marked
    ``genome_wide_strong_hit`` and appended after the other records, so that phylo, ownership and
    the coverage demotion see them; every count of the result is that of the other records, which
    are exactly as without the answers."""
    answers = [a for a in annotations if a.get("goi_class") == synvoy_label.GENOME_WIDE_CLASS]
    goi_dedup = dedupe_goi_annotations([a for a in annotations
                                        if a.get("goi_class") != synvoy_label.GENOME_WIDE_CLASS])
    if answers:
        extra = dedupe_goi_annotations(answers)["records"]
        for r in extra:
            r["genome_wide_strong_hit"] = True
        goi_dedup["records"].extend(extra)
    return goi_dedup


def genome_wide_answer_records(records, demoted_classes):
    """The answers of the genome-wide path that stand: ``genome_wide_strong_hit`` records not
    demoted by ownership (a demoted one stays a record of its class and is no answer)."""
    return [r for r in records or () if r.get("genome_wide_strong_hit")
            and r.get("goi_class") not in demoted_classes]


def drop_genome_wide_strong_hits(goi_dedup, demoted_classes):
    """The report's rule for the answers of the genome-wide path (genome_wide_answer_records): an
    answer is dropped, with its reason, in a genome that holds a counted call
    (``genome_has_counted_call``, _is_counted), where it overlaps (same sequence, any base) a HIGH,
    MEDIUM or AMBIGUOUS record of any other source, another home locus or a rescue pass
    (``overlaps_record``; a LOW record holds a poorer model of the same DNA and is no record here),
    and where another answer of the genome ranks first (``not_best_in_genome``, _record_rank, then
    sequence and start): at most one answer per genome, and only where nothing is counted. Runs
    after ownership and the coverage demotion, before drop_gene_order_only. The dropped records are
    listed (``genome_wide_strong_hits_dropped_records``); apply_report_drops.py leaves their rows
    out of the figures. Mutates ``goi_dedup``; returns the dropped records."""
    records = goi_dedup.get("records", [])
    counted = counted_genomes(records, demoted_classes)
    others = defaultdict(list)
    for r in records:
        if not r.get("genome_wide_strong_hit") and r.get("confidence") in ("HIGH", "MEDIUM", "AMBIGUOUS"):
            others[(r["genome"], r["chrom"])].append(r)
    answers = genome_wide_answer_records(records, demoted_classes)
    reasons = {}
    for r in answers:
        if r["genome"] in counted:
            reasons[id(r)] = "genome_has_counted_call"
        elif any(int(o["start"]) <= int(r["end"]) and int(r["start"]) <= int(o["end"])
                 for o in others.get((r["genome"], r["chrom"]), ())):
            reasons[id(r)] = "overlaps_record"
    standing = defaultdict(list)
    for r in answers:
        if id(r) not in reasons:
            standing[r["genome"]].append(r)
    for genome in sorted(standing):
        rs = standing[genome]
        best = min(rs, key=lambda r: (_record_rank(r), r["chrom"], int(r["start"]), r.get("mrna_id", "")))
        for r in rs:
            if r is not best:
                reasons[id(r)] = "not_best_in_genome"
    dropped = [r for r in answers if id(r) in reasons]
    goi_dedup["records"] = [r for r in records if id(r) not in reasons]
    goi_dedup["genome_wide_strong_hits_dropped"] = len(dropped)
    goi_dedup["genome_wide_strong_hits_dropped_records"] = [
        {"genome": r["genome"], "chrom": r["chrom"], "start": r["start"], "end": r["end"],
         "mrna_id": r.get("mrna_id", ""), "provenance": r.get("provenance", []),
         "reason": reasons[id(r)]}
        for r in dropped]
    return dropped


def drop_gene_order_only(goi_dedup, demoted_classes):
    """A weak match where the gene should be (Basis=gene_order_only) is offered for curation
    only in a genome without a counted call of the gene from any home locus or rescue pass
    (plan 3.4): where the gene was found, it is a second, weaker candidate. Counted = HIGH or
    MEDIUM and not demoted by the ownership or the coverage rule, so this runs after them
    (review fix RF-3: it ran inside the dedup, and a genome whose only call was then demoted
    had lost its candidate as well); an answer of the genome-wide path is not counted
    (_is_counted). The dropped records are listed
    (``gene_order_only_dropped_records``): the per-genome GFFs and the figures drop the same
    rows (bin/apply_report_drops.py). Mutates ``goi_dedup``; returns the dropped records."""
    records = goi_dedup.get("records", [])
    counted = counted_genomes(records, demoted_classes)
    dropped = [r for r in records if r.get("basis") == "gene_order_only" and r["genome"] in counted]
    keep = {id(r) for r in records} - {id(r) for r in dropped}
    goi_dedup["records"] = [r for r in records if id(r) in keep]
    goi_dedup["gene_order_only_dropped"] = len(dropped)
    goi_dedup["gene_order_only_dropped_records"] = [
        {"genome": r["genome"], "chrom": r["chrom"], "start": r["start"], "end": r["end"],
         "mrna_id": r.get("mrna_id", ""), "provenance": r.get("provenance", [])}
        for r in dropped]
    main = [r for r in goi_dedup["records"] if not r.get("genome_wide_strong_hit")]
    goi_dedup["ambiguous_goi_deduped"] = sum(1 for r in main if r.get("confidence") == "AMBIGUOUS")
    goi_dedup["post_dedup_records"] = len(main)
    return dropped


GOI_LOW_RECORDS_NOTE = ("LOW models: not orthology calls, counted nowhere in this report; one entry per "
                        "distinct model over all home loci and rescue passes. basis '' = not judged "
                        "(the strong-synteny rescue's models, rows the judge could not label).")


def goi_low_records(goi_annotations, records, demoted_classes, ownership_rows=()):
    """The LOW GOI models of the run, one entry per distinct model: rows of the same genome,
    sequence, start, end and strand (one model found from several home loci or rescue passes)
    are one entry; its representative is the row with CDS rows, else the first by (source, model
    id), and ``provenance`` lists the sources. Each entry names the best record of ``records``
    (the report's records after ownership, the coverage rule and the drops) on the same sequence
    that overlaps it (``overlapping_record``, null when none), whether its genome holds a counted
    call (``genome_has_counted_call``) and, when an ownership row (reciprocal overlap above 0.5,
    best by bit score, then overlap) covers it, the ownership check's best panel entry, its score
    and the gap to the second (raw Smith-Waterman scores). No count, flag or label reads them."""
    groups = defaultdict(list)
    for a in goi_annotations or ():
        if a.get("confidence") == "LOW":
            groups[(a["genome"], a["chrom"], a["start"], a["end"], a.get("strand", ""))].append(a)
    counted = counted_genomes(records, demoted_classes)
    recs_by_gc = defaultdict(list)
    for r in records or ():
        recs_by_gc[(r.get("genome"), r.get("chrom"))].append(r)
    own_by_gc = defaultdict(list)
    for row in ownership_rows or ():
        own_by_gc[(canonical_genome_id(row.get("genome", "")), row.get("chrom", ""))].append(row)
    out = []
    for key in sorted(groups):
        genome, chrom, start, end, strand = key
        members = groups[key]
        rep = min(members, key=lambda a: (0 if a.get("cds") else 1, a.get("source", ""), a.get("mrna_id", "")))
        sources = sorted({m.get("source", "") for m in members})
        over = [r for r in recs_by_gc.get((genome, chrom), ())
                if int(r["start"]) <= end and start <= int(r["end"])]
        best = min(over, key=lambda r: (_record_rank(r), int(r["start"]), r.get("mrna_id", ""))) if over else None
        entry = {
            "genome": genome, "chrom": chrom, "start": start, "end": end, "strand": strand,
            "confidence": "LOW",
            "basis": rep.get("basis", ""),
            "eeff": rep.get("eeff", ""),
            "place": rep.get("place", ""),
            "block_esyn": rep.get("block_esyn", ""),
            "bits": rep.get("bits", ""),
            "identity": rep.get("identity", ""),
            "query_cov": rep.get("query_cov", ""),
            **({"core_cov": rep["core_cov"]} if rep.get("core_cov") else {}),
            "evidence": rep.get("evidence", ""),
            "goi_class": rep.get("goi_class", ""),
            "model_status": rep.get("model_status", ""),
            "mrna_id": rep.get("mrna_id", ""),
            "cds": [[cs, ce] for cs, ce in sorted(rep.get("cds") or ())],
            "provenance": sources,
            "n_source_loci": len(sources),
            "target_gene": rep.get("target_gene", ""),
            "overlapping_record": ({"confidence": best.get("confidence", ""), "goi_class": best.get("goi_class", ""),
                                    "start": best["start"], "end": best["end"],
                                    "counted": _is_counted(best, demoted_classes)}
                                   if best is not None else None),
            "genome_has_counted_call": genome in counted,
        }
        scored = []
        for row in own_by_gc.get((genome, chrom), ()):
            ov = _reciprocal_overlap({"start": start, "end": end}, {"start": row["start_i"], "end": row["end_i"]})
            if ov > 0.5:
                scored.append((row["best_bit_f"], ov, row))
        if scored:
            row = max(scored, key=lambda t: (t[0], t[1]))[2]
            entry["ownership_best_paralog"] = row.get("best_paralog", "")
            entry["ownership_bit"] = row["best_bit_f"]
            entry["ownership_gap"] = row["gap_f"]
        out.append(entry)
    return {"note": GOI_LOW_RECORDS_NOTE, "count": len(out), "records": out}


def _panel_scores(field):
    """``name=score;name=score`` of reciprocal_best_paralog_check.py as {name: score}. A name may
    hold ';' or '=' (a GFF id with %3B): the number ends each item."""
    return {name: float(score) for name, score in
            re.findall(r"(.+?)=(-?\d+(?:\.\d+)?)(?:;|$)", field or "")}


def _is_goi_for_call(label, scores, panel_meta):
    """A panel entry counts as the GOI for one call: the GOI's own entry, or an own copy (found_by
    own_copy: near-identical to the GOI, build_home_paralog_panel) that the call does not favour
    over the GOI's own entry by the bar at which the search's p* tells them apart
    (synvoy_specificity.SEED_MIN_MARGIN_BITS). Scores are raw Smith-Waterman scores."""
    meta = panel_meta.get(label, {})
    if not meta.get("is_goi"):
        return False
    if meta.get("goi_found_by") != "own_copy":
        return True
    own = [bit for name, bit in scores.items() if panel_meta.get(name, {}).get("is_goi")
           and panel_meta.get(name, {}).get("goi_found_by") != "own_copy"]
    if not own or label not in scores:
        return True
    return (synvoy_specificity.raw_to_bits(scores[label]) - synvoy_specificity.raw_to_bits(max(own))
            < synvoy_specificity.SEED_MIN_MARGIN_BITS)


def _ownership_scores(row):
    """Every panel entry's score in an ownership row. A table written before panel_scores names only
    the two best entries."""
    return row.get("panel_scores_d") or {
        label: bit for label, bit in ((row.get("best_paralog"), row.get("best_bit_f")),
                                      (row.get("second_paralog"), row.get("second_bit_f"))) if label}


def load_paralog_check_rows(paralog_check_dir, suffix=".locus_ownership.tsv"):
    """Read every ``<locus>__<genome><suffix>`` TSV of reciprocal_best_paralog_check.py: the §1m
    locus-ownership tables of ASSIGN_LOCUS_OWNERSHIP (locus_ownership/), whose ``best_paralog``
    names a ``locus_<id>|<gene>`` entry of the home-paralog panel.

    Returns a list of dicts keyed by the TSV header. Skips header-only files and silently
    tolerates missing dirs."""
    rows = []
    if not paralog_check_dir or not os.path.isdir(paralog_check_dir):
        return rows
    with os.scandir(paralog_check_dir) as it:
        for entry in sorted(it, key=lambda e: e.name):
            if not entry.is_file() or not entry.name.endswith(suffix):
                continue
            try:
                with open(entry.path) as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        if not row.get("best_paralog"):
                            continue
                        try:
                            row["best_bit_f"] = float(row.get("best_bit") or 0)
                            row["second_bit_f"] = float(row.get("second_bit") or 0)
                            row["gap_f"] = float(row.get("bitscore_gap") or 0)
                            row["start_i"] = int(row.get("start") or 0)
                            row["end_i"] = int(row.get("end") or 0)
                            row["panel_scores_d"] = _panel_scores(row.get("panel_scores"))
                        except ValueError:
                            continue
                        rows.append(row)
            except OSError:
                continue
    return rows


# A step that failed and was recorded by bin/run_step.py (fix-phase d5): one line in its output.
STEP_FAILED_MARKER = "# SYNVOY_STEP_FAILED"
_STEP_FAILURE_SENTENCE = {
    "ASSIGN_LOCUS_OWNERSHIP": "its calls were not checked against the other home genes, so a call there may be "
                              "another home gene's (ownership_not_checked)",
    "BUILD_HOME_PARALOG_PANEL": "no call was checked against the other home genes, so any call may be another home "
                                "gene's (ownership_not_checked)",
    "RESCUE_GOI_HULL": "a gene only the hull rescue would model is missing there",
    "RESCUE_STRONG_SYNTENY": "a gene only the strong-synteny rescue would model is missing there",
}


def collect_step_failures(results_dir):
    """The steps that failed and were recorded instead of being ignored (bin/run_step.py,
    fix-phase d5): the marker line in a staged region GFF (the rescues), ownership table or the
    panel's meta table. Sorted list of {step, locus, genome, exit, reason}."""
    found = set()
    for sub in ("regions", "locus_ownership"):
        folder = os.path.join(results_dir, sub)
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            path = os.path.join(folder, name)
            if not os.path.isfile(path):
                continue
            try:
                with open(path, errors="replace") as fh:
                    head = [fh.readline() for _ in range(3)]
            except OSError:
                continue
            for line in head:
                if not line.startswith(STEP_FAILED_MARKER):
                    continue
                f = dict(kv.split("=", 1) for kv in line.rstrip("\n").split("\t")[1:] if "=" in kv)
                genome = canonical_genome_id(f["genome"]) if f.get("genome") else ""
                # The ownership task of a hull rescue's models runs under <genome>.hull_rescue: it
                # checked only those models.
                part = "hull_rescue" if (f.get("genome") or "").endswith(".hull_rescue") else ""
                found.add((f.get("step", ""), f.get("locus", ""), genome, f.get("exit", ""), f.get("reason", ""), part))
    out = []
    for t in sorted(found):
        row = dict(zip(("step", "locus", "genome", "exit", "reason"), t[:5]))
        if t[5]:
            row["only"] = t[5]
        out.append(row)
    return out


def mark_ownership_not_checked(goi_dedup, failures):
    """``ownership_not_checked: true`` on every record a failed ownership step would have
    checked: the records of that genome recovered from that locus, or all of them when the
    panel failed. Returns how many records were marked."""
    pairs = [(f["locus"], f["genome"], f.get("only")) for f in failures if f["step"] == "ASSIGN_LOCUS_OWNERSHIP"]
    panel_failed = any(f["step"] == "BUILD_HOME_PARALOG_PANEL" for f in failures)
    n = 0
    for rec in goi_dedup.get("records", []):
        hull = "synteny_hull_rescue" in (rec.get("goi_class"), rec.get("original_goi_class"))
        if panel_failed or any(rec.get("genome") == g and (not l or l in (rec.get("provenance") or []))
                               and (only != "hull_rescue" or hull)
                               for l, g, only in pairs):
            rec["ownership_not_checked"] = True
            n += 1
    return n


def step_failure_warnings(failures):
    """One sentence per recorded failure: the step, the genome and locus, the exit and its last
    stderr line, and what is missing from the report because of it."""
    out = []
    for f in failures:
        where = ", ".join(x for x in (f["genome"], f["locus"]) if x)
        out.append(f"{f['step']} failed" + (f" for {where}" if where else "")
                   + f" (exit {f['exit']}" + (f": {f['reason']}" if f["reason"] else "") + "): "
                   + _STEP_FAILURE_SENTENCE.get(f["step"], "its result is missing") + " (failed_steps).")
    return out


def count_home_loci(gff_files):
    """The home loci that produced region files: their staging tags (``locus_N__``). Files without
    a tag (an older run, a direct call) say nothing about loci, so they add none (fix item 38: the
    advisory counted (genome, staging source, chromosome) cells with a HIGH flanking gene, which
    gave '9 home loci' for a gene with one)."""
    tags = set()
    for p in gff_files or ():
        m = _STAGING_PREFIX_RE.match(os.path.basename(p))
        if m and m.group(1).startswith("locus_"):
            tags.add(m.group(1))
    return len(tags)


def explained_by_other_home_gene(rec):
    """A record that another home gene explains at least as well: the paralog check gave it to
    another home gene (``paralog_not_goi``), or the judge found another home protein within the
    specificity margin (AMBIGUOUS, Basis=not_specific)."""
    return (rec.get("goi_class") == "paralog_not_goi"
            or (rec.get("confidence") == "AMBIGUOUS" and rec.get("basis") == "not_specific"))


def build_gene_family_advisory(records, n_home_loci=0,
                               family_ratio=0.4, family_min_records=5,
                               family_min_loci=3):
    """Detect a gene-family / multi-paralog query and advise honest interpretation.
    Purely data-driven from signals the report already has — no search-side change.
    SynVoy is validated for divergent LOW-COPY genes via conserved flanking synteny; inside
    a large paralog family ortholog-vs-paralog assignment is unreliable, so we say so
    instead of asserting the GOI name.

    ``records``: the GOI records after the dedup, ownership, the coverage rule and the drops
    (the answers of the genome-wide path are left out). Triggers (either):
      - at least ``family_min_records`` records, and at least ``family_ratio`` of them, are
        explained by another home gene at least as well (explained_by_other_home_gene: the
        paralog check, or specificity), or
      - the query fanned out to many home loci (paralog cluster in the home genome).
    LOW rows are no records and do not count: a LOW model is a hit without a significant
    match, which says nothing about a family.
    """
    recs = [r for r in records or () if not r.get("genome_wide_strong_hit")]
    other = sum(1 for r in recs if explained_by_other_home_gene(r))
    total = len(recs)
    reasons = []
    if total and other >= family_min_records and other >= family_ratio * total:
        pct = round(100 * other / total)
        reasons.append(
            f"{other} of {total} GOI records ({pct}%) are explained by another home gene at "
            "least as well (paralog check or specificity)")
    if n_home_loci >= family_min_loci:
        reasons.append(f"the query resolved to {n_home_loci} home loci (paralog cluster)")
    is_family = bool(reasons)
    return {
        "is_gene_family_query": is_family,
        "reasons": reasons,
        "advice": (
            "Gene-family / multi-paralog query detected. SynVoy is validated for "
            "divergent LOW-COPY genes via conserved flanking synteny; within a large "
            "paralog family, ortholog-vs-paralog assignment is unreliable. Treat the GOI "
            "calls as candidates that need curation, and read each record's locus "
            "ownership (goi_class paralog_not_goi, owning_gene) and the report's "
            "self_consistency flags." if is_family else ""
        ),
    }


def load_panel_meta(panel_meta_path):
    """Read ``home_paralog_panel.faa.meta.tsv`` (build_home_paralog_panel.py, §1m).

    Returns {panel_id: {locus, gene, chrom, is_goi, goi_found_by}} so the report can
    map an ownership row's ``best_paralog`` (a ``locus_<id>|<gene>`` panel label)
    back to its home locus / gene / GOI-ness. ``goi_found_by`` (GOI entry only, since
    2026-10-28): how the GOI entry was found, 'goi_info', 'identity_rule' or 'query';
    '' in a meta file written before. Tolerates a missing file (ownership becomes
    label-only / advisory)."""
    meta = {}
    if not panel_meta_path or not os.path.isfile(panel_meta_path):
        return meta
    try:
        with open(panel_meta_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                pid = row.get("panel_id")
                if not pid:
                    continue
                meta[pid] = {
                    "locus": row.get("locus_id", ""),
                    "gene": row.get("gene", ""),
                    "chrom": row.get("chrom", ""),
                    "is_goi": _is_true(row.get("is_goi_gene")),
                    "goi_found_by": row.get("goi_found_by") or "",
                }
    except OSError:
        pass
    return meta


def goi_panel_entry(panel_meta):
    """The panel's GOI entry as the report records it (``locus_ownership_summary.goi_entry``):
    {gene, panel_ids, found_by}; found_by is 'goi_info' (the home gene the home step took),
    'identity_rule' (the home protein near-identical to the query), 'query' (no home gene is
    the query's own, so the query itself), or 'unrecorded' (a panel built before
    2026-10-28). None without a panel or without a GOI entry."""
    goi = sorted(pid for pid, m in (panel_meta or {}).items() if m.get("is_goi"))
    if not goi:
        return None
    # Own copies (found_by own_copy: home copies near-identical to the query) are listed in
    # panel_ids; the entry itself is the gene the rule above names.
    first = panel_meta[next((pid for pid in goi if panel_meta[pid].get("goi_found_by") != "own_copy"), goi[0])]
    return {"gene": first.get("gene", ""), "panel_ids": goi,
            "found_by": first.get("goi_found_by") or "unrecorded"}


def _panel_locus(panel_label):
    """``locus_4|gene-OMD`` -> ``locus_4``. Gene names may contain '-' but never
    '|', so splitting on the first pipe is safe."""
    return panel_label.split("|", 1)[0] if "|" in panel_label else panel_label


def load_phylo_placement_rows(phylo_dir):
    """Read the per-locus TSVs from PHYLO_PLACEMENT_CHECK (F9).

    Returns {(locus, model_id): row}, the locus from the file's name (<locus>.phylo_placement.tsv):
    model ids restart per home locus, so keyed by the id alone the last locus's file won and a verdict
    could land on another model of the same id (fix item 28). Header-only files (the no-op case, and
    the sentinel) are skipped, and a missing directory is not an error -- the check is optional.
    """
    rows = {}
    if not phylo_dir or not os.path.isdir(phylo_dir):
        return rows
    with os.scandir(phylo_dir) as it:
        for entry in sorted(it, key=lambda e: e.name):
            if not entry.is_file() or not entry.name.endswith(".tsv"):
                continue
            locus = entry.name.split(".", 1)[0]
            try:
                with open(entry.path) as fh:
                    for row in csv.DictReader(fh, delimiter="\t"):
                        mid = (row.get("model_id") or "").strip()
                        if mid:
                            rows[(locus, mid)] = row
            except OSError:
                continue
    return rows


def _phylo_row_of(rec, phylo_rows):
    """The F9 row of a record: the row of its own model (its representative's home locus and id), or
    of a member model at the record's span (the same model found from another home locus); a fitted
    verdict before insufficient_data, the representative's first; None without one. A model id names
    another model in another home locus (ids restart per locus), so a record never takes the row of
    an id it merely shares (fix item 28: a cow chrX record took decorin's verdict that way). A record
    without such a row gets none (F9 sees only HIGH and MEDIUM models)."""
    rep = tuple(rec.get("_rep") or ())
    span = (int(rec["start"]), int(rec["end"]))
    keys = [rep] if len(rep) == 2 else []
    keys += sorted({(m[0], m[1]) for m in rec.get("_members") or ()
                    if (int(m[2]), int(m[3])) == span and (m[0], m[1]) != rep})
    found = [phylo_rows[k] for k in keys if k in phylo_rows]
    return next((r for r in found if (r.get("verdict") or "").strip() != "insufficient_data"),
                found[0] if found else None)


def apply_phylo_placement(goi_dedup, phylo_rows):
    """Attach F9 phylogenetic-placement verdicts to the deduped GOI records.

    Why this check and not another threshold: identity, coverage, flanking support and
    home-paralog RBH were each measured on the melittin benchmark and none separates a
    true ortholog from a lineage-specific paralog -- the highest-identity call in the run
    (93.8 %, Bombus terrestris, melittin known lost) beats every true positive. Those are
    all PER-CALL signals, and one call in isolation cannot carry the answer. Whether a
    call's divergence tracks its SPECIES' divergence is a property of the whole set.

    ``phylo_discordant`` -- the call is more diverged than orthology predicts for its
    species distance. Recorded on the record and surfaced as a self-consistency flag. No
    label moves. (Until 2026-11-10 an option could raise a concordant AMBIGUOUS call to
    MEDIUM; since e3 a record takes the verdict of its own model, and F9 fits only the models
    the search labelled HIGH or MEDIUM, so no AMBIGUOUS record had a verdict and the option
    could not act: removed, main.nf removedParams.)

    Mutates records in place; returns (flags, summary).
    """
    records = goi_dedup.get("records", [])
    summary = {"evaluated": 0, "concordant": 0, "discordant": 0,
               "insufficient_data": 0}
    if not records or not phylo_rows:
        return [], summary

    flags = []
    for rec in records:
        row = _phylo_row_of(rec, phylo_rows)
        if not row:
            continue
        verdict = (row.get("verdict") or "").strip()
        rec["phylo_verdict"] = verdict
        rec["phylo_z"] = row.get("z_score", "")
        summary["evaluated"] += 1
        if verdict == "phylo_discordant":
            summary["discordant"] += 1
            flags.append({
                "type": "phylo_discordant",
                "genome": rec["genome"], "chrom": rec["chrom"], "start": rec["start"],
                "confidence": rec["confidence"], "identity": rec.get("identity", ""),
                "z_score": row.get("z_score", ""),
                "gene_divergence": row.get("gene_divergence", ""),
                "expected_divergence": row.get("expected_divergence", ""),
                "species_distance": row.get("species_distance", ""),
                "detail": ("sequence is more diverged than orthology predicts for this "
                           "species distance — consistent with a lineage-specific "
                           "paralog rather than the ortholog"),
            })
        elif verdict == "phylo_concordant":
            summary["concordant"] += 1
        else:
            summary["insufficient_data"] += 1
    return flags, summary


def build_locus_ownership(goi_dedup, ownership_rows, panel_meta, flanking_per_locus,
                          *, tiebreak_gap=10.0):
    """Re-attribute each deduped GOI ortholog to the home locus it *truly* belongs
    to (docs/TODO.md §1m), using reciprocal-best alignment against the home-paralog
    panel (the RBH leg) with a flanking-synteny tiebreak when the panel scores are
    near-tied.

    Why: home loci search independently, and cross-paralog flanking can anchor a
    block over the wrong gene — so the chr9 osteomodulin locus modelled real decorin
    on cow BTA5 while the chr12 decorin locus produced nothing there. §1h dedup keeps
    the best copy but attributes it to whichever locus *recovered* it (osteomodulin),
    and a true paralog (asporin) can keep the GOI label. This pass fixes both:

      * ``owning_locus`` — the panel locus the recovered protein aligns best to.
      * ``locus_reattributed`` flag — owning_locus differs from the recovered-from
        provenance. When owning_locus is *absent* from provenance entirely and the
        owner is the GOI's own gene, that's the §1m "owner search fumble" (the true
        home locus failed to model its own gene).
      * ``paralog_misassignment`` flag + ``inferred_paralog`` — the call aligns best
        to a non-GOI home paralog (asporin, not decorin); relabel before it's read
        as an ortholog of the GOI.

    A non-GOI gene owns a call only when it beats the GOI's best panel entry by at least
    ``tiebreak_gap``: below that the two cannot be told apart, and the call stays the GOI's
    (``ownership_near_tie``). An own copy counts as the GOI unless the call favours it over the
    GOI's own entry (_is_goi_for_call). A demo had 30 of 30 true orthologs relabelled when another
    home copy was identical to the GOI, a margin of 0 (fix-phase item 10). Neither case says the
    GOI's locus failed to recover the call, so neither raises locus_reattributed or
    goi_owner_search_fumble; nor does an answer of the genome-wide path (genome_wide_strong_hit),
    which is never counted and is dropped where its genome has a counted call. A paralog that
    owns an answer still demotes it (paralog_not_goi), so it is no answer.

    Mutates the records in ``goi_dedup['records']`` in place (adds owning_* keys) and
    returns (flags, summary)."""
    records = goi_dedup.get("records", [])
    if not records or not ownership_rows:
        return [], {"n_reattributed": 0, "n_paralog_misassignment": 0,
                    "n_owner_search_fumble": 0, "n_tiebreak_applied": 0, "n_near_tie": 0, "evaluated": 0,
                    "goi_entry": goi_panel_entry(panel_meta)}

    # Index ownership rows by genome+chrom for coordinate matching.
    # canonical_genome_id() on the row's genome is load-bearing: the rescue-pass
    # ownership tasks are tagged "<genome>.hull_rescue" so their output filename cannot
    # collide with the seeded task's, and that tag is written verbatim into the TSV's
    # `genome` column. The dedup records on the other side of this lookup are ALREADY
    # canonicalized, so without stripping it here the key never matches and every
    # rescue-derived call silently skips the RBH paralog check — i.e. the guard misses
    # exactly the calls most likely to need it (most HIGH calls now come from rescue).
    rows_by_gc = defaultdict(list)
    for r in ownership_rows:
        rows_by_gc[(canonical_genome_id(r["genome"]), r["chrom"])].append(r)

    flags = []
    n_reattributed = n_paralog = n_fumble = n_tiebreak = n_eval = n_near_tie = 0

    for rec in records:
        cands = rows_by_gc.get((rec["genome"], rec["chrom"]), [])
        if not cands:
            continue
        # Best-aligning ownership row overlapping this rep's span (multiple home loci
        # may have recovered the same gene; take the highest-scoring model).
        scored = []
        for r in cands:
            ov = _reciprocal_overlap(rec, {"start": r["start_i"], "end": r["end_i"]})
            if ov > 0.5:
                scored.append((r["best_bit_f"], ov, r))
        if not scored:
            continue
        scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
        row = scored[0][2]
        n_eval += 1

        best_label = row["best_paralog"]
        second_label = row.get("second_paralog", "")
        gap = row["gap_f"]
        owning_locus = _panel_locus(best_label)
        owning_gene = best_label.split("|", 1)[1] if "|" in best_label else best_label
        scores = _ownership_scores(row)
        owning_is_goi = _is_goi_for_call(best_label, scores, panel_meta)
        second_is_goi = bool(second_label) and _is_goi_for_call(second_label, scores, panel_meta)
        goi_entries = [(bit, label) for label, bit in scores.items()
                       if _is_goi_for_call(label, scores, panel_meta)]
        owner_label = best_label

        # A non-GOI gene on top, closer to the GOI's best entry than the gap: the GOI keeps the call.
        near_tie = False
        if not owning_is_goi and goi_entries and row["best_bit_f"] - max(goi_entries)[0] < tiebreak_gap:
            owner_label = max(goi_entries)[1]
            owning_locus = _panel_locus(owner_label)
            owning_gene = owner_label.split("|", 1)[1] if "|" in owner_label else owner_label
            owning_is_goi = near_tie = True
            n_near_tie += 1

        # --- Synteny tiebreak: when the top-2 panel genes are near-tied and neither or both are the
        # GOI, defer to which candidate locus's HIGH flanking neighbourhood is present on this target
        # chromosome (true local synteny beats a coin-flip alignment gap).
        tiebroken = False
        if not near_tie and second_label and gap < tiebreak_gap and owning_is_goi == second_is_goi:
            second_locus = _panel_locus(second_label)
            if second_locus != owning_locus:
                f_best = flanking_per_locus.get((rec["genome"], owning_locus, rec["chrom"]), 0)
                f_second = flanking_per_locus.get((rec["genome"], second_locus, rec["chrom"]), 0)
                if f_second > f_best:
                    owner_label = second_label
                    owning_locus = second_locus
                    owning_gene = second_label.split("|", 1)[1] if "|" in second_label else second_label
                    tiebroken = True
                    n_tiebreak += 1

        # An own copy owns the call, by a near tie or not: the call stays at a searched locus, not at
        # the copy's pseudo-locus home_<gene>, from which no search started.
        own_copy_owner = owning_is_goi and panel_meta.get(owner_label, {}).get("goi_found_by") == "own_copy"
        if own_copy_owner and owning_locus.startswith("home_"):
            own_entries = [(bit, label) for label, bit in scores.items() if panel_meta.get(label, {}).get("is_goi")
                           and panel_meta.get(label, {}).get("goi_found_by") != "own_copy"]
            if own_entries:
                owning_locus = _panel_locus(max(own_entries)[1])

        provenance = rec.get("provenance", [])
        # Re-attribution only makes sense when the owning gene IS the GOI (moving a true
        # GOI ortholog to its correct home locus). When the best match is a home PARALOG
        # (owning_is_goi=False — e.g. a homology-family entry like home_gene-BGN), that's a
        # paralog_misassignment, not a reattribution to a searched locus. A near tie or an own copy
        # says nothing about which home locus the call belongs to. Nor does an answer of the
        # genome-wide path: it lies outside synteny by construction, and the report drops it where
        # the genome has a counted call (drop_genome_wide_strong_hits), so it says nothing about
        # which home locus recovered the gene. Its paralog demotion below stays.
        cannot_tell = near_tie or own_copy_owner or bool(rec.get("genome_wide_strong_hit"))
        reattributed = owning_is_goi and not cannot_tell and (owning_locus not in provenance or (
            len(provenance) > 1 and provenance[0] != owning_locus
        ))
        owner_absent = owning_locus not in provenance and not cannot_tell

        rec["owning_locus"] = owning_locus
        rec["owning_gene"] = owning_gene
        rec["owning_is_goi_gene"] = owning_is_goi
        rec["ownership_bit"] = round(row["best_bit_f"], 1)
        rec["ownership_gap"] = round(gap, 1)
        rec["ownership_tiebroken"] = tiebroken
        if near_tie:
            rec["ownership_near_tie"] = True
        rec["locus_reattributed"] = reattributed
        if not owning_is_goi:
            rec["inferred_paralog"] = owning_gene
            # Don't destroy the original; mark the family-member reassignment.
            rec["goi_class"] = "paralog_not_goi"

        if reattributed:
            n_reattributed += 1
            flags.append({
                "type": "locus_reattributed",
                "genome": rec["genome"],
                "chrom": rec["chrom"],
                "start": rec["start"],
                "end": rec["end"],
                "recovered_from": provenance,
                "owning_locus": owning_locus,
                "owning_gene": owning_gene,
                "ownership_bit": round(row["best_bit_f"], 1),
                "ownership_gap": round(gap, 1),
                "tiebroken_by_synteny": tiebroken,
                "advice": (
                    f"Ortholog recovered from {provenance or '[]'} but aligns best to "
                    f"'{best_label}' — re-attributed to {owning_locus}."
                ),
            })
        if owner_absent and owning_is_goi:
            n_fumble += 1
            flags.append({
                "type": "goi_owner_search_fumble",
                "genome": rec["genome"],
                "chrom": rec["chrom"],
                "start": rec["start"],
                "end": rec["end"],
                "owning_locus": owning_locus,
                "owning_gene": owning_gene,
                "recovered_from": provenance,
                "advice": (
                    f"The GOI's own home locus ({owning_locus}/{owning_gene}) did NOT "
                    f"recover this ortholog — it was only found via {provenance}. The "
                    "true locus's flanking is likely rearranged in this genome so no "
                    "GOI-search block formed there."
                ),
            })
        if not owning_is_goi:
            n_paralog += 1
            flags.append({
                "type": "paralog_misassignment",
                "genome": rec["genome"],
                "chrom": rec["chrom"],
                "start": rec["start"],
                "end": rec["end"],
                "labeled_as": "GOI",
                "inferred_paralog": owning_gene,
                "owning_locus": owning_locus,
                "ownership_bit": round(row["best_bit_f"], 1),
                "ownership_gap": round(gap, 1),
                "advice": (
                    f"Call labeled as the GOI aligns best to home paralog '{owning_gene}' "
                    f"(bit={row['best_bit_f']:.1f}, gap={gap:.1f}) — it is most likely an "
                    f"ortholog of {owning_gene}, not the GOI. Relabeled goi_class=paralog_not_goi."
                ),
            })

    summary = {
        "evaluated": n_eval,
        "n_reattributed": n_reattributed,
        "n_paralog_misassignment": n_paralog,
        "n_owner_search_fumble": n_fumble,
        "n_tiebreak_applied": n_tiebreak,
        "n_near_tie": n_near_tie,
        "goi_entry": goi_panel_entry(panel_meta),
    }
    return flags, summary


def build_self_consistency(goi_annotations, goi_dedup, flanking_per_locus,
                           *, strong_flanking_min=20,
                           ownership_flags=None, ownership_summary=None,
                           phylo_flags=None, phylo_summary=None,
                           identity_decoupled_min_identity=50.0,
                           identity_decoupled_max_qcov=0.35,
                           core_completeness=False):
    """End-of-run sanity checks (docs/TODO.md §1j).

    Currently emits two flag types:

    * ``strong_synteny_no_goi`` — a (genome, locus, chrom) cell with >=
      ``strong_flanking_min`` HIGH-confidence flanking genes but zero HIGH-confidence
      GOI calls. That's the orthologous neighbourhood with a GOI too divergent for the
      current miniprot/classify thresholds; the advice asks the user to lower
      ``--miniprot_min_identity`` or inspect the block. (The cluster_grs side of this,
      §1e, surfaces the *block*; this flags the *genome+locus* as worth user attention.)

    * ``cross_locus_duplicate`` — same target gene found from >=2 home loci with
      bit-for-bit identical identity (already collapsed by §1h's
      ``dedupe_goi_annotations``; surfaced here so it is visible in the report
      alongside other self-consistency outputs).

    The reciprocal-best paralog check's ``paralog_confusion`` flag went with its step
    (RECIPROCAL_BEST_PARALOG, removed on 2026-11-16: it always received one query sequence and
    could flag nothing); locus ownership compares each call with the home genes
    (``ownership_flags``).
    """
    # Count HIGH-conf GOI per (genome, locus, chrom) for the identity-decay check.
    high_goi_count = defaultdict(int)
    for a in goi_annotations:
        if a["confidence"] == "HIGH":
            high_goi_count[(a["genome"], a["source"], a["chrom"])] += 1

    flags = []
    for (genome, locus, chrom), flank_n in flanking_per_locus.items():
        if flank_n >= strong_flanking_min and high_goi_count.get((genome, locus, chrom), 0) == 0:
            flags.append({
                "type": "strong_synteny_no_goi",
                "genome": genome,
                "locus": locus,
                "chrom": chrom,
                "high_flanking_count": flank_n,
                "high_goi_count": 0,
                "advice": (
                    "Strong synteny but no HIGH-confidence GOI call — inspect this block: "
                    "the orthologous neighbourhood is there, but no gene model in it is "
                    "significant, specific and complete enough for HIGH. Look at the "
                    "AMBIGUOUS (needs curation) and LOW calls of this genome; the gene may "
                    "be too divergent, broken in the assembly, or lost."
                ),
            })

    for rec in goi_dedup.get("records", []):
        if rec.get("cross_locus_duplicate"):
            flags.append({
                "type": "cross_locus_duplicate",
                "genome": rec["genome"],
                "chrom": rec["chrom"],
                "start": rec["start"],
                "end": rec["end"],
                "loci": rec["provenance"],
                "identity": rec["identity"],
                "note": "same gene identically modeled from multiple home loci; collapsed into one record",
            })

    # QW3 (docs/TODO_JUN.md): identity<->coverage decoupling. A GOI call with high
    # %identity but missing/low query coverage is a short high-identity local window
    # masquerading as a confident call — e.g. Vollenhovia "79% melittin" that is
    # really ~23% over the full query (the tandem-copy pident is a local-window
    # identity; QW2 now records its true QueryCoverage). Generic + advisory: no
    # demotion, no biology. Tetramorium's OV788322 ortholog (40% id) is below the
    # identity floor so it never trips here.
    for a in goi_annotations:
        ident = _safe_float(a.get("identity"), 0.0)
        if ident < identity_decoupled_min_identity:
            continue
        qc = _judged_coverage(a, core_completeness)
        if qc is None or qc < identity_decoupled_max_qcov:
            flags.append({
                "type": "identity_coverage_decoupled",
                "genome": a.get("genome"),
                "chrom": a.get("chrom"),
                "start": a.get("start"),
                "end": a.get("end"),
                "mrna_id": a.get("mrna_id"),
                "identity": ident,
                "query_coverage": qc,  # None = coverage not recorded for this call
                "coverage_recorded": qc is not None,
                "confidence": a.get("confidence"),
                "advice": (
                    "High %identity over a small fraction of the query — a short "
                    "high-identity local window, not a full-length match. This call is "
                    "DEMOTED out of the headline counts."
                    if qc is not None else
                    "High %identity but query coverage was never recorded for this call, "
                    "so it cannot be verified as full-length. Advisory only — not demoted; "
                    "inspect the alignment span manually."
                ),
            })

    # §1m locus-ownership flags (locus_reattributed / goi_owner_search_fumble /
    # paralog_misassignment) computed in build_locus_ownership and threaded through.
    ownership_flags = ownership_flags or []
    phylo_flags = phylo_flags or []
    ownership_summary = ownership_summary or {}
    flags.extend(ownership_flags)
    flags.extend(phylo_flags)

    checks_performed = ["strong_synteny_no_goi", "cross_locus_duplicate",
                        "identity_coverage_decoupled"]
    deferred_checks = []
    if ownership_summary.get("evaluated"):
        checks_performed.append("locus_ownership")
    else:
        deferred_checks.append("locus_ownership")
    if (phylo_summary or {}).get("evaluated"):
        checks_performed.append("phylo_placement")
    else:
        deferred_checks.append("phylo_placement")

    return {
        "checks_performed": checks_performed,
        "deferred_checks": deferred_checks,
        "thresholds": {
            "strong_flanking_min": strong_flanking_min,
            "identity_decoupled_min_identity": identity_decoupled_min_identity,
            "identity_decoupled_max_qcov": identity_decoupled_max_qcov,
        },
        "flags": flags,
        "locus_ownership_summary": ownership_summary,
        "phylo_placement_summary": phylo_summary or {},
        "summary": {
            "n_strong_synteny_no_goi": sum(1 for f in flags if f["type"] == "strong_synteny_no_goi"),
            "n_cross_locus_duplicates": sum(1 for f in flags if f["type"] == "cross_locus_duplicate"),
            "n_locus_reattributed": sum(1 for f in flags if f["type"] == "locus_reattributed"),
            "n_goi_owner_search_fumble": sum(1 for f in flags if f["type"] == "goi_owner_search_fumble"),
            "n_paralog_misassignment": sum(1 for f in flags if f["type"] == "paralog_misassignment"),
            "n_identity_coverage_decoupled": sum(1 for f in flags if f["type"] == "identity_coverage_decoupled"),
            "total_flags": len(flags),
        },
    }


def _judged_coverage(rec, core_completeness=False):
    """The coverage a record is tested on (constraint C6, as the judge does): its QueryCoverage,
    under core_completeness the larger of that and its CoreCoverage. None when neither is
    recorded."""
    query_raw = str(rec.get("query_cov") or "").strip()
    core_raw = str(rec.get("core_cov") or "").strip()
    return completeness_coverage(_safe_float(query_raw, None) if query_raw else None,
                                 _safe_float(core_raw, None) if core_raw else None, core_completeness)


def apply_coverage_demotion(goi_dedup, min_identity=50.0, max_qcov=0.35, core_completeness=False):
    """Demote a high-identity / low-coverage GOI call out of the headline counts.

    A short high-identity local window is not an ortholog call: 100 % over 10 aa of a
    70 aa query says almost nothing. ``build_self_consistency`` has flagged this as
    ``identity_coverage_decoupled`` since QW3, but only advisorily — the headline counted
    the call anyway. This makes the existing, already-calibrated detector actually demote.

    **Only demotes when coverage is KNOWN and low.** A call whose coverage was never
    recorded stays where it is and remains advisory: absence of evidence is not evidence
    of low coverage, and demoting on it would take out good full-length rescue models
    (before ``rescue_goi_hull.py`` learned to emit ``QueryCoverage``, that was every
    single one of them — including the true Apis cerana melittin).

    Mutates ``goi_dedup["records"]`` in place, mirroring how §1m relabels
    ``paralog_not_goi``. Returns a summary dict.
    """
    demoted, unknown = 0, 0
    for rec in goi_dedup.get("records", []):
        if rec.get("confidence") not in ("HIGH", "MEDIUM"):
            continue
        if _safe_float(rec.get("identity"), 0.0) < min_identity:
            continue
        qc = _judged_coverage(rec, core_completeness)
        if qc is None:
            unknown += 1
            continue
        if qc >= max_qcov:
            continue
        rec.setdefault("original_goi_class", rec.get("goi_class"))
        rec["goi_class"] = "identity_coverage_decoupled"
        rec["coverage_demoted"] = True
        rec["query_coverage"] = qc
        demoted += 1
    return {
        "demoted": demoted,
        "coverage_not_recorded": unknown,
        "min_identity": min_identity,
        "max_qcov": max_qcov,
    }


def build_report(results_dir, qc_json=None, qc_policy=None,
                 locus_ownership_tiebreak_gap=10.0, coverage_demotion=True,
                 identity_decoupled_min_identity=50.0, identity_decoupled_max_qcov=0.35,
                 phylo_placement_dir=None,
                 strong_flanking_min=20, core_completeness=False, tree_infos=None):
    qc_records = load_qc_records(qc_json)
    qc_summary = summarize_qc(qc_records)

    if phylo_placement_dir is None:
        phylo_placement_dir = os.path.join(results_dir, "phylo_placement")

    regions_dir = os.path.join(results_dir, "regions")
    hits_dir = os.path.join(results_dir, "hits")
    scores_dir = os.path.join(results_dir, "scores")

    # Single scandir() per directory, classify by suffix in one pass (docs/TODO.md §1c).
    regions_scan = scan_dir_by_suffix(regions_dir, {
        "fasta": (".faa", ".fna"),
        "gff": (".gff", ".gff3"),
    })
    scores_scan = scan_dir_by_suffix(scores_dir, {"scores": (".scores.tsv",)})
    # The search's per-genome tables ride along in hits/, so they need no channel of their
    # own; generate_report.nf stages each by name (the raw search hits are not staged).
    hits_scan = scan_dir_by_suffix(
        hits_dir, {"nonsyntenic": (NONSYNTENIC_SUFFIX,),
                   "not_searched": (NOT_SEARCHED_SUFFIX,),
                   "target_policy": (TARGET_POLICY_SUFFIX,),
                   "genome_wide": (GENOME_WIDE_SUFFIX,)})
    not_searched = collect_not_searched(hits_scan["not_searched"])
    target_policy = collect_target_policy(hits_scan["target_policy"])
    fasta_files = regions_scan["fasta"]
    gff_files = regions_scan["gff"]
    score_files = scores_scan["scores"]
    nonsyntenic_files = hits_scan["nonsyntenic"]
    nonsyntenic_dropped = []
    nonsyntenic_records, nonsyntenic_by_genome = collect_nonsyntenic_candidates(
        nonsyntenic_files, dropped_out=nonsyntenic_dropped)

    staging_diagnostics = build_staging_diagnostics(
        results_dir,
        dir_patterns={
            "regions_fasta": (regions_dir, ["*.faa", "*.fna"]),
            "regions_gff": (regions_dir, ["*.gff", "*.gff3"]),
            "scores": (scores_dir, ["*.scores.tsv"]),
        },
        match_counts={
            "fasta_files": len(fasta_files),
            "gff_files": len(gff_files),
            "score_files": len(score_files),
        },
    )

    genes_added_per_genome = summarize_fasta_outputs(fasta_files)
    # One read-pass over every GFF instead of three (docs/TODO.md §1c). The three
    # public functions (summarize_annotations / collect_goi_annotations /
    # collect_high_flanking_per_locus) still work for tests + external callers; they
    # just throw away the two outputs they don't need each.
    annotation_summary, goi_annotations, flanking_per_locus = _process_gffs_unified(gff_files)
    region_summary = summarize_region_scores(score_files)

    downstream_genomes = set(genes_added_per_genome)
    downstream_genomes |= {row["genome"] for row in annotation_summary["per_genome"]}
    downstream_genomes |= {row["genome"] for row in region_summary["per_genome"]}

    failed_downstream = []
    for rec in qc_summary["failed_genomes"]:
        genome = rec["genome"]
        if genome in downstream_genomes:
            failed_downstream.append(genome)

    # GOI ortholog headline metrics, from the GOI models of the region GFFs (one row per
    # model). The raw search hits do not reach the report.
    goi_conf = annotation_summary.get("goi_confidence_counts", {})
    goi_high = goi_conf.get("HIGH", 0)
    goi_medium = goi_conf.get("MEDIUM", 0)
    goi_low = goi_conf.get("LOW", 0)
    confident_goi = annotation_summary.get("goi_class_counts", {}).get("confident_goi", 0)
    resolved_goi = sum(row.get("resolved_goi_annotations", 0) for row in annotation_summary["per_genome"])
    genomes_with_goi = sum(1 for row in annotation_summary["per_genome"] if row.get("goi_annotations", 0) > 0)

    # Cross-locus dedup: the same target ortholog is reported once per home locus it was
    # found from (e.g. the luciferase rerun's Aquatica gene appears identically via locus_3
    # and locus_4). Collapse those so the headline counts distinct genes, not seeds, and
    # flag bit-for-bit-identical multi-locus hits as cross_locus_duplicate. docs/TODO.md §1h.
    # ``goi_annotations`` + ``flanking_per_locus`` came out of the unified GFF pass above.
    # The answers of the genome-wide path are deduped apart and appended (dedupe_goi_records).
    goi_dedup = dedupe_goi_records(goi_annotations)
    # F9 phylogenetic placement: each record gets its own model's verdict (advisory, no label
    # moves; the counts below are read after it, as when an option could promote, until 2026-11-10).
    phylo_rows = load_phylo_placement_rows(phylo_placement_dir)
    phylo_flags, phylo_summary = apply_phylo_placement(
        goi_dedup, phylo_rows)
    if phylo_flags or phylo_summary["evaluated"]:
        recount = [r for r in goi_dedup.get("records", []) if not r.get("genome_wide_strong_hit")]
        goi_dedup["high_confidence_goi_deduped"] = sum(
            1 for r in recount if r["confidence"] == "HIGH")
        goi_dedup["medium_confidence_goi_deduped"] = sum(
            1 for r in recount if r["confidence"] == "MEDIUM")
        goi_dedup["ambiguous_goi_deduped"] = sum(
            1 for r in recount if r["confidence"] == "AMBIGUOUS")

    dedup_high = goi_dedup["high_confidence_goi_deduped"]
    dedup_medium = goi_dedup["medium_confidence_goi_deduped"]

    # §1m locus ownership: reciprocal-best of each deduped GOI ortholog against the
    # home-paralog panel (built by build_home_paralog_panel.py, RBH run by
    # ASSIGN_LOCUS_OWNERSHIP and staged under locus_ownership/). Re-attributes
    # orthologs the §1h dedup filed under the wrong home locus and relabels paralogs
    # carried under the GOI name. Mutates goi_dedup records in place.
    ownership_rows = load_paralog_check_rows(
        os.path.join(results_dir, "locus_ownership"),
        suffix=".locus_ownership.tsv",
    )
    panel_meta = load_panel_meta(
        os.path.join(results_dir, "locus_ownership", "home_paralog_panel.faa.meta.tsv")
    )
    ownership_flags, ownership_summary = build_locus_ownership(
        goi_dedup, ownership_rows, panel_meta, flanking_per_locus,
        tiebreak_gap=locus_ownership_tiebreak_gap,
    )
    # Steps that failed and were recorded (fix-phase d5): a failed ownership check marks the
    # calls it would have checked; every failure is a warning below.
    step_failures = collect_step_failures(results_dir)
    ownership_not_checked = mark_ownership_not_checked(goi_dedup, step_failures)

    # §1m: ownership may relabel a call goi_class=paralog_not_goi (its best home match is a
    # paralog, not the GOI — e.g. a 55% biglycan hit mislabeled GOI_DCN). Such a call is NOT
    # a GOI ortholog, so drop it from the headline HIGH/MEDIUM counts. It stays in
    # goi_dedup["records"] with its flag for transparency; the pre-ownership counts are kept
    # as *_pre_ownership. This is what makes the RBH check actually demote the false positive.
    # A high-identity call over a tiny slice of the query is not an ortholog call either.
    # Same treatment as paralog_not_goi: relabel, keep the record, drop it from the
    # headline. Only fires on KNOWN-low coverage (see apply_coverage_demotion).
    coverage_demotion_summary = (
        apply_coverage_demotion(goi_dedup,
                                min_identity=identity_decoupled_min_identity,
                                max_qcov=identity_decoupled_max_qcov,
                                core_completeness=core_completeness)
        if coverage_demotion else
        {"demoted": 0, "coverage_not_recorded": 0, "disabled": True}
    )

    _DEMOTED_CLASSES = {"paralog_not_goi", "identity_coverage_decoupled"}
    # The answers of the genome-wide path: one per genome, only where no call is counted.
    _dropped_answers = drop_genome_wide_strong_hits(goi_dedup, _DEMOTED_CLASSES)
    _dropped_gene_order = drop_gene_order_only(goi_dedup, _DEMOTED_CLASSES)
    # (record, why) of every record the two drops removed: a region they leave without a record
    # still names its model (goi_regions_without_call).
    _report_dropped = (
        list(zip(_dropped_answers,
                 (d["reason"] for d in goi_dedup.get("genome_wide_strong_hits_dropped_records", []))))
        + [(r, "gene_order_only") for r in _dropped_gene_order])
    absence = genome_absence_categories(goi_dedup.get("records", []), target_policy,
                                        _DEMOTED_CLASSES)
    dedup_high_pre_ownership, dedup_medium_pre_ownership = dedup_high, dedup_medium
    _own_recs = goi_dedup.get("records", [])
    _counted_recs = [r for r in _own_recs if _is_counted(r, _DEMOTED_CLASSES)]
    # The LOW models, listed apart (never counted): one entry per distinct model.
    low_block = goi_low_records(goi_annotations, _own_recs, _DEMOTED_CLASSES, ownership_rows)
    # Gene-bearing regions of MEDIUM or better (and goi_hit_no_call regions) without a counted
    # call, with what they hold and why.
    region_summary["goi_regions_without_call"] = goi_regions_without_call(
        collect_region_rows(score_files), goi_annotations, _own_recs, nonsyntenic_records,
        nonsyntenic_dropped, _DEMOTED_CLASSES, report_dropped=_report_dropped)
    dedup_high = sum(1 for r in _counted_recs if r.get("confidence") == "HIGH")
    dedup_medium = sum(1 for r in _counted_recs if r.get("confidence") == "MEDIUM")
    # The answers of the genome-wide path the report keeps (a reading of their own, never
    # counted): outside every searched window, or near a block (inside its window).
    _answers = genome_wide_answer_records(_own_recs, _DEMOTED_CLASSES)
    _near_block_answers = [r for r in _answers if str(r.get("near_block") or "").strip()]
    _no_synteny_answers = [r for r in _answers if not str(r.get("near_block") or "").strip()]
    # §1x: the refused hits per place; a place that a counted call of the run overlaps (from any
    # home locus or rescue pass) is "modelled", not "found, not syntenic".
    nonsyntenic_block = summarize_nonsyntenic_candidates(
        nonsyntenic_records, nonsyntenic_by_genome, counted_records=_counted_recs,
        dropped=nonsyntenic_dropped,
        dropped_hits=sum_target_policy_column(hits_scan["target_policy"], "rescue_dropped_outside_zones"))
    # Calls found only because the dispersed-GOI rescue made the search look there.
    dispersed_rescue_by_rule = dict(sorted(Counter(
        r.get("rescue_rule") or "unknown" for r in _counted_recs if r.get("dispersed_rescue")).items()))
    # AMBIGUOUS calls: shown as "needs curation", never counted, split by why (Basis). One that
    # ownership gave to another home gene (paralog_not_goi) is that gene's call, named apart.
    _curation = [r for r in _own_recs if r.get("confidence") == "AMBIGUOUS"
                 and r.get("goi_class") != "paralog_not_goi" and not r.get("genome_wide_strong_hit")]
    dedup_ambiguous = len(_curation)
    ambiguous_other_home_gene = sum(1 for r in _own_recs if r.get("confidence") == "AMBIGUOUS"
                                    and r.get("goi_class") == "paralog_not_goi"
                                    and not r.get("genome_wide_strong_hit"))
    needs_curation_by_basis = dict(Counter(r.get("basis") or "unlabelled" for r in _curation))

    # End-of-run self-consistency checks (docs/TODO.md §1j): identity-decay sanity +
    # cross-locus duplicate visibility + locus ownership (§1m).
    # The flanking-only-block surfacing is the §1e piece that lives in cluster_grs.py.
    goi_dedup["coverage_demotion"] = coverage_demotion_summary
    self_consistency = build_self_consistency(
        goi_annotations, goi_dedup, flanking_per_locus,
        strong_flanking_min=strong_flanking_min,
        ownership_flags=ownership_flags,
        phylo_flags=phylo_flags,
        phylo_summary=phylo_summary,
        ownership_summary=ownership_summary,
        identity_decoupled_min_identity=identity_decoupled_min_identity,
        identity_decoupled_max_qcov=identity_decoupled_max_qcov,
        core_completeness=core_completeness,
    )

    # Detect a gene-family / multi-paralog query and advise treating GOI calls as
    # candidates (SynVoy is validated for divergent low-copy genes, not large paralog
    # families). Report-only, driven by the share of records another home gene explains as
    # well and by the home-locus count.
    gene_family_advisory = build_gene_family_advisory(
        goi_dedup.get("records", []),
        n_home_loci=count_home_loci(gff_files),
    )

    # A genome that failed the quality check and was searched anyway (the user chose it):
    # every call in it carries the flag, and the report says so at the top.
    low_quality_ids, warnings = low_quality_warnings(
        qc_summary["failed_genomes"], failed_downstream)
    warnings.extend(step_failure_warnings(step_failures))
    if not_searched:
        warnings.append(
            f"{len(not_searched)} genome(s) could not be searched and have no result: "
            f"{', '.join(sorted(not_searched))}. They are not GOI-absent "
            "(summary.genomes_not_searched).")
        too_small = sorted(g for g, reasons in not_searched.items()
                           if "mmseqs_split_memory_limit_too_small" in reasons)
        if too_small:
            warnings.append(
                f"MMseqs2 could not fit the database of {len(too_small)} genome(s) into its split limit "
                f"({', '.join(too_small)}; reason mmseqs_split_memory_limit_too_small): raise "
                f"--mmseqs_split_memory_limit, and --iterative_search_memory when it is set to a size "
                f"(as in the laptop_safe and low_mem profiles), and run again with -resume.")
    for rec in goi_dedup.get("records", []):
        if rec.get("genome") in low_quality_ids:
            rec["low_quality_assembly"] = True
    # Genomes where a search found no flanking gene at all: the GOI's own hits seeded every
    # block there, with no neighbour test.
    no_flanking_genomes = genomes_without_flanking_hits(hits_scan["target_policy"])
    if no_flanking_genomes:
        warnings.append(
            f"{len(no_flanking_genomes)} genome(s) had no hit of any flanking gene in at least one "
            f"home locus's search: {', '.join(no_flanking_genomes)}. Every cluster of GOI hits there "
            "seeded a block of its own, so a call there rests on the protein match alone "
            "(DispersedRescue=1, RescueRule=no_flanking_hits; summary.genomes_without_flanking_hits).")
    genome_wide_records = collect_genome_wide_places(hits_scan["genome_wide"])
    _no_synteny_genomes = sorted({r["genome"] for r in _no_synteny_answers})
    _near_block_genomes = sorted({r["genome"] for r in _near_block_answers})
    if _no_synteny_genomes:
        warnings.append(
            f"{len(_no_synteny_genomes)} genome(s) without a counted call hold a strong hit without "
            "synteny (genome-wide path; MEDIUM, not counted): summary.genome_wide_strong_hit_genomes.")
    if _near_block_genomes:
        warnings.append(
            f"{len(_near_block_genomes)} genome(s) without a counted call hold a strong hit near a "
            "block, outside its GOI zone (genome-wide path; MEDIUM, not counted): "
            "summary.near_block_strong_hit_genomes.")

    report = {
        # Things a reader must know before reading any number below. Empty on a clean run.
        "warnings": warnings,
        "genome_qc": qc_records,
        "qc_summary": {
            **qc_summary,
            "qc_fail_policy": qc_policy or "unspecified",
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
        },
        "gene_family_advisory": gene_family_advisory,
        "synteny_results": {
            "genes_discovered": genes_added_per_genome,
        },
        "annotations": annotation_summary,
        "regions": region_summary,
        "goi_dedup": goi_dedup,
        # The LOW models (not orthology calls, counted nowhere), one entry per distinct model.
        "goi_low_records": low_block,
        # §1x — the run's known false negatives: strong GOI hits it refused on synteny
        # grounds and used to discard in silence.
        "rejected_candidates": nonsyntenic_block,
        # Whenever the genome-wide path ran (its tables were staged), also without a place.
        **({"genome_wide_path": summarize_genome_wide_places(genome_wide_records)}
           if genome_wide_records or hits_scan["genome_wide"] else {}),
        # Per target genome: the median flanking identity, the blocks searched and whether
        # the best block's window ran past a scaffold end (hits/<genome>.target_policy.tsv).
        "target_policy": target_policy,
        "self_consistency": self_consistency,
        "staging_diagnostics": staging_diagnostics,
        "summary": {
            # --- Headline: the number(s) a user actually cares about ---
            # Post-dedup (docs/TODO.md §1h): distinct HIGH/MEDIUM ortholog genes, not the
            # per-home-locus seed count, so the same gene found from N loci reads as one.
            "headline": _format_goi_headline(dedup_high, dedup_medium, goi_low,
                                             genomes_with_goi, ambiguous=dedup_ambiguous,
                                             curation_by_basis=needs_curation_by_basis,
                                             other_home_gene=ambiguous_other_home_gene,
                                             genome_wide=len(_no_synteny_answers),
                                             near_block=len(_near_block_answers)),
            # How many calls the coverage rule removed from the headline above.
            "coverage_demoted_calls": coverage_demotion_summary.get("demoted", 0),
            "coverage_not_recorded_calls": coverage_demotion_summary.get("coverage_not_recorded", 0),
            "headline_metric": dedup_high,  # distinct high-confidence GOI orthologs (post-dedup)
            "high_confidence_goi": dedup_high,
            "medium_confidence_goi": dedup_medium,
            # AMBIGUOUS calls, reported as "needs curation" and never folded into the ortholog
            # counts (see _format_goi_headline), and why: gene_order_only = a weak match in
            # the slot where the gene was expected; not_specific = another home gene explains
            # the match as well. ambiguous_goi is the old name of the same number (kept for
            # one release). ambiguous_other_home_gene_goi: AMBIGUOUS calls that ownership
            # gave to another home gene (paralog_not_goi); they need no curation for this GOI.
            "needs_curation_goi": dedup_ambiguous,
            "needs_curation_by_basis": needs_curation_by_basis,
            "ambiguous_other_home_gene_goi": ambiguous_other_home_gene,
            "ambiguous_goi": dedup_ambiguous,
            # Calls (AMBIGUOUS or better) whose specificity left an own copy of q out of p*
            # (record field own_copy_set_aside; fix-phase d2).
            "own_copy_set_aside_goi": sum(1 for r in goi_dedup.get("records", [])
                                          if r.get("own_copy_set_aside")
                                          and not r.get("genome_wide_strong_hit")),
            # The answers of the genome-wide path the report kept (one per genome, only where no
            # call is counted; never counted, not in the headline's ortholog count): outside every
            # searched window, and near a block (inside its window, outside its GOI zone); and
            # how many the report dropped (goi_dedup.genome_wide_strong_hits_dropped_records).
            "genome_wide_strong_hits": len(_no_synteny_answers),
            "genome_wide_strong_hit_genomes": _no_synteny_genomes,
            "near_block_strong_hits": len(_near_block_answers),
            "near_block_strong_hit_genomes": _near_block_genomes,
            "genome_wide_strong_hits_dropped": goi_dedup.get("genome_wide_strong_hits_dropped", 0),
            "gene_order_only_dropped": goi_dedup.get("gene_order_only_dropped", 0),
            "phylo_discordant_calls": phylo_summary["discordant"],
            "high_confidence_goi_pre_dedup": goi_high,
            "medium_confidence_goi_pre_dedup": goi_medium,
            # Pre-ownership = before §1m relabeled paralog_not_goi calls out of the headline.
            "high_confidence_goi_pre_ownership": dedup_high_pre_ownership,
            "medium_confidence_goi_pre_ownership": dedup_medium_pre_ownership,
            "cross_locus_duplicate_goi": goi_dedup["cross_locus_duplicates"],
            "goi_hits_collapsed_by_dedup": goi_dedup["hits_collapsed_by_dedup"],
            "self_consistency_flag_count": self_consistency["summary"]["total_flags"],
            "strong_synteny_no_goi_flags": self_consistency["summary"]["n_strong_synteny_no_goi"],
            "locus_reattributed_orthologs": self_consistency["summary"]["n_locus_reattributed"],
            "goi_owner_search_fumbles": self_consistency["summary"]["n_goi_owner_search_fumble"],
            "paralog_misassignments": self_consistency["summary"]["n_paralog_misassignment"],
            "gene_family_query": gene_family_advisory["is_gene_family_query"],
            "low_confidence_goi": goi_low,
            "confident_goi_annotations": confident_goi,
            # Class-based, deprecated (kept for this release): GOI rows of every class but the
            # LOW class ambiguous_goi_family_member, counted before the dedup and ownership.
            "resolved_goi_annotations": resolved_goi,
            "genomes_with_goi_annotations": genomes_with_goi,
            "total_new_genes": sum(genes_added_per_genome.values()),
            "genomes_with_annotations": len(annotation_summary["per_genome"]),
            "total_annotations": annotation_summary["total_annotations"],
            "total_goi_annotations": annotation_summary["role_counts"].get("goi", 0),
            # Class-based, deprecated (kept for this release): the rows of the LOW class
            # ambiguous_goi_family_member; low_confidence_goi counts LOW rows by label.
            "ambiguous_goi_annotations": annotation_summary["goi_class_counts"].get("ambiguous_goi_family_member", 0),
            "fallback_goi_annotations": annotation_summary["fallback_goi_annotations"],
            "low_confidence_regions": region_summary["confidence_counts"].get("LOW", 0),
            # Gene-bearing regions of MEDIUM or better without a counted call
            # (regions.goi_regions_without_call, with the reason for each).
            "goi_regions_without_call": len(region_summary["goi_regions_without_call"]),
            # §1x: "we found nothing" and "we found something and could not place it"
            # are opposite claims. goi_absent_genomes is now the former ONLY; a genome
            # holding a strong-but-unplaceable hit moves to goi_found_but_not_syntenic.
            # Derived from the genomes that produced no GOI call, which is NOT the same
            # set as genomes_without_goi: a genome whose search emitted no GFF at all
            # never reaches that list, yet is exactly the case where a refused hit is
            # the only thing we know about it.
            "goi_absent_genomes": sorted(
                g for g in annotation_summary["genomes_without_goi"]
                if g not in nonsyntenic_by_genome
                and g not in {r["genome"] for r in _answers}),
            # Genomes without a counted call (records after ownership, the coverage rule and
            # the drops above) that hold a refused place no counted call of the run overlaps
            # (another home locus's, a rescue pass's). Such a genome may hold LOW or AMBIGUOUS
            # models, or a call ownership gave to another home gene. A genome holding a kept
            # answer of the genome-wide path is in the strong_hit lists.
            "goi_found_but_not_syntenic": sorted(
                set(nonsyntenic_block["unmodelled_genomes"])
                - counted_genomes(_own_recs, _DEMOTED_CLASSES) - {r["genome"] for r in _answers}),
            "nonsyntenic_candidate_count": nonsyntenic_block["total"],
            "nonsyntenic_place_count": nonsyntenic_block["place_count"],
            "nonsyntenic_unmodelled_place_count": nonsyntenic_block["unmodelled_place_count"],
            # Counted calls (HIGH or MEDIUM, not demoted) found only because the dispersed-GOI
            # rescue made the search look there (DispersedRescue=1), and by its rule.
            "dispersed_rescue_goi": sum(dispersed_rescue_by_rule.values()),
            "dispersed_rescue_by_rule": dispersed_rescue_by_rule,
            # Counted calls whose model is a fragment (ModelStatus=fragment), and those that
            # are hit chains the fallback gate refused, written because their block is
            # significant (RefusalLifted; they never seed and never reach the tree).
            "counted_fragment_calls": sum(1 for r in _counted_recs if r.get("model_status") == "fragment"),
            "refusal_lifted_calls": sum(1 for r in _counted_recs if r.get("refusal_lifted")),
            # Counted calls of a tandem-copy lineage (record field copy_lineage): their identity
            # and coverage refer to the copy; home_identity / home_query_cov to the home protein.
            "counted_copy_lineage_calls": sum(1 for r in _counted_recs if r.get("copy_lineage")),
            "genomes_without_flanking_hits": no_flanking_genomes,
            # Genomes the search skipped (e.g. MMseqs resource failure): no result at all.
            "genomes_not_searched": {g: r for g, r in sorted(not_searched.items())},
            # Genomes without a counted call, by what the search saw there
            # (genome_absence_categories): the best block complete, its window cut at a
            # scaffold end ("unknown", not "absent"), or no block at all; and genomes whose
            # best call is a model another home gene explains as well.
            "absent_context_complete": absence["absent_context_complete"],
            "unknown_context_truncated": absence["unknown_context_truncated"],
            "absent_no_block": absence["absent_no_block"],
            "strong_hit_without_synteny": absence["strong_hit_without_synteny"],
            "strong_hit_near_block": absence["strong_hit_near_block"],
            "nearest_paralog_only": absence["nearest_paralog_only"],
            "nearest_paralog_only_genomes": len(absence["nearest_paralog_only"]),
            # Class-based count of the LOW class (deprecated; kept for this release):
            # goi_low_only_genomes is the same reading by label.
            "goi_ambiguous_only_genomes": annotation_summary["genomes_with_only_ambiguous_goi"],
            "goi_low_only_genomes": annotation_summary["genomes_with_only_low_goi"],
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
            "low_quality_genomes_searched": low_quality_ids,
            "staging_empty": staging_diagnostics["empty"],
        },
    }
    # Decoy mode (the low_complexity module or --decoy_mode): the run's decoy rate.
    decoy_rate = collect_decoy_rate(hits_scan["target_policy"])
    if decoy_rate is not None:
        report["decoy_rate"] = decoy_rate
    # What each limit of the search dropped (counts; they change no call).
    limits = collect_limits(hits_scan["target_policy"])
    if limits is not None:
        report["limits"] = limits
        # MMseqs2 on a lower rung of its low-memory ladder (fix-phase d6): fewer hits possible.
        for key, what in (("mmseqs_genome_lower_rung", "genome search(es)"),
                          ("mmseqs_window_lower_rung", "window search(es)")):
            c = limits["counts"].get(key) or {}
            if c.get("total"):
                report["warnings"].append(
                    f"MMseqs2 ran on a lower rung of its low-memory ladder (sensitivity 7.0 or 6.0, one "
                    f"thread) after a memory failure in {c['total']} {what} ({', '.join(c['genomes'])}): "
                    f"a hit, and so a call, may be missing there (limits.{key}; the genome's GFF says "
                    f"which, ##synvoy-reduced-search). Give MMseqs2 more memory.")
        # The model cap of a GOI locus stopped its copy rounds with a copy possibly left.
        c = limits["counts"].get("goi_model_cap_stopped") or {}
        if c.get("total"):
            report["warnings"].append(
                f"The model cap of a GOI locus stopped the copy rounds {c['total']} time(s) "
                f"({', '.join(c['genomes'])}) while an alignment that would have been kept, or a window "
                f"hit that no model explains, was left: a copy of the gene may be missing there "
                f"(limits.goi_model_cap_stopped). Raise --goi_max_models_per_locus.")
        # MMseqs2 keeps at most --mmseqs_max_seqs targets per query.
        maxseqs = [(key, limits["counts"].get(key) or {})
                   for key in ("mmseqs_maxseqs_genome", "mmseqs_maxseqs_window")]
        if any(c.get("total") for _key, c in maxseqs):
            where = "; ".join(f"{key} {c['total']} ({', '.join(c['genomes'])})"
                              for key, c in maxseqs if c.get("total"))
            report["warnings"].append(
                f"A query reached MMseqs2's limit of targets per query (--mmseqs_max_seqs) in a search: "
                f"its further hits there may be missing ({where}; limits). The counts are a lower "
                f"bound: they count queries with at least that many result rows, and a query whose "
                f"prefilter kept that many targets but fewer of them passed the E-value is not "
                f"counted. Raise --mmseqs_max_seqs.")
    # What ran the Smith-Waterman search (fix-phase item 18); main() moves it into run_settings.
    sw_backend = collect_sw_backend(hits_scan["target_policy"])
    if sw_backend is not None:
        report["smith_waterman"] = sw_backend
    # Steps that failed and were recorded (fix-phase d5); absent on a run where none failed.
    if step_failures:
        report["failed_steps"] = step_failures
        report["summary"]["ownership_not_checked_calls"] = ownership_not_checked
    # The figures (item 36): per home locus with a GFF of the search, the outcome of each GOI model
    # the report merged, dropped or does not count, and the other loci that found a record drawn
    # there. apply_report_drops.py marks the per-genome GFFs the figures are drawn from with them, so
    # the figures draw each record once. Written for every locus, with or without a tree; empty lists
    # mean that the report has nothing to mark there.
    report["figures"] = {
        locus: figure_marks(goi_annotations, goi_dedup.get("records", []), locus, _DEMOTED_CLASSES,
                            dropped=_report_dropped)
        for locus in sorted({staging_source_label(p) for p in gff_files
                             if _STAGING_PREFIX_RE.match(os.path.basename(p)) and not rescue_pass_of(p)})}
    # The trees (COMPUTE_TREE's infos): the leaves the report does not count (fix item 46). Only where
    # a tree was inferred: a run with --skip_tree has no trees block.
    if tree_infos:
        report["trees"] = tree_infos
        for locus, info in sorted(tree_infos.items()):
            if not isinstance(info, dict) or not info.get("leaves"):
                continue
            info["not_counted"] = tree_leaves_not_counted(goi_annotations, goi_dedup.get("records", []),
                                                          locus, info["leaves"])
            if info["not_counted"]:
                models = {a.get("mrna_id") for a in goi_annotations if a.get("source") == locus}
                n_models = sum(1 for leaf in info["leaves"] if leaf in models)
                n_out = len(info["not_counted"])
                verb, pron = ("is a model", "it") if n_out == 1 else ("are models", "them")
                report["warnings"].append(
                    f"{n_out} of the {n_models} model leaves of the {locus} tree {verb} the report does not "
                    f"count (merged into another model's call, overruled, or another home gene's): the tree "
                    f"page draws {pron} hollow (trees.{locus}.not_counted).")
    for rec in goi_dedup.get("records", []):
        for key in ("_rep", "_members", "_overruled"):
            rec.pop(key, None)
    return report


def format_empty_staging_message(diagnostics):
    lines = [
        "ERROR: generate_report found zero annotation and zero region files under staged_results.",
        "This usually means ITERATIVE_SEARCH produced no hits, or the Nextflow channel wiring is broken.",
        f"  results_dir: {diagnostics['results_dir']} (exists={diagnostics['results_dir_exists']})",
    ]
    for name, info in diagnostics["dirs"].items():
        sample = ", ".join(info["sample_entries"]) if info["sample_entries"] else "(none)"
        lines.append(
            f"  {name}: path={info['path']} exists={info['exists']} "
            f"entries={info['entry_count']} matches={info['pattern_match_count']} "
            f"patterns={info['patterns']} sample={sample}"
        )
    lines.append(
        "Try: inspect work/<hash>/iterative_results/ for the ITERATIVE_SEARCH process; "
        "check logs/iterative_search/ for zero-hit reports; "
        "rerun with --allow-empty if zero-hit is expected."
    )
    return "\n".join(lines)


RUN_SETTINGS_KEYS = ("modules_mode", "modules_applied", "modules_would_fire", "preset_applied",
                     "apply_reason", "legacy_preset", "user_set_kept", "settings", "source",
                     "features", "home_copies", "profile", "derivation", "constraints")


def copies_in_slot_by_genome(records, home_copies):
    """Per genome, the copies of the gene family at home and in the target's slot (preset
    plan of 2026-10-15, PR-6; only with --report_copies_in_slot, whose CopiesInSlot
    attribute this reads): the slot is the one of the block that holds the genome's best
    record with a slot (the judge's order). {} when no record carries the attribute."""
    by_genome = {}
    for r in records:
        if str(r.get("copies_in_slot") or "").strip():
            by_genome.setdefault(r.get("genome"), []).append(r)
    out = {}
    for genome in sorted(g for g in by_genome if g):
        best = min(by_genome[genome], key=_record_rank)
        try:
            m = int(str(best["copies_in_slot"]).strip())
        except ValueError:
            continue
        n = home_copies if isinstance(home_copies, int) else None
        out[genome] = {
            "home_copies": n,
            "copies_in_slot": m,
            "unit": f"{n if n is not None else '?'}:{m}",
            "slot_of": {"chrom": best.get("chrom"), "start": best.get("start"),
                        "end": best.get("end"), "confidence": best.get("confidence"),
                        "slot_rank": best.get("slot_rank") or None},
        }
    return out


def load_sort_provenance(folder):
    """The sort_provenance.json files of PHYLO_SORT (one per home locus; fix-phase d7): per
    file its strategy and the NCBI lookups that did not answer, unique, sorted. [] without any."""
    seen = set()
    for root, _dirs, files in os.walk(folder or ""):
        for name in files:
            if name != "sort_provenance.json":
                continue
            try:
                with open(os.path.join(root, name)) as fh:
                    p = json.load(fh)
            except (OSError, ValueError):
                continue
            if isinstance(p, dict):
                seen.add((p.get("strategy") or "", tuple(p.get("network_failures") or ())))
    return [{"strategy": s, **({"network_failures": list(f)} if f else {})} for s, f in sorted(seen)]


def load_download_record(path):
    """Rows of easy mode's download_record.tsv (fix-phase d7): accession, species, genome, gff."""
    if not path or not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def _member_overlap(m, r):
    """The overlap of a record's member (locus, model id, start, end) with the record's span, as a
    share of the stretch the two cover together."""
    s, e = int(m[2]), int(m[3])
    return (max(0, min(e, int(r["end"])) - max(s, int(r["start"])) + 1)
            / max(1, max(e, int(r["end"])) - min(s, int(r["start"])) + 1))


def _shown_member(r, locus, ids):
    """The model id of home locus ``locus`` that draws the record ``r`` there, in its tree and in its
    figures: among the record's members of that locus whose id is in ``ids``, the record's
    representative, else the member at the record's span (the same model found from another home
    locus), else the member that overlaps the record most, then the smallest id. None when no member
    of the record is in ``ids``."""
    mine = [m for m in r.get("_members") or () if m[0] == locus and m[1] in ids]
    if not mine:
        return None
    rep, span = tuple(r.get("_rep") or ()), (int(r["start"]), int(r["end"]))
    return min(mine, key=lambda m: ((m[0], m[1]) != rep, (int(m[2]), int(m[3])) != span,
                                    -_member_overlap(m, r), m[1]))[1]


def tree_leaves_not_counted(goi_annotations, records, locus, leaves,
                            not_counted_classes=("paralog_not_goi", "identity_coverage_decoupled")):
    """The leaves of ``locus``'s tree that the report does not count, each with its reason (fix item
    46). The tree is inferred before the report, from every HIGH or MEDIUM model of the search, so
    a model the report merged into another, relabelled or overruled stayed a leaf, drawn as a copy
    of its species. A leaf is a model id of this locus (ids repeat across loci); its record is the
    one it is a member of (the dedup's members, (locus, model id)). A counted record (HIGH or
    MEDIUM, not of a class in ``not_counted_classes``) is drawn once in each tree it has leaves in:
    by its representative where that is a model of this locus, else by its member at the record's
    span (the same model found from another home locus), else by the member that overlaps the
    record most. Reasons: 'merged' (another leaf of this tree is the gene's), 'overruled_piece' (a
    piece the record's whole gene overruled), the record's class (paralog_not_goi,
    identity_coverage_decoupled), 'not_in_report' (no record holds it, e.g. one the report dropped),
    else the record's confidence (no case known). Leaves without a GOI row of this locus (the home
    protein, the query) are left out: [{'leaf', 'reason', 'record'}, ...] in the leaves' order."""
    rows = {}
    for a in goi_annotations:
        if a.get("source") == locus and a.get("mrna_id"):
            rows.setdefault(a["mrna_id"], a)
    member_of = {}
    for r in records:
        for m in r.get("_members") or ():
            if m[0] == locus and m[1]:
                member_of.setdefault(m[1], r)

    def counted(r):
        return r.get("confidence") in ("HIGH", "MEDIUM") and r.get("goi_class") not in not_counted_classes

    leaf_set = set(leaves)
    shown = {}
    for r in records:
        if not counted(r):
            continue
        member = _shown_member(r, locus, leaf_set)
        if member is not None:
            shown[id(r)] = member
    out = []
    for leaf in leaves:
        if leaf not in rows:
            continue
        rec = member_of.get(leaf)
        if rec is None:
            reason = "not_in_report"
        elif counted(rec):
            if shown.get(id(rec)) == leaf:
                continue
            reason = "merged"
        elif (locus, leaf) in {tuple(p) for p in rec.get("_overruled") or ()}:
            reason = "overruled_piece"
        elif rec.get("goi_class") in not_counted_classes:
            reason = rec["goi_class"]
        else:
            reason = rec.get("confidence") or "not_counted"
        out.append({"leaf": leaf, "reason": reason,
                    "record": "" if rec is None else f"{rec['genome']} {rec['chrom']}:{rec['start']}-{rec['end']}"})
    return out


# The GOI rows of the figures that the report's records are made of: HIGH, MEDIUM and AMBIGUOUS (a
# LOW model is in no record).
FIGURE_RECORD_LABELS = ("HIGH", "MEDIUM", "AMBIGUOUS")


def figure_marks(goi_annotations, records, locus, demoted_classes, dropped=()):
    """What the figures of home locus ``locus`` need to draw the report's records once each (item 36):
    the outcome of each GOI model of that locus's per-genome GFFs (the search's rows labelled HIGH,
    MEDIUM or AMBIGUOUS; the rescue passes' rows are no input of the figures) where the report has
    one, and the other home loci that found a record drawn there.

    A record is drawn by one model in each home locus whose search modelled it, its shown member
    (_shown_member, the rule of the trees), taken among the record's models of that locus that are
    no fragment (ModelStatus=fragment, which the default figures hide) where it has one: a counted
    record whose representative is a fragment inside a wider model of the same gene would
    otherwise be in no default figure, its other models merged and the fragment hidden.
    Outcomes: 'merged' (the record's other models of this locus: a second model of the same gene,
    a wider model that a compact one replaced, the other half of one reading frame, a piece of the
    whole gene that holds the record), 'dropped' (a model of a record the report dropped,
    ``dropped``: (record, why) pairs of drop_genome_wide_strong_hits and drop_gene_order_only); for
    a shown member whose record the report does not count: 'overruled_piece' (the record is a
    whole gene and the shown member its piece), the record's demoted class (``demoted_classes``:
    paralog_not_goi, identity_coverage_decoupled), 'genome_wide_strong_hit' (a kept answer of the
    genome-wide path); 'not_in_report' for a model no record holds. The figures leave out 'merged'
    and 'dropped' (bin/plot_synteny.py, REPORT_LEFT_OUT) and draw the others by their label.
    Returns {'models': [{'genome', 'model', 'outcome', 'record', 'drawn_by'}], 'also_found_from':
    [{'genome', 'model', 'loci'}]}, each sorted by genome and model id; ``record`` is the span
    ('chrom:start-end') of the model's record, '' for not_in_report; ``drawn_by`` is the model that
    draws the record of a merged model (its shown member), '' for every other outcome."""
    ids, fragments = defaultdict(set), defaultdict(set)
    for a in goi_annotations or ():
        if (a.get("source") == locus and a.get("mrna_id") and not a.get("rescue_pass")
                and a.get("confidence") in FIGURE_RECORD_LABELS):
            ids[a["genome"]].add(a["mrna_id"])
            if (a.get("model_status") or "").strip().lower() == "fragment":
                fragments[a["genome"]].add(a["mrna_id"])

    def span(r):
        return f"{r['chrom']}:{r['start']}-{r['end']}"
    marks, also, held = {}, [], set()
    for r, _why in dropped or ():
        for m in r.get("_members") or ():
            if m[0] == locus and m[1] in ids.get(r["genome"], ()):
                marks[(r["genome"], m[1])] = ("dropped", span(r), "")
                held.add((r["genome"], m[1]))
    for r in records or ():
        mine = ids.get(r["genome"], set())
        held |= {(r["genome"], m[1]) for m in r.get("_members") or () if m[0] == locus and m[1] in mine}
        shown = (_shown_member(r, locus, mine - fragments.get(r["genome"], set()))
                 or _shown_member(r, locus, mine))
        if shown is None:
            continue
        for m in r.get("_members") or ():
            if m[0] == locus and m[1] in mine and m[1] != shown:
                marks[(r["genome"], m[1])] = ("merged", span(r), shown)
        if (locus, shown) in {tuple(p) for p in r.get("_overruled") or ()}:
            outcome = "overruled_piece"
        elif r.get("genome_wide_strong_hit") and r.get("goi_class") not in demoted_classes:
            outcome = "genome_wide_strong_hit"
        elif r.get("goi_class") in demoted_classes:
            outcome = r["goi_class"]
        else:
            outcome = ""
        if outcome:
            marks[(r["genome"], shown)] = (outcome, span(r), "")
        others = sorted(set(r.get("provenance") or ()) - {locus})
        if others:
            also.append({"genome": r["genome"], "model": shown, "loci": others})
    for genome in ids:
        for mid in ids[genome]:
            if (genome, mid) not in held:
                marks[(genome, mid)] = ("not_in_report", "", "")
    return {"models": [{"genome": g, "model": mid, "outcome": o, "record": where, "drawn_by": by}
                       for (g, mid), (o, where, by) in sorted(marks.items())],
            "also_found_from": sorted(also, key=lambda e: (e["genome"], e["model"]))}


def load_tree_infos(folder):
    """{locus: info} of COMPUTE_TREE's <locus>_tree.info.json (fix-phase d6): seed, threads,
    sequences and the bootstrap (status ultrafast, none_fewer_than_4_sequences,
    dropped_after_failure, placeholder_tree; replicates; reason). {} without any."""
    out = {}
    if not folder or not os.path.isdir(folder):
        return out
    for name in sorted(os.listdir(folder)):
        if not name.endswith("_tree.info.json"):
            continue
        try:
            with open(os.path.join(folder, name)) as fh:
                info = json.load(fh)
        except (OSError, ValueError):
            continue
        if isinstance(info, dict):
            out[info.get("locus") or name[:-len("_tree.info.json")]] = info
    return out


def load_home_goi(path):
    """{query_id, gene_id, transcript_id, method, found_by} of ANNOTATE_GOI's goi_info.json (the keys
    it holds), or None without a file or a gene."""
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            info = json.load(fh)
    except (OSError, ValueError):
        return None
    if not isinstance(info, dict) or not info.get("gene_id"):
        return None
    return {k: info[k] for k in ("query_id", "gene_id", "transcript_id", "method", "found_by") if info.get(k)}


def load_goi_info_block(path, key):
    """The block ``key`` of ANNOTATE_GOI's goi_info.json, or None (no file, no such block, or a
    file written before 2026-11-02)."""
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            info = json.load(fh)
    except (OSError, ValueError):
        return None
    block = info.get(key) if isinstance(info, dict) else None
    return dict(block) if isinstance(block, dict) else None


def load_home_annotation(folder):
    """Where the home annotation of a run without a home GFF came from (borrow_annotations.py,
    batch HG-5b: home_annotation.json in the home protein database folder), with a warning:
    the home proteome then holds only the predicted genes near the query's hits and the gene
    models transferred from one annotated target, so the specificity terms compare a model
    with these alone. None with a home GFF (no such file)."""
    path = os.path.join(folder, "home_annotation.json") if folder else ""
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            info = json.load(fh)
    except (OSError, ValueError):
        return None
    if not isinstance(info, dict):
        return None
    out = {k: info.get(k) for k in ("source", "partial", "target", "target_choice", "target_region",
                                    "borrowed", "predicted", "predictions_dropped", "proteins")}
    out["genes_left_out_over_goi"] = len(info.get("left_out_over_goi") or [])
    borrowed = (f" and {out['borrowed']} gene model(s) transferred from {out['target']}"
                if out.get("borrowed") else "")
    out["warning"] = ("The home genome has no annotation: its proteome is partial, the "
                      f"{out.get('predicted') or 0} predicted gene(s) near the query's hits{borrowed}. "
                      "BestOtherHome and RealCompetitor compare the calls with these proteins only.")
    return out


def load_run_settings(path):
    """What decided the run's settings (effective_params.json of RESOLVE_EFFECTIVE_PARAMS):
    the preset or the preset modules that were applied and why, the triggers that fire, the
    resulting settings and the features the triggers read, the design quantities the module
    rules read (profile), how each value that is not the default came about (derivation) and
    the constraint checks. None without the file (or for a sentinel)."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path) as fh:
            text = fh.read().strip()
        if not text or text.lower().startswith("sentinel"):
            return None
        data = json.loads(text)
    except (OSError, ValueError):
        return None
    out = {k: data.get(k) for k in RUN_SETTINGS_KEYS}
    # A record without modules_mode predates --modules (2026-10-15): such a run was legacy. The default
    # of a run is auto since 2026-11-10, and every record since 2026-10-15 names its mode.
    out["modules_mode"] = out["modules_mode"] or "legacy"
    out["modules_applied"] = out["modules_applied"] or []
    out["modules_would_fire"] = out["modules_would_fire"] or []
    # What the module triggers say about the query and its home (reporting only; the
    # trigger is the one definition, select_modules.py).
    feats = out.get("features") or {}
    fired = {d.get("id") for d in out["modules_would_fire"]}
    contig_end = "C2" in fired
    out["flags"] = {
        "domain_shared": "A3" in fired,
        "dispersed_family_loci": (feats.get("loci") or []) if "A2" in fired else [],
        "home_contig_end": contig_end,
        "home_subtelomeric": contig_end and bool(feats.get("telomeric_repeat_at_end")),
        "dist_to_scaffold_end": feats.get("dist_to_scaffold_end"),
    }
    # The resources the run used (fix-phase d6): absent in a file written before 2026-11-02.
    if data.get("run_resources"):
        out["run_resources"] = data["run_resources"]
    # The tools' versions (fix-phase d8): absent before 2026-11-02.
    if data.get("tool_versions"):
        out["tool_versions"] = data["tool_versions"]
    return out


def main():
    parser = argparse.ArgumentParser(description="Generate SynVoy final evidence report")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--qc_json", help="Path to QC summary JSON")
    parser.add_argument("--qc_policy", default="unspecified", help="QC handling policy used in the workflow")
    parser.add_argument("--output", required=True, help="Report JSON")
    parser.add_argument(
        "--locus_ownership_tiebreak_gap", type=float, default=10.0,
        help="Panel score margin a paralog needs over the GOI's best entry to take a call; below "
             "it the call stays the GOI's, and between two entries that are both the GOI or both "
             "not, the flanking-synteny tiebreak decides",
    )
    parser.add_argument(
        "--phylo_placement_dir", default=None,
        help="Directory of PHYLO_PLACEMENT_CHECK TSVs (the phylogenetic placement check, "
             "docs/ALGORITHM.md 10). Defaults to <results_dir>/phylo_placement.",
    )
    parser.add_argument(
        "--disable_coverage_demotion", action="store_true",
        help="Keep high-identity/low-coverage GOI calls in the headline counts "
             "(listed only, as in runs made before 2026-08-22).",
    )
    parser.add_argument(
        "--identity_decoupled_min_identity", type=float, default=50.0,
        help="Identity at or above which a low-coverage GOI call is demoted.",
    )
    parser.add_argument(
        "--identity_decoupled_max_qcov", type=float, default=0.35,
        help="Query coverage below which a high-identity GOI call is demoted.",
    )
    parser.add_argument(
        "--strong_flanking_min", type=int, default=20,
        help="HIGH flanking models of a (genome, locus, scaffold) with no HIGH GOI call at "
             "which self_consistency flags strong_synteny_no_goi (the pipeline passes "
             "strong_synteny_no_goi_min_flanking: 2 x n_flanking_genes by default).",
    )
    parser.add_argument(
        "--effective_params", default=None,
        help="effective_params.json of RESOLVE_EFFECTIVE_PARAMS, written into the report as "
             "run_settings: which preset or preset modules ran, why, and on which features.",
    )
    parser.add_argument(
        "--home_db_dir", default=None,
        help="The home protein database folder of the search. Without a home annotation it "
             "holds home_annotation.json (borrow_annotations.py): the report says where the "
             "home annotation and proteome came from and that the proteome is partial.",
    )
    parser.add_argument(
        "--goi_info", default=None,
        help="goi_info.json of ANNOTATE_GOI: what it records about the home genome and the GOI's "
             "home gene goes into the report (docs/OUTPUT.md).",
    )
    parser.add_argument(
        "--tree_info_dir", default=None,
        help="The <locus>_tree.info.json files of COMPUTE_TREE: how each tree was made (seed, "
             "threads, the bootstrap, also one dropped after a failure), written as trees.",
    )
    parser.add_argument(
        "--sort_provenance_dir", default=None,
        help="PHYLO_SORT's sort_provenance.json files (any depth): how the target genomes were "
             "ordered, and the NCBI lookups that did not answer, written as genome_order.",
    )
    parser.add_argument(
        "--download_record", default=None,
        help="download_record.tsv of easy mode's genome download: the genomes and GFFs that did "
             "not arrive are named in warnings (downloads).",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Do not exit non-zero when the staging directories contain no GFFs, scores, or hits.",
    )
    args = parser.parse_args()

    report = build_report(
        args.results_dir,
        qc_json=args.qc_json,
        qc_policy=args.qc_policy,
        locus_ownership_tiebreak_gap=args.locus_ownership_tiebreak_gap,
        coverage_demotion=not args.disable_coverage_demotion,
        phylo_placement_dir=args.phylo_placement_dir,
        identity_decoupled_min_identity=args.identity_decoupled_min_identity,
        identity_decoupled_max_qcov=args.identity_decoupled_max_qcov,
        strong_flanking_min=args.strong_flanking_min,
        core_completeness=synvoy_settings.pick(
            None, synvoy_settings.load_settings(args.effective_params), "core_completeness", False,
            synvoy_settings.flag)[0],
        tree_infos=load_tree_infos(args.tree_info_dir),
    )
    report["run_settings"] = load_run_settings(args.effective_params)
    sw_backend = report.pop("smith_waterman", None)
    if isinstance(report["run_settings"], dict) and sw_backend is not None:
        report["run_settings"]["smith_waterman"] = sw_backend
    home_annotation = load_home_annotation(args.home_db_dir)
    if home_annotation:
        report["home_annotation"] = home_annotation
        report.setdefault("warnings", []).append(home_annotation["warning"])
    home_filter = load_goi_info_block(args.goi_info, "home_filter")
    if home_filter is not None:
        report["home_genome_filter"] = home_filter
    order = load_sort_provenance(args.sort_provenance_dir)
    if order:
        report["genome_order"] = order
        for o in order:
            if o.get("network_failures"):
                report.setdefault("warnings", []).append(
                    f"NCBI did not answer while the target genomes were ordered ({o['network_failures'][0]}"
                    + (f"; {len(o['network_failures']) - 1} more" if len(o["network_failures"]) > 1 else "")
                    + f"): they were ordered by {o.get('strategy')}. The order of the search waves, and "
                    "so the seeds, can differ from a run where NCBI answered (genome_order).")
    downloads = load_download_record(args.download_record)
    if downloads:
        report["downloads"] = downloads
        failed = [d["accession"] for d in downloads if d["genome"] == "failed"]
        no_gff = [d["accession"] for d in downloads if d["genome"] == "ok" and d["gff"] == "none"]
        if failed:
            report.setdefault("warnings", []).append(
                f"{len(failed)} genome(s) the run chose could not be downloaded and were not searched: "
                f"{', '.join(failed)} (downloads).")
        if no_gff:
            report.setdefault("warnings", []).append(
                f"{len(no_gff)} downloaded genome(s) came without a GFF: {', '.join(no_gff)}; their genes "
                "carry no target annotation (downloads).")
    trees = report.get("trees") or {}
    if trees:
        for locus, t in trees.items():
            boot = t.get("bootstrap") or {}
            if boot.get("status") == "dropped_after_failure":
                report.setdefault("warnings", []).append(
                    f"The {locus} tree has no branch support: {boot.get('reason', '')} (trees).")
            threads = t.get("threads")
            if isinstance(threads, int) and threads > 1 and boot.get("status") != "placeholder_tree":
                report.setdefault("warnings", []).append(
                    f"The {locus} tree was inferred with {threads} IQ-TREE threads: a rerun can give "
                    "other branch lengths and supports (--compute_tree_cpus 1 gives the same tree) (trees).")
    uniprot = load_goi_info_block(args.goi_info, "uniprot_lookup")
    if uniprot and uniprot.get("status") == "failed":
        # Its names can decide which home gene is the GOI (fix-phase d7).
        report.setdefault("warnings", []).append(
            f"UniProt did not answer for {uniprot.get('query_id')} ({uniprot.get('error', '')}): the home "
            "gene of the GOI was chosen without the query's UniProt names (goi_info.json, "
            "uniprot_lookup); with an answer another home gene can be chosen.")
    # The home gene ANNOTATE_GOI matched to the query: the figures break a tie between home genes
    # (or transcripts) under the query's hits with it (fix item 28).
    home_goi = load_home_goi(args.goi_info)
    if home_goi:
        report["home_goi"] = home_goi
    unmatched = load_goi_info_block(args.goi_info, "home_gff_unmatched")
    if unmatched:
        report["home_gff_unmatched"] = unmatched
        report.setdefault("warnings", []).append(
            f"{unmatched.get('seqids')} of the home GFF's {unmatched.get('of')} sequences with genes are "
            f"not in the home FASTA (e.g. {', '.join(unmatched.get('examples', [])[:3])}): their genes "
            "are no flanking genes and no competitors of a call (home_gff_unmatched).")
    slot_units = copies_in_slot_by_genome(report.get("goi_dedup", {}).get("records", []),
                                          (report["run_settings"] or {}).get("home_copies"))
    if slot_units:
        report["copies_in_slot"] = slot_units
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report generated: {args.output}")
    # Surface the GOI headline so the run's terminal/log shows the real result rather
    # than only the (often-zero) raw-hit diagnostic. See docs/TODO.md §1k.
    print(f"  {report['summary']['headline']}")
    decoy = report.get("decoy_rate")
    if decoy and decoy["cases"]:
        med = decoy["medium_or_better"]
        print(f"  Decoy calls MEDIUM or better: {med['calls']} of {decoy['cases']} decoy cases "
              f"(upper end of the exact 95 % interval: {100 * med['ci95_upper']:.2f} %)")
    gw_decoy = (decoy or {}).get("genome_wide") or {}
    if gw_decoy.get("cases"):
        ans = gw_decoy["answers"]
        print(f"  Decoy answers of the genome-wide path: {ans['calls']} of {gw_decoy['cases']} decoy "
              f"cases (upper end of the exact 95 % interval: {100 * ans['ci95_upper']:.2f} %)")
    for warning in report.get("warnings", []):
        print(f"  [!] {warning}")

    # §1.1 honesty: loudly caveat a gene-family query so a user doesn't read paralog
    # calls as confident orthologs.
    advisory = report.get("gene_family_advisory", {})
    if advisory.get("is_gene_family_query"):
        print("  [!] GENE-FAMILY QUERY — " + advisory["advice"])
        for reason in advisory["reasons"]:
            print(f"      · {reason}")

    if report["staging_diagnostics"]["empty"]:
        msg = format_empty_staging_message(report["staging_diagnostics"])
        if args.allow_empty:
            print(f"Warning (--allow-empty set): {msg}", file=sys.stderr)
        else:
            print(msg, file=sys.stderr)
            sys.exit(2)


if __name__ == "__main__":
    main()
