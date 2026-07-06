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


def summarize_hits(hit_files):
    hits_per_genome = {}
    for hit_path in hit_files:
        genome = canonical_genome_id(hit_path)
        count = 0
        try:
            with open(hit_path) as hit_fh:
                count = sum(1 for line in hit_fh if line.strip())
        except Exception:
            count = 0
        hits_per_genome[genome] = hits_per_genome.get(genome, 0) + count
    return hits_per_genome


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
    # ── collect_goi_annotations accumulator ───────────────────────────────
    goi_records = []
    # ── collect_high_flanking_per_locus accumulator ───────────────────────
    flanking_seen = defaultdict(set)  # (genome, locus, chrom) -> set of gene names

    for gff_path in gff_files:
        genome = canonical_genome_id(gff_path)
        source = staging_source_label(gff_path)
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
            },
        )

        try:
            with open(gff_path) as fh:
                for line in fh:
                    # NB: rstrip("\n") matches the collect_goi_annotations behaviour
                    # which trimmed only the newline (not surrounding whitespace).
                    # summarize_annotations used .strip() but parts of the GFF aren't
                    # whitespace-sensitive past the 9-tab field split, so the looser
                    # rstrip preserves identical semantics for the line.split("\t")
                    # downstream.
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 9 or parts[2] not in {"mRNA", "gene"}:
                        continue
                    attrs = _parse_gff_attrs(parts[8])
                    model_id = attrs.get("ID", "")
                    role = (attrs.get("SynVoyRole") or "").strip().lower()
                    if not role:
                        role = "goi" if model_id.startswith("GOI_") or attrs.get("SynVoy_Parent", "").startswith("GOI_") else "flanking"
                    confidence = (attrs.get("Confidence", "") or "UNKNOWN").upper()
                    goi_class = attrs.get("GOIClass", "")
                    evidence_type = attrs.get("EvidenceType", attrs.get("Type", "")) or "unknown"
                    model_status = attrs.get("ModelStatus", "")

                    # — summarize_annotations bookkeeping —
                    ann_stats["total_annotations"] += 1
                    ann_stats["role_counts"][role] += 1
                    ann_stats["evidence_type_counts"][evidence_type] += 1
                    if model_status:
                        ann_stats["model_status_counts"][model_status] += 1
                        model_status_counts[model_status] += 1
                    ann_stats["goi_confidence_counts"][confidence] += 1 if role == "goi" else 0
                    if role == "goi" and goi_class:
                        ann_stats["goi_class_counts"][goi_class] += 1

                    total_annotations += 1
                    role_counts[role] += 1
                    evidence_type_counts[evidence_type] += 1

                    if role == "goi":
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
                            "chrom": parts[0],
                            "start": start,
                            "end": end,
                            "confidence": confidence,
                            "goi_class": goi_class,
                            "identity": (attrs.get("Identity", "") or "").strip(),
                            "query_cov": (attrs.get("QueryCoverage", "") or "").strip(),
                            "target_gene": attrs.get("TargetGene", ""),
                            "model_status": model_status,
                            "mrna_id": model_id,
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
        }
        per_genome_list.append(row)
        if stats["goi_annotations"] == 0:
            genomes_without_goi.append(genome)
        elif stats["resolved_goi_annotations"] == 0:
            genomes_with_only_ambiguous_goi.append(genome)

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
        "genomes_with_only_ambiguous_goi": genomes_with_only_ambiguous_goi,
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
            matches.extend(glob.glob(os.path.join(dir_path, pat)))
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


def _format_goi_headline(high, medium, low, n_genomes):
    """One-line, unambiguous summary of the GOI ortholog yield, for the summary block.

    Replaces the old reliance on `total_hits` (raw .m8 count, often 0) as the at-a-glance
    signal — that made successful runs look empty. See docs/TODO.md §1k.
    """
    if not (high or medium or low):
        return "No GOI ortholog annotations were produced."
    confident = []
    if high:
        confident.append(f"{high} high-confidence")
    if medium:
        confident.append(f"{medium} medium-confidence")
    if confident:
        lead = " + ".join(confident)
        tail = f" (+{low} low-confidence/ambiguous)" if low else ""
    else:
        lead, tail = f"{low} low-confidence/ambiguous", ""
    return f"{lead} GOI ortholog annotation(s){tail} across {n_genomes} genome(s)."


def _confidence_rank(confidence):
    return {"HIGH": 2, "MEDIUM": 1}.get((confidence or "").upper(), 0)


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


def dedupe_goi_annotations(annotations, min_overlap=0.8):
    """Collapse HIGH/MEDIUM GOI hits that are the same target gene found from
    multiple home loci (docs/TODO.md §1h).

    Key: ``(genome, target_chrom)`` plus reciprocal coordinate overlap > ``min_overlap``.
    The kept record is the highest-confidence / highest-identity member, with a
    ``provenance`` list of the source home loci it was recovered from. When the
    identity is bit-for-bit identical across >=2 loci (a near-certain "one gene,
    many seeds" signal) the record is reclassified ``goi_class=cross_locus_duplicate``.

    LOW-confidence hits are left untouched (they are noisy fallback spans, not
    ortholog calls), matching the headline metric which only counts HIGH/MEDIUM.
    """
    considered = [a for a in annotations if a["confidence"] in {"HIGH", "MEDIUM"}]
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
                if _reciprocal_overlap(a, cl["rep"]) > min_overlap:
                    cl["members"].append(a)
                    better = (_confidence_rank(a["confidence"]), _safe_float(a["identity"]))
                    current = (_confidence_rank(cl["rep"]["confidence"]), _safe_float(cl["rep"]["identity"]))
                    if better > current:
                        cl["rep"] = a
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
                "goi_class": "cross_locus_duplicate" if is_cross_dup else rep["goi_class"],
                "original_goi_class": rep["goi_class"],
                "target_gene": rep["target_gene"],
                "model_status": rep["model_status"],
                "provenance": sources,
                "n_source_loci": len(sources),
                "n_merged_hits": len(members),
                "cross_locus_duplicate": is_cross_dup,
            })

    records.sort(key=lambda r: (-_confidence_rank(r["confidence"]), -_safe_float(r["identity"]), r["genome"], r["chrom"]))
    return {
        "records": records,
        "high_confidence_goi_deduped": sum(1 for r in records if r["confidence"] == "HIGH"),
        "medium_confidence_goi_deduped": sum(1 for r in records if r["confidence"] == "MEDIUM"),
        "cross_locus_duplicates": sum(1 for r in records if r["cross_locus_duplicate"]),
        "hits_collapsed_by_dedup": len(considered) - len(records),
        "pre_dedup_high_medium": len(considered),
        "post_dedup_records": len(records),
    }


def load_paralog_check_rows(paralog_check_dir, suffix=".paralog_check.tsv"):
    """Read every ``<genome><suffix>`` TSV produced by RECIPROCAL_BEST_PARALOG
    (docs/TODO.md §1j Phase B) — also reused for the §1m locus-ownership TSVs,
    which share the same column layout but carry ``locus_<id>|<gene>`` panel
    labels in ``best_paralog`` (pass ``suffix='.locus_ownership.tsv'``).

    Returns a list of dicts keyed by the TSV header. Skips header-only files (the
    rescue/single-paralog no-op case) and silently tolerates missing dirs."""
    rows = []
    if not paralog_check_dir or not os.path.isdir(paralog_check_dir):
        return rows
    with os.scandir(paralog_check_dir) as it:
        for entry in it:
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
                        except ValueError:
                            continue
                        rows.append(row)
            except OSError:
                continue
    return rows


def build_gene_family_advisory(goi_class_counts, n_home_loci=0,
                               family_ratio=0.4, family_min_ambiguous=5,
                               family_min_loci=3):
    """Detect a gene-family / multi-paralog query and advise honest interpretation
    (docs/TODO.md §1.1). Purely data-driven from signals the report already has —
    no search-side change. SynVoy is validated for divergent LOW-COPY genes via
    conserved flanking synteny; inside a large paralog family ortholog-vs-paralog
    assignment is unreliable, so we say so instead of asserting the GOI name.

    Triggers (either):
      - a large share of GOI calls are ``ambiguous_goi_family_member`` (the class the
        classifier already uses for "this is in the family but I can't pin the ortholog"), or
      - the query fanned out to many home loci (paralog cluster in the home genome).
    """
    counts = goi_class_counts or {}
    ambiguous = counts.get("ambiguous_goi_family_member", 0)
    total = sum(counts.values())
    reasons = []
    if total and ambiguous >= family_min_ambiguous and ambiguous >= family_ratio * total:
        pct = round(100 * ambiguous / total)
        reasons.append(
            f"{ambiguous}/{total} GOI calls ({pct}%) are ambiguous family members")
    if n_home_loci >= family_min_loci:
        reasons.append(f"the query resolved to {n_home_loci} home loci (paralog cluster)")
    is_family = bool(reasons)
    return {
        "is_gene_family_query": is_family,
        "reasons": reasons,
        "advice": (
            "Gene-family / multi-paralog query detected. SynVoy is validated for "
            "divergent LOW-COPY genes via conserved flanking synteny; within a large "
            "paralog family, ortholog-vs-paralog assignment is unreliable. Treat GOI "
            "calls as CANDIDATES needing manual curation and weight the §1m ownership / "
            "§1j self-consistency flags over raw identity." if is_family else ""
        ),
    }


def build_paralog_confusion_flags(paralog_rows, min_gap):
    """Emit a ``paralog_confusion`` flag for each call whose best home paralog
    differs from the modal best paralog at the same (genome, locus), provided
    the bitscore gap to the runner-up exceeds ``min_gap`` (so near-ties don't
    fire).

    Per-call summary (``best_paralog_per_call``) is also returned for transparency
    so the report carries the alignment evidence behind every flag.
    """
    if not paralog_rows:
        return [], [], {}

    # Modal best-paralog per (genome, locus).
    by_locus = defaultdict(Counter)
    for r in paralog_rows:
        by_locus[(r["genome"], r["locus"])][r["best_paralog"]] += 1
    modal_per_locus = {key: ctr.most_common(1)[0][0] for key, ctr in by_locus.items()}

    flags = []
    per_call_summary = []
    for r in paralog_rows:
        key = (r["genome"], r["locus"])
        modal = modal_per_locus.get(key, "")
        is_confused = (r["best_paralog"] != modal) and (r["gap_f"] >= min_gap)
        per_call_summary.append({
            "genome": r["genome"],
            "locus": r["locus"],
            "chrom": r["chrom"],
            "start": r["start_i"],
            "end": r["end_i"],
            "mrna_id": r.get("mrna_id", ""),
            "best_paralog": r["best_paralog"],
            "best_bit": r["best_bit_f"],
            "second_paralog": r.get("second_paralog", ""),
            "second_bit": r["second_bit_f"],
            "bitscore_gap": r["gap_f"],
            "locus_modal_paralog": modal,
            "confused": is_confused,
        })
        if is_confused:
            flags.append({
                "type": "paralog_confusion",
                "genome": r["genome"],
                "locus": r["locus"],
                "chrom": r["chrom"],
                "start": r["start_i"],
                "end": r["end_i"],
                "mrna_id": r.get("mrna_id", ""),
                "assigned_paralog": modal,           # the locus's "own" paralog (modal best)
                "best_match_paralog": r["best_paralog"],
                "bitscore_gap": r["gap_f"],
                "advice": (
                    f"Recovered protein aligns best to '{r['best_paralog']}' (bit={r['best_bit_f']:.1f}) "
                    f"but the modal best paralog at this locus is '{modal}' "
                    f"(bit-gap {r['gap_f']:.1f} >= {min_gap}). Likely a sibling paralog "
                    "dragged in by MEDIUM-confidence iterative search — re-check the call's "
                    "family assignment before reporting it as an ortholog of the locus paralog."
                ),
            })

    return flags, per_call_summary, modal_per_locus


def load_panel_meta(panel_meta_path):
    """Read ``home_paralog_panel.faa.meta.tsv`` (build_home_paralog_panel.py, §1m).

    Returns {panel_id: {locus, gene, chrom, start, end, is_goi}} so the report can
    map an ownership row's ``best_paralog`` (a ``locus_<id>|<gene>`` panel label)
    back to its home locus / gene / GOI-ness. Tolerates a missing file (ownership
    becomes label-only / advisory)."""
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
                }
    except OSError:
        pass
    return meta


def _panel_locus(panel_label):
    """``locus_4|gene-OMD`` -> ``locus_4``. Gene names may contain '-' but never
    '|', so splitting on the first pipe is safe."""
    return panel_label.split("|", 1)[0] if "|" in panel_label else panel_label


def build_locus_ownership(goi_dedup, ownership_rows, panel_meta, flanking_per_locus,
                          *, tiebreak_gap=15.0):
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

    Mutates the records in ``goi_dedup['records']`` in place (adds owning_* keys) and
    returns (flags, summary)."""
    records = goi_dedup.get("records", [])
    if not records or not ownership_rows:
        return [], {"n_reattributed": 0, "n_paralog_misassignment": 0,
                    "n_owner_search_fumble": 0, "n_tiebreak_applied": 0, "evaluated": 0}

    # Index ownership rows by genome+chrom for coordinate matching.
    rows_by_gc = defaultdict(list)
    for r in ownership_rows:
        rows_by_gc[(r["genome"], r["chrom"])].append(r)

    flags = []
    n_reattributed = n_paralog = n_fumble = n_tiebreak = n_eval = 0

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
        owning_is_goi = panel_meta.get(best_label, {}).get("is_goi", False)

        # --- Synteny tiebreak: when the top-2 panel genes are near-tied, defer to
        # which candidate locus's HIGH flanking neighbourhood is present on this
        # target chromosome (true local synteny beats a coin-flip alignment gap).
        tiebroken = False
        if second_label and gap < tiebreak_gap:
            second_locus = _panel_locus(second_label)
            if second_locus != owning_locus:
                f_best = flanking_per_locus.get((rec["genome"], owning_locus, rec["chrom"]), 0)
                f_second = flanking_per_locus.get((rec["genome"], second_locus, rec["chrom"]), 0)
                if f_second > f_best:
                    owning_locus = second_locus
                    owning_gene = second_label.split("|", 1)[1] if "|" in second_label else second_label
                    owning_is_goi = panel_meta.get(second_label, {}).get("is_goi", False)
                    tiebroken = True
                    n_tiebreak += 1

        provenance = rec.get("provenance", [])
        # Re-attribution only makes sense when the owning gene IS the GOI (moving a true
        # GOI ortholog to its correct home locus). When the best match is a home PARALOG
        # (owning_is_goi=False — e.g. a homology-family entry like home_gene-BGN), that's a
        # paralog_misassignment, not a reattribution to a searched locus.
        reattributed = owning_is_goi and (owning_locus not in provenance or (
            len(provenance) > 1 and provenance[0] != owning_locus
        ))
        owner_absent = owning_locus not in provenance

        rec["owning_locus"] = owning_locus
        rec["owning_gene"] = owning_gene
        rec["owning_is_goi_gene"] = owning_is_goi
        rec["ownership_bit"] = round(row["best_bit_f"], 1)
        rec["ownership_gap"] = round(gap, 1)
        rec["ownership_tiebroken"] = tiebroken
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
                    "GOI-search block formed there (docs/TODO.md §1m search-core fix)."
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
    }
    return flags, summary


def build_self_consistency(goi_annotations, goi_dedup, flanking_per_locus,
                           *, strong_flanking_min=20,
                           paralog_flags=None, paralog_per_call=None,
                           paralog_modal_per_locus=None,
                           paralog_min_gap=5.0,
                           ownership_flags=None, ownership_summary=None,
                           identity_decoupled_min_identity=50.0,
                           identity_decoupled_max_qcov=0.35):
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

    * ``paralog_confusion`` (§1j Phase B, landed) — a call whose best home paralog
      (parasail SW vs every sequence in the multi-FASTA query) differs from the
      modal best at the same (genome, locus) and the bitscore gap to the runner-up
      exceeds ``paralog_min_gap``. Catches TP53/63/73-style sibling leakage at
      MEDIUM confidence (the residual issue from docs/TP53_PARALOG_ASSESSMENT.md).

    For single-paralog queries the §1j check is a no-op and ``paralog_flags`` is
    empty; for multi-paralog runs the RECIPROCAL_BEST_PARALOG process produces a
    TSV per (locus, genome) that ``load_paralog_check_rows`` then aggregates.
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
                    "Strong synteny but no HIGH-confidence GOI call — try a lower "
                    "--classify_high_min_identity / --classify_medium_min_identity, or "
                    "inspect this block manually (the orthologous neighbourhood likely "
                    "exists but the GOI is too divergent for current thresholds)."
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
                "note": "same gene identically modeled from multiple home loci; collapsed into one record (§1h)",
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
        qc_raw = (a.get("query_cov") or "").strip()
        qc = _safe_float(qc_raw, None) if qc_raw else None
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
                "confidence": a.get("confidence"),
                "advice": (
                    "High %identity but missing/low query coverage — likely a short "
                    "high-identity local window, not a full-length match. Treat the "
                    "identity as unreliable and inspect the alignment span/coverage."
                ),
            })

    # §1j Phase B paralog confusion flags (computed in build_paralog_confusion_flags
    # and threaded through; we only attach them here so the existing flag-summary
    # plumbing carries them).
    paralog_flags = paralog_flags or []
    paralog_per_call = paralog_per_call or []
    paralog_modal_per_locus = paralog_modal_per_locus or {}
    flags.extend(paralog_flags)

    # §1m locus-ownership flags (locus_reattributed / goi_owner_search_fumble /
    # paralog_misassignment) computed in build_locus_ownership and threaded through.
    ownership_flags = ownership_flags or []
    ownership_summary = ownership_summary or {}
    flags.extend(ownership_flags)

    checks_performed = ["strong_synteny_no_goi", "cross_locus_duplicate",
                        "identity_coverage_decoupled"]
    deferred_checks = []
    if paralog_per_call:
        checks_performed.append("paralog_confusion")
    else:
        deferred_checks.append("paralog_confusion")
    if ownership_summary.get("evaluated"):
        checks_performed.append("locus_ownership")
    else:
        deferred_checks.append("locus_ownership")

    return {
        "checks_performed": checks_performed,
        "deferred_checks": deferred_checks,
        "thresholds": {
            "strong_flanking_min": strong_flanking_min,
            "paralog_confusion_min_gap": paralog_min_gap,
            "identity_decoupled_min_identity": identity_decoupled_min_identity,
            "identity_decoupled_max_qcov": identity_decoupled_max_qcov,
        },
        "flags": flags,
        # Carry the full per-call paralog alignment table so the paper can
        # cite exact numbers without re-running anything.
        "paralog_alignment_table": paralog_per_call,
        "modal_paralog_per_locus": [
            {"genome": g, "locus": l, "paralog": p}
            for (g, l), p in sorted(paralog_modal_per_locus.items())
        ],
        "locus_ownership_summary": ownership_summary,
        "summary": {
            "n_strong_synteny_no_goi": sum(1 for f in flags if f["type"] == "strong_synteny_no_goi"),
            "n_cross_locus_duplicates": sum(1 for f in flags if f["type"] == "cross_locus_duplicate"),
            "n_paralog_confusion": sum(1 for f in flags if f["type"] == "paralog_confusion"),
            "n_locus_reattributed": sum(1 for f in flags if f["type"] == "locus_reattributed"),
            "n_goi_owner_search_fumble": sum(1 for f in flags if f["type"] == "goi_owner_search_fumble"),
            "n_paralog_misassignment": sum(1 for f in flags if f["type"] == "paralog_misassignment"),
            "n_identity_coverage_decoupled": sum(1 for f in flags if f["type"] == "identity_coverage_decoupled"),
            "total_flags": len(flags),
        },
    }


def build_report(results_dir, qc_json=None, qc_policy=None, paralog_confusion_min_gap=5.0,
                 locus_ownership_tiebreak_gap=10.0):
    qc_records = load_qc_records(qc_json)
    qc_summary = summarize_qc(qc_records)

    regions_dir = os.path.join(results_dir, "regions")
    hits_dir = os.path.join(results_dir, "hits")
    scores_dir = os.path.join(results_dir, "scores")

    # Single scandir() per directory, classify by suffix in one pass (docs/TODO.md §1c).
    regions_scan = scan_dir_by_suffix(regions_dir, {
        "fasta": (".faa", ".fna"),
        "gff": (".gff", ".gff3"),
    })
    scores_scan = scan_dir_by_suffix(scores_dir, {"scores": (".scores.tsv",)})
    hits_scan = scan_dir_by_suffix(hits_dir, {"hits": (".m8",)})
    fasta_files = regions_scan["fasta"]
    gff_files = regions_scan["gff"]
    score_files = scores_scan["scores"]
    hit_files = hits_scan["hits"]

    staging_diagnostics = build_staging_diagnostics(
        results_dir,
        dir_patterns={
            "regions_fasta": (regions_dir, ["*.faa", "*.fna"]),
            "regions_gff": (regions_dir, ["*.gff", "*.gff3"]),
            "scores": (scores_dir, ["*.scores.tsv"]),
            "hits": (hits_dir, ["*.m8"]),
        },
        match_counts={
            "fasta_files": len(fasta_files),
            "gff_files": len(gff_files),
            "score_files": len(score_files),
            "hit_files": len(hit_files),
        },
    )

    genes_added_per_genome = summarize_fasta_outputs(fasta_files)
    hits_per_genome = summarize_hits(hit_files)
    # One read-pass over every GFF instead of three (docs/TODO.md §1c). The three
    # public functions (summarize_annotations / collect_goi_annotations /
    # collect_high_flanking_per_locus) still work for tests + external callers; they
    # just throw away the two outputs they don't need each.
    annotation_summary, goi_annotations, flanking_per_locus = _process_gffs_unified(gff_files)
    region_summary = summarize_region_scores(score_files)

    downstream_genomes = set(genes_added_per_genome) | set(hits_per_genome)
    downstream_genomes |= {row["genome"] for row in annotation_summary["per_genome"]}
    downstream_genomes |= {row["genome"] for row in region_summary["per_genome"]}

    failed_downstream = []
    for rec in qc_summary["failed_genomes"]:
        genome = rec["genome"]
        if genome in downstream_genomes:
            failed_downstream.append(genome)

    # GOI ortholog headline metrics — what a user actually wants to read. The historical
    # total_hits/genomes_with_hits counted raw .m8 search hits staged under hits/, which is
    # frequently 0 even on successful runs (annotations come from the region GFFs, not the
    # m8 staging) and was being misread as "the run found nothing". See docs/TODO.md §1k.
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
    goi_dedup = dedupe_goi_annotations(goi_annotations)
    dedup_high = goi_dedup["high_confidence_goi_deduped"]
    dedup_medium = goi_dedup["medium_confidence_goi_deduped"]

    # §1j Phase B: per-call reciprocal-best paralog check (run by RECIPROCAL_BEST_PARALOG
    # per (locus, genome) and staged under paralog_check/). Aggregated here so the
    # self-consistency block can emit paralog_confusion flags for any call whose best
    # home paralog differs from the modal best at the same locus.
    paralog_check_rows = load_paralog_check_rows(
        os.path.join(results_dir, "paralog_check")
    )
    paralog_flags, paralog_per_call, paralog_modal = build_paralog_confusion_flags(
        paralog_check_rows, min_gap=paralog_confusion_min_gap,
    )

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

    # §1m: ownership may relabel a call goi_class=paralog_not_goi (its best home match is a
    # paralog, not the GOI — e.g. a 55% biglycan hit mislabeled GOI_DCN). Such a call is NOT
    # a GOI ortholog, so drop it from the headline HIGH/MEDIUM counts. It stays in
    # goi_dedup["records"] with its flag for transparency; the pre-ownership counts are kept
    # as *_pre_ownership. This is what makes the RBH check actually demote the false positive.
    dedup_high_pre_ownership, dedup_medium_pre_ownership = dedup_high, dedup_medium
    _own_recs = goi_dedup.get("records", [])
    dedup_high = sum(1 for r in _own_recs
                     if r.get("confidence") == "HIGH" and r.get("goi_class") != "paralog_not_goi")
    dedup_medium = sum(1 for r in _own_recs
                       if r.get("confidence") == "MEDIUM" and r.get("goi_class") != "paralog_not_goi")

    # End-of-run self-consistency checks (docs/TODO.md §1j): identity-decay sanity +
    # cross-locus duplicate visibility + paralog confusion + locus ownership (§1m).
    # The flanking-only-block surfacing is the §1e piece that lives in cluster_grs.py.
    self_consistency = build_self_consistency(
        goi_annotations, goi_dedup, flanking_per_locus,
        paralog_flags=paralog_flags,
        paralog_per_call=paralog_per_call,
        paralog_modal_per_locus=paralog_modal,
        paralog_min_gap=paralog_confusion_min_gap,
        ownership_flags=ownership_flags,
        ownership_summary=ownership_summary,
    )

    # §1.1 honesty: detect a gene-family / multi-paralog query and advise treating GOI
    # calls as candidates (SynVoy is validated for divergent low-copy genes, not large
    # paralog families). Report-only, driven by the ambiguous-call share + home-locus count.
    gene_family_advisory = build_gene_family_advisory(
        annotation_summary.get("goi_class_counts", {}),
        n_home_loci=len(flanking_per_locus or {}),
    )

    report = {
        "genome_qc": qc_records,
        "qc_summary": {
            **qc_summary,
            "qc_fail_policy": qc_policy or "unspecified",
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
        },
        "gene_family_advisory": gene_family_advisory,
        "synteny_results": {
            "genes_discovered": genes_added_per_genome,
            "synteny_hits_count": hits_per_genome,
        },
        "annotations": annotation_summary,
        "regions": region_summary,
        "goi_dedup": goi_dedup,
        "self_consistency": self_consistency,
        "staging_diagnostics": staging_diagnostics,
        "summary": {
            # --- Headline: the number(s) a user actually cares about ---
            # Post-dedup (docs/TODO.md §1h): distinct HIGH/MEDIUM ortholog genes, not the
            # per-home-locus seed count, so the same gene found from N loci reads as one.
            "headline": _format_goi_headline(dedup_high, dedup_medium, goi_low, genomes_with_goi),
            "headline_metric": dedup_high,  # distinct high-confidence GOI orthologs (post-dedup)
            "high_confidence_goi": dedup_high,
            "medium_confidence_goi": dedup_medium,
            "high_confidence_goi_pre_dedup": goi_high,
            "medium_confidence_goi_pre_dedup": goi_medium,
            # Pre-ownership = before §1m relabeled paralog_not_goi calls out of the headline.
            "high_confidence_goi_pre_ownership": dedup_high_pre_ownership,
            "medium_confidence_goi_pre_ownership": dedup_medium_pre_ownership,
            "cross_locus_duplicate_goi": goi_dedup["cross_locus_duplicates"],
            "goi_hits_collapsed_by_dedup": goi_dedup["hits_collapsed_by_dedup"],
            "self_consistency_flag_count": self_consistency["summary"]["total_flags"],
            "strong_synteny_no_goi_flags": self_consistency["summary"]["n_strong_synteny_no_goi"],
            "paralog_confusion_flags": self_consistency["summary"]["n_paralog_confusion"],
            "locus_reattributed_orthologs": self_consistency["summary"]["n_locus_reattributed"],
            "goi_owner_search_fumbles": self_consistency["summary"]["n_goi_owner_search_fumble"],
            "paralog_misassignments": self_consistency["summary"]["n_paralog_misassignment"],
            "gene_family_query": gene_family_advisory["is_gene_family_query"],
            "low_confidence_goi": goi_low,
            "confident_goi_annotations": confident_goi,
            "resolved_goi_annotations": resolved_goi,
            "genomes_with_goi_annotations": genomes_with_goi,
            "total_new_genes": sum(genes_added_per_genome.values()),
            "genomes_with_annotations": len(annotation_summary["per_genome"]),
            "total_annotations": annotation_summary["total_annotations"],
            "total_goi_annotations": annotation_summary["role_counts"].get("goi", 0),
            "ambiguous_goi_annotations": annotation_summary["goi_class_counts"].get("ambiguous_goi_family_member", 0),
            "fallback_goi_annotations": annotation_summary["fallback_goi_annotations"],
            "low_confidence_regions": region_summary["confidence_counts"].get("LOW", 0),
            "goi_absent_genomes": annotation_summary["genomes_without_goi"],
            "goi_ambiguous_only_genomes": annotation_summary["genomes_with_only_ambiguous_goi"],
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
            "staging_empty": staging_diagnostics["empty"],
            # --- Low-level diagnostic (NOT the result count) ---
            # Raw MMseqs/BLAST .m8 hits staged under hits/; often 0 even when annotations
            # were produced. Renamed from total_hits/genomes_with_hits so it stops reading
            # as "the run found nothing" (docs/TODO.md §1k).
            "total_raw_search_hits": sum(hits_per_genome.values()),
            "genomes_with_raw_search_hits": len(hits_per_genome),
        },
    }
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


def main():
    parser = argparse.ArgumentParser(description="Generate SynVoy final evidence report")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--qc_json", help="Path to QC summary JSON")
    parser.add_argument("--qc_policy", default="unspecified", help="QC handling policy used in the workflow")
    parser.add_argument("--output", required=True, help="Report JSON")
    parser.add_argument(
        "--paralog_confusion_min_gap", type=float, default=5.0,
        help="Bitscore gap required to flag a paralog_confusion (docs/TODO.md §1j Phase B)",
    )
    parser.add_argument(
        "--locus_ownership_tiebreak_gap", type=float, default=15.0,
        help="Panel bitscore gap below which locus ownership defers to the "
             "flanking-synteny tiebreak (docs/TODO.md §1m)",
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
        paralog_confusion_min_gap=args.paralog_confusion_min_gap,
        locus_ownership_tiebreak_gap=args.locus_ownership_tiebreak_gap,
    )
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report generated: {args.output}")
    # Surface the GOI headline so the run's terminal/log shows the real result rather
    # than only the (often-zero) raw-hit diagnostic. See docs/TODO.md §1k.
    print(f"  {report['summary']['headline']}")

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
