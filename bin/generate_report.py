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
    ".m8",
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


def summarize_annotations(gff_files):
    per_genome = {}
    role_counts = Counter()
    goi_confidence_counts = Counter()
    goi_class_counts = Counter()
    evidence_type_counts = Counter()
    model_status_counts = Counter()
    goi_evidence_counts = Counter()
    fallback_goi_annotations = 0
    total_annotations = 0

    for gff_path in gff_files:
        genome = canonical_genome_id(gff_path)
        stats = per_genome.setdefault(
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
                    line = line.strip()
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

                    stats["total_annotations"] += 1
                    stats["role_counts"][role] += 1
                    stats["evidence_type_counts"][evidence_type] += 1
                    if model_status:
                        stats["model_status_counts"][model_status] += 1
                        model_status_counts[model_status] += 1
                    stats["goi_confidence_counts"][confidence] += 1 if role == "goi" else 0
                    if role == "goi" and goi_class:
                        stats["goi_class_counts"][goi_class] += 1

                    total_annotations += 1
                    role_counts[role] += 1
                    evidence_type_counts[evidence_type] += 1

                    if role == "goi":
                        stats["goi_annotations"] += 1
                        goi_confidence_counts[confidence] += 1
                        if goi_class:
                            goi_class_counts[goi_class] += 1
                        goi_evidence_counts[evidence_type] += 1
                        if goi_class == "ambiguous_goi_family_member":
                            stats["ambiguous_goi_annotations"] += 1
                        else:
                            stats["resolved_goi_annotations"] += 1
                        if evidence_type in {"fallback_hit_span", "raw_hit", "rescued_exon"}:
                            fallback_goi_annotations += 1
        except Exception as exc:
            print(f"Warning: Could not parse GFF {gff_path}: {exc}")

    per_genome_list = []
    genomes_without_goi = []
    genomes_with_only_ambiguous_goi = []
    for genome in sorted(per_genome):
        stats = per_genome[genome]
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

    return {
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

    The genome is canonicalised (staging prefix stripped) so the same target genome
    found from several home loci groups together; the staging tag is kept as the
    per-locus source for dedup provenance (docs/TODO.md §1h).
    """
    annotations = []
    for gff_path in gff_files:
        genome = canonical_genome_id(gff_path)
        source = staging_source_label(gff_path)
        try:
            with open(gff_path) as fh:
                for line in fh:
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
                    if role != "goi":
                        continue
                    try:
                        start, end = int(parts[3]), int(parts[4])
                    except ValueError:
                        continue
                    if end < start:
                        start, end = end, start
                    annotations.append({
                        "genome": genome,
                        "source": source,
                        "chrom": parts[0],
                        "start": start,
                        "end": end,
                        "confidence": (attrs.get("Confidence", "") or "UNKNOWN").upper(),
                        "goi_class": attrs.get("GOIClass", ""),
                        "identity": (attrs.get("Identity", "") or "").strip(),
                        "target_gene": attrs.get("TargetGene", ""),
                        "model_status": attrs.get("ModelStatus", ""),
                        "mrna_id": model_id,
                    })
        except Exception as exc:
            print(f"Warning: Could not parse GFF for dedup {gff_path}: {exc}")
    return annotations


def collect_high_flanking_per_locus(gff_files):
    """Count distinct HIGH-confidence flanking genes per (genome, locus, chrom) from
    the staged per-locus GFFs. Used by build_self_consistency to spot blocks with
    strong conserved synteny but no GOI model (docs/TODO.md §1j identity-decay sanity).
    De-duplicated by gene Name so a multi-exon annotation and its coarse hit-span
    don't double-count one flanking gene.
    """
    seen = defaultdict(set)  # (genome, locus, chrom) -> set of gene names
    for gff_path in gff_files:
        genome = canonical_genome_id(gff_path)
        locus = staging_source_label(gff_path)
        try:
            with open(gff_path) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 9 or parts[2] not in {"mRNA", "gene"}:
                        continue
                    attrs = _parse_gff_attrs(parts[8])
                    role = (attrs.get("SynVoyRole") or "").strip().lower()
                    if role != "flanking":
                        continue
                    if (attrs.get("Confidence") or "").strip().upper() != "HIGH":
                        continue
                    gene = (attrs.get("Name") or attrs.get("SynVoy_Parent")
                            or attrs.get("ID") or f"{parts[0]}:{parts[3]}")
                    seen[(genome, locus, parts[0])].add(gene)
        except Exception as exc:
            print(f"Warning: Could not parse GFF for self-consistency {gff_path}: {exc}")
    return {key: len(genes) for key, genes in seen.items()}


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


def build_self_consistency(goi_annotations, goi_dedup, flanking_per_locus,
                           *, strong_flanking_min=20):
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

    Reciprocal-best protein check across home paralogs (TODO §1j first bullet) is
    deferred — it needs an aligner and the per-locus home queries staged at report
    time, but for paralogs on different chromosomes (TP53/63/73, the motivating case)
    synteny anchoring already provides the locus assignment. Same-coord cross-locus
    cases are caught by §1h's coord-overlap dedup.
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

    return {
        "checks_performed": ["strong_synteny_no_goi", "cross_locus_duplicate"],
        "deferred_checks": ["reciprocal_best_paralog"],
        "thresholds": {"strong_flanking_min": strong_flanking_min},
        "flags": flags,
        "summary": {
            "n_strong_synteny_no_goi": sum(1 for f in flags if f["type"] == "strong_synteny_no_goi"),
            "n_cross_locus_duplicates": sum(1 for f in flags if f["type"] == "cross_locus_duplicate"),
            "total_flags": len(flags),
        },
    }


def build_report(results_dir, qc_json=None, qc_policy=None):
    qc_records = load_qc_records(qc_json)
    qc_summary = summarize_qc(qc_records)

    regions_dir = os.path.join(results_dir, "regions")
    hits_dir = os.path.join(results_dir, "hits")
    scores_dir = os.path.join(results_dir, "scores")

    fasta_files = glob.glob(os.path.join(regions_dir, "*.faa")) + glob.glob(os.path.join(regions_dir, "*.fna"))
    gff_files = glob.glob(os.path.join(regions_dir, "*.gff")) + glob.glob(os.path.join(regions_dir, "*.gff3"))
    score_files = glob.glob(os.path.join(scores_dir, "*.scores.tsv"))
    hit_files = glob.glob(os.path.join(hits_dir, "*.m8"))

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
    annotation_summary = summarize_annotations(gff_files)
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
    goi_annotations = collect_goi_annotations(gff_files)
    goi_dedup = dedupe_goi_annotations(goi_annotations)
    dedup_high = goi_dedup["high_confidence_goi_deduped"]
    dedup_medium = goi_dedup["medium_confidence_goi_deduped"]

    # End-of-run self-consistency checks (docs/TODO.md §1j): identity-decay sanity +
    # cross-locus duplicate visibility. The flanking-only-block surfacing is the §1e
    # piece that already lives in cluster_grs.py.
    flanking_per_locus = collect_high_flanking_per_locus(gff_files)
    self_consistency = build_self_consistency(goi_annotations, goi_dedup, flanking_per_locus)

    report = {
        "genome_qc": qc_records,
        "qc_summary": {
            **qc_summary,
            "qc_fail_policy": qc_policy or "unspecified",
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
        },
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
            "cross_locus_duplicate_goi": goi_dedup["cross_locus_duplicates"],
            "goi_hits_collapsed_by_dedup": goi_dedup["hits_collapsed_by_dedup"],
            "self_consistency_flag_count": self_consistency["summary"]["total_flags"],
            "strong_synteny_no_goi_flags": self_consistency["summary"]["n_strong_synteny_no_goi"],
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
        "--allow-empty",
        action="store_true",
        help="Do not exit non-zero when the staging directories contain no GFFs, scores, or hits.",
    )
    args = parser.parse_args()

    report = build_report(args.results_dir, qc_json=args.qc_json, qc_policy=args.qc_policy)
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report generated: {args.output}")
    # Surface the GOI headline so the run's terminal/log shows the real result rather
    # than only the (often-zero) raw-hit diagnostic. See docs/TODO.md §1k.
    print(f"  {report['summary']['headline']}")

    if report["staging_diagnostics"]["empty"]:
        msg = format_empty_staging_message(report["staging_diagnostics"])
        if args.allow_empty:
            print(f"Warning (--allow-empty set): {msg}", file=sys.stderr)
        else:
            print(msg, file=sys.stderr)
            sys.exit(2)


if __name__ == "__main__":
    main()
