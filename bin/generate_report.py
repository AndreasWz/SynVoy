#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import re
from collections import Counter, defaultdict


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


def _parse_gff_attrs(attr_field):
    attrs = {}
    for kv in (attr_field or "").split(";"):
        if "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def _safe_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def _is_true(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def canonical_genome_id(path_or_name):
    name = os.path.basename(path_or_name or "")
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
                    role = (attrs.get("SynTerraRole") or "").strip().lower()
                    if not role:
                        role = "goi" if model_id.startswith("GOI_") or attrs.get("SynTerra_Parent", "").startswith("GOI_") else "flanking"
                    confidence = (attrs.get("Confidence", "") or "UNKNOWN").upper()
                    goi_class = attrs.get("GOIClass", "")
                    evidence_type = attrs.get("EvidenceType", attrs.get("Type", "")) or "unknown"

                    stats["total_annotations"] += 1
                    stats["role_counts"][role] += 1
                    stats["evidence_type_counts"][evidence_type] += 1
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
        "summary": {
            "total_new_genes": sum(genes_added_per_genome.values()),
            "genomes_with_hits": len(hits_per_genome),
            "total_hits": sum(hits_per_genome.values()),
            "genomes_with_annotations": len(annotation_summary["per_genome"]),
            "total_annotations": annotation_summary["total_annotations"],
            "total_goi_annotations": annotation_summary["role_counts"].get("goi", 0),
            "ambiguous_goi_annotations": annotation_summary["goi_class_counts"].get("ambiguous_goi_family_member", 0),
            "fallback_goi_annotations": annotation_summary["fallback_goi_annotations"],
            "low_confidence_regions": region_summary["confidence_counts"].get("LOW", 0),
            "goi_absent_genomes": annotation_summary["genomes_without_goi"],
            "goi_ambiguous_only_genomes": annotation_summary["genomes_with_only_ambiguous_goi"],
            "failed_qc_genomes_with_downstream_results": sorted(set(failed_downstream)),
        },
    }
    return report


def main():
    parser = argparse.ArgumentParser(description="Generate SynTerra final evidence report")
    parser.add_argument("--results_dir", required=True)
    parser.add_argument("--qc_json", help="Path to QC summary JSON")
    parser.add_argument("--qc_policy", default="unspecified", help="QC handling policy used in the workflow")
    parser.add_argument("--output", required=True, help="Report JSON")
    args = parser.parse_args()

    report = build_report(args.results_dir, qc_json=args.qc_json, qc_policy=args.qc_policy)
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report generated: {args.output}")


if __name__ == "__main__":
    main()
