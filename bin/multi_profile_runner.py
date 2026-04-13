#!/usr/bin/env python3
"""
multi_profile_runner.py — Run SynVoy with multiple parameter profiles and select best.

For small searches (where total_work = loci × targets × profiles ≤ max_jobs),
generates several parameter profile variants from the LLM-estimated baseline,
scores the results of each, and selects the best-performing profile.

Scoring criteria:
  1. GOI hit rate (% of targets with HIGH/MEDIUM confidence GOI)
  2. Average synteny score across regions
  3. Gene model completeness (complete vs fragment ratio)
  4. Noise ratio (fewer LOW-confidence results = cleaner)
"""

import argparse
import copy
import json
import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Profile generation
# ---------------------------------------------------------------------------

# Adjustments relative to the LLM-estimated baseline.
# Each dict maps param_name → (operation, operand).
# Operations: "mul" (multiply), "add", "set", "min_of", "max_of"
PROFILE_ADJUSTMENTS = {
    "sensitive": {
        "mmseqs_sensitivity": ("add", 1.0),
        "min_hit_identity": ("mul", 0.5),
        "min_synteny_score": ("mul", 0.6),
        "search_evalue": ("mul", 100.0),
        "sw_min_score": ("mul", 0.6),
        "min_gene_identity": ("mul", 0.6),
        "enable_smith_waterman": ("set", True),
        "min_hit_length": ("max_of", 5),
    },
    "stringent": {
        "min_synteny_score": ("mul", 1.3),
        "min_hit_identity": ("mul", 1.5),
        "search_evalue": ("mul", 0.01),
        "min_gene_identity": ("mul", 1.3),
        "min_hit_length": ("mul", 1.5),
    },
}

# Parameter specs for clamping (subset of llm_param_advisor.ESTIMABLE_PARAMS)
PARAM_RANGES = {
    "mmseqs_sensitivity": (1.0, 12.0),
    "min_hit_identity": (0.0, 100.0),
    "min_synteny_score": (0.1, 1.0),
    "search_evalue": (1e-10, 100.0),
    "sw_min_score": (5.0, 200.0),
    "min_gene_identity": (5.0, 100.0),
    "min_hit_length": (5, 500),
}


def _apply_adjustment(base_value: Any, operation: str, operand: Any) -> Any:
    """Apply a single profile adjustment to a parameter value."""
    if operation == "set":
        return operand
    if operation == "mul":
        return type(base_value)(base_value * operand)
    if operation == "add":
        return type(base_value)(base_value + operand)
    if operation == "min_of":
        return min(base_value, operand)
    if operation == "max_of":
        return max(base_value, operand)
    return base_value


def generate_profiles(
    base_params: Dict[str, Any], defaults: Dict[str, Any]
) -> Dict[str, Dict[str, Any]]:
    """
    Generate parameter profiles from LLM-estimated baseline.

    Returns dict mapping profile_name → full param set (base + adjustments).
    The 'balanced' profile is the LLM estimate unchanged.
    """
    profiles = {"balanced": dict(base_params)}

    for profile_name, adjustments in PROFILE_ADJUSTMENTS.items():
        profile = dict(base_params)
        for param, (op, operand) in adjustments.items():
            current = profile.get(param, defaults.get(param, 0))
            adjusted = _apply_adjustment(current, op, operand)

            # Clamp to valid range
            if param in PARAM_RANGES:
                lo, hi = PARAM_RANGES[param]
                if isinstance(adjusted, (int, float)):
                    adjusted = max(lo, min(hi, adjusted))
                    if isinstance(current, int):
                        adjusted = int(adjusted)

            profile[param] = adjusted

        profiles[profile_name] = profile

    return profiles


# ---------------------------------------------------------------------------
# Result scoring
# ---------------------------------------------------------------------------


def _parse_gff_confidence_stats(gff_dir: str) -> Dict[str, int]:
    """
    Parse GFF files to count confidence levels of GOI models.
    Returns {"HIGH": n, "MEDIUM": n, "LOW": n, "total": n}.
    """
    counts = defaultdict(int)
    if not os.path.isdir(gff_dir):
        return dict(counts)

    for gff_file in Path(gff_dir).glob("**/*.gff"):
        try:
            with open(gff_file) as fh:
                for line in fh:
                    if line.startswith("#") or "\t" not in line:
                        continue
                    if "GOI_" not in line and "SynVoy_Parent=GOI" not in line:
                        continue
                    # Extract Confidence attribute
                    m = re.search(r"Confidence=(\w+)", line)
                    if m:
                        counts[m.group(1)] += 1
                        counts["total"] += 1
        except Exception:
            pass

    return dict(counts)


def _parse_model_status_stats(gff_dir: str) -> Dict[str, int]:
    """Count ModelStatus values (complete/partial/fragment) from GFF files."""
    counts = defaultdict(int)
    if not os.path.isdir(gff_dir):
        return dict(counts)

    for gff_file in Path(gff_dir).glob("**/*.gff"):
        try:
            with open(gff_file) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    m = re.search(r"ModelStatus=(\w+)", line)
                    if m:
                        counts[m.group(1)] += 1
        except Exception:
            pass

    return dict(counts)


def _parse_region_scores(regions_dir: str) -> List[float]:
    """Extract synteny scores from region BED files."""
    scores = []
    if not os.path.isdir(regions_dir):
        return scores

    for bed_file in Path(regions_dir).glob("**/*.bed"):
        try:
            with open(bed_file) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        try:
                            scores.append(float(parts[4]))
                        except ValueError:
                            # Score might be in the name field (Reg1_G5_CHIGH_S0.85)
                            name = parts[3] if len(parts) > 3 else ""
                            m = re.search(r"_S([\d.]+)", name)
                            if m:
                                scores.append(float(m.group(1)))
        except Exception:
            pass

    return scores


def score_result(result_dir: str, n_targets: int) -> Dict[str, Any]:
    """
    Score a pipeline result set by multiple quality metrics.

    Returns a dict with individual scores and a composite score.
    """
    metrics = {
        "goi_hit_rate": 0.0,
        "avg_synteny_score": 0.0,
        "completeness_ratio": 0.0,
        "noise_ratio": 1.0,
        "composite_score": 0.0,
    }

    if not os.path.isdir(result_dir):
        return metrics

    # GFF confidence stats
    gff_dir = os.path.join(result_dir, "iterative_results", "regions")
    conf_stats = _parse_gff_confidence_stats(gff_dir)
    high = conf_stats.get("HIGH", 0)
    medium = conf_stats.get("MEDIUM", 0)
    low = conf_stats.get("LOW", 0)
    total_goi = conf_stats.get("total", 0)

    # GOI hit rate: fraction of targets with HIGH or MEDIUM GOI
    if n_targets > 0:
        # Approximate: count unique GFF files with GOI hits
        goi_genomes = set()
        if os.path.isdir(gff_dir):
            for gff_file in Path(gff_dir).glob("*.gff"):
                try:
                    with open(gff_file) as fh:
                        for line in fh:
                            if "GOI_" in line and ("Confidence=HIGH" in line or "Confidence=MEDIUM" in line):
                                goi_genomes.add(gff_file.stem)
                                break
                except Exception:
                    pass
        metrics["goi_hit_rate"] = len(goi_genomes) / n_targets

    # Average synteny score
    regions_dir = os.path.join(result_dir, "regions")
    if not os.path.isdir(regions_dir):
        regions_dir = os.path.join(result_dir, "iterative_results", "regions")
    scores = _parse_region_scores(regions_dir)
    if scores:
        metrics["avg_synteny_score"] = sum(scores) / len(scores)

    # Model completeness
    model_stats = _parse_model_status_stats(gff_dir)
    complete = model_stats.get("complete", 0)
    partial = model_stats.get("partial", 0)
    fragment = model_stats.get("fragment", 0)
    total_models = complete + partial + fragment
    if total_models > 0:
        metrics["completeness_ratio"] = (complete + 0.5 * partial) / total_models

    # Noise ratio (lower is better → we invert it)
    if total_goi > 0:
        signal = high + medium
        metrics["noise_ratio"] = low / total_goi
    else:
        metrics["noise_ratio"] = 1.0

    # Composite score (weighted combination)
    metrics["composite_score"] = (
        0.40 * metrics["goi_hit_rate"]
        + 0.25 * metrics["avg_synteny_score"]
        + 0.20 * metrics["completeness_ratio"]
        + 0.15 * (1.0 - metrics["noise_ratio"])
    )

    return metrics


# ---------------------------------------------------------------------------
# Decision logic
# ---------------------------------------------------------------------------


def should_multi_profile(
    n_loci: int,
    n_targets: int,
    n_profiles: int = 3,
    max_total_jobs: int = 30,
) -> bool:
    """
    Determine whether multi-profile mode should be used.
    Based on total work = loci × targets × profiles.
    """
    total_work = n_loci * n_targets * n_profiles
    decision = total_work <= max_total_jobs
    logger.info(
        f"Multi-profile check: {n_loci} loci × {n_targets} targets × {n_profiles} profiles "
        f"= {total_work} jobs (max={max_total_jobs}) → {'YES' if decision else 'NO'}"
    )
    return decision


def select_best_profile(
    profile_scores: Dict[str, Dict[str, Any]],
) -> Tuple[str, Dict[str, Any]]:
    """
    Select the best-performing profile based on composite scores.
    Returns (profile_name, metrics).
    """
    if not profile_scores:
        return "balanced", {}

    best_name = max(profile_scores, key=lambda k: profile_scores[k].get("composite_score", 0))
    best_metrics = profile_scores[best_name]

    logger.info("Profile comparison:")
    for name, metrics in sorted(profile_scores.items()):
        marker = " ★" if name == best_name else ""
        logger.info(
            f"  {name}: composite={metrics.get('composite_score', 0):.3f} "
            f"(GOI={metrics.get('goi_hit_rate', 0):.2f}, "
            f"synteny={metrics.get('avg_synteny_score', 0):.2f}, "
            f"complete={metrics.get('completeness_ratio', 0):.2f}, "
            f"noise={metrics.get('noise_ratio', 0):.2f}){marker}"
        )

    return best_name, best_metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Generate multi-profile parameter sets or score completed results"
    )
    subparsers = parser.add_subparsers(dest="command")

    # Generate profiles
    gen = subparsers.add_parser("generate", help="Generate profile param sets from baseline")
    gen.add_argument("--baseline", required=True, help="Path to estimated_params.json")
    gen.add_argument("--output_dir", required=True, help="Directory to write profile JSONs")

    # Score results
    score = subparsers.add_parser("score", help="Score completed profile results")
    score.add_argument("--results_dir", required=True, help="Base directory with profile results")
    score.add_argument("--n_targets", type=int, required=True, help="Number of target genomes")
    score.add_argument("--output", required=True, help="Output JSON with scores and best profile")

    # Check if multi-profile should be used
    check = subparsers.add_parser("check", help="Check if multi-profile is appropriate")
    check.add_argument("--n_loci", type=int, required=True)
    check.add_argument("--n_targets", type=int, required=True)
    check.add_argument("--max_jobs", type=int, default=30)

    args = parser.parse_args()

    if args.command == "generate":
        with open(args.baseline) as fh:
            baseline = json.load(fh)

        base_params = baseline.get("parameters", {})
        # Use ESTIMABLE_PARAMS defaults
        try:
            from llm_param_advisor import ESTIMABLE_PARAMS
            defaults = {k: v["default"] for k, v in ESTIMABLE_PARAMS.items()}
        except ImportError:
            defaults = {}

        profiles = generate_profiles(base_params, defaults)

        outdir = Path(args.output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        for name, params in profiles.items():
            out_path = outdir / f"profile_{name}.json"
            with open(out_path, "w") as fh:
                json.dump({"profile": name, "parameters": params}, fh, indent=2)
            logger.info(f"Written profile '{name}' → {out_path}")

    elif args.command == "score":
        profile_scores = {}
        results_base = Path(args.results_dir)

        for profile_dir in results_base.iterdir():
            if profile_dir.is_dir() and profile_dir.name.startswith("profile_"):
                profile_name = profile_dir.name.replace("profile_", "")
                metrics = score_result(str(profile_dir), args.n_targets)
                profile_scores[profile_name] = metrics

        best_name, best_metrics = select_best_profile(profile_scores)

        result = {
            "best_profile": best_name,
            "best_metrics": best_metrics,
            "all_scores": profile_scores,
        }

        with open(args.output, "w") as fh:
            json.dump(result, fh, indent=2)

        logger.info(f"Best profile: {best_name} (composite={best_metrics.get('composite_score', 0):.3f})")

    elif args.command == "check":
        result = should_multi_profile(args.n_loci, args.n_targets, max_total_jobs=args.max_jobs)
        print(json.dumps({"should_multi_profile": result}))

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
