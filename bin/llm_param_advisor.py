#!/usr/bin/env python3
"""
llm_param_advisor.py — LLM-powered parameter estimation for SynVoy.

Uses Gemma 4 (via Ollama local API or Google Cloud Gemini API) to analyze
biological context and estimate optimal pipeline parameters.  Falls back
to deterministic heuristics when no LLM backend is available.

All estimated parameters are validated against allowed ranges and checked
for breaking combinations before being emitted.
"""

import argparse
import json
import logging
import os
import platform
import re
import socket
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

DEFAULT_OLLAMA_TIMEOUT = 300

# ---------------------------------------------------------------------------
# Parameter specifications: allowed ranges, types, and defaults
# ---------------------------------------------------------------------------

# Only Tier 1 + Tier 2 parameters are LLM-estimable.
# Tier 3 (classification thresholds, scoring weights, viz) are never touched.
ESTIMABLE_PARAMS: Dict[str, Dict[str, Any]] = {
    "max_intron": {"type": int, "min": 0, "max": 500000, "default": 20000},
    "cluster_distance": {"type": int, "min": 1000, "max": 2000000, "default": 150000},
    "n_flanking_genes": {"type": int, "min": 2, "max": 30, "default": 10},
    "min_synteny_score": {"type": float, "min": 0.1, "max": 1.0, "default": 0.6},
    "region_padding": {"type": int, "min": 5000, "max": 1000000, "default": 150000},
    "padding_min": {"type": int, "min": 5000, "max": 500000, "default": 50000},
    "padding_max": {"type": int, "min": 10000, "max": 1000000, "default": 200000},
    "search_evalue": {"type": float, "min": 1e-10, "max": 100.0, "default": 0.01},
    "min_hit_identity": {"type": float, "min": 0.0, "max": 100.0, "default": 10.0},
    "min_hit_length": {"type": int, "min": 5, "max": 500, "default": 10},
    "mmseqs_sensitivity": {"type": float, "min": 1.0, "max": 12.0, "default": 9.5},
    "max_flanking_goi_similarity": {"type": float, "min": 10.0, "max": 100.0, "default": 35.0},
    "max_flanking_distance": {"type": int, "min": 0, "max": 5000000, "default": 0},
    "expand_goi_similar": {"type": bool, "default": True},
    "expand_goi_similar_distance": {"type": int, "min": 10000, "max": 2000000, "default": 300000},
    "min_gene_identity": {"type": float, "min": 5.0, "max": 100.0, "default": 30.0},
    "enable_smith_waterman": {"type": bool, "default": True},
    "sw_min_score": {"type": float, "min": 5.0, "max": 200.0, "default": 20.0},
    "sw_min_identity": {"type": float, "min": 0.0, "max": 100.0, "default": 10.0},
    "enable_plm_search": {"type": bool, "default": False},
    "enable_structural_search": {"type": bool, "default": False},
    "max_blocks_per_genome": {"type": int, "min": 5, "max": 500, "default": 80},
    "min_block_genes": {"type": int, "min": 1, "max": 10, "default": 2},
    "max_consecutive_empty_blocks": {"type": int, "min": 5, "max": 100, "default": 25},
    "aug_relaxed_evalue_mult": {"type": float, "min": 10, "max": 100000, "default": 1000},
    "gap_search_window": {"type": int, "min": 5000, "max": 500000, "default": 50000},
    "prefer_large_genes": {"type": bool, "default": True},
    "min_flanking_size": {"type": int, "min": 100, "max": 5000, "default": 500},
    "exon_level_search": {"type": bool, "default": True},
}


# ---------------------------------------------------------------------------
# System prompt for Gemma 4
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """\
You are a genomic synteny analysis expert configuring SynVoy, a pipeline that \
finds orthologous genes across genomes.

HOW SYNVOY WORKS:
1. The user provides a query gene (GOI = Gene of Interest) and a home genome.
2. The pipeline locates the GOI in the home genome and extracts N flanking genes \
on each side as "synteny anchors" — these neighbouring genes form a conserved \
gene neighbourhood (synteny block).
3. For each target genome (searched in phylogenetic order, closest first): \
the flanking gene proteins are searched against the target using MMseqs2 \
(sequence homology). Hits that cluster within cluster_distance bp form a \
candidate synteny block.
4. min_synteny_score is the fraction of flanking anchors that must have hits \
in a candidate block. If too strict, real orthologs in rearranged genomes \
are missed. If too loose, false positive regions accumulate.
5. Within each candidate block, the pipeline runs gene prediction \
(Augustus for eukaryotes, Prodigal for prokaryotes) to find ORFs, then \
searches for the GOI itself via tblastn, Smith-Waterman, and optionally \
PLM embeddings (ProtT5) and structural comparison (ESMFold + Foldseek).
6. max_intron controls gene model construction — too small and multi-exon \
genes in large genomes are fragmented; too large and spurious gene merges occur.
7. cluster_distance and region_padding must scale with genome architecture: \
compact bacterial operons need small values, large vertebrate/plant genomes \
with huge intergenic regions need large values.

KEY TRADE-OFFS:
- Strict params (high identity, high synteny score) → high precision, low recall. \
Good for close species. Misses distant orthologs.
- Relaxed params (low identity, low synteny, high sensitivity) → high recall, \
lower precision. Needed for distant species. Slower, more false positives.
- Small peptides (<100 aa) produce weak alignment scores by nature — you MUST \
lower sw_min_score and min_hit_length, or the GOI will be missed entirely.
- Gene families (paralogs) require careful flanking gene filtering \
(max_flanking_goi_similarity) to avoid using paralogs as anchors.

CONTEXT JSON you will receive has these fields:
- query.length_aa: protein length — small (<100), normal (100-800), large (>800)
- query.estimated_exon_count: 1 for single-exon, >1 for multi-exon genes
- query.has_signal_peptide / is_secreted: secreted peptides are often small/fast-evolving
- query.gene_family_size_estimate: "small", "medium", or "large"
- query.domain_families: known protein domains (empty = novel/unknown)
- home_species.kingdom: Bacteria, Archaea, Fungi, Plantae, Animalia, or Unknown
- home_species.genome_size_mb: 0 means lookup failed — infer from kingdom
- home_species.gene_count: total annotated genes
- target_context.max_evolutionary_distance_mya: estimated divergence to most distant target
- target_context.kingdoms_represented: kingdoms in the target set
- target_context.target_count: number of target genomes

PARAMETER REFERENCE (only override parameters that should differ from defaults):

Genome Architecture:
- max_intron (int, default=20000): Max intron size in bp for gene models.
  Bacteria=0-500, Fungi=100-500, Insects=500-30000, Vertebrates=5000-100000, Plants=3000-100000+
- cluster_distance (int, default=150000): Max gap in bp to merge flanking gene hits into a block.
  Bacteria=10000-30000, Fungi=20000-50000, Compact animals=50000-150000, Vertebrates=150000-500000, Plants=300000-1000000
- n_flanking_genes (int, default=10): Synteny anchors per side of GOI.
  Gene-dense genomes=5-8, Normal=10, Gene-sparse or rearranged=12-20
- region_padding (int, default=150000): Extra bp around candidate blocks for gene prediction.
  Should be ~0.5-1× cluster_distance.
- padding_min/padding_max (int, default 50000/200000): Adaptive padding bounds.

Search Sensitivity:
- min_synteny_score (float, default=0.6): Fraction of flanking anchors required to call a block.
  Close species (same genus/family)=0.6-0.8, Moderate=0.4-0.6, Distant/rearranged=0.2-0.4
- search_evalue (float, default=0.01): MMseqs2 e-value cutoff.
  Close=0.001, Moderate=0.01-0.1, Distant=1.0-10
- min_hit_identity (float, default=10): Min alignment identity % for initial hits.
  Close=20-30, Moderate=10, Distant=5
- min_hit_length (int, default=10): Min alignment length (aa).
  Large proteins=30, Normal=10, Small peptides (<100aa)=8
- mmseqs_sensitivity (float, default=9.5): MMseqs2 sensitivity (1-12, higher=slower+more sensitive).
  Close=7-8, Default=9.5, Distant=10-11, Very distant (>500 Mya)=11-12
- min_gene_identity (float, default=30): Min identity % for flanking gene reciprocal best hits.
  Close=30-40, Moderate=20-30, Distant=10-20

Gene Family:
- max_flanking_goi_similarity (float, default=35): Exclude flanking genes >X% similar to GOI.
  Large gene families=20-25 (strict filter), Normal=35, Unique genes=50-100 (relaxed)
- expand_goi_similar (bool, default=true): Use GOI-like neighbors as extra search queries.
- max_flanking_distance (int, default=0): Max bp from GOI to walk for flanking genes. 0=unlimited.

Advanced Search:
- enable_plm_search (bool, default=false): ProtT5 protein language model embedding search. \
Finds remote homologs missed by sequence alignment. Enable for distance >400 Mya.
- enable_structural_search (bool, default=false): ESMFold structure prediction + Foldseek 3Di search. \
Finds structural homologs with no sequence similarity. Enable for distance >600 Mya. Requires GPU.
- sw_min_score (float, default=20): Min Smith-Waterman alignment score for GOI candidates.
  Small peptides=10-15, Normal=20, Large proteins=30
- sw_min_identity (float, default=10): Min Smith-Waterman identity %.

CRITICAL RULES:
1. Plants ALWAYS need: max_intron≥50000, cluster_distance≥300000, padding_max≥400000
2. Bacteria/Archaea: max_intron≤500, cluster_distance≤30000, n_flanking_genes≤8
3. Small peptides (<100 aa): MUST set sw_min_score≤15 and min_hit_length≤8
4. Cross-phylum searches (>500 Mya): mmseqs_sensitivity≥10, enable_plm_search=true
5. Never set min_synteny_score below 0.15 (generates too many false positives)
6. padding_max must be ≥ padding_min
7. When genome_size_mb=0 (lookup failed), infer genome architecture from kingdom: \
Bacteria ~5Mb, Fungi ~30Mb, Insects ~200-400Mb, Vertebrates ~1000-3000Mb, Plants ~500-5000Mb

REASONING: Think step-by-step before outputting JSON:
1. What kingdom? → sets genome architecture params (introns, distances)
2. How big is the query protein? → sets alignment thresholds
3. How far apart are the species? → sets search sensitivity and advanced methods
4. Gene family concerns? → sets flanking gene filtering
Then output ONLY the JSON override object.

OUTPUT FORMAT: A single JSON object with parameter overrides. Include ONLY \
parameters that differ from defaults. No explanations outside the JSON.

EXAMPLES:

Honeybee melittin (70aa secreted peptide) searching in bumblebees (~80 Mya, same family):
{"max_intron": 15000, "cluster_distance": 100000, "sw_min_score": 10, "min_hit_length": 8, "min_hit_identity": 8}

Human p53 (393aa tumor suppressor) searching across vertebrates (~450 Mya):
{"max_intron": 50000, "cluster_distance": 300000, "region_padding": 200000, "padding_max": 350000, "mmseqs_sensitivity": 10, "min_gene_identity": 15, "search_evalue": 0.1, "min_synteny_score": 0.4, "enable_plm_search": true}

Arabidopsis defensin (80aa) searching across angiosperms (~150 Mya):
{"max_intron": 80000, "cluster_distance": 400000, "region_padding": 250000, "padding_min": 100000, "padding_max": 500000, "sw_min_score": 12, "min_hit_length": 8, "min_synteny_score": 0.4, "n_flanking_genes": 12}

E. coli beta-lactamase (286aa) in Enterobacteriaceae (~300 Mya):
{"max_intron": 0, "cluster_distance": 20000, "region_padding": 15000, "padding_min": 5000, "padding_max": 25000, "n_flanking_genes": 5, "min_flanking_size": 200, "mmseqs_sensitivity": 10, "min_gene_identity": 15, "search_evalue": 0.1, "min_synteny_score": 0.4}
"""


# ---------------------------------------------------------------------------
# Backend: Ollama (local)
# ---------------------------------------------------------------------------


def _detect_system_resources() -> Dict[str, Any]:
    """Detect available system resources for model selection."""
    import psutil  # safe — installed with most Python envs

    ram_gb = psutil.virtual_memory().total / (1024**3)

    # Attempt GPU VRAM detection
    vram_gb = 0
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.total", "--format=csv,noheader,nounits"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            vram_gb = sum(int(x.strip()) for x in result.stdout.strip().split("\n") if x.strip()) / 1024
    except Exception:
        pass

    return {"ram_gb": round(ram_gb, 1), "vram_gb": round(vram_gb, 1)}


def _auto_select_model(resources: Optional[Dict] = None) -> str:
    """Select the best Gemma 4 model based on available resources."""
    if resources is None:
        try:
            resources = _detect_system_resources()
        except ImportError:
            # psutil not available; assume modest hardware
            return "gemma4:e4b"

    ram = resources.get("ram_gb", 8)
    vram = resources.get("vram_gb", 0)

    if vram >= 20 or ram >= 48:
        model = "gemma4:31b"
    elif vram >= 10 or ram >= 24:
        model = "gemma4:26b"
    else:
        model = "gemma4:e4b"

    logger.info(
        f"Auto-selected model: {model} "
        f"(RAM={ram:.0f}GB, VRAM={vram:.0f}GB)"
    )
    return model


def _check_ollama_available(ollama_url: str, model: str = "") -> bool:
    """Check if Ollama server is reachable and the model is available."""
    try:
        req = urllib.request.Request(f"{ollama_url}/api/version")
        with urllib.request.urlopen(req, timeout=3) as resp:
            if resp.status != 200:
                return False
    except Exception:
        return False

    # If a model name is provided, verify it's actually pulled
    if model:
        try:
            req = urllib.request.Request(f"{ollama_url}/api/tags")
            with urllib.request.urlopen(req, timeout=5) as resp:
                data = json.loads(resp.read().decode())
                available = {m.get("name", "") for m in data.get("models", [])}
                # Exact match first; if model has no tag, allow any variant of that base name
                if model not in available:
                    has_tag = ":" in model
                    if has_tag or not any(a.startswith(model + ":") or a == model for a in available):
                        logger.warning(
                            f"Ollama server reachable but model '{model}' not found. "
                            f"Available: {sorted(available) or 'none'}. "
                            f"Run 'ollama pull {model}' to download it."
                        )
                        return False
        except Exception as exc:
            logger.warning(f"Could not verify model availability: {exc}")
            # Server is reachable, proceed anyway — the model call will fail fast

    return True


def _call_ollama(
    context: Dict,
    model: str,
    ollama_url: str,
    timeout: int = DEFAULT_OLLAMA_TIMEOUT,
) -> Optional[Dict]:
    """Call Gemma 4 via Ollama's native /api/chat endpoint."""
    if not _check_ollama_available(ollama_url, model=model):
        logger.warning("Ollama server not reachable or model not available")
        return None

    api_url = f"{ollama_url}/api/chat"
    user_msg = (
        f"Analyze this biological context and estimate optimal SynVoy parameters.\n\n"
        f"Context:\n{json.dumps(context, indent=2)}"
    )

    payload = json.dumps(
        {
            "model": model,
            "messages": [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": user_msg},
            ],
            "stream": False,
            "think": False,  # Disable reasoning trace — critical for CPU speed
            "options": {
                "temperature": 0.1,
                "num_predict": 1024,
            },
        }
    ).encode()

    headers = {"Content-Type": "application/json"}

    req = urllib.request.Request(api_url, data=payload, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            result = json.loads(resp.read().decode())
            content = result.get("message", {}).get("content", "")
            eval_count = result.get("eval_count", 0)
            total_s = round(result.get("total_duration", 0) / 1e9, 1)
            logger.info(f"Ollama response: {eval_count} tokens in {total_s}s")
            return _parse_llm_json(content)
    except (TimeoutError, socket.timeout) as exc:
        logger.error(f"Ollama call timed out after {timeout}s: {exc}")
        return None
    except Exception as exc:
        logger.error(f"Ollama call failed: {exc}")
        return None


# ---------------------------------------------------------------------------
# Backend: Google Cloud Gemini API
# ---------------------------------------------------------------------------


def _call_google_cloud(
    context: Dict,
    api_key: str,
    timeout: int = 30,
) -> Optional[Dict]:
    """Call Google Gemini API for parameter estimation."""
    if not api_key:
        return None

    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-lite:generateContent?key={api_key}"

    user_msg = (
        f"Analyze this biological context and estimate optimal SynVoy parameters.\n\n"
        f"Context:\n{json.dumps(context, indent=2)}"
    )

    payload = json.dumps(
        {
            "system_instruction": {"parts": [{"text": SYSTEM_PROMPT}]},
            "contents": [{"parts": [{"text": user_msg}]}],
            "generationConfig": {
                "temperature": 0.1,
                "maxOutputTokens": 1024,
                "responseMimeType": "application/json",
            },
        }
    ).encode()

    headers = {"Content-Type": "application/json"}

    # Retry on 429 (rate limit) with backoff — free tier has low rpm limits
    max_retries = 3
    for attempt in range(max_retries):
        req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                result = json.loads(resp.read().decode())
                candidates = result.get("candidates", [])
                if candidates:
                    content = candidates[0].get("content", {}).get("parts", [{}])[0].get("text", "")
                    return _parse_llm_json(content)
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 503) and attempt < max_retries - 1:
                wait = 15 * (attempt + 1)
                logger.warning(f"Gemini HTTP {exc.code}, retrying in {wait}s ({attempt+1}/{max_retries})...")
                import time
                time.sleep(wait)
                continue
            logger.error(f"Google Cloud API call failed: {exc}")
        except Exception as exc:
            logger.error(f"Google Cloud API call failed: {exc}")
            break
    return None


# ---------------------------------------------------------------------------
# LLM output parsing and validation
# ---------------------------------------------------------------------------


def _parse_llm_json(raw_text: str) -> Optional[Dict]:
    """Extract JSON from LLM output, handling markdown fences etc."""
    if not raw_text:
        return None

    text = raw_text.strip()

    # Strip markdown code fences
    text = re.sub(r"^```(?:json)?\s*\n?", "", text)
    text = re.sub(r"\n?```\s*$", "", text)
    text = text.strip()

    # Try parsing directly
    try:
        result = json.loads(text)
        if isinstance(result, dict):
            return result
    except json.JSONDecodeError:
        pass

    # Try extracting JSON object from surrounding text
    match = re.search(r"\{[^{}]*\}", text, re.DOTALL)
    if match:
        try:
            result = json.loads(match.group())
            if isinstance(result, dict):
                return result
        except json.JSONDecodeError:
            pass

    # Try nested braces
    match = re.search(r"\{(?:[^{}]|\{[^{}]*\})*\}", text, re.DOTALL)
    if match:
        try:
            result = json.loads(match.group())
            if isinstance(result, dict):
                return result
        except json.JSONDecodeError:
            pass

    logger.error(f"Could not parse LLM output as JSON: {text[:200]}")
    return None


def validate_and_clamp(raw_params: Dict) -> Tuple[Dict, List[str]]:
    """
    Validate LLM-estimated parameters against allowed ranges.
    Returns (clamped_params, list_of_warnings).
    """
    validated = {}
    warnings = []

    for key, value in raw_params.items():
        if key not in ESTIMABLE_PARAMS:
            warnings.append(f"Ignoring unknown parameter: {key}")
            continue

        spec = ESTIMABLE_PARAMS[key]
        expected_type = spec["type"]

        # Type coercion
        try:
            if expected_type == bool:
                if isinstance(value, bool):
                    coerced = value
                elif isinstance(value, str):
                    coerced = value.strip().lower() in ("true", "1", "yes")
                else:
                    coerced = bool(value)
            elif expected_type == int:
                coerced = int(float(value))
            elif expected_type == float:
                coerced = float(value)
            else:
                coerced = value
        except (ValueError, TypeError):
            warnings.append(f"Cannot convert {key}={value} to {expected_type.__name__}, skipping")
            continue

        # Range clamping for numeric types
        if expected_type in (int, float) and "min" in spec and "max" in spec:
            original = coerced
            coerced = max(spec["min"], min(spec["max"], coerced))
            if coerced != original:
                warnings.append(
                    f"Clamped {key}: {original} → {coerced} "
                    f"(range [{spec['min']}, {spec['max']}])"
                )

        validated[key] = coerced

    return validated, warnings


def detect_breaking_combos(params: Dict) -> Tuple[Dict, List[str]]:
    """
    Detect parameter combinations that are likely to break the pipeline
    or produce garbage results.  Returns (fixed_params, issues) — the
    input dict is NOT mutated.
    """
    params = dict(params)  # shallow copy to avoid side effects
    issues = []

    # Padding consistency
    pmin = params.get("padding_min", 50000)
    pmax = params.get("padding_max", 200000)
    if pmax < pmin:
        issues.append(
            f"BREAKING: padding_max ({pmax}) < padding_min ({pmin}). "
            f"Auto-fixing: setting padding_max = padding_min + 50000"
        )
        params["padding_max"] = pmin + 50000

    # Overly strict synteny + few flanking genes
    min_score = params.get("min_synteny_score", 0.6)
    n_flank = params.get("n_flanking_genes", 10)
    if min_score > 0.8 and n_flank <= 4:
        issues.append(
            f"WARNING: min_synteny_score={min_score} with n_flanking_genes={n_flank} "
            f"means ≥{int(min_score * n_flank)} of {n_flank} genes must match. "
            f"This is extremely strict and may produce zero results."
        )

    # Relaxed e-value with strict identity
    evalue = params.get("search_evalue", 0.01)
    identity = params.get("min_hit_identity", 10)
    if evalue > 1.0 and identity > 30:
        issues.append(
            f"WARNING: search_evalue={evalue} is very relaxed but "
            f"min_hit_identity={identity} is strict. These conflict — "
            f"the relaxed e-value finds distant hits that the identity filter discards."
        )

    # PLM/structural enabled without Smith-Waterman
    if params.get("enable_plm_search") and not params.get("enable_smith_waterman", True):
        issues.append(
            "WARNING: PLM search enabled but Smith-Waterman disabled. "
            "PLM-discovered ORFs benefit from SW validation."
        )

    return params, issues


# ---------------------------------------------------------------------------
# Heuristic fallback (deterministic, no LLM needed)
# ---------------------------------------------------------------------------


def heuristic_estimate(context: Dict) -> Dict[str, Any]:
    """
    Rule-based parameter estimation as fallback when no LLM is available.
    Encodes the same biological reasoning as the system prompt.
    """
    params: Dict[str, Any] = {}
    query = context.get("query", {})
    home = context.get("home_species", {})
    targets = context.get("target_context", {})

    kingdom = home.get("kingdom", "Animalia")
    genome_mb = home.get("genome_size_mb", 0)
    query_length = query.get("length_aa", 300)
    max_dist_mya = targets.get("max_evolutionary_distance_mya", 0)
    gene_family = query.get("gene_family_size_estimate", "unknown")

    # ---- Fallback for unknown kingdom ----
    # When NCBI lookup failed, guess kingdom from genome size or default to Animalia
    if kingdom == "Unknown":
        if genome_mb > 0:
            if genome_mb < 20:
                kingdom = "Bacteria"
            elif genome_mb < 100:
                kingdom = "Fungi"
            elif genome_mb < 1000:
                kingdom = "Animalia"
            else:
                kingdom = "Plantae"
        else:
            kingdom = "Animalia"  # safe default — moderate introns/padding
        logger.info(f"Unknown kingdom, falling back to '{kingdom}' (genome_size={genome_mb}Mb)")

    # ---- Kingdom-based genome architecture ----
    if kingdom in ("Bacteria", "Archaea"):
        params["max_intron"] = 0 if kingdom == "Bacteria" else 200
        params["cluster_distance"] = 20000
        params["region_padding"] = 20000
        params["padding_min"] = 5000
        params["padding_max"] = 30000
        params["n_flanking_genes"] = 5
        params["min_flanking_size"] = 200
    elif kingdom == "Fungi":
        params["max_intron"] = 500
        params["cluster_distance"] = 40000
        params["region_padding"] = 30000
        params["padding_max"] = 60000
        params["n_flanking_genes"] = 8
    elif kingdom == "Plantae":
        params["max_intron"] = 80000
        params["cluster_distance"] = 400000
        params["region_padding"] = 300000
        params["padding_min"] = 100000
        params["padding_max"] = 500000
        if genome_mb > 3000:
            params["max_intron"] = 100000
            params["cluster_distance"] = 600000
            params["region_padding"] = 400000
            params["padding_max"] = 700000
        params["min_synteny_score"] = 0.4  # plants are more rearranged
        params["n_flanking_genes"] = 12
    elif kingdom == "Animalia":
        if genome_mb > 2000:
            # Large vertebrate genome
            params["max_intron"] = 50000
            params["cluster_distance"] = 300000
            params["region_padding"] = 200000
            params["padding_max"] = 350000
        elif genome_mb < 300:
            # Compact genome (insects, nematodes)
            params["max_intron"] = 15000
            params["cluster_distance"] = 100000

    # ---- Query-size adaptive ----
    if query_length < 80:
        # Small peptide (melittin, defensin)
        params["sw_min_score"] = 10
        params["min_hit_length"] = 8
        params["min_hit_identity"] = 8
    elif query_length < 150:
        params["sw_min_score"] = 15
        params["min_hit_length"] = 8
    elif query_length > 1000:
        # Large multi-domain protein
        params["sw_min_score"] = 30
        params["max_blocks_per_genome"] = 120

    # ---- Gene family adaptive ----
    if gene_family == "large" or query.get("domain_families"):
        params["max_flanking_goi_similarity"] = 22
        params["expand_goi_similar"] = True
        if not params.get("max_flanking_distance"):
            params["max_flanking_distance"] = 400000
    elif gene_family == "medium":
        params["max_flanking_goi_similarity"] = 28

    # ---- Evolutionary distance adaptive ----
    if max_dist_mya > 700:
        params["mmseqs_sensitivity"] = 11
        params["min_hit_identity"] = 5
        params["min_gene_identity"] = 10
        params["search_evalue"] = 1.0
        params["min_synteny_score"] = min(params.get("min_synteny_score", 0.6), 0.3)
        params["enable_plm_search"] = True
        params["enable_structural_search"] = True
    elif max_dist_mya > 400:
        params["mmseqs_sensitivity"] = 10.5
        params["min_hit_identity"] = 8
        params["min_gene_identity"] = 15
        params["search_evalue"] = 0.1
        params["min_synteny_score"] = min(params.get("min_synteny_score", 0.6), 0.4)
        params["enable_plm_search"] = True
    elif max_dist_mya > 200:
        params["mmseqs_sensitivity"] = 10
        params["min_gene_identity"] = 20
        params["search_evalue"] = 0.05

    # ---- Remove entries that match defaults ----
    cleaned = {}
    for key, value in params.items():
        if key in ESTIMABLE_PARAMS:
            if value != ESTIMABLE_PARAMS[key]["default"]:
                cleaned[key] = value

    return cleaned


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------


def estimate_params(
    context: Dict,
    model: str = "auto",
    ollama_url: str = "http://localhost:11434",
    ollama_timeout: int = DEFAULT_OLLAMA_TIMEOUT,
    google_api_key: str = "",
) -> Dict[str, Any]:
    """
    Estimate pipeline parameters using the best available backend.
    Priority: Ollama → Google Cloud → Heuristic.
    """
    raw_params: Optional[Dict] = None
    backend_used = "none"

    # --- Try Ollama first ---
    if model == "auto":
        try:
            model = _auto_select_model()
        except Exception:
            model = "gemma4:e4b"

    if _check_ollama_available(ollama_url):
        logger.info(f"Attempting LLM estimation via Ollama ({model}, timeout={ollama_timeout}s)...")
        raw_params = _call_ollama(context, model, ollama_url, timeout=ollama_timeout)
        if raw_params is not None:
            backend_used = f"ollama:{model}"
            logger.info(f"LLM estimation successful ({backend_used})")

    # --- Try Google Cloud as fallback ---
    if raw_params is None and google_api_key:
        logger.info("Attempting LLM estimation via Google Cloud Gemini API...")
        raw_params = _call_google_cloud(context, google_api_key)
        if raw_params is not None:
            backend_used = "google_cloud:gemini"
            logger.info(f"LLM estimation successful ({backend_used})")

    # --- Heuristic fallback ---
    if raw_params is None:
        logger.info(
            "No LLM backend available. Using heuristic parameter estimation. "
            "Install Ollama (https://ollama.com) with 'ollama pull gemma4:e4b' "
            "for LLM-powered estimation."
        )
        raw_params = heuristic_estimate(context)
        backend_used = "heuristic"

    # Validate and clamp
    validated_params, warnings = validate_and_clamp(raw_params)
    for w in warnings:
        logger.warning(f"Validation: {w}")

    # Detect breaking combinations
    validated_params, issues = detect_breaking_combos(validated_params)
    for issue in issues:
        if issue.startswith("BREAKING"):
            logger.error(issue)
        else:
            logger.warning(issue)

    # Build result with metadata
    result = {
        "backend": backend_used,
        "parameters": validated_params,
        "warnings": warnings,
        "issues": issues,
        "context_summary": {
            "kingdom": context.get("home_species", {}).get("kingdom", "Unknown"),
            "genome_size_mb": context.get("home_species", {}).get("genome_size_mb", 0),
            "query_length_aa": context.get("query", {}).get("length_aa", 0),
            "max_distance_mya": context.get("target_context", {}).get(
                "max_evolutionary_distance_mya", 0
            ),
        },
    }

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="LLM-powered parameter estimation for SynVoy"
    )
    parser.add_argument(
        "--context", required=True, help="Path to context.json from build_llm_context.py"
    )
    parser.add_argument(
        "--model",
        default="auto",
        help="Ollama model name or 'auto' for resource-based selection (default: auto)",
    )
    parser.add_argument(
        "--ollama_url",
        default="http://localhost:11434",
        help="Ollama server URL (default: http://localhost:11434)",
    )
    parser.add_argument(
        "--ollama_timeout",
        type=int,
        default=DEFAULT_OLLAMA_TIMEOUT,
        help=f"Ollama request timeout in seconds (default: {DEFAULT_OLLAMA_TIMEOUT})",
    )
    parser.add_argument(
        "--google_api_key",
        default="",
        help="Google Cloud Gemini API key (optional fallback)",
    )
    parser.add_argument("--output", required=True, help="Output JSON path")

    args = parser.parse_args()

    # Load context
    with open(args.context) as fh:
        context = json.load(fh)

    # Get Google API key from env if not in CLI
    google_key = args.google_api_key or os.environ.get("GOOGLE_API_KEY", "")

    # Estimate
    result = estimate_params(
        context=context,
        model=args.model,
        ollama_url=args.ollama_url,
        ollama_timeout=max(args.ollama_timeout, 1),
        google_api_key=google_key,
    )

    # Write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    # Log summary
    params = result["parameters"]
    logger.info(f"Parameter estimation complete (backend: {result['backend']})")
    logger.info(f"Overriding {len(params)} parameter(s):")
    for key, value in sorted(params.items()):
        default = ESTIMABLE_PARAMS.get(key, {}).get("default", "?")
        logger.info(f"  {key}: {default} → {value}")

    if result["warnings"]:
        logger.info(f"Warnings: {len(result['warnings'])}")
    if result["issues"]:
        logger.info(f"Issues: {len(result['issues'])}")

    # Also print the flat params to stdout for easy consumption
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
