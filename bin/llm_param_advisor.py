#!/usr/bin/env python3
"""
llm_param_advisor.py — LLM-powered parameter estimation for SynVoy.

Uses a cloud LLM API (Google Gemini or OpenAI-compatible) to analyze biological
context and estimate optimal pipeline parameters.  Falls back to deterministic
heuristics when no API key is available.

All estimated parameters are validated against allowed ranges and checked
for breaking combinations before being emitted.
"""

import argparse
import json
import logging
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sequence_utils import setup_logging  # noqa: E402

logger = setup_logging(name=__name__)

# ---------------------------------------------------------------------------
# Parameter specifications: allowed ranges, types, and defaults
# ---------------------------------------------------------------------------

# Only Tier 1 + Tier 2 parameters are LLM-estimable.
# Tier 3 (classification thresholds, scoring weights, viz) are never touched.
# The distances of the search (block gap, window margins, intron limit) are not
# estimated here any more: since 2026-10-08 the search derives them from the home
# neighbourhood and from each target genome (bin/synvoy_scale.py), which an estimate by
# kingdom and genome size cannot improve on.
# Also not estimated since 2026-10-09: the hit floors (min_hit_identity, min_hit_length,
# sw_min_score, sw_min_identity) -- lowering them changes nothing, because the E-value
# filters of the search are stricter (a 20-aa query in a 2 kb window needs a raw
# Smith-Waterman score of 31 for E <= 1; the floor is 20) -- and the two family options
# (expand_goi_similar*, max_flanking_distance): the first lowered recall where it was
# measured, the second is a distance in bp.
ESTIMABLE_PARAMS: Dict[str, Dict[str, Any]] = {
    "n_flanking_genes": {"type": int, "min": 2, "max": 30, "default": 10},
    "search_evalue": {"type": float, "min": 1e-10, "max": 100.0, "default": 0.01},
    "mmseqs_sensitivity": {"type": float, "min": 1.0, "max": 12.0, "default": 9.5},
    "enable_smith_waterman": {"type": bool, "default": True},
    "enable_plm_search": {"type": bool, "default": False},
    "enable_structural_search": {"type": bool, "default": False},
    "max_blocks_per_genome": {"type": int, "min": 5, "max": 500, "default": 80},
    "min_block_genes": {"type": int, "min": 1, "max": 10, "default": 2},
    "min_flanking_size": {"type": int, "min": 100, "max": 5000, "default": 500},
    "exon_level_search": {"type": bool, "default": True},
}


# ---------------------------------------------------------------------------
# System prompt
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
(sequence homology). Hits that lie close together form a candidate synteny block; \
how close is derived by the pipeline from the gene spacing of the home neighbourhood \
and of each target genome.
4. A candidate block needs at least min_block_genes distinct flanking genes; \
at most max_blocks_per_genome blocks are searched per genome.
5. Within each candidate block, the pipeline runs gene prediction \
(Augustus for eukaryotes, Prodigal for prokaryotes) to find ORFs, then \
searches for the GOI itself via tblastn, Smith-Waterman, and optionally \
PLM embeddings (ProtT5) and structural comparison (ESMFold + Foldseek).
6. All distances of the search (block gap, window margins, intron limits) are derived \
by the pipeline itself. Do NOT output distance parameters.

KEY TRADE-OFFS:
- Strict params (low E-value, low sensitivity) → high precision, low recall. \
Good for close species. Misses distant orthologs.
- Relaxed params (high E-value, high sensitivity) → high recall, \
lower precision. Needed for distant species. Slower, more false positives.
- Small peptides (<80 aa): alignments are judged by their E-value inside the search \
window, so no score, identity or length threshold has to be lowered; nothing to set.
- Gene families (paralogs): the pipeline keeps relatives of the GOI out of the \
anchors by an alignment E-value; nothing to set.

CONTEXT JSON you will receive has these fields:
- query.length_aa: protein length — very small (<80), small (80-150), normal (150-800), large (>800)
- query.estimated_exon_count: 1 for single-exon, >1 for multi-exon genes
- query.has_signal_peptide / is_secreted: secreted peptides are often small/fast-evolving
- query.gene_family_size_estimate: "small", "medium", or "large"
- query.domain_families: known protein domains (empty = novel/unknown)
- home_species.kingdom: Bacteria, Archaea, Fungi, Plantae, Animalia, or Unknown
- home_species.genome_size_mb: 0 means lookup failed — infer from kingdom
- home_species.gene_count: total annotated genes
- home_species.avg_intron_length_bp: average intron length (for information; intron \
limits are derived by the pipeline)
- home_species.avg_gene_density_per_mb: genes per Mb — high density (>100) means fewer \
flanking genes needed; low density (<20) means more needed
- target_context.max_evolutionary_distance_mya: estimated divergence to most distant target
- target_context.kingdoms_represented: kingdoms in the target set
- target_context.target_count: number of target genomes

PARAMETER REFERENCE (only override parameters that should differ from defaults):

Genome Architecture:
- n_flanking_genes (int, default=10): Synteny anchors per side of GOI.
  Gene-dense genomes (>100 genes/Mb)=5-8, Normal=10, Gene-sparse or rearranged=12-20
- min_flanking_size (int, default=500): Min flanking gene size in bp; used only with \
exon_level_search=false. Bacteria=200-300, Eukaryotes=500 (default).
- exon_level_search (bool, default=true): Search at exon level for divergent gene detection. \
Keep true for cross-order+ searches. May disable for within-genus searches to save time.

Search Sensitivity:
- search_evalue (float, default=0.01): MMseqs2 e-value cutoff.
  Close=0.001, Moderate=0.01-0.1, Distant=1.0-10
- mmseqs_sensitivity (float, default=9.5): MMseqs2 sensitivity (1-12, higher=slower+more sensitive).
  Close=7-8, Default=9.5, Distant=10-11, Very distant (>500 Mya)=11-12

Smith-Waterman:
- enable_smith_waterman (bool, default=true): Use rigorous SW alignment. Keep true unless \
searching only very close species where tblastn alone suffices.

Advanced Search:
- enable_plm_search (bool, default=false): ProtT5 protein language model embedding search. \
Finds remote homologs missed by sequence alignment. Enable for distance >400 Mya.
- enable_structural_search (bool, default=false): ESMFold structure prediction + Foldseek 3Di search. \
Finds structural homologs with no sequence similarity. Enable for distance >600 Mya. Requires GPU.

Block Filtering:
- max_blocks_per_genome (int, default=80): Safety cap on candidate blocks per target genome. \
Increase to 120-200 for large gene families; decrease to 20-30 for quick runs.
- min_block_genes (int, default=2): Min flanking gene hits in a block to keep it. \
Keep at 2 for most cases. Set to 1 only for extremely distant searches with few anchors expected.

CRITICAL RULES:
1. Never output a distance in base pairs (cluster_distance, region_padding, max_intron, \
...): the pipeline derives them.
2. Bacteria/Archaea: n_flanking_genes≤8, min_flanking_size≤300
3. Cross-phylum searches (>500 Mya): mmseqs_sensitivity≥10, enable_plm_search=true
4. When genome_size_mb=0 (lookup failed), infer genome architecture from kingdom: \
Bacteria ~5Mb, Fungi ~30Mb, Insects ~200-400Mb, Vertebrates ~1000-3000Mb, Plants ~500-5000Mb
5. Only output parameters you are confident about. If uncertain, omit the parameter \
and let the pipeline default handle it.

REASONING: Think step-by-step before outputting JSON:
1. What kingdom? → sets the number of flanking genes (and their minimum size for prokaryotes)
2. How far apart are the species? → sets search sensitivity and advanced methods
3. Any special genome features (gene density)? → fine-tune the number of flanking genes
Then output ONLY the JSON override object.

OUTPUT FORMAT: A single JSON object with parameter overrides. Include ONLY \
parameters that differ from defaults. No explanations outside the JSON.

EXAMPLES:

Honeybee melittin (70aa secreted peptide) searching in bumblebees (~80 Mya, same family):
{}

Human p53 (393aa tumor suppressor) searching across vertebrates (~450 Mya):
{"mmseqs_sensitivity": 10, "search_evalue": 0.1, "enable_plm_search": true}

Arabidopsis defensin (80aa) searching across angiosperms (~150 Mya):
{"n_flanking_genes": 12}

E. coli beta-lactamase (286aa) in Enterobacteriaceae (~300 Mya):
{"n_flanking_genes": 5, "min_flanking_size": 200, "mmseqs_sensitivity": 10, "search_evalue": 0.1}
"""


# ---------------------------------------------------------------------------
# Backend: Google Gemini API
# ---------------------------------------------------------------------------


def _call_google(
    context: Dict,
    api_key: str,
    model: str = "gemini-2.5-flash-lite",
    timeout: int = 30,
) -> Optional[Dict]:
    """Call Google Gemini API for parameter estimation."""
    if not api_key:
        return None

    url = (
        f"https://generativelanguage.googleapis.com/v1beta/models/"
        f"{model}:generateContent?key={api_key}"
    )

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

    max_retries = 3
    for attempt in range(max_retries):
        req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                result = json.loads(resp.read().decode())
                candidates = result.get("candidates", [])
                if candidates:
                    content = (
                        candidates[0]
                        .get("content", {})
                        .get("parts", [{}])[0]
                        .get("text", "")
                    )
                    return _parse_llm_json(content)
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 503) and attempt < max_retries - 1:
                import time
                wait = 15 * (attempt + 1)
                logger.warning(
                    f"Gemini HTTP {exc.code}, retrying in {wait}s "
                    f"({attempt + 1}/{max_retries})..."
                )
                time.sleep(wait)
                continue
            logger.error(f"Google Gemini API call failed: HTTP {exc.code}")
            body = exc.read().decode(errors="replace")
            logger.debug(f"Response body: {body[:300]}")
        except Exception as exc:
            logger.error(f"Google Gemini API call failed: {exc}")
            break
    return None


# ---------------------------------------------------------------------------
# Backend: OpenAI-compatible API (OpenAI, Together, Groq, LM Studio, etc.)
# ---------------------------------------------------------------------------


def _call_openai(
    context: Dict,
    api_key: str,
    model: str = "gpt-4o-mini",
    api_base_url: str = "https://api.openai.com",
    timeout: int = 30,
) -> Optional[Dict]:
    """Call any OpenAI-compatible chat completions API."""
    if not api_key:
        return None

    url = f"{api_base_url.rstrip('/')}/v1/chat/completions"
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
            "temperature": 0.1,
            "max_tokens": 1024,
            "response_format": {"type": "json_object"},
        }
    ).encode()

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
    }

    max_retries = 3
    for attempt in range(max_retries):
        req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                result = json.loads(resp.read().decode())
                content = (
                    result.get("choices", [{}])[0]
                    .get("message", {})
                    .get("content", "")
                )
                return _parse_llm_json(content)
        except urllib.error.HTTPError as exc:
            if exc.code == 429 and attempt < max_retries - 1:
                import time
                wait = 15 * (attempt + 1)
                logger.warning(
                    f"OpenAI API rate-limited, retrying in {wait}s "
                    f"({attempt + 1}/{max_retries})..."
                )
                time.sleep(wait)
                continue
            logger.error(f"OpenAI-compatible API call failed: HTTP {exc.code}")
            body = exc.read().decode(errors="replace")
            logger.debug(f"Response body: {body[:300]}")
        except Exception as exc:
            logger.error(f"OpenAI-compatible API call failed: {exc}")
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
            warnings.append(
                f"Cannot convert {key}={value} to {expected_type.__name__}, skipping"
            )
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
    home = context.get("home_species", {})
    targets = context.get("target_context", {})

    kingdom = home.get("kingdom", "Animalia")
    genome_mb = home.get("genome_size_mb", 0)
    max_dist_mya = targets.get("max_evolutionary_distance_mya", 0)

    # Fallback for unknown kingdom
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
            kingdom = "Animalia"
        logger.info(f"Unknown kingdom, falling back to '{kingdom}' (genome_size={genome_mb}Mb)")

    # Kingdom-based choices. Distances (block gap, margins, intron limit) are not set
    # here: the search derives them from the home neighbourhood and each target genome.
    if kingdom in ("Bacteria", "Archaea"):
        params["n_flanking_genes"] = 5
        params["min_flanking_size"] = 200
    elif kingdom == "Fungi":
        params["n_flanking_genes"] = 8
    elif kingdom == "Plantae":
        params["n_flanking_genes"] = 12

    # Nothing is set by query length or family size any more (2026-10-09): the score and
    # length floors do not bind below their defaults (the E-value filters are stricter),
    # and relatives of the GOI are kept out of the anchors by an alignment E-value.

    # Evolutionary distance adaptive
    if max_dist_mya > 700:
        params["mmseqs_sensitivity"] = 11
        params["search_evalue"] = 1.0
        params["enable_plm_search"] = True
        params["enable_structural_search"] = True
    elif max_dist_mya > 400:
        params["mmseqs_sensitivity"] = 10.5
        params["search_evalue"] = 0.1
        params["enable_plm_search"] = True
    elif max_dist_mya > 200:
        params["mmseqs_sensitivity"] = 10
        params["search_evalue"] = 0.05

    # Remove entries that match defaults
    cleaned = {}
    for key, value in params.items():
        if key in ESTIMABLE_PARAMS:
            if value != ESTIMABLE_PARAMS[key]["default"]:
                cleaned[key] = value

    return cleaned


# ---------------------------------------------------------------------------
# Default models per provider
# ---------------------------------------------------------------------------

DEFAULT_MODELS = {
    "google": "gemini-2.5-flash-lite",
    "openai": "gpt-4o-mini",
}


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------


def estimate_params(
    context: Dict,
    provider: str = "google",
    api_key: str = "",
    api_base_url: str = "",
    model: str = "",
) -> Dict[str, Any]:
    """
    Estimate pipeline parameters using the best available backend.
    Priority: API key → heuristic.

    provider: 'google' or 'openai' (or any OpenAI-compatible endpoint via api_base_url)
    api_key: provider API key
    api_base_url: custom base URL for OpenAI-compatible providers (e.g. Together, Groq)
    model: model name override; if empty, uses provider default
    """
    raw_params: Optional[Dict] = None
    backend_used = "none"

    resolved_model = model or DEFAULT_MODELS.get(provider, "")

    if api_key:
        if provider == "google":
            logger.info(f"Attempting LLM estimation via Google Gemini ({resolved_model})...")
            raw_params = _call_google(context, api_key, model=resolved_model)
            if raw_params is not None:
                backend_used = f"google:{resolved_model}"
        else:
            # openai or any OpenAI-compatible provider
            base = api_base_url or "https://api.openai.com"
            provider_label = api_base_url if api_base_url else "openai"
            logger.info(
                f"Attempting LLM estimation via OpenAI-compatible API "
                f"({provider_label}, model={resolved_model})..."
            )
            raw_params = _call_openai(
                context, api_key, model=resolved_model, api_base_url=base
            )
            if raw_params is not None:
                backend_used = f"openai:{resolved_model}"

        if raw_params is not None:
            logger.info(f"LLM estimation successful (backend: {backend_used})")
        else:
            logger.warning(
                "LLM API call failed or returned unparseable output. "
                "Falling back to heuristic estimation."
            )
    else:
        logger.info(
            "No LLM API key provided. Using heuristic parameter estimation.\n"
            "To enable LLM-quality estimation, set --auto_params true and provide\n"
            "an API key via --llm_api_key or the LLM_API_KEY environment variable.\n"
            "Supported providers: Google Gemini (--llm_provider google) and\n"
            "OpenAI / compatible APIs (--llm_provider openai)."
        )

    # Heuristic fallback
    if raw_params is None:
        raw_params = heuristic_estimate(context)
        backend_used = "heuristic"
        logger.info("Using heuristic parameter estimation.")

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
        description="LLM-powered parameter estimation for SynVoy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Google Gemini (recommended)
  llm_param_advisor.py --context ctx.json --provider google --api_key $GOOGLE_API_KEY --output params.json

  # OpenAI
  llm_param_advisor.py --context ctx.json --provider openai --api_key $OPENAI_API_KEY --output params.json

  # OpenAI-compatible (e.g. Together, Groq, local LM Studio)
  llm_param_advisor.py --context ctx.json --provider openai \\
      --api_key $TOGETHER_API_KEY --api_base_url https://api.together.xyz \\
      --model meta-llama/Llama-3.1-8B-Instruct-Turbo --output params.json

  # No API key: pure heuristic estimation
  llm_param_advisor.py --context ctx.json --output params.json
""",
    )
    parser.add_argument(
        "--context", required=True, help="Path to context.json from build_llm_context.py"
    )
    parser.add_argument(
        "--provider",
        default="google",
        choices=["google", "openai"],
        help="LLM provider: 'google' (Gemini) or 'openai' (OpenAI or compatible). Default: google",
    )
    parser.add_argument(
        "--api_key",
        default="",
        help=(
            "API key for the chosen provider. "
            "Also read from LLM_API_KEY, GOOGLE_API_KEY (for google), "
            "or OPENAI_API_KEY (for openai) env vars."
        ),
    )
    parser.add_argument(
        "--api_base_url",
        default="",
        help=(
            "Custom API base URL for OpenAI-compatible providers "
            "(e.g. https://api.together.xyz, https://api.groq.com/openai). "
            "Ignored for the google provider."
        ),
    )
    parser.add_argument(
        "--model",
        default="",
        help=(
            "Model name override. Defaults: google=gemini-2.5-flash-lite, "
            "openai=gpt-4o-mini. Use this for custom models on compatible endpoints."
        ),
    )
    parser.add_argument("--output", required=True, help="Output JSON path")

    args = parser.parse_args()

    # Load context
    with open(args.context) as fh:
        context = json.load(fh)

    # Resolve API key: CLI arg > LLM_API_KEY > provider-specific env var
    api_key = args.api_key
    if not api_key:
        api_key = os.environ.get("LLM_API_KEY", "")
    if not api_key and args.provider == "google":
        api_key = os.environ.get("GOOGLE_API_KEY", "")
    if not api_key and args.provider == "openai":
        api_key = os.environ.get("OPENAI_API_KEY", "")

    result = estimate_params(
        context=context,
        provider=args.provider,
        api_key=api_key,
        api_base_url=args.api_base_url,
        model=args.model,
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

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
