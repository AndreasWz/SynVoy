#!/usr/bin/env python3
"""
build_llm_context.py — Gather biological context for LLM parameter estimation.

Scrapes NCBI Datasets API for genome statistics (genome size, gene count,
scaffold N50) and combines with resolved query metadata (protein length,
organism, domains) to produce a structured context JSON for the LLM advisor.

Results are cached to ~/.synvoy/species_cache/ to avoid redundant API calls.
"""

import argparse
import hashlib
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CACHE_DIR = Path.home() / ".synvoy" / "species_cache"
CACHE_TTL_DAYS = 30

NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_TAXONOMY_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

# Rough kingdom→intron/gene-density defaults when NCBI lookup fails.
KINGDOM_DEFAULTS = {
    "Animalia": {
        "avg_intron_length_bp": 5000,
        "avg_gene_density_per_mb": 8.0,
        "is_polyploid": False,
    },
    "Plantae": {
        "avg_intron_length_bp": 4000,
        "avg_gene_density_per_mb": 5.0,
        "is_polyploid": False,  # conservative; many are, but default safe
    },
    "Fungi": {
        "avg_intron_length_bp": 250,
        "avg_gene_density_per_mb": 30.0,
        "is_polyploid": False,
    },
    "Bacteria": {
        "avg_intron_length_bp": 0,
        "avg_gene_density_per_mb": 900.0,
        "is_polyploid": False,
    },
    "Archaea": {
        "avg_intron_length_bp": 0,
        "avg_gene_density_per_mb": 800.0,
        "is_polyploid": False,
    },
    "Protista": {
        "avg_intron_length_bp": 300,
        "avg_gene_density_per_mb": 50.0,
        "is_polyploid": False,
    },
}

# Taxonomic rank → approximate divergence time (Mya) from common ancestor.
# Used as fallback when proper phylogenetic distance is unavailable.
RANK_DISTANCE_MYA = {
    "species": 5,
    "genus": 30,
    "family": 80,
    "order": 150,
    "class": 300,
    "phylum": 550,
    "kingdom": 900,
    "superkingdom": 1500,
}

RANKED_LINEAGE_ORDER = (
    "superkingdom",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)

# ---------------------------------------------------------------------------
# Caching helpers
# ---------------------------------------------------------------------------


def _cache_key(species_name: str) -> str:
    """Stable filesystem-safe cache key from species name."""
    normalized = species_name.strip().lower().replace(" ", "_")
    h = hashlib.md5(normalized.encode()).hexdigest()[:8]
    safe = re.sub(r"[^a-z0-9_]", "", normalized)[:60]
    return f"{safe}_{h}"


def _coerce_float(value: Any, default: float = 0.0) -> float:
    """Best-effort float coercion for mixed string/int cache fields."""
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _lineage_is_trivial(lineage: Any) -> bool:
    """Return True for empty or species-only lineage snapshots."""
    if not isinstance(lineage, list) or not lineage:
        return True

    valid_entries = [entry for entry in lineage if isinstance(entry, dict)]
    if not valid_entries:
        return True
    if len(valid_entries) == 1:
        return valid_entries[0].get("rank", "").lower() == "species"
    return False


def _is_degraded_cache_entry(data: Dict[str, Any]) -> bool:
    """
    Detect caches created from fallback-only species lookups.

    We treat a cache as degraded when taxonomy stayed unknown, genome stats are
    all zero-like, and the lineage is empty or only contains the terminal
    species node. Those snapshots are safe to refetch because they do not add
    meaningful biological context for parameter estimation.
    """
    kingdom = str(data.get("kingdom", "Unknown"))
    phylum = str(data.get("phylum", "Unknown"))
    genome_size_mb = _coerce_float(data.get("genome_size_mb"))
    gene_count = _coerce_int(data.get("gene_count"))
    scaffold_n50 = _coerce_int(data.get("scaffold_n50"))
    lineage = data.get("lineage", [])

    unknown_taxonomy = kingdom == "Unknown" and phylum == "Unknown"
    no_genome_stats = genome_size_mb <= 0 and gene_count <= 0 and scaffold_n50 <= 0
    trivial_lineage = _lineage_is_trivial(lineage)
    return unknown_taxonomy and no_genome_stats and trivial_lineage


def _read_cache(species_name: str) -> Optional[Dict]:
    key = _cache_key(species_name)
    cache_file = CACHE_DIR / f"{key}.json"
    if not cache_file.exists():
        return None
    try:
        data = json.loads(cache_file.read_text())
        cached_at = data.get("_cached_at", 0)
        if time.time() - cached_at > CACHE_TTL_DAYS * 86400:
            return None  # expired
        if _is_degraded_cache_entry(data):
            logger.info(f"Ignoring degraded species cache for '{species_name}'")
            return None
        return data
    except Exception:
        return None


def _write_cache(species_name: str, data: Dict) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        key = _cache_key(species_name)
        cache_file = CACHE_DIR / f"{key}.json"
        data["_cached_at"] = time.time()
        cache_file.write_text(json.dumps(data, indent=2))
    except Exception as exc:
        logger.warning(f"Could not write species cache: {exc}")


# ---------------------------------------------------------------------------
# NCBI API helpers
# ---------------------------------------------------------------------------


def _ncbi_api_key_header() -> Dict[str, str]:
    """Return NCBI API key header if set in environment."""
    key = os.environ.get("NCBI_API_KEY", "")
    if key:
        return {"api-key": key}
    return {}


def _fetch_json(url: str, timeout: int = 20) -> Optional[Dict]:
    """Fetch JSON from a URL, return None on failure."""
    headers = {"Accept": "application/json"}
    headers.update(_ncbi_api_key_header())
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError, Exception) as exc:
        logger.warning(f"API request failed ({url}): {exc}")
        return None


def _extract_taxon_name(value: Any) -> str:
    """Return a taxon/scientific name from a few common NCBI response shapes."""
    if isinstance(value, str):
        return value
    if not isinstance(value, dict):
        return ""

    current_name = value.get("currentScientificName")
    if isinstance(current_name, dict) and current_name.get("name"):
        return current_name["name"]

    for key in ("name", "scientific_name", "organism_name", "organismName", "value"):
        raw = value.get(key)
        if isinstance(raw, str) and raw:
            return raw

    return ""


def _normalize_lineage_entry(rank: str, value: Any) -> Optional[Dict[str, str]]:
    """Normalize one taxonomy/classification entry into the local lineage shape."""
    if not rank:
        return None

    name = _extract_taxon_name(value)
    if not name:
        return None

    normalized_rank = rank.lower()
    if normalized_rank == "domain":
        # NCBI commonly uses "domain" in classification where SynVoy expects superkingdom.
        normalized_rank = "superkingdom"

    return {"name": name, "rank": normalized_rank}


def _resolve_lineage_taxids(taxids: List[int]) -> List[Dict[str, str]]:
    """
    Batch-resolve a list of integer taxids into ranked lineage entries.

    The NCBI Datasets v2 taxonomy endpoint returns lineage as integer taxid
    arrays.  This function queries the batch endpoint to get rank+name for
    each, keeping only the major ranks SynVoy uses.
    """
    if not taxids:
        return []

    # Only resolve major-rank taxids — filter unranked ancestors later
    # Send up to ~30 taxids per batch (API supports comma-separated)
    wanted_ranks = set(RANKED_LINEAGE_ORDER)
    wanted_ranks.add("domain")  # NCBI sometimes uses "domain" instead of "superkingdom"
    results: List[Dict[str, str]] = []

    batch_size = 30
    for i in range(0, len(taxids), batch_size):
        chunk = taxids[i : i + batch_size]
        ids_str = ",".join(str(t) for t in chunk)
        url = f"{NCBI_TAXONOMY_BASE}/taxon/{ids_str}"
        data = _fetch_json(url, timeout=15)
        if not data:
            continue
        for tnode in data.get("taxonomy_nodes", []):
            tax = tnode.get("taxonomy", {})
            rank = (tax.get("rank") or "").lower()
            name = tax.get("organism_name") or ""
            if rank and name and rank in wanted_ranks:
                normalized = _normalize_lineage_entry(rank, name)
                if normalized:
                    results.append(normalized)

    return results


def _extract_ranked_lineage(node: Dict[str, Any]) -> List[Dict[str, str]]:
    """
    Extract lineage entries from an NCBI taxonomy node.

    The Datasets v2 API returns ``lineage`` as an array of integer taxids.
    When that happens we batch-resolve them via a second API call.  If the
    response already contains a ``classification`` dict with named ranks
    (older/alternative response shape), we use that directly.
    """
    lineage: List[Dict[str, str]] = []
    seen_ranks: set = set()

    # Path A: classification dict (named ranks) — may be present in some responses
    classification = node.get("classification", {})
    if isinstance(classification, dict):
        for raw_rank in RANKED_LINEAGE_ORDER:
            normalized = _normalize_lineage_entry(raw_rank, classification.get(raw_rank))
            if normalized and normalized["rank"] not in seen_ranks:
                lineage.append(normalized)
                seen_ranks.add(normalized["rank"])

    # Path B: lineage as list — may be integer taxids or dicts
    raw_lineage = node.get("lineage", [])
    if isinstance(raw_lineage, list) and raw_lineage:
        if all(isinstance(e, int) for e in raw_lineage):
            # Integer taxid array — batch-resolve to get rank+name
            resolved = _resolve_lineage_taxids(raw_lineage)
            for entry in resolved:
                if entry["rank"] not in seen_ranks:
                    lineage.append(entry)
                    seen_ranks.add(entry["rank"])
        else:
            # Dict-shaped entries (older API or pre-processed)
            for entry in raw_lineage:
                if not isinstance(entry, dict):
                    continue
                rank = entry.get("rank") or entry.get("taxon_rank") or entry.get("rank_type")
                normalized = _normalize_lineage_entry(rank, entry)
                if normalized and normalized["rank"] not in seen_ranks:
                    lineage.append(normalized)
                    seen_ranks.add(normalized["rank"])

    # Always include the current node itself
    current_rank = node.get("rank", "")
    current_name = _extract_taxon_name(node.get("currentScientificName")) or _extract_taxon_name(node)
    normalized_current = _normalize_lineage_entry(current_rank, current_name)
    if normalized_current and normalized_current["rank"] not in seen_ranks:
        lineage.append(normalized_current)

    return lineage


def _coerce_int(value: Any, default: int = 0) -> int:
    """Best-effort integer coercion for mixed string/int API fields."""
    try:
        if value in (None, ""):
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Species / genome info from NCBI Datasets v2
# ---------------------------------------------------------------------------


def _classify_kingdom(lineage: List[Dict]) -> str:
    """
    Determine kingdom from NCBI taxonomy lineage.
    lineage: list of dicts with 'name' and 'rank' keys.
    """
    rank_names = {entry.get("rank", "").lower(): entry.get("name", "") for entry in lineage}

    superkingdom = rank_names.get("superkingdom", "")
    kingdom = rank_names.get("kingdom", "")
    # Normalize clade-level entries
    clade_names = [
        entry.get("name", "")
        for entry in lineage
        if entry.get("rank", "").lower() in ("clade", "no rank", "")
    ]

    if superkingdom == "Bacteria":
        return "Bacteria"
    if superkingdom == "Archaea":
        return "Archaea"
    if superkingdom == "Viruses":
        return "Viruses"

    if kingdom == "Metazoa" or "Metazoa" in clade_names:
        return "Animalia"
    if kingdom in ("Viridiplantae", "Streptophyta") or "Viridiplantae" in clade_names:
        return "Plantae"
    if kingdom == "Fungi" or "Fungi" in clade_names:
        return "Fungi"

    # Broad eukaryote fallback
    if superkingdom == "Eukaryota":
        return "Protista"

    return "Unknown"


def _estimate_avg_intron(kingdom: str, genome_size_mb: float) -> int:
    """Rough intron length estimate from kingdom + genome size."""
    defaults = KINGDOM_DEFAULTS.get(kingdom, KINGDOM_DEFAULTS["Animalia"])
    base = defaults["avg_intron_length_bp"]

    if kingdom == "Animalia":
        # Larger animal genomes tend to have longer introns
        if genome_size_mb > 2000:
            return max(base, 6000)
        if genome_size_mb > 500:
            return max(base, 3000)
        if genome_size_mb < 200:
            return min(base, 1500)
    elif kingdom == "Plantae":
        # Plant introns scale with genome size but plateau
        if genome_size_mb > 5000:
            return 8000
        if genome_size_mb > 1000:
            return 5000
        return 3000
    return base


def fetch_species_info(species_name: str) -> Dict[str, Any]:
    """
    Fetch genome statistics for a species from NCBI Datasets API.
    Returns a dict with genome_size_mb, gene_count, kingdom, etc.
    Falls back to heuristics if API is unavailable.
    """
    # Check cache first
    cached = _read_cache(species_name)
    if cached:
        logger.info(f"Species info for '{species_name}' loaded from cache")
        cached.pop("_cached_at", None)
        return cached

    result = {
        "name": species_name,
        "kingdom": "Unknown",
        "phylum": "Unknown",
        "genome_size_mb": 0,
        "gene_count": 0,
        "scaffold_n50": 0,
        "avg_intron_length_bp": 5000,
        "avg_gene_density_per_mb": 8.0,
        "is_polyploid": False,
        "lineage": [],
    }

    # Step 1: Resolve taxonomy to get taxid + lineage
    tax_url = f"{NCBI_TAXONOMY_BASE}/taxon/{urllib.parse.quote(species_name)}"
    tax_data = _fetch_json(tax_url)

    taxid = None
    lineage_items: List[Dict] = []

    if tax_data:
        # Navigate NCBI Datasets v2 taxonomy response
        taxonomy = tax_data.get("taxonomy_nodes", [])
        if taxonomy:
            node = taxonomy[0].get("taxonomy", {})
            taxid = node.get("tax_id") or node.get("taxId")
            lineage_items = _extract_ranked_lineage(node)
            result["lineage"] = list(lineage_items)
            result["phylum"] = next(
                (
                    li.get("name", "Unknown")
                    for li in lineage_items
                    if li.get("rank", "").lower() == "phylum"
                ),
                "Unknown",
            )
    else:
        logger.warning(f"Could not resolve taxonomy for '{species_name}'")

    # Classify kingdom from lineage
    if lineage_items:
        result["kingdom"] = _classify_kingdom(lineage_items)
    elif tax_data is None:
        # Pure heuristic from name patterns
        lower = species_name.lower()
        if any(kw in lower for kw in ("arabidopsis", "oryza", "zea", "solanum", "nicotiana")):
            result["kingdom"] = "Plantae"
        elif any(kw in lower for kw in ("escherichia", "salmonella", "bacillus", "streptococcus")):
            result["kingdom"] = "Bacteria"
        elif any(kw in lower for kw in ("saccharomyces", "aspergillus", "neurospora")):
            result["kingdom"] = "Fungi"

    # Step 2: Fetch genome assembly stats
    if taxid:
        genome_url = (
            f"{NCBI_DATASETS_BASE}/genome/taxon/{taxid}/dataset_report"
            f"?filters.reference_only=true&page_size=1"
        )
        genome_data = _fetch_json(genome_url)
        if genome_data:
            reports = genome_data.get("reports", [])
            if reports:
                asm = reports[0].get("assembly_stats") or reports[0].get("assemblyStats") or {}
                ann = reports[0].get("annotation_info") or reports[0].get("annotationInfo") or {}

                total_length = _coerce_int(
                    asm.get("total_sequence_length") or asm.get("totalSequenceLength")
                )
                result["genome_size_mb"] = round(total_length / 1e6, 1) if total_length else 0
                result["scaffold_n50"] = _coerce_int(
                    asm.get("scaffold_n50") or asm.get("scaffoldN50")
                )
                result["gene_count"] = _coerce_int(
                    ann.get("gene_count_total")
                    or ann.get("geneCountTotal")
                    or ann.get("stats", {}).get("gene_counts", {}).get("total", 0)
                    or ann.get("stats", {}).get("geneCounts", {}).get("total", 0)
                )
        else:
            logger.warning(f"Could not fetch genome stats for taxid {taxid}")

    # Step 3: Derive secondary stats
    kingdom = result["kingdom"]
    genome_mb = result["genome_size_mb"]

    result["avg_intron_length_bp"] = _estimate_avg_intron(kingdom, genome_mb)

    if result["gene_count"] and genome_mb > 0:
        result["avg_gene_density_per_mb"] = round(result["gene_count"] / genome_mb, 2)
    else:
        defaults = KINGDOM_DEFAULTS.get(kingdom, KINGDOM_DEFAULTS["Animalia"])
        result["avg_gene_density_per_mb"] = defaults["avg_gene_density_per_mb"]

    # Cache only when the fetch produced usable biological context.
    if _is_degraded_cache_entry(result):
        logger.info(f"Skipping degraded species cache for '{species_name}'")
    else:
        _write_cache(species_name, dict(result))
    logger.info(
        f"Species info for '{species_name}': kingdom={result['kingdom']}, "
        f"genome={result['genome_size_mb']}Mb, genes={result['gene_count']}, "
        f"avg_intron={result['avg_intron_length_bp']}bp"
    )

    return result


# ---------------------------------------------------------------------------
# Query protein analysis
# ---------------------------------------------------------------------------


def _parse_fasta_length(fasta_path: str) -> int:
    """Return the length of the first protein sequence in a FASTA file."""
    total = 0
    started = False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if started:
                    break  # only first sequence
                started = True
                continue
            total += len(line)
    return total


def _fetch_uniprot_features(uniprot_id: str) -> Dict[str, Any]:
    """Fetch domain/family/signal peptide info from UniProt."""
    result = {
        "has_signal_peptide": False,
        "is_secreted": False,
        "domain_families": [],
        "gene_family_size_estimate": "unknown",
        "protein_name": "",
    }
    if not uniprot_id:
        return result

    # Clean accession
    clean_id = uniprot_id.strip().split(".")[0]
    url = f"{UNIPROT_BASE}/{clean_id}.json"
    data = _fetch_json(url, timeout=15)
    if not data:
        return result

    # Protein name
    prot_desc = data.get("proteinDescription", {})
    rec = prot_desc.get("recommendedName", {})
    result["protein_name"] = rec.get("fullName", {}).get("value", "")

    # Features: signal peptides, domains
    features = data.get("features", [])
    for feat in features:
        ftype = feat.get("type", "")
        if ftype == "Signal":
            result["has_signal_peptide"] = True
            result["is_secreted"] = True
        if ftype in ("Domain", "Region"):
            desc = feat.get("description", "")
            if desc and desc not in result["domain_families"]:
                result["domain_families"].append(desc)

    # Subcellular location
    comments = data.get("comments", [])
    for comment in comments:
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locations = comment.get("subcellularLocations", [])
            for loc in locations:
                loc_val = loc.get("location", {}).get("value", "").lower()
                if "secret" in loc_val or "extracellular" in loc_val:
                    result["is_secreted"] = True

    # Gene family size estimate from cross-references
    xrefs = data.get("uniProtKBCrossReferences", [])
    panther_count = sum(1 for x in xrefs if x.get("database") == "PANTHER")
    if panther_count > 10:
        result["gene_family_size_estimate"] = "large"
    elif panther_count > 3:
        result["gene_family_size_estimate"] = "medium"
    else:
        result["gene_family_size_estimate"] = "small"

    return result


def analyze_query(resolved_json_path: str, fasta_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Build query context from resolved_input.json and optional FASTA.
    """
    result = {
        "length_aa": 0,
        "estimated_exon_count": 1,
        "has_signal_peptide": False,
        "is_secreted": False,
        "domain_families": [],
        "gene_family_size_estimate": "unknown",
        "protein_name": "",
    }

    # Load resolved_input.json
    resolved = {}
    if os.path.exists(resolved_json_path):
        try:
            with open(resolved_json_path) as fh:
                resolved = json.load(fh)
        except Exception:
            pass

    # Get protein length from FASTA
    fasta = fasta_path or resolved.get("fasta_path", "")
    if fasta and os.path.exists(fasta):
        result["length_aa"] = _parse_fasta_length(fasta)

    # Estimate exon count from protein length (rough heuristic)
    length = result["length_aa"]
    if length > 0:
        if length < 80:
            result["estimated_exon_count"] = 1
        elif length < 200:
            result["estimated_exon_count"] = 2
        elif length < 500:
            result["estimated_exon_count"] = 4
        else:
            result["estimated_exon_count"] = max(5, length // 100)

    # Fetch UniProt features if we have an accession
    source = resolved.get("source", "")
    input_id = resolved.get("input_id", "")
    if source == "uniprot" and input_id:
        features = _fetch_uniprot_features(input_id)
        result.update(features)
    elif resolved.get("protein_name"):
        result["protein_name"] = resolved["protein_name"]

    return result


# ---------------------------------------------------------------------------
# Evolutionary distance estimation
# ---------------------------------------------------------------------------


def estimate_max_distance(
    home_lineage: List[Dict], target_species: List[str]
) -> Tuple[float, List[str], List[str]]:
    """
    Estimate max evolutionary distance in Mya between home species and targets.
    Returns (max_distance_mya, kingdoms_represented, phyla_represented).
    """
    if not target_species:
        return 0, [], []

    home_ranks = {}
    for entry in home_lineage:
        rank = entry.get("rank", "").lower()
        name = entry.get("name", "")
        if rank and name:
            home_ranks[rank] = name

    max_dist = 0
    kingdoms = set()
    phyla = set()

    for sp in target_species:
        sp_info = fetch_species_info(sp)
        kingdoms.add(sp_info.get("kingdom", "Unknown"))
        phyla.add(sp_info.get("phylum", "Unknown"))

        # Find lowest common rank
        sp_ranks = {}
        for entry in sp_info.get("lineage", []):
            rank = entry.get("rank", "").lower()
            name = entry.get("name", "")
            if rank and name:
                sp_ranks[rank] = name

        # Walk from superkingdom down to find the lowest common rank (LCA).
        # The divergence happened one rank below the LCA.
        lca_rank = None
        for rank in ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            home_val = home_ranks.get(rank, "")
            sp_val = sp_ranks.get(rank, "")
            if home_val and sp_val and home_val == sp_val:
                lca_rank = rank  # they still agree at this rank
            elif home_val and sp_val:
                break  # first disagreement — LCA is the previous rank

        # Map LCA rank to approximate divergence time.
        # If they agree at genus, divergence ~ genus-level Mya.
        # If no agreement at any rank, assume superkingdom-level distance.
        divergence_rank = lca_rank if lca_rank else "superkingdom"

        dist = RANK_DISTANCE_MYA.get(divergence_rank, 1000)
        max_dist = max(max_dist, dist)

    return max_dist, sorted(kingdoms), sorted(phyla)


# ---------------------------------------------------------------------------
# Main context builder
# ---------------------------------------------------------------------------


def build_context(
    resolved_json: str,
    home_species: str,
    target_species_str: str,
    fasta_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Build the full context JSON for LLM parameter estimation.
    """
    logger.info("Building LLM context...")

    # 1. Analyze query protein
    query_info = analyze_query(resolved_json, fasta_path)
    logger.info(
        f"Query: {query_info['protein_name'] or 'unknown'}, "
        f"{query_info['length_aa']} aa, "
        f"~{query_info['estimated_exon_count']} exon(s)"
    )

    # 2. Fetch home species info
    home_info = fetch_species_info(home_species)

    # 3. Parse target species list
    targets = []
    if target_species_str:
        targets = [s.strip() for s in target_species_str.split(",") if s.strip()]

    # 4. Estimate evolutionary distance
    max_dist_mya, kingdoms_repr, phyla_repr = estimate_max_distance(
        home_info.get("lineage", []), targets
    )

    # Build target context
    target_context = {
        "kingdoms_represented": kingdoms_repr,
        "phyla_represented": phyla_repr,
        "max_evolutionary_distance_mya": max_dist_mya,
        "includes_plants": "Plantae" in kingdoms_repr,
        "includes_bacteria": "Bacteria" in kingdoms_repr,
        "target_count": len(targets),
    }

    # Remove internal fields from home_info
    home_info_clean = {k: v for k, v in home_info.items() if not k.startswith("_")}
    home_info_clean.pop("lineage", None)

    context = {
        "query": query_info,
        "home_species": home_info_clean,
        "target_context": target_context,
    }

    return context


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Build biological context for LLM parameter estimation"
    )
    parser.add_argument(
        "--resolved_query",
        required=True,
        help="Path to resolved_input.json from RESOLVE_GENE_INPUT",
    )
    parser.add_argument("--home_species", required=True, help="Home species name")
    parser.add_argument(
        "--target_species",
        default="",
        help="Comma-separated target species names (empty = auto-detected)",
    )
    parser.add_argument(
        "--fasta", default=None, help="Path to query FASTA (optional, for length measurement)"
    )
    parser.add_argument("--output", required=True, help="Output context JSON path")
    args = parser.parse_args()

    context = build_context(
        args.resolved_query,
        args.home_species,
        args.target_species,
        args.fasta,
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(context, fh, indent=2)

    logger.info(f"Context written to {args.output}")
    print(json.dumps(context, indent=2))


if __name__ == "__main__":
    main()
