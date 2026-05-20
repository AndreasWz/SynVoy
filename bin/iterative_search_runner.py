#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
import glob
import re
import concurrent.futures
import uuid
import sys
import json
import logging
import contextlib
import time
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
from urllib.parse import unquote

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def has_parasail_available() -> bool:
    """
    Fast runtime probe for parasail availability in the current Python env.
    """
    try:
        import parasail  # type: ignore
        return True
    except Exception:
        return False


def _tail_lines(text: str, n: int = 20) -> str:
    if not text:
        return ""
    return "\n".join(text.strip().splitlines()[-n:])


def _format_mmseqs_failure(proc: subprocess.CompletedProcess) -> str:
    """
    Return compact failure details from an mmseqs subprocess result.
    """
    stderr_tail = _tail_lines(proc.stderr or "", n=20)
    stdout_tail = _tail_lines(proc.stdout or "", n=20)
    return stderr_tail or stdout_tail or "no tool output captured"


def _is_mmseqs_resource_failure(details: str) -> bool:
    """
    Detect common low-memory / prefilter-failure signatures from MMseqs.
    """
    txt = (details or "").lower()
    patterns = (
        "prefilter died",
        "search step died",
        "search died",
        "out of memory",
        "cannot allocate memory",
        "std::bad_alloc",
        "killed",
    )
    return any(p in txt for p in patterns)


# Shared helpers
def str2bool(value: str) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    val = value.strip().lower()
    if val in {"true", "1", "yes", "y", "t"}:
        return True
    if val in {"false", "0", "no", "n", "f"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

# Use our own sequence utilities (no BioPython dependency)
try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )
    from annotate_goi_exons import annotate_exons_from_hit_list, MINIPROT_AVAILABLE
except ImportError:
    # Fallback if not in path - add bin directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        parse_gff, get_feature_id, load_genome, reverse_complement, translate
    )
    from annotate_goi_exons import annotate_exons_from_hit_list, MINIPROT_AVAILABLE

# Import fragment utilities if available
try:
    from fragment_query import generate_fragments, parse_fragment_id, merge_fragment_hits
    FRAGMENT_SUPPORT = True
except ImportError:
    FRAGMENT_SUPPORT = False

# Import PLM (Protein Language Model) search if available
try:
    from plm_search import (
        check_plm_available,
        plm_search_region,
        compute_candidate_similarities,
        precompute_goi_embeddings,
        load_embeddings,
        embed_proteins,
    )
    PLM_IMPORT_OK = True
except ImportError:
    PLM_IMPORT_OK = False

    def check_plm_available():
        return False

# Import structural search (Foldseek + ESMFold) if available
try:
    from structural_search import (
        check_structural_search_available,
        check_esmfold_available,
        check_foldseek_available,
        structural_search_region,
        compute_candidate_structural_similarities,
        prefold_goi_structures,
        save_structure_index,
        load_structure_index,
    )
    STRUCTURAL_IMPORT_OK = True
except ImportError:
    STRUCTURAL_IMPORT_OK = False

    def check_structural_search_available():
        return False

    def check_esmfold_available():
        return False

    def check_foldseek_available():
        return False

# Classification thresholds for GOI confidence assignment.  Populated from
# CLI args in main() so that they can be tuned per-run via nextflow.config.
CLASSIFY_THRESHOLDS = {
    # exon_annotation -> HIGH
    "high_min_exons": 2,
    "high_min_identity": 50.0,
    "high_min_flanking": 2,
    # exon_annotation -> MEDIUM
    "medium_min_identity": 35.0,
    "medium_min_flanking": 1,
    "medium_min_qcov": 0.65,
    # fallback_hit_span -> MEDIUM (flanking-supported)
    "fallback_med_min_flanking": 2,
    "fallback_med_min_qcov": 0.75,
    "fallback_med_min_identity": 60.0,
    # fallback_hit_span -> MEDIUM (strong flanking context)
    "fallback_strong_min_flanking": 5,
    "fallback_strong_min_qcov": 0.25,
    "fallback_strong_min_identity": 35.0,
    # tandem_copy -> MEDIUM vs LOW
    "tandem_min_identity": 40.0,
    # model_status thresholds
    "fragment_max_qcov": 0.4,
    "complete_min_qcov": 0.7,
}

# Family-consistency config (Change A in TP53_IMPLEMENTATION_FIX_PLAN.md).
# `tokens` is a set of uppercase-alphanum tokens derived from the query gene name.
# When populated, every GOI annotation receives a GoiFamilyConsistent attribute.
# When `strict` is True, family-inconsistent fallback/rescued_exon/raw_hit calls
# are downgraded to LOW/ambiguous so they do not masquerade as probable GOI.
FAMILY_CONFIG: Dict[str, Any] = {
    "strict": False,
    "tokens": set(),
    "strict_evidence_types": {"fallback_hit_span", "rescued_exon", "raw_hit"},
}

# Weak-GOI emission filter. Drops LOW-confidence fallback/rescued_exon/raw_hit
# calls whose (identity/100)*qcov is below `min_id_x_qcov`. Motivated by TP53
# runs where ~70% of LOW fallback rows had id*qcov < 0.10 (pure noise) while
# genuine MEDIUM/probable melittin fallback hits never fall below 0.12.
# Default 0.05 drops obvious garbage without touching any validated toxin call.
# Set to 0.0 to disable entirely.
EMISSION_CONFIG: Dict[str, Any] = {
    "min_id_x_qcov": 0.0,
    "weak_evidence_types": {"fallback_hit_span", "rescued_exon", "raw_hit"},
    "_dropped_counter": 0,
}


def _should_skip_weak_goi_emission(
    evidence_type: str,
    identity: float,
    query_cov: Optional[float],
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
) -> bool:
    """Return True if this GOI emission is weak enough to drop before GFF write.

    Only fires on LOW-confidence fallback/rescued_exon/raw_hit. A hit that
    would classify as MEDIUM/HIGH (strong flanking, PLM rescue, structural
    rescue, or adequate qcov+identity) is always kept.
    """
    threshold = float(EMISSION_CONFIG.get("min_id_x_qcov") or 0.0)
    if threshold <= 0.0:
        return False
    if evidence_type not in EMISSION_CONFIG.get("weak_evidence_types", set()):
        return False
    qcov = float(query_cov or 0.0)
    ident = float(identity or 0.0)
    # Safety: if qcov is missing (0), the emission lacks coverage info — keep.
    # Prevents accidentally dropping tandem-copy-like rows that slip through.
    if qcov <= 0.0 or ident <= 0.0:
        return False
    conf, _cls, _reason = _classify_goi_evidence(
        evidence_type=evidence_type,
        identity=ident,
        exon_count=1,
        query_cov=qcov,
        flanking_support=flanking_support,
        embedding_similarity=embedding_similarity,
        structural_similarity=structural_similarity,
    )
    if conf != "LOW":
        return False
    return (ident / 100.0) * qcov < threshold


def _normalize_family_token(value: Any) -> str:
    if not value:
        return ""
    return re.sub(r"[^A-Za-z0-9]", "", str(value)).upper()


def _auto_derive_family_tokens(fasta_path: Optional[str]) -> set:
    """Extract likely gene-name tokens from a UniProt/NCBI FASTA header.

    Looks for ``GN=XYZ`` (UniProt) or falls back to the accession's gene suffix.
    Returns a set of normalized tokens. Empty set if nothing derivable.
    """
    tokens: set = set()
    if not fasta_path or not os.path.exists(fasta_path):
        return tokens
    try:
        with open(fasta_path, "r") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                m = re.search(r"GN=([A-Za-z0-9_\-]+)", header)
                if m:
                    tokens.add(_normalize_family_token(m.group(1)))
                # UniProt sp|ACC|NAME_SPECIES pattern
                m2 = re.match(r"sp\|[^|]+\|([A-Za-z0-9]+)_", header)
                if m2:
                    tokens.add(_normalize_family_token(m2.group(1)))
                break
    except Exception:
        return tokens
    tokens.discard("")
    return tokens


def _check_family_consistency(target_gene: str, target_product: str) -> Tuple[bool, str]:
    tokens = FAMILY_CONFIG.get("tokens") or set()
    if not tokens:
        return True, "family_tokens_unset"
    gene_norm = _normalize_family_token(target_gene)
    product_norm = _normalize_family_token(target_product)
    if not gene_norm and not product_norm:
        return False, "no_target_annotation"
    for tok in tokens:
        if not tok:
            continue
        if gene_norm and tok in gene_norm:
            return True, f"matched_gene:{tok}"
        if product_norm and tok in product_norm:
            return True, f"matched_product:{tok}"
    return False, "no_family_match"


def run_command(cmd):
    subprocess.check_call(cmd)


@contextlib.contextmanager
def maybe_quiet_streams(quiet: bool):
    """
    Suppress noisy third-party stdout/stderr (e.g. miniprot diagnostics)
    when running in low-noise mode.
    """
    if not quiet:
        yield
        return
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield

def normalize_coordinates(start: int, end: int) -> Tuple[int, int]:
    return min(start, end), max(start, end)

def filter_exon_hits(hits: List[Dict[str, Any]], query_len: int,
                     min_query_cov: float, min_alnlen: int) -> List[Dict[str, Any]]:
    """Filter exon hits by query coverage and alignment length."""
    if not hits or query_len <= 0:
        return hits
    kept = []
    for h in hits:
        qstart = h.get('qstart')
        qend = h.get('qend')
        if qstart is None or qend is None:
            continue
        qspan = abs(qend - qstart) + 1
        alnlen = h.get('alnlen', qspan)
        qcov = qspan / query_len if query_len > 0 else 0
        
        # We comment out the strict constraints here because they destroy fragmented exon hits
        # before the more robust miniprot/fallback logic can assemble them.
        # if min_alnlen and alnlen < min_alnlen:
        #    continue
        # if min_query_cov and qcov < min_query_cov:
        #    continue
        kept.append(h)
    return kept

def parse_hits(
    hits_file: str,
    min_identity: float,
    min_length: int,
    evalue_thresh: float,
    query_lengths: Optional[Dict[str, int]] = None,
    short_query_frac: float = 0.60,
    short_query_min: int = 12,
) -> List[Dict[str, Any]]:
    """
    Parse MMseqs2 hits and return a list of hit dictionaries.
    Filters by basic quality metrics.
    Preserves qstart/qend (query protein positions) and strand for exon annotation.
    """
    hits = []
    if not os.path.exists(hits_file):
        return hits

    skipped_lines = 0
    total_lines = 0
    first_error = None
    with open(hits_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            total_lines += 1
            try:
                # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
                # 0      1       2       3       4         5        6       7     8       9     10      11
                if len(parts) < 11:
                    skipped_lines += 1
                    continue

                pident = float(parts[2])
                alnlen = int(parts[3])
                evalue = float(parts[10])

                q_id = parts[0]
                q_len = None
                if query_lengths:
                    q_len = query_lengths.get(q_id)
                    if q_len is None:
                        q_base = extract_base_gene_id(q_id)
                        q_len = query_lengths.get(q_base)

                effective_min_length = min_length
                if q_len and q_len > 0:
                    adaptive_min = max(short_query_min, int(round(q_len * short_query_frac)))
                    effective_min_length = min(min_length, adaptive_min)

                if (evalue <= evalue_thresh and
                    pident >= min_identity and
                    alnlen >= effective_min_length):

                    t_start = int(parts[8])
                    t_end = int(parts[9])
                    start, end = normalize_coordinates(t_start, t_end)
                    # Convert 1-based mmseqs/BLAST coordinates to 0-based half-open
                    # for Python slicing: start-1 becomes 0-based, end stays (exclusive)
                    start -= 1
                    strand = '+' if t_start <= t_end else '-'

                    q_start = int(parts[6])
                    q_end = int(parts[7])

                    bits = float(parts[11]) if len(parts) > 11 else 0.0

                    hits.append({
                        'query': q_id,
                        'target': parts[1], # Chromosome/Scaffold
                        'chrom': parts[1],
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'qstart': min(q_start, q_end),
                        'qend': max(q_start, q_end),
                        'evalue': evalue,
                        'pident': pident,
                        'alnlen': alnlen,
                        'bits': bits
                    })
            except (ValueError, IndexError) as e:
                skipped_lines += 1
                if first_error is None:
                    first_error = f"line {total_lines}: {e}"
                continue

    if skipped_lines > 0:
        logger.warning(
            f"parse_hits({os.path.basename(hits_file)}): skipped {skipped_lines}/{total_lines} "
            f"malformed lines (first: {first_error})"
        )
    return hits

def create_locus_object(query_id, hits):
    chrom = hits[0]['chrom']
    start = min(h['start'] for h in hits)
    end = max(h['end'] for h in hits)
    return {
        'query': query_id,
        'chrom': chrom,
        'start': start,
        'end': end,
        'hits': hits
    }

def extract_base_gene_id(query_id: str) -> str:
    """
    Extract the base gene ID from a query ID that may contain exon info.
    
    Handles formats:
    - "gene-LOC726866" -> "gene-LOC726866"
    - "gene-LOC726866|exon_1" -> "gene-LOC726866"
    - "gene-LOC726866|var1" -> "gene-LOC726866"
    - "gene-LOC726866|GCA_xxx_MP000001" -> "gene-LOC726866"
    """
    # Split on | and take the first part
    parts = query_id.split('|')
    base_id = parts[0]
    
    # Also handle cases where gene ID itself contains underscores but not exon info
    # The exon suffix is always "|exon_N"
    
    return base_id


def is_goi_query_id(query_id: str) -> bool:
    """
    Identify GOI-derived queries in the expanding DB.

    GOI IDs produced by annotate_goi_exons.py are prefixed with `GOI_`
    (full protein, exon entries, tandem copies).
    """
    if not query_id:
        return False
    base_id = extract_base_gene_id(query_id)
    # Backward compatibility: older runs could emit bare "exon_N" IDs
    # for GOI-derived exons (missing GOI_ prefix).
    bare_legacy_goi_exon = bool(re.fullmatch(r'exon_\d+', query_id))
    return (
        query_id.startswith('GOI_') or
        base_id.startswith('GOI_') or
        query_id.startswith('GOI_copy_') or
        bare_legacy_goi_exon
    )


def is_full_length_goi_query_id(query_id: str) -> bool:
    """
    Full-length GOI-like queries used for final gene model inference.
    Excludes exon and synthetic fragment entries.
    """
    if not is_goi_query_id(query_id):
        return False
    if '|exon_' in query_id:
        return False
    if '|frag_' in query_id:
        return False
    return True


def split_hits_into_loci(hits: List[Dict[str, Any]], max_gap: int) -> List[List[Dict[str, Any]]]:
    """
    Split hits into local loci so distant paralogs are not forced into one model.
    """
    if not hits:
        return []

    sorted_hits = sorted(hits, key=lambda h: (h.get('chrom', ''), h.get('gstart', h.get('start', 0))))
    loci: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = [sorted_hits[0]]

    for h in sorted_hits[1:]:
        prev = current[-1]
        same_chrom = h.get('chrom') == prev.get('chrom')
        prev_end = prev.get('gend', prev.get('end', 0))
        curr_start = h.get('gstart', h.get('start', 0))

        if same_chrom and (curr_start - prev_end) <= max_gap:
            current.append(h)
        else:
            loci.append(current)
            current = [h]

    loci.append(current)
    return loci


def build_flanking_query_by_parent(
    db_sequences: Dict[str, Dict[str, str]],
    parent_ids: set
) -> Dict[str, Dict[str, str]]:
    """
    Build one representative flanking query sequence per parent gene.

    Flanking DB entries are often exon-level (`gene|exon_N`). For these, we
    reconstruct a pseudo full-length protein by concatenating exon proteins in
    transcript order (reverse exon order on minus strand).

    Returns dict mapping parent_id -> {id, seq, header, exon_offsets}.
    exon_offsets maps exon_num -> (start_offset, length) in the reconstructed
    protein, enabling downstream qstart/qend remapping.
    """
    if not parent_ids:
        return {}

    # Deduplicate aliased entries (base IDs and explicit IDs can point to same record).
    unique_records: Dict[str, Dict[str, str]] = {}
    for rec in db_sequences.values():
        rid = rec.get('id', '')
        if rid and rid not in unique_records:
            unique_records[rid] = rec

    by_parent: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for rid, rec in unique_records.items():
        parent = extract_base_gene_id(rid)
        if parent in parent_ids and not is_goi_query_id(rid):
            by_parent[parent].append(rec)

    result: Dict[str, Dict[str, str]] = {}
    for parent, recs in by_parent.items():
        # Prefer true full-length record if present.
        full_records = [
            r for r in recs
            if r.get('id', '') == parent and '|exon_' not in r.get('id', '')
        ]
        if full_records:
            best = max(full_records, key=lambda r: len(r.get('seq', '')))
            if best.get('seq'):
                result[parent] = {
                    'id': parent,
                    'seq': best.get('seq', ''),
                    'header': best.get('header', parent),
                    'exon_offsets': {},
                }
                continue

        exon_parts = []
        strand_hint = '+'
        for r in recs:
            rid = r.get('id', '')
            m = re.search(r'\|exon_(\d+)$', rid)
            if not m:
                continue
            exon_num = int(m.group(1))
            seq = r.get('seq', '')
            if not seq:
                continue
            header = r.get('header', '')
            sm = re.search(r'\bstrand=([+-])\b', header)
            if sm:
                strand_hint = sm.group(1)
            exon_parts.append((exon_num, seq))

        if exon_parts:
            exon_parts.sort(key=lambda x: x[0], reverse=(strand_hint == '-'))
            # Build offset map: exon_num -> (start_offset_in_reconstructed, length)
            exon_offsets = {}
            offset = 0
            for exon_num, seq in exon_parts:
                exon_offsets[exon_num] = (offset, len(seq))
                offset += len(seq)
            recon_seq = ''.join(seq for _, seq in exon_parts)
            if recon_seq:
                result[parent] = {
                    'id': parent,
                    'seq': recon_seq,
                    'header': f"{parent}|reconstructed_exons={len(exon_parts)}|strand={strand_hint}",
                    'exon_offsets': exon_offsets,
                }
                continue

        # Last fallback: use longest available record.
        best = max(recs, key=lambda r: len(r.get('seq', '')))
        if best.get('seq'):
            result[parent] = {
                'id': parent,
                'seq': best.get('seq', ''),
                'header': best.get('header', parent),
                'exon_offsets': {},
            }

    return result


def _parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for kv in (attr_field or "").split(";"):
        if not kv or "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def _select_parent_id(attrs: Dict[str, str], model_id: str) -> str:
    """
    Select a stable parent identifier from heterogeneous GFF attributes.

    Supports SynVoy, NCBI/Ensembl-style, and custom tags.
    """
    for key in [
        "SynVoy_Parent",
        "ParentProtein",
        "protein_id",
        "gene_id",
        "gene",
        "locus_tag",
        "Name",
    ]:
        value = (attrs.get(key) or "").strip()
        if not value:
            continue
        # If multiple values are provided, keep the first token.
        value = value.split(",")[0].strip()
        if value:
            return value
    return extract_base_gene_id(model_id)


def _is_generic_label(label: str) -> bool:
    """Heuristic for non-informative locus-tag style labels."""
    if not label:
        return True
    txt = str(label).strip()
    if not txt:
        return True
    if txt.startswith("gene-"):
        txt = txt[5:]
    if re.match(r"^[A-Za-z]+\d*_\d+$", txt):
        return True
    if re.match(r"^[A-Z]{2,6}\d{0,3}_\d+$", txt):
        return True
    if re.match(r"^LOC\d+$", txt, re.IGNORECASE):
        return True
    return False


def _safe_gff_value(value: Any) -> str:
    """Sanitize attribute values so GFF parsing remains robust."""
    txt = str(value) if value is not None else ""
    txt = txt.replace(";", ",").replace("\t", " ").replace("\n", " ").strip()
    return txt


def _compose_gff_attrs(base_attrs: Dict[str, Any],
                       native_annot: Optional[Dict[str, str]] = None) -> str:
    """
    Compose a GFF attribute field from base attrs plus optional native-annotation attrs.
    """
    merged: Dict[str, Any] = {}
    for key, value in base_attrs.items():
        if value is None:
            continue
        sval = _safe_gff_value(value)
        if sval:
            merged[key] = sval

    label = ""
    product = ""
    if native_annot:
        if native_annot.get("label"):
            merged["TargetGene"] = _safe_gff_value(native_annot["label"])
            label = native_annot["label"]
        if native_annot.get("product"):
            merged["TargetProduct"] = _safe_gff_value(native_annot["product"])
            product = native_annot["product"]
        if native_annot.get("feature_id"):
            merged["TargetID"] = _safe_gff_value(native_annot["feature_id"])

    if base_attrs.get("SynVoyRole") == "goi" and FAMILY_CONFIG.get("tokens"):
        consistent, reason = _check_family_consistency(label, product)
        merged["GoiFamilyConsistent"] = "true" if consistent else "false"
        merged["GoiFamilyReason"] = reason
        # Strict downgrade only fires on informative mismatches (annotated locus
        # whose name disagrees with family tokens). Unannotated regions — common
        # for novel toxin/venom discovery where orthologs are de novo — are NOT
        # downgraded: synteny + sequence evidence must stand on their own there.
        if (FAMILY_CONFIG.get("strict")
                and not consistent
                and reason == "no_family_match"):
            ev = base_attrs.get("EvidenceType") or base_attrs.get("Type") or ""
            if ev in FAMILY_CONFIG.get("strict_evidence_types", set()):
                merged["Confidence"] = "LOW"
                merged["GOIClass"] = "ambiguous_goi_family_member"
                prev = merged.get("InferenceReason", "")
                merged["InferenceReason"] = (
                    f"{prev}|strict_family_downgrade" if prev else "strict_family_downgrade"
                )

    return ";".join(f"{k}={v}" for k, v in merged.items())


def _format_attr_float(value: Optional[float], digits: int = 3) -> Optional[str]:
    if value is None:
        return None
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return None


def _synteny_context_label(flanking_support: int) -> str:
    flanking_support = max(0, int(flanking_support or 0))
    if flanking_support >= 3:
        return "strong_flanking_support"
    if flanking_support >= 1:
        return "weak_flanking_support"
    return "no_flanking_support"


def _model_status(
    query_cov: float,
    exon_count: int = 1,
    evidence_type: str = "",
) -> str:
    """
    Label a gene model's completeness independent of confidence.

    Returns one of:
      - ``complete``  — query coverage >= complete_min_qcov and multi-exon
        (or single-exon gene that is expected to be single-exon, e.g.
        tandem_copy)
      - ``partial``   — between fragment and complete
      - ``fragment``  — query coverage < fragment_max_qcov or rescued_exon /
        raw_hit
    """
    ct = CLASSIFY_THRESHOLDS
    qcov = float(query_cov or 0.0)
    exons = max(1, int(exon_count or 1))

    if evidence_type in {"rescued_exon", "raw_hit"}:
        return "fragment"
    if qcov < ct["fragment_max_qcov"]:
        return "fragment"
    if qcov >= ct["complete_min_qcov"] and (exons >= 2 or evidence_type == "tandem_copy"):
        return "complete"
    return "partial"


def _classify_goi_evidence(
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
) -> Tuple[str, str, str]:
    """
    Assign a conservative confidence/class label to GOI-derived candidates.

    The goal is not to prove orthology here, but to prevent fallback-heavy
    output from masquerading as confident GOI evidence downstream.

    embedding_similarity (0-1 cosine) from ProtT5 can BOOST confidence when
    the structural/functional signal is strong even if sequence identity is low.
    structural_similarity (0-1 TM-score) from ESMFold+Foldseek can BOOST
    confidence when 3D structure matches despite extreme sequence divergence.
    Neither signal ever reduces confidence — sequence methods already confirmed
    a match.
    """
    identity = float(identity or 0.0)
    exon_count = max(1, int(exon_count or 1))
    qcov = float(query_cov or 0.0)
    context = _synteny_context_label(flanking_support)

    ct = CLASSIFY_THRESHOLDS

    if evidence_type == "tandem_copy":
        if identity >= ct["tandem_min_identity"]:
            return "MEDIUM", "tandem_goi_copy", "goi_tandem_copy_detected"
        return "LOW", "tandem_goi_copy", "goi_tandem_copy_low_identity"

    if evidence_type == "exon_annotation":
        if (exon_count >= ct["high_min_exons"]
                and identity >= ct["high_min_identity"]
                and flanking_support >= ct["high_min_flanking"]):
            return "HIGH", "confident_goi", "multi_exon_model_with_flanking_support"
        if identity >= ct["medium_min_identity"] and (
                flanking_support >= ct["medium_min_flanking"]
                or qcov >= ct["medium_min_qcov"]):
            confidence = "MEDIUM"
            reason = "modeled_goi_with_partial_support"
            # PLM boost: MEDIUM → HIGH when embedding strongly agrees
            if (embedding_similarity is not None
                    and embedding_similarity >= ct.get("plm_high_threshold", 0.85)
                    and flanking_support >= ct["high_min_flanking"]):
                confidence = "HIGH"
                reason += "_plm_boosted"
            # Structural boost: MEDIUM → HIGH when TM-score strongly agrees
            if (confidence == "MEDIUM"
                    and structural_similarity is not None
                    and structural_similarity >= ct.get("structural_high_threshold", 0.7)
                    and flanking_support >= ct["high_min_flanking"]):
                confidence = "HIGH"
                reason += "_structural_boosted"
            return confidence, "probable_goi", reason
        # LOW from sequence — but PLM may rescue
        confidence = "LOW"
        goi_class = "ambiguous_goi_family_member"
        reason = "modeled_goi_but_family_context_is_weak"
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            confidence = "MEDIUM"
            goi_class = "probable_goi"
            reason = "modeled_goi_plm_rescued"
        # Structural rescue: LOW → MEDIUM when TM-score confirms structural match
        if (confidence == "LOW"
                and structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            confidence = "MEDIUM"
            goi_class = "probable_goi"
            reason = "modeled_goi_structural_rescued"
        return confidence, goi_class, reason

    if evidence_type == "fallback_hit_span":
        if (flanking_support >= ct["fallback_med_min_flanking"]
                and qcov >= ct["fallback_med_min_qcov"]
                and identity >= ct["fallback_med_min_identity"]):
            return "MEDIUM", "probable_goi", "fallback_span_supported_by_flanking_context"
        # Strong syntenic evidence upgrades even low-qcov hits.
        # Short secreted peptides (e.g. melittin, ~70 aa) can only reach ~45% qcov
        # because the signal peptide exon is too diverged to align; the flanking
        # gene context provides orthology evidence beyond per-alignment statistics.
        if (flanking_support >= ct["fallback_strong_min_flanking"]
                and (qcov >= ct["fallback_strong_min_qcov"]
                     or identity >= ct["fallback_strong_min_identity"])):
            return "MEDIUM", "probable_goi", "fallback_span_with_strong_flanking_support"
        # PLM rescue for fallback hits
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "fallback_span_plm_rescued"
        # Structural rescue for fallback hits
        if (structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "fallback_span_structural_rescued"
        return "LOW", "ambiguous_goi_family_member", "fallback_span_only"

    if evidence_type == "rescued_exon":
        # PLM can rescue isolated exons when embedding strongly matches
        if (embedding_similarity is not None
                and embedding_similarity >= ct.get("plm_medium_threshold", 0.7)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "rescued_exon_plm_boosted"
        # Structural rescue for isolated exons
        if (structural_similarity is not None
                and structural_similarity >= ct.get("structural_medium_threshold", 0.5)
                and flanking_support >= 1):
            return "MEDIUM", "probable_goi", "rescued_exon_structural_boosted"
        return "LOW", "ambiguous_goi_family_member", "isolated_rescued_exon"

    if evidence_type == "raw_hit":
        return "LOW", "ambiguous_goi_family_member", "single_raw_hit_only"

    return "LOW", "ambiguous_goi_family_member", f"unclassified_{evidence_type or 'goi'}"


def _classify_flanking_evidence(
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
) -> Tuple[str, str]:
    identity = float(identity or 0.0)
    exon_count = max(1, int(exon_count or 1))
    qcov = float(query_cov or 0.0)

    if evidence_type in {"flanking_miniprot", "rearranged_flanking"}:
        if exon_count >= 2 or identity >= 55.0:
            return "HIGH", "modeled_flanking_gene"
        return "MEDIUM", "single_exon_flanking_model"

    if evidence_type in {"flanking_hit_span", "rearranged_flanking_fallback"}:
        if qcov >= 0.65 and identity >= 55.0:
            return "MEDIUM", "coarse_flanking_span_with_good_support"
        return "LOW", "coarse_flanking_span_only"

    return "MEDIUM", f"unclassified_{evidence_type or 'flanking'}"


def _goi_feature_attrs(
    base_attrs: Dict[str, Any],
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    flanking_support: int = 0,
    embedding_similarity: Optional[float] = None,
    structural_similarity: Optional[float] = None,
) -> Dict[str, Any]:
    confidence, goi_class, reason = _classify_goi_evidence(
        evidence_type=evidence_type,
        identity=identity,
        exon_count=exon_count,
        query_cov=query_cov,
        flanking_support=flanking_support,
        embedding_similarity=embedding_similarity,
        structural_similarity=structural_similarity,
    )
    attrs = dict(base_attrs)
    attrs.setdefault("Identity", f"{float(identity or 0.0):.1f}")
    attrs["SynVoyRole"] = "goi"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
    attrs["GOIClass"] = goi_class
    attrs["ModelStatus"] = _model_status(query_cov, exon_count, evidence_type)
    attrs["SyntenyContext"] = _synteny_context_label(flanking_support)
    attrs["BlockFlankingSupport"] = str(max(0, int(flanking_support or 0)))
    attrs["InferenceReason"] = reason
    if query_cov is not None:
        attrs["QueryCoverage"] = _format_attr_float(query_cov)
    if exon_count:
        attrs.setdefault("Exons", str(int(exon_count)))
    if embedding_similarity is not None:
        attrs["EmbeddingSimilarity"] = f"{embedding_similarity:.3f}"
    if structural_similarity is not None:
        attrs["StructuralSimilarity"] = f"{structural_similarity:.3f}"
    return attrs


def _flanking_feature_attrs(
    base_attrs: Dict[str, Any],
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    context: Optional[str] = None,
) -> Dict[str, Any]:
    confidence, reason = _classify_flanking_evidence(
        evidence_type=evidence_type,
        identity=identity,
        exon_count=exon_count,
        query_cov=query_cov,
    )
    attrs = dict(base_attrs)
    attrs.setdefault("Identity", f"{float(identity or 0.0):.1f}")
    attrs["SynVoyRole"] = "flanking"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
    attrs["ModelStatus"] = _model_status(query_cov, exon_count, evidence_type)
    attrs["InferenceReason"] = reason
    if context:
        attrs["SyntenyContext"] = context
    if query_cov is not None:
        attrs["QueryCoverage"] = _format_attr_float(query_cov)
    if exon_count:
        attrs.setdefault("Exons", str(int(exon_count)))
    return attrs


def find_native_annotation_path(genome_path: str) -> Optional[str]:
    """
    Locate the native annotation file next to a genome FASTA.
    Supports .gff/.gff3 (+ gz variants).
    """
    p = Path(genome_path)
    candidates = [
        p.with_suffix(".gff"),
        p.with_suffix(".gff3"),
        Path(str(p) + ".gff"),
        Path(str(p) + ".gff3"),
        p.with_suffix(".gff.gz"),
        p.with_suffix(".gff3.gz"),
        Path(str(p) + ".gff.gz"),
        Path(str(p) + ".gff3.gz"),
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


def load_native_annotation_index(gff_path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """
    Build a per-chromosome interval index from native target annotations.
    """
    if not gff_path or not os.path.exists(gff_path):
        return {}

    entries_by_chrom: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    accepted_types = {"gene", "mRNA", "mrna", "transcript"}

    for feat in parse_gff(gff_path):
        ftype = feat.get("type")
        if ftype not in accepted_types:
            continue

        attrs = feat.get("attributes", {}) or {}
        label = (
            attrs.get("gene")
            or attrs.get("Name")
            or attrs.get("locus_tag")
            or attrs.get("gene_id")
            or attrs.get("ID")
            or ""
        )
        product = attrs.get("product", "") or ""
        feature_id = attrs.get("ID") or attrs.get("locus_tag") or ""

        label = unquote(str(label)).strip()
        product = unquote(str(product)).strip()
        feature_id = unquote(str(feature_id)).strip()

        if not label and not product and not feature_id:
            continue

        try:
            s = int(feat.get("start", 0))
            e = int(feat.get("end", 0))
        except (ValueError, TypeError):
            continue
        if e < s:
            s, e = e, s

        entries_by_chrom[str(feat.get("seqid", ""))].append({
            "start": s,
            "end": e,
            "strand": feat.get("strand", "."),
            "label": label,
            "product": product,
            "feature_id": feature_id,
        })

    index: Dict[str, Dict[str, Any]] = {}
    for chrom, entries in entries_by_chrom.items():
        if not chrom:
            continue
        entries.sort(key=lambda x: (x["start"], x["end"]))
        index[chrom] = {
            "entries": entries,
            "starts": [e["start"] for e in entries],
        }
    return index


def lookup_native_annotation(native_index: Dict[str, Dict[str, Any]],
                             chrom: str,
                             start: int,
                             end: int,
                             strand: Optional[str] = None,
                             max_distance: int = 5000) -> Optional[Dict[str, str]]:
    """
    Find the best native annotation overlapping (or very near) an interval.
    Coordinates are expected as 1-based inclusive.
    """
    if not native_index or not chrom:
        return None
    chrom_data = native_index.get(chrom)
    if not chrom_data:
        return None

    s, e = int(start), int(end)
    if e < s:
        s, e = e, s

    entries = chrom_data["entries"]
    starts = chrom_data["starts"]
    i = bisect_right(starts, e + max_distance)
    if i <= 0:
        return None

    best = None
    best_score = None
    j = i - 1
    while j >= 0:
        ent = entries[j]
        if ent["end"] < s - max_distance:
            break
        ov = max(0, min(e, ent["end"]) - max(s, ent["start"]) + 1)
        if ov > 0:
            dist = 0
        else:
            dist = min(abs(s - ent["end"]), abs(e - ent["start"]))
            if dist > max_distance:
                j -= 1
                continue

        strand_score = 1 if strand and strand in {"+", "-"} and ent.get("strand") == strand else 0
        label_quality = 0 if _is_generic_label(ent.get("label", "")) else 1
        score = (
            1 if ov > 0 else 0,
            ov,
            -dist,
            strand_score,
            label_quality,
            (ent["end"] - ent["start"] + 1),
        )
        if best_score is None or score > best_score:
            best_score = score
            best = ent
        j -= 1

    if not best:
        return None

    return {
        "label": best.get("label", ""),
        "product": best.get("product", ""),
        "feature_id": best.get("feature_id", ""),
    }


def _flanking_chain_consistent(cds_intervals: List[Tuple[int, int]], strand: str) -> bool:
    """
    Validate exon chain direction from emitted CDS order for one transcript.

    '+' strand expects left-to-right exon centers.
    '-' strand expects right-to-left exon centers.
    """
    if len(cds_intervals) <= 1:
        return True

    centers = [0.5 * (s + e) for s, e in cds_intervals]
    if strand == "-":
        return all(centers[i] >= centers[i + 1] for i in range(len(centers) - 1))
    return all(centers[i] <= centers[i + 1] for i in range(len(centers) - 1))


def _split_models_into_loci(
    model_infos: List[Dict[str, Any]],
    locus_gap_bp: int = 50000
) -> List[List[Dict[str, Any]]]:
    """
    Split models into genomic loci for a single parent ID.

    Locus definition:
    - same chromosome
    - same strand
    - model start within `locus_gap_bp` of current locus end
    """
    if not model_infos:
        return []

    ordered = sorted(
        model_infos,
        key=lambda m: (
            m.get("chrom", ""),
            m.get("strand", "+"),
            m.get("start", 0),
            m.get("end", 0),
        ),
    )

    loci: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = []
    current_chrom = ""
    current_strand = "+"
    current_end = -1

    for info in ordered:
        chrom = info.get("chrom", "")
        strand = info.get("strand", "+")
        start = info.get("start", 0)
        end = info.get("end", 0)

        if not current:
            current = [info]
            current_chrom = chrom
            current_strand = strand
            current_end = end
            continue

        same_locus = (
            chrom == current_chrom and
            strand == current_strand and
            start <= (current_end + locus_gap_bp)
        )
        if same_locus:
            current.append(info)
            current_end = max(current_end, end)
        else:
            loci.append(current)
            current = [info]
            current_chrom = chrom
            current_strand = strand
            current_end = end

    if current:
        loci.append(current)
    return loci


def _longest_monotonic_query_chain(
    ordered_hits: List[Dict[str, Any]],
    strand: str,
) -> List[Dict[str, Any]]:
    """
    Keep the longest genomic-order/query-order-consistent hit chain.

    `ordered_hits` must already be sorted by genomic start.
    For '+' strands, query centers should increase with genomic position.
    For '-' strands, query centers should decrease with genomic position.
    """
    n = len(ordered_hits)
    if n <= 1:
        return ordered_hits

    qcenters = []
    for h in ordered_hits:
        qs = h.get("qstart", 0)
        qe = h.get("qend", 0)
        qcenters.append(0.5 * (min(qs, qe) + max(qs, qe)))

    dp = [1] * n
    prev = [-1] * n
    best_idx = 0

    for i in range(n):
        for j in range(i):
            ok = qcenters[j] <= qcenters[i] if strand != "-" else qcenters[j] >= qcenters[i]
            if ok and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j
        if dp[i] > dp[best_idx]:
            best_idx = i

    chain_idx = []
    cur = best_idx
    while cur != -1:
        chain_idx.append(cur)
        cur = prev[cur]
    chain_idx.reverse()
    return [ordered_hits[i] for i in chain_idx]


def deduplicate_flanking_models(
    gene_records: List[Dict[str, Any]],
    gff_lines: List[str],
    genome_name: str = "",
    locus_gap_bp: int = 50000,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """
    Keep one best flanking model per parent ID per genomic locus.

    GOI annotations are untouched. Flanking best model is chosen by:
    1) strand-consistent exon chain
    2) miniprot model over hit-span fallback
    3) more exons
    4) longer total CDS
    5) longer transcript span
    6) higher identity
    """
    if not gff_lines:
        return gene_records, gff_lines

    flanking_sources = {"flanking_annotation", "flanking_hits"}
    models: Dict[str, Dict[str, Any]] = {}
    cds_by_model: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        source = parts[1]
        ftype = parts[2]
        if source not in flanking_sources:
            continue
        attrs = _parse_gff_attributes(parts[8])

        if ftype == "mRNA":
            model_id = attrs.get("ID")
            if not model_id:
                continue
            parent_id = _select_parent_id(attrs, model_id)
            try:
                identity = float(parts[5]) if parts[5] != "." else 0.0
            except (ValueError, TypeError):
                identity = 0.0
            try:
                start = int(parts[3])
                end = int(parts[4])
            except (ValueError, TypeError):
                continue
            models[model_id] = {
                "id": model_id,
                "parent": parent_id,
                "source": source,
                "chrom": parts[0],
                "strand": parts[6],
                "start": start,
                "end": end,
                "identity": identity,
            }
        elif ftype == "CDS":
            parent_model = attrs.get("Parent")
            if not parent_model:
                continue
            try:
                cds_by_model[parent_model].append((int(parts[3]), int(parts[4])))
            except (ValueError, TypeError):
                continue

    if not models:
        return gene_records, gff_lines

    models_by_parent: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for info in models.values():
        models_by_parent[info["parent"]].append(info)

    keep_model_ids = set()
    for parent_id, parent_models in models_by_parent.items():
        loci = _split_models_into_loci(parent_models, locus_gap_bp=max(0, int(locus_gap_bp)))
        for locus_models in loci:
            best = None
            best_score = None
            for info in locus_models:
                model_id = info["id"]
                cds_intervals = cds_by_model.get(model_id, [])
                exon_count = len(cds_intervals)
                cds_total_nt = sum(max(0, e - s + 1) for s, e in cds_intervals)
                span_nt = max(0, info["end"] - info["start"] + 1)
                consistent = _flanking_chain_consistent(cds_intervals, info.get("strand", "+"))
                source_rank = 1 if info.get("source") == "flanking_annotation" else 0

                score = (
                    1 if consistent else 0,
                    source_rank,
                    exon_count,
                    cds_total_nt,
                    span_nt,
                    info.get("identity", 0.0),
                )
                if best_score is None or score > best_score:
                    best_score = score
                    best = model_id

            if best:
                keep_model_ids.add(best)

    if len(keep_model_ids) == len(models):
        return gene_records, gff_lines

    filtered_gff: List[str] = []
    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            filtered_gff.append(line)
            continue
        source = parts[1]
        ftype = parts[2]
        if source not in flanking_sources:
            filtered_gff.append(line)
            continue
        attrs = _parse_gff_attributes(parts[8])
        if ftype == "mRNA":
            model_id = attrs.get("ID", "")
            if model_id in keep_model_ids:
                filtered_gff.append(line)
        elif ftype == "CDS":
            model_id = attrs.get("Parent", "")
            if model_id in keep_model_ids:
                filtered_gff.append(line)
        else:
            filtered_gff.append(line)

    filtered_gene_records = []
    for rec in gene_records:
        rec_id = rec.get("id", "")
        if rec_id in models and rec_id not in keep_model_ids:
            continue
        filtered_gene_records.append(rec)

    dropped = len(models) - len(keep_model_ids)
    if dropped > 0:
        gname = genome_name or "genome"
        logger.info(
            f"[{gname}] Flanking model dedup (per-locus): kept {len(keep_model_ids)}/{len(models)} "
            f"models across {len(models_by_parent)} parent IDs (dropped {dropped} duplicates)."
        )

    return filtered_gene_records, filtered_gff


def collapse_flanking_cds_to_gene_span(gff_lines: List[str]) -> List[str]:
    """
    Replace flanking multi-CDS models with a single span CDS per transcript.

    This keeps flanking annotations gene-centric for downstream plotting while
    preserving GOI CDS detail.
    """
    if not gff_lines:
        return gff_lines

    flanking_sources = {"flanking_annotation", "flanking_hits"}
    mRNAs = []
    output = []

    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            output.append(line)
            continue
        source = parts[1]
        ftype = parts[2]
        attrs = _parse_gff_attributes(parts[8])
        if source in flanking_sources and ftype == "mRNA":
            model_id = attrs.get("ID")
            if model_id:
                try:
                    mRNAs.append({
                        "chrom": parts[0],
                        "source": source,
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "id": model_id,
                    })
                except (ValueError, IndexError) as e:
                    logger.warning(f"Skipping malformed GFF mRNA line for {model_id}: {e}")
            output.append(line)
        elif source in flanking_sources and ftype == "CDS":
            # Drop detailed flanking CDS entries; replaced by one span CDS below.
            continue
        else:
            output.append(line)

    for m in mRNAs:
        output.append(
            f"{m['chrom']}\t{m['source']}\tCDS\t{m['start']}\t{m['end']}\t.\t{m['strand']}\t0\t"
            f"ID={m['id']}_CDS1;Parent={m['id']};Type=flanking_gene_span"
        )

    return output


def identify_synteny_blocks(hits, max_intron=20000, cluster_distance=50000):
    """
    Identify all synteny blocks from hits.
    
    Key Logic:
    - Group hits by query gene (handling multi-exon genes AND exon-level queries)
    - Cluster loci that are close together (likely same gene/region)
    - Score blocks by number of unique flanking genes found
    
    Updated to handle exon-level queries like "gene-LOC726866|exon_1"
    
    Args:
        hits: List of hit dictionaries
        max_intron: Maximum distance between exons of same gene (bp)
        cluster_distance: Maximum distance to cluster genes into synteny block (bp)
    
    Returns:
        List of dictionaries (blocks), sorted by score (descending).
    """
    if not hits:
        return []
        
    # --- Step 1: Group hits by Query Gene (Base ID) ---
    # Multiple hits from same query = different exons or duplicates
    # With exon_mode, query may be "gene-LOC726866|exon_1" -> base is "gene-LOC726866"
    hits_by_query = defaultdict(list)
    for h in hits:
        # Extract base query ID (handles both |var and |exon_ suffixes)
        base_query = extract_base_gene_id(h['query'])
        hits_by_query[base_query].append(h)
        
    # --- Step 2: Define Gene Loci per Query ---
    # Each query gene may have multiple loci (paralogs/duplications)
    # Hits close together (<max_intron) = same locus (multi-exon gene)
    all_loci = []
    for query_id, q_hits in hits_by_query.items():
        # Sort hits by genomic position
        q_hits.sort(key=lambda x: (x['chrom'], x['start']))
        
        current_locus_hits = []
        for h in q_hits:
            if not current_locus_hits:
                current_locus_hits.append(h)
                continue
            
            last_hit = current_locus_hits[-1]
            
            # Same chromosome and close enough = same locus (exons)
            if (h['chrom'] == last_hit['chrom'] and 
                h['start'] - last_hit['end'] < max_intron):
                current_locus_hits.append(h)
            else:
                # Start new locus
                all_loci.append(create_locus_object(query_id, current_locus_hits))
                current_locus_hits = [h]
        
        # Don't forget last locus
        if current_locus_hits:
            all_loci.append(create_locus_object(query_id, current_locus_hits))

    # --- Step 3: Cluster Loci into Synteny Blocks ---
    # Loci from different genes that are close = synteny block
    all_loci.sort(key=lambda x: (x['chrom'], x['start']))
    if not all_loci: 
        return []
        
    synteny_blocks = []
    current_block = [all_loci[0]]
    
    for locus in all_loci[1:]:
        last_locus = current_block[-1]
        
        # Same chromosome and within clustering distance
        if (locus['chrom'] == last_locus['chrom'] and 
            locus['start'] - last_locus['end'] < cluster_distance):
            current_block.append(locus)
        else:
            synteny_blocks.append(current_block)
            current_block = [locus]
    synteny_blocks.append(current_block)
    
    # --- Step 4: Format and Sort Blocks ---
    final_blocks = []
    for block in synteny_blocks:
        chrom = block[0]['chrom']
        start = min(l['start'] for l in block)
        end = max(l['end'] for l in block)
        genes_list = list(set(extract_base_gene_id(l['query']) for l in block))
        
        final_blocks.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'genes_count': len(genes_list),
            'loci_count': len(block),
            'genes': genes_list
        })

    # Sort by 'genes_count' descending
    final_blocks.sort(key=lambda x: x['genes_count'], reverse=True)
    
    return final_blocks


def identify_best_synteny_block(hits, max_intron=20000, cluster_distance=50000):
    """
    Backwards-compatible wrapper kept for tests/importers.
    """
    blocks = identify_synteny_blocks(hits, max_intron=max_intron, cluster_distance=cluster_distance)
    return blocks[0] if blocks else None

def calculate_adaptive_padding(hits: List[Dict[str, Any]], best_region: Dict[str, Any],
                               default: int = 100000, min_pad: int = 50000, max_pad: int = 200000) -> int:
    """
    Calculate region padding based on gene spacing in hits.
    Returns padding distance in base pairs.
    """
    # Filter hits to the best region's chromosome
    region_hits = [h for h in hits if h['chrom'] == best_region['chrom']]
    
    if len(region_hits) < 2:
        return default
    
    # Sort by position
    sorted_hits = sorted(region_hits, key=lambda h: h['start'])
    
    # Calculate inter-gene gaps
    gaps = []
    for i in range(len(sorted_hits) - 1):
        gap = sorted_hits[i+1]['start'] - sorted_hits[i]['end']
        if gap > 0:  # Only positive gaps
            gaps.append(gap)
    
    if not gaps:
        return default
    
    # Average gap * 2 (to cover one gene on each side)
    avg_gap = sum(gaps) / len(gaps)
    adaptive_padding = int(avg_gap * 2)
    
    # Clamp to reasonable range
    final_padding = max(min_pad, min(max_pad, adaptive_padding))
    
    return final_padding


def estimate_cluster_distance(genome_file: str, gff_file: Optional[str] = None, default_dist: int = 50000) -> int:
    """
    Estimate gene density to adjust cluster_distance intelligently.
    
    Strategy:
    1. If GFF provided: Calculate actual inter-gene distances
    2. Else: Use genome size heuristic (improved)
    3. Return 2-3x median inter-gene distance as cluster threshold
    """
    
    # Method 1: Use GFF if available (most accurate)
    if gff_file and os.path.exists(gff_file) and gff_file != "NO_GFF":
        try:
            # parse_gff returns a generator; consume it into a list and
            # filter to gene-level features only for density estimation.
            genes = [f for f in parse_gff(gff_file, feature_types=['gene'])]
            if len(genes) > 10:  # Need reasonable sample size
                # Sort by chromosome and position
                by_chrom = defaultdict(list)
                for gene in genes:
                    by_chrom[gene['seqid']].append(gene['start'])
                
                # Calculate inter-gene distances per chromosome
                all_distances = []
                for chrom, positions in by_chrom.items():
                    sorted_pos = sorted(positions)
                    for i in range(len(sorted_pos) - 1):
                        dist = sorted_pos[i+1] - sorted_pos[i]
                        if dist > 0:  # Skip overlapping genes
                            all_distances.append(dist)
                
                if all_distances:
                    # Use median distance * 2.5 as clustering threshold
                    all_distances.sort()
                    median_dist = all_distances[len(all_distances) // 2]
                    cluster_distance = int(median_dist * 2.5)
                    # Clamp to reasonable range
                    cluster_distance = max(10000, min(200000, cluster_distance))
                    print(f"Estimated cluster distance from GFF: {cluster_distance} bp "
                          f"(median inter-gene: {median_dist} bp)", file=sys.stderr)
                    return cluster_distance
        except Exception as e:
            print(f"Warning: Could not parse GFF for gene density: {e}", file=sys.stderr)
    
    # Method 2: Improved genome size heuristic
    try:
        size = os.path.getsize(genome_file)

        # Gzipped FASTA is ~25% of real genome size → scale up to avoid wrong bin
        if genome_file.endswith('.gz') or genome_file.endswith('.gzip'):
            size = size * 4
        
        # More refined heuristics based on typical genomes
        if size < 5_000_000:  # < 5MB: Bacteria/Archaea
            return 15000  # Dense gene packing
        elif size < 20_000_000:  # 5-20MB: Large bacteria, fungi
            return 25000
        elif size < 100_000_000:  # 20-100MB: Small eukaryotes
            return 40000
        elif size < 500_000_000:  # 100-500MB: Insects, small vertebrates
            return 70000
        elif size < 2_000_000_000:  # 0.5-2GB: Mammals, birds
            return 100000
        else:  # > 2GB: Plants, large genomes
            return 150000
    except Exception as e:
        logger.warning(f"Could not auto-detect cluster distance from genome size ({genome_file}): {e}. Using default {default_dist} bp.")

    return default_dist

def run_augmented_mmseqs_with_retries(
    variants_fasta: str,
    region_fasta: str,
    aug_hits_file: str,
    aug_tmp_dir: str,
    args,
    threads: int,
    genome_name: str,
) -> Tuple[bool, str]:
    """
    Run augmented MMseqs search with low-memory retries.

    Unlike the primary whole-genome search, this runs on per-block regions, but
    still needs retry logic because strict split-memory limits (e.g. 1G) can
    fail on translated ORF DB construction for certain blocks.
    """
    base_attempts = [
        {
            "label": "aug_primary",
            "threads": max(1, int(threads)),
            "sens": float(args.mmseqs_sens),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "aug_lowmem_auto",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "aug_lowmem_2g",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": "2G",
        },
        {
            "label": "aug_lowmem_verysafe",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 6.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
    ]

    attempts = []
    seen = set()
    for a in base_attempts:
        key = (a["threads"], a["sens"], a["split"])
        if key in seen:
            continue
        seen.add(key)
        attempts.append(a)

    os.makedirs(aug_tmp_dir, exist_ok=True)
    last_details = ""

    # Create query and target databases ONCE.
    query_db = os.path.join(aug_tmp_dir, "queryDB")
    target_db = os.path.join(aug_tmp_dir, "targetDB")
    fmt_output = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    proc_qdb = subprocess.run(
        ["mmseqs", "createdb", variants_fasta, query_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_qdb.returncode != 0:
        return False, _format_mmseqs_failure(proc_qdb)

    proc_tdb = subprocess.run(
        ["mmseqs", "createdb", region_fasta, target_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_tdb.returncode != 0:
        return False, _format_mmseqs_failure(proc_tdb)

    for i, attempt in enumerate(attempts, start=1):
        if os.path.exists(aug_hits_file):
            try:
                os.remove(aug_hits_file)
            except OSError:
                pass

        attempt_tmp = os.path.join(aug_tmp_dir, f"attempt_{i}")
        os.makedirs(attempt_tmp, exist_ok=True)
        result_db = os.path.join(attempt_tmp, "resultDB")

        search_cmd = [
            "mmseqs", "search",
            query_db, target_db, result_db, attempt_tmp,
            "--search-type", "2",
            "--threads", str(attempt["threads"]),
            "-s", str(attempt["sens"]),
            "-e", str(min(args.aug_relaxed_evalue_cap, args.evalue * args.aug_relaxed_evalue_mult)),
            "--min-seq-id", "0.0",
            "--split-memory-limit", str(attempt["split"]),
            "-v", str(args.mmseqs_verbosity),
        ]

        proc = subprocess.run(
            search_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if proc.returncode != 0:
            details = _format_mmseqs_failure(proc)
            last_details = details
            resource_fail = _is_mmseqs_resource_failure(details)
            logger.warning(
                f"[{genome_name}] Augmented MMseqs search attempt {i}/{len(attempts)} failed "
                f"({attempt['label']}): {details}"
            )
            if not resource_fail:
                break
            continue

        # convertalis
        conv_cmd = [
            "mmseqs", "convertalis",
            query_db, target_db, result_db, aug_hits_file,
            "--format-output", fmt_output,
            "-v", str(args.mmseqs_verbosity),
        ]
        conv_proc = subprocess.run(
            conv_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        if conv_proc.returncode != 0:
            last_details = _format_mmseqs_failure(conv_proc)
            break

        if i > 1:
            logger.warning(
                f"[{genome_name}] Augmented MMseqs succeeded on retry {i}/{len(attempts)} "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']})."
            )
        return True, ""

    return False, last_details or "unknown augmented MMseqs failure"

def run_augmented_search(region_fasta: str, goi_queries: List[Dict[str, str]], 
                        genome_name: str, args, unique_id: str, threads: int) -> List[Dict[str, Any]]:
    """
    Run augmented search (MMseqs2 + tblastn + Smith-Waterman) for GOI queries.
    
    Uses both methods for maximum sensitivity:
    1. MMseqs2 with query fragments (fast, good for similar sequences)
    2. tblastn over the region (sensitive translated search)
    3. Smith-Waterman via parasail/ssearch36 (slower, better for divergent sequences)
    
    Args:
        region_fasta: Path to extracted region FASTA
        goi_queries: List of GOI query dicts with 'id' and 'seq'
        genome_name: Name of current genome
        args: Command line arguments
        unique_id: Unique ID for temp files
        threads: Number of threads to use
        
    Returns:
        List of hit dictionaries combining MMseqs2 and Smith-Waterman results
    """
    all_hits = []
    
    try:
        # Generate variants for each GOI query
        variants_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_variants.faa"
        query_fasta = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_goi_query.faa"
        all_variants = []
        
        # Write full query sequences for Smith-Waterman
        write_fasta([(q['id'], q['seq']) for q in goi_queries], query_fasta)
        goi_query_lengths = {
            q['id']: len(q['seq'])
            for q in goi_queries
            if q.get('id') and q.get('seq')
        }
        
        if not FRAGMENT_SUPPORT:
            print(f"[{genome_name}] Warning: fragment_query module not available, using full sequences only", flush=True)
            # Just use the original sequences
            all_variants = [(q['id'], q['seq']) for q in goi_queries]
        else:
            for query in goi_queries:
                # Generate fragments (halves, thirds, quarters)
                fragments = generate_fragments(query['seq'], query['id'], min_size=15)
                all_variants.extend([(f[0], f[1]) for f in fragments])
        
        write_fasta(all_variants, variants_fasta)
        variant_query_lengths = {
            vid: len(vseq)
            for vid, vseq in all_variants
            if vid and vseq
        }
        
        # Shared relaxed thresholds for augmented search.
        relaxed_evalue = min(args.aug_relaxed_evalue_cap, args.evalue * args.aug_relaxed_evalue_mult)
        relaxed_identity = max(args.aug_relaxed_identity_min, args.min_identity * args.aug_relaxed_identity_factor)
        relaxed_length = max(args.aug_relaxed_length_min, int(args.min_length // args.aug_relaxed_length_div))

        # ========== 1. MMseqs2 Search ==========
        aug_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_hits.m8"
        aug_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_aug_mmseqs"
        
        os.makedirs(aug_tmp_dir, exist_ok=True)
        
        mmseqs_ok, mmseqs_details = run_augmented_mmseqs_with_retries(
            variants_fasta=variants_fasta,
            region_fasta=region_fasta,
            aug_hits_file=aug_hits_file,
            aug_tmp_dir=aug_tmp_dir,
            args=args,
            threads=threads,
            genome_name=genome_name,
        )

        if mmseqs_ok:
            # Parse hits - use relaxed thresholds for augmented search
            relaxed_identity = max(args.aug_relaxed_identity_min, args.min_identity * args.aug_relaxed_identity_factor)
            relaxed_length = max(args.aug_relaxed_length_min, int(args.min_length // args.aug_relaxed_length_div))
            
            mmseqs_hits = parse_hits(
                aug_hits_file,
                relaxed_identity,
                relaxed_length,
                args.evalue * args.aug_relaxed_parse_evalue_mult,
                query_lengths=variant_query_lengths
            )
            if mmseqs_hits:
                print(f"[{genome_name}] MMseqs2 augmented search found {len(mmseqs_hits)} hits.", flush=True)
                all_hits.extend(mmseqs_hits)
        else:
            print(f"[{genome_name}] Augmented MMseqs failed after retries: {mmseqs_details}", flush=True)

        # ========== 2. tblastn Search ==========
        tblastn_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_tblastn_hits.m8"
        blast_db_prefix = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_tblastn_db"
        try:
            makedb_result = subprocess.run(
                ["makeblastdb", "-in", region_fasta, "-dbtype", "nucl", "-out", blast_db_prefix],
                check=True,
                capture_output=True, text=True
            )
            tblastn_result = subprocess.run(
                [
                    "tblastn",
                    "-query", query_fasta,
                    "-db", blast_db_prefix,
                    "-evalue", str(relaxed_evalue),
                    "-seg", "no",
                    "-max_target_seqs", "5000",
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                    "-out", tblastn_hits_file
                ],
                check=True,
                capture_output=True, text=True
            )

            tblastn_hits = parse_hits(
                tblastn_hits_file,
                relaxed_identity,
                relaxed_length,
                args.evalue * args.aug_relaxed_parse_evalue_mult,
                query_lengths=goi_query_lengths
            )
            if tblastn_hits:
                print(f"[{genome_name}] tblastn augmented search found {len(tblastn_hits)} hits.", flush=True)
                for hit in tblastn_hits:
                    hit['method'] = 'tblastn'
                all_hits.extend(tblastn_hits)
        except FileNotFoundError:
            logger.warning(f"[{genome_name}] tblastn/makeblastdb not found on PATH. Install BLAST+ to enable tblastn search.")
        except subprocess.CalledProcessError as blast_err:
            stderr_snippet = (blast_err.stderr or '')[:300].strip()
            logger.warning(
                f"[{genome_name}] tblastn failed (exit {blast_err.returncode}), continuing with MMseqs/SW only."
                + (f" stderr: {stderr_snippet}" if stderr_snippet else "")
            )
        
        # ========== 3. Smith-Waterman Search ==========
        # Use Smith-Waterman for very divergent sequences (more sensitive than MMseqs2)
        if args.enable_smith_waterman:
            sw_hits_file = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_sw_hits.m8"
            
            try:
                # Use smith_waterman_search.py script
                sw_cmd = [
                    "python3", os.path.join(os.path.dirname(__file__), "smith_waterman_search.py"),
                    "--query", query_fasta,
                    "--target", region_fasta,
                    "--output", sw_hits_file,
                    "--min_score", str(int(round(args.sw_min_score))),
                    "--min_identity", str(args.sw_min_identity),
                    "--threads", str(threads),
                    "--method", str(args.sw_method)
                ]
                
                result = subprocess.run(sw_cmd, capture_output=True, text=True, timeout=args.sw_timeout_seconds)
                
                if result.returncode == 0 and os.path.exists(sw_hits_file):
                    # Parse Smith-Waterman hits (BLAST m8 format)
                    sw_hits = parse_hits(
                        sw_hits_file,
                        args.sw_min_identity,
                        10,    # min 10 aa (was 2 – way too permissive)
                        1.0,   # evalue <= 1.0 (was 20000 – effectively no filter)
                        query_lengths=goi_query_lengths
                    )
                    if sw_hits:
                        print(f"[{genome_name}] Smith-Waterman found {len(sw_hits)} additional hits.", flush=True)
                        # Mark hits as from Smith-Waterman
                        for hit in sw_hits:
                            hit['method'] = 'smith_waterman'
                        all_hits.extend(sw_hits)
                elif result.stderr:
                    print(f"[{genome_name}] Smith-Waterman warning: {result.stderr[:200]}", flush=True)
                    
            except subprocess.TimeoutExpired:
                print(f"[{genome_name}] Smith-Waterman timed out after {args.sw_timeout_seconds} seconds, using MMseqs2 only.", flush=True)
            except FileNotFoundError:
                print(f"[{genome_name}] Smith-Waterman script not found, using MMseqs2 only.", flush=True)
            except Exception as sw_err:
                print(f"[{genome_name}] Smith-Waterman failed: {sw_err}, using MMseqs2 only.", flush=True)
        
        # ========== 4. PLM Embedding Search ==========
        # Use ProtT5 protein language model to find ORFs that are functionally
        # similar to the GOI even when sequence identity is undetectable (<15%).
        if getattr(args, 'enable_plm_search', False) and PLM_IMPORT_OK and check_plm_available():
            try:
                goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
                if os.path.exists(goi_emb_path):
                    goi_embeddings = load_embeddings(goi_emb_path)
                    if goi_embeddings:
                        plm_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_plm"
                        os.makedirs(plm_tmp_dir, exist_ok=True)

                        # Short queries produce noisy mean-pooled embeddings —
                        # raise the similarity threshold to suppress false hits.
                        max_goi_len = max(goi_query_lengths.values()) if goi_query_lengths else 0
                        plm_thresh = getattr(args, 'plm_similarity_threshold', 0.5)
                        if 0 < max_goi_len < 100:
                            plm_thresh = max(plm_thresh, 0.75)

                        plm_hits = plm_search_region(
                            goi_embeddings=goi_embeddings,
                            region_fasta=region_fasta,
                            output_dir=plm_tmp_dir,
                            similarity_threshold=plm_thresh,
                            device=getattr(args, 'plm_device', 'cpu'),
                            predictor='augustus',  # Enforce Augustus for eukaryotic ORF prediction; no Prodigal fallback (loses multi-exon genes)
                            augustus_species=getattr(args, 'augustus_species', 'fly'),
                        )

                        if plm_hits:
                            print(
                                f"[{genome_name}] PLM embedding search found "
                                f"{len(plm_hits)} candidate ORF(s).",
                                flush=True,
                            )
                            all_hits.extend(plm_hits)

                        if os.path.exists(plm_tmp_dir):
                            shutil.rmtree(plm_tmp_dir, ignore_errors=True)
            except Exception as plm_err:
                logger.debug(f"[{genome_name}] PLM search failed: {plm_err}")

        # ========== 5. Foldseek Structural Search ==========
        # Use ESMFold to predict ORF structures and Foldseek 3Di alphabet to
        # find structural homologs that have diverged beyond all sequence methods.
        # Skip structural search for very short queries — ESMFold produces
        # degenerate structures below ~50 aa and TM-score is meaningless for
        # single short helices (e.g. melittin, defensins).
        _max_goi_len = max(goi_query_lengths.values()) if goi_query_lengths else 0
        _struct_too_short = 0 < _max_goi_len < 50
        if (getattr(args, 'enable_structural_search', False)
                and STRUCTURAL_IMPORT_OK and check_structural_search_available()
                and _struct_too_short):
            print(
                f"[{genome_name}] Skipping structural search: "
                f"longest GOI query is {_max_goi_len} aa (<50 aa threshold).",
                flush=True,
            )

        if (getattr(args, 'enable_structural_search', False)
                and STRUCTURAL_IMPORT_OK and check_structural_search_available()
                and not _struct_too_short):
            try:
                goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
                if os.path.exists(goi_struct_index):
                    goi_structures = load_structure_index(goi_struct_index)
                    if goi_structures:
                        struct_tmp_dir = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_struct"
                        os.makedirs(struct_tmp_dir, exist_ok=True)

                        struct_hits = structural_search_region(
                            goi_structures=goi_structures,
                            region_fasta=region_fasta,
                            output_dir=struct_tmp_dir,
                            tm_threshold=getattr(
                                args, 'structural_tm_threshold', 0.3
                            ),
                            device=getattr(args, 'structural_device', 'cpu'),
                            max_length=getattr(args, 'structural_max_length', 700),
                            threads=getattr(args, 'threads', 1),
                            predictor='augustus',  # Enforce Augustus for eukaryotic ORF prediction; no Prodigal fallback (loses multi-exon genes)
                            augustus_species=getattr(args, 'augustus_species', 'fly'),
                        )

                        if struct_hits:
                            print(
                                f"[{genome_name}] Foldseek structural search found "
                                f"{len(struct_hits)} candidate ORF(s).",
                                flush=True,
                            )
                            all_hits.extend(struct_hits)

                        if os.path.exists(struct_tmp_dir):
                            shutil.rmtree(struct_tmp_dir, ignore_errors=True)
            except Exception as struct_err:
                logger.debug(f"[{genome_name}] Structural search failed: {struct_err}")

        # Clean up temp files
        files_to_remove = [variants_fasta, query_fasta, aug_hits_file]
        if 'sw_hits_file' in locals():
             files_to_remove.append(sw_hits_file)
        files_to_remove.append(tblastn_hits_file)
            
        for f in files_to_remove:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(aug_tmp_dir):
            shutil.rmtree(aug_tmp_dir, ignore_errors=True)
        for f in glob.glob(f"{blast_db_prefix}*"):
            if os.path.exists(f):
                os.remove(f)
        
        # Deduplicate hits by position (Greedy Spatial Clustering)
        if all_hits:
            # Sort by bitscore desc, then evalue asc
            # Ensure bits is float
            for h in all_hits:
                if 'bits' not in h: h['bits'] = 0.0
                h['bits'] = float(h['bits'])
            
            all_hits.sort(key=lambda x: (-x['bits'], x['evalue']))
            
            kept_hits = []
            for hit in all_hits:
                # Check overlap with kept hits
                is_overlapped = False
                h_start = hit['start']
                h_end = hit['end']
                h_len = h_end - h_start
                
                for kept in kept_hits:
                    if hit['chrom'] != kept['chrom']: continue
                    
                    k_start = kept['start']
                    k_end = kept['end']
                    
                    # Calculate overlap
                    o_start = max(h_start, k_start)
                    o_end = min(h_end, k_end)
                    overlap = max(0, o_end - o_start)
                    
                    # If overlap covers > 50% of THIS hit, drop it (we already have a better one)
                    if overlap > 0.5 * h_len:
                        is_overlapped = True
                        break
                
                if not is_overlapped:
                    kept_hits.append(hit)
            
            all_hits = kept_hits
            print(f"[{genome_name}] Combined augmented search: {len(all_hits)} unique hits (spatially filtered).", flush=True)
        
        return all_hits
        
    except Exception as e:
        import traceback
        logger.warning(
            f"[{genome_name}] Augmented search failed: {e}\n"
            f"  Traceback: {traceback.format_exc().strip().splitlines()[-3:]}"
        )
        return []

def batch_rbh_check(
    candidates,
    home_db,
    cand_map,
    threads=1,
    evalue=1e-5,
    min_coverage=0.5,
    min_identity=25.0
):
    """
    Perform Reciprocal Best Hit check for multiple candidates at once.
    
    Enhanced validation:
    1. RBH to home genome (traditional)
    2. Coverage check: alignment must cover >50% of both query and target
    3. Identity must be reasonable (>25%)
    
    candidates: list of dicts with 'id' and 'seq' keys (or BioPython SeqRecord)
    """
    if not candidates:
        return set()
    
    unique_id = uuid.uuid4().hex
    query_fasta = f"batch_candidates_{unique_id}.fasta"
    rbh_out = f"batch_rbh_{unique_id}.m8"
    tmp_subdir = f"tmp_rbh_batch_{unique_id}"
    
    valid_ids = set()

    try:
        def _ids_match(parent_name: str, target_name: str) -> bool:
            parent_base = extract_base_gene_id(parent_name).strip()
            target_base = extract_base_gene_id(target_name).strip()
            return (
                parent_base == target_base or
                parent_name == target_name or
                target_name == parent_name or
                parent_base in target_base or
                target_base in parent_base
            )

        def _best_ungapped_metrics(query_seq: str, target_seq: str) -> Tuple[float, float, float]:
            """
            Best ungapped overlap identity between two protein sequences.
            Returns (identity_pct, qcov_pct, tcov_pct).
            """
            q = query_seq or ""
            t = target_seq or ""
            if not q or not t:
                return 0.0, 0.0, 0.0

            best_identity = 0.0
            best_qcov = 0.0
            best_tcov = 0.0
            best_overlap = 0

            # offset < 0 means target starts before query
            for offset in range(-len(t) + 1, len(q)):
                q_start = max(0, offset)
                t_start = max(0, -offset)
                overlap = min(len(q) - q_start, len(t) - t_start)
                if overlap <= 0:
                    continue

                matches = 0
                for idx in range(overlap):
                    if q[q_start + idx] == t[t_start + idx]:
                        matches += 1

                identity = (matches / overlap) * 100.0
                qcov = (overlap / len(q)) * 100.0
                tcov = (overlap / len(t)) * 100.0

                if (
                    identity > best_identity or
                    (identity == best_identity and overlap > best_overlap)
                ):
                    best_identity = identity
                    best_qcov = qcov
                    best_tcov = tcov
                    best_overlap = overlap

            return best_identity, best_qcov, best_tcov

        def _fallback_rbh_local(records_local, db_for_convert, mapping) -> set:
            """
            Fallback RBH path for short proteins where mmseqs prefilter cannot run.
            """
            fallback_valid = set()
            tmp_home_faa = f"home_export_{uuid.uuid4().hex}.faa"
            try:
                subprocess.run(
                    ["mmseqs", "convert2fasta", db_for_convert, tmp_home_faa],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                home_records = [(hid, hseq) for _, hid, hseq in parse_fasta(tmp_home_faa)]
                if not home_records:
                    return fallback_valid

                for cand_id, cand_seq in records_local:
                    if cand_id not in mapping:
                        continue

                    best_target = None
                    best_identity = -1.0
                    best_qcov = 0.0
                    best_tcov = 0.0

                    for target_id, target_seq in home_records:
                        ident, qcov, tcov = _best_ungapped_metrics(cand_seq, target_seq)
                        if (
                            ident > best_identity or
                            (ident == best_identity and (qcov + tcov) > (best_qcov + best_tcov))
                        ):
                            best_identity = ident
                            best_target = target_id
                            best_qcov = qcov
                            best_tcov = tcov

                    if not best_target:
                        continue

                    parent = mapping[cand_id]
                    coverage_ok = (
                        best_qcov >= min_coverage * 100 and
                        best_tcov >= min_coverage * 100
                    )
                    identity_ok = best_identity >= min_identity
                    if _ids_match(parent, best_target) and coverage_ok and identity_ok:
                        fallback_valid.add(cand_id)

            finally:
                if os.path.exists(tmp_home_faa):
                    os.remove(tmp_home_faa)

            return fallback_valid

        def _candidate_to_pair(candidate):
            # Dict form used by pipeline internals
            if isinstance(candidate, dict):
                cid = candidate.get('id')
                cseq = candidate.get('seq')
                return cid, str(cseq) if cseq is not None else None

            # SeqRecord-like form used by unit tests
            cid = getattr(candidate, 'id', None)
            cseq_obj = getattr(candidate, 'seq', None)
            cseq = str(cseq_obj) if cseq_obj is not None else None
            return cid, cseq

        # Write FASTA using our utility
        records = []
        for c in candidates:
            cid, cseq = _candidate_to_pair(c)
            if cid and cseq:
                records.append((cid, cseq))

        if not records:
            return set()

        write_fasta(records, query_fasta)
            
        db_path = home_db
        if os.path.isdir(home_db):
            db_path = os.path.join(home_db, "db")
            
        cmd = [
            "mmseqs", "easy-search",
            query_fasta, db_path, rbh_out, tmp_subdir,
            "-e", str(evalue),
            "--format-output", "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen",
            "--max-seqs", "1", # Top hit only
            "--split-memory-limit", str(args.mmseqs_split_memory_limit),
            "-v", str(args.mmseqs_verbosity),
            "--threads", str(threads)
        ]

        mmseqs_ok = True
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            mmseqs_ok = False
            logger.warning(
                "RBH mmseqs easy-search failed; using local fallback "
                "(expected for very short peptides)."
            )
            valid_ids.update(_fallback_rbh_local(records, db_path, cand_map))

        if mmseqs_ok and os.path.exists(rbh_out):
            with open(rbh_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 9: continue
                    
                    cand_id = parts[0]
                    target_id = parts[1]
                    pident = float(parts[2])
                    qcov = float(parts[3]) if len(parts) > 3 else 100
                    tcov = float(parts[4]) if len(parts) > 4 else 100

                    # MMseqs coverage can be reported as fraction [0,1] or percent [0,100]
                    if qcov <= 1.0:
                        qcov *= 100.0
                    if tcov <= 1.0:
                        tcov *= 100.0
                    
                    # Enhanced validation
                    if cand_id not in cand_map:
                        continue
                    
                    parent = cand_map[cand_id]
                    
                    # Check 1: ID matching (exact or close)
                    ids_match = _ids_match(parent, target_id)
                    
                    # Check 2: Coverage (both query and target must be well-covered)
                    coverage_ok = (qcov >= min_coverage * 100 and 
                                  tcov >= min_coverage * 100)
                    
                    # Check 3: Identity must be reasonable
                    identity_ok = pident >= min_identity
                    
                    if ids_match and coverage_ok and identity_ok:
                        valid_ids.add(cand_id)
                    elif ids_match and not coverage_ok:
                        print(f"RBH: {cand_id} matches {target_id} but low coverage "
                              f"(qcov={qcov:.0f}%, tcov={tcov:.0f}%). Likely fragment/paralog.",
                              file=sys.stderr)
                    elif ids_match and not identity_ok:
                        logger.debug(f"RBH: {cand_id} matches {target_id} but very low identity "
                              f"({pident:.1f}%). Possible pseudogene.")
                              
    except Exception as e:
        logger.error(f"RBH check failed: {e}")
        return set()
    finally:
        # Cleanup
        if os.path.exists(query_fasta): os.remove(query_fasta)
        if os.path.exists(rbh_out): os.remove(rbh_out)
        if os.path.exists(tmp_subdir):
             shutil.rmtree(tmp_subdir, ignore_errors=True)
                          
    return valid_ids

def process_region_block(block_idx, block, hits, genome_seqs, db_sequences, genome_name, args, unique_id,
                         threads_per_job, native_annot_index=None):
    """
    Process a single synteny block: extract region, identify queries, run augmented search, and annotate.
    Returns: (list_of_new_genes, list_of_gff_lines)
    """
    chrom = block['chrom']
    if chrom not in genome_seqs:
        return [], []

    slen = len(genome_seqs[chrom])
    padding = calculate_adaptive_padding(
        hits,
        block,
        default=args.region_padding,
        min_pad=args.padding_min,
        max_pad=args.padding_max
    )

    w_start = max(0, block['start'] - padding)
    w_end = min(slen, block['end'] + padding)
    subseq = genome_seqs[chrom][w_start:w_end]
    # Keep final GOI calls anchored near the flanking-defined block.
    core_start = max(0, block['start'] - args.gap_search_window)
    core_end = min(slen, block['end'] + args.gap_search_window)

    temp_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_block{block_idx}.fasta"
    query_mini_fa = f"{args.output_dir}/tmp_{unique_id}_{genome_name}_block{block_idx}_query.faa"
    write_fasta([("region_seq", subseq)], temp_fa)

    relevant_hits = [
        h for h in hits
        if h['chrom'] == chrom and h['end'] > w_start and h['start'] < w_end
    ]
    block_flanking_support = len({
        extract_base_gene_id(h.get('query', ''))
        for h in relevant_hits
        if h.get('query') and not is_goi_query_id(h.get('query', ''))
    })
    unique_queries = set(extract_base_gene_id(h['query']) for h in relevant_hits)

    found_queries = []

    # Always include GOI-derived queries so blocks remain anchored to the target family.
    goi_ids = [k for k in db_sequences.keys() if is_goi_query_id(k)]
    full_length_goi = {gid for gid in goi_ids if '|frag_' not in gid}
    if full_length_goi:
        unique_queries.update(full_length_goi)
    else:
        unique_queries.update(goi_ids)

    for query_id in unique_queries:
        if query_id in db_sequences:
            found_queries.append(db_sequences[query_id])
            continue
        base = extract_base_gene_id(query_id)
        if base in db_sequences:
            found_queries.append(db_sequences[base])

    if not found_queries:
        if os.path.exists(temp_fa):
            os.remove(temp_fa)
        return [], []

    # Region annotation should be driven by full-length GOI queries.
    # Flanking genes are used to define candidate regions, not to produce GOI models.
    full_goi_queries = [q for q in found_queries if is_full_length_goi_query_id(q.get('id', ''))]
    goi_queries = full_goi_queries if full_goi_queries else [q for q in found_queries if is_goi_query_id(q.get('id', ''))]
    if not goi_queries:
        if os.path.exists(temp_fa):
            os.remove(temp_fa)
        return [], []

    # Collapse duplicated GOI entries by parent ID and keep the longest representative.
    goi_query_by_parent = {}
    for q in goi_queries:
        parent = extract_base_gene_id(q.get('id', ''))
        if not parent:
            continue
        if parent not in goi_query_by_parent or len(q.get('seq', '')) > len(goi_query_by_parent[parent].get('seq', '')):
            goi_query_by_parent[parent] = q
    goi_queries = list(goi_query_by_parent.values())

    write_fasta([(q['id'], q['seq']) for q in goi_queries], query_mini_fa)

    annotated_records_raw = []
    valid_gff_lines = []
    raw_candidates = []
    flanking_candidates = []
    clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')

    def _mRNA_attrs(base_attrs: Dict[str, Any], g_start: int, g_end: int, g_strand: str) -> str:
        native = lookup_native_annotation(
            native_annot_index or {},
            chrom=chrom,
            start=g_start,
            end=g_end,
            strand=g_strand,
        )
        return _compose_gff_attrs(base_attrs, native)

    try:
        augmented_hits_mmseqs = run_augmented_search(
            temp_fa,
            goi_queries,
            f"{genome_name}_b{block_idx}",
            args,
            f"{unique_id}_b{block_idx}",
            threads_per_job
        )

        all_search_hits = list(augmented_hits_mmseqs) if augmented_hits_mmseqs else []

        # Include original block-supporting hits (converted to local coords).
        for hit in relevant_hits:
            if not is_goi_query_id(hit.get('query', '')):
                continue
            region_hit = dict(hit)
            local_start = max(0, hit['start'] - w_start)
            local_end = min(w_end - w_start, hit['end'] - w_start)
            if local_end <= local_start:
                continue
            region_hit['start'] = local_start
            region_hit['end'] = local_end
            all_search_hits.append(region_hit)

        if all_search_hits:
            hits_by_gene = defaultdict(list)
            for hit in all_search_hits:
                local_start = max(0, int(hit.get('start', 0)))
                local_end = min(len(subseq), int(hit.get('end', 0)))
                if local_end <= local_start:
                    continue

                global_start = w_start + local_start
                global_end = w_start + local_end
                if global_end <= core_start or global_start >= core_end:
                    continue

                query_id = hit.get('query', '')
                if not query_id:
                    continue

                parent = query_id.split('|frag_')[0] if '|frag_' in query_id else query_id
                parent = extract_base_gene_id(parent)

                # Calculate offset for fragment or exon queries so qstart/qend
                # are expressed in full-protein coordinates.
                q_offset = 0
                if '|frag_' in query_id and FRAGMENT_SUPPORT:
                    try:
                        frag_info = parse_fragment_id(query_id)
                        q_offset = frag_info['start'] - 1
                    except Exception:
                        q_offset = 0
                elif '|exon_' in query_id:
                    # Build exon offset map from DB sequences on first encounter.
                    if not hasattr(run_augmented_search, '_goi_exon_offsets'):
                        run_augmented_search._goi_exon_offsets = {}
                    if parent not in run_augmented_search._goi_exon_offsets:
                        exon_seqs = []
                        for rid, rec in db_sequences.items():
                            if extract_base_gene_id(rid) == parent:
                                m = re.search(r'\|exon_(\d+)$', rid)
                                if m and rec.get('seq'):
                                    exon_seqs.append((int(m.group(1)), len(rec['seq'])))
                        exon_seqs.sort()
                        offmap = {}
                        off = 0
                        for enum, elen in exon_seqs:
                            offmap[enum] = off
                            off += elen
                        run_augmented_search._goi_exon_offsets[parent] = offmap
                    offmap = run_augmented_search._goi_exon_offsets.get(parent, {})
                    m = re.search(r'\|exon_(\d+)$', query_id)
                    if m:
                        enum = int(m.group(1))
                        q_offset = offmap.get(enum, 0)

                hits_by_gene[parent].append({
                    'qstart': hit.get('qstart', 1) + q_offset,
                    'qend': hit.get('qend', 100) + q_offset,
                    'gstart': local_start,
                    'gend': local_end,
                    'evalue': hit.get('evalue', 1),
                    'pident': hit.get('pident', 0),
                    'alnlen': hit.get('alnlen', 0),
                    'bits': hit.get('bits', 0),
                    'strand': hit.get('strand', '+'),
                    'chrom': chrom
                })

            exon_parents = set()
            for q in found_queries:
                qid = q.get('id', '')
                if '|exon_' in qid:
                    exon_parents.add(extract_base_gene_id(qid))

            locus_gap = max(5000, min(100000, args.max_intron * 2))

            for parent_id, parent_hits in hits_by_gene.items():
                is_goi_parent = is_goi_query_id(parent_id)
                if not is_goi_parent:
                    continue

                parent_query = goi_query_by_parent.get(parent_id)
                if not parent_query:
                    continue
                parent_query_seq = parent_query.get('seq', '')
                if not parent_query_seq:
                    continue

                parent_loci = split_hits_into_loci(parent_hits, max_gap=locus_gap)
                print(f"[DEBUG PARENT] parent_id={parent_id} is_goi={is_goi_parent} num_loci={len(parent_loci)} exon_parents={exon_parents}", flush=True)
                if not parent_loci:
                    continue

                def _locus_rank(locus_hits):
                    qmin = min(h.get('qstart', 1) for h in locus_hits)
                    qmax = max(h.get('qend', 1) for h in locus_hits)
                    qcov = ((qmax - qmin + 1) / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                    best_bits = max(h.get('bits', 0) for h in locus_hits)
                    return (qcov, best_bits, len(locus_hits))

                # Keep only the strongest loci per parent to avoid noisy paralog floods.
                parent_loci = sorted(parent_loci, key=_locus_rank, reverse=True)[:2]

                for locus_idx, gene_hits in enumerate(parent_loci, start=1):
                    work_hits = gene_hits
                    if parent_id in exon_parents:
                        print(f"[DEBUG EXON] {parent_id} starting with {len(work_hits)} hits", flush=True)
                        work_hits = filter_exon_hits(
                            work_hits,
                            len(parent_query_seq),
                            args.min_exon_query_cov,
                            args.min_exon_alnlen
                        )
                        print(f"[DEBUG EXON] {parent_id} after filter: {len(work_hits)} hits", flush=True)
                        if not work_hits:
                            print(f"[DEBUG EXON] {parent_id} dropped entirely by filter_exon_hits!", flush=True)
                            continue

                    exons = []
                    try:
                        is_tandem = False
                        tandem_copies = []
                        if is_goi_parent:
                            from annotate_goi_exons import detect_tandem_duplications
                            is_tandem, tandem_copies = detect_tandem_duplications(
                                work_hits, parent_query_seq, subseq, chrom
                            )

                        if is_tandem and tandem_copies:
                            exons = tandem_copies
                        else:
                            with maybe_quiet_streams(args.quiet_subtools):
                                exons, _ = annotate_exons_from_hit_list(
                                    work_hits,
                                    parent_query_seq,
                                    subseq,
                                    chrom,
                                    search_missing=True,
                                    gap_min_size=args.gap_min_size,
                                    gap_search_window=args.gap_search_window,
                                    gap_evalue=args.gap_evalue,
                                    gap_min_identity=args.gap_min_identity,
                                    gap_min_alnlen=args.gap_min_alnlen,
                                    gap_max_hits=args.gap_max_hits,
                                    exon_query_mode=(parent_id in exon_parents),
                                    min_exon_query_cov=args.min_exon_query_cov,
                                    min_exon_alnlen=args.min_exon_alnlen,
                                    sensitive=is_goi_parent
                                )
                    except Exception as e:
                        import traceback
                        print(f"[ERROR] miniprot annotation failed for {parent_id} at {chrom}:{w_start}-{w_end}: {e}")
                        traceback.print_exc()
                        exons = []

                    if not exons:
                        if work_hits:
                            strand_votes = {'+': 0.0, '-': 0.0}
                            for h in work_hits:
                                strand_votes[h.get('strand', '+')] += float(h.get('bits', 1))
                            strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                            ordered_hits = sorted(
                                [h for h in work_hits if h.get('strand', '+') == strand],
                                key=lambda h: h.get('gstart', 0)
                            )
                            if ordered_hits:
                                ordered_hits = _longest_monotonic_query_chain(ordered_hits, strand)
                                if ordered_hits:
                                    qmin = min(min(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                    qmax = max(max(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                    query_len = max(1, len(parent_query_seq))
                                    qcov = (qmax - qmin + 1) / query_len
                                    aln_total = sum(max(0, int(h.get('alnlen', 0))) for h in ordered_hits)
                                    best_bits = max(float(h.get('bits', 0)) for h in ordered_hits)
                                    print(f"[DEBUG FALLBACK] GOI={parent_id} qcov={qcov:.2f} aln_total={aln_total} bits={best_bits:.1f} len(hits)={len(ordered_hits)}", flush=True)

                                    valid_fallback = True
                                    if len(ordered_hits) == 1:
                                        # Less strict coverage constraint for GOI fragments matching strongly
                                        if qcov < 0.25 or int(ordered_hits[0].get('alnlen', 0)) < 25:
                                            valid_fallback = False
                                    else:
                                        if qcov < 0.25 and aln_total < 35 and best_bits < 60.0:
                                            valid_fallback = False

                                    # Span sanity check: reject fallback hits whose
                                    # genomic span is wildly disproportionate to the
                                    # query protein length. For short peptides
                                    # (<150 aa) the intron:exon ratio is dominated
                                    # by introns, so a flat 30x multiplier rejects
                                    # valid multi-exon loci (e.g. melittin at 70 aa
                                    # can span tens of kb in bee genomes).
                                    if valid_fallback:
                                        # Short queries (<150 aa, e.g. melittin/toxin peptides):
                                        #   genomic span is dominated by introns; allow up to
                                        #   600x the CDS length so real multi-exon loci pass.
                                        # Long queries (>=150 aa, e.g. TP53/TP63/TP73):
                                        #   raised from 30x -> 50x after TP53 runs showed real
                                        #   TP63/TP73 loci spanning >100 kb being rejected.
                                        if query_len < 150:
                                            max_span_nt = max(100_000, query_len * 3 * 200)
                                        else:
                                            max_span_nt = max(5000, query_len * 3 * 50)
                                        gspan_min = min(h.get('gstart', 0) for h in ordered_hits)
                                        gspan_max = max(h.get('gend', 0) for h in ordered_hits)
                                        actual_span = gspan_max - gspan_min
                                        if actual_span > max_span_nt:
                                            print(
                                                f"[DEBUG FALLBACK] REJECTED span={actual_span} "
                                                f"> max_span={max_span_nt} for {parent_id}",
                                                flush=True,
                                            )
                                            valid_fallback = False

                                    print(f"[DEBUG FALLBACK] valid_fallback={valid_fallback}", flush=True)

                                    if valid_fallback:
                                        coding_frags = []
                                        cds_intervals = []
                                        for h in ordered_hits:
                                            hs = h.get('gstart', 0)
                                            he = h.get('gend', 0)
                                            if he <= hs: continue
                                            dna = subseq[hs:he]
                                            if strand == '-': dna = reverse_complement(dna)
                                            dna = dna[:len(dna) - (len(dna) % 3)]
                                            if len(dna) < 9: continue
                                            aa = translate(dna).replace('*', '')
                                            if not aa: continue
                                            coding_frags.append(aa)
                                            cds_intervals.append((hs, he))

                                        if coding_frags and cds_intervals:
                                            if strand == '-': coding_frags = list(reversed(coding_frags))
                                            flank_protein = ''.join(coding_frags)
                                            global_start = w_start + min(s for s, _ in cds_intervals) + 1
                                            global_end = w_start + max(e for _, e in cds_intervals)
                                            avg_pident = sum(h.get('pident', 0) for h in ordered_hits) / len(ordered_hits)
                                            
                                            copy_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_fallback"
                                            
                                            cand_gff = [
                                                (
                                                    f"{chrom}\tfallback_hits\tmRNA\t{global_start}\t{global_end}\t"
                                                    f"{avg_pident:.1f}\t{strand}\t.\t"
                                                    f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'fallback_hit_span'}, evidence_type='fallback_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=qcov, flanking_support=block_flanking_support), global_start, global_end, strand)}"
                                                )
                                            ]
                                            for eidx, (hs, he) in enumerate(cds_intervals, 1):
                                                exon_gs = w_start + hs + 1
                                                exon_ge = w_start + he
                                                cand_gff.append(
                                                    f"{chrom}\tfallback_hits\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\tID={copy_id}_CDS{eidx};Parent={copy_id}"
                                                )
                                            
                                            raw_candidates.append({
                                                'start': global_start,
                                                'end': global_end,
                                                'score': avg_pident,
                                                'record': {
                                                    'id': copy_id,
                                                    'seq': flank_protein,
                                                    'description': f"coords:{global_start}-{global_end} parent:{parent_id} type:fallback identity:{avg_pident:.1f}"
                                                },
                                                'gff': cand_gff
                                            })

                    if exons:
                        is_tandem_result = any(e.get('id', '').startswith('GOI_copy_') for e in exons)
                        if is_tandem_result:
                            for copy in exons:
                                global_start = w_start + copy['gstart'] + 1
                                global_end = w_start + copy['gend']
                                strand = copy.get('strand', '+')
                                copy_id = f"{copy['id']}|{clean_gname}_b{block_idx}_l{locus_idx}"
                                raw_candidates.append({
                                    'start': global_start,
                                    'end': global_end,
                                    'score': copy.get('pident', 0),
                                    'record': {
                                        'id': copy_id,
                                        'seq': copy['seq'],
                                        'description': (
                                            f"coords:{global_start}-{global_end} "
                                            f"parent:{parent_id} tandem_copy "
                                            f"identity:{copy.get('pident', 0):.1f}"
                                        )
                                    },
                                    'gff': [
                                        f"{chrom}\ttandem_copy\tgene\t{global_start}\t{global_end}\t"
                                        f"{copy.get('pident', 0):.1f}\t{strand}\t.\t"
                                        f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': copy['id'], 'SynVoy_Parent': parent_id, 'Type': 'tandem_copy'}, evidence_type='tandem_copy', identity=copy.get('pident', 0), exon_count=1, query_cov=None, flanking_support=block_flanking_support), global_start, global_end, strand)}"
                                    ]
                                })
                        else:
                            exons.sort(key=lambda e: e.get('qstart', 0))
                            exon_protein = ''.join(e['seq'] for e in exons)
                            strand = exons[0].get('strand', '+')
                            avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                            model_qcov = None
                            if parent_query_seq:
                                qmin_model = min(e.get('qstart', 1) for e in exons)
                                qmax_model = max(e.get('qend', 1) for e in exons)
                                model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))

                            global_start = w_start + min(e['gstart'] for e in exons) + 1
                            global_end = w_start + max(e['gend'] for e in exons)
                            new_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_exon_ann"

                            cand_gff = [
                                (
                                    f"{chrom}\texon_annotation\tmRNA\t{global_start}\t{global_end}\t"
                                    f"{avg_pident:.1f}\t{strand}\t.\t"
                                    f"{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='exon_annotation', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, flanking_support=block_flanking_support), global_start, global_end, strand)}"
                                )
                            ]
                            for eidx, exon in enumerate(exons, 1):
                                exon_gs = w_start + exon['gstart'] + 1
                                exon_ge = w_start + exon['gend']
                                attrs = f"ID={new_id}_CDS{eidx};Parent={new_id}"
                                if exon.get('splice_acceptor'):
                                    attrs += f";SpliceAcceptor={exon['splice_acceptor']}"
                                if exon.get('splice_donor'):
                                    attrs += f";SpliceDonor={exon['splice_donor']}"
                                if exon.get('has_start_codon'):
                                    attrs += ";StartCodon=ATG"
                                if exon.get('has_stop_codon'):
                                    attrs += ";StopCodon=yes"
                                cand_gff.append(
                                    f"{chrom}\texon_annotation\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\t{attrs}"
                                )

                            raw_candidates.append({
                                'start': global_start,
                                'end': global_end,
                                'score': avg_pident,
                                'record': {
                                    'id': new_id,
                                    'seq': exon_protein,
                                    'description': (
                                        f"coords:{global_start}-{global_end} "
                                        f"parent:{parent_id} exons:{len(exons)} "
                                        f"identity:{avg_pident:.1f}"
                                    )
                                },
                                'gff': cand_gff
                            })

                            # Rescue only for GOI-derived parents to prevent fallback explosions.
                            if is_goi_parent:
                                covered_intervals = [(e['gstart'], e['gend']) for e in exons]
                                sorted_raw = sorted(work_hits, key=lambda h: h.get('gstart', h.get('start', 0)))
                                extra_count = 0
                                for hit in sorted_raw:
                                    h_start = hit.get('gstart', hit.get('start', 0))
                                    h_end = hit.get('gend', hit.get('end', 0))
                                    h_len = h_end - h_start
                                    if h_len < 1:
                                        continue

                                    is_covered = False
                                    for es, ee in covered_intervals:
                                        overlap = max(0, min(h_end, ee) - max(h_start, es))
                                        if overlap > 0.5 * h_len:
                                            is_covered = True
                                            break
                                    if is_covered:
                                        continue

                                    extra_count += 1
                                    region_dna = subseq[h_start:h_end]
                                    if hit.get('strand', '+') == '-':
                                        region_dna = reverse_complement(region_dna)
                                    region_dna = region_dna[:len(region_dna) - len(region_dna) % 3]
                                    if len(region_dna) < 6:
                                        continue

                                    hit_prot = translate(region_dna).replace('*', '')
                                    if not hit_prot:
                                        continue

                                    gs = w_start + h_start + 1
                                    ge = w_start + h_end
                                    hit_qspan = abs(hit.get('qend', 0) - hit.get('qstart', 0)) + 1
                                    hit_qcov = (hit_qspan / len(parent_query_seq)) if parent_query_seq else None
                                    raw_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_extra{extra_count}"
                                    raw_candidates.append({
                                        'start': gs,
                                        'end': ge,
                                        'score': hit.get('pident', 0),
                                        'record': {
                                            'id': raw_id,
                                            'seq': hit_prot,
                                            'description': f"coords:{gs}-{ge} parent:{parent_id} type:rescued_exon"
                                        },
                                        'gff': [
                                            f"{chrom}\trescued_exon\tmRNA\t{gs}\t{ge}\t{hit.get('pident',0):.1f}\t{hit.get('strand','+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': raw_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='rescued_exon', identity=hit.get('pident', 0), exon_count=1, query_cov=hit_qcov, flanking_support=block_flanking_support), gs, ge, hit.get('strand','+'))}",
                                            f"{chrom}\trescued_exon\tCDS\t{gs}\t{ge}\t.\t{hit.get('strand','+')}\t0\tID={raw_id}_CDS1;Parent={raw_id}"
                                        ]
                                    })

                    # Structure inference failed: keep only a very strict GOI fallback set.
                    elif is_goi_parent:
                        locus_qmin = min(h.get('qstart', 1) for h in work_hits)
                        locus_qmax = max(h.get('qend', 1) for h in work_hits)
                        locus_qcov = ((locus_qmax - locus_qmin + 1) / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                        best_bits = max(h.get('bits', 0) for h in work_hits) if work_hits else 0.0

                        # If miniprot is available, only allow fallback on very strong evidence.
                        if MINIPROT_AVAILABLE and (len(work_hits) < 2 or locus_qcov < 0.75 or best_bits < 60):
                            continue

                        # Fallback should be conservative; otherwise short/noisy hits explode.
                        min_fallback_identity = max(args.gap_min_identity, 90.0)
                        min_fallback_alnlen = max(
                            args.min_exon_alnlen,
                            int(max(15, len(parent_query_seq) * 0.6))
                        )
                        sorted_hits = sorted(
                            work_hits,
                            key=lambda h: (h.get('pident', 0), h.get('alnlen', 0), h.get('bits', 0)),
                            reverse=True
                        )[:max(1, args.gap_max_hits)]

                        accepted_fallback = 0
                        for h_idx, hit in enumerate(sorted_hits, 1):
                            if hit.get('pident', 0) < min_fallback_identity:
                                continue
                            if hit.get('alnlen', 0) < min_fallback_alnlen:
                                continue
                            qspan = abs(hit.get('qend', 0) - hit.get('qstart', 0)) + 1
                            qcov = (qspan / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                            if qcov < 0.6:
                                continue

                            g_s, g_e = hit['gstart'], hit['gend']
                            region_dna = subseq[g_s:g_e]
                            if hit.get('strand', '+') == '-':
                                region_dna = reverse_complement(region_dna)
                            region_dna = region_dna[:len(region_dna) - len(region_dna) % 3]
                            if len(region_dna) < 9:
                                continue

                            hit_protein = translate(region_dna).replace('*', '')
                            if not hit_protein:
                                continue

                            nt_s = w_start + g_s + 1
                            nt_e = w_start + g_e
                            qcov = (qspan / len(parent_query_seq)) if len(parent_query_seq) > 0 else 0.0
                            new_id = f"{parent_id}|{clean_gname}_b{block_idx}_l{locus_idx}_fallback_{h_idx}"
                            raw_candidates.append({
                                'start': nt_s,
                                'end': nt_e,
                                'score': hit.get('pident', 0),
                                'record': {
                                    'id': new_id,
                                    'seq': hit_protein,
                                    'description': f"coords:{nt_s}-{nt_e} parent:{parent_id} identity:{hit.get('pident', 0):.1f}"
                                },
                                'gff': [
                                    f"{chrom}\traw_hit\tmRNA\t{nt_s}\t{nt_e}\t{hit.get('pident', 0):.1f}\t{hit.get('strand', '+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id}, evidence_type='raw_hit', identity=hit.get('pident', 0), exon_count=1, query_cov=qcov, flanking_support=block_flanking_support), nt_s, nt_e, hit.get('strand', '+'))}",
                                    f"{chrom}\traw_hit\tCDS\t{nt_s}\t{nt_e}\t.\t{hit.get('strand', '+')}\t0\tID={new_id}_CDS1;Parent={new_id}"
                                ]
                            })
                            accepted_fallback += 1
                            if accepted_fallback >= 1:
                                break

        # Flanking gene annotation pass (kept separate from GOI logic).
        flanking_parent_ids = {
            extract_base_gene_id(h.get('query', ''))
            for h in relevant_hits
            if h.get('query') and not is_goi_query_id(h.get('query', ''))
        }
        flanking_query_by_parent = build_flanking_query_by_parent(db_sequences, flanking_parent_ids)
        flanking_hits_by_parent = defaultdict(list)
        flanking_locus_gap = max(5000, min(100000, args.max_intron * 2))

        for hit in relevant_hits:
            query_id = hit.get('query', '')
            if not query_id or is_goi_query_id(query_id):
                continue

            parent_id = extract_base_gene_id(query_id)
            if parent_id not in flanking_query_by_parent:
                continue

            local_start = max(0, int(hit.get('start', 0)) - w_start)
            local_end = min(len(subseq), int(hit.get('end', 0)) - w_start)
            if local_end <= local_start:
                continue

            global_start = w_start + local_start
            global_end = w_start + local_end
            if global_end <= core_start or global_start >= core_end:
                continue

            flanking_hits_by_parent[parent_id].append({
                'query': query_id,
                'qstart': hit.get('qstart', 1),
                'qend': hit.get('qend', 1),
                'gstart': local_start,
                'gend': local_end,
                'evalue': hit.get('evalue', 1),
                'pident': hit.get('pident', 0),
                'alnlen': hit.get('alnlen', 0),
                'bits': hit.get('bits', 0),
                'strand': hit.get('strand', '+'),
                'chrom': chrom
            })

        # Remap per-exon qstart/qend to reconstructed full-protein coordinates.
        # Without this, a hit from gene-TOP1MT|exon_5 (qstart=1,qend=37) would
        # be compared against the 601aa full protein → qcov=6% → wrongly fails
        # coverage thresholds.  After remapping, it becomes qstart=170,qend=207
        # which correctly reflects its position in the reconstructed protein.
        for parent_id, parent_hits in flanking_hits_by_parent.items():
            query_rec = flanking_query_by_parent.get(parent_id)
            if not query_rec:
                continue
            exon_offsets = query_rec.get('exon_offsets', {})
            if not exon_offsets:
                continue
            for h in parent_hits:
                qid = h.get('query', '')
                m = re.search(r'\|exon_(\d+)$', qid)
                if not m:
                    continue
                exon_num = int(m.group(1))
                if exon_num not in exon_offsets:
                    continue
                offset, _exon_len = exon_offsets[exon_num]
                # Shift from exon-local to full-protein coordinates (1-based)
                h['qstart'] = h.get('qstart', 1) + offset
                h['qend'] = h.get('qend', 1) + offset

        for parent_id, parent_hits in flanking_hits_by_parent.items():
            query_rec = flanking_query_by_parent.get(parent_id)
            if not query_rec:
                continue
            parent_query_seq = query_rec.get('seq', '')
            if not parent_query_seq:
                continue

            parent_loci = split_hits_into_loci(parent_hits, max_gap=flanking_locus_gap)
            if not parent_loci:
                continue

            # Keep strongest loci first; flanking is anchor support, not exhaustive paralog mining.
            parent_loci = sorted(
                parent_loci,
                key=lambda hs: (max(h.get('bits', 0) for h in hs), len(hs)),
                reverse=True
            )[:2]

            for locus_idx, work_hits in enumerate(parent_loci, start=1):
                try:
                    with maybe_quiet_streams(args.quiet_subtools):
                        exons, _ = annotate_exons_from_hit_list(
                            work_hits,
                            parent_query_seq,
                            subseq,
                            chrom,
                            search_missing=True,
                            gap_min_size=args.gap_min_size,
                            gap_search_window=args.gap_search_window,
                            gap_evalue=args.gap_evalue,
                            gap_min_identity=args.gap_min_identity,
                            gap_min_alnlen=args.gap_min_alnlen,
                            gap_max_hits=args.gap_max_hits,
                            exon_query_mode=False,
                            min_exon_query_cov=args.min_exon_query_cov,
                            min_exon_alnlen=args.min_exon_alnlen,
                            sensitive=False
                        )
                except Exception as e:
                    logger.warning(f"[{genome_name}] Flanking exon annotation failed for {parent_id} locus {locus_idx}: {e}")
                    exons = []

                if not exons:
                    # Conservative fallback: use hit-defined exon blocks from first to last hit.
                    if not work_hits:
                        continue

                    strand_votes = {'+': 0, '-': 0}
                    for h in work_hits:
                        strand_votes[h.get('strand', '+')] += 1
                    strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                    ordered_hits = sorted(
                        [h for h in work_hits if h.get('strand', '+') == strand],
                        key=lambda h: h.get('gstart', 0)
                    )
                    if not ordered_hits:
                        continue

                    # Keep only the longest chain consistent with query-order direction.
                    ordered_hits = _longest_monotonic_query_chain(ordered_hits, strand)
                    if not ordered_hits:
                        continue

                    # Require minimum support before emitting coarse flanking fallback models.
                    qmin = min(min(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                    qmax = max(max(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                    query_len = max(1, len(parent_query_seq))
                    qcov = (qmax - qmin + 1) / query_len
                    aln_total = sum(max(0, int(h.get('alnlen', 0))) for h in ordered_hits)
                    best_bits = max(float(h.get('bits', 0)) for h in ordered_hits)

                    # Single-hit fallback must look like substantial coverage.
                    if len(ordered_hits) == 1:
                        min_single_qcov = 0.45
                        min_single_aln = max(30, int(0.35 * query_len))
                        if qcov < min_single_qcov or int(ordered_hits[0].get('alnlen', 0)) < min_single_aln:
                            continue
                    else:
                        # Multi-hit chain fallback: require meaningful overall support.
                        if qcov < 0.35 and aln_total < max(45, int(0.45 * query_len)) and best_bits < 80.0:
                            continue

                    coding_frags = []
                    cds_intervals = []
                    for h in ordered_hits:
                        hs = h.get('gstart', 0)
                        he = h.get('gend', 0)
                        if he <= hs:
                            continue
                        dna = subseq[hs:he]
                        if strand == '-':
                            dna = reverse_complement(dna)
                        dna = dna[:len(dna) - (len(dna) % 3)]
                        if len(dna) < 9:
                            continue
                        aa = translate(dna).replace('*', '')
                        if not aa:
                            continue
                        coding_frags.append(aa)
                        cds_intervals.append((hs, he))

                    if not coding_frags or not cds_intervals:
                        continue

                    if strand == '-':
                        coding_frags = list(reversed(coding_frags))

                    flank_protein = ''.join(coding_frags)
                    global_start = w_start + min(s for s, _ in cds_intervals) + 1
                    global_end = w_start + max(e for _, e in cds_intervals)
                    avg_pident = sum(h.get('pident', 0) for h in ordered_hits) / len(ordered_hits)
                    new_id = f"{parent_id}|{clean_gname}_b{block_idx}_fl{locus_idx}_flank_hits"

                    cand_gff = [
                        (
                            f"{chrom}\tflanking_hits\tmRNA\t{global_start}\t{global_end}\t"
                            f"{avg_pident:.1f}\t{strand}\t.\t"
                            f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'flanking_hit_span'}, evidence_type='flanking_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
                        )
                    ]
                    for eidx, (hs, he) in enumerate(cds_intervals, 1):
                        exon_gs = w_start + hs + 1
                        exon_ge = w_start + he
                        cand_gff.append(
                            f"{chrom}\tflanking_hits\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\tID={new_id}_CDS{eidx};Parent={new_id}"
                        )

                    flanking_candidates.append({
                        'start': global_start,
                        'end': global_end,
                        'score': avg_pident,
                        'record': {
                            'id': new_id,
                            'seq': flank_protein,
                            'description': (
                                f"coords:{global_start}-{global_end} parent:{parent_id} "
                                f"type:flanking_hit_span identity:{avg_pident:.1f}"
                            )
                        },
                        'gff': cand_gff
                    })
                    continue

                exons.sort(key=lambda e: e.get('qstart', 0))
                flank_protein = ''.join(e.get('seq', '') for e in exons if e.get('seq'))
                if not flank_protein:
                    continue

                strand = exons[0].get('strand', '+')
                avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                model_qcov = None
                if parent_query_seq:
                    qmin_model = min(e.get('qstart', 1) for e in exons)
                    qmax_model = max(e.get('qend', 1) for e in exons)
                    model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))
                global_start = w_start + min(e['gstart'] for e in exons) + 1
                global_end = w_start + max(e['gend'] for e in exons)
                new_id = f"{parent_id}|{clean_gname}_b{block_idx}_fl{locus_idx}_flank_ann"

                cand_gff = [
                    (
                        f"{chrom}\tflanking_annotation\tmRNA\t{global_start}\t{global_end}\t"
                        f"{avg_pident:.1f}\t{strand}\t.\t"
                        f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'flanking_miniprot'}, evidence_type='flanking_miniprot', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
                    )
                ]
                for eidx, exon in enumerate(exons, 1):
                    exon_gs = w_start + exon['gstart'] + 1
                    exon_ge = w_start + exon['gend']
                    cand_gff.append(
                        f"{chrom}\tflanking_annotation\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\tID={new_id}_CDS{eidx};Parent={new_id}"
                    )

                flanking_candidates.append({
                    'start': global_start,
                    'end': global_end,
                    'score': avg_pident,
                    'record': {
                        'id': new_id,
                        'seq': flank_protein,
                        'description': (
                            f"coords:{global_start}-{global_end} parent:{parent_id} "
                            f"exons:{len(exons)} identity:{avg_pident:.1f}"
                        )
                    },
                    'gff': cand_gff
                })

    except Exception as e:
        import traceback
        logger.warning(
            f"[{genome_name}] Block {block_idx} processing failed "
            f"({chrom}:{block.get('start','?')}-{block.get('end','?')}): {e}\n"
            f"  {traceback.format_exc().strip().splitlines()[-1]}"
        )

    if os.path.exists(temp_fa):
        os.remove(temp_fa)
    if os.path.exists(query_mini_fa):
        os.remove(query_mini_fa)

    # ── PLM re-ranking: batch-embed GOI model proteins, compute embedding
    #    similarity, re-run classification where PLM boost applies, and
    #    update GFF attributes accordingly. ──
    if (getattr(args, 'enable_plm_search', False)
            and PLM_IMPORT_OK and check_plm_available()
            and raw_candidates):
        try:
            goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
            if os.path.exists(goi_emb_path):
                goi_embeddings = load_embeddings(goi_emb_path)
                goi_cand_seqs = [
                    (cand['record']['id'], cand['record']['seq'])
                    for cand in raw_candidates
                    if cand['record'].get('seq') and len(cand['record']['seq']) >= 5
                ]
                if goi_cand_seqs and goi_embeddings:
                    cand_sims = compute_candidate_similarities(
                        goi_cand_seqs, goi_embeddings,
                        device=getattr(args, 'plm_device', 'cpu'),
                    )
                    boosted = 0
                    for cand in raw_candidates:
                        cid = cand['record']['id']
                        sim = cand_sims.get(cid)
                        if sim is None:
                            continue
                        cand['embedding_similarity'] = sim

                        if not cand.get('gff'):
                            continue
                        mRNA_line = cand['gff'][0]

                        # Parse existing classification from GFF attrs to re-run
                        # with embedding_similarity.
                        attrs_str = mRNA_line.split('\t')[-1] if '\t' in mRNA_line else ""
                        attr_kv = {}
                        for pair in attrs_str.split(';'):
                            if '=' in pair:
                                k, v = pair.split('=', 1)
                                attr_kv[k.strip()] = v.strip()

                        old_conf = attr_kv.get('Confidence', '')
                        ev_type = attr_kv.get('EvidenceType', '')
                        try:
                            ident = float(attr_kv.get('Identity', 0))
                        except (ValueError, TypeError):
                            ident = 0.0
                        try:
                            exons = int(attr_kv.get('Exons', 1))
                        except (ValueError, TypeError):
                            exons = 1
                        try:
                            qcov = float(attr_kv.get('QueryCoverage', 0))
                        except (ValueError, TypeError):
                            qcov = None
                        try:
                            fl_sup = int(attr_kv.get('BlockFlankingSupport', 0))
                        except (ValueError, TypeError):
                            fl_sup = 0

                        # Re-classify with embedding similarity
                        new_conf, new_class, new_reason = _classify_goi_evidence(
                            evidence_type=ev_type,
                            identity=ident,
                            exon_count=exons,
                            query_cov=qcov,
                            flanking_support=fl_sup,
                            embedding_similarity=sim,
                        )

                        # Update GFF attrs if confidence was boosted
                        updated_attrs = attrs_str + f";EmbeddingSimilarity={sim:.3f}"
                        if new_conf != old_conf and old_conf in ('LOW', 'MEDIUM'):
                            updated_attrs = updated_attrs.replace(
                                f"Confidence={old_conf}", f"Confidence={new_conf}"
                            )
                            old_reason = attr_kv.get('InferenceReason', '')
                            if old_reason:
                                updated_attrs = updated_attrs.replace(
                                    f"InferenceReason={old_reason}",
                                    f"InferenceReason={new_reason}",
                                )
                            old_class = attr_kv.get('GOIClass', '')
                            if old_class:
                                updated_attrs = updated_attrs.replace(
                                    f"GOIClass={old_class}",
                                    f"GOIClass={new_class}",
                                )
                            boosted += 1

                        # Rebuild the mRNA line with updated attributes
                        cols = mRNA_line.split('\t')
                        if len(cols) >= 9:
                            cols[8] = updated_attrs
                            cand['gff'][0] = '\t'.join(cols)

                    if boosted:
                        print(
                            f"[{genome_name}] PLM re-ranking: boosted {boosted}/{len(raw_candidates)} "
                            f"candidate(s) via embedding similarity.",
                            flush=True,
                        )
        except Exception as plm_rerank_err:
            logger.debug(
                f"[{genome_name}] PLM re-ranking failed: {plm_rerank_err}"
            )

    # ── Structural re-ranking: fold GOI candidate proteins with ESMFold,
    #    compare against GOI structures via Foldseek, re-run classification
    #    where structural boost applies, and update GFF attributes. ──
    if (getattr(args, 'enable_structural_search', False)
            and STRUCTURAL_IMPORT_OK and check_structural_search_available()
            and raw_candidates):
        try:
            goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
            if os.path.exists(goi_struct_index):
                goi_structures = load_structure_index(goi_struct_index)
                goi_cand_seqs = [
                    (cand['record']['id'], cand['record']['seq'])
                    for cand in raw_candidates
                    if cand['record'].get('seq') and len(cand['record']['seq']) >= 10
                ]
                if goi_cand_seqs and goi_structures:
                    struct_rerank_dir = os.path.join(
                        args.output_dir, f"tmp_struct_rerank_{genome_name}"
                    )
                    os.makedirs(struct_rerank_dir, exist_ok=True)

                    cand_tm_scores = compute_candidate_structural_similarities(
                        goi_cand_seqs, goi_structures,
                        output_dir=struct_rerank_dir,
                        device=getattr(args, 'structural_device', 'cpu'),
                        max_length=getattr(args, 'structural_max_length', 700),
                        threads=getattr(args, 'threads', 1),
                    )

                    boosted = 0
                    for cand in raw_candidates:
                        cid = cand['record']['id']
                        tm = cand_tm_scores.get(cid)
                        if tm is None:
                            continue
                        cand['structural_similarity'] = tm

                        if not cand.get('gff'):
                            continue
                        mRNA_line = cand['gff'][0]

                        attrs_str = mRNA_line.split('\t')[-1] if '\t' in mRNA_line else ""
                        attr_kv = {}
                        for pair in attrs_str.split(';'):
                            if '=' in pair:
                                k, v = pair.split('=', 1)
                                attr_kv[k.strip()] = v.strip()

                        old_conf = attr_kv.get('Confidence', '')
                        ev_type = attr_kv.get('EvidenceType', '')
                        try:
                            ident = float(attr_kv.get('Identity', 0))
                        except (ValueError, TypeError):
                            ident = 0.0
                        try:
                            exons = int(attr_kv.get('Exons', 1))
                        except (ValueError, TypeError):
                            exons = 1
                        try:
                            qcov = float(attr_kv.get('QueryCoverage', 0))
                        except (ValueError, TypeError):
                            qcov = None
                        try:
                            fl_sup = int(attr_kv.get('BlockFlankingSupport', 0))
                        except (ValueError, TypeError):
                            fl_sup = 0
                        try:
                            emb_sim = float(attr_kv.get('EmbeddingSimilarity', 0))
                        except (ValueError, TypeError):
                            emb_sim = None

                        new_conf, new_class, new_reason = _classify_goi_evidence(
                            evidence_type=ev_type,
                            identity=ident,
                            exon_count=exons,
                            query_cov=qcov,
                            flanking_support=fl_sup,
                            embedding_similarity=emb_sim if emb_sim else None,
                            structural_similarity=tm,
                        )

                        updated_attrs = attrs_str + f";StructuralSimilarity={tm:.3f}"
                        if new_conf != old_conf and old_conf in ('LOW', 'MEDIUM'):
                            updated_attrs = updated_attrs.replace(
                                f"Confidence={old_conf}", f"Confidence={new_conf}"
                            )
                            old_reason = attr_kv.get('InferenceReason', '')
                            if old_reason:
                                updated_attrs = updated_attrs.replace(
                                    f"InferenceReason={old_reason}",
                                    f"InferenceReason={new_reason}",
                                )
                            old_class = attr_kv.get('GOIClass', '')
                            if old_class:
                                updated_attrs = updated_attrs.replace(
                                    f"GOIClass={old_class}",
                                    f"GOIClass={new_class}",
                                )
                            boosted += 1

                        cols = mRNA_line.split('\t')
                        if len(cols) >= 9:
                            cols[8] = updated_attrs
                            cand['gff'][0] = '\t'.join(cols)

                    if boosted:
                        print(
                            f"[{genome_name}] Structural re-ranking: boosted {boosted}/{len(raw_candidates)} "
                            f"candidate(s) via TM-score.",
                            flush=True,
                        )

                    if os.path.exists(struct_rerank_dir):
                        shutil.rmtree(struct_rerank_dir, ignore_errors=True)
        except Exception as struct_rerank_err:
            logger.debug(
                f"[{genome_name}] Structural re-ranking failed: {struct_rerank_err}"
            )

    found_ids = set()

    # Greedy non-overlapping selection for GOI/tandem candidates.
    if raw_candidates:
        raw_candidates.sort(
            key=lambda x: (x.get('score', 0), x['end'] - x['start']),
            reverse=True
        )

        final_candidates = []
        accepted_intervals = []

        for cand in raw_candidates:
            c_start, c_end = cand['start'], cand['end']
            c_len = c_end - c_start
            if c_len <= 0:
                continue

            is_overlapping = False
            for a_start, a_end in accepted_intervals:
                overlap = max(0, min(c_end, a_end) - max(c_start, a_start))
                if overlap > 0.5 * c_len:
                    is_overlapping = True
                    break

            if not is_overlapping:
                final_candidates.append(cand)
                accepted_intervals.append((c_start, c_end))

        for cand in final_candidates:
            rec = cand['record']
            if rec['id'] in found_ids:
                continue
            found_ids.add(rec['id'])
            annotated_records_raw.append(rec)
            valid_gff_lines.extend(cand['gff'])

    # Add flanking candidates without interfering with GOI selection.
    if flanking_candidates:
        flanking_candidates.sort(
            key=lambda x: (x.get('score', 0), x['end'] - x['start']),
            reverse=True
        )
        flanking_intervals = []
        for cand in flanking_candidates:
            rec = cand['record']
            if rec['id'] in found_ids:
                continue

            c_start, c_end = cand['start'], cand['end']
            c_len = c_end - c_start
            if c_len <= 0:
                continue

            is_overlapping = False
            for a_start, a_end in flanking_intervals:
                overlap = max(0, min(c_end, a_end) - max(c_start, a_start))
                if overlap > 0.5 * c_len:
                    is_overlapping = True
                    break
            if is_overlapping:
                continue

            flanking_intervals.append((c_start, c_end))
            found_ids.add(rec['id'])
            annotated_records_raw.append(rec)
            valid_gff_lines.extend(cand['gff'])

    return annotated_records_raw, valid_gff_lines


def merge_synteny_blocks(blocks, padding):
    """
    Merge overlapping synteny blocks to prevent redundant searches.
    
    If the search windows (block +/- padding) overlap, the blocks are merged.
    """
    if not blocks:
        return []

    # Sort by chrom, then start
    # Ensure all blocks have 'chrom' and 'start'
    valid_blocks = [b for b in blocks if 'chrom' in b and 'start' in b]
    if len(valid_blocks) < len(blocks):
        print(f"Warning: Dropped {len(blocks) - len(valid_blocks)} blocks missing coordinates during merge.")
    
    valid_blocks.sort(key=lambda b: (b['chrom'], b['start']))

    merged = []
    if not valid_blocks:
        return []
        
    current_block = valid_blocks[0]
    current_chrom = current_block['chrom']
    # Effective end of the search region for the current block
    current_search_end = current_block['end'] + padding

    for i in range(1, len(valid_blocks)):
        next_block = valid_blocks[i]
        next_search_start = next_block['start'] - padding
        
        # Check for overlap
        if next_block['chrom'] == current_chrom and next_search_start <= current_search_end:
            # Merge!
            # 1. Update core coordinates to cover both blocks
            current_block['start'] = min(current_block['start'], next_block['start'])
            current_block['end'] = max(current_block['end'], next_block['end'])
            
            # 2. Update stats
            # Sum gene counts (approximate, but safe for filtering)
            current_block['genes_count'] = current_block.get('genes_count', 0) + next_block.get('genes_count', 0)
            # Max score (best local evidence)
            current_block['score'] = max(current_block.get('score', 0), next_block.get('score', 0))
            
            # 3. Combine hits if present
            if 'hits' in next_block:
                if 'hits' not in current_block:
                    current_block['hits'] = []
                current_block['hits'].extend(next_block['hits'])
            
            # 4. Extend the search horizon
            current_search_end = max(current_search_end, next_block['end'] + padding)
            
            # Note: We keep modifying 'current_block' in place
        else:
            # No overlap, finalize current and move to next
            merged.append(current_block)
            current_block = next_block
            current_chrom = current_block['chrom']
            current_search_end = current_block['end'] + padding
            
    merged.append(current_block)
    return merged


def run_mmseqs_easy_search_with_retries(
    db_path: str,
    genome_path: str,
    hits_file: str,
    tmp_dir: str,
    args,
    threads_per_job: int,
    genome_name: str,
) -> Tuple[bool, str, bool]:
    """
    Run mmseqs easy-search with conservative retries for resource failures.

    Returns:
      (success, details, resource_related_failure)
    """
    base_attempts = [
        {
            "label": "primary",
            "threads": max(1, int(threads_per_job)),
            "sens": float(args.mmseqs_sens),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "lowmem_retry_1",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 7.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
        {
            "label": "lowmem_retry_2",
            "threads": 1,
            "sens": min(float(args.mmseqs_sens), 6.0),
            "split": str(args.mmseqs_split_memory_limit),
        },
    ]

    attempts = []
    seen = set()
    for a in base_attempts:
        key = (a["threads"], a["sens"], a["split"])
        if key in seen:
            continue
        seen.add(key)
        attempts.append(a)

    last_details = ""
    last_resource_fail = False
    os.makedirs(tmp_dir, exist_ok=True)

    # Create query and target databases ONCE (outside retry loop).
    query_db = os.path.join(tmp_dir, "queryDB")
    target_db = os.path.join(tmp_dir, "targetDB")
    fmt_output = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

    # -- createdb for query (protein FASTA) --
    proc_qdb = subprocess.run(
        ["mmseqs", "createdb", db_path, query_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_qdb.returncode != 0:
        details = _format_mmseqs_failure(proc_qdb)
        logger.warning(f"[{genome_name}] mmseqs createdb (query) failed: {details}")
        return False, details, False

    # -- createdb for target (nucleotide genome) --
    proc_tdb = subprocess.run(
        ["mmseqs", "createdb", genome_path, target_db],
        check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if proc_tdb.returncode != 0:
        details = _format_mmseqs_failure(proc_tdb)
        logger.warning(f"[{genome_name}] mmseqs createdb (target) failed: {details}")
        return False, details, False

    for i, attempt in enumerate(attempts, start=1):
        if os.path.exists(hits_file):
            try:
                os.remove(hits_file)
            except OSError:
                pass

        attempt_tmp = os.path.join(tmp_dir, f"attempt_{i}")
        os.makedirs(attempt_tmp, exist_ok=True)
        result_db = os.path.join(attempt_tmp, "resultDB")

        # -- mmseqs search (protein query vs translated nucleotide target) --
        search_cmd = [
            "mmseqs", "search",
            query_db, target_db, result_db, attempt_tmp,
            "--search-type", "2",
            "--threads", str(attempt["threads"]),
            "-s", str(attempt["sens"]),
            "-e", str(args.evalue),
            "--split-memory-limit", str(attempt["split"]),
            "-v", str(args.mmseqs_verbosity),
        ]

        search_proc = subprocess.run(
            search_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if search_proc.returncode != 0:
            details = _format_mmseqs_failure(search_proc)
            resource_fail = _is_mmseqs_resource_failure(details)
            last_details = details
            last_resource_fail = resource_fail
            logger.warning(
                f"[{genome_name}] MMseqs search attempt {i}/{len(attempts)} failed "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']}): {details}"
            )
            if not resource_fail:
                break
            continue

        # -- mmseqs convertalis → m8 output --
        conv_cmd = [
            "mmseqs", "convertalis",
            query_db, target_db, result_db, hits_file,
            "--format-output", fmt_output,
            "-v", str(args.mmseqs_verbosity),
        ]

        conv_proc = subprocess.run(
            conv_cmd,
            check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

        if conv_proc.returncode != 0:
            details = _format_mmseqs_failure(conv_proc)
            last_details = details
            logger.warning(f"[{genome_name}] mmseqs convertalis failed: {details}")
            break

        if i > 1:
            logger.warning(
                f"[{genome_name}] MMseqs succeeded on retry {i}/{len(attempts)} "
                f"({attempt['label']}: threads={attempt['threads']}, "
                f"s={attempt['sens']}, split={attempt['split']})."
            )
        return True, "", False

    return False, last_details, last_resource_fail


_SEED_QUALIFYING_CLASSES = {"confident_goi", "probable_goi"}
_TREE_QUALIFYING_CONFIDENCE = {"HIGH", "MEDIUM"}


def _classify_goi_for_seed_and_tree(all_genes, feature_meta, is_goi_query_id_fn):
    """Split GOI-role genes into (seed, tree_extra, suppressed_count).

    Seed genes feed `expanded_db.faa` (used as the next wave's MMseqs2 query DB
    AND historically as the tree input). Strict filter: HIGH/MEDIUM confidence
    AND goi_class in {confident_goi, probable_goi}.

    Tree-extra genes feed `goi_for_tree.faa` *in addition to* the seed list:
    HIGH/MEDIUM-confidence GOI hits whose goi_class disqualified them from
    seeding (typically `tandem_goi_copy`). They appear in MAFFT/IQ-TREE so
    species like Apis cerana — which only carry tandem-duplicate evidence —
    are not silently dropped from the phylogeny.

    Suppressed_count is the total of GOI-role genes that did NOT qualify for
    seeding (used for logging).
    """
    seed = []
    tree_extra = []
    suppressed = 0
    for g in all_genes:
        gid = g.get("id", "")
        meta = feature_meta.get(gid, {})
        role = meta.get("role", "goi" if is_goi_query_id_fn(gid) else "flanking")
        if role != "goi":
            continue
        confidence = meta.get("confidence", "")
        goi_class = meta.get("goi_class", "")
        if (confidence in _TREE_QUALIFYING_CONFIDENCE
                and goi_class in _SEED_QUALIFYING_CLASSES):
            seed.append(g)
        else:
            suppressed += 1
            if confidence in _TREE_QUALIFYING_CONFIDENCE:
                tree_extra.append(g)
    return seed, tree_extra, suppressed


def process_single_genome(genome_path, db_path, args, home_db_dir, prefix, threads_per_job):
    """
    Worker function to search a single genome.
    Returns: (genome_name, list_of_new_genes, list_of_tree_extra_genes,
             error_message_or_none)
    """
    genome_name = os.path.basename(genome_path)
    if not os.path.exists(genome_path):
        msg = "Genome file not found"
        logger.error(f"[{genome_name}] {msg}.")
        return genome_name, [], [], msg
    
    unique_id = uuid.uuid4().hex
    hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
    tmp_dir = f"{args.output_dir}/tmp_mmseqs_{unique_id}_{genome_name}"
    
    new_genes = []
    tree_extra_genes = []
    error_message = None

    try:
        # Parse DB once so we can adapt hit-length filtering for short proteins.
        base_db_sequences = {}
        query_lengths = {}
        for header, clean_id, seq in parse_fasta(db_path):
            base_db_sequences[clean_id] = {'id': clean_id, 'seq': seq, 'header': header}
            query_lengths[clean_id] = len(seq)
            base = extract_base_gene_id(clean_id)
            if base != clean_id and base not in base_db_sequences:
                base_db_sequences[base] = {'id': clean_id, 'seq': seq, 'header': header}
            query_lengths.setdefault(base, len(seq))

        # Identify flanking proteins that are proxies for the GOI query
        # (e.g., gene-LY6E is a proxy for GOI_LY6E — same gene, different prefix).
        # Proxy blocks seeded by these genes will bypass the min_block_genes filter.
        _goi_bare_names: set = set()
        for _rid in list(base_db_sequences.keys()):
            if is_goi_query_id(_rid) and '|' not in _rid:
                _bare = _rid
                if _bare.startswith('GOI_copy_'):
                    _bare = _bare[9:]
                elif _bare.startswith('GOI_'):
                    _bare = _bare[4:]
                _goi_bare_names.add(_bare)
        goi_proxy_flanking_parents: set = set()
        for _rid in list(base_db_sequences.keys()):
            if is_goi_query_id(_rid) or '|' in _rid:
                continue
            _bare = _rid[5:] if _rid.startswith('gene-') else _rid
            if _bare in _goi_bare_names:
                goi_proxy_flanking_parents.add(_rid)
        if goi_proxy_flanking_parents:
            logger.info(
                f"[{genome_name}] GOI proxy flanking parents: "
                f"{sorted(goi_proxy_flanking_parents)}"
            )

        c_dist = args.cluster_distance
        if c_dist <= 0:
            c_dist = estimate_cluster_distance(genome_path)
            
        # 1. Search (MMseqs) with low-memory retries for fragmented genomes.
        mmseqs_ok, mmseqs_details, resource_fail = run_mmseqs_easy_search_with_retries(
            db_path=db_path,
            genome_path=genome_path,
            hits_file=hits_file,
            tmp_dir=tmp_dir,
            args=args,
            threads_per_job=threads_per_job,
            genome_name=genome_name,
        )
        if not mmseqs_ok:
            if resource_fail:
                logger.warning(
                    f"[{genome_name}] Skipping genome after MMseqs resource failures "
                    "(prefilter/search crash persisted across retries)."
                )
                return genome_name, [], [], None
            raise RuntimeError(
                f"MMseqs easy-search failed for {genome_name}: {mmseqs_details}"
            )

        hits = parse_hits(
            hits_file,
            args.min_identity,
            args.min_length,
            args.evalue,
            query_lengths=query_lengths,
        )
        if not hits:
            logger.info(f"[{genome_name}] No hits found in MMseqs output.")
            return genome_name, [], [], None
            
        logger.info(f"[{genome_name}] Parsed {len(hits)} hits.")

        # 2. Identify Synteny (Multiple Blocks)
        # Anchor blocks with non-GOI (flanking) hits first, then search GOI inside them.
        # This follows the intended flow: flanking -> region -> GOI modeling.
        flanking_seed_hits = [h for h in hits if not is_goi_query_id(h.get('query', ''))]
        seed_hits = flanking_seed_hits if flanking_seed_hits else hits
        if flanking_seed_hits:
            logger.info(f"[{genome_name}] Using {len(flanking_seed_hits)} flanking-anchor hits for block seeding.")
        else:
            logger.info(f"[{genome_name}] No flanking anchors found; falling back to GOI/mixed hits for block seeding.")

        synteny_blocks = identify_synteny_blocks(
            seed_hits,
            max_intron=args.max_intron,
            cluster_distance=c_dist
        )
        
        if not synteny_blocks:
            logger.info(f"[{genome_name}] No valid syntenic region found.")
            return genome_name, [], [], None
            
        # Optimization: Merge overlapping search regions
        pre_merge_count = len(synteny_blocks)
        synteny_blocks = merge_synteny_blocks(synteny_blocks, args.region_padding)
        if len(synteny_blocks) < pre_merge_count:
            logger.info(f"[{genome_name}] Merged {pre_merge_count} blocks into {len(synteny_blocks)} discrete search regions.")
        else:
            logger.info(f"[{genome_name}] Found {len(synteny_blocks)} discrete syntenic blocks.")

        # Keep only blocks with enough anchor genes to avoid spending hours
        # on singleton/noise loci in fragmented genomes.
        # Exception: single-gene blocks seeded by a GOI proxy (e.g., gene-LY6E)
        # are retained so the full GOI search runs on the correct scaffold.
        if args.min_block_genes > 1:
            pre_filter_count = len(synteny_blocks)
            filtered_blocks = []
            proxy_blocks = []
            for _b in synteny_blocks:
                if _b.get('genes_count', 0) >= args.min_block_genes:
                    filtered_blocks.append(_b)
                elif goi_proxy_flanking_parents and any(
                    _g in goi_proxy_flanking_parents for _g in _b.get('genes', [])
                ):
                    proxy_blocks.append(_b)
            synteny_blocks = filtered_blocks + proxy_blocks
            logger.info(
                f"[{genome_name}] Block filter (min_block_genes={args.min_block_genes}): "
                f"{len(filtered_blocks)} standard + {len(proxy_blocks)} GOI-proxy "
                f"= {len(synteny_blocks)}/{pre_filter_count} retained."
            )

        # Fallback: if filtering removed everything, still keep top blocks by score.
        if not synteny_blocks:
            logger.info(
                f"[{genome_name}] No blocks after filtering. "
                f"Falling back to top {args.max_blocks_per_genome} unfiltered blocks."
            )
            synteny_blocks = identify_synteny_blocks(
                seed_hits,
                max_intron=args.max_intron,
                cluster_distance=c_dist
            )
            synteny_blocks = merge_synteny_blocks(synteny_blocks, args.region_padding)

        # Hard cap per genome to prevent pathological runtimes.
        if args.max_blocks_per_genome > 0 and len(synteny_blocks) > args.max_blocks_per_genome:
            logger.info(
                f"[{genome_name}] Capping blocks: {len(synteny_blocks)} -> "
                f"{args.max_blocks_per_genome} (max_blocks_per_genome)."
            )
            synteny_blocks = synteny_blocks[:args.max_blocks_per_genome]
        
        all_genes = []
        all_gff_lines = []
        
        genome_seqs = load_genome(genome_path)
        native_gff_path = find_native_annotation_path(genome_path)
        native_annot_index = load_native_annotation_index(native_gff_path)
        if native_gff_path and native_annot_index:
            logger.info(
                f"[{genome_name}] Native annotation index loaded from "
                f"{os.path.basename(native_gff_path)} ({len(native_annot_index)} contigs)."
            )
        elif native_gff_path:
            logger.info(
                f"[{genome_name}] Native annotation file found ({os.path.basename(native_gff_path)}) "
                f"but no usable gene/transcript features parsed."
            )
        else:
            logger.info(f"[{genome_name}] No native annotation GFF found next to genome FASTA.")

        def _target_mrna_attrs(
            chrom_name: str,
            g_start: int,
            g_end: int,
            g_strand: str,
            base_attrs: Dict[str, Any],
        ) -> str:
            native = lookup_native_annotation(
                native_annot_index or {},
                chrom=chrom_name,
                start=g_start,
                end=g_end,
                strand=g_strand,
            )
            return _compose_gff_attrs(base_attrs, native)

        empty_block_streak = 0
        for i, block in enumerate(synteny_blocks):
            # Limit number of blocks to process? (e.g. top 50?)
            # For now process all valid blocks.
            block_genes, block_gff = process_region_block(
                i,
                block,
                hits,
                genome_seqs,
                base_db_sequences,
                genome_name,
                args,
                unique_id,
                threads_per_job,
                native_annot_index=native_annot_index,
            )
            all_genes.extend(block_genes)
            all_gff_lines.extend(block_gff)

            block_goi = sum(1 for g in block_genes if is_goi_query_id(g.get('id', '')))
            if block_goi > 0:
                empty_block_streak = 0
            else:
                empty_block_streak += 1

            if (i + 1) % 10 == 0:
                logger.info(
                    f"[{genome_name}] Block progress: {i + 1}/{len(synteny_blocks)} "
                    f"(empty_streak={empty_block_streak})."
                )

            if (
                args.max_consecutive_empty_blocks > 0 and
                empty_block_streak >= args.max_consecutive_empty_blocks
            ):
                logger.info(
                    f"[{genome_name}] Early stop logic disabled to ensure complete block evaluation (streak={empty_block_streak})."
                )
                # break

        # Flanking-only post-pass: keep one best chain per flanking parent per locus across blocks.
        all_genes, all_gff_lines = deduplicate_flanking_models(
            all_genes,
            all_gff_lines,
            genome_name=genome_name,
            locus_gap_bp=max(5000, int(args.cluster_distance)),
        )
        all_gff_lines = collapse_flanking_cds_to_gene_span(all_gff_lines)

        # --- Cross-chromosome flanking recovery ---
        # Some flanking genes may have translocated to a different chromosome.
        # After standard block processing, check which flanking parents are still
        # missing, look for strong hits on other chromosomes, and annotate them
        # with a "rearranged" flag so the plot can show them.
        all_db_flanking_parents = {
            extract_base_gene_id(rid)
            for rid in base_db_sequences
            if not is_goi_query_id(rid) and '|' not in rid
        }
        found_flanking_parents = set()
        for g in all_genes:
            gid = g.get('id', '')
            if not is_goi_query_id(gid):
                found_flanking_parents.add(extract_base_gene_id(gid))

        missing_flanking = all_db_flanking_parents - found_flanking_parents
        block_chroms = {b['chrom'] for b in synteny_blocks}

        if missing_flanking and hits:
            # Gather off-block-chrom hits for missing parents
            off_chrom_hits_by_parent = defaultdict(list)
            for h in hits:
                parent = extract_base_gene_id(h.get('query', ''))
                if parent not in missing_flanking:
                    continue
                if h['chrom'] in block_chroms:
                    continue  # Already considered during block processing
                off_chrom_hits_by_parent[parent].append(h)

            if off_chrom_hits_by_parent:
                logger.info(
                    f"[{genome_name}] Cross-chromosome flanking recovery: "
                    f"{len(off_chrom_hits_by_parent)} missing parents have off-block hits."
                )
                flanking_query_map = build_flanking_query_by_parent(
                    base_db_sequences, set(off_chrom_hits_by_parent.keys())
                )
                clean_gname = genome_name.replace('.', '_').replace('-', '_').replace(' ', '_')
                flanking_locus_gap = max(5000, min(100000, args.max_intron * 2))

                for parent_id, parent_hits in off_chrom_hits_by_parent.items():
                    query_rec = flanking_query_map.get(parent_id)
                    if not query_rec:
                        continue
                    parent_query_seq = query_rec.get('seq', '')
                    if not parent_query_seq:
                        continue

                    # Remap exon-level qstart/qend
                    exon_offsets = query_rec.get('exon_offsets', {})
                    if exon_offsets:
                        for h in parent_hits:
                            m = re.search(r'\|exon_(\d+)$', h.get('query', ''))
                            if m:
                                exon_num = int(m.group(1))
                                if exon_num in exon_offsets:
                                    offset, _ = exon_offsets[exon_num]
                                    h['qstart'] = h.get('qstart', 1) + offset
                                    h['qend'] = h.get('qend', 1) + offset

                    parent_loci = split_hits_into_loci(parent_hits, max_gap=flanking_locus_gap)
                    if not parent_loci:
                        continue

                    # Only process the single strongest off-chrom locus
                    best_locus = max(
                        parent_loci,
                        key=lambda hs: (max(h.get('bits', 0) for h in hs), len(hs))
                    )
                    best_bits = max(h.get('bits', 0) for h in best_locus)
                    if best_bits < 40:
                        continue  # Too weak to be a real ortholog

                    off_chrom = best_locus[0].get('chrom', '')
                    if off_chrom not in genome_seqs:
                        continue

                    # Build a local region around the off-chrom hits
                    off_start = min(h.get('start', 0) for h in best_locus)
                    off_end = max(h.get('end', 0) for h in best_locus)

                    # If this off-chrom hit is a GOI proxy (e.g., gene-LY6E),
                    # run a full GOI block search instead of just annotating it
                    # as rearranged_flanking.
                    if parent_id in goi_proxy_flanking_parents:
                        _proxy_block = {
                            'chrom': off_chrom,
                            'start': off_start,
                            'end': off_end,
                            'genes_count': 1,
                            'loci_count': len(best_locus),
                            'genes': [parent_id],
                        }
                        logger.info(
                            f"[{genome_name}] Cross-chrom GOI proxy {parent_id} on "
                            f"{off_chrom}:{off_start}-{off_end} "
                            f"({best_bits:.0f} bits) — running full GOI search."
                        )
                        _proxy_idx = len(synteny_blocks) + 1000 + len(all_genes)
                        _p_genes, _p_gff = process_region_block(
                            _proxy_idx,
                            _proxy_block,
                            hits,
                            genome_seqs,
                            base_db_sequences,
                            genome_name,
                            args,
                            unique_id,
                            threads_per_job,
                            native_annot_index=native_annot_index,
                        )
                        all_genes.extend(_p_genes)
                        all_gff_lines.extend(_p_gff)
                        found_flanking_parents.add(parent_id)
                        continue

                    off_slen = len(genome_seqs[off_chrom])
                    off_pad = min(args.region_padding, 50000)
                    off_w_start = max(0, off_start - off_pad)
                    off_w_end = min(off_slen, off_end + off_pad)
                    off_subseq = genome_seqs[off_chrom][off_w_start:off_w_end]

                    # Remap hits to local coordinates for annotation
                    local_hits = []
                    for h in best_locus:
                        local_start = max(0, int(h.get('start', 0)) - off_w_start)
                        local_end = min(len(off_subseq), int(h.get('end', 0)) - off_w_start)
                        if local_end <= local_start:
                            continue
                        local_hits.append({
                            'query': h.get('query', ''),
                            'qstart': h.get('qstart', 1),
                            'qend': h.get('qend', 1),
                            'gstart': local_start,
                            'gend': local_end,
                            'evalue': h.get('evalue', 1),
                            'pident': h.get('pident', 0),
                            'alnlen': h.get('alnlen', 0),
                            'bits': h.get('bits', 0),
                            'strand': h.get('strand', '+'),
                            'chrom': off_chrom,
                        })

                    if not local_hits:
                        continue

                    # Attempt miniprot-style exon annotation
                    exons = []
                    try:
                        with maybe_quiet_streams(args.quiet_subtools):
                            exons, _ = annotate_exons_from_hit_list(
                                local_hits,
                                parent_query_seq,
                                off_subseq,
                                off_chrom,
                                search_missing=True,
                                gap_min_size=args.gap_min_size,
                                gap_search_window=args.gap_search_window,
                                gap_evalue=args.gap_evalue,
                                gap_min_identity=args.gap_min_identity,
                                gap_min_alnlen=args.gap_min_alnlen,
                                gap_max_hits=args.gap_max_hits,
                                exon_query_mode=False,
                                min_exon_query_cov=args.min_exon_query_cov,
                                min_exon_alnlen=args.min_exon_alnlen,
                                sensitive=False
                            )
                    except Exception as e:
                        logger.warning(f"[{genome_name}] Off-block exon annotation failed at {off_chrom}:{off_w_start}-{off_w_end}: {e}")
                        exons = []

                    if exons:
                        exons.sort(key=lambda e: e.get('qstart', 0))
                        exon_protein = ''.join(e['seq'] for e in exons)
                        strand = exons[0].get('strand', '+')
                        avg_pident = sum(e.get('pident', 0) for e in exons) / len(exons)
                        model_qcov = None
                        if parent_query_seq:
                            qmin_model = min(e.get('qstart', 1) for e in exons)
                            qmax_model = max(e.get('qend', 1) for e in exons)
                            model_qcov = (qmax_model - qmin_model + 1) / max(1, len(parent_query_seq))
                        global_start = off_w_start + min(e['gstart'] for e in exons) + 1
                        global_end = off_w_start + max(e['gend'] for e in exons)
                        new_id = f"{parent_id}|{clean_gname}_rearranged"

                        rearr_gff = [
                            f"{off_chrom}\trearranged_flanking\tmRNA\t{global_start}\t{global_end}\t"
                            f"{avg_pident:.1f}\t{strand}\t.\t"
                            f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'rearranged_flanking', 'Rearranged_from': ','.join(block_chroms)}, evidence_type='rearranged_flanking', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
                        ]
                        for eidx, e in enumerate(exons, 1):
                            exon_gs = off_w_start + e['gstart'] + 1
                            exon_ge = off_w_start + e['gend']
                            rearr_gff.append(
                                f"{off_chrom}\trearranged_flanking\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\t"
                                f"ID={new_id}_CDS{eidx};Parent={new_id}"
                            )

                        all_genes.append({
                            'id': new_id,
                            'seq': exon_protein,
                            'description': (
                                f"coords:{off_chrom}:{global_start}-{global_end} "
                                f"parent:{parent_id} type:rearranged_flanking "
                                f"identity:{avg_pident:.1f}"
                            )
                        })
                        all_gff_lines.extend(rearr_gff)
                        found_flanking_parents.add(parent_id)
                        logger.info(
                            f"[{genome_name}] Recovered rearranged flanking gene "
                            f"{parent_id} on {off_chrom} ({avg_pident:.1f}% identity, "
                            f"{len(exons)} exons)."
                        )
                    else:
                        # Fallback: use concatenated hit-derived protein
                        strand_votes = {'+': 0, '-': 0}
                        for h in local_hits:
                            strand_votes[h.get('strand', '+')] += 1
                        strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

                        ordered_hits = sorted(
                            [h for h in local_hits if h.get('strand', '+') == strand],
                            key=lambda h: h.get('gstart', 0)
                        )
                        ordered_hits = _longest_monotonic_query_chain(ordered_hits, strand) if ordered_hits else []

                        if ordered_hits:
                            coding_frags = []
                            cds_intervals = []
                            for h in ordered_hits:
                                hs = h.get('gstart', 0)
                                he = h.get('gend', 0)
                                if he <= hs:
                                    continue
                                dna = off_subseq[hs:he]
                                if strand == '-':
                                    dna = reverse_complement(dna)
                                dna = dna[:len(dna) - (len(dna) % 3)]
                                if len(dna) < 9:
                                    continue
                                aa = translate(dna).replace('*', '')
                                if not aa:
                                    continue
                                coding_frags.append(aa)
                                cds_intervals.append((hs, he))

                            if coding_frags and cds_intervals:
                                if strand == '-':
                                    coding_frags = list(reversed(coding_frags))
                                flank_protein = ''.join(coding_frags)
                                global_start = off_w_start + min(s for s, _ in cds_intervals) + 1
                                global_end = off_w_start + max(e for _, e in cds_intervals)
                                avg_pident = sum(h.get('pident', 0) for h in ordered_hits) / len(ordered_hits)
                                qmin = min(min(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                qmax = max(max(h.get('qstart', 0), h.get('qend', 0)) for h in ordered_hits)
                                model_qcov = ((qmax - qmin + 1) / len(parent_query_seq)) if parent_query_seq else None
                                new_id = f"{parent_id}|{clean_gname}_rearranged_fallback"

                                rearr_gff = [
                                    f"{off_chrom}\trearranged_flanking\tmRNA\t{global_start}\t{global_end}\t"
                                    f"{avg_pident:.1f}\t{strand}\t.\t"
                                    f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynVoy_Parent': parent_id, 'Type': 'rearranged_flanking_fallback', 'Rearranged_from': ','.join(block_chroms)}, evidence_type='rearranged_flanking_fallback', identity=avg_pident, exon_count=len(cds_intervals), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
                                ]
                                for cidx, (hs, he) in enumerate(cds_intervals, 1):
                                    exon_gs = off_w_start + hs + 1
                                    exon_ge = off_w_start + he
                                    rearr_gff.append(
                                        f"{off_chrom}\trearranged_flanking\tCDS\t{exon_gs}\t{exon_ge}\t.\t{strand}\t0\t"
                                        f"ID={new_id}_CDS{cidx};Parent={new_id}"
                                    )

                                all_genes.append({
                                    'id': new_id,
                                    'seq': flank_protein,
                                    'description': (
                                        f"coords:{off_chrom}:{global_start}-{global_end} "
                                        f"parent:{parent_id} type:rearranged_flanking_fallback "
                                        f"identity:{avg_pident:.1f}"
                                    )
                                })
                                all_gff_lines.extend(rearr_gff)
                                found_flanking_parents.add(parent_id)
                                logger.info(
                                    f"[{genome_name}] Recovered rearranged flanking gene "
                                    f"{parent_id} on {off_chrom} (fallback, {avg_pident:.1f}% identity)."
                                )

        # --- GOI proxy sweep on block chromosomes (Step 3b) ---
        # Some GOI proxies (e.g., gene-LY6E) may have hits on block chromosomes
        # that fall OUTSIDE the padded regions already processed.  This happens
        # when the proxy hit is on the same large scaffold as another flanking
        # gene (e.g., CYP11B on the same NW_ contig) but far from that anchor.
        # We sweep all such loci and run a full GOI block search on each one.
        if goi_proxy_flanking_parents:
            _proxy_gap = max(5000, min(100000, args.max_intron * 2))
            # Build rough coverage set from already-processed blocks
            _covered: list = []
            for _b in synteny_blocks:
                _pad = args.region_padding
                _covered.append((_b['chrom'], _b['start'] - _pad, _b['end'] + _pad))

            for _proxy_id in goi_proxy_flanking_parents:
                if _proxy_id in found_flanking_parents:
                    continue
                _proxy_on_block = [
                    h for h in hits
                    if extract_base_gene_id(h.get('query', '')) == _proxy_id
                    and h['chrom'] in block_chroms
                ]
                if not _proxy_on_block:
                    continue
                _proxy_loci = split_hits_into_loci(_proxy_on_block, max_gap=_proxy_gap)
                for _lh in _proxy_loci:
                    _lbits = max(h.get('bits', 0) for h in _lh)
                    if _lbits < 40:
                        continue
                    _lchrom = _lh[0]['chrom']
                    _lstart = min(h.get('start', 0) for h in _lh)
                    _lend = max(h.get('end', 0) for h in _lh)
                    # Skip if already covered by a processed block
                    if any(
                        _cc == _lchrom and _lstart < _ce and _lend > _cs
                        for _cc, _cs, _ce in _covered
                    ):
                        continue
                    _sweep_block = {
                        'chrom': _lchrom,
                        'start': _lstart,
                        'end': _lend,
                        'genes_count': 1,
                        'loci_count': len(_lh),
                        'genes': [_proxy_id],
                    }
                    logger.info(
                        f"[{genome_name}] On-block GOI proxy sweep: {_proxy_id} on "
                        f"{_lchrom}:{_lstart}-{_lend} ({_lbits:.0f} bits) outside "
                        f"processed regions — running full GOI search."
                    )
                    _sw_idx = len(synteny_blocks) + 2000 + len(all_genes)
                    _sw_genes, _sw_gff = process_region_block(
                        _sw_idx,
                        _sweep_block,
                        hits,
                        genome_seqs,
                        base_db_sequences,
                        genome_name,
                        args,
                        unique_id,
                        threads_per_job,
                        native_annot_index=native_annot_index,
                    )
                    all_genes.extend(_sw_genes)
                    all_gff_lines.extend(_sw_gff)
                    if any(is_goi_query_id(_g.get('id', '')) for _g in _sw_genes):
                        found_flanking_parents.add(_proxy_id)
                    # Add to covered so we don't re-run overlapping loci
                    _covered.append((_lchrom, _lstart - args.region_padding,
                                     _lend + args.region_padding))

        # --- Flanking gene tracking summary (Step 4) ---
        still_missing = all_db_flanking_parents - found_flanking_parents
        if all_db_flanking_parents:
            found_count = len(found_flanking_parents & all_db_flanking_parents)
            total_count = len(all_db_flanking_parents)
            logger.info(
                f"[{genome_name}] Flanking gene summary: "
                f"{found_count}/{total_count} parents found."
            )
            if still_missing:
                # For each missing parent, report best hit info from raw hits
                missing_details = []
                for parent_id in sorted(still_missing):
                    parent_hits = [
                        h for h in hits
                        if extract_base_gene_id(h.get('query', '')) == parent_id
                    ]
                    if parent_hits:
                        best = max(parent_hits, key=lambda h: h.get('bits', 0))
                        missing_details.append(
                            f"{parent_id}(best_hit:{best['chrom']},"
                            f"bits={best.get('bits', 0):.1f},"
                            f"id={best.get('pident', 0):.1f}%)"
                        )
                    else:
                        missing_details.append(f"{parent_id}(no_hits)")
                logger.info(
                    f"[{genome_name}] Missing flanking parents: "
                    f"{'; '.join(missing_details)}"
                )

        feature_meta = {}

        # Write aggregated GFF
        if all_gff_lines:
             gff_out = f"{args.output_dir}/regions/{genome_name}.gff"
             faa_out = f"{args.output_dir}/regions/{genome_name}.faa"
             tsv_out = f"{args.output_dir}/regions/{genome_name}.homology.tsv"
             
             with open(gff_out, 'w') as gf:
                 gf.write("##gff-version 3\n")
                 for gl in all_gff_lines:
                     gf.write(gl + "\n")

             write_fasta([(g['id'], g['seq']) for g in all_genes], faa_out)

             for gl in all_gff_lines:
                 parts = gl.split("\t")
                 if len(parts) < 9 or parts[2] not in {"mRNA", "gene"}:
                     continue
                 attrs = _parse_gff_attributes(parts[8])
                 model_id = attrs.get("ID")
                 if not model_id:
                     continue
                 feature_meta[model_id] = {
                     "parent": _select_parent_id(attrs, model_id),
                     "role": attrs.get("SynVoyRole", "goi" if is_goi_query_id(model_id) else "flanking"),
                     "confidence": attrs.get("Confidence", ""),
                     "goi_class": attrs.get("GOIClass", ""),
                     "model_status": attrs.get("ModelStatus", ""),
                     "evidence_type": attrs.get("EvidenceType", attrs.get("Type", "")),
                     "identity": attrs.get("Identity", ""),
                     "n_exons": attrs.get("Exons", ""),
                     "synteny_context": attrs.get("SyntenyContext", ""),
                     "block_flanking_support": attrs.get("BlockFlankingSupport", ""),
                     "query_coverage": attrs.get("QueryCoverage", ""),
                     "target_gene": attrs.get("TargetGene", ""),
                     "target_product": attrs.get("TargetProduct", ""),
                     "embedding_similarity": attrs.get("EmbeddingSimilarity", ""),
                     "structural_similarity": attrs.get("StructuralSimilarity", ""),
                 }

             with open(tsv_out, 'w') as tf:
                 tf.write(
                     "target_id\thome_id\trole\tconfidence\tgoi_class\tmodel_status\t"
                     "evidence_type\tidentity\tn_exons\tsynteny_context\t"
                     "block_flanking_support\tquery_coverage\ttarget_gene\ttarget_product\t"
                     "embedding_similarity\tstructural_similarity\n"
                 )
                 for rec in all_genes:
                     meta = feature_meta.get(rec['id'], {})
                     parent = meta.get("parent") or extract_base_gene_id(rec['id'])
                     tf.write(
                         "\t".join([
                             rec['id'],
                             parent,
                             meta.get("role", "goi" if is_goi_query_id(rec['id']) else "flanking"),
                             meta.get("confidence", ""),
                             meta.get("goi_class", ""),
                             meta.get("model_status", ""),
                             meta.get("evidence_type", ""),
                             meta.get("identity", ""),
                             meta.get("n_exons", ""),
                             meta.get("synteny_context", ""),
                             meta.get("block_flanking_support", ""),
                             meta.get("query_coverage", ""),
                             meta.get("target_gene", ""),
                             meta.get("target_product", ""),
                             meta.get("embedding_similarity", ""),
                             meta.get("structural_similarity", ""),
                         ]) + "\n"
                     )
        
        # Expansion DB should only be augmented with GOI-derived models that
        # survived confidence/ambiguity triage. Low-confidence ambiguous/tandem
        # calls are still reported in GFF/plots but should not recursively seed
        # later waves.
        new_genes, tree_extra_genes, suppressed_seed_count = (
            _classify_goi_for_seed_and_tree(all_genes, feature_meta, is_goi_query_id)
        )
        if all_genes:
            logger.info(
                f"[{genome_name}] Expansion payload: {len(new_genes)} GOI-derived / "
                f"{len(all_genes)} total annotations"
                + (
                    f" ({suppressed_seed_count} GOI-like annotations withheld from seeding"
                    + (f", {len(tree_extra_genes)} of which still feed the tree"
                       if tree_extra_genes else "")
                    + ")."
                    if suppressed_seed_count
                    else "."
                )
            )

    except Exception as e:
        error_message = str(e)
        logger.error(f"[{genome_name}] Error processing: {error_message}")
    finally:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)

    return genome_name, new_genes, tree_extra_genes, error_message
def main():
    parser = argparse.ArgumentParser(description="Iterative Genome Search Runner (Wavefront Parallel)")
    parser.add_argument("--initial_db", required=True)
    parser.add_argument("--sorted_genomes", required=True, 
                        help="Tab-separated file: genome_path\\tdistance. Genomes sorted by distance.")
    parser.add_argument("--genomes_dir", help="Directory containing genome files (if paths in sorted_genomes are relative)")
    parser.add_argument("--home_db_dir", help="Home Proteome MMseqs DB for RBH")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--min_identity", type=float, default=40.0)
    parser.add_argument("--min_length", type=int, default=50)
    parser.add_argument("--max_intron", type=int, default=20000)
    parser.add_argument("--threads", type=int, default=4, help="Total threads available for parallel processing")
    parser.add_argument("--cluster_distance", type=int, default=-1, help="Auto-detect if -1")
    parser.add_argument("--mmseqs_sens", type=float, default=7.5, help="MMseqs2 sensitivity (higher = more sensitive but slower)")
    parser.add_argument("--mmseqs_split_memory_limit", default="0", help="MMseqs split memory limit (e.g. 3G, 8000M, 0=auto)")
    parser.add_argument("--mmseqs_verbosity", type=int, default=1, help="MMseqs verbosity (0-3)")
    parser.add_argument("--min_gene_identity", type=float, default=25.0, help="Minimum identity for RBH validation")
    parser.add_argument("--region_padding", type=int, default=100000, help="Default padding around synteny block")
    parser.add_argument("--padding_min", type=int, default=50000, help="Minimum adaptive padding")
    parser.add_argument("--padding_max", type=int, default=200000, help="Maximum adaptive padding")
    parser.add_argument("--enable_smith_waterman", type=str2bool, default=True, help="Enable Smith-Waterman search")
    parser.add_argument("--sw_method", type=str, default="auto", help="Smith-Waterman backend (auto, parasail, ssearch36)")
    parser.add_argument("--sw_min_score", type=float, default=50.0)
    parser.add_argument("--sw_min_identity", type=float, default=20.0)
    parser.add_argument("--sw_timeout_seconds", type=int, default=300)
    parser.add_argument("--aug_relaxed_evalue_mult", type=float, default=1000.0)
    parser.add_argument("--aug_relaxed_evalue_cap", type=float, default=10.0)
    parser.add_argument("--aug_relaxed_parse_evalue_mult", type=float, default=10.0)
    parser.add_argument("--aug_relaxed_identity_factor", type=float, default=0.6)
    parser.add_argument("--aug_relaxed_identity_min", type=float, default=25.0)
    parser.add_argument("--aug_relaxed_length_div", type=float, default=2.0)
    parser.add_argument("--aug_relaxed_length_min", type=int, default=15)
    parser.add_argument("--aug_dedup_bin_bp", type=int, default=100)
    parser.add_argument("--gap_search_window", type=int, default=50000)
    parser.add_argument("--gap_min_size", type=int, default=10)
    parser.add_argument("--gap_evalue", type=float, default=10.0)
    parser.add_argument("--gap_min_identity", type=float, default=25.0)
    parser.add_argument("--gap_min_alnlen", type=int, default=10)
    parser.add_argument("--gap_max_hits", type=int, default=5)
    parser.add_argument("--min_exon_query_cov", type=float, default=0.25,
                        help="Minimum query coverage for exon hits")
    parser.add_argument("--min_exon_alnlen", type=int, default=30,
                        help="Minimum alignment length for exon hits")
    parser.add_argument(
        "--max_blocks_per_genome",
        type=int,
        default=80,
        help="Maximum synteny blocks to process per genome (0 = no cap)"
    )
    parser.add_argument(
        "--min_block_genes",
        type=int,
        default=2,
        help="Minimum unique anchor genes required for a synteny block"
    )
    parser.add_argument(
        "--max_consecutive_empty_blocks",
        type=int,
        default=25,
        help="Stop processing more blocks after this many GOI-empty blocks in a row (0 = disable)"
    )
    parser.add_argument(
        "--quiet_subtools",
        type=str2bool,
        default=True,
        help="Suppress noisy third-party tool output (miniprot/tblastn diagnostics)"
    )
    parser.add_argument("--prefix", default="", help="Prefix for output files (e.g. locus ID)")
    parser.add_argument("--resume", action="store_true", help="Resume from previous checkpoint if available")

    # Classification thresholds (GOI confidence / model status)
    parser.add_argument("--classify_high_min_identity", type=float, default=50.0,
                        help="Min identity for HIGH confidence exon_annotation "
                             "(lowered from 60 to accommodate cross-vertebrate orthologs)")
    parser.add_argument("--classify_medium_min_identity", type=float, default=35.0,
                        help="Min identity for MEDIUM confidence exon_annotation "
                             "(lowered from 45 to accommodate divergent orthologs)")
    parser.add_argument("--classify_tandem_min_identity", type=float, default=40.0,
                        help="Min identity for MEDIUM confidence tandem_copy (below = LOW)")
    parser.add_argument("--classify_fragment_max_qcov", type=float, default=0.4,
                        help="Query coverage below this marks model as fragment")
    parser.add_argument("--classify_complete_min_qcov", type=float, default=0.7,
                        help="Query coverage above this (with multi-exon) marks model as complete")
    parser.add_argument("--strict_goi_family", type=str2bool, default=False,
                        help="Downgrade fallback/rescued_exon/raw_hit GOI calls whose annotated "
                             "TargetGene/TargetProduct does not match the family tokens")
    parser.add_argument("--goi_family_tokens", type=str, default="",
                        help="Comma-separated family name tokens for strict mode "
                             "(default: auto-derive from query FASTA header GN=...)")

    # PLM (Protein Language Model) embedding search
    parser.add_argument("--enable_plm_search", type=str2bool, default=False,
                        help="Enable ProtT5 embedding search for remote homolog detection")
    parser.add_argument("--plm_device", type=str, default="cpu",
                        help="Device for PLM inference (cpu or cuda)")
    parser.add_argument("--plm_similarity_threshold", type=float, default=0.5,
                        help="Minimum cosine similarity for PLM-discovered ORFs (0-1)")
    parser.add_argument("--plm_medium_threshold", type=float, default=0.7,
                        help="Embedding similarity above this can boost LOW → MEDIUM")
    parser.add_argument("--plm_high_threshold", type=float, default=0.85,
                        help="Embedding similarity above this can boost MEDIUM → HIGH")

    # Gene predictor selection (Augustus for eukaryotes, Prodigal for prokaryotes)
    parser.add_argument("--gene_predictor", type=str, default="auto",
                        help="Gene predictor: auto (prefer Augustus), augustus, or prodigal")
    parser.add_argument("--augustus_species", type=str, default="fly",
                        help="Augustus species model (e.g. fly, human, honeybee1)")

    # Structural search (ESMFold + Foldseek)
    parser.add_argument("--enable_structural_search", type=str2bool, default=False,
                        help="Enable ESMFold + Foldseek structural search for remote homologs")
    parser.add_argument("--structural_device", type=str, default="cpu",
                        help="Device for ESMFold inference (cpu or cuda)")
    parser.add_argument("--structural_tm_threshold", type=float, default=0.3,
                        help="Minimum TM-score for Foldseek-discovered ORFs (0-1)")
    parser.add_argument("--structural_medium_threshold", type=float, default=0.5,
                        help="TM-score above this can boost LOW → MEDIUM")
    parser.add_argument("--structural_high_threshold", type=float, default=0.7,
                        help="TM-score above this can boost MEDIUM → HIGH")
    parser.add_argument("--structural_max_length", type=int, default=700,
                        help="Max sequence length for ESMFold (VRAM safety, default 700)")

    args = parser.parse_args()

    # Populate classification thresholds from CLI args
    CLASSIFY_THRESHOLDS["high_min_identity"] = args.classify_high_min_identity
    CLASSIFY_THRESHOLDS["medium_min_identity"] = args.classify_medium_min_identity
    CLASSIFY_THRESHOLDS["tandem_min_identity"] = args.classify_tandem_min_identity
    CLASSIFY_THRESHOLDS["fragment_max_qcov"] = args.classify_fragment_max_qcov
    CLASSIFY_THRESHOLDS["complete_min_qcov"] = args.classify_complete_min_qcov
    CLASSIFY_THRESHOLDS["plm_medium_threshold"] = args.plm_medium_threshold
    CLASSIFY_THRESHOLDS["plm_high_threshold"] = args.plm_high_threshold
    CLASSIFY_THRESHOLDS["structural_medium_threshold"] = args.structural_medium_threshold
    CLASSIFY_THRESHOLDS["structural_high_threshold"] = args.structural_high_threshold

    # Family-consistency config (strict GOI mode)
    FAMILY_CONFIG["strict"] = bool(args.strict_goi_family)
    override_tokens = {
        _normalize_family_token(t) for t in (args.goi_family_tokens or "").split(",")
    }
    override_tokens.discard("")
    if override_tokens:
        FAMILY_CONFIG["tokens"] = override_tokens
    else:
        FAMILY_CONFIG["tokens"] = _auto_derive_family_tokens(args.initial_db)
    if FAMILY_CONFIG["strict"] and not FAMILY_CONFIG["tokens"]:
        logger.warning(
            "strict_goi_family=true but no family tokens derived from query FASTA; "
            "strict downgrade will be disabled. Provide --goi_family_tokens explicitly."
        )
        FAMILY_CONFIG["strict"] = False
    if FAMILY_CONFIG["tokens"]:
        logger.info(
            f"GOI family tokens: {sorted(FAMILY_CONFIG['tokens'])} "
            f"(strict={FAMILY_CONFIG['strict']})"
        )

    # INPUT VALIDATION
    # 1. Validate required files exist
    if not os.path.exists(args.initial_db):
        logger.error(
            f"Initial database file not found: {args.initial_db}. "
            f"This FASTA is produced by EXTRACT_FLANKING_GENES and contains "
            f"the flanking-gene proteins used as synteny parents. "
            f"Try: (1) rerun without -resume to regenerate it; "
            f"(2) check the Nextflow work dir for a failed EXTRACT_FLANKING_GENES "
            f"process; (3) ensure the home GFF has protein-coding genes "
            f"flanking the GOI."
        )
        sys.exit(1)

    if not os.path.exists(args.sorted_genomes):
        logger.error(
            f"Sorted genomes file not found: {args.sorted_genomes}. "
            f"This is the ranked target-genome list produced by "
            f"RANK_BY_SIMILARITY. "
            f"Try: rerun without -resume; if it still fails, check the "
            f"upstream RANK_BY_SIMILARITY process logs in Nextflow's work/ dir."
        )
        sys.exit(1)

    # 2. Validate initial_db is not empty
    if os.path.getsize(args.initial_db) == 0:
        logger.error(
            f"Initial database file is empty: {args.initial_db}. "
            f"EXTRACT_FLANKING_GENES produced a zero-byte FASTA, which means "
            f"no flanking genes were found. "
            f"Try: (1) verify the home GFF has 'gene' and 'CDS' features "
            f"around the GOI coordinates; "
            f"(2) increase --n_flanking_genes; "
            f"(3) check for GOI-vs-flanking similarity filtering being too "
            f"aggressive (see EXTRACT_FLANKING_GENES logs)."
        )
        sys.exit(1)
    
    # 3. Validate parameters are in valid ranges
    if args.min_identity < 0 or args.min_identity > 100:
        logger.error(f"Invalid min_identity: {args.min_identity}. Must be between 0 and 100")
        sys.exit(1)
    
    if args.min_length < 1:
        logger.error(f"Invalid min_length: {args.min_length}. Must be >= 1")
        sys.exit(1)
    
    if args.evalue <= 0:
        logger.error(f"Invalid evalue: {args.evalue}. Must be > 0")
        sys.exit(1)
    
    if args.threads < 1:
        logger.error(f"Invalid threads: {args.threads}. Must be >= 1")
        sys.exit(1)
    
    if args.mmseqs_sens < 1 or args.mmseqs_sens > 9:
        logger.warning(f"MMseqs sensitivity {args.mmseqs_sens} outside typical range (1-9)")
    if args.mmseqs_verbosity < 0 or args.mmseqs_verbosity > 3:
        logger.error("mmseqs_verbosity must be between 0 and 3")
        sys.exit(1)
    if args.padding_min < 0 or args.padding_max < 0 or args.region_padding < 0:
        logger.error("Padding values must be non-negative")
        sys.exit(1)
    if args.padding_max < args.padding_min:
        logger.error("padding_max must be >= padding_min")
        sys.exit(1)
    if args.max_intron < 0:
        logger.error("max_intron must be >= 0")
        sys.exit(1)
    if args.gap_min_size < 1:
        logger.error("gap_min_size must be >= 1")
        sys.exit(1)
    if args.gap_max_hits < 1:
        logger.error("gap_max_hits must be >= 1")
        sys.exit(1)
    if args.min_exon_query_cov < 0 or args.min_exon_query_cov > 1:
        logger.error("min_exon_query_cov must be between 0 and 1")
        sys.exit(1)
    if args.min_exon_alnlen < 1:
        logger.error("min_exon_alnlen must be >= 1")
        sys.exit(1)
    if args.max_blocks_per_genome < 0:
        logger.error("max_blocks_per_genome must be >= 0")
        sys.exit(1)
    if args.min_block_genes < 1:
        logger.error("min_block_genes must be >= 1")
        sys.exit(1)
    if args.max_consecutive_empty_blocks < 0:
        logger.error("max_consecutive_empty_blocks must be >= 0")
        sys.exit(1)
    if args.aug_relaxed_length_div <= 0:
        logger.error("aug_relaxed_length_div must be > 0")
        sys.exit(1)

    # Classification threshold validation
    for name in ['classify_high_min_identity', 'classify_medium_min_identity', 'classify_tandem_min_identity']:
        val = getattr(args, name)
        if val < 0 or val > 100:
            logger.error(f"{name} must be between 0 and 100 (got {val})")
            sys.exit(1)
    for name in ['classify_fragment_max_qcov', 'classify_complete_min_qcov']:
        val = getattr(args, name)
        if val < 0 or val > 1:
            logger.error(f"{name} must be between 0 and 1 (got {val})")
            sys.exit(1)
    if args.classify_fragment_max_qcov >= args.classify_complete_min_qcov:
        logger.warning(
            f"classify_fragment_max_qcov ({args.classify_fragment_max_qcov}) >= "
            f"classify_complete_min_qcov ({args.classify_complete_min_qcov}); "
            f"fragment and complete model status ranges overlap"
        )

    # Smith-Waterman validation
    if args.sw_min_score < 0:
        logger.error(f"sw_min_score must be >= 0 (got {args.sw_min_score})")
        sys.exit(1)
    if args.sw_min_identity < 0 or args.sw_min_identity > 100:
        logger.error(f"sw_min_identity must be between 0 and 100 (got {args.sw_min_identity})")
        sys.exit(1)
    if args.sw_timeout_seconds < 1:
        logger.error(f"sw_timeout_seconds must be >= 1 (got {args.sw_timeout_seconds})")
        sys.exit(1)

    # In auto mode, only keep SW enabled when parasail is available.
    # This avoids per-block fallback overhead on very large block sets.
    if args.enable_smith_waterman and args.sw_method == "auto" and not has_parasail_available():
        logger.warning(
            "parasail not available in auto mode; disabling Smith-Waterman for this run. "
            "Set --sw_method ssearch36 to force the slower fallback."
        )
        args.enable_smith_waterman = False

    # Gene predictor validation
    if args.gene_predictor not in ("auto", "augustus", "prodigal"):
        logger.error(f"gene_predictor must be 'auto', 'augustus', or 'prodigal' (got '{args.gene_predictor}')")
        sys.exit(1)

    # PLM search validation & availability check
    if args.enable_plm_search:
        if not PLM_IMPORT_OK:
            logger.error(
                "PLM search enabled but plm_search module could not be imported.\n"
                "  Ensure bin/plm_search.py exists and torch/transformers are installed:\n"
                "  pip install torch transformers sentencepiece"
            )
            sys.exit(1)
        if not check_plm_available():
            logger.error(
                "PLM search enabled but PyTorch/Transformers not installed.\n"
                "  Install with: pip install torch transformers sentencepiece"
            )
            sys.exit(1)
        if args.plm_similarity_threshold < 0 or args.plm_similarity_threshold > 1:
            logger.error("plm_similarity_threshold must be between 0 and 1")
            sys.exit(1)
        if args.plm_device not in ("cpu", "cuda"):
            logger.error("plm_device must be 'cpu' or 'cuda'")
            sys.exit(1)

    # Structural search (ESMFold + Foldseek) validation & availability check
    if args.enable_structural_search:
        if not STRUCTURAL_IMPORT_OK:
            logger.error(
                "Structural search enabled but structural_search module could not be imported.\n"
                "  Ensure bin/structural_search.py exists and dependencies are installed:\n"
                "  pip install torch transformers\n"
                "  conda install -c bioconda foldseek"
            )
            sys.exit(1)
        if not check_esmfold_available():
            logger.error(
                "Structural search enabled but ESMFold not available.\n"
                "  Install with: pip install torch transformers"
            )
            sys.exit(1)
        if not check_foldseek_available():
            logger.error(
                "Structural search enabled but Foldseek binary not found.\n"
                "  Install with: conda install -c bioconda foldseek"
            )
            sys.exit(1)
        if args.structural_tm_threshold < 0 or args.structural_tm_threshold > 1:
            logger.error("structural_tm_threshold must be between 0 and 1")
            sys.exit(1)
        if args.structural_device not in ("cpu", "cuda"):
            logger.error("structural_device must be 'cpu' or 'cuda'")
            sys.exit(1)
        if args.structural_max_length < 10:
            logger.error("structural_max_length must be >= 10")
            sys.exit(1)

    logger.info("")
    logger.info("═══ Iterative Search Configuration ═══")
    logger.info(f"  Threads:          {args.threads}")
    logger.info(f"  Min identity:     {args.min_identity}%")
    logger.info(f"  Min length:       {args.min_length}")
    logger.info(f"  E-value cutoff:   {args.evalue}")
    logger.info(f"  MMseqs sens:      {args.mmseqs_sens}")
    logger.info(f"  MMseqs mem limit: {args.mmseqs_split_memory_limit}")
    logger.info(f"  Gene predictor:   {args.gene_predictor}" + (f" (species={args.augustus_species})" if args.gene_predictor != "prodigal" else ""))
    logger.info(f"  Smith-Waterman:   {'enabled (' + args.sw_method + ')' if args.enable_smith_waterman else 'disabled'}")
    if args.enable_plm_search:
        logger.info(f"  PLM search:       enabled (ProtT5, device={args.plm_device})")
        logger.info(f"  PLM threshold:    {args.plm_similarity_threshold}")
        logger.info(f"  PLM boost:        MEDIUM≥{args.plm_medium_threshold}, HIGH≥{args.plm_high_threshold}")
    else:
        logger.info(f"  PLM search:       disabled")
    if args.enable_structural_search:
        logger.info(f"  Structural search: enabled (ESMFold+Foldseek, device={args.structural_device})")
        logger.info(f"  Struct TM thresh: {args.structural_tm_threshold}")
        logger.info(f"  Struct boost:     MEDIUM≥{args.structural_medium_threshold}, HIGH≥{args.structural_high_threshold}")
        logger.info(f"  Struct max len:   {args.structural_max_length}")
    else:
        logger.info(f"  Structural search: disabled")
    logger.info(f"  Max blocks/genome:{args.max_blocks_per_genome}")
    logger.info(f"  Region padding:   {args.padding_min}-{args.padding_max} bp")
    if args.prefix:
        logger.info(f"  Locus prefix:     {args.prefix}")
    logger.info("══════════════════════════════════════")
    
    prefix = f"{args.prefix}_" if args.prefix else ""
    
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(f"{args.output_dir}/hits", exist_ok=True)
    os.makedirs(f"{args.output_dir}/regions", exist_ok=True)
    
    # CHECKPOINTING: Check for resume
    checkpoint_file = f"{args.output_dir}/.checkpoint"
    start_wave = 0
    current_db = f"{args.output_dir}/current_db.faa"
    
    if args.resume and os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as cf:
            checkpoint_data = json.loads(cf.read())
            start_wave = checkpoint_data.get('completed_waves', 0)
            total_saved_waves = checkpoint_data.get('total_waves', '?')
            last_db = checkpoint_data.get('last_db', None)

            if last_db and os.path.exists(last_db):
                logger.info(
                    f"Resuming from checkpoint: {start_wave}/{total_saved_waves} waves already completed. "
                    f"Continuing from wave {start_wave + 1}."
                )
                current_db = last_db
            else:
                logger.warning(
                    f"Checkpoint found ({start_wave} waves completed) but database file is missing "
                    f"(expected: {last_db}). Starting from beginning."
                )
                start_wave = 0
                shutil.copyfile(args.initial_db, current_db)
    else:
        if args.resume:
            logger.info("No checkpoint found; starting fresh.")
        shutil.copyfile(args.initial_db, current_db)
    
    # Parse Genomes and Distances
    genome_entries = []
    parse_warnings = []
    with open(args.sorted_genomes, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            gname = parts[0]
            if len(parts) < 2:
                parse_warnings.append(
                    f"  Line {line_num}: missing distance column (using 0.0): {line!r}"
                )
            try:
                dist = float(parts[1]) if len(parts) > 1 else 0.0
            except ValueError:
                logger.error(
                    f"sorted_genomes line {line_num}: cannot parse distance "
                    f"from {parts[1]!r}. Expected tab-separated: genome_path<TAB>distance"
                )
                sys.exit(1)

            gpath = gname
            if args.genomes_dir:
                if not os.path.isabs(gname):
                    gpath = os.path.join(args.genomes_dir, os.path.basename(gname))

            genome_entries.append({'name': gname, 'path': gpath, 'dist': dist})

    if parse_warnings:
        for w in parse_warnings:
            logger.warning(w)

    # Validate we loaded genomes
    if not genome_entries:
        logger.error(
            f"No genomes found in sorted_genomes file: {args.sorted_genomes}\n"
            f"  The file should contain one genome per line: genome_path<TAB>distance\n"
            f"  Check that the upstream PHYLO_SORT step completed successfully."
        )
        sys.exit(1)

    # Pre-flight: verify all genome files are accessible before starting the search
    missing_genomes = []
    for entry in genome_entries:
        if not os.path.exists(entry['path']):
            missing_genomes.append(f"  {entry['name']} -> {entry['path']}")
    if missing_genomes:
        logger.error(
            f"{len(missing_genomes)} of {len(genome_entries)} genome file(s) not found:\n"
            + "\n".join(missing_genomes[:10])
            + (f"\n  ... and {len(missing_genomes) - 10} more" if len(missing_genomes) > 10 else "")
            + f"\n  genomes_dir={args.genomes_dir or '(not set)'}"
            + "\n  Check that target genomes were downloaded/staged correctly."
        )
        sys.exit(1)

    logger.info(f"Loaded {len(genome_entries)} genomes (all files accessible).")

    # Normalize phylogenetic distances if they are not in [0,1]
    finite_dists = [g['dist'] for g in genome_entries if g.get('dist') not in [None, float('inf')]]
    if finite_dists:
        max_dist = max(finite_dists)
        if max_dist > 1.0:
            for g in genome_entries:
                if g['dist'] == float('inf'):
                    g['dist'] = 1.0
                else:
                    g['dist'] = g['dist'] / max_dist
            logger.info(f"Normalized phylogenetic distances by max {max_dist:.3f}")
    
    # Define Waves
    waves = []

    # IMPROVED WAVEFRONT STRATEGY:
    # - Closest genomes (dist < 0.05): Process strictly serially for maximum sensitivity
    # - Medium distance (0.05 - 0.15): Small waves (2-3 genomes)
    # - Distant genomes (> 0.15): Larger waves (can parallelize more)
    
    i = 0
    while i < len(genome_entries):
        curr = genome_entries[i]
        
        if curr['dist'] < 0.05:
            # Very close: Serial processing (wave of 1)
            waves.append([curr])
            i += 1
        elif curr['dist'] < 0.15:
            # Medium distance: Small waves of 2-3 genomes with similar distance
            wave = [curr]
            i += 1
            while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.01:
                wave.append(genome_entries[i])
                i += 1
                if len(wave) >= 3:  # Max 3 per wave for medium distance
                    break
            waves.append(wave)
        else:
            # Distant: Can parallelize more (waves of up to 5)
            wave = [curr]
            i += 1
            while i < len(genome_entries) and abs(genome_entries[i]['dist'] - curr['dist']) < 0.02:
                wave.append(genome_entries[i])
                i += 1
                if len(wave) >= 5:  # Max 5 per wave for distant genomes
                    break
            waves.append(wave)
    
    logger.info(f"Defined {len(waves)} waves of execution.")

    # ── PLM: Pre-compute GOI embeddings before waves start ──
    if args.enable_plm_search:
        goi_emb_path = os.path.join(args.output_dir, "goi_embeddings.npz")
        if not os.path.exists(goi_emb_path):
            logger.info("Pre-computing GOI embeddings with ProtT5 ...")
            goi_embs = precompute_goi_embeddings(
                db_fasta=args.initial_db,
                output_path=goi_emb_path,
                device=args.plm_device,
            )
            if not goi_embs:
                logger.warning(
                    "No GOI sequences found for PLM embedding. "
                    "PLM search will be skipped."
                )
                args.enable_plm_search = False
            else:
                logger.info(
                    f"GOI embeddings ready: {len(goi_embs)} sequence(s) → {goi_emb_path}"
                )
        else:
            logger.info(f"Using cached GOI embeddings: {goi_emb_path}")

    # ── Structural: Pre-fold GOI structures before waves start ──
    if args.enable_structural_search:
        goi_struct_index = os.path.join(args.output_dir, "goi_structure_index.tsv")
        if not os.path.exists(goi_struct_index):
            logger.info("Pre-folding GOI structures with ESMFold ...")
            goi_structs = prefold_goi_structures(
                db_fasta=args.initial_db,
                output_dir=args.output_dir,
                device=args.structural_device,
                max_length=args.structural_max_length,
            )
            if not goi_structs:
                logger.warning(
                    "No GOI sequences found for structural folding. "
                    "Structural search will be skipped."
                )
                args.enable_structural_search = False
            else:
                save_structure_index(goi_structs, goi_struct_index)
                logger.info(
                    f"GOI structures ready: {len(goi_structs)} structure(s) → {goi_struct_index}"
                )
        else:
            logger.info(f"Using cached GOI structures: {goi_struct_index}")

    # Compute cumulative genome index for progress tracking
    total_genomes = len(genome_entries)
    cumulative_genome_idx = 0  # genomes completed so far

    latest_db = current_db
    total_new_genes = 0
    genomes_with_hits = 0
    genomes_without_hits = 0
    # Tree-only GOI hits (e.g. tandem_goi_copy) — accumulated across waves
    # and written to goi_for_tree.faa so that COMPUTE_TREE sees species like
    # Apis cerana whose only GOI evidence is a tandem duplicate.
    tree_extras_total = []
    tree_extras_seen_ids = set()

    for i, wave in enumerate(waves):
        # Skip already completed waves
        if i < start_wave:
            cumulative_genome_idx += len(wave)
            logger.info(f"Skipping wave {i+1}/{len(waves)} (already completed, resuming)")
            continue

        wave_genome_names = [e['name'] for e in wave]
        logger.info(
            f"═══ Wave {i+1}/{len(waves)} "
            f"[genomes {cumulative_genome_idx+1}-{cumulative_genome_idx+len(wave)}/{total_genomes}] "
            f"({len(wave)} genome(s), dist≈{wave[0]['dist']:.3f}) ═══"
        )
        for entry in wave:
            gbase = os.path.basename(entry['name'])
            logger.info(f"  • {gbase} (dist={entry['dist']:.4f})")

        # Parallel Execution
        max_workers = min(len(wave), args.threads)
        threads_per_job = max(1, args.threads // max_workers)

        wave_results = []
        wave_errors = []
        wave_start_time = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for entry in wave:
                future = executor.submit(
                    process_single_genome,
                    entry['path'], latest_db, args, args.home_db_dir, prefix, threads_per_job
                )
                futures[future] = entry

            completed = 0
            for future in concurrent.futures.as_completed(futures):
                entry = futures[future]
                completed += 1
                cumulative_genome_idx += 1
                gbase = os.path.basename(entry['name'])
                try:
                    gname, new_genes, tree_extras, genome_error = future.result()
                    if genome_error:
                        wave_errors.append((gname, genome_error))
                        logger.error(
                            f"  ✗ {gbase} FAILED ({completed}/{len(wave)} in wave, "
                            f"{cumulative_genome_idx}/{total_genomes} overall): {genome_error}"
                        )
                        continue
                    gene_count = len(new_genes) if new_genes else 0
                    extra_count = len(tree_extras) if tree_extras else 0
                    if gene_count > 0 or extra_count > 0:
                        genomes_with_hits += 1
                        msg = f"  ✓ {gbase}: {gene_count} new gene(s)"
                        if extra_count:
                            msg += f" + {extra_count} tree-only hit(s)"
                        logger.info(
                            f"{msg} ({cumulative_genome_idx}/{total_genomes} overall)"
                        )
                    else:
                        genomes_without_hits += 1
                        logger.info(
                            f"  – {gbase}: no new genes "
                            f"({cumulative_genome_idx}/{total_genomes} overall)"
                        )
                    if new_genes:
                        wave_results.extend(new_genes)
                    if tree_extras:
                        for g in tree_extras:
                            gid = g.get('id', '')
                            if gid and gid not in tree_extras_seen_ids:
                                tree_extras_seen_ids.add(gid)
                                tree_extras_total.append(g)
                except Exception as exc:
                    wave_errors.append((entry.get('name', 'unknown'), str(exc)))
                    logger.error(
                        f"  ✗ {gbase} CRASHED ({cumulative_genome_idx}/{total_genomes}): {exc}"
                    )

        wave_elapsed = time.time() - wave_start_time

        if wave_errors:
            logger.error(f"Wave {i+1} had {len(wave_errors)} error(s):")
            for gname, msg in wave_errors:
                logger.error(f"  • {os.path.basename(gname)}: {msg}")
            # Classify errors to give actionable advice
            oom_errors = [m for _, m in wave_errors if 'memory' in m.lower() or 'killed' in m.lower() or 'oom' in m.lower()]
            if oom_errors:
                logger.error(
                    "  Hint: some failures look like out-of-memory (OOM). Try:\n"
                    "    --mmseqs_split_memory_limit <lower value>  or\n"
                    "    increase the process memory in nextflow.config"
                )
            logger.error("Aborting iterative search due to per-genome processing errors.")
            sys.exit(1)

        # Update DB after Wave
        total_new_genes += len(wave_results)
        if wave_results:
            logger.info(
                f"Wave {i+1} complete: {len(wave_results)} new gene(s) in {wave_elapsed:.0f}s. "
                f"Running total: {total_new_genes} gene(s) from {genomes_with_hits} genome(s)."
            )

            new_genes_fasta = f"{args.output_dir}/iter_{i+1}_new_genes.faa"
            write_fasta([(g['id'], g['seq']) for g in wave_results], new_genes_fasta)

            next_db = f"{args.output_dir}/db_iter_{i+1}.faa"
            with open(next_db, 'w') as ndb:
                with open(latest_db, 'r') as old_db:
                    shutil.copyfileobj(old_db, ndb)
                with open(new_genes_fasta, 'r') as new_g:
                    shutil.copyfileobj(new_g, ndb)

            # Clean up previous DB if it's not the initial one
            if i > 0 and latest_db != current_db:
                try:
                    os.remove(latest_db)
                except OSError as e:
                    logger.warning(f"Could not remove old DB file {latest_db}: {e}")

            latest_db = next_db
        else:
            logger.info(f"Wave {i+1} complete: no new genes ({wave_elapsed:.0f}s).")

        # CHECKPOINT: Save progress after each wave
        with open(checkpoint_file, 'w') as cf:
            json.dump({
                'completed_waves': i + 1,
                'last_db': latest_db,
                'total_waves': len(waves)
            }, cf)

    expanded_db = f"{args.output_dir}/expanded_db.faa"
    if os.path.exists(latest_db):
        shutil.move(latest_db, expanded_db)

    # Build a tree-input FASTA = expanded_db + tree-only GOI hits (tandem_copy
    # and other MEDIUM/HIGH GOI hits that were withheld from wave seeding).
    # This is what COMPUTE_TREE consumes — keeping seeding strict but the tree
    # complete.
    goi_for_tree = f"{args.output_dir}/goi_for_tree.faa"
    if os.path.exists(expanded_db):
        shutil.copyfile(expanded_db, goi_for_tree)
    else:
        # Defensive: emit an empty file so downstream wiring doesn't trip.
        open(goi_for_tree, 'w').close()
    if tree_extras_total:
        existing_ids = set()
        if os.path.exists(goi_for_tree):
            with open(goi_for_tree) as fh:
                for line in fh:
                    if line.startswith('>'):
                        existing_ids.add(line[1:].split()[0].strip())
        appended = 0
        with open(goi_for_tree, 'a') as fh:
            for g in tree_extras_total:
                gid = g.get('id', '')
                seq = g.get('seq', '')
                if not gid or not seq or gid in existing_ids:
                    continue
                fh.write(f">{gid}\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i+80] + "\n")
                existing_ids.add(gid)
                appended += 1
        logger.info(
            f"  Tree-input FASTA: {goi_for_tree} (added {appended} tree-only hit(s) "
            f"on top of expanded_db.faa)"
        )

    # Final summary
    logger.info("")
    logger.info("═══ Iterative Search Summary ═══")
    logger.info(f"  Genomes searched:  {total_genomes}")
    logger.info(f"  Genomes with hits: {genomes_with_hits}")
    logger.info(f"  Genomes no hits:   {genomes_without_hits}")
    logger.info(f"  Total new genes:   {total_new_genes}")
    logger.info(f"  Final database:    {expanded_db}")
    logger.info("════════════════════════════════")

if __name__ == "__main__":
    main()
