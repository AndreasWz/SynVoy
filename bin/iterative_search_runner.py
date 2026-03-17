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
        
    with open(hits_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            try:
                # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
                # 0      1       2       3       4         5        6       7     8       9     10      11
                if len(parts) < 11: continue
                
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
            except Exception as e:
                continue
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

    Supports SynTerra, NCBI/Ensembl-style, and custom tags.
    """
    for key in [
        "SynTerra_Parent",
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

    if native_annot:
        if native_annot.get("label"):
            merged["TargetGene"] = _safe_gff_value(native_annot["label"])
        if native_annot.get("product"):
            merged["TargetProduct"] = _safe_gff_value(native_annot["product"])
        if native_annot.get("feature_id"):
            merged["TargetID"] = _safe_gff_value(native_annot["feature_id"])

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


def _classify_goi_evidence(
    evidence_type: str,
    identity: float = 0.0,
    exon_count: int = 1,
    query_cov: Optional[float] = None,
    flanking_support: int = 0,
) -> Tuple[str, str, str]:
    """
    Assign a conservative confidence/class label to GOI-derived candidates.

    The goal is not to prove orthology here, but to prevent fallback-heavy
    output from masquerading as confident GOI evidence downstream.
    """
    identity = float(identity or 0.0)
    exon_count = max(1, int(exon_count or 1))
    qcov = float(query_cov or 0.0)
    context = _synteny_context_label(flanking_support)

    if evidence_type == "tandem_copy":
        return "MEDIUM", "tandem_goi_copy", "goi_tandem_copy_detected"

    if evidence_type == "exon_annotation":
        if exon_count >= 2 and identity >= 60.0 and flanking_support >= 2:
            return "HIGH", "confident_goi", "multi_exon_model_with_flanking_support"
        if identity >= 45.0 and (flanking_support >= 1 or qcov >= 0.65):
            return "MEDIUM", "probable_goi", "modeled_goi_with_partial_support"
        return "LOW", "ambiguous_goi_family_member", "modeled_goi_but_family_context_is_weak"

    if evidence_type == "fallback_hit_span":
        if flanking_support >= 2 and qcov >= 0.75 and identity >= 60.0:
            return "MEDIUM", "probable_goi", "fallback_span_supported_by_flanking_context"
        return "LOW", "ambiguous_goi_family_member", "fallback_span_only"

    if evidence_type == "rescued_exon":
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
) -> Dict[str, Any]:
    confidence, goi_class, reason = _classify_goi_evidence(
        evidence_type=evidence_type,
        identity=identity,
        exon_count=exon_count,
        query_cov=query_cov,
        flanking_support=flanking_support,
    )
    attrs = dict(base_attrs)
    attrs.setdefault("Identity", f"{float(identity or 0.0):.1f}")
    attrs["SynTerraRole"] = "goi"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
    attrs["GOIClass"] = goi_class
    attrs["SyntenyContext"] = _synteny_context_label(flanking_support)
    attrs["BlockFlankingSupport"] = str(max(0, int(flanking_support or 0)))
    attrs["InferenceReason"] = reason
    if query_cov is not None:
        attrs["QueryCoverage"] = _format_attr_float(query_cov)
    if exon_count:
        attrs.setdefault("Exons", str(int(exon_count)))
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
    attrs["SynTerraRole"] = "flanking"
    attrs["EvidenceType"] = evidence_type
    attrs["Confidence"] = confidence
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
        except Exception:
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
            except Exception:
                identity = 0.0
            try:
                start = int(parts[3])
                end = int(parts[4])
            except Exception:
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
            except Exception:
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
                except Exception:
                    pass
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
    except:
        pass
    
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
            subprocess.run(
                ["makeblastdb", "-in", region_fasta, "-dbtype", "nucl", "-out", blast_db_prefix],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            subprocess.run(
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
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
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
            print(f"[{genome_name}] tblastn/makeblastdb not available, skipping tblastn step.", flush=True)
        except subprocess.CalledProcessError:
            # Keep MMseqs/SW path active even when BLAST fails for a block.
            print(f"[{genome_name}] tblastn step failed for this block, continuing.", flush=True)
        
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
                        2,
                        20000.0,
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
        print(f"[{genome_name}] Augmented search failed: {e}", flush=True)
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
                            strand_votes = {'+': 0, '-': 0}
                            for h in work_hits:
                                strand_votes[h.get('strand', '+')] += 1
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
                                    # query protein length.  A gene spanning 50 kb for
                                    # a 70-aa protein is clearly wrong.  Allow up to
                                    # 30x the expected coding length (accounts for
                                    # large introns in some lineages).
                                    if valid_fallback:
                                        max_span_nt = max(3000, query_len * 3 * 30)
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
                                                    f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': parent_id, 'SynTerra_Parent': parent_id, 'Type': 'fallback_hit_span'}, evidence_type='fallback_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=qcov, flanking_support=block_flanking_support), global_start, global_end, strand)}"
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
                                        f"{_mRNA_attrs(_goi_feature_attrs({'ID': copy_id, 'Name': copy['id'], 'SynTerra_Parent': parent_id, 'Type': 'tandem_copy'}, evidence_type='tandem_copy', identity=copy.get('pident', 0), exon_count=1, query_cov=None, flanking_support=block_flanking_support), global_start, global_end, strand)}"
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
                                    f"{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id}, evidence_type='exon_annotation', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, flanking_support=block_flanking_support), global_start, global_end, strand)}"
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
                                            f"{chrom}\trescued_exon\tmRNA\t{gs}\t{ge}\t{hit.get('pident',0):.1f}\t{hit.get('strand','+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': raw_id, 'Name': parent_id, 'SynTerra_Parent': parent_id}, evidence_type='rescued_exon', identity=hit.get('pident', 0), exon_count=1, query_cov=hit_qcov, flanking_support=block_flanking_support), gs, ge, hit.get('strand','+'))}",
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
                                    f"{chrom}\traw_hit\tmRNA\t{nt_s}\t{nt_e}\t{hit.get('pident', 0):.1f}\t{hit.get('strand', '+')}\t.\t{_mRNA_attrs(_goi_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id}, evidence_type='raw_hit', identity=hit.get('pident', 0), exon_count=1, query_cov=qcov, flanking_support=block_flanking_support), nt_s, nt_e, hit.get('strand', '+'))}",
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
                except Exception:
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
                            f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id, 'Type': 'flanking_hit_span'}, evidence_type='flanking_hit_span', identity=avg_pident, exon_count=len(cds_intervals), query_cov=qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
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
                        f"{_mRNA_attrs(_flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id, 'Type': 'flanking_miniprot'}, evidence_type='flanking_miniprot', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='candidate_region_anchor'), global_start, global_end, strand)}"
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
        print(f"Block {block_idx} failed: {e}")

    if os.path.exists(temp_fa):
        os.remove(temp_fa)
    if os.path.exists(query_mini_fa):
        os.remove(query_mini_fa)

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


def process_single_genome(genome_path, db_path, args, home_db_dir, prefix, threads_per_job):
    """
    Worker function to search a single genome.
    Returns: (genome_name, list_of_new_genes, error_message_or_none)
    """
    genome_name = os.path.basename(genome_path)
    if not os.path.exists(genome_path):
        msg = "Genome file not found"
        logger.error(f"[{genome_name}] {msg}.")
        return genome_name, [], msg
    
    unique_id = uuid.uuid4().hex
    hits_file = f"{args.output_dir}/hits/{prefix}{genome_name}.m8"
    tmp_dir = f"{args.output_dir}/tmp_mmseqs_{unique_id}_{genome_name}"
    
    new_genes = []
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
                return genome_name, [], None
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
            return genome_name, [], None
            
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
            return genome_name, [], None
            
        # Optimization: Merge overlapping search regions
        pre_merge_count = len(synteny_blocks)
        synteny_blocks = merge_synteny_blocks(synteny_blocks, args.region_padding)
        if len(synteny_blocks) < pre_merge_count:
            logger.info(f"[{genome_name}] Merged {pre_merge_count} blocks into {len(synteny_blocks)} discrete search regions.")
        else:
            logger.info(f"[{genome_name}] Found {len(synteny_blocks)} discrete syntenic blocks.")

        # Keep only blocks with enough anchor genes to avoid spending hours
        # on singleton/noise loci in fragmented genomes.
        if args.min_block_genes > 1:
            pre_filter_count = len(synteny_blocks)
            synteny_blocks = [b for b in synteny_blocks if b.get('genes_count', 0) >= args.min_block_genes]
            logger.info(
                f"[{genome_name}] Block filter (min_block_genes={args.min_block_genes}): "
                f"{len(synteny_blocks)}/{pre_filter_count} blocks retained."
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
                    except Exception:
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
                            f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id, 'Type': 'rearranged_flanking', 'Rearranged_from': ','.join(block_chroms)}, evidence_type='rearranged_flanking', identity=avg_pident, exon_count=len(exons), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
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
                                    f"{_target_mrna_attrs(off_chrom, global_start, global_end, strand, _flanking_feature_attrs({'ID': new_id, 'Name': parent_id, 'SynTerra_Parent': parent_id, 'Type': 'rearranged_flanking_fallback', 'Rearranged_from': ','.join(block_chroms)}, evidence_type='rearranged_flanking_fallback', identity=avg_pident, exon_count=len(cds_intervals), query_cov=model_qcov, context='cross_chromosome_rearranged'))}"
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
                     "role": attrs.get("SynTerraRole", "goi" if is_goi_query_id(model_id) else "flanking"),
                     "confidence": attrs.get("Confidence", ""),
                     "goi_class": attrs.get("GOIClass", ""),
                     "evidence_type": attrs.get("EvidenceType", attrs.get("Type", "")),
                     "identity": attrs.get("Identity", ""),
                     "n_exons": attrs.get("Exons", ""),
                     "synteny_context": attrs.get("SyntenyContext", ""),
                     "block_flanking_support": attrs.get("BlockFlankingSupport", ""),
                     "query_coverage": attrs.get("QueryCoverage", ""),
                     "target_gene": attrs.get("TargetGene", ""),
                     "target_product": attrs.get("TargetProduct", ""),
                 }

             with open(tsv_out, 'w') as tf:
                 tf.write(
                     "target_id\thome_id\trole\tconfidence\tgoi_class\tevidence_type\t"
                     "identity\tn_exons\tsynteny_context\tblock_flanking_support\t"
                     "query_coverage\ttarget_gene\ttarget_product\n"
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
                             meta.get("evidence_type", ""),
                             meta.get("identity", ""),
                             meta.get("n_exons", ""),
                             meta.get("synteny_context", ""),
                             meta.get("block_flanking_support", ""),
                             meta.get("query_coverage", ""),
                             meta.get("target_gene", ""),
                             meta.get("target_product", ""),
                         ]) + "\n"
                     )
        
        # Expansion DB should only be augmented with GOI-derived models that
        # survived confidence/ambiguity triage. Low-confidence ambiguous/tandem
        # calls are still reported in GFF/plots but should not recursively seed
        # later waves.
        new_genes = []
        suppressed_seed_count = 0
        for g in all_genes:
            gid = g.get('id', '')
            meta = feature_meta.get(gid, {})
            role = meta.get("role", "goi" if is_goi_query_id(gid) else "flanking")
            confidence = meta.get("confidence", "")
            goi_class = meta.get("goi_class", "")
            if role != "goi":
                continue
            if confidence in {"HIGH", "MEDIUM"} and goi_class in {"confident_goi", "probable_goi"}:
                new_genes.append(g)
            else:
                suppressed_seed_count += 1
        if all_genes:
            logger.info(
                f"[{genome_name}] Expansion payload: {len(new_genes)} GOI-derived / "
                f"{len(all_genes)} total annotations"
                + (
                    f" ({suppressed_seed_count} GOI-like annotations withheld from seeding)."
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
            
    return genome_name, new_genes, error_message
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
    
    args = parser.parse_args()
    
    # INPUT VALIDATION
    # 1. Validate required files exist
    if not os.path.exists(args.initial_db):
        logger.error(f"Initial database file not found: {args.initial_db}")
        sys.exit(1)
    
    if not os.path.exists(args.sorted_genomes):
        logger.error(f"Sorted genomes file not found: {args.sorted_genomes}")
        sys.exit(1)
    
    # 2. Validate initial_db is not empty
    if os.path.getsize(args.initial_db) == 0:
        logger.error("Initial database file is empty")
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

    # In auto mode, only keep SW enabled when parasail is available.
    # This avoids per-block fallback overhead on very large block sets.
    if args.enable_smith_waterman and args.sw_method == "auto" and not has_parasail_available():
        logger.warning(
            "parasail not available in auto mode; disabling Smith-Waterman for this run. "
            "Set --sw_method ssearch36 to force the slower fallback."
        )
        args.enable_smith_waterman = False
    
    logger.info(f"Starting iterative search with {args.threads} threads")
    logger.info(f"Parameters: identity>={args.min_identity}%, length>={args.min_length}, evalue<={args.evalue}")
    logger.info(
        f"MMseqs controls: split_memory_limit={args.mmseqs_split_memory_limit}, "
        f"verbosity={args.mmseqs_verbosity}, quiet_subtools={args.quiet_subtools}"
    )
    
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
            last_db = checkpoint_data.get('last_db', None)
            
            if last_db and os.path.exists(last_db):
                logger.info(f"Resuming from wave {start_wave + 1}, using DB: {last_db}")
                current_db = last_db
            else:
                logger.warning("Checkpoint found but DB missing, starting from beginning")
                start_wave = 0
                shutil.copyfile(args.initial_db, current_db)
    else:
        shutil.copyfile(args.initial_db, current_db)
    
    # Parse Genomes and Distances
    genome_entries = []
    with open(args.sorted_genomes, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            gname = parts[0]
            dist = float(parts[1]) if len(parts) > 1 else 0.0
            
            gpath = gname
            if args.genomes_dir:
                if not os.path.isabs(gname): # If gname is not an absolute path, assume it's relative to genomes_dir
                    gpath = os.path.join(args.genomes_dir, os.path.basename(gname))
            
            genome_entries.append({'name': gname, 'path': gpath, 'dist': dist})
    
    # Validate we loaded genomes
    if not genome_entries:
        logger.error("No genomes found in sorted_genomes file")
        sys.exit(1)
            
    logger.info(f"Loaded {len(genome_entries)} genomes.")

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
    
    latest_db = current_db
    
    for i, wave in enumerate(waves):
        # Skip already completed waves
        if i < start_wave:
            logger.info(f"Skipping wave {i+1}/{len(waves)} (already completed)")
            continue
            
        logger.info(f"=== Starting Wave {i+1}/{len(waves)} ({len(wave)} genomes, dist={wave[0]['dist']:.3f}) ===")
        
        # Parallel Execution
        max_workers = min(len(wave), args.threads)
        threads_per_job = max(1, args.threads // max_workers)
        
        logger.info(f"  Running {len(wave)} jobs in parallel with {max_workers} workers, each using {threads_per_job} threads.")
        
        wave_results = []
        wave_errors = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for entry in wave:
                futures.append(
                    executor.submit(process_single_genome, 
                                    entry['path'], latest_db, args, args.home_db_dir, prefix, threads_per_job)
                )
            
            for future in concurrent.futures.as_completed(futures):
                try:
                    gname, new_genes, genome_error = future.result()
                    if genome_error:
                        wave_errors.append((gname, genome_error))
                        continue
                    if new_genes:
                        wave_results.extend(new_genes)
                except Exception as exc:
                    wave_errors.append(("unknown", str(exc)))

        if wave_errors:
            for gname, msg in wave_errors:
                logger.error(f"Wave {i+1} failed for genome '{gname}': {msg}")
            logger.error("Aborting iterative search due to per-genome processing errors.")
            sys.exit(1)

        # Update DB after Wave
        if wave_results:
            logger.info(f"Wave {i+1} completed. Found {len(wave_results)} new genes. Updating DB.")
            
            new_genes_fasta = f"{args.output_dir}/iter_{i+1}_new_genes.faa"
            # wave_results is list of {'id': ..., 'seq': ...}
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
            logger.info(f"Wave {i+1} completed. No new genes found.")
        
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
        
    logger.info(f"Iterative wavefront search complete. Final DB: {expanded_db}")

if __name__ == "__main__":
    main()
