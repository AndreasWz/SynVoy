#!/usr/bin/env python3
"""
Utilities to normalize flanking-gene query FASTA records.

Flanking queries can come in two forms:
1. Explicit exon entries: `gene_id|exon_N` (+ metadata in header).
2. Repeated IDs with many fragment-like records.

This module collapses those inputs to one protein query per parent gene.
"""

import re
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple


EXON_SUFFIX_RE = re.compile(r"^(?P<parent>.+)\|exon_(?P<exon>\d+)$")
PARENT_FIELD_RE = re.compile(r"\bparent=([^\s;]+)")
STRAND_FIELD_RE = re.compile(r"\bstrand=([+-])\b")
NON_AA_RE = re.compile(r"[^A-Za-z*]")


def _clean_protein_sequence(seq: str) -> str:
    """Uppercase and strip non-AA chars and stop codons from query proteins."""
    if not seq:
        return ""
    cleaned = NON_AA_RE.sub("", str(seq).upper())
    return cleaned.replace("*", "")


def _pick_best_sequence(seqs: List[str]) -> str:
    """
    Pick the best representative sequence from duplicates/fragments.

    Priorities:
    1. Longer sequence
    2. Fewer ambiguous X residues
    3. Higher duplicate support
    """
    if not seqs:
        return ""
    counts = Counter(seqs)
    best_seq, _ = max(
        counts.items(),
        key=lambda item: (len(item[0]), -item[0].count("X"), item[1]),
    )
    return best_seq


def collapse_flanking_query_records(
    records: Iterable[Tuple[str, str, str]]
) -> Tuple[List[Tuple[str, str]], Dict[str, int]]:
    """
    Collapse flanking records to one protein per parent gene.

    Args:
        records: iterable of (raw_header, clean_id, sequence)

    Returns:
        (collapsed_records, stats)
        collapsed_records format: [(parent_id, protein_seq), ...]
    """
    grouped = {}
    parent_order: List[str] = []
    stats = {
        "input_records": 0,
        "parents": 0,
        "exon_reconstructed": 0,
        "full_selected": 0,
        "fragment_collapsed": 0,
        "dropped_empty": 0,
    }

    for raw_header, clean_id, seq in records:
        stats["input_records"] += 1
        seq_clean = _clean_protein_sequence(seq)
        if not seq_clean:
            stats["dropped_empty"] += 1
            continue

        exon_match = EXON_SUFFIX_RE.match(clean_id or "")
        if exon_match:
            parent_id = exon_match.group("parent")
            exon_num = int(exon_match.group("exon"))
        else:
            parent_match = PARENT_FIELD_RE.search(raw_header or "")
            parent_id = parent_match.group(1) if parent_match else clean_id
            exon_num = None

        if not parent_id:
            stats["dropped_empty"] += 1
            continue

        if parent_id not in grouped:
            grouped[parent_id] = {
                "full": [],
                "exons": defaultdict(list),
                "strand": "+",
            }
            parent_order.append(parent_id)

        strand_match = STRAND_FIELD_RE.search(raw_header or "")
        if strand_match:
            grouped[parent_id]["strand"] = strand_match.group(1)

        if exon_num is not None:
            grouped[parent_id]["exons"][exon_num].append(seq_clean)
        else:
            grouped[parent_id]["full"].append(seq_clean)

    stats["parents"] = len(grouped)

    collapsed: List[Tuple[str, str]] = []
    for parent_id in parent_order:
        info = grouped[parent_id]
        exons = info["exons"]
        full = info["full"]
        out_seq = ""

        if exons:
            exon_numbers = sorted(exons.keys(), reverse=(info["strand"] == "-"))
            parts = []
            for exon_num in exon_numbers:
                exon_seq = _pick_best_sequence(exons[exon_num])
                if exon_seq:
                    parts.append(exon_seq)
            out_seq = "".join(parts)
            if out_seq:
                stats["exon_reconstructed"] += 1

        if not out_seq and full:
            out_seq = _pick_best_sequence(full)
            if out_seq:
                unique_full = len(set(full))
                if unique_full > 1 or len(full) > 1:
                    stats["fragment_collapsed"] += 1
                else:
                    stats["full_selected"] += 1

        if out_seq:
            collapsed.append((parent_id, out_seq))

    return collapsed, stats

