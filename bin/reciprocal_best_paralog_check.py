#!/usr/bin/env python3
"""Reciprocal-best paralog check per (locus, genome) — docs/TODO.md §1j Phase B.

For paralog-family queries (multi-sequence home FASTA like TP53 + TP63 + TP73),
the iterative search at MEDIUM confidence can drag in sibling paralogs at the
"wrong" locus. This script aligns each recovered GOI target protein against
every home paralog with Smith–Waterman (parasail) and emits a TSV recording
the best-matching paralog + bitscore gap to the runner-up.

generate_report.py reads these TSVs, joins them onto the GOI annotation list,
and flags a call as ``paralog_confusion`` when its best home paralog differs
from the *modal* best-paralog at the same locus (with a bitscore gap >= 5).
The flag goes into ``self_consistency.flags`` in synvoy_report.json so the
paper text can cite an exact confusion count.

Single-paralog queries (one sequence in --home_query) are a no-op: there's
nothing to be reciprocal against.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Iterator

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import parse_gff_attributes  # noqa: E402

try:
    import parasail  # type: ignore
except ImportError:
    parasail = None


HEADER = [
    "genome", "locus", "chrom", "start", "end",
    "mrna_id", "target_len_aa",
    "best_paralog", "best_bit", "best_identity",
    "second_paralog", "second_bit", "bitscore_gap",
    "n_paralogs_compared",
]


def _iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    name = ""
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        yield name, "".join(parts)


def _find_goi_models(gff_path: str) -> list[dict]:
    """Each GOI mRNA in the per-locus region GFF, keyed enough to look up its
    protein sequence in the matching .faa later."""
    out = []
    try:
        with open(gff_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "mRNA":
                    continue
                attrs = parse_gff_attributes(parts[8])
                role = (attrs.get("SynVoyRole") or "").strip().lower()
                # Treat ID-prefix as a fallback when SynVoyRole is missing, matching
                # the upstream classifier in generate_report.summarize_annotations.
                if role != "goi":
                    model_id = attrs.get("ID", "")
                    parent = attrs.get("SynVoy_Parent", "")
                    if not (model_id.startswith("GOI_") or parent.startswith("GOI_")):
                        continue
                try:
                    s, e = int(parts[3]), int(parts[4])
                except ValueError:
                    continue
                if e < s:
                    s, e = e, s
                out.append({
                    "mrna_id": attrs.get("ID", ""),
                    "chrom": parts[0],
                    "start": s,
                    "end": e,
                    "confidence": (attrs.get("Confidence") or "UNKNOWN").upper(),
                })
    except OSError:
        pass
    return out


def _align(query_seq: str, paralog_seq: str) -> tuple[float, float]:
    """Return (bitscore_proxy, percent_identity) for query_seq vs paralog_seq.

    parasail returns the raw SW score; we treat that as a bitscore proxy (no need
    for an evalue-based calibration since we only use it to RANK paralogs against
    each other for the same query). Identity is computed from the traceback CIGAR.
    """
    if parasail is None:
        raise RuntimeError("parasail is required for §1j reciprocal-best paralog check")
    res = parasail.sw_trace_striped_16(query_seq, paralog_seq, 11, 1,
                                       parasail.blosum62)
    score = float(res.score)
    cigar = res.cigar
    # cigar.seq is a numpy array in newer parasail; check length, not truthiness.
    cigar_ops = getattr(cigar, "seq", None)
    if cigar_ops is None or len(cigar_ops) == 0:
        return score, 0.0
    matches = 0
    aligned = 0
    for op in cigar_ops:
        length = op >> 4
        code = op & 0xF
        # parasail cigar ops: 0=M (match-or-mismatch), 7==, 8=X, 1=I, 2=D
        if code in (0, 7, 8):
            aligned += length
            if code == 7:
                matches += length
    if aligned == 0 or matches == 0:
        # cigar didn't expose match/mismatch (alphabet=blosum) — approximate from
        # the aligned span length as a fraction of the query, capped at 1.0.
        ident = min(1.0, score / max(1, len(query_seq))) * 100
    else:
        ident = 100.0 * matches / aligned
    return score, ident


def _best_two_paralogs(query_seq: str, paralogs: dict[str, str]) -> tuple[str, float, float, str, float]:
    """Return (best_name, best_bit, best_identity, second_name, second_bit).
    ``second_name``/``second_bit`` are '' / 0.0 when only one paralog is available."""
    scored = []
    for name, seq in paralogs.items():
        if not seq:
            continue
        bit, ident = _align(query_seq, seq)
        scored.append((bit, ident, name))
    scored.sort(reverse=True)
    if not scored:
        return "", 0.0, 0.0, "", 0.0
    best_bit, best_ident, best_name = scored[0]
    if len(scored) == 1:
        return best_name, best_bit, best_ident, "", 0.0
    second_bit, _, second_name = scored[1]
    return best_name, best_bit, best_ident, second_name, second_bit


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--home_query", required=True,
                    help="Home GOI multi-FASTA — each sequence is a paralog")
    ap.add_argument("--target_faa", required=True,
                    help="Per-locus per-genome region .faa")
    ap.add_argument("--target_gff", required=True,
                    help="Per-locus per-genome region GFF (to filter to GOI mRNAs)")
    ap.add_argument("--locus_id", required=True)
    ap.add_argument("--genome_name", required=True)
    ap.add_argument("--output", required=True, help="Output TSV (HEADER above)")
    args = ap.parse_args()

    # 1) Read the home paralog panel.
    paralogs = {name: seq for name, seq in _iter_fasta(args.home_query) if seq}
    with open(args.output, "w") as out:
        out.write("\t".join(HEADER) + "\n")

        # 2) Single-paralog query — no reciprocal-best to check; emit the header
        # only so the downstream collect() still has a file.
        if len(paralogs) < 2:
            return

        # 3) Find the recovered GOI models from the per-locus region GFF.
        goi_models = _find_goi_models(args.target_gff)
        if not goi_models:
            return
        models_by_id = {m["mrna_id"]: m for m in goi_models if m["mrna_id"]}

        # 4) Read target proteins from the per-locus .faa; only keep ones whose
        # ID matches a GOI mRNA.
        for name, seq in _iter_fasta(args.target_faa):
            model = models_by_id.get(name)
            if model is None or not seq:
                continue
            best, best_bit, best_ident, second, second_bit = \
                _best_two_paralogs(seq, paralogs)
            if not best:
                continue
            gap = best_bit - second_bit
            row = [
                args.genome_name, args.locus_id, model["chrom"],
                str(model["start"]), str(model["end"]),
                name, str(len(seq)),
                best, f"{best_bit:.1f}", f"{best_ident:.1f}",
                second, f"{second_bit:.1f}", f"{gap:.1f}",
                str(len(paralogs)),
            ]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
