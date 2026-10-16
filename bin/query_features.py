#!/usr/bin/env python3
"""Features of the query protein that the preset modules test (bin/select_modules.py).

Computed from the query sequence alone, cheaply, before any search:

  query_len, self_bits  length; bit score of the query against itself (BLOSUM62, BLAST 11/1)
  seg_frac              share of residues inside a 12-residue window of Shannon entropy
                        <= 2.2 bits (an approximation of BLAST's -seg 12 2.2 2.5)
  repeat_frac           share of residues covered by an exact 4-mer that recurs at an
                        offset >= 4
  sp, sp_len            an N-terminal signal peptide within residues 1-40 (n-, h- and
                        c-region, von Heijne's rules; no SignalP) and its cleavage site
  gpi, gpi_len          a C-terminal GPI-anchor signal in the last 40 residues (omega site,
                        spacer, hydrophobic tail; Eisenhaber, Bork & Eisenhaber 2001) and the
                        number of residues from omega+1 to the end
  n_cys, cys_max90      cysteines; the most cysteines in any 90-residue window

The signal-peptide and GPI calls are heuristics with assumed constants. They decide
which preset module a run uses, never a label.

    query_features.py --query query.faa [--output features.json]
"""
import argparse
import json
import os
import sys
from typing import Dict, Optional, Tuple

try:
    from sequence_utils import (low_complexity_mask, mean_hydropathy, parse_fasta,
                                repeat_kmer_mask, sw_bits)
except ImportError:  # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (low_complexity_mask, mean_hydropathy, parse_fasta,
                                repeat_kmer_mask, sw_bits)

# Signal peptide (von Heijne 1983, 1986): a short n-region, a hydrophobic h-region of 7-15
# residues, a polar c-region of 3-8 residues; at the cleavage site a small residue at -1
# and a small or aliphatic one at -3 (melittin's precursor is cleaved after ...YIYA: -3 = I).
SP_MAX_END = 40
SP_H_LEN = (7, 15)
SP_H_MIN_KD = 1.6
SP_N_MAX = 15                  # the h-region starts within the first 15 residues after Met
SP_C_LEN = (3, 8)
SP_C_MAX_KD = 1.0              # the c-region is polar: its mean hydropathy stays below this
SP_MINUS1 = frozenset("AGSCT")
SP_MINUS3 = frozenset("AGSCTVIL")

# GPI-anchor signal (Eisenhaber et al. 2001): a small residue at omega and omega+2, a
# moderately polar spacer of about 7 residues, then a hydrophobic stretch of >= 12
# residues that runs to the C-terminus (a few residues may follow it).
GPI_REGION = 40
GPI_TAIL_MIN = 12
GPI_TAIL_KD = 1.6
GPI_TAIL_WINDOW = 7            # windows of this length mark the hydrophobic stretch
GPI_MAX_TRAILING = 5
GPI_SPACER = (4, 12)
GPI_SPACER_ORDER = (7, 6, 8, 5, 9, 4, 10, 11, 12)
GPI_SPACER_MAX_KD = 1.0
GPI_OMEGA = frozenset("GASCDN")
GPI_OMEGA2 = frozenset("GASCT")

SEG_WINDOW, SEG_MAX_ENTROPY = 12, 2.2
REPEAT_K, REPEAT_MIN_OFFSET = 4, 4
CYS_WINDOW = 90


def signal_peptide(seq: str) -> Tuple[bool, int]:
    """(found, cleavage site = signal peptide length) of an N-terminal signal peptide.

    An h-region (7-15 residues, mean Kyte-Doolittle >= 1.6) that starts within the first
    15 residues after the Met, followed by a polar c-region of 3-8 residues whose last
    residue (-1) is small and whose -3 residue is small or aliphatic; the cleavage site at
    or before residue 40. An h-region alone also matches a signal anchor or a
    transmembrane segment, so the c-region and the -1/-3 rule are required. The earliest
    cleavage site of the most hydrophobic h-region wins."""
    s = seq.upper()
    n = min(len(s), SP_MAX_END)
    best: Optional[Tuple[float, int, int]] = None      # (-mean kd, h end, cleavage site)
    for hs in range(1, SP_N_MAX + 1):
        for w in range(SP_H_LEN[0], SP_H_LEN[1] + 1):
            he = hs + w
            if he + SP_C_LEN[0] > n:
                break
            kd = mean_hydropathy(s[hs:he])
            if kd < SP_H_MIN_KD:
                continue
            for c in range(he + SP_C_LEN[0], min(he + SP_C_LEN[1], n) + 1):
                if s[c - 1] not in SP_MINUS1 or s[c - 3] not in SP_MINUS3:
                    continue
                if mean_hydropathy(s[he:c]) > SP_C_MAX_KD:
                    continue
                cand = (-round(kd, 6), c, he)
                if best is None or cand < best:
                    best = cand
                break
    if best is None:
        return False, 0
    return True, best[1]


def gpi_signal(seq: str) -> Tuple[bool, int]:
    """(found, residues from omega+1 to the end) of a C-terminal GPI-anchor signal.

    Searched in the last 40 residues: a hydrophobic stretch (the run of 7-residue windows
    with mean Kyte-Doolittle >= 1.6, together >= 12 residues) that ends at most 5 residues
    before the C-terminus; before it a spacer of 4-12 residues (about 7) of moderate
    polarity; before that the omega site, with small residues at omega and omega+2."""
    s = seq.upper()
    n = len(s)
    lo = max(0, n - GPI_REGION)
    if n - lo < GPI_TAIL_MIN + GPI_SPACER[0] + 3:
        return False, 0
    w = GPI_TAIL_WINDOW
    hydro = [mean_hydropathy(s[i:i + w]) >= GPI_TAIL_KD for i in range(lo, n - w + 1)]
    # the last run of hydrophobic windows
    end_idx = None
    for k in range(len(hydro) - 1, -1, -1):
        if hydro[k]:
            end_idx = k
            break
    if end_idx is None:
        return False, 0
    start_idx = end_idx
    while start_idx > 0 and hydro[start_idx - 1]:
        start_idx -= 1
    t_start, t_end = lo + start_idx, lo + end_idx + w        # [t_start, t_end)
    if t_end - t_start < GPI_TAIL_MIN or n - t_end > GPI_MAX_TRAILING:
        return False, 0
    for spacer in GPI_SPACER_ORDER:
        omega = t_start - spacer - 3                           # 0-based position of omega
        if omega < lo or omega + 2 >= t_start:
            continue
        if s[omega] not in GPI_OMEGA or s[omega + 2] not in GPI_OMEGA2:
            continue
        if mean_hydropathy(s[omega + 3:t_start]) > GPI_SPACER_MAX_KD:
            continue
        return True, n - (omega + 1)
    return False, 0


def query_features(seq: str) -> Dict[str, object]:
    """The feature dict of one query protein (module docstring). Deterministic: fixed
    rounding, so the dict can sit in a cached Nextflow task's output."""
    s = "".join(c for c in seq.upper() if c.isalpha())
    n = len(s)
    sp, sp_len = signal_peptide(s)
    gpi, gpi_len = gpi_signal(s)
    try:
        self_bits = round(sw_bits(s, s), 1) if s else 0.0
    except RuntimeError:            # no parasail
        self_bits = None
    cys = [i for i, c in enumerate(s) if c == "C"]
    cys_max = 0
    for i, p in enumerate(cys):
        k = i
        while k < len(cys) and cys[k] < p + CYS_WINDOW:
            k += 1
        cys_max = max(cys_max, k - i)
    return {
        "query_len": n,
        "self_bits": self_bits,
        "seg_frac": round(sum(low_complexity_mask(s, SEG_WINDOW, SEG_MAX_ENTROPY)) / n, 3) if n else 0.0,
        "repeat_frac": round(sum(repeat_kmer_mask(s, REPEAT_K, REPEAT_MIN_OFFSET)) / n, 3) if n else 0.0,
        "sp": sp,
        "sp_len": sp_len,
        "gpi": gpi,
        "gpi_len": gpi_len,
        "n_cys": len(cys),
        "cys_max90": cys_max,
    }


def first_sequence(path: str) -> str:
    """The first record of a FASTA file ('' when there is none)."""
    for _h, _id, seq in parse_fasta(path):
        return seq
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description="Features of a query protein (preset modules)")
    ap.add_argument("--query", required=True, help="Query protein FASTA (first record)")
    ap.add_argument("--output", help="Write the features as JSON here (default: stdout)")
    args = ap.parse_args()
    feats = query_features(first_sequence(args.query))
    text = json.dumps(feats, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(text + "\n")
    else:
        print(text)


if __name__ == "__main__":
    main()
