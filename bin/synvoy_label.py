#!/usr/bin/env python3
"""The label of a GOI candidate (import-only helper; label rework, step 4).

One function judges every GOI candidate wherever one is judged: the choice in a search
window, its label, the report's dedup, the hull rescue, the figures and the benchmark's
first answer, the seed gate. Its inputs have a statistical meaning, its thresholds are
significance levels (tmp/formula_check_2026-10-14/FINAL_for_vscode.md, section 3):

  Eeff        E-value of the home protein against the candidate in the smallest expected
              place that holds it (slot, core, window in a block with E_syn <= 0.1, else
              the genome), times the number of significant blocks (synvoy_evalue)
  specificity no other home protein explains the candidate as well (synvoy_specificity):
              no real competitor, or a margin >= 3.3 bits; within +-3.3 bits the place
              decides for the gene (slot or core of a significant block), unless the best
              other home protein lies in the gene's own home neighbourhood
  complete    the model covers >= 70 % of the query

  HIGH       significant_specific  Eeff <= 1e-6, margin >= 10 bits or no real competitor,
                                   complete
  MEDIUM     significant_specific  Eeff <= 1e-3 and specific
  AMBIGUOUS  not_specific          Eeff <= 1e-3, specificity fails ("another home gene fits
                                   as well"; shown as needs curation, not counted)
  AMBIGUOUS  gene_order_only       1e-3 < Eeff <= 1 in the slot of a significant block, and
                                   specific ("weak match where the gene should be")
  LOW        not_significant       everything else
"""
import math
from typing import Any, Dict, Optional, Tuple

HIGH_EVALUE = 1e-6
MEDIUM_EVALUE = 1e-3
GENE_ORDER_EVALUE = 1.0
HIGH_MARGIN_BITS = 10.0
SPECIFIC_MARGIN_BITS = 3.3
COMPLETE_QCOV = 0.7
EVALUE_TIE = 1e-10          # E-values below this are equal for ranking

LABEL_RANK = {"HIGH": 4, "MEDIUM": 3, "AMBIGUOUS": 2, "LOW": 1}
BASIS_TEXT = {
    "significant_specific": "ortholog",
    "not_specific": "needs curation: another home gene fits as well",
    "gene_order_only": "needs curation: weak match where the gene should be",
    "not_significant": "no significant match",
}


def is_specific(margin: Optional[float], real_competitor: Optional[bool], place: str,
                block_significant: bool, best_other_in_neighbourhood: Optional[bool]) -> bool:
    """No other home protein explains the candidate as well. Without specificity terms
    (no home proteome) the candidate is not held back."""
    if real_competitor is None or not real_competitor:
        return True
    if margin is None:
        return True
    if margin >= SPECIFIC_MARGIN_BITS:
        return True
    if margin > -SPECIFIC_MARGIN_BITS and place in ("slot", "core") and block_significant:
        # positional tie-break, refused where position cannot tell the two apart: a best
        # other protein in the gene's own home neighbourhood (or unknown positions)
        return best_other_in_neighbourhood is False
    return False


def judge(e_eff: Optional[float], place: str, block_significant: bool,
          margin: Optional[float], real_competitor: Optional[bool],
          best_other_in_neighbourhood: Optional[bool], query_cov: Optional[float]) -> Tuple[str, str]:
    """(confidence, basis) of a candidate (module docstring)."""
    e = math.inf if e_eff is None or (isinstance(e_eff, float) and math.isnan(e_eff)) else float(e_eff)
    specific = is_specific(margin, real_competitor, place, block_significant, best_other_in_neighbourhood)
    if e <= MEDIUM_EVALUE:
        if not specific:
            return "AMBIGUOUS", "not_specific"
        strong = (not real_competitor) or margin is None or margin >= HIGH_MARGIN_BITS
        complete = query_cov is not None and query_cov >= COMPLETE_QCOV
        if e <= HIGH_EVALUE and strong and complete:
            return "HIGH", "significant_specific"
        return "MEDIUM", "significant_specific"
    if e <= GENE_ORDER_EVALUE and place == "slot" and block_significant and specific:
        return "AMBIGUOUS", "gene_order_only"
    return "LOW", "not_significant"


def rank_key(confidence: str, e_eff: Optional[float], block_esyn: Optional[float],
             bits: Optional[float], span: int, position: int, ident: str = "") -> Tuple[Any, ...]:
    """Sort key, best first, used everywhere a choice is made: label, then E (below 1e-10
    all equal), then the stronger block (E_syn), then bits, then the shorter span, then
    position (and an id for a total order)."""
    e = math.inf if e_eff is None else max(float(e_eff), EVALUE_TIE)
    esyn = math.inf if block_esyn is None else float(block_esyn)
    return (-LABEL_RANK.get(confidence, 0), e, esyn, -float(bits or 0.0), int(span), int(position), ident)


# ------------------------------------------------------------------ frame guard (3.1)
def _translate(dna: str) -> str:
    try:
        from sequence_utils import translate
    except ImportError:  # pragma: no cover
        import os
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from sequence_utils import translate
    return translate(dna)


def frame_guard(cds: str) -> Tuple[str, str]:
    """(protein, status) of a model's CDS (exons joined in transcript order).

    status ``ok``: the length is a multiple of 3 and the translation has no stop before
    its end; ``fixed``: one of the two other frames over the same exons is clean (it is
    used); ``broken``: none is (a pseudogene is as likely as a frameshift; judge the hit
    regions only, the protein returned is that of frame 0 with its stops)."""
    cds = cds.upper()
    def clean(frame):
        part = cds[frame:]
        part = part[:len(part) - len(part) % 3]
        prot = _translate(part)
        body = prot[:-1] if prot.endswith("*") else prot
        return prot, ("*" not in body) and len(body) > 0
    prot0, ok0 = clean(0)
    if ok0 and len(cds) % 3 == 0:
        return prot0.rstrip("*"), "ok"
    for frame in (1, 2):
        prot, ok = clean(frame)
        if ok:
            return prot.rstrip("*"), "fixed"
    return prot0, "broken"
