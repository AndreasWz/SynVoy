#!/usr/bin/env python3
"""Specificity of a GOI model against the home proteome (import-only helper).

A model x of the gene of interest q, found in a target genome, is specific to q when no
other protein of the home genome explains it as well. Used by iterative_search_runner.py
(label rework of 2026-10-14, step 2: the seed rule; the labels follow in a later step).

  p*             the best OTHER protein of the home proteome for x: Smith-Waterman of x
                 against every home protein (one per gene, the longest isoform, as
                 gff_to_faa.py writes the proteome), BLOSUM62, BLAST gaps 11/1. Exhaustive
                 on purpose: on the Candida pilot a prefilter (MMseqs2 top 15) missed two
                 real home paralogs that exhaustive Smith-Waterman found.
  q's own gene   the one home protein that is q itself (near-identical to q); never p*.
  margin         bits(q, x) - bits(p*, x).
  real           p* is a real competitor when it aligns to x at E <= 0.01 over the whole
  competitor     proteome, or when it is a home family member of q (q and p* align at
                 E <= 1e-3 with composition-based statistics; see FAMILY_TEST_SEG for what
                 that does and does not do against low-complexity repeats).
  neighbourhood  every protein-coding home gene on the scaffold of this home locus between
                 its outermost flanking genes, the relatives the flanking walk skipped
                 included (extract_flanking_genes.py drops genes that align to q). A p* in
                 there is a gene position cannot tell from q (a home tandem array).

The proteome and its coordinates travel inside the home MMseqs2 database folder that the
search already receives (modules/prepare_home.nf writes home_proteome.faa and
home_genes.tsv there), so no new pipeline input is needed.
"""
import json
import math
import os
import shutil
import subprocess
import tempfile
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    from sequence_utils import SW_K_11_1, SW_LAMBDA_11_1, parse_fasta, sw_align_span
except ImportError:  # pragma: no cover - when imported from outside bin/
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import SW_K_11_1, SW_LAMBDA_11_1, parse_fasta, sw_align_span

PROTEOME_FILE = "home_proteome.faa"
GENES_FILE = "home_genes.tsv"

REAL_COMPETITOR_EVALUE = 0.01   # p* aligns to x at this E over the whole home proteome
FAMILY_MEMBER_EVALUE = 1e-3     # q and p* are one family (pairwise, composition-corrected)
# SEG masking in the family test. The plan (2026-10-14, correction 9) asks for
# composition-based statistics alone, "so repeats do not count". Measured the same day on the
# honeybee proteome: melittin against its two best repeat look-alikes (PEPEP stretch, query
# residues 21-45) gives E 5.3e-4 and 4.5e-4 with -comp_based_stats 2 -seg no (plain
# Karlin-Altschul 8.4e-7 and 6.4e-6), still family members at 1e-3; with -seg yes 2.6 and 5.3.
# Kept as the plan says until that is decided (tmp/label_rework_2026-10-14/STATE.md).
FAMILY_TEST_SEG = "no"
SEED_MIN_MARGIN_BITS = 10.0     # specificity bar of a seed: margin >= 10 bits or no real competitor
OWN_GENE_MIN_IDENTITY = 95.0    # a home protein is q itself at this identity ...
OWN_GENE_MIN_COVERAGE = 0.9     # ... over this share of q

_LN2 = math.log(2.0)


def raw_to_bits(score: float) -> float:
    """Bit score of a raw BLOSUM62 score at BLAST gaps 11/1."""
    return (SW_LAMBDA_11_1 * score - math.log(SW_K_11_1)) / _LN2


class HomeProteome:
    """The home proteome with gene coordinates, and the searches against it."""

    def __init__(self, records: List[Tuple[str, str]], coords: Dict[str, Tuple[str, int, int]]):
        # Sorted by id: among equal scores the first id wins, whatever the file order.
        self.records = sorted((rid, seq) for rid, seq in records if rid and seq)
        self._by_id = dict(self.records)
        self.coords = coords
        self.total_length = sum(len(s) for _r, s in self.records)
        self._own: Dict[str, Optional[str]] = {}
        self._family: Dict[Tuple[str, str], float] = {}

    @classmethod
    def load(cls, folder: Optional[str]) -> Optional["HomeProteome"]:
        """The proteome written into ``folder`` (the home MMseqs2 database folder), or None."""
        if not folder:
            return None
        faa = os.path.join(folder, PROTEOME_FILE)
        if not os.path.exists(faa) or os.path.getsize(faa) == 0:
            return None
        records = [(rid, seq.upper().rstrip("*")) for _h, rid, seq in parse_fasta(faa)]
        coords: Dict[str, Tuple[str, int, int]] = {}
        tsv = os.path.join(folder, GENES_FILE)
        if os.path.exists(tsv):
            with open(tsv) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 4 and f[0] != "gene_id":
                        try:
                            coords[f[0]] = (f[1], int(f[2]), int(f[3]))
                        except ValueError:
                            continue
        return cls(records, coords)

    # ------------------------------------------------------------------ searches
    def scores(self, seq: str) -> List[Tuple[float, str]]:
        """Raw Smith-Waterman score of ``seq`` against every home protein, best first."""
        import parasail  # type: ignore
        if not seq:
            return []
        profile = parasail.profile_create_sat(seq, parasail.blosum62)
        out = [(float(parasail.sw_striped_profile_sat(profile, p, 12, 1).score), rid)
               for rid, p in self.records]
        out.sort(key=lambda t: (-t[0], t[1]))
        return out

    def own_gene(self, q_id: str, q_seq: str) -> Optional[str]:
        """The home protein that is q itself: the best-scoring one, when it is near-identical
        to q (>= 95 % identity over >= 90 % of q). None when q is not in the proteome (a
        query from another species, or an unannotated ORF)."""
        if q_id not in self._own:
            own = None
            for _score, rid in self.scores(q_seq)[:3]:
                aln = sw_align_span(q_seq, self._seq(rid))
                if aln["identity"] >= OWN_GENE_MIN_IDENTITY and aln["query_cov"] >= OWN_GENE_MIN_COVERAGE:
                    own = rid
                    break
            self._own[q_id] = own
        return self._own[q_id]

    def _seq(self, rid: str) -> str:
        return self._by_id.get(rid, "")

    def best_other(self, x_seq: str, exclude: Iterable[str]) -> Tuple[str, float]:
        """(p*, raw score of p* against x), or ("", 0.0) when nothing aligns."""
        skip = set(exclude)
        for score, rid in self.scores(x_seq):
            if rid not in skip:
                return (rid, score) if score > 0 else ("", 0.0)
        return "", 0.0

    def proteome_evalue(self, raw_score: float, x_len: int) -> float:
        """E-value of an alignment of x with this score in a search of the whole proteome."""
        return SW_K_11_1 * max(1, x_len) * max(1, self.total_length) * math.exp(-SW_LAMBDA_11_1 * raw_score)

    def family_evalue(self, q_id: str, q_seq: str, p_id: str) -> float:
        """Pairwise E-value of q and a home protein with composition-based statistics
        (blastp -comp_based_stats 2, read at full precision from its JSON output). Plain
        Karlin-Altschul when blastp is missing."""
        key = (q_id, p_id)
        if key not in self._family:
            self._family[key] = pairwise_evalue(q_seq, self._seq(p_id))
        return self._family[key]

    # ------------------------------------------------------------------ positions
    def neighbourhood(self, flanking_ids: Iterable[str]) -> Optional[Set[str]]:
        """Home genes between the outermost flanking genes of a home locus, on the scaffold
        that holds most of them. None without coordinates."""
        spans = [self.coords[g] for g in flanking_ids if g in self.coords]
        if not spans:
            return None
        chroms: Dict[str, int] = {}
        for c, _s, _e in spans:
            chroms[c] = chroms.get(c, 0) + 1
        chrom = sorted(chroms.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        lo = min(s for c, s, _e in spans if c == chrom)
        hi = max(e for c, _s, e in spans if c == chrom)
        return {g for g, (c, s, e) in self.coords.items() if c == chrom and e >= lo and s <= hi}


_WARNED: Set[str] = set()


def _blastp() -> Optional[str]:
    """blastp on PATH, else next to the running interpreter (a conda env's bin/)."""
    found = shutil.which("blastp")
    if found:
        return found
    import sys
    sibling = os.path.join(os.path.dirname(sys.executable), "blastp")
    return sibling if os.access(sibling, os.X_OK) else None


def pairwise_evalue(a: str, b: str) -> float:
    """E-value of the best alignment of two proteins, composition-corrected (blastp
    -comp_based_stats 2, -seg no). Without blastp: Karlin-Altschul at BLAST 11/1, with a
    warning, because that one calls low-complexity look-alikes a family."""
    if not a or not b:
        return float("inf")
    blastp = _blastp()
    if not blastp and "blastp" not in _WARNED:
        _WARNED.add("blastp")
        import logging
        logging.getLogger(__name__).warning(
            "blastp not found: the home-family test of the specificity check uses plain "
            "Karlin-Altschul statistics, which calls low-complexity look-alikes a family.")
    if blastp:
        with tempfile.TemporaryDirectory(prefix="synvoy_spec_") as tmp:
            qa, sb = os.path.join(tmp, "a.faa"), os.path.join(tmp, "b.faa")
            with open(qa, "w") as fh:
                fh.write(f">a\n{a}\n")
            with open(sb, "w") as fh:
                fh.write(f">b\n{b}\n")
            proc = subprocess.run([blastp, "-query", qa, "-subject", sb, "-comp_based_stats", "2",
                                   "-seg", FAMILY_TEST_SEG, "-outfmt", "15"],
                                  capture_output=True, text=True, check=False)
            if proc.returncode == 0 and proc.stdout.strip():
                try:
                    report = json.loads(proc.stdout)["BlastOutput2"][0]["report"]["results"]
                    hits = (report.get("bl2seq") or [{}])[0].get("hits", [])
                    evalues = [h["evalue"] for hit in hits for h in hit.get("hsps", [])]
                    return min(evalues) if evalues else float("inf")
                except (KeyError, IndexError, ValueError, TypeError):
                    pass
    import parasail  # type: ignore
    score = parasail.sw_striped_sat(a, b, 12, 1, parasail.blosum62).score
    return SW_K_11_1 * len(a) * len(b) * math.exp(-SW_LAMBDA_11_1 * score)


def specificity(proteome: HomeProteome, q_id: str, q_seq: str, x_seq: str, bits_qx: float,
                neighbourhood: Optional[Set[str]]) -> Dict[str, object]:
    """Specificity terms of model x of GOI q (see the module docstring).

    Returns ``best_other`` (p*, "" when no other home protein aligns), ``best_other_bits``,
    ``margin``, ``real_competitor`` (bool), ``in_neighbourhood`` (bool, None without
    coordinates) and ``seed_specific`` (margin >= SEED_MIN_MARGIN_BITS or no real competitor).
    """
    own = proteome.own_gene(q_id, q_seq)
    p_id, p_raw = proteome.best_other(x_seq, [own] if own else [])
    p_bits = raw_to_bits(p_raw) if p_id else 0.0
    margin = bits_qx - p_bits
    real = False
    if p_id:
        real = (proteome.proteome_evalue(p_raw, len(x_seq)) <= REAL_COMPETITOR_EVALUE
                or proteome.family_evalue(q_id, q_seq, p_id) <= FAMILY_MEMBER_EVALUE)
    return {
        "best_other": p_id,
        "best_other_bits": p_bits,
        "margin": margin,
        "real_competitor": real,
        "in_neighbourhood": (None if neighbourhood is None or not p_id else p_id in neighbourhood),
        "seed_specific": (not real) or margin >= SEED_MIN_MARGIN_BITS,
    }


def gff_attrs(spec: Dict[str, object]) -> str:
    """The terms as GFF attributes (appended to a model's mRNA/gene row)."""
    near = spec.get("in_neighbourhood")
    out = (f"BestOtherHome={spec['best_other']};SpecMargin={float(spec['margin']):.1f};"
           f"RealCompetitor={1 if spec['real_competitor'] else 0}")
    if near is not None:
        out += f";BestOtherInNeighbourhood={1 if near else 0}"
    return out
