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
  q's own gene   the home gene that is q itself; never p*. It is the gene the home step took
                 as q's (goi_info.json: gene_id, found by name or by position; review fix
                 RF-2b, own_gene_id_from_goi_info). Without such a gene (a query from another
                 species that no annotated gene matches, an unannotated ORF) or when the home
                 step flagged its choice as ambiguous, the home protein near-identical to q
                 (>= 95 % identity over >= 90 % of q). That identity rule alone failed where
                 the home annotation of q differs from q (an extra exon, a partial model: the
                 own gene became p*) and took a >= 95 % paralog for q when q is unannotated.
  margin         bits(q, x) - bits(p*, x).
  real           p* is a real competitor when it aligns to x at E <= 0.01 over the whole
  competitor     proteome, or when it is a home family member of q (q and p* align at
                 E <= 1e-3 with composition-based statistics, low-complexity stretches of
                 both masked: see FAMILY_TEST_SEG).
  neighbourhood  every protein-coding home gene on the scaffold of this home locus between
                 its outermost flanking genes, the relatives the flanking walk skipped
                 included (extract_flanking_genes.py drops genes that align to q), the hull
                 extended to the home locus itself; and q's own home locus: q's own gene and
                 its n_flanking_genes neighbours on each side, also in a run from a secondary
                 home locus (review fix RF-2). A p* in there is a gene position cannot tell
                 from q (a home tandem array).
  direction      the family test is not symmetric: q is the query, p the subject. E(q->p) and
                 E(p->q) differ (55 of 949 home pairs disagree at 1e-3 under -seg no, 34 under
                 -seg yes, label session 2026-10-15); a symmetric rule is a v2 question.

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
# SEG in the family test, so that a shared low-complexity stretch does not make two proteins
# one family. In a pairwise blastp run (-subject), -seg masks the QUERY only; with
# -comp_based_stats 1-3 blastp masks the SUBJECT itself whatever -seg says. With -seg no the
# test was one-sided: melittin as the query against its two best look-alikes in the honeybee
# proteome (the PEPEP stretch, query residues 21-45) gave blastp E 5.3e-4 and 4.5e-4, family
# members at 1e-3, while each look-alike as the query gave blastp E 2.6 and 5.3. With -seg yes
# both sides are masked: blastp E 2.6 and 5.3 in both directions (BLAST+ 2.16 and 2.17).
FAMILY_TEST_SEG = "yes"
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
        self._own: Dict[Tuple[str, str], Optional[str]] = {}
        self._own_hood: Dict[Tuple[str, int, str], Optional[Set[str]]] = {}
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

    def own_gene(self, q_id: str, q_seq: str, gene_id: Optional[str] = None) -> Optional[str]:
        """The home protein that is q itself. ``gene_id`` (own_gene_id_from_goi_info) names
        it when the proteome holds that gene. Else the identity rule: the best-scoring home
        protein, when it is near-identical to q (>= 95 % identity over >= 90 % of q); None
        when there is none (a query from another species, or an unannotated ORF)."""
        key = (q_id, gene_id or "")
        if key not in self._own:
            own = None
            if gene_id and gene_id in self._by_id:
                own = gene_id
            else:
                if gene_id:
                    _warn_once(f"own:{gene_id}",
                               f"q's home gene {gene_id} is not in the home proteome: q's own "
                               f"gene is the home protein near-identical to q, if any.")
                for _score, rid in self.scores(q_seq)[:3]:
                    aln = sw_align_span(q_seq, self._seq(rid))
                    if aln["identity"] >= OWN_GENE_MIN_IDENTITY and aln["query_cov"] >= OWN_GENE_MIN_COVERAGE:
                        own = rid
                        break
            self._own[key] = own
        return self._own[key]

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
        (blastp -comp_based_stats 2 -seg yes, read at full precision from its JSON output).
        Plain Karlin-Altschul when blastp is missing."""
        key = (q_id, p_id)
        if key not in self._family:
            self._family[key] = pairwise_evalue(q_seq, self._seq(p_id))
        return self._family[key]

    # ------------------------------------------------------------------ positions
    def neighbourhood(self, flanking_ids: Iterable[str],
                      locus_span: Optional[Tuple[int, int]] = None) -> Optional[Set[str]]:
        """Home genes between the outermost flanking genes of a home locus, on the scaffold
        that holds most of them, the hull extended to ``locus_span`` (start, end of the home
        locus on that scaffold): with flanking genes on one side only, the hull ended before
        q and its adjacent relatives (review fix RF-2). None without coordinates."""
        spans = [self.coords[g] for g in flanking_ids if g in self.coords]
        if not spans:
            return None
        chroms: Dict[str, int] = {}
        for c, _s, _e in spans:
            chroms[c] = chroms.get(c, 0) + 1
        chrom = sorted(chroms.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        lo = min(s for c, s, _e in spans if c == chrom)
        hi = max(e for c, _s, e in spans if c == chrom)
        if locus_span:
            lo, hi = min(lo, int(locus_span[0])), max(hi, int(locus_span[1]))
        return {g for g, (c, s, e) in self.coords.items() if c == chrom and e >= lo and s <= hi}

    def own_neighbourhood(self, q_id: str, q_seq: str, n: int,
                          gene_id: Optional[str] = None) -> Optional[Set[str]]:
        """q's own home locus: q's own gene (own_gene, with ``gene_id``) and the ``n`` home
        genes on each side of it on its scaffold (the genes the flanking walk of the primary
        locus starts from). A run from a secondary home locus knew only its own locus, so a p*
        next to q's own gene was not in its neighbourhood (review fix RF-2). None when q's
        own gene is not found."""
        key = (q_id, int(n), gene_id or "")
        if key not in self._own_hood:
            own = self.own_gene(q_id, q_seq, gene_id)
            hood = None
            if own and own in self.coords and n > 0:
                chrom = self.coords[own][0]
                order = sorted((s, e, g) for g, (c, s, e) in self.coords.items() if c == chrom)
                i = next(k for k, (_s, _e, g) in enumerate(order) if g == own)
                hood = {g for _s, _e, g in order[max(0, i - n): i + n + 1]}
            self._own_hood[key] = hood
        return self._own_hood[key]


_WARNED: Set[str] = set()


def _warn_once(key: str, message: str) -> None:
    if key not in _WARNED:
        _WARNED.add(key)
        import logging
        logging.getLogger(__name__).warning(message)


def own_gene_id_from_goi_info(goi_info) -> Optional[str]:
    """q's own home gene as the home step recorded it (goi_info.json of ANNOTATE_GOI, batch
    HG-3; a path or the loaded dict): its ``gene_id`` when the gene was found by name or by
    position. None without one (a miniprot model, no annotated gene, a goi_info.json written
    before HG-3), and None with a warning when the home step flagged its choice as ambiguous:
    the gene it took may then be a relative with q's own gene as the runner-up, so the
    identity rule decides (review fix RF-2b)."""
    info = goi_info
    if isinstance(goi_info, str):
        try:
            with open(goi_info) as fh:
                info = json.load(fh)
        except (OSError, ValueError):
            return None
    if not isinstance(info, dict):
        return None
    gene_id = info.get("gene_id")
    if not gene_id or info.get("found_by") not in ("name", "position"):
        return None
    if info.get("home_gene_ambiguous"):
        _warn_once(f"ambiguous:{gene_id}",
                   f"The home step flagged its gene {gene_id} as ambiguous (runner-up "
                   f"{info.get('runner_up_gene_id')} at score fraction "
                   f"{info.get('runner_up_score_fraction')}): q's own gene for the specificity "
                   f"terms is the home protein near-identical to q, if any.")
        return None
    return str(gene_id)


# Failed blastp calls of the family test in this process (review fix RF-2): each one is logged,
# and the pair counts as no family (E = inf), as for the place E-value (synvoy_evalue, RF-1c).
FAMILY_TEST_FAILURES = {"count": 0}


def _family_test_failed(reason: str) -> float:
    import logging
    FAMILY_TEST_FAILURES["count"] += 1
    logging.getLogger(__name__).warning(
        "blastp of the home-family test failed (%s); the pair counts as no family (E = inf) "
        "(failure %d of this process)", reason, FAMILY_TEST_FAILURES["count"])
    return float("inf")


def _blastp() -> Optional[str]:
    """blastp on PATH, else next to the running interpreter (a conda env's bin/)."""
    found = shutil.which("blastp")
    if found:
        return found
    import sys
    sibling = os.path.join(os.path.dirname(sys.executable), "blastp")
    return sibling if os.access(sibling, os.X_OK) else None


def pairwise_evalue(a: str, b: str) -> float:
    """E-value of the best alignment of two proteins, composition-corrected, with SEG on
    both proteins (blastp -comp_based_stats 2 -seg yes): ``a`` is the query (q in the family
    test), ``b`` the subject; the two directions differ (module docstring). Read like
    synvoy_evalue.blastp_evalue (``bl2seq`` or ``search`` blocks). A failed call (nonzero exit,
    no output, unparsable JSON or blocks missing) is logged and gives E = inf (review fix
    RF-2; it fell through to plain Karlin-Altschul in silence). Without blastp: Karlin-Altschul
    at BLAST 11/1, with a warning, because that one calls low-complexity look-alikes a
    family (the search stops at the start without blastp)."""
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
        if proc.returncode != 0:
            return _family_test_failed(f"exit {proc.returncode}: {proc.stderr.strip()[-300:]}")
        if not proc.stdout.strip():
            return _family_test_failed("no output")
        try:
            data = json.loads(proc.stdout)
            reports = data["BlastOutput2"]
            evalues = []
            for rep in reports:
                results = rep["report"]["results"]
                if "bl2seq" not in results and "search" not in results:
                    raise KeyError("bl2seq")
                for search in (results["bl2seq"] if "bl2seq" in results else [results["search"]]) or []:
                    for hit in search.get("hits", []):
                        evalues += [float(h["evalue"]) for h in hit.get("hsps", [])]
        except (KeyError, IndexError, ValueError, TypeError) as exc:
            return _family_test_failed(f"unparsable JSON output ({type(exc).__name__}: {exc})")
        return min(evalues) if evalues else float("inf")
    import parasail  # type: ignore
    score = parasail.sw_striped_sat(a, b, 12, 1, parasail.blosum62).score
    return SW_K_11_1 * len(a) * len(b) * math.exp(-SW_LAMBDA_11_1 * score)


def specificity(proteome: HomeProteome, q_id: str, q_seq: str, x_seq: str, bits_qx: float,
                neighbourhood: Optional[Set[str]], n_flanking: int = 0,
                own_gene_id: Optional[str] = None) -> Dict[str, object]:
    """Specificity terms of model x of GOI q (see the module docstring).

    ``neighbourhood`` is that of this home locus (HomeProteome.neighbourhood); with
    ``n_flanking`` > 0 q's own home locus (own_neighbourhood) is added. ``own_gene_id``:
    q's own home gene from the home step (own_gene_id_from_goi_info); None = the identity
    rule. Returns
    ``best_other`` (p*, "" when no other home protein aligns), ``best_other_bits``,
    ``margin``, ``real_competitor`` (bool), ``in_neighbourhood`` (bool, None without
    coordinates) and ``seed_specific`` (margin >= SEED_MIN_MARGIN_BITS or no real competitor).
    """
    if n_flanking:
        own_hood = proteome.own_neighbourhood(q_id, q_seq, n_flanking, own_gene_id)
        if own_hood is not None:
            neighbourhood = (neighbourhood or set()) | own_hood
    own = proteome.own_gene(q_id, q_seq, own_gene_id)
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
