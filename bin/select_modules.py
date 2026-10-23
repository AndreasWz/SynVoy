#!/usr/bin/env python3
"""Preset modules: which kinds of gene this query is, and what each kind changes.

The four presets (preset_short_peptide, _tandem_family, _paralog_discrimination,
_single_copy) were chosen by the first matching rule from the home hits alone. A gene is
not one kind: it can be short, cysteine-rich, secreted and tandem at once (proposal of
2026-10-15, cloud_claude_input/proposal_151026.md). Modules replace that choice, and one
profile sets every parameter they touch (profile handoff of 2026-10-20, stage 1):

* every module tests its own trigger on the run's features (hit_profile.json with its
  query_features, goi_info.json, home_context.json); they stack, none comes first;
* the features give a few design quantities, once (design_profile: home copies, the reach
  of the home array, query and exon lengths; PROFILE_SYMBOLS);
* every parameter a module may change has ONE rule (RULES): the modules that switch it on
  (its owners), one formula of the design quantities, a cap, the evidence and a line of
  doc. A rule applies when any of its owners is applied, so two modules can never disagree
  about a key and the order of the modules never matters;
* after the rules, check_constraints tests how the settings must relate to each other; a
  violation that only the rules made is a bug and stops the run before the search;
* no rule writes a label threshold: one judge (synvoy_label) labels every call, and the
  judge's constants are not parameters at all. Switches, n_flanking_genes and the module
  mechanism itself are forbidden to the rules (FORBIDDEN_FOR_MODULES).

`--modules` (nextflow.config):
  legacy   the default until the measured switch: the old chooser's preset is applied
           verbatim (resolve_effective_params.PRESET_OVERRIDES); the triggers are only
           evaluated and recorded (modules_would_fire)
  auto     every module whose trigger fires
  none     no module
  a list   exactly these modules, whether their trigger fires or not:
           "short_query,local_family" (names or ids, A1 ...)
  deltas   "auto,-dispersed_family,+cys_framework"; "none,+short_query" = a list

    select_modules.py --hit_profile hit_profile.json [--goi_info goi_info.json]
                      [--home_context home_context.json] [--modules auto]
                      [--n_flanking_genes 10] [--window_margin_genes 7] [--output x.json]
"""
import argparse
import json
import math
import os
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

# ----------------------------------------------------------------------------- constants
# Every constant here is assumed (proposal §3.1) or today's value; none is tuned.
SHORT_QUERY_MAX = 100          # B1: today's short-peptide rule
LOW_COMPLEXITY_SEG = 0.10      # B2 (assumed)
LOW_COMPLEXITY_REPEAT = 0.25   # B2 (assumed)
TANDEM_MIN_COPIES = 1.5        # A1: today's tandem rule
CYS_FRAMEWORK_MIN = 6          # B4 (assumed): an LU domain has 10 Cys in ~80 aa, insect defensins 6
CONTIG_END_BP = 40_000         # C2 (assumed)
CONTIG_END_MIN_SCAFFOLD = 1_000_000
TYPICAL_WINDOW_BP = 100_000    # B1 headroom, logged only (assumed)
HEADROOM_EVALUE = 1e-3

# The old presets as module sets ("legacy" names them; --preset_override uses them).
LEGACY_ALIAS = {
    "preset_short_peptide": ("B1",),
    "preset_tandem_family": ("A1",),
    "preset_paralog_discrimination": ("A2",),
    "preset_single_copy": (),
}

# Switches: a whole step on or off, a benchmark control or a run-level choice. A switch
# becomes a default after its own gate, never per gene (profile handoff, section 3).
SWITCHES = frozenset({
    "exon_level_search", "enable_smith_waterman", "sw_method", "allow_missing_smith_waterman",
    "deterministic_goi_search", "tandem_copy_parent", "enable_dispersed_goi_rescue",
    "report_nonsyntenic_candidates", "disable_goi_hull_rescue", "disable_strong_synteny_rescue",
    "disable_distant_synteny_rescue", "goi_genome_wide_path", "goi_genome_wide_without_blocks",
    "disable_locus_ownership", "disable_paralog_check", "require_paralog_panel",
    "rank_wave_binning", "disable_wavefront", "disable_anchor_resolution",
    "disable_relative_query", "disable_synteny_collinearity", "synteny_bridge_two_sided",
    "goi_chain_query_aware", "tandem_query_restart", "goi_first_exon_guard",
    "goi_join_leftover_hits", "decoy_mode",
})
# Keys no rule may write (constraint C10): the home side (anchors, home loci), what a
# measurement banned, the module mechanism itself, and every switch.
# expand_goi_similar: with it on, LY6 scaffold recall fell from 8/9 to 0/9.
# decoy_mode: a calibration switch of the run; B2 fires on 47 % of the Candida genes, so as
# a module key it would add the decoy searches to about half of all runs.
FORBIDDEN_FOR_MODULES = frozenset({
    "n_flanking_genes", "expand_goi_similar", "expand_goi_similar_distance",
    "modules", "preset_override", "auto_apply_preset",
}) | SWITCHES


# ----------------------------------------------------------------------------- features
def _real_exons(goi_info: Optional[dict]) -> Optional[List[dict]]:
    """The home gene's own exons: none of a tandem-duplication or whole-protein model, no
    GOI_copy_ rows. None when goi_info holds no such gene model."""
    gi = goi_info or {}
    if not gi or gi.get("method") in ("tandem_duplication", "full_protein_only"):
        return None
    return [e for e in gi.get("exons") or [] if not str(e.get("id", "")).startswith("GOI_copy_")]


def flat_features(hit_profile: Optional[dict], goi_info: Optional[dict],
                  home_context: Optional[dict]) -> Dict[str, Any]:
    """One flat dict of everything a trigger may read (missing inputs give None values)."""
    hp = hit_profile or {}
    qf = hp.get("query_features") or {}
    gi = goi_info or {}
    hc = home_context or {}
    real = _real_exons(gi)
    f = {
        # home hits (profile_hits.py)
        "n_hits": hp.get("n_hits"),
        "n_independent_loci": hp.get("n_independent_loci"),
        "n_strong_loci": hp.get("n_strong_loci"),
        "bit_ratio": hp.get("bit_ratio"),
        "primary_copies": hp.get("primary_copies"),
        "secondary_full_length": hp.get("secondary_full_length"),
        "shared_domain_only": hp.get("shared_domain_only"),
        "loci": hp.get("loci"),
        "locate_retry": hp.get("locate_retry"),
        # the query (query_features.py)
        "query_len": qf.get("query_len", hp.get("query_len")),
        "self_bits": qf.get("self_bits"),
        "seg_frac": qf.get("seg_frac"),
        "repeat_frac": qf.get("repeat_frac"),
        "sp": qf.get("sp"),
        "sp_len": qf.get("sp_len"),
        "gpi": qf.get("gpi"),
        "gpi_len": qf.get("gpi_len"),
        "n_cys": qf.get("n_cys"),
        "cys_max90": qf.get("cys_max90"),
        # the home gene model (annotate_goi_exons.py)
        "goi_method": gi.get("method"),
        "goi_found_by": gi.get("found_by"),
        "goi_home_gene_ambiguous": gi.get("home_gene_ambiguous"),
        "goi_real_exons": None if real is None else len(real),
        # the home neighbourhood (home_context.py)
        "goi_gene": hc.get("goi_gene"),
        "relatives_nearby": hc.get("relatives_nearby"),
        "relative_adjacent": hc.get("relative_adjacent"),
        "relative_proximal": hc.get("relative_proximal"),
        "family_statistic": hc.get("statistic"),
        "scaffold_len": hc.get("scaffold_len"),
        "dist_to_scaffold_end": hc.get("dist_to_scaffold_end"),
        "home_genome_bp": hc.get("home_genome_bp"),
        "telomeric_repeat_at_end": hc.get("telomeric_repeat_at_end"),
    }
    return f


def short_query_headroom(f: Dict[str, Any]) -> Optional[float]:
    """Self bits of the query over the bits an alignment needs for E <= 1e-3 in a typical
    window (100 kb, six frames as protein residues): logged with B1, decides nothing."""
    n, bits = f.get("query_len"), f.get("self_bits")
    if not n or not bits:
        return None
    needed = math.log2(n * (TYPICAL_WINDOW_BP // 3) / HEADROOM_EVALUE)
    return round(bits / needed, 2)


def home_copies(f: Dict[str, Any]) -> Optional[int]:
    """Copies of the gene family at the home locus: the times the query is covered at its
    strongest home locus (rounded half up), or the home gene plus its family members within
    the flank walk, whichever is more. None when neither is known (no hits, no GFF)."""
    copies, rel = f.get("primary_copies"), f.get("relatives_nearby")
    if copies is None and rel is None:
        return None
    return max(int(math.floor(float(copies or 1) + 0.5)), 1 + len(rel or []))


# ----------------------------------------------------------------------------- the profile
# The design quantities, computed once per run and written to effective_params.json as
# `profile`. A rule reads these and nothing else (Rule.inputs names the ones it reads).
PROFILE_SYMBOLS: Dict[str, str] = {
    "c": "home copies: the larger of primary_copies (rounded half up) and 1 + the family "
         "members within the flank walk (home_copies); unknown without hits and GFF",
    "R": "reach of the home array in genes: the largest |rank| of a family member within the "
         "flank walk; 0 without one",
    "Fam": "home loci of the query's family (n_independent_loci)",
    "Lq": "query length (aa)",
    "Lx": "longest real exon of the home gene (aa); real as in goi_real_exons",
    "Lpiece": "min(Lq, Lx), Lq when Lx is unknown: the longest piece one true hit can be, "
              "since a hit does not cross an intron",
    "Lcore": "Lq without the signal peptide and the GPI tail",
    "nx": "real exons of the home gene (goi_real_exons)",
    "n": "n_flanking_genes of the run",
    "W": "window_margin_genes of the run: the neighbourhood around a block in gene pitches; "
         "the DNA searched ends there, so it caps the GOI margin (C1)",
}

# The shipped window_margin_genes, for a profile made without the run's value (W unknown).
WINDOW_MARGIN_GENES_DEFAULT = 7


def _whole(v: Any) -> Optional[float]:
    """A number from a run value: an int when it is whole (7, not 7.0); None when unset."""
    if v in (None, ""):
        return None
    x = float(v)
    return int(x) if x.is_integer() else x


def longest_real_exon_aa(goi_info: Optional[dict]) -> Optional[int]:
    """Lx: the longest real exon of the home gene in aa; None when unknown."""
    lengths = [int(e["length_aa"]) for e in _real_exons(goi_info) or [] if e.get("length_aa")]
    return max(lengths) if lengths else None


def design_profile(f: Dict[str, Any], goi_info: Optional[dict] = None,
                   n_flanking_genes: Optional[Any] = None,
                   window_margin_genes: Optional[Any] = None) -> Dict[str, Any]:
    """The design quantities of PROFILE_SYMBOLS from the flat features, the home gene model
    (for Lx) and the run's n_flanking_genes and window_margin_genes; an unknown quantity is
    None."""
    lq, lx = f.get("query_len"), longest_real_exon_aa(goi_info)
    lpiece = lq if lx is None else (lx if lq is None else min(int(lq), lx))
    lcore = None
    if lq is not None:
        sp = int(f.get("sp_len") or 0) if f.get("sp") else 0
        gpi = int(f.get("gpi_len") or 0) if f.get("gpi") else 0
        lcore = max(0, int(lq) - sp - gpi)
    ranks = [abs(int(r["rank"])) for r in f.get("relatives_nearby") or [] if r.get("rank") is not None]
    return {
        "c": home_copies(f),
        "R": max(ranks) if ranks else 0,
        "Fam": f.get("n_independent_loci"),
        "Lq": lq,
        "Lx": lx,
        "Lpiece": lpiece,
        "Lcore": lcore,
        "nx": f.get("goi_real_exons"),
        "n": None if n_flanking_genes in (None, "") else int(n_flanking_genes),
        "W": _whole(window_margin_genes),
    }


# ----------------------------------------------------------------------------- triggers
def _t_local_family(f):
    copies = f.get("primary_copies")
    rel = f.get("relatives_nearby") or []
    if copies is not None and copies >= TANDEM_MIN_COPIES:
        return True, f"query covered {copies:g} times at its home locus"
    if rel:
        r = rel[0]
        kind = "tandem" if f.get("relative_adjacent") else "proximal"
        return True, (f"{len(rel)} family member(s) within the flank walk, nearest "
                      f"{r['gene_id']} at rank {r['rank']:+d} (E {r['evalue_seg_no']:.2g}; {kind})")
    return False, ""


def _t_dispersed_family(f):
    if (f.get("n_strong_loci") or 0) >= 2 and f.get("secondary_full_length"):
        return True, f"{f['n_strong_loci']} strong home loci, a secondary one covers >= 0.5 of the query"
    return False, ""


def _t_shared_domain(f):
    if f.get("shared_domain_only"):
        return True, f"{f.get('n_strong_loci')} strong home loci, every secondary one covers < 0.5 of the query"
    return False, ""


def _t_short_query(f):
    n = f.get("query_len")
    if n is not None and n < SHORT_QUERY_MAX:
        head = short_query_headroom(f)
        return True, f"query {n} aa" + (f", headroom {head}" if head is not None else "")
    return False, ""


def _t_low_complexity(f):
    seg, rep = f.get("seg_frac") or 0.0, f.get("repeat_frac") or 0.0
    if seg >= LOW_COMPLEXITY_SEG or rep >= LOW_COMPLEXITY_REPEAT:
        return True, f"low complexity {seg:.2f}, repeats {rep:.2f}"
    return False, ""


def _t_signal_anchor(f):
    parts = []
    if f.get("sp"):
        parts.append(f"signal peptide 1-{f.get('sp_len')}")
    if f.get("gpi"):
        parts.append(f"GPI tail of {f.get('gpi_len')} aa")
    return (True, ", ".join(parts)) if parts else (False, "")


def _t_cys_framework(f):
    c = f.get("cys_max90")
    if c is not None and c >= CYS_FRAMEWORK_MIN:
        return True, f"{c} Cys in 90 aa"
    return False, ""


def _t_exon_structure(f):
    n = f.get("goi_real_exons")
    if n is not None and n >= 2:
        return True, f"home gene with {n} exons"
    return False, ""


def _t_home_context(f):
    d, slen = f.get("dist_to_scaffold_end"), f.get("scaffold_len")
    if d is not None and slen is not None and slen >= CONTIG_END_MIN_SCAFFOLD and d < CONTIG_END_BP:
        tel = " (telomeric repeat at that end)" if f.get("telomeric_repeat_at_end") else ""
        return True, f"home gene {d} bp from the end of a {slen} bp scaffold{tel}"
    return False, ""


# ----------------------------------------------------------------------------- modules
class Module:
    """One kind of gene: its trigger, and the switch that enables the rules it owns."""

    def __init__(self, mid: str, name: str, axis: str, tier: str,
                 trigger: Callable[[Dict[str, Any]], Tuple[bool, str]],
                 legacy_alias: Optional[str] = None):
        self.id, self.name, self.axis, self.tier = mid, name, axis, tier
        self.trigger, self.legacy_alias = trigger, legacy_alias


MODULES: List[Module] = [
    Module("A1", "local_family", "copy architecture", "after benchmark", _t_local_family,
           "preset_tandem_family"),
    Module("A2", "dispersed_family", "copy architecture", "now", _t_dispersed_family,
           "preset_paralog_discrimination"),
    Module("A3", "shared_domain", "copy architecture", "now", _t_shared_domain),
    Module("B1", "short_query", "sequence signal", "after benchmark", _t_short_query,
           "preset_short_peptide"),
    Module("B2", "low_complexity", "sequence signal", "after benchmark", _t_low_complexity),
    Module("B3", "signal_anchor", "sequence signal", "after benchmark", _t_signal_anchor),
    Module("B4", "cys_framework", "sequence signal", "now (advisory)", _t_cys_framework),
    Module("C1", "exon_structure", "structure and locus at home", "now (advisory)", _t_exon_structure),
    Module("C2", "home_context", "structure and locus at home", "now (advisory)", _t_home_context),
]
BY_KEY = {m.id.lower(): m for m in MODULES}
BY_KEY.update({m.name.lower(): m for m in MODULES})


def module(name: str) -> Module:
    try:
        return BY_KEY[name.strip().lower()]
    except KeyError:
        known = ", ".join(f"{m.name} ({m.id})" for m in MODULES)
        raise ValueError(f"unknown preset module '{name}' (known: {known})") from None


# ----------------------------------------------------------------------------- the rules
# measured: a run or a count showed it; found in the code: a limit read off the code;
# inferred: argued, not measured yet; report only: an attribute that decides nothing.
EVIDENCE = ("measured", "found in the code", "inferred", "report only")


class Rule:
    """The one rule of a parameter. It applies when any of its owners (module ids) is
    applied; formula(profile, applied module ids) gives the value, which `cap` bounds."""

    def __init__(self, key: str, owners: Sequence[str],
                 formula: Callable[[Dict[str, Any], FrozenSet[str]], Any], text: str,
                 inputs: Sequence[str] = (), cap: Optional[float] = None,
                 evidence: str = "inferred", doc: str = ""):
        self.key, self.owners, self.formula, self.text = key, tuple(owners), formula, text
        self.inputs, self.cap, self.evidence, self.doc = tuple(inputs), cap, evidence, doc


def _const(value: Any) -> Callable[[Dict[str, Any], FrozenSet[str]], Any]:
    return lambda _p, _on: value


def _goi_margin(p: Dict[str, Any], _on: FrozenSet[str]) -> Any:
    w = p.get("W")
    return min(WINDOW_MARGIN_GENES_DEFAULT if w is None else w, max(2, (p.get("R") or 0) + 1))


def _piece_floor(p: Dict[str, Any], _on: FrozenSet[str]) -> int:
    """C7: a per-hit length floor that one true hit can clear: 0.6 of the longest piece,
    between 10 and 50; 50 (the default) when the piece length is unknown."""
    lp = p.get("Lpiece")
    return 50 if lp is None else int(min(50, max(10, round(0.6 * float(lp)))))


def _sw_alignments(p: Dict[str, Any], _on: FrozenSet[str]) -> int:
    return 3 if (p.get("c") or 0) >= 2 or (p.get("nx") or 0) >= 2 else 1


# Stage 1 of the profile handoff: the 15 keys the modules set on 2026-10-15, each as a
# rule with the value it had. Stage 3 adds the search's rules (the model and region caps,
# the GOI margin, Smith-Waterman alignments per piece, the short-piece length floors) and
# raises the regions cap from 10 to 40. Registry order = documentation order.
RULES: List[Rule] = [
    Rule("goi_relative_evalue", ("A1",), _const(1e-3), "1e-3",
         doc="family members next to the GOI: weaker relatives are kept out of the anchors too "
             "(the old preset_tandem_family)"),
    Rule("goi_candidate_regions_per_query", ("A1",),
         lambda p, _on: max(2, (p["c"] or 1) + 1), "max(2, c + 1)", inputs=("c",), cap=40,
         doc="GOI loci kept per query in a search window: each home copy, plus one (an array in "
             "4 sub-clusters lost 15 of 31 copies at the old cap of 10: magic-numbers report #6)"),
    Rule("goi_max_models_per_locus", ("A1",),
         lambda p, _on: max(10, (p["c"] or 1) + 1), "max(10, c + 1)", inputs=("c",), cap=60,
         doc="gene models per GOI locus (masked miniprot re-runs): each home copy, plus one; "
             "copies beyond the cap became rescued_exon fragments (magic-numbers report #3)"),
    Rule("goi_margin_genes", ("A1",), _goi_margin, "min(W, max(2, R + 1))", inputs=("R", "W"),
         doc="the GOI core reaches the farthest home family member next to the GOI, plus one "
             "gene; never beyond the window margin W, where the DNA searched ends (C1)"),
    Rule("sw_max_alignments_per_piece", ("A1", "C1"), _sw_alignments,
         "3 if c >= 2 or nx >= 2, else 1", inputs=("c", "nx"),
         doc="Smith-Waterman reports one alignment per frame: a second copy or exon in the same "
             "frame of a window piece is found only with more (Round 2 could not separate this "
             "switch from --tandem_copy_parent home in arm A, so it is a rule, not a default)"),
    Rule("report_copies_in_slot", ("A1",), _const(True), "on", evidence="report only",
         doc="CopiesInSlot / SlotRank on the GOI models, copies_in_slot in the report"),
    Rule("adaptive_score_floor_frac", ("A2",), _const(0.45), "0.45",
         doc="figure regions, relatives elsewhere in the genome: a higher relative score floor "
             "(the old preset_paralog_discrimination)"),
    Rule("adaptive_score_floor_abs", ("A2",), _const(0.08), "0.08",
         doc="figure regions: a higher absolute score floor (the old preset_paralog_discrimination)"),
    Rule("adaptive_max_regions", ("A2",), _const(4), "4",
         doc="figure regions: at most 4 per genome (the old preset_paralog_discrimination)"),
    Rule("classify_complete_min_qcov", ("B1",), _const(0.5), "0.5",
         doc="ModelStatus only: one missed exon is a large share of a short protein "
             "(the old preset_short_peptide)"),
    Rule("classify_fragment_max_qcov", ("B1",), _const(0.3), "0.3",
         doc="ModelStatus only (the old preset_short_peptide)"),
    Rule("window_tblastn_word_size", ("B1",), _const(2), "2",
         doc="finding only: a short query has few words of 3"),
    Rule("window_tblastn_comp_stats", ("B1",), _const(0), "0",
         doc="finding only: a short query loses score to the composition adjustment; the run "
             "with tblastn's defaults still gives EeffTblastn. Also for a low-complexity short "
             "query: keeping 2 there (C8, stage 3) cost real tblastn window hits on curated "
             "SPIN copies and shortened a counted call, with decoys at 0 MEDIUM+ either way"),
    Rule("window_mmseqs_min_length", ("B1",), _const(10), "10",
         doc="window ORFs from 10 codons, so an exon of a short peptide between two stops is "
             "searched"),
    Rule("miniprot_short_orf_max_len", ("B1",), _const(100), "100",
         doc="every query under 100 aa gets miniprot -L 10 in the windows"),
    Rule("aug_relaxed_length_min", ("B1",), _const(10), "10", evidence="found in the code",
         doc="window hits from 10 aa: B1 opens 10-codon ORFs, but the window search dropped "
             "every GOI hit under 15 aa (C7)"),
    Rule("fallback_short_min_aln_aa", ("B1",), _const(10), "10", evidence="found in the code",
         doc="the short-query fallback judges the same short hits as the window floor, with "
             "the same floor (C7); its 30-bit floor stays"),
    Rule("dispersed_goi_min_alnlen", ("B1", "C1"), _piece_floor,
         "min(50, max(10, round(0.6 x Lpiece)))", inputs=("Lpiece",), evidence="found in the code",
         doc="the dispersed rescue judged each hit against 50 aa, so no piece of a gene whose "
             "exons are shorter could be rescued; its E <= 1e-10 and 40 % identity stay (C7)"),
    Rule("goi_relative_statistic", ("B2",), _const("blastp_cbs"), "blastp_cbs",
         doc="the flank walk's relative test with composition statistics: a low-complexity "
             "look-alike is no relative (a constructed PEPEP protein: SW E 4e-8, blastp_cbs no hit)"),
    Rule("report_low_complexity", ("B2",), _const(True), "on", evidence="report only",
         doc="LowComplexity= on the GOI models"),
    Rule("core_completeness", ("B3",), _const(True), "on",
         doc="completeness for HIGH also on the core of the home protein (no signal peptide, "
             "no GPI tail)"),
]


def validate_registry(rules: Sequence[Rule], modules: Sequence[Module]) -> Dict[str, Rule]:
    """{key: rule}; ValueError for two rules of one key, a forbidden key, an unknown or
    missing owner, an unknown profile input or evidence, or a rule without text or doc."""
    ids = {m.id for m in modules}
    by_key: Dict[str, Rule] = {}
    for r in rules:
        if r.key in by_key:
            raise ValueError(f"two rules for {r.key}")
        if r.key in FORBIDDEN_FOR_MODULES:
            raise ValueError(f"a rule writes the forbidden key {r.key}")
        if not r.owners or set(r.owners) - ids:
            raise ValueError(f"rule {r.key}: owners {r.owners} are not module ids")
        if set(r.inputs) - set(PROFILE_SYMBOLS):
            raise ValueError(f"rule {r.key} reads {sorted(set(r.inputs) - set(PROFILE_SYMBOLS))}, "
                             "which the profile does not have")
        if r.evidence not in EVIDENCE:
            raise ValueError(f"rule {r.key}: evidence '{r.evidence}' is not one of {EVIDENCE}")
        if not r.text or not r.doc:
            raise ValueError(f"rule {r.key} needs its formula as text and a line of doc")
        by_key[r.key] = r
    return by_key


RULES_BY_KEY = validate_registry(RULES, MODULES)


def rule_keys(mid: str) -> List[str]:
    """The parameters set by the rules a module (an id) owns."""
    return [r.key for r in RULES if mid in r.owners]


def _shown(v: Any) -> str:
    if v is None:
        return "unknown"
    if isinstance(v, bool):
        return "on" if v else "off"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def apply_rules(on: Iterable[str], p: Dict[str, Any],
                rules: Optional[Sequence[Rule]] = None) -> Tuple[Dict[str, Any], Dict[str, str], Dict[str, dict]]:
    """(settings, source, derivation) of every rule with an applied owner. source is
    "module:<owner ids>"; a derivation entry gives the value, the formula, the owners, the
    profile inputs with their values, the evidence, the cap and the uncapped value when the
    cap binds, and one readable line ("goi_candidate_regions_per_query 4 = max(2, c + 1),
    c = 3 (A1)")."""
    applied = frozenset(on)
    settings: Dict[str, Any] = {}
    source: Dict[str, str] = {}
    derivation: Dict[str, dict] = {}
    for r in RULES if rules is None else rules:
        owners = sorted(set(r.owners) & applied)
        if not owners:
            continue
        want = r.formula(p, applied)
        capped = r.cap is not None and want > r.cap
        value = r.cap if capped else want
        settings[r.key] = value
        source[r.key] = "module:" + "+".join(owners)
        entry = {"value": value, "rule": r.text, "owners": owners,
                 "inputs": {s: p.get(s) for s in r.inputs}, "evidence": r.evidence}
        line = f"{r.key} {_shown(value)}"
        if r.inputs:
            line += f" = {r.text}, " + ", ".join(f"{s} = {_shown(p.get(s))}" for s in r.inputs)
        if capped:
            entry.update(cap=r.cap, uncapped=want)
            line += f", capped at {_shown(r.cap)} (wanted {_shown(want)})"
        entry["line"] = line + " (" + "+".join(owners) + ")"
        derivation[r.key] = entry
    return settings, source, derivation


def merge(chosen: Sequence[Module], f: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """(settings, source) of these modules' rules on the profile of the features alone (no
    exon lengths, n unknown): apply_rules for callers that hold Module objects."""
    settings, source, _derivation = apply_rules({m.id for m in chosen}, design_profile(f))
    return settings, source


# ----------------------------------------------------------------------------- constraints
# How the settings must relate to each other, checked after the rules and the values the
# user set (profile handoff, section 4). The later stages add theirs with their rules
# (C2 is a code fix of the dispersed rescue, C4/C5 the panels and the figure copies, C6 one
# coverage basis, C9 counts tied to n). C8 (composition statistics kept for a low-complexity
# short query) was tried in stage 3 and taken out again; see window_tblastn_comp_stats. C10,
# no rule writes a forbidden key, is enforced when this module loads (validate_registry).
def _num(v: Any) -> float:
    return float(v)


def _c1_goi_margin_within_window(settings, _source, p, _on):
    key = "goi_margin_genes"
    if key not in settings:
        return None
    w = p.get("W")
    w = WINDOW_MARGIN_GENES_DEFAULT if w is None else w
    return (_num(settings[key]) <= _num(w), (key,),
            f"{key} {settings[key]}, window margin {w} (no DNA is searched beyond it)")


def _c3_caps_reach_the_copies(settings, _source, p, on):
    keys = [k for k in ("goi_candidate_regions_per_query", "goi_max_models_per_locus")
            if k in settings]
    if "A1" not in on or not keys:
        return None
    ok, parts = True, []
    for key in keys:
        cap = RULES_BY_KEY[key].cap
        need = (p.get("c") or 1) + 1
        part = f"{key} {settings[key]}, home copies + 1 = {need}"
        if cap is not None and need > cap:
            part += f": the cap of {cap} binds"
            need = cap
        ok = ok and _num(settings[key]) >= need
        parts.append(part)
    return ok, tuple(keys), "; ".join(parts)


C7_KEYS = ("aug_relaxed_length_min", "fallback_short_min_aln_aa", "dispersed_goi_min_alnlen")


def _c7_one_length_rule(settings, source, p, _on):
    keys = [k for k in C7_KEYS if k in settings and str(source.get(k, "")).startswith("module:")]
    if not keys:
        return None
    lp = p.get("Lpiece")
    bound = None if lp is None else max(10, int(round(0.6 * float(lp))))
    over = [k for k in keys if _num(settings[k]) < 10 or (bound is not None and _num(settings[k]) > bound)]
    detail = (", ".join(f"{k} {settings[k]}" for k in keys)
              + (f"; bound max(10, round(0.6 x Lpiece)) = {bound}" if bound is not None
                 else "; Lpiece unknown"))
    return not over, tuple(keys), detail


def _c11_caps(settings, source, _p, _on):
    keys = [r.key for r in RULES if r.cap is not None and r.key in settings
            and str(source.get(r.key, "")).startswith("module:")]
    if not keys:
        return None
    over = [k for k in keys if _num(settings[k]) > RULES_BY_KEY[k].cap]
    detail = ("above the cap: " + ", ".join(over)) if over else "within the caps: " + ", ".join(keys)
    return not over, tuple(keys), detail


def _c12_model_status_order(settings, _source, _p, _on):
    lo, hi = "classify_fragment_max_qcov", "classify_complete_min_qcov"
    if lo not in settings or hi not in settings:
        return None
    return (_num(settings[lo]) < _num(settings[hi]), (lo, hi),
            f"fragment below {settings[lo]}, complete from {settings[hi]}")


CONSTRAINTS: List[Tuple[str, str, Callable]] = [
    ("C1", "the GOI margin lies within the window margin, where the DNA searched ends",
     _c1_goi_margin_within_window),
    ("C3", "with A1, the GOI loci kept per query and the models per GOI locus reach the home "
           "copies + 1, within their caps (a cap that binds is named)", _c3_caps_reach_the_copies),
    ("C7", "a per-hit length floor a rule sets is at least 10 and at most max(10, 0.6 x Lpiece)",
     _c7_one_length_rule),
    ("C11", "a value set by a rule stays within the rule's cap", _c11_caps),
    ("C12", "ModelStatus: the fragment bound lies below the complete bound", _c12_model_status_order),
]


def check_constraints(settings: Dict[str, Any], source: Dict[str, str], p: Dict[str, Any],
                      on: Iterable[str]) -> List[Dict[str, Any]]:
    """Every constraint: applies, ok, on_violation and detail. `on` = the modules whose
    rules were applied (none under legacy). on_violation is "stop" when a rule set one of
    the constraint's keys and the user set none of them (a violation is then a bug, and the
    run stops before the search), else "warn" (the value the user set wins, as everywhere;
    defaults and presets alone never stop a run)."""
    applied = frozenset(on)
    out = []
    for cid, doc, fn in CONSTRAINTS:
        try:
            res = fn(settings, source, p, applied)
        except (TypeError, ValueError) as exc:
            # only a value set by hand can be no number; the step that reads it decides
            out.append({"id": cid, "doc": doc, "applies": True, "ok": False, "on_violation": "warn",
                        "keys": [], "detail": f"not checked: {exc}"})
            continue
        if res is None:
            out.append({"id": cid, "doc": doc, "applies": False})
            continue
        ok, keys, detail = res
        srcs = [str(source.get(k, "default")) for k in keys]
        stop = (not any(s in ("user_set", "user_cli") for s in srcs)
                and any(s.startswith("module:") for s in srcs))
        out.append({"id": cid, "doc": doc, "applies": True, "ok": bool(ok),
                    "on_violation": "stop" if stop else "warn", "keys": list(keys),
                    "detail": detail})
    return out


# ----------------------------------------------------------------------------- --modules
def parse_spec(spec: Any) -> Tuple[str, Set[str], Set[str]]:
    """(mode, ids forced on, ids forced off) of a --modules value; mode is legacy, auto,
    none or list. Raises ValueError for anything else."""
    if isinstance(spec, bool) or str(spec).strip().lower() in ("true", "false"):
        raise ValueError("--modules got a boolean: write the value in quotes, e.g. "
                         "--modules 'auto,-dispersed_family' (a bare -x is read by Nextflow as true)")
    text = str(spec if spec is not None else "").strip()
    if not text:
        return "legacy", set(), set()
    parts = [p.strip() for p in text.split(",") if p.strip()]
    base = parts[0].lower()
    if base in ("legacy", "auto", "none"):
        rest = parts[1:]
        if base == "legacy" and rest:
            raise ValueError("--modules legacy takes no module list (it applies the old preset)")
    else:
        base, rest = "list", parts
    on, off = set(), set()
    for p in rest:
        if p.startswith("-"):
            if base != "auto":
                raise ValueError(f"--modules: '{p}' removes a module, which only follows 'auto'")
            off.add(module(p[1:]).id)
        else:
            on.add(module(p[1:] if p.startswith("+") else p).id)
    if base == "none" and on:
        base = "list"
    return base, on, off


def evaluate(f: Dict[str, Any]) -> List[Dict[str, str]]:
    """Every module whose trigger fires on the features, in registry order."""
    out = []
    for m in MODULES:
        fired, reason = m.trigger(f)
        if fired:
            out.append({"id": m.id, "name": m.name, "reason": reason})
    return out


def choose(spec: Any, f: Dict[str, Any]) -> Tuple[str, List[Dict[str, str]], List[Dict[str, str]]]:
    """(mode, modules applied, modules whose trigger fires) for a --modules value. In
    legacy and none mode nothing is applied here."""
    mode, on, off = parse_spec(spec)
    fired = evaluate(f)
    fired_ids = {d["id"] for d in fired}
    if mode == "auto":
        ids = (fired_ids - off) | on
    elif mode == "list":
        ids = on
    else:
        ids = set()
    reasons = {d["id"]: d["reason"] for d in fired}
    applied = [{"id": m.id, "name": m.name,
                "reason": reasons.get(m.id) or "named in --modules (trigger did not fire)"}
               for m in MODULES if m.id in ids]
    return mode, applied, fired


def _read(path: Optional[str]) -> Optional[dict]:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        with open(path) as fh:
            text = fh.read().strip()
        return None if not text or text.lower().startswith("sentinel") else json.loads(text)
    except (OSError, ValueError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Which preset modules fire for a run (trigger audit)")
    ap.add_argument("--hit_profile", required=True, help="hit_profile.json of LOCATE_GENE")
    ap.add_argument("--goi_info", help="goi_info.json of ANNOTATE_GOI")
    ap.add_argument("--home_context", help="home_context.json of HOME_CONTEXT")
    ap.add_argument("--modules", default="auto", help="--modules value to evaluate (default auto)")
    ap.add_argument("--n_flanking_genes", type=int, default=None,
                    help="n_flanking_genes of the run, for the profile (default: unknown)")
    ap.add_argument("--window_margin_genes", type=float, default=None,
                    help="window_margin_genes of the run, for the profile (default: unknown, "
                         "the GOI margin rule then caps at the shipped 7)")
    ap.add_argument("--output", help="Write {features, profile, modules_would_fire, modules_applied, "
                                     "settings, source, derivation} here")
    args = ap.parse_args()
    goi_info = _read(args.goi_info)
    f = flat_features(_read(args.hit_profile), goi_info, _read(args.home_context))
    p = design_profile(f, goi_info, args.n_flanking_genes, args.window_margin_genes)
    mode, applied, fired = choose(args.modules, f)
    settings, source, derivation = apply_rules([d["id"] for d in applied], p)
    out = {"modules_mode": mode, "modules_applied": applied, "modules_would_fire": fired,
           "settings": settings, "source": source, "features": f, "profile": p,
           "derivation": derivation}
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
            fh.write("\n")
    print("modules would fire: " + (", ".join(f"{d['id']} {d['name']} ({d['reason']})" for d in fired) or "none"))


if __name__ == "__main__":
    main()
