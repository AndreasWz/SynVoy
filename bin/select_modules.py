#!/usr/bin/env python3
"""Preset modules: which kinds of gene this query is, and what each kind changes.

The four presets (preset_short_peptide, _tandem_family, _paralog_discrimination,
_single_copy) were chosen by the first matching rule from the home hits alone. A gene is
not one kind: it can be short, cysteine-rich, secreted and tandem at once (proposal of
2026-10-15, cloud_claude_input/proposal_151026.md). Modules replace that choice:

* every module tests its own trigger on the run's features (hit_profile.json with its
  query_features, goi_info.json, home_context.json); they stack, none comes first;
* a module may only write keys of KEY_CLASS, and each key class has a merge rule, so the
  order of the modules never matters (MERGE below);
* no module writes a label threshold: one judge (synvoy_label) labels every call, and
  the judge's constants are not parameters at all. Modules widen the search, tighten the
  inputs of the specificity test, change how completeness is measured, or add report
  labels.

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
                      [--home_context home_context.json] [--modules auto] [--output x.json]
"""
import argparse
import json
import math
import os
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

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

# Keys a module may write, with the class that decides how two modules' values merge.
#   anchor_hygiene  the largest E (excludes the most relatives from the anchors)
#   search_breadth  the most permissive value (per key: the larger or the smaller one),
#                   then RESOURCE_CAPS
#   figure_regions  one module only (the region caps of the figures)
#   model_status    one module only (ModelStatus coverage bounds)
#   report          on when any module turns it on (attributes that decide nothing)
#   specificity_guard  the stricter value (STRICTER_ORDER)
#   completeness_basis one module only (B3: completeness on the core of the protein)
KEY_CLASS: Dict[str, Tuple[str, str]] = {
    "goi_relative_evalue": ("anchor_hygiene", "max"),
    "goi_candidate_regions_per_query": ("search_breadth", "max"),
    "window_tblastn_word_size": ("search_breadth", "min"),
    "window_tblastn_comp_stats": ("search_breadth", "min"),
    "window_mmseqs_min_length": ("search_breadth", "min"),
    "miniprot_short_orf_max_len": ("search_breadth", "max"),
    "adaptive_score_floor_frac": ("figure_regions", "single"),
    "adaptive_score_floor_abs": ("figure_regions", "single"),
    "adaptive_max_regions": ("figure_regions", "single"),
    "classify_complete_min_qcov": ("model_status", "single"),
    "classify_fragment_max_qcov": ("model_status", "single"),
    "report_copies_in_slot": ("report", "any"),
    "report_low_complexity": ("report", "any"),
    "goi_relative_statistic": ("specificity_guard", "stricter"),
    "core_completeness": ("completeness_basis", "single"),
}
# specificity_guard: the stricter value wins; booleans: on when any module turns them on,
# other keys by this order (last = strictest).
STRICTER_ORDER: Dict[str, Tuple[str, ...]] = {"goi_relative_statistic": ("sw", "blastp_cbs")}
# Upper bounds of search_breadth keys, whatever the modules want.
RESOURCE_CAPS: Dict[str, float] = {"goi_candidate_regions_per_query": 10}
# Keys that change the home side (anchors, home loci), that a measurement banned, or that
# stay run-level options.
# expand_goi_similar: with it on, LY6 scaffold recall fell from 8/9 to 0/9.
# decoy_mode: a calibration switch of the run; B2 fires on 47 % of the Candida genes, so as
# a module key it would add the decoy searches to about half of all runs.
FORBIDDEN_FOR_MODULES = frozenset({"n_flanking_genes", "expand_goi_similar", "decoy_mode"})


# ----------------------------------------------------------------------------- features
def flat_features(hit_profile: Optional[dict], goi_info: Optional[dict],
                  home_context: Optional[dict]) -> Dict[str, Any]:
    """One flat dict of everything a trigger may read (missing inputs give None values)."""
    hp = hit_profile or {}
    qf = hp.get("query_features") or {}
    gi = goi_info or {}
    hc = home_context or {}
    method = gi.get("method")
    exons = gi.get("exons") or []
    real_exons = None
    if gi and method not in ("tandem_duplication", "full_protein_only"):
        real_exons = sum(1 for e in exons if not str(e.get("id", "")).startswith("GOI_copy_"))
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
        "goi_method": method,
        "goi_found_by": gi.get("found_by"),
        "goi_home_gene_ambiguous": gi.get("home_gene_ambiguous"),
        "goi_real_exons": real_exons,
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


def family_candidate_regions(f: Dict[str, Any]) -> int:
    """A1: GOI loci kept per query in a search window, home copies + 1, between 2 and
    RESOURCE_CAPS (merge applies the cap)."""
    return max(2, (home_copies(f) or 1) + 1)


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


# ----------------------------------------------------------------------------- registry
class Module:
    def __init__(self, mid: str, name: str, axis: str, tier: str,
                 trigger: Callable[[Dict[str, Any]], Tuple[bool, str]],
                 overrides: Callable[[Dict[str, Any]], Dict[str, Any]],
                 legacy_alias: Optional[str] = None):
        self.id, self.name, self.axis, self.tier = mid, name, axis, tier
        self.trigger, self.overrides, self.legacy_alias = trigger, overrides, legacy_alias


def _none(_f):
    return {}


MODULES: List[Module] = [
    Module("A1", "local_family", "copy architecture", "after benchmark", _t_local_family,
           lambda f: {"goi_relative_evalue": 1e-3,
                      "goi_candidate_regions_per_query": family_candidate_regions(f),
                      "report_copies_in_slot": True}, "preset_tandem_family"),
    Module("A2", "dispersed_family", "copy architecture", "now", _t_dispersed_family,
           lambda f: {"adaptive_score_floor_frac": 0.45, "adaptive_score_floor_abs": 0.08,
                      "adaptive_max_regions": 4}, "preset_paralog_discrimination"),
    Module("A3", "shared_domain", "copy architecture", "now", _t_shared_domain, _none),
    Module("B1", "short_query", "sequence signal", "after benchmark", _t_short_query,
           lambda f: {"classify_complete_min_qcov": 0.5, "classify_fragment_max_qcov": 0.3,
                      # finding only: the judge's own blastp is untouched
                      "window_tblastn_word_size": 2, "window_tblastn_comp_stats": 0,
                      "window_mmseqs_min_length": 10, "miniprot_short_orf_max_len": 100},
           "preset_short_peptide"),
    Module("B2", "low_complexity", "sequence signal", "after benchmark", _t_low_complexity,
           lambda f: {"goi_relative_statistic": "blastp_cbs", "report_low_complexity": True}),
    Module("B3", "signal_anchor", "sequence signal", "after benchmark", _t_signal_anchor,
           lambda f: {"core_completeness": True}),
    Module("B4", "cys_framework", "sequence signal", "now (advisory)", _t_cys_framework, _none),
    Module("C1", "exon_structure", "structure and locus at home", "now (advisory)",
           _t_exon_structure, _none),
    Module("C2", "home_context", "structure and locus at home", "now (advisory)",
           _t_home_context, _none),
]
BY_KEY = {m.id.lower(): m for m in MODULES}
BY_KEY.update({m.name.lower(): m for m in MODULES})


def module(name: str) -> Module:
    try:
        return BY_KEY[name.strip().lower()]
    except KeyError:
        known = ", ".join(f"{m.name} ({m.id})" for m in MODULES)
        raise ValueError(f"unknown preset module '{name}' (known: {known})") from None


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


# ----------------------------------------------------------------------------- merge
def merge(chosen: Sequence[Module], f: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """(settings, source) of the chosen modules' overrides, merged by key class.
    Raises ValueError for a key outside KEY_CLASS, a forbidden key, or two different
    values for a key only one module may write (a registry error)."""
    wanted: Dict[str, List[Tuple[str, Any]]] = {}
    for m in chosen:
        for k, v in m.overrides(f).items():
            if k in FORBIDDEN_FOR_MODULES:
                raise ValueError(f"module {m.id} writes the forbidden key {k}")
            if k not in KEY_CLASS:
                raise ValueError(f"module {m.id} writes {k}, which has no key class")
            wanted.setdefault(k, []).append((m.id, v))
    settings, source = {}, {}
    for k in sorted(wanted):
        _cls, rule = KEY_CLASS[k]
        vals = sorted(wanted[k], key=lambda t: t[0])
        if rule == "max":
            mid, v = max(vals, key=lambda t: (float(t[1]), t[0]))
        elif rule == "min":
            mid, v = min(vals, key=lambda t: (float(t[1]), t[0]))
        elif rule == "single":
            if len({json.dumps(v) for _m, v in vals}) > 1:
                raise ValueError(f"key {k} is written by {[m for m, _ in vals]} with different values")
            mid, v = vals[0]
        elif rule == "any":
            v = any(bool(x) for _m, x in vals)
        elif rule == "stricter":
            order = STRICTER_ORDER.get(k)
            if order is None:
                v = any(bool(x) for _m, x in vals)
            else:
                v = max((x for _m, x in vals), key=order.index)
        else:  # pragma: no cover - a rule without code is a registry error
            raise ValueError(f"key {k}: unknown merge rule {rule}")
        if k in RESOURCE_CAPS:
            v = min(v, RESOURCE_CAPS[k])
        settings[k] = v
        source[k] = "module:" + "+".join(m for m, _ in vals)
    return settings, source


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
    ap.add_argument("--output", help="Write {features, modules_would_fire, modules_applied, settings} here")
    args = ap.parse_args()
    f = flat_features(_read(args.hit_profile), _read(args.goi_info), _read(args.home_context))
    mode, applied, fired = choose(args.modules, f)
    by_id = {m.id: m for m in MODULES}
    settings, source = merge([by_id[d["id"]] for d in applied], f)
    out = {"modules_mode": mode, "modules_applied": applied, "modules_would_fire": fired,
           "settings": settings, "source": source, "features": f}
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
            fh.write("\n")
    print("modules would fire: " + (", ".join(f"{d['id']} {d['name']} ({d['reason']})" for d in fired) or "none"))


if __name__ == "__main__":
    main()
