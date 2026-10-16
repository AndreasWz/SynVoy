"""B2 low_complexity, guards (preset plan of 2026-10-15, PR-8): stricter only.

* the relative test of the flank walk can use the label's family test (blastp with
  composition-based statistics) instead of plain Smith-Waterman: goi_relative_statistic;
  the gene at the GOI's own position is always tested with Smith-Waterman;
* LowComplexity= on GOI models (report_low_complexity);
* decoy mode stays a run-level option (no module may set it: B2 fires on 47 % of the Candida
  genes); in a decoy run the report carries the decoy calls of the SEARCHED cases with the
  upper end of the exact 95 % interval.
legacy never applies them.
"""
from __future__ import annotations

import csv
import random
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import extract_flanking_genes as efg  # noqa: E402
import generate_report  # noqa: E402
import iterative_search_runner as isr  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402
import sequence_utils as su  # noqa: E402

MELITTIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
AA = "ACDEFGHIKLMNPQRSTVWY"
DEFAULTS = {"goi_relative_statistic": "sw", "report_low_complexity": False}


# ============================================================ module and merge
def test_b2_turns_the_guards_on_and_legacy_does_not():
    f = {"seg_frac": 0.3, "repeat_frac": 0.0}
    auto = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="auto", features=f)
    assert auto["settings"] == {"goi_relative_statistic": "blastp_cbs", "report_low_complexity": True}
    legacy = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_single_copy"},
                         modules_spec="legacy", features=f)
    assert legacy["settings"] == DEFAULTS


def test_the_stricter_value_wins(monkeypatch):
    a = sm.Module("Z1", "z1", "", "", lambda f: (True, ""),
                  lambda f: {"goi_relative_statistic": "sw", "zz_guard": False})
    b = sm.Module("Z2", "z2", "", "", lambda f: (True, ""),
                  lambda f: {"goi_relative_statistic": "blastp_cbs", "zz_guard": True})
    monkeypatch.setitem(sm.KEY_CLASS, "zz_guard", ("specificity_guard", "stricter"))  # a boolean guard
    for order in ([a, b], [b, a]):
        assert sm.merge(order, {})[0] == {"goi_relative_statistic": "blastp_cbs", "zz_guard": True}
    assert sm.merge([a], {})[0] == {"goi_relative_statistic": "sw", "zz_guard": False}


def test_decoy_mode_is_a_run_level_option():
    """B2 fires on 47 % of the Candida genes: as a module key, decoy mode would add the decoy
    searches to about half of all runs (the plan's author, 2026-10-15)."""
    assert "decoy_mode" not in sm.KEY_CLASS and "decoy_mode" in sm.FORBIDDEN_FOR_MODULES
    assert all("decoy_mode" not in m.overrides({"query_len": 50, "seg_frac": 0.3, "sp": True,
                                                  "primary_copies": 3.0}) for m in sm.MODULES)
    rogue = sm.Module("Z3", "z3", "", "", lambda f: (True, ""), lambda f: {"decoy_mode": True})
    with pytest.raises(ValueError):
        sm.merge([rogue], {})
    assert "--decoy_mode ${params.decoy_mode}" in (ROOT / "modules" / "iterative_search.nf").read_text()
    keys = (ROOT / "main.nf").read_text().split("def presetDefaultKeys()")[1].split("}")[0]
    assert "decoy_mode" not in keys


# ============================================================ the relative test
@pytest.mark.skipif(shutil.which("blastp") is None, reason="blastp not on PATH")
def test_a_low_complexity_look_alike_is_a_relative_only_for_plain_smith_waterman():
    """A protein with a PEPEP stretch like melittin's pro-region: Smith-Waterman E ~4e-8,
    blastp with composition statistics no hit. A real relative stays one either way."""
    rng = random.Random(7)
    look_alike = (''.join(rng.choice(AA) for _ in range(120)) + "APEPEPAPEPEAEPEPEPAPEPEAEAEPEPEPAPEAEPEPEAEP"
                  + ''.join(rng.choice(AA) for _ in range(150)))
    assert efg.is_goi_relative(look_alike, [MELITTIN], 1e-5, "sw")
    assert not efg.is_goi_relative(look_alike, [MELITTIN], 1e-5, "blastp_cbs")
    relative = list(MELITTIN)
    for i in range(0, len(relative), 7):
        relative[i] = "A" if relative[i] != "A" else "S"
    relative = "MSTQ" + "".join(relative) + "GKR"
    assert efg.is_goi_relative(relative, [MELITTIN], 1e-5, "blastp_cbs")


def test_the_gene_at_the_goi_position_is_always_tested_with_smith_waterman():
    src = (ROOT / "bin" / "extract_flanking_genes.py").read_text()
    i = src.index("def _is_goi_itself")
    body = src[i:src.index("# Walk outward", i)]
    assert "goi_relative_statistic" not in body
    j = src.index("def _is_goi_similar")
    assert "args.goi_relative_statistic" in src[j:i]


# ============================================================ decoys
def _decoy_row(parent, conf):
    return (f"c1\texon_annotation\tmRNA\t1\t300\t.\t+\t.\tID={parent}|g_b0_l1_exon_ann;"
            f"SynVoy_Parent={parent};SynVoyRole=goi;Confidence={conf}")


def test_decoy_cases_count_every_kind_also_without_a_model():
    lines = [_decoy_row("GOI_q~rev", "LOW"), _decoy_row("GOI_q~shuf1", "MEDIUM"),
             _decoy_row("GOI_q~shuf1", "HIGH"), "c1\tx\tCDS\t1\t300\t.\t+\t0\tParent=a"]
    assert isr.decoy_case_counts(lines, ["GOI_q"]) == {
        "decoy_cases": 3, "decoy_high": 1, "decoy_medium_plus": 1, "decoy_ambiguous_plus": 1}
    assert isr.decoy_case_counts([], ["GOI_a", "GOI_b"])["decoy_cases"] == 6
    assert isr.decoy_case_counts([], [])["decoy_cases"] == 0          # no window: nothing searched


def _policy(tmp_path, decoy, ctx=None):
    args = SimpleNamespace(output_dir=str(tmp_path), decoy_mode=decoy, home_goi_ids={"GOI_q", "GOI_q|exon_1"})
    isr._write_target_policy(args, "locus_1_", "gA", 55.0, 8, ctx)
    with open(tmp_path / "hits" / "locus_1_gA.target_policy.tsv") as fh:
        return next(csv.DictReader(fh, delimiter="\t"))


def test_the_target_policy_carries_decoy_counts_only_in_decoy_mode(tmp_path):
    assert "decoy_cases" not in _policy(tmp_path, False)
    row = _policy(tmp_path, True)                         # before the genome end / no window
    assert (row["decoy_cases"], row["decoy_medium_plus"]) == ("0", "0")
    ctx = {"decoy_counts": {"decoy_cases": 3, "decoy_high": 0, "decoy_medium_plus": 1, "decoy_ambiguous_plus": 2}}
    assert _policy(tmp_path, True, ctx)["decoy_ambiguous_plus"] == "2"


def test_the_report_sums_decoy_cases_with_exact_intervals(tmp_path):
    def write(name, rows, decoy=True):
        cols = ["genome", "median_flank_identity"] + (list(generate_report.DECOY_COUNT_COLUMNS) if decoy else [])
        with open(tmp_path / name, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for r in rows:
                w.writerow(r)
    write("a.target_policy.tsv", [["gA", "50", 135, 0, 1, 2]])
    write("b.target_policy.tsv", [["gB", "60", 135, 0, 1, 1]])
    out = generate_report.collect_decoy_rate([tmp_path / "a.target_policy.tsv", tmp_path / "b.target_policy.tsv"])
    assert out["cases"] == 270 and out["genome_searches"] == 2
    assert out["medium_or_better"] == {"calls": 2, "ci95_upper": pytest.approx(0.0265, abs=1e-4)}
    assert out["high"] == {"calls": 0, "ci95_upper": pytest.approx(0.0136, abs=1e-4)}  # never a bare 0 %
    write("c.target_policy.tsv", [["gC", "50"]], decoy=False)
    assert generate_report.collect_decoy_rate([tmp_path / "c.target_policy.tsv"]) is None


@pytest.mark.parametrize("k,n,lo,hi", [(0, 270, 0.0, 0.013570), (2, 270, 0.000898, 0.026500),
                                       (5, 40, 0.041860, 0.268033)])
def test_clopper_pearson_matches_the_beta_quantiles(k, n, lo, hi):
    assert su.clopper_pearson(k, n) == pytest.approx((lo, hi), abs=2e-6)


# ============================================================ LowComplexity
def test_low_complexity_of_the_model_protein(tmp_path):
    protein = "MKT" + "PEPEPEPEPEPEPEPEPEPE" + "WYCHRKLMNF"
    cand = {"gff": ["c1\texon_annotation\tmRNA\t1\t99\t.\t+\t.\tID=GOI_q|g_b0_l1_exon_ann;"
                    "SynVoy_Parent=GOI_q;EvidenceType=exon_annotation"],
            "record": {"id": "GOI_q|g_b0_l1_exon_ann", "seq": protein}}
    off = {"gff": list(cand["gff"]), "record": dict(cand["record"])}
    isr.annotate_model_extras([cand], {}, None, "", low_complexity=True)
    isr.annotate_model_extras([off], {}, None, "")
    value = float(cand["gff"][0].split("LowComplexity=")[1].split(";")[0])
    assert value == pytest.approx(sum(su.low_complexity_mask(protein)) / len(protein), abs=1e-3)
    assert 0.5 < value < 1.0
    assert "LowComplexity" not in off["gff"][0]
