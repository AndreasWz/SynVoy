"""Preset modules (2026-10-15, PR-1 of the preset plan): features, registry, merge, resolution.

* query_features.py   features of the query sequence alone
* profile_hits.py     per-locus union coverage, loci, secondary_full_length, shared_domain_only
* home_context.py     family members next to the home gene (the label's family test)
* select_modules.py   triggers, key classes, an order-free merge, the --modules value
* resolve_effective_params.py  --modules legacy reproduces the old preset exactly; auto,
                      lists and deltas apply modules below the preset layer, and a value the
                      user set always wins
"""
from __future__ import annotations

import itertools
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import auto_select_preset as asp  # noqa: E402
import generate_report  # noqa: E402
import home_context as hc  # noqa: E402
import profile_hits  # noqa: E402
import query_features as qf  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
import select_modules as sm  # noqa: E402
import synvoy_specificity  # noqa: E402

# UniProt P01501 (melittin precursor) and Q16553 (LY6E), public sequences.
MELITTIN = "MKFLVNVALVFMVVYISYIYAAPEPEPAPEPEAEADAEADPEAGIGAVLKVLTTGLPALISWIKRKRQQG"
LY6E = ("MKIFLPVLLAALLGVERASSLMCFSCLNQKSNLYCLKPTICSDQDNYCVTVSASAGIGNLVTFGHSLSKTCSPACPIPEGVNVG"
        "VASMGISCCQSFLCNFSAADGGLRASVTLLGAGLLLSLLPALLRFGP")

DEFAULTS = {
    "sw_min_score": 20, "sw_min_identity": 10.0, "min_hit_identity": 10, "min_hit_length": 10,
    "classify_fragment_max_qcov": 0.4, "classify_complete_min_qcov": 0.7,
    "expand_goi_similar": False, "goi_family_tokens": "", "goi_relative_evalue": 1e-5,
    "adaptive_score_floor_frac": 0.3, "adaptive_score_floor_abs": 0.03,
    "adaptive_max_regions": 6, "adaptive_unique_gene_floor": 3,
}


def _profile(query_len=300, n_hits=4, n_strong=1, n_loci=1, ratio=None, copies=1.0,
             secondary_qcov=None, seq=None):
    loci = [{"chrom": "c1", "start": 1, "end": 900, "bit": 400.0, "qcov": 1.0}]
    if secondary_qcov is not None:
        loci.append({"chrom": "c2", "start": 1, "end": 900, "bit": 300.0, "qcov": secondary_qcov})
    return {
        "query_len": query_len, "n_hits": n_hits, "n_hits_raw": n_hits,
        "n_independent_loci": n_loci, "n_strong_loci": n_strong, "top_bit": 400.0,
        "second_bit": 300.0 if n_strong >= 2 else None, "bit_ratio": ratio,
        "primary_copies": copies, "loci": loci,
        "secondary_full_length": bool(n_strong >= 2 and (secondary_qcov or 0) >= 0.5),
        "shared_domain_only": bool(n_strong >= 2 and (secondary_qcov or 0) < 0.5),
        "query_features": qf.query_features(seq or ("MKT" + "AILVG" * (query_len // 5))[:query_len]),
    }


# ============================================================ query features
def test_ly6e_has_a_signal_peptide_a_gpi_tail_and_ten_cysteines():
    f = qf.query_features(LY6E)
    assert (f["sp"], f["sp_len"]) == (True, 20)            # UniProt: signal 1-20
    assert f["gpi"] and 25 <= f["gpi_len"] <= 35            # UniProt: omega at 101 of 131
    assert f["n_cys"] == 10 and f["cys_max90"] == 10
    assert LY6E.endswith("GGLRASVTLLGAGLLLSLLPALLRFGP")


def test_melittin_precursor_has_its_signal_peptide_and_low_complexity():
    f = qf.query_features(MELITTIN)
    assert (f["sp"], f["sp_len"]) == (True, 21)            # cleaved after ...YIYA
    assert not f["gpi"] and f["n_cys"] == 0
    assert f["seg_frac"] >= 0.10                            # the PEPEP pro-region


def test_low_complexity_and_repeats_are_measured():
    polyq = "MKT" + "Q" * 30 + "LVNVALVFMVVYISYIYAWERTYIPASDFGHKLCVNM"
    assert qf.query_features(polyq)["seg_frac"] >= 0.10
    assert qf.query_features("MKTAPEPAPEPAPEPAPEPAPEPLVNW")["repeat_frac"] >= 0.5
    plain = qf.query_features("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    assert not plain["sp"] and not plain["gpi"]


def test_query_features_are_deterministic_and_rounded():
    assert qf.query_features(LY6E) == qf.query_features(LY6E)
    assert json.dumps(qf.query_features(LY6E), sort_keys=True) == json.dumps(
        qf.query_features(LY6E.lower()), sort_keys=True)


# ============================================================ profile_hits
def _hit(chrom, qs, qe, ts, te, bits, ev=1e-30, pid=90.0):
    return {"chrom": chrom, "qstart": qs, "qend": qe, "tstart": ts, "tend": te,
            "pident": pid, "bits": bits, "evalue": ev}


def test_locus_coverage_is_a_union_not_a_sum():
    hits = [_hit("c1", 1, 60, 1000, 1180, 100.0), _hit("c1", 40, 100, 2000, 2183, 90.0)]
    p = profile_hits.build_profile(hits, 100, 50_000, 0.5)
    assert p["loci"][0]["qcov"] == 1.0                     # 1-100, not 60 + 61
    assert p["primary_copies"] == 1.21                     # the copy number keeps the sum


def test_profile_keeps_every_old_key_and_adds_the_new_ones():
    hits = [_hit("c1", 1, 100, 1000, 1300, 200.0), _hit("c2", 1, 40, 5000, 5120, 150.0)]
    p = profile_hits.build_profile(hits, 100, 50_000, 0.5, query_seq=MELITTIN)
    old = {"query_len", "n_hits", "n_hits_raw", "n_independent_loci", "n_strong_loci", "top_bit",
           "second_bit", "bit_ratio", "primary_copies", "identity_min", "identity_q1",
           "identity_median", "identity_q3", "identity_max", "identity_iqr"}
    assert old <= set(p)
    assert (p["n_strong_loci"], p["secondary_full_length"], p["shared_domain_only"]) == (2, False, True)
    assert p["query_features"]["sp_len"] == 21
    empty = profile_hits.build_profile([], 70, 50_000, 0.5, query_seq=MELITTIN)
    assert empty["loci"] == [] and empty["secondary_full_length"] is False
    assert empty["query_features"]["query_len"] == 70


# ============================================================ registry and merge
def test_every_key_a_module_writes_has_a_class_a_default_and_a_reader():
    keys = set(re.findall(r"'([a-z_]+)'", re.search(
        r"def presetDefaultKeys\(\) \{\s*return \[(.*?)\]", (ROOT / "main.nf").read_text(), re.S).group(1)))
    readers = "".join(p.read_text() for p in (ROOT / "modules").glob("*.nf"))
    big = sm.flat_features(_profile(query_len=50, copies=3.0, n_strong=3, secondary_qcov=0.9), None, None)
    for m in sm.MODULES:
        for k in m.overrides(big):
            assert k in sm.KEY_CLASS, (m.id, k)
            assert k in keys, (m.id, k)
            assert f"settings.{k}" in readers, (m.id, k)
            assert k not in sm.FORBIDDEN_FOR_MODULES, (m.id, k)
    assert not set(sm.KEY_CLASS) & sm.FORBIDDEN_FOR_MODULES


def test_merge_does_not_depend_on_the_order_of_the_modules():
    f = sm.flat_features(_profile(query_len=50, copies=3.0, n_strong=3, secondary_qcov=0.9), None, None)
    chosen = [m for m in sm.MODULES if m.id in ("A1", "A2", "B1", "B2", "C1")]
    first = sm.merge(chosen, f)
    for perm in itertools.permutations(chosen):
        assert sm.merge(list(perm), f) == first


def test_the_old_presets_are_module_sets_with_the_same_values():
    by_id = {m.id: m for m in sm.MODULES}
    f = sm.flat_features(_profile(), None, None)
    for preset, ids in sm.LEGACY_ALIAS.items():
        settings, _src = sm.merge([by_id[i] for i in ids], f)
        assert settings == rep.PRESET_OVERRIDES[preset], preset
    assert set(sm.LEGACY_ALIAS) == set(rep.PRESET_OVERRIDES)


def test_groovy_module_list_matches_the_registry():
    text = (ROOT / "main.nf").read_text()
    block = re.search(r"def presetModules\(\) \{\s*return \[(.*?)\]", text, re.S).group(1)
    groovy = dict(re.findall(r"([A-C]\d): '([a-z_]+)'", block))
    assert groovy == {m.id: m.name for m in sm.MODULES}


@pytest.mark.parametrize("spec,expected", [
    ("legacy", ("legacy", set(), set())),
    ("", ("legacy", set(), set())),
    ("auto", ("auto", set(), set())),
    ("none", ("none", set(), set())),
    ("short_query,A1", ("list", {"B1", "A1"}, set())),
    ("auto,-dispersed_family,+cys_framework", ("auto", {"B4"}, {"A2"})),
    ("none,+short_query", ("list", {"B1"}, set())),
])
def test_modules_spec(spec, expected):
    assert sm.parse_spec(spec) == expected


@pytest.mark.parametrize("spec", [True, "true", "bogus", "legacy,short_query", "none,-B1", "-B1"])
def test_bad_modules_specs_are_refused(spec):
    with pytest.raises(ValueError):
        sm.parse_spec(spec)


# ============================================================ triggers
def test_a_family_member_next_to_the_gene_fires_local_family_under_auto_only():
    """LY6E-like: one copy at home, SLURP1-like relative four genes to the left."""
    hp = _profile(query_len=131, seq=LY6E)
    home = {"relatives_nearby": [{"gene_id": "gene-X", "rank": -4, "distance_bp": 276000,
                                  "evalue_seg_no": 2.8e-7, "evalue_seg_yes": 3e-7}],
            "relative_adjacent": False, "relative_proximal": True, "statistic": "blastp"}
    f = sm.flat_features(hp, None, home)
    fired = {d["id"] for d in sm.evaluate(f)}
    assert {"A1", "B3", "B4"} <= fired and "B1" not in fired
    auto = rep.resolve(DEFAULTS, auto_preset={"preset": asp.select_preset(hp)[0]},
                       modules_spec="auto", features=f)
    assert auto["settings"]["goi_relative_evalue"] == 1e-3
    assert auto["source"]["goi_relative_evalue"] == "module:A1"
    legacy = rep.resolve(DEFAULTS, auto_preset={"preset": asp.select_preset(hp)[0]},
                         modules_spec="legacy", features=f)
    assert legacy["preset_applied"] == "preset_single_copy"
    assert legacy["settings"] == DEFAULTS
    assert "A1" in {d["id"] for d in legacy["modules_would_fire"]}


def test_shared_domain_and_dispersed_family_are_told_apart():
    luciferase = sm.flat_features(_profile(n_strong=4, n_loci=50, ratio=1.3, secondary_qcov=0.3), None, None)
    fired = {d["id"] for d in sm.evaluate(luciferase)}
    assert "A3" in fired and "A2" not in fired
    tp53 = sm.flat_features(_profile(n_strong=3, n_loci=3, ratio=1.2, secondary_qcov=0.8), None, None)
    fired = {d["id"] for d in sm.evaluate(tp53)}
    assert "A2" in fired and "A3" not in fired


def test_exon_structure_counts_real_exons_only():
    exons = [{"id": "GOI_q|exon_1"}, {"id": "GOI_q|exon_2"}, {"id": "GOI_copy_1"}]
    f = sm.flat_features(None, {"method": "gff_annotation", "exons": exons}, None)
    assert f["goi_real_exons"] == 2 and sm._t_exon_structure(f)[0]
    f = sm.flat_features(None, {"method": "tandem_duplication", "exons": exons}, None)
    assert f["goi_real_exons"] is None and not sm._t_exon_structure(f)[0]


def test_contig_end_needs_a_long_scaffold():
    near = {"dist_to_scaffold_end": 12_000, "scaffold_len": 3_000_000}
    assert sm._t_home_context(sm.flat_features(None, None, near))[0]
    assert not sm._t_home_context(sm.flat_features(None, None, dict(near, scaffold_len=200_000)))[0]
    assert not sm._t_home_context(sm.flat_features(None, None, dict(near, dist_to_scaffold_end=90_000)))[0]


# ============================================================ resolution
GRID = [
    dict(n_hits=0, query_len=300),
    dict(query_len=50), dict(query_len=99), dict(query_len=100), dict(query_len=101),
    dict(n_strong=2, ratio=1.2, secondary_qcov=0.9), dict(n_strong=2, ratio=2.0, secondary_qcov=0.9),
    dict(copies=1.4), dict(copies=1.5), dict(copies=2.0), dict(),
]


@pytest.mark.parametrize("kw", GRID)
def test_legacy_reproduces_the_old_preset_exactly(kw):
    """Same settings, preset and per-key source as the resolver before the modules."""
    hp = _profile(**kw)
    preset = asp.select_preset(hp)[0]
    res = rep.resolve(DEFAULTS, auto_preset={"preset": preset}, modules_spec="legacy",
                      features=sm.flat_features(hp, None, None))
    expected = dict(DEFAULTS, **rep.PRESET_OVERRIDES[preset])
    assert res["preset_applied"] == preset and res["settings"] == expected
    for k in DEFAULTS:
        want = f"preset:{preset}" if k in rep.PRESET_OVERRIDES[preset] else "default"
        assert res["source"][k] == want, (k, res["source"][k])
    assert [d["id"] for d in res["modules_applied"]] == list(sm.LEGACY_ALIAS[preset])


def test_none_applies_nothing_and_a_list_forces_its_modules():
    f = sm.flat_features(_profile(), None, None)
    none = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"}, modules_spec="none",
                       features=f)
    assert none["settings"] == DEFAULTS and none["preset_applied"] is None
    forced = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="short_query", features=f)
    assert forced["settings"]["classify_complete_min_qcov"] == 0.5
    assert forced["modules_applied"][0]["reason"].startswith("named in --modules")


def test_preset_override_sits_above_the_modules():
    f = sm.flat_features(_profile(copies=2.0), None, None)
    res = rep.resolve(DEFAULTS, auto_preset=None, modules_spec="auto", features=f,
                      preset_override="preset_short_peptide")
    assert res["settings"]["goi_relative_evalue"] == 1e-3                 # A1 fired
    assert res["source"]["classify_fragment_max_qcov"] == "preset:preset_short_peptide"
    assert res["preset_applied"] == "preset_short_peptide"


def test_a_value_the_user_set_survives_a_module_and_is_reported():
    f = sm.flat_features(_profile(copies=2.0), None, None)
    res = rep.resolve(dict(DEFAULTS, goi_relative_evalue="1e-6"), auto_preset=None,
                      modules_spec="auto", features=f, user_set={"goi_relative_evalue": "1e-6"})
    assert res["settings"]["goi_relative_evalue"] == "1e-6"
    assert res["source"]["goi_relative_evalue"] == "user_set"
    assert res["user_set_kept"] == ["goi_relative_evalue"]


def test_a_profile_preset_enters_as_user_set():
    """-profile preset_X sets params at launch: the resolver sees them as user-set values."""
    f = sm.flat_features(_profile(), None, None)
    preset_vals = rep.PRESET_OVERRIDES["preset_paralog_discrimination"]
    res = rep.resolve(dict(DEFAULTS, **preset_vals), auto_preset={"preset": "preset_single_copy"},
                      modules_spec="legacy", features=f, user_set=dict(preset_vals))
    for k, v in preset_vals.items():
        assert res["settings"][k] == v and res["source"][k] == "user_set"


def test_switching_auto_apply_off_applies_nothing_but_records_the_triggers():
    f = sm.flat_features(_profile(query_len=50), None, None)
    res = rep.resolve(DEFAULTS, auto_preset={"preset": "preset_short_peptide"}, modules_spec="auto",
                      features=f, disable_auto_apply=True)
    assert res["settings"] == DEFAULTS and res["modules_applied"] == []
    assert "B1" in {d["id"] for d in res["modules_would_fire"]}


def test_resolver_cli_with_sentinels_and_modules(tmp_path):
    (tmp_path / "defaults.json").write_text(json.dumps(DEFAULTS))
    (tmp_path / "hp.json").write_text(json.dumps(_profile(query_len=50)))
    sentinel = tmp_path / "NO_HOME_CONTEXT"
    sentinel.write_text("sentinel: none\n")
    out = tmp_path / "eff.json"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--auto_preset", str(tmp_path / "NO_AUTO_PRESET"), "--hit_profile", str(tmp_path / "hp.json"),
         "--goi_info", str(tmp_path / "missing.json"), "--home_context", str(sentinel),
         "--modules", "auto", "--defaults", str(tmp_path / "defaults.json"), "--output", str(out)],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    res = json.loads(out.read_text())
    assert res["modules_mode"] == "auto" and res["settings"]["classify_fragment_max_qcov"] == 0.3
    assert res["features"]["relatives_nearby"] is None
    assert "MODULES: mode auto; applied: B1 short_query" in proc.stdout
    bad = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
         "--modules", "bogus", "--defaults", str(tmp_path / "defaults.json"), "--output", str(out)],
        capture_output=True, text=True)
    assert bad.returncode == 2 and "unknown preset module" in bad.stderr


def test_report_reads_the_run_settings(tmp_path):
    eff = {"modules_mode": "legacy", "modules_applied": [], "modules_would_fire": [{"id": "B3"}],
           "preset_applied": "preset_single_copy", "apply_reason": "x", "settings": {"a": 1},
           "features": {"sp": True}}
    p = tmp_path / "effective_params.json"
    p.write_text(json.dumps(eff))
    rs = generate_report.load_run_settings(str(p))
    assert rs["modules_would_fire"] == [{"id": "B3"}] and rs["features"] == {"sp": True}
    (tmp_path / "NO_EFFECTIVE_PARAMS").write_text("sentinel: none\n")
    assert generate_report.load_run_settings(str(tmp_path / "NO_EFFECTIVE_PARAMS")) is None


# ============================================================ home context
CODON = {"M": "ATG", "K": "AAA", "I": "ATT", "F": "TTT", "L": "CTG", "P": "CCG", "V": "GTG",
         "A": "GCG", "G": "GGC", "E": "GAA", "R": "CGT", "S": "AGC", "C": "TGC", "N": "AAC",
         "Q": "CAG", "Y": "TAT", "T": "ACC", "D": "GAT", "H": "CAT", "W": "TGG"}


def _gene_dna(protein):
    return "".join(CODON[a] for a in protein) + "TAA"


@pytest.mark.skipif(shutil.which("blastp") is None and synvoy_specificity._blastp() is None,
                    reason="blastp not available")
def test_home_context_finds_a_family_member_among_the_neighbours(tmp_path):
    import random
    rnd = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    relative = "".join(c if rnd.random() < 0.6 else rnd.choice(aa) for c in LY6E)
    others = ["M" + "".join(rnd.choice(aa) for _ in range(150)) for _ in range(5)]
    proteins = others[:2] + [relative] + [LY6E] + others[2:]     # relative at rank -1
    spacer = "".join(rnd.choice("ACGT") for _ in range(3000))
    seq, gff, pos = spacer, ["##gff-version 3"], len(spacer)
    for i, prot in enumerate(proteins):
        dna = _gene_dna(prot)
        s, e = pos + 1, pos + len(dna)
        gff += [f"chr1\tt\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-g{i}",
                f"chr1\tt\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-g{i};Parent=gene-g{i}",
                f"chr1\tt\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-g{i};Parent=rna-g{i}"]
        seq += dna + spacer
        pos = len(seq)
    (tmp_path / "g.fna").write_text(">chr1\n" + seq + "\n")
    (tmp_path / "g.gff").write_text("\n".join(gff) + "\n")
    goi_line = [l for l in gff if l.startswith("chr1\tt\tCDS") and "cds-g3" in l][0].split("\t")
    (tmp_path / "goi.bed").write_text(f"chr1\t{int(goi_line[3]) - 1}\t{goi_line[4]}\tGOI_x|exon_1\t100\t+\n")
    shipped = synvoy_specificity.FAMILY_TEST_SEG
    ctx = hc.home_context(LY6E, str(tmp_path / "g.fna"), str(tmp_path / "g.gff"),
                          str(tmp_path / "goi.bed"), 10)
    assert ctx["goi_gene"] == "gene-g3" and ctx["n_neighbours_tested"] == 6
    assert [r["gene_id"] for r in ctx["relatives_nearby"]] == ["gene-g2"]
    assert ctx["relatives_nearby"][0]["rank"] == -1 and ctx["relative_adjacent"] is True
    assert ctx["scaffold_len"] == len(seq) and ctx["home_genome_bp"] == len(seq)
    assert synvoy_specificity.FAMILY_TEST_SEG == shipped       # restored after the test calls
    no_gff = hc.home_context(LY6E, str(tmp_path / "g.fna"), "NO_GFF", str(tmp_path / "goi.bed"), 10)
    assert no_gff["relatives_nearby"] is None and no_gff["dist_to_scaffold_end"] is not None


def test_home_context_without_blastp_reports_no_relatives(tmp_path, monkeypatch):
    (tmp_path / "g.fna").write_text(">chr1\n" + "ACGT" * 3000 + "\n")
    (tmp_path / "g.gff").write_text("chr1\tt\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=r1\n")
    (tmp_path / "goi.bed").write_text("chr1\t99\t300\tGOI_x|exon_1\t100\t+\n")
    monkeypatch.setattr(synvoy_specificity, "_blastp", lambda: None)
    ctx = hc.home_context("MKT", str(tmp_path / "g.fna"), str(tmp_path / "g.gff"),
                          str(tmp_path / "goi.bed"), 10)
    assert ctx["statistic"] == "fallback" and ctx["relatives_nearby"] is None
