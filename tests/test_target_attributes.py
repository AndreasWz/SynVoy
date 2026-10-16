"""Attributes and report categories that change no call (preset plan of 2026-10-15, PR-2).

* ContextTruncated  a block's search window ran past a scaffold end (left/right/both)
* ExonsKept         exons of a gene model / exons of the home gene
* CysKept           home-protein cysteines aligned to a cysteine of the model protein
* hits/<genome>.target_policy.tsv  median flanking identity, blocks, the best block's context
* report: absent_context_complete / unknown_context_truncated / absent_no_block,
  nearest_paralog_only, target_policy, run_settings.flags
"""
from __future__ import annotations

import csv
import json
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import generate_report as gr  # noqa: E402
import iterative_search_runner as isr  # noqa: E402

LY6E = ("MKIFLPVLLAALLGVERASSLMCFSCLNQKSNLYCLKPTICSDQDNYCVTVSASAGIGNLVTFGHSLSKTCSPACPIPEGVNVG"
        "VASMGISCCQSFLCNFSAADGGLRASVTLLGAGLLLSLLPALLRFGP")


# ============================================================ ContextTruncated
@pytest.mark.parametrize("start,end,pad,slen,want", [
    (50_000, 60_000, 20_000, 1_000_000, ""),
    (10_000, 20_000, 20_000, 1_000_000, "left"),
    (980_000, 990_000, 20_000, 1_000_000, "right"),
    (5_000, 15_000, 20_000, 30_000, "both"),
    (20_000, 980_000, 20_000, 1_000_000, ""),         # exactly at the ends: not cut
])
def test_context_truncation(start, end, pad, slen, want):
    assert isr.context_truncation(start, end, pad, slen) == want


def _cand(evidence, parent="GOI_LY6E", exons="3", protein=LY6E):
    attrs = (f"ID={parent}|G_b0_l1;SynVoy_Parent={parent};SynVoyRole=goi;EvidenceType={evidence};"
             f"Exons={exons};Confidence=LOW")
    return {"gff": [f"chr1\tSynVoy\tmRNA\t100\t900\t.\t+\t.\t{attrs}"],
            "record": {"id": f"{parent}|G_b0_l1", "seq": protein}, "start": 100, "end": 900}


def _attrs(cand):
    return isr._parse_gff_attr_string(cand["gff"][0].split("\t")[8])


def test_model_extras_on_a_gene_model_in_a_cut_window():
    by_parent = {"GOI_LY6E": {"seq": LY6E}}
    mutated = LY6E.replace("CFSC", "SFSC", 1)            # one cysteine lost
    model, chain = _cand("exon_annotation", protein=mutated), _cand("fallback_hit_span")
    isr.annotate_model_extras([model, chain], by_parent, 5, "right")
    a = _attrs(model)
    assert (a["ContextTruncated"], a["ExonsKept"], a["CysKept"]) == ("right", "3/5", "9/10")
    b = _attrs(chain)
    assert b["ContextTruncated"] == "right" and "ExonsKept" not in b and "CysKept" not in b


def test_model_extras_add_nothing_in_a_complete_window_without_home_exons():
    model = _cand("exon_annotation", protein="MKTAYIAKQRQISFVKSHFSRQ")
    before = model["gff"][0]
    isr.annotate_model_extras([model], {"GOI_LY6E": {"seq": "MKTAYIAKQRQISFVKSHFSRQ"}}, None, "")
    assert model["gff"][0] == before                      # no Cys at home, no exon count, no cut


def test_home_exon_count_ignores_tandem_copies_of_the_home_model(tmp_path):
    info = {"method": "hit_annotation",
            "exons": [{"id": "GOI_q|exon_1"}, {"id": "GOI_q|exon_2"}, {"id": "GOI_copy_1"}]}
    (tmp_path / "goi_info.json").write_text(json.dumps(info))
    assert isr.home_exon_terms(*isr.read_home_exons(str(tmp_path / "goi_info.json"))) == (2, None)
    for method in ("tandem_duplication", "full_protein_only"):          # no gene model at home
        assert isr.home_exon_terms(method, [40, 50])[0] is None
    assert isr.home_exon_terms(*isr.read_home_exons(str(tmp_path / "missing.json"))) == (None, None)


# ============================================================ target_policy.tsv
def _policy(tmp_path, ctx):
    args = types.SimpleNamespace(output_dir=str(tmp_path))
    isr._write_target_policy(args, "locus_1_", "G.fna", 61.25, 14, ctx)
    with open(tmp_path / "hits" / "locus_1_G.fna.target_policy.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_target_policy_names_the_best_block_and_its_context(tmp_path):
    rows = _policy(tmp_path, {"block_esyn": {0: 0.5, 1: 2e-6, 2: 0.03},
                              "block_context": {0: "", 1: "left", 2: ""}})
    assert len(rows) == 1
    r = rows[0]
    assert (r["median_flank_identity"], r["n_flanking_genes"], r["blocks_searched"]) == ("61.2", "14", "3")
    assert (r["significant_blocks"], r["best_block_esyn"], r["best_block_context"]) == ("2", "2e-06", "left")
    assert r["truncated_blocks"] == "1"


def test_target_policy_exists_without_blocks(tmp_path):
    args = types.SimpleNamespace(output_dir=str(tmp_path))
    isr._write_target_policy(args, "", "G.fna", None, 0)
    with open(tmp_path / "hits" / "G.fna.target_policy.tsv", newline="") as fh:
        (r,) = list(csv.DictReader(fh, delimiter="\t"))
    assert (r["median_flank_identity"], r["best_block_context"], r["blocks_searched"]) == ("NA", "NA", "0")


# ============================================================ report
def _write_policy(path, genome, ctx, esyn="1e-08"):
    with open(path, "w") as fh:
        fh.write("\t".join(isr.TARGET_POLICY_COLUMNS) + "\n")
        fh.write("\t".join([genome, "55.0", "12", "4", "2", esyn, ctx, "1" if ctx not in ("complete", "NA") else "0"]) + "\n")


def test_report_splits_genomes_without_a_counted_call(tmp_path):
    for i, (g, ctx) in enumerate([("A.fna", "complete"), ("B.fna", "left"), ("C.fna", "NA"),
                                  ("D.fna", "complete")]):
        _write_policy(tmp_path / f"locus_1__{g}{gr.TARGET_POLICY_SUFFIX}", g, ctx)
    policy = gr.collect_target_policy(sorted(str(p) for p in tmp_path.iterdir()))
    assert set(policy) == {"A", "B", "C", "D"}
    records = [
        {"genome": "D", "confidence": "MEDIUM", "goi_class": "probable_goi", "eeff": "1e-5"},
        {"genome": "A", "confidence": "AMBIGUOUS", "basis": "not_specific", "eeff": "1e-6",
         "goi_class": "syntenic_candidate_unconfirmed"},
        {"genome": "A", "confidence": "LOW", "basis": "not_significant", "eeff": "3"},
        {"genome": "B", "confidence": "HIGH", "goi_class": "paralog_not_goi", "eeff": "1e-20"},
    ]
    for i, r in enumerate(records):
        r.update(chrom="c1", start=1000 * (i + 1), end=1000 * (i + 1) + 600)
    cat = gr.genome_absence_categories(records, policy, {"paralog_not_goi"})
    assert cat == {"nearest_paralog_only": ["A"], "absent_context_complete": ["A"],
                   "unknown_context_truncated": ["B"], "absent_no_block": ["C"]}


def test_target_policy_of_several_loci_keeps_the_complete_view(tmp_path):
    _write_policy(tmp_path / f"locus_1__X.fna{gr.TARGET_POLICY_SUFFIX}", "X.fna", "right", "1e-9")
    _write_policy(tmp_path / f"locus_2__X.fna{gr.TARGET_POLICY_SUFFIX}", "X.fna", "complete", "0.04")
    policy = gr.collect_target_policy(sorted(str(p) for p in tmp_path.iterdir()))
    assert policy["X"]["best_block_context"] == "complete"


def test_run_settings_carry_the_module_flags(tmp_path):
    eff = {"modules_mode": "legacy",
           "modules_would_fire": [{"id": "A2"}, {"id": "C2"}],
           "features": {"loci": [{"chrom": "c1", "qcov": 1.0}, {"chrom": "c2", "qcov": 0.8}],
                        "dist_to_scaffold_end": 12000, "telomeric_repeat_at_end": True}}
    p = tmp_path / "effective_params.json"
    p.write_text(json.dumps(eff))
    flags = gr.load_run_settings(str(p))["flags"]
    assert flags["dispersed_family_loci"][1]["qcov"] == 0.8 and not flags["domain_shared"]
    assert flags["home_contig_end"] and flags["home_subtelomeric"]
    assert flags["dist_to_scaffold_end"] == 12000


def test_report_module_stages_the_target_policy():
    text = (ROOT / "modules" / "generate_report.nf").read_text()
    assert "*.target_policy.tsv" in text
