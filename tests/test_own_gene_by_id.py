"""Review fix RF-2b (batch of 2026-10-15, S2-F2): q's own home gene in the specificity terms is
the gene the home step took as q's (goi_info.json gene_id, found by name or by position; batch
HG-3), no longer only the home protein near-identical to q (>= 95 % identity over >= 90 % of q).
That identity rule failed where the home annotation of q differs from q: the own gene became
p*. It stays the fallback without a gene id, and when the home step flagged its choice as
ambiguous (the gene it took may be a relative with q's own gene as the runner-up). The runner
names the gene for the home GOI only; the hull rescue reads goi_info.json too. All synthetic."""
import json
import logging
import os
import random
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
import synvoy_specificity as spec  # noqa: E402
from sequence_utils import sw_bits  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _protein(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def _mutate(rng, seq, frac):
    return "".join(rng.choice(AA) if rng.random() < frac else c for c in seq)


def _case(seed=11):
    """q (200 aa) whose home model is partial (its first 150 residues: 75 % of q, below the
    identity rule's 90 %), a relative at about 60 % identity, a close relative at about 92 %,
    unrelated genes, and a target ortholog x of q."""
    rng = random.Random(seed)
    q = _protein(rng, 200)
    records = [("gene-Q", q[:150]), ("gene-P", _mutate(rng, q, 0.4)),
               ("gene-C", _mutate(rng, q, 0.08))] + [(f"gene-{i}", _protein(rng, 200)) for i in range(6)]
    x = _mutate(rng, q, 0.1)
    return q, records, x


def _write(folder, records):
    with open(os.path.join(folder, spec.PROTEOME_FILE), "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


# ------------------------------------------------------------------ the module
def test_a_partial_home_model_is_still_q_s_own_gene():
    q, records, x = _case()
    records = [r for r in records if r[0] != "gene-C"]
    hp = spec.HomeProteome(records, {})
    assert hp.own_gene("GOI_q", q) is None                    # the identity rule misses it ...
    before = spec.specificity(hp, "GOI_q", q, x, sw_bits(q, x), None)
    assert before["best_other"] == "gene-Q"                   # ... so the own gene is p*
    after = spec.specificity(hp, "GOI_q", q, x, sw_bits(q, x), None, own_gene_id="gene-Q")
    assert hp.own_gene("GOI_q", q, "gene-Q") == "gene-Q"
    assert after["best_other"] == "gene-P" and after["margin"] > before["margin"] + 50


def test_an_ambiguous_choice_keeps_the_identity_rule_and_warns(caplog):
    # q's own gene is annotated whole; the home step took the close relative gene-C by position
    # and flagged it (q's own gene as the runner-up within 0.05).
    q, records, x = _case(12)
    records = [("gene-Q", q) if rid == "gene-Q" else (rid, s) for rid, s in records]
    hp = spec.HomeProteome(records, {})
    info = {"found_by": "position", "gene_id": "gene-C", "home_gene_ambiguous": True,
            "runner_up_gene_id": "gene-Q", "runner_up_score_fraction": 0.97}
    spec._WARNED.discard("ambiguous:gene-C")
    with caplog.at_level(logging.WARNING, logger=spec.__name__):
        assert spec.own_gene_id_from_goi_info(info) is None
    assert "ambiguous" in caplog.text and "gene-Q" in caplog.text
    kept = spec.specificity(hp, "GOI_q", q, x, sw_bits(q, x), None, own_gene_id=None)
    assert hp.own_gene("GOI_q", q) == "gene-Q" and kept["best_other"] == "gene-C"
    # Taken blindly, the flagged gene would have made q's own gene p*.
    blind = spec.specificity(hp, "GOI_q", q, x, sw_bits(q, x), None, own_gene_id="gene-C")
    assert blind["best_other"] == "gene-Q" and blind["margin"] < kept["margin"]


def test_the_gene_comes_only_from_a_name_or_position_match(tmp_path):
    get = spec.own_gene_id_from_goi_info
    assert get({"found_by": "name", "gene_id": "gene-Q", "home_gene_ambiguous": False}) == "gene-Q"
    assert get({"found_by": "position", "gene_id": "gene-Q"}) == "gene-Q"
    for info in ({"found_by": "miniprot", "gene_id": None}, {"found_by": "none"},
                 {"found_by": "position_3mer", "gene_id": "gene-Q"},
                 {"method": "gff_annotation", "num_exons": 2},          # written before HG-3
                 {"found_by": "name", "gene_id": ""}, None, [], "not a file"):
        assert get(info) is None
    path = tmp_path / "goi_info.json"
    path.write_text(json.dumps({"found_by": "name", "gene_id": "gene-Q"}))
    assert get(str(path)) == "gene-Q"
    (tmp_path / "bad.json").write_text("{")
    assert get(str(tmp_path / "bad.json")) is None


def test_a_gene_missing_from_the_proteome_falls_back_with_a_warning(caplog):
    q, records, _x = _case(13)
    records = [("gene-Q2", q) if rid == "gene-Q" else (rid, s) for rid, s in records]
    hp = spec.HomeProteome(records, {})
    spec._WARNED.discard("own:gene-Q")
    with caplog.at_level(logging.WARNING, logger=spec.__name__):
        assert hp.own_gene("GOI_q", q, "gene-Q") == "gene-Q2"      # the identity rule
    assert "not in the home proteome" in caplog.text


def test_q_s_own_locus_follows_the_named_gene():
    q, records, _x = _case(14)
    records = [r for r in records if r[0] != "gene-C"]
    coords = {"gene-0": ("c", 100, 200), "gene-Q": ("c", 1000, 1600), "gene-1": ("c", 3000, 3500),
              "gene-P": ("d", 100, 900), "gene-2": ("d", 2000, 2500)}
    hp = spec.HomeProteome(records, coords)
    assert hp.own_neighbourhood("GOI_q", q, 1) is None             # identity rule: no own gene
    assert hp.own_neighbourhood("GOI_q", q, 1, "gene-Q") == {"gene-0", "gene-Q", "gene-1"}


# ------------------------------------------------------------------ the callers
def _cand(mid, seq, parent):
    attrs = f"ID={parent}|g_b0_{mid};Name={parent};SynVoy_Parent={parent};SynVoyRole=goi;Confidence=MEDIUM"
    return {"start": 100, "end": 900, "record": {"id": f"{parent}|g_b0_{mid}", "seq": seq},
            "gff": [f"chr1\texon_annotation\tmRNA\t100\t900\t50\t+\t.\t{attrs}"]}


def test_the_runner_names_the_own_gene_for_the_home_goi_only(tmp_path, monkeypatch):
    q, records, x = _case(15)
    records = [r for r in records if r[0] != "gene-C"]
    _write(str(tmp_path), records)
    seen = []
    real = spec.specificity

    def spy(*a, **kw):
        seen.append((a[1], kw.get("own_gene_id")))
        return real(*a, **kw)
    monkeypatch.setattr(isr.synvoy_specificity, "specificity", spy)
    decoy = q[::-1]
    by_parent = {"GOI_q": {"id": "GOI_q", "seq": q}, "GOI_q~rev": {"id": "GOI_q~rev", "seq": decoy}}

    def run(gene_id):
        isr._HOME_PROTEOME_CACHE.clear()
        args = types.SimpleNamespace(home_db_dir=str(tmp_path), home_flanking_ids=[],
                                     home_goi_ids={"GOI_q", "GOI_q|exon_1"}, home_goi_gene_id=gene_id)
        cands = [_cand("l1", x, "GOI_q"), _cand("l2", _mutate(random.Random(1), decoy, 0.1), "GOI_q~rev")]
        isr.score_goi_candidates(cands, by_parent)
        assert isr.annotate_goi_specificity(cands, by_parent, args) == 2
        return [isr._parse_gff_attr_string(c["gff"][0].split("\t")[8]) for c in cands]

    with_id = run("gene-Q")
    assert seen == [("GOI_q", "gene-Q"), ("GOI_q~rev", None)]   # the decoy keeps the identity rule
    assert with_id[0]["BestOtherHome"] == "gene-P"
    seen.clear()
    without = run(None)
    assert seen == [("GOI_q", None), ("GOI_q~rev", None)]
    assert without[0]["BestOtherHome"] == "gene-Q"                # the old result: the own gene as p*
    assert float(with_id[0]["SpecMargin"]) > float(without[0]["SpecMargin"])
    isr._HOME_PROTEOME_CACHE.clear()


def test_the_hull_judge_passes_the_own_gene(monkeypatch):
    q, records, x = _case(16)
    records = [r for r in records if r[0] != "gene-C"]
    hp = spec.HomeProteome(records, {})
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", lambda qq, s, sp: (1e-40, 300.0))
    rows = [["w", "miniprot", "CDS", "1", str(3 * len(x)), ".", "+", "0", "Parent=m"]]
    _c, _k, attrs_id, _b = hull.judge_hull_model("q", q, x, rows, 10 ** 7, hp, coverage=0.95,
                                                 own_gene_id="gene-Q")
    _c, _k, attrs_rule, _b = hull.judge_hull_model("q", q, x, rows, 10 ** 7, hp, coverage=0.95)
    assert "BestOtherHome=gene-P" in attrs_id and "BestOtherHome=gene-Q" in attrs_rule
