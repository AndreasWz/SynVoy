"""Item 37 (fix phase, 2 Nov): on equal bits the first subject is the judged protein.

score_goi_candidates took ``max((sw_bits(home, s), s) for s in subjects)``, so between two
subjects of equal bits the protein's letters decided which one Bits and the specificity read.
Since the same-frame join (fix-phase item 1) a joined piece ties with its stronger hit, and on
a constructed ORF one residue flips the label between HIGH and AMBIGUOUS. Now the first subject
in ``label_subjects`` order wins: the label protein's pieces in transcript order, then the hit
regions in theirs. rescue_goi_hull.judge_hull_model had the same line over the pieces of one
protein; the same rule there.
"""
import os
import random
import sys

import pytest

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

pytest.importorskip("parasail")

import iterative_search_runner as isr  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
import synvoy_label  # noqa: E402
import synvoy_specificity as spec  # noqa: E402
from sequence_utils import sw_bits  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _random(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice(AA) for _ in range(n))


def _mutate(seq, frac, seed):
    rng = random.Random(seed)
    out = list(seq)
    for i in rng.sample(range(len(seq)), int(len(seq) * frac)):
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return "".join(out)


Q = _random(220, 1)                          # the home GOI
P = _random(260, 7)                          # another home gene, unrelated to q
HIT = "M" + _mutate(Q[1:160], 0.35, 3)       # a diverged hit of q (first letter between A and Y)
LINKER = _random(40, 11)


def _joined(first):
    """The ORF's joined piece: a domain close to gene-p, a linker, then the hit of q. Only its
    first residue differs between the two cases."""
    return first + _mutate(P, 0.2, 5)[1:] + LINKER + HIT


def _proteome(folder):
    recs = [("gene-q", Q), ("gene-p", P)] + [(f"gene-u{i}", _random(300, 100 + i)) for i in range(20)]
    with open(os.path.join(folder, spec.PROTEOME_FILE), "w") as fh:
        for rid, seq in recs:
            fh.write(f">{rid}\n{seq}\n")
    return spec.HomeProteome.load(folder)


def _cand(subjects):
    mrna = ("c1\tfallback_hits\tmRNA\t1\t1500\t60.0\t+\t.\tID=GOI_q|g_b0_l1_fallback;"
            "SynVoy_Parent=GOI_q;SynVoyRole=goi;EvidenceType=fallback_hit_span")
    return {"start": 1, "end": 1500, "gff": [mrna], "record": {"seq": subjects[0]},
            "label_subjects": list(subjects)}


def test_the_case_is_a_tie():
    for first in "AY":
        assert sw_bits(Q, _joined(first)) == sw_bits(Q, HIT)


def test_under_the_old_rule_the_letter_decides():
    picked = {first: max((sw_bits(Q, s), s) for s in (_joined(first), HIT))[1] for first in "AY"}
    assert picked["A"] == HIT and picked["Y"] == _joined("Y")


@pytest.mark.parametrize("first", ["A", "Y"])
def test_the_first_subject_is_judged_whatever_its_letters(first):
    cand = _cand([_joined(first), HIT])
    assert isr.score_goi_candidates([cand], {"GOI_q": {"seq": Q}}) == 1
    assert cand["judged_protein"] == _joined(first)
    assert cand["bits"] == sw_bits(Q, HIT)


def test_a_hit_region_listed_first_is_judged_first():
    cand = _cand([HIT, _joined("Y")])
    isr.score_goi_candidates([cand], {"GOI_q": {"seq": Q}})
    assert cand["judged_protein"] == HIT


def test_a_strictly_better_subject_still_wins_from_behind():
    cand = _cand([HIT, Q[:200]])
    isr.score_goi_candidates([cand], {"GOI_q": {"seq": Q}})
    assert cand["judged_protein"] == Q[:200]


def test_the_label_no_longer_depends_on_one_residue(tmp_path):
    """Through the real specificity and judge: under the old rule 'A' read HIGH (the hit is
    judged; gene-p is nowhere near it) and 'Y' AMBIGUOUS not_specific (the joined piece is
    judged; gene-p explains it far better than q). Now both read the joined piece's label."""
    proteome = _proteome(str(tmp_path))
    labels = {}
    for first in "AY":
        cand = _cand([_joined(first), HIT])
        isr.score_goi_candidates([cand], {"GOI_q": {"seq": Q}})
        s = spec.specificity(proteome, "GOI_q", Q, cand["judged_protein"], cand["bits"], None)
        labels[first] = synvoy_label.judge(1e-30, "slot", True, s["margin"], s["real_competitor"], None, 0.9)
    assert labels["A"] == labels["Y"] == ("AMBIGUOUS", "not_specific")
    old = {}
    for first in "AY":
        bits, x = max((sw_bits(Q, s), s) for s in (_joined(first), HIT))
        s = spec.specificity(proteome, "GOI_q", Q, x, bits, None)
        old[first] = synvoy_label.judge(1e-30, "slot", True, s["margin"], s["real_competitor"], None, 0.9)[0]
    assert old == {"A": "HIGH", "Y": "AMBIGUOUS"}


@pytest.mark.parametrize("letters", ["AY", "YA"])
def test_the_hull_judges_its_first_piece_on_a_tie(monkeypatch, letters):
    """Two pieces of one hull protein with equal bits: the first in transcript order is judged."""
    seen = {}
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", lambda q, s, sp: (1e-30, 200.0))

    def fake_spec(proteome, qid, q, x, bits, hood, n_flanking=0, own_gene_id=None):
        seen["x"] = x
        return {"margin": 50.0, "real_competitor": False}
    monkeypatch.setattr(hull.synvoy_specificity, "specificity", fake_spec)
    monkeypatch.setattr(hull.synvoy_specificity, "gff_attrs", lambda spec_: "SpecMargin=50")
    one, two = letters[0] + LINKER + HIT, letters[1] + LINKER + HIT
    assert sw_bits(Q, one) == sw_bits(Q, two)
    rows = [["w", "miniprot", "CDS", "1", "900", ".", "+", "0", "Parent=m"]]
    hull.judge_hull_model("q", Q, one + "X" + two, rows, 10 ** 7, object(),
                          label_protein=one + "*" + two, coverage=0.9)
    assert seen["x"] == one


def test_the_hull_without_pieces_judges_nothing():
    conf, _key, _attrs, bits = hull.judge_hull_model("q", Q, "", [], 10 ** 7, None, label_protein="")
    assert bits == 0.0 and conf == "LOW"


def test_the_genome_wide_path_takes_the_first_piece_on_a_tie():
    """The same rule in genome_wide_goi_pass (opt-in): its subjects were a sorted set, and max over (bits, piece) let
    the letters decide. Now the place's hits in order, each hit's pieces in transcript order, each piece once."""
    src = open(os.path.join(BIN, "iterative_search_runner.py")).read()
    assert "subjects = list(dict.fromkeys(piece for h in p['hit_list'] for piece in _hit_pieces(genome_seqs, h)))" in src
    assert "bits, best = first_best_subject(q, subjects)" in src
    assert "key=lambda s: (sw_bits(q, s), s)" not in src
    for first in "AY":
        assert isr.first_best_subject(Q, [_joined(first), HIT]) == (sw_bits(Q, HIT), _joined(first))
    assert isr.first_best_subject(Q, [HIT, Q[:200]])[1] == Q[:200]
