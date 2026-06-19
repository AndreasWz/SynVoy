"""Gap-bounded fallback chaining (_longest_monotonic_query_chain).

Regression for the "two seqs far apart called one gene" failure: a melittin
(~70 aa) run produced GOI_Melt fallback_hit_span mRNAs spanning 27-60 kb that
were really two short, weak, unrelated hits (e.g. 17 aa @ 3,003,077 + 12 aa @
3,063,163, 60 kb apart, on a "protein sprouty" gene) stitched into one model.

The chain must only join two hits as consecutive exons of one gene when the
genomic gap between them is within max_gap (the max-intron distance). When the
gap bound leaves only single-hit chains it must return the single best hit.
"""
import os
import sys

BIN = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402

chain = isr._longest_monotonic_query_chain


def _hit(gstart, gend, qstart, qend, bits, alnlen, strand="+"):
    return {"gstart": gstart, "gend": gend, "qstart": qstart, "qend": qend,
            "bits": bits, "alnlen": alnlen, "strand": strand}


# --- the actual melittin b18_l1 failure: two short hits 60 kb apart --------- #
MELITTIN_FAR = [
    _hit(3003077, 3003128, 1, 17, 30.0, 17),
    _hit(3063163, 3063199, 55, 67, 22.0, 12),
]


def test_no_max_gap_preserves_legacy_behaviour():
    # Without a bound the two far-apart hits still chain (old behaviour).
    out = chain(MELITTIN_FAR, "+")
    assert len(out) == 2


def test_far_apart_hits_not_chained():
    out = chain(MELITTIN_FAR, "+", max_gap=20000)
    assert len(out) == 1, "60 kb gap must not be stitched into one gene"


def test_broken_chain_keeps_best_single_hit():
    # The stronger hit (bits 30) must win over the lower-gstart-but-weaker one
    # only if it were weaker — here the first hit is also the strongest.
    out = chain(MELITTIN_FAR, "+", max_gap=20000)
    assert out[0]["bits"] == 30.0
    # Flip the strengths: best single hit should be the downstream one now.
    flipped = [_hit(3003077, 3003128, 1, 17, 10.0, 17),
               _hit(3063163, 3063199, 55, 67, 40.0, 30)]
    out2 = chain(flipped, "+", max_gap=20000)
    assert len(out2) == 1 and out2[0]["bits"] == 40.0


def test_compact_multi_exon_still_chains():
    # A genuine 2-exon locus with a 300 bp intron must survive the bound.
    compact = [_hit(1000, 1051, 1, 17, 30.0, 17),
               _hit(1351, 1450, 18, 50, 40.0, 33)]
    out = chain(compact, "+", max_gap=20000)
    assert len(out) == 2


def test_gap_just_under_bound_chains():
    # Two exons exactly max_gap apart are allowed (boundary inclusive).
    pair = [_hit(1000, 1100, 1, 20, 30.0, 20),
            _hit(1100 + 20000, 1100 + 20100, 21, 40, 30.0, 20)]
    out = chain(pair, "+", max_gap=20000)
    assert len(out) == 2


def test_minus_strand_far_apart_not_chained():
    # On '-' strand query centers decrease with genomic position; gap bound
    # must apply the same way.
    minus = [_hit(8343568, 8343628, 55, 67, 30.0, 20, "-"),
             _hit(8370962, 8371016, 1, 18, 22.0, 18, "-")]
    out = chain(minus, "-", max_gap=20000)
    assert len(out) == 1, "27 kb gap on - strand must not be stitched"


def test_single_hit_unchanged():
    one = [_hit(1000, 1100, 1, 20, 30.0, 20)]
    assert chain(one, "+", max_gap=20000) == one


# --- data-driven max-gap: derive from the home GOI's own largest intron ----- #
import json
import tempfile

derive = isr.derive_goi_fallback_max_gap


def _info(**kw):
    fh = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(kw, fh)
    fh.close()
    return fh.name


def test_derive_melittin_floored():
    # Home melittin: one ~1.15 kb intron. 1148 x 5 = 5740 -> floored to 10 kb.
    # The real run's 27/60 kb spurious joins are well above this and die.
    gap, reason = derive(_info(max_intron_bp=1148, genomic_span_bp=1361), 5.0, 10000, 20000)
    assert gap == 10000, reason


def test_derive_large_gene_scales_up():
    # A gene with a 15 kb home intron must be allowed big introns in targets.
    gap, _ = derive(_info(max_intron_bp=15000), 5.0, 10000, 20000)
    assert gap == 75000


def test_derive_single_exon_uses_floor():
    gap, _ = derive(_info(max_intron_bp=0), 5.0, 10000, 20000)
    assert gap == 10000


def test_derive_missing_file_falls_back_to_max_intron():
    gap, _ = derive("/no/such/file.json", 5.0, 10000, 20000)
    assert gap == 20000


def test_derive_null_field_falls_back():
    gap, _ = derive(_info(max_intron_bp=None), 5.0, 10000, 20000)
    assert gap == 20000


def test_derive_none_path_falls_back():
    gap, _ = derive(None, 5.0, 10000, 20000)
    assert gap == 20000


# --- home GOI footprint (annotate_goi_exons.compute_genomic_footprint) ------ #
import importlib

age = importlib.import_module("annotate_goi_exons")


def test_footprint_two_exon_melittin():
    # Real home melittin exon coords (NC_037641.1):
    exons = [
        {"gstart": 12398900, "gend": 12399017},
        {"gstart": 12397656, "gend": 12397752},
    ]
    span, max_intron = age.compute_genomic_footprint(exons)
    assert span == 12399017 - 12397656  # 1361
    assert max_intron == 12398900 - 12397752  # 1148


def test_footprint_single_exon_zero_intron():
    span, max_intron = age.compute_genomic_footprint([{"gstart": 100, "gend": 400}])
    assert span == 300 and max_intron == 0


def test_footprint_coords_fallback_and_no_coords():
    span, max_intron = age.compute_genomic_footprint([{"coords": (10, 60)}, {"coords": (200, 250)}])
    assert span == 240 and max_intron == 140
    assert age.compute_genomic_footprint([{"id": "x"}]) == (None, None)
