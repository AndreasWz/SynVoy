"""Fix-phase item 8: the hit-chain fallback's span cap (iterative_search_runner.fallback_span_limit).

The old cap was max(100 kb, 600 bp per aa) below 150 aa and max(5 kb, 150 bp per aa) from 150 aa on, so it fell
from 100 kb to 22.5 kb at 150 aa, the 150 was a literal, and the cap ignored the GOI's own span at home.
The counter that shows a refusal was 0 on the Candida pilot and the family seeds, so the rule rests on these
tests.
"""
import os
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import iterative_search_runner as isr  # noqa: E402
import synvoy_scale  # noqa: E402


def _old_cap(query_len):
    if query_len < 150:
        return max(100_000, query_len * 3 * 200)
    return max(5000, query_len * 3 * 50)


def test_the_old_cap_jumped():
    assert _old_cap(149) == 100_000 and _old_cap(150) == 22_500


def test_the_cap_is_continuous_in_length():
    """No step larger than one residue's 150 bp, at any length."""
    caps = [isr.fallback_span_limit(n) for n in range(1, 5001)]
    assert max(b - a for a, b in zip(caps, caps[1:])) <= 150
    assert all(b >= a for a, b in zip(caps, caps[1:]))


def test_no_chain_the_old_cap_kept_is_refused_now():
    """Only refusals can lift: a chain the old cap kept goes to the judge as before."""
    for n in range(1, 5001):
        assert isr.fallback_span_limit(n) >= _old_cap(n), n


@pytest.mark.parametrize("n,cap", [(30, 100_000), (149, 100_000), (150, 100_000), (666, 100_000),
                                   (667, 100_050), (1000, 150_000), (5000, 750_000)])
def test_values(n, cap):
    assert isr.fallback_span_limit(n) == cap


def test_the_home_span_times_the_expansion_is_the_floor():
    """As an anchor: 1.5 x home span x expansion (synvoy_scale.gene_gap_limit)."""
    assert isr.fallback_span_limit(200, home_span=300_000, expansion=2.0) == 900_000
    assert isr.fallback_span_limit(200, home_span=300_000, expansion=0.5) == 225_000
    # A compact gene in a compact target leaves the length arm.
    assert isr.fallback_span_limit(200, home_span=2_000, expansion=1.0) == 100_000
    # The same rule as an anchor with this span.
    geometry = {"genes": {"g": {"span": 300_000}}}
    assert synvoy_scale.gene_gap_limit(geometry, "g", 2.0, 100_000) == isr.fallback_span_limit(200, 300_000, 2.0)


@pytest.mark.parametrize("span,expansion", [(0, 1.0), (None, 1.0), (0, None), (5_000, None)])
def test_without_a_home_span_or_expansion_the_length_arm_stands(span, expansion):
    assert isr.fallback_span_limit(400, home_span=span, expansion=expansion) == 100_000


def _hits(a, b):
    return [{"gstart": a, "gend": a + 100}, {"gstart": b - 100, "gend": b}]


def _args(**kw):
    counts = {}
    return types.SimpleNamespace(limit_counts=counts, **kw), counts


def test_too_wide_reads_the_targets_span_and_counts(monkeypatch):
    calls = []
    monkeypatch.setattr(isr, "count_limit", lambda args, name, n=1: calls.append(name))
    args, _ = _args(goi_home_span=100_000, expansion=2.0)
    assert not isr.fallback_span_too_wide(args, "GOI_x", _hits(0, 300_000), 200)
    assert calls == []
    assert isr.fallback_span_too_wide(args, "GOI_x", _hits(0, 300_001), 200)
    assert calls == ["fallback_span_cap"]


def test_a_decoy_has_the_same_cap_and_is_not_counted(monkeypatch):
    calls = []
    monkeypatch.setattr(isr, "count_limit", lambda args, name, n=1: calls.append(name))
    args, _ = _args(goi_home_span=100_000, expansion=2.0)
    assert isr.fallback_span_too_wide(args, "GOI_x~rev", _hits(0, 300_001), 200)
    assert not isr.fallback_span_too_wide(args, "GOI_x~shuf1", _hits(0, 300_000), 200)
    assert calls == []


def test_args_without_a_span_keep_the_length_arm(monkeypatch):
    monkeypatch.setattr(isr, "count_limit", lambda *a, **k: None)
    args = types.SimpleNamespace()
    assert not isr.fallback_span_too_wide(args, "GOI_x", _hits(0, 100_000), 200)
    assert isr.fallback_span_too_wide(args, "GOI_x", _hits(0, 100_001), 200)
