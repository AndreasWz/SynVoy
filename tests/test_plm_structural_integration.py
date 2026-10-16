#!/usr/bin/env python3
"""
test_plm_structural_integration.py — guards for the experimental ML search
layers (PLM embedding search + ESMFold/Foldseek structural search) and, in
particular, the soundness fix that stops a SYNTHETIC ML similarity (ProtT5
cosine / Foldseek fident scaled to look like sequence stats) from being read by
the identity/coverage classifier.

These tests are numpy-only — they do NOT require torch, transformers, or
foldseek, so they run in plain synvoy_env. The torch-dependent code paths
(embedding, folding) are exercised separately on the cluster.

Run with the synvoy_env interpreter (see CLAUDE.md):
    /home/faw/miniforge3/envs/synvoy_env/bin/python3 -m pytest tests/test_plm_structural_integration.py -v
"""
import os
import sys

import numpy as np
import pytest

BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")
if BIN not in sys.path:
    sys.path.insert(0, BIN)

import iterative_search_runner as isr  # noqa: E402
import plm_search  # noqa: E402
import structural_search  # noqa: E402


# ---------------------------------------------------------------------------
# Soundness fix: _resolve_fallback_signals must keep synthetic ML stats out of
# the identity channel.
# ---------------------------------------------------------------------------

def test_resolve_signals_noop_when_no_ml_hits():
    """With only real sequence hits (PLM off), identity == mean(pident) and the
    locus is never flagged ml_only — i.e. a perfect no-op vs the old behavior."""
    hits = [{"pident": 80.0}, {"pident": 60.0, "method": None}]
    seq_id, emb, struct, ml_only = isr._resolve_fallback_signals(hits)
    assert seq_id == pytest.approx(70.0)
    assert emb is None and struct is None
    assert ml_only is False


def test_resolve_signals_ml_only_zeroes_identity():
    """An ML-only locus must NOT pass its synthetic cosine (85) as identity; the
    real signal is surfaced as embedding_similarity and ml_only is True."""
    hits = [{"pident": 85.0, "method": "plm_embedding", "embedding_similarity": 0.72}]
    seq_id, emb, struct, ml_only = isr._resolve_fallback_signals(hits)
    assert seq_id == 0.0
    assert emb == pytest.approx(0.72)
    assert struct is None
    assert ml_only is True


def test_resolve_signals_mixed_keeps_real_identity_plus_booster():
    """A mixed locus keeps the honest real identity and carries the ML signal
    as a booster (not folded into identity)."""
    hits = [
        {"pident": 40.0, "method": None},
        {"pident": 90.0, "method": "foldseek_structural", "structural_similarity": 0.61},
    ]
    seq_id, emb, struct, ml_only = isr._resolve_fallback_signals(hits)
    assert seq_id == pytest.approx(40.0)  # synthetic 90 excluded
    assert struct == pytest.approx(0.61)
    assert emb is None
    assert ml_only is False


def test_resolve_signals_takes_best_ml_signal():
    hits = [
        {"method": "plm_embedding", "embedding_similarity": 0.55},
        {"method": "plm_embedding", "embedding_similarity": 0.81},
    ]
    _, emb, _, ml_only = isr._resolve_fallback_signals(hits)
    assert emb == pytest.approx(0.81)
    assert ml_only is True


# ---------------------------------------------------------------------------
# ML relabel (rows the judge did not label): synthetic identity alone can't reach
# MEDIUM; the ML rescue requires real similarity AND flanking (synteny) corroboration.
# ---------------------------------------------------------------------------

def test_ml_only_locus_without_a_score_stays_low():
    """An ML-only locus carries identity 0 and no coverage; without an ML score above the
    threshold nothing lifts its provisional LOW."""
    assert isr._provisional_goi_label("fallback_hit_span")[0] == "LOW"
    assert isr._ml_relabel("LOW", "fallback_hit_span", flanking_support=10) == ("LOW", None, None)


def test_ml_rescue_requires_flanking():
    # similarity over threshold but no flanking -> stays LOW
    assert isr._ml_relabel("LOW", "fallback_hit_span", 0, embedding_similarity=0.75)[0] == "LOW"
    # similarity over threshold + flanking -> MEDIUM via the rescue
    conf, cls, reason = isr._ml_relabel("LOW", "fallback_hit_span", 2, embedding_similarity=0.75)
    assert (conf, cls, reason) == ("MEDIUM", "probable_goi", "fallback_span_plm_rescued")


def test_structural_rescue_path():
    conf, _, reason = isr._ml_relabel("LOW", "fallback_hit_span", 2, structural_similarity=0.55)
    assert conf == "MEDIUM"
    assert reason == "fallback_span_structural_rescued"


def test_boost_never_lowers_confidence():
    """ML signals only ever raise confidence, by one step; only a gene model goes from
    MEDIUM to HIGH, and only in a block with enough flanking genes."""
    assert isr._ml_relabel("HIGH", "exon_annotation", 3, embedding_similarity=0.95)[0] == "HIGH"
    assert isr._ml_relabel("MEDIUM", "exon_annotation", 3, embedding_similarity=0.1)[0] == "MEDIUM"
    assert isr._ml_relabel("MEDIUM", "exon_annotation", 2, embedding_similarity=0.9) == (
        "HIGH", "confident_goi", "modeled_goi_plm_boosted")
    assert isr._ml_relabel("MEDIUM", "exon_annotation", 1, embedding_similarity=0.9)[0] == "MEDIUM"
    assert isr._ml_relabel("MEDIUM", "fallback_hit_span", 5, embedding_similarity=0.99)[0] == "MEDIUM"
    assert isr._ml_relabel("LOW", "exon_annotation", 5, embedding_similarity=0.99)[0] == "MEDIUM"
    assert isr._ml_relabel("LOW", "raw_hit", 5, embedding_similarity=0.99)[0] == "LOW"


# ---------------------------------------------------------------------------
# plm_search numpy-only helpers
# ---------------------------------------------------------------------------

def test_cosine_similarity_identical_vectors():
    q = {"GOI": np.array([1.0, 0.0, 0.0], dtype=np.float32)}
    t = {"a": np.array([2.0, 0.0, 0.0], dtype=np.float32),   # same direction
         "b": np.array([0.0, 1.0, 0.0], dtype=np.float32)}   # orthogonal
    sims = plm_search.cosine_similarity_matrix(q, t)
    assert sims[("GOI", "a")] == pytest.approx(1.0, abs=1e-5)
    assert sims[("GOI", "b")] == pytest.approx(0.0, abs=1e-5)


def test_best_similarities_picks_max():
    q = {"q1": np.array([1.0, 0.0], dtype=np.float32),
         "q2": np.array([0.0, 1.0], dtype=np.float32)}
    t = {"t1": np.array([1.0, 1.0], dtype=np.float32)}
    best = plm_search.best_similarities(q, t)
    assert best["t1"] == pytest.approx(np.cos(np.pi / 4), abs=1e-5)


def test_chunk_sequence_overlap():
    seq = "A" * 250
    chunks = plm_search._chunk_sequence(seq, chunk_size=100, overlap=20)
    assert all(len(c) <= 100 for c in chunks)
    assert chunks[0] == seq[:100]


def test_embeddings_roundtrip(tmp_path):
    embs = {"GOI_x": np.arange(8, dtype=np.float32),
            "GOI_y": np.ones(8, dtype=np.float32)}
    path = str(tmp_path / "e.npz")
    plm_search.save_embeddings(embs, path)
    loaded = plm_search.load_embeddings(path)
    assert set(loaded) == set(embs)
    np.testing.assert_allclose(loaded["GOI_x"], embs["GOI_x"])


def test_check_plm_available_returns_bool():
    # Must not raise whether or not torch is installed.
    assert isinstance(plm_search.check_plm_available(), bool)


# ---------------------------------------------------------------------------
# structural_search numpy-only helpers + availability
# ---------------------------------------------------------------------------

def test_structure_index_roundtrip(tmp_path):
    idx = {"GOI_a": "/x/a.pdb", "GOI_b": "/x/b.pdb"}
    path = str(tmp_path / "idx.tsv")
    structural_search.save_structure_index(idx, path)
    assert structural_search.load_structure_index(path) == idx


def test_load_missing_structure_index_is_empty(tmp_path):
    assert structural_search.load_structure_index(str(tmp_path / "nope.tsv")) == {}


def test_structural_availability_returns_bool():
    assert isinstance(structural_search.check_esmfold_available(), bool)
    assert isinstance(structural_search.check_foldseek_available(), bool)
    assert isinstance(structural_search.check_structural_search_available(), bool)


def test_safe_filename_sanitises():
    assert structural_search._safe_filename("GOI|x/y z") == "GOI_x_y_z"


@pytest.mark.parametrize("vram,expect_cap", [(4, 150), (8, 300), (16, 400), (40, 10_000)])
def test_vram_tier_caps(vram, expect_cap):
    cap, chunk = structural_search._vram_tier_caps(vram)
    assert cap == expect_cap
    assert chunk in (32, 48, 64)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
