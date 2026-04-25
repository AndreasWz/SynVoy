#!/usr/bin/env python3
"""Tests for the VRAM tier probe in structural_search.py (plan item 4f).

These cover only the pure-logic helpers. The CUDA probe is exercised at
runtime on real hardware — see LOCAL_PLAN.md §8 for the deferred
hardware test on the GTX 1650.
"""

import os
import sys
import unittest
from unittest import mock


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))


def _import_ss():
    """Import lazily so the test file doesn't try to pull in torch at import."""
    import structural_search as ss  # noqa: E402
    return ss


class TestVramTierCaps(unittest.TestCase):
    def test_gtx_1650_class_under_6gb(self):
        ss = _import_ss()
        cap, chunk = ss._vram_tier_caps(4.0)
        self.assertEqual(cap, 150)
        self.assertEqual(chunk, 32)

    def test_rtx_3060_class_under_10gb(self):
        ss = _import_ss()
        cap, chunk = ss._vram_tier_caps(8.0)
        self.assertEqual(cap, 300)
        self.assertEqual(chunk, 48)

    def test_rtx_3090_class_under_20gb(self):
        ss = _import_ss()
        cap, chunk = ss._vram_tier_caps(16.0)
        self.assertEqual(cap, 400)
        self.assertEqual(chunk, 64)

    def test_a100_class_20gb_and_up(self):
        ss = _import_ss()
        cap, chunk = ss._vram_tier_caps(24.0)
        self.assertGreaterEqual(cap, 1000)  # effectively uncapped
        self.assertEqual(chunk, 64)
        cap, chunk = ss._vram_tier_caps(80.0)
        self.assertGreaterEqual(cap, 1000)
        self.assertEqual(chunk, 64)

    def test_tier_boundaries_are_strict_less_than(self):
        ss = _import_ss()
        self.assertEqual(ss._vram_tier_caps(5.99)[0], 150)
        self.assertEqual(ss._vram_tier_caps(6.0)[0], 300)
        self.assertEqual(ss._vram_tier_caps(9.99)[0], 300)
        self.assertEqual(ss._vram_tier_caps(10.0)[0], 400)
        self.assertEqual(ss._vram_tier_caps(19.99)[0], 400)
        self.assertEqual(ss._vram_tier_caps(20.0)[0], 10_000)


class TestEffectiveMaxLength(unittest.TestCase):
    def test_cpu_device_honors_request(self):
        ss = _import_ss()
        self.assertEqual(ss._effective_max_length("cpu", 700), 700)
        self.assertEqual(ss._effective_max_length("cpu", 50), 50)

    def test_cuda_but_no_vram_probe_honors_request(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=None):
            self.assertEqual(ss._effective_max_length("cuda", 700), 700)

    def test_cuda_4gb_caps_to_150(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=4.0):
            self.assertEqual(ss._effective_max_length("cuda", 700), 150)
            self.assertEqual(ss._effective_max_length("cuda", 100), 100)

    def test_cuda_16gb_caps_to_400(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=16.0):
            self.assertEqual(ss._effective_max_length("cuda", 700), 400)
            self.assertEqual(ss._effective_max_length("cuda", 250), 250)

    def test_cuda_24gb_uncapped(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=24.0):
            self.assertEqual(ss._effective_max_length("cuda", 700), 700)


class TestRecommendedChunkSize(unittest.TestCase):
    def test_cpu_returns_default(self):
        ss = _import_ss()
        self.assertEqual(ss._recommended_chunk_size("cpu", default=64), 64)

    def test_cuda_4gb_returns_32(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=4.0):
            self.assertEqual(ss._recommended_chunk_size("cuda"), 32)

    def test_cuda_8gb_returns_48(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=8.0):
            self.assertEqual(ss._recommended_chunk_size("cuda"), 48)

    def test_cuda_no_probe_keeps_default(self):
        ss = _import_ss()
        with mock.patch.object(ss, "_probe_vram_gb", return_value=None):
            self.assertEqual(ss._recommended_chunk_size("cuda", default=64), 64)


if __name__ == "__main__":
    unittest.main()
