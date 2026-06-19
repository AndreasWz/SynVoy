#!/usr/bin/env python3
"""
A4 rank/quantile wave-binning regression tests (TODO_JUN §A4).

The legacy absolute-distance binning collapses for a tight clade whose phylo
distances saturate (e.g. all 19 bees land at dist≈0.99 after max-normalization),
producing ~4 undifferentiated waves with no close→far gradient. `assign_waves_by_rank`
bins by RANK instead, so the closest genomes go in small (serial) waves — the
precondition for §A3's closest-first seeding — and distant genomes in larger waves.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from iterative_search_runner import assign_waves_by_rank


def entries(dists):
    """Build genome_entries pre-sorted closest-first."""
    return [{'name': f'g{i}.fna', 'path': f'g{i}.fna', 'dist': d}
            for i, d in enumerate(dists)]


class TestRankWaveBinning(unittest.TestCase):

    def test_saturated_distances_still_grade(self):
        # The melittin pathology: 19 bees all at dist≈0.99 (0.7% spread).
        sat = entries([0.993 + i * 0.0004 for i in range(19)])
        waves = assign_waves_by_rank(sat)
        # Must produce a graded set of waves, not ~4 lumps.
        self.assertGreaterEqual(len(waves), 5)
        # Every genome is placed exactly once, order preserved.
        flat = [e for w in waves for e in w]
        self.assertEqual(len(flat), 19)
        self.assertEqual([e['name'] for e in flat], [e['name'] for e in sat])

    def test_closest_genomes_are_serial(self):
        # The closest tier (rank-quantile < 0.10) must be waves of 1 so a close
        # relative's GOI seeds the DB before divergent genomes are searched.
        waves = assign_waves_by_rank(entries([0.99] * 19))
        self.assertEqual(len(waves[0]), 1)
        self.assertEqual(len(waves[1]), 1)

    def test_distant_genomes_parallelize(self):
        # The farthest tier (rank-quantile >= 0.70) groups into larger waves.
        waves = assign_waves_by_rank(entries([0.99] * 19))
        self.assertEqual(max(len(w) for w in waves), 5)
        # The largest wave is among the last (farthest) waves.
        self.assertEqual(len(waves[-1]), 5)

    def test_wave_size_grows_monotonically_nondecreasing(self):
        waves = assign_waves_by_rank(entries([0.99] * 30))
        sizes = [len(w) for w in waves]
        # All but the final wave grow non-decreasing close→far; the last wave is
        # whatever genomes remain (may be a partial chunk, hence smaller).
        head = sizes[:-1]
        self.assertEqual(head, sorted(head),
                         f"wave sizes should be non-decreasing close→far (excl. remainder): {sizes}")
        self.assertLessEqual(sizes[-1], max(sizes))

    def test_deterministic(self):
        e = entries([0.99] * 19)
        a = assign_waves_by_rank(e)
        b = assign_waves_by_rank(e)
        self.assertEqual([[x['name'] for x in w] for w in a],
                         [[x['name'] for x in w] for w in b])

    def test_single_genome(self):
        waves = assign_waves_by_rank(entries([0.5]))
        self.assertEqual(len(waves), 1)
        self.assertEqual(len(waves[0]), 1)

    def test_empty(self):
        self.assertEqual(assign_waves_by_rank([]), [])

    def test_partition_is_complete_for_various_n(self):
        for n in (2, 3, 7, 11, 20, 50):
            e = entries([0.99] * n)
            waves = assign_waves_by_rank(e)
            flat = [x for w in waves for x in w]
            self.assertEqual(len(flat), n, f"n={n} lost/duplicated genomes")
            self.assertEqual([x['name'] for x in flat], [x['name'] for x in e])


if __name__ == '__main__':
    unittest.main()
