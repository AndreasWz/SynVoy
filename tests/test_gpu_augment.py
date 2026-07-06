#!/usr/bin/env python3
"""Tests for bin/gpu_augment.py — the main-process GPU augmentation engine.

These cover the pure orchestration/merge logic with the GPU model calls mocked,
so they run on any machine (no torch / CUDA / foldseek). The real fold/embed is
exercised on GPU hardware (LRZ). The point here is that the engine:
  * loads manifests / templates correctly,
  * is FAIL-LOUD when a requested layer is unavailable (the whole reason this
    rework exists — the old code silently produced no scores),
  * assembles re-rank scores and discovery hits from embed/fold results,
  * raises the silent-no-op guard when every fold fails.
"""

import json
import os
import sys
import unittest
from unittest import mock

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "bin"))

import gpu_augment as g  # noqa: E402


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


class TestManifestAndTemplates(unittest.TestCase):
    def test_load_manifest_skips_bad_and_empty(self):
        import tempfile
        d = tempfile.mkdtemp()
        mf = os.path.join(d, "m.jsonl")
        _write(mf, "\n".join([
            json.dumps({"kind": "candidate", "id": "c1", "seq": "MKAYIL", "genome": "g1"}),
            "{ this is not json",
            json.dumps({"kind": "orf", "id": "o1", "seq": "MEEP", "genome": "g1",
                        "chrom": "c", "start": 5, "end": 20}),
            json.dumps({"kind": "candidate", "id": "c2", "seq": ""}),   # empty seq -> skip
            json.dumps({"kind": "candidate", "seq": "MMMM"}),            # no id -> skip
        ]) + "\n")
        items = g.load_manifest(mf)
        self.assertEqual([i.wid for i in items], ["c1", "o1"])
        orf = [i for i in items if i.kind == "orf"][0]
        self.assertEqual(orf.meta["start"], 5)
        self.assertEqual(orf.meta["chrom"], "c")

    def test_load_manifest_missing_file(self):
        self.assertEqual(g.load_manifest("/no/such/file"), [])

    def test_load_goi_templates_prefix_filter(self):
        import tempfile
        d = tempfile.mkdtemp()
        db = os.path.join(d, "db.faa")
        _write(db, ">GOI_osk\nMKAYILGG\n>flank_1\nMMMM\n>GOI_osk_exon_1\nMKAY\n")
        tpls = g.load_goi_templates_from_fasta(db)
        self.assertEqual(sorted(t[0] for t in tpls), ["GOI_osk", "GOI_osk_exon_1"])


class TestFailLoudContract(unittest.TestCase):
    """The core reason for the rework: never a silent no-op."""

    def _db(self):
        import tempfile
        d = tempfile.mkdtemp()
        db = os.path.join(d, "db.faa")
        _write(db, ">GOI_osk\nMKAYILGGVQ\n")
        return d, db

    def test_require_raises_when_structural_unavailable(self):
        d, db = self._db()
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1")]
        with mock.patch.object(g, "structural_available", return_value=False):
            with self.assertRaises(RuntimeError):
                g.run_augmentation(db, items, os.path.join(d, "w"),
                                   do_structural=True, require=True)

    def test_require_raises_when_plm_unavailable(self):
        d, db = self._db()
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1")]
        with mock.patch.object(g, "plm_available", return_value=False):
            with self.assertRaises(RuntimeError):
                g.run_augmentation(db, items, os.path.join(d, "w"),
                                   do_plm=True, require=True)

    def test_no_require_skips_gracefully(self):
        d, db = self._db()
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1")]
        with mock.patch.object(g, "structural_available", return_value=False):
            res = g.run_augmentation(db, items, os.path.join(d, "w"),
                                     do_structural=True, require=False)
        self.assertEqual(res["scores"], {})
        self.assertEqual(res["diagnostics"]["structural"]["skipped"], "unavailable")

    def test_no_templates_raises(self):
        import tempfile
        d = tempfile.mkdtemp()
        db = os.path.join(d, "db.faa")
        _write(db, ">flank_only\nMMMM\n")   # no GOI_ prefix
        with self.assertRaises(RuntimeError):
            g.run_augmentation(db, [], os.path.join(d, "w"), do_plm=False, do_structural=False)


class TestScoreAssembly(unittest.TestCase):
    """With embed/fold mocked, verify scores + discovery hits are assembled."""

    def _run(self, emb_best, tm_best, items, **kw):
        import tempfile
        d = tempfile.mkdtemp()
        db = os.path.join(d, "db.faa")
        _write(db, ">GOI_osk\nMKAYILGGVQ\n")
        with mock.patch.object(g, "plm_available", return_value=True), \
             mock.patch.object(g, "structural_available", return_value=True), \
             mock.patch.object(g, "_embed_all", return_value=(emb_best, {"embedded": len(emb_best)})), \
             mock.patch.object(g, "_fold_all", return_value=({"GOI_osk": "x.pdb"}, {}, {"folded_items": 1})), \
             mock.patch.object(g, "_foldseek_best_tm", return_value=tm_best):
            return g.run_augmentation(db, items, os.path.join(d, "w"),
                                      do_plm=True, do_structural=True, require=True, **kw)

    def test_rerank_scores_merge_emb_and_tm(self):
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1"),
                 g.WorkItem("candidate", "c2", "MKAYILGG", "g1")]
        res = self._run({"c1": 0.9}, {"c1": 0.8, "c2": 0.4}, items)
        self.assertAlmostEqual(res["scores"]["c1"]["embedding_similarity"], 0.9)
        self.assertAlmostEqual(res["scores"]["c1"]["structural_similarity"], 0.8)
        # c2 had only a TM score, no embedding
        self.assertIsNone(res["scores"]["c2"]["embedding_similarity"])
        self.assertAlmostEqual(res["scores"]["c2"]["structural_similarity"], 0.4)

    def test_discovery_hits_above_threshold(self):
        items = [g.WorkItem("orf", "o1", "MKAYILGG", "g1", {"chrom": "c", "start": 1, "end": 30}),
                 g.WorkItem("orf", "o2", "MKAYILGG", "g1")]
        res = self._run({"o1": 0.8}, {"o1": 0.6, "o2": 0.1},
                        items, plm_threshold=0.5, structural_tm_threshold=0.3)
        hit_ids = {h["id"] for h in res["discovery_hits"]}
        self.assertIn("o1", hit_ids)          # emb 0.8 >= 0.5 and tm 0.6 >= 0.3
        self.assertNotIn("o2", hit_ids)       # emb None, tm 0.1 < 0.3
        o1 = [h for h in res["discovery_hits"] if h["id"] == "o1"][0]
        self.assertEqual(o1["method"], "foldseek_structural")
        self.assertEqual(o1["chrom"], "c")

    def test_candidates_never_emitted_as_discovery(self):
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1")]
        res = self._run({"c1": 0.99}, {"c1": 0.99}, items)
        self.assertEqual(res["discovery_hits"], [])


class TestSilentNoOpGuard(unittest.TestCase):
    def test_fold_all_raises_when_all_items_fail(self):
        """_fold_all must raise if a non-empty item set yields zero structures
        — the exact silent no-op the fork-CUDA bug caused."""
        items = [g.WorkItem("candidate", "c1", "MKAYILGG", "g1")]
        fake_ss = mock.MagicMock()
        # GOI template folds fine (1 window), but every item fold returns nothing.
        fake_ss.fold_protein_windows.side_effect = (
            lambda sid, seq, d, **k: [(sid, "goi.pdb")] if sid.startswith("GOI_") else []
        )
        with mock.patch.object(g, "_struct_mod", return_value=fake_ss):
            with self.assertRaises(RuntimeError):
                g._fold_all([("GOI_osk", "MKAYILGGVQ")], items, "/tmp/x", "cuda", 380, 60)
        fake_ss.unload_esmfold.assert_called_once()


if __name__ == "__main__":
    unittest.main()
