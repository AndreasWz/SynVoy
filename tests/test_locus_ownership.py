#!/usr/bin/env python3
"""Tests for §1m locus ownership (docs/TODO.md §1m).

Two layers:
  1. build_home_paralog_panel.py — builds >locus_<id>|<gene> panel + meta from the
     home GFF, home proteome, and per-locus split BEDs; flags the GOI's own gene.
  2. generate_report.build_locus_ownership — re-attributes each deduped GOI ortholog
     to the home locus its protein aligns best to (RBH), with a flanking-synteny
     tiebreak, emitting locus_reattributed / goi_owner_search_fumble /
     paralog_misassignment flags.

Grounded in the real decorin (P07585) failure: the chr9 osteomodulin locus modelled
real decorin on cow BTA5 (recovered_from = the OMD locus) while the chr12 DCN locus
produced nothing there; a 49.5 % asporin was carried under the GOI_DCN label.
"""

import os
import subprocess
import sys
import tempfile
import unittest

BIN_DIR = os.path.join(os.path.dirname(__file__), "..", "bin")
sys.path.insert(0, BIN_DIR)

from build_home_paralog_panel import (  # noqa: E402
    locus_id_from_path,
    genes_for_locus,
    locus_span,
    find_family_paralogs,
)
from generate_report import build_locus_ownership  # noqa: E402

# A realistic base protein + a ~half-mutated paralog + an unrelated sequence, for the
# homology-family augmentation tests (the decorin/biglycan fix).
_BASE = "MKATIILLLLAQVSWAGPTSDQALVNLAEKLYESGDFEKAVPLLKEAVRLNPNDAEAWNLLG"
_PARALOG = "".join(c if i % 3 else "G" for i, c in enumerate(_BASE))   # ~2/3 identity
_UNRELATED = "WYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWYWY"

PANEL_SCRIPT = os.path.join(BIN_DIR, "build_home_paralog_panel.py")
PY = sys.executable


def _ownership_row(genome, chrom, start, end, best_label, best_bit,
                   second_label="", second_bit=0.0, locus="locus_x", mrna="m1"):
    """Mimic a post-load_paralog_check_rows row (adds the *_f/_i numeric keys)."""
    return {
        "genome": genome, "locus": locus, "chrom": chrom,
        "start": str(start), "end": str(end), "mrna_id": mrna,
        "best_paralog": best_label, "best_bit": str(best_bit),
        "second_paralog": second_label, "second_bit": str(second_bit),
        "bitscore_gap": str(best_bit - second_bit),
        "best_bit_f": float(best_bit), "second_bit_f": float(second_bit),
        "gap_f": float(best_bit - second_bit),
        "start_i": int(start), "end_i": int(end),
    }


def _dedup_record(genome, chrom, start, end, provenance, goi_class="confident_goi",
                  confidence="HIGH", identity="90.0"):
    return {
        "genome": genome, "chrom": chrom, "start": start, "end": end,
        "confidence": confidence, "identity": identity, "goi_class": goi_class,
        "target_gene": "", "model_status": "complete",
        "provenance": list(provenance), "n_source_loci": len(provenance),
    }


class TestPanelBuilderUnits(unittest.TestCase):
    def test_locus_id_from_path(self):
        self.assertEqual(locus_id_from_path("/x/y/locus_4.bed"), "locus_4")
        self.assertEqual(locus_id_from_path("intermediate/split_loci/locus_12.bed"), "locus_12")

    def test_genes_for_locus_overlap_ranks_by_overlap(self):
        idx = {"c1": [(100, 200, "gene-A"), (150, 400, "gene-B"), (9000, 9100, "gene-FAR")]}
        # Locus span 140-210 overlaps A and B; B has the larger overlap.
        got = genes_for_locus("c1", 140, 210, idx, pad=0, fallback_window=1000, max_genes=5)
        self.assertEqual(set(got), {"gene-A", "gene-B"})
        self.assertNotIn("gene-FAR", got)

    def test_genes_for_locus_fallback_nearest(self):
        idx = {"c1": [(1000, 2000, "gene-N")]}
        # No overlap; nearest gene within window is returned.
        got = genes_for_locus("c1", 2500, 2600, idx, pad=0, fallback_window=1000, max_genes=5)
        self.assertEqual(got, ["gene-N"])
        # Outside the fallback window -> nothing.
        got2 = genes_for_locus("c1", 50000, 50100, idx, pad=0, fallback_window=1000, max_genes=5)
        self.assertEqual(got2, [])


class TestPanelBuilderEndToEnd(unittest.TestCase):
    def test_panel_and_meta(self):
        with tempfile.TemporaryDirectory() as d:
            gff = os.path.join(d, "home.gff")
            with open(gff, "w") as fh:
                fh.write("\n".join([
                    "NC_12\tsrc\tgene\t1000\t2000\t.\t-\t.\tID=gene-DCN;Name=DCN",
                    "NC_9\tsrc\tgene\t5000\t6000\t.\t+\t.\tID=gene-OGN;Name=OGN",
                    "NC_9\tsrc\tgene\t6200\t7000\t.\t+\t.\tID=gene-OMD;Name=OMD",
                    "NC_9\tsrc\tgene\t7200\t8000\t.\t+\t.\tID=gene-ASPN;Name=ASPN",
                    "NC_9\tsrc\tgene\t90000\t91000\t.\t+\t.\tID=gene-FAR;Name=FAR",
                ]) + "\n")
            prot = os.path.join(d, "home_proteome.faa")
            with open(prot, "w") as fh:
                fh.write(
                    ">gene-DCN\nMKATIILLLLAQVSWADCNDCNDCNDCNDCNDCNPPPWWWQQQEEERRRTTT\n"
                    ">gene-OGN\nMKKKOGNOGNOGNOGNOGNGGGHHHJJJKKKLLLZZZXXXCCCVVV\n"
                    ">gene-OMD\nMOMDOMDOMDOMDOMDAAASSSDDDFFFGGGHHHJJJKKK\n"
                    ">gene-ASPN\nMASPNASPNASPNASPNQQQWWWEEERRRTTTYYYUUU\n"
                    ">gene-FAR\nMFARFARFARFARFAR\n"
                )
            l1 = os.path.join(d, "locus_1.bed")
            l4 = os.path.join(d, "locus_4.bed")
            with open(l1, "w") as fh:
                fh.write("NC_12\t1100\t1900\tgene_loc\t1e-30\t-\n")
            with open(l4, "w") as fh:
                fh.write("NC_9\t5500\t7500\tgene_loc\t1e-20\t+\n")
            query = os.path.join(d, "query.faa")
            with open(query, "w") as fh:
                fh.write(">P07585 decorin\nMKATIILLLLAQVSWADCNDCNDCNDCNDCNDCNPPPWWWQQQEEERRRTTT\n")
            out = os.path.join(d, "panel.faa")

            r = subprocess.run(
                [PY, PANEL_SCRIPT, "--home_gff", gff, "--home_proteome", prot,
                 "--locus_beds", l1, l4, "--query", query, "--output", out],
                capture_output=True, text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr)

            labels = {ln[1:].split()[0] for ln in open(out) if ln.startswith(">")}
            self.assertIn("locus_1|gene-DCN", labels)
            self.assertEqual({"locus_4|gene-OGN", "locus_4|gene-OMD", "locus_4|gene-ASPN"},
                             labels - {"locus_1|gene-DCN"})
            self.assertNotIn("gene-FAR", "".join(labels))  # outside any locus span

            meta = {}
            with open(out + ".meta.tsv") as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    row = dict(zip(header, line.rstrip("\n").split("\t")))
                    meta[row["panel_id"]] = row
            self.assertEqual(meta["locus_1|gene-DCN"]["is_goi_gene"], "1")
            self.assertEqual(meta["locus_4|gene-ASPN"]["is_goi_gene"], "0")
            self.assertEqual(meta["locus_4|gene-OMD"]["locus_id"], "locus_4")


class TestFamilyAugmentation(unittest.TestCase):
    """§1m fix: the panel must source the GOI's paralog FAMILY from homology, not from
    split_loci spans (which capture LRR cross-hits, missing the real BGN/ASPN/OMD)."""

    def test_find_family_paralogs_threshold_and_exclusion(self):
        proteome = {
            "gene-DCN": _BASE,        # the GOI's own gene (would be 'already' present)
            "gene-BGN": _PARALOG,     # a real paralog ~2/3 identity
            "gene-XYZ": _UNRELATED,   # below threshold
        }
        fam = find_family_paralogs(_BASE, proteome, already={"gene-DCN"},
                                   min_identity=30.0, max_n=15)
        gids = [g for g, _s, _i in fam]
        self.assertIn("gene-BGN", gids)        # paralog included
        self.assertNotIn("gene-DCN", gids)     # already-present excluded
        self.assertNotIn("gene-XYZ", gids)     # below threshold excluded

    def test_find_family_paralogs_caps_and_ranks(self):
        proteome = {f"g{i}": _PARALOG for i in range(10)}
        proteome["gene-IDENT"] = _BASE  # 100% — should rank first
        fam = find_family_paralogs(_BASE, proteome, already=set(),
                                   min_identity=30.0, max_n=3)
        self.assertEqual(len(fam), 3)               # capped
        self.assertEqual(fam[0][0], "gene-IDENT")   # highest identity first

    def test_panel_includes_homology_family_member(self):
        """A paralog NOT overlapping any locus BED is still added (home_<gene>|<gene>)."""
        with tempfile.TemporaryDirectory() as d:
            gff = os.path.join(d, "home.gff")
            with open(gff, "w") as fh:
                fh.write("\n".join([
                    "NC_12\tsrc\tgene\t1000\t2000\t.\t-\t.\tID=gene-DCN;Name=DCN",
                    "NC_X\tsrc\tgene\t5000\t6000\t.\t+\t.\tID=gene-BGN;Name=BGN",
                ]) + "\n")
            prot = os.path.join(d, "home_proteome.faa")
            with open(prot, "w") as fh:
                fh.write(f">gene-DCN\n{_BASE}\n>gene-BGN\n{_PARALOG}\n")
            l1 = os.path.join(d, "locus_1.bed")
            with open(l1, "w") as fh:
                fh.write("NC_12\t1100\t1900\tgene_loc\t1e-30\t-\n")
            query = os.path.join(d, "query.faa")
            with open(query, "w") as fh:
                fh.write(f">P07585\n{_BASE}\n")
            out = os.path.join(d, "panel.faa")
            r = subprocess.run(
                [PY, PANEL_SCRIPT, "--home_gff", gff, "--home_proteome", prot,
                 "--locus_beds", l1, "--query", query, "--output", out,
                 "--paralog_min_identity", "30"],
                capture_output=True, text=True,
            )
            self.assertEqual(r.returncode, 0, r.stderr)
            labels = {ln[1:].split()[0] for ln in open(out) if ln.startswith(">")}
            self.assertIn("locus_1|gene-DCN", labels)            # GOI's own gene (locus)
            self.assertIn("home_gene-BGN|gene-BGN", labels)      # paralog via homology
            meta = {}
            with open(out + ".meta.tsv") as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    row = dict(zip(header, line.rstrip("\n").split("\t")))
                    meta[row["panel_id"]] = row
            self.assertEqual(meta["locus_1|gene-DCN"]["is_goi_gene"], "1")
            self.assertEqual(meta["home_gene-BGN|gene-BGN"]["is_goi_gene"], "0")
            self.assertEqual(meta["home_gene-BGN|gene-BGN"]["chrom"], "NC_X")


class TestBuildLocusOwnership(unittest.TestCase):
    def _panel_meta(self):
        return {
            "locus_1|gene-DCN":  {"locus": "locus_1", "gene": "gene-DCN", "is_goi": True},
            "locus_4|gene-OMD":  {"locus": "locus_4", "gene": "gene-OMD", "is_goi": False},
            "locus_4|gene-ASPN": {"locus": "locus_4", "gene": "gene-ASPN", "is_goi": False},
        }

    def test_reattribution_and_owner_fumble(self):
        # Real decorin recovered ONLY by the chr9 OMD locus (locus_4); aligns best to
        # the GOI's own gene at locus_1 -> re-attributed + owner-search-fumble.
        rec = _dedup_record("cow", "NC_t5", 100, 500, provenance=["locus_4"])
        goi_dedup = {"records": [rec]}
        rows = [_ownership_row("cow", "NC_t5", 100, 500,
                               best_label="locus_1|gene-DCN", best_bit=375.0,
                               second_label="locus_4|gene-ASPN", second_bit=104.0,
                               locus="locus_4")]
        flags, summary = build_locus_ownership(goi_dedup, rows, self._panel_meta(), {})
        self.assertEqual(rec["owning_locus"], "locus_1")
        self.assertTrue(rec["locus_reattributed"])
        types = {f["type"] for f in flags}
        self.assertIn("locus_reattributed", types)
        self.assertIn("goi_owner_search_fumble", types)
        self.assertEqual(summary["n_reattributed"], 1)
        self.assertEqual(summary["n_owner_search_fumble"], 1)
        self.assertEqual(summary["n_paralog_misassignment"], 0)

    def test_paralog_misassignment_relabel(self):
        # 49.5 % asporin carried under GOI_DCN; aligns best to gene-ASPN (non-GOI).
        rec = _dedup_record("cow", "NC_t8", 900, 1300, provenance=["locus_4"],
                            confidence="MEDIUM", goi_class="probable_goi", identity="49.5")
        goi_dedup = {"records": [rec]}
        rows = [_ownership_row("cow", "NC_t8", 900, 1300,
                               best_label="locus_4|gene-ASPN", best_bit=248.0,
                               second_label="locus_1|gene-DCN", second_bit=92.0,
                               locus="locus_4")]
        flags, summary = build_locus_ownership(goi_dedup, rows, self._panel_meta(), {})
        self.assertEqual(rec["owning_locus"], "locus_4")
        self.assertEqual(rec["goi_class"], "paralog_not_goi")
        self.assertEqual(rec["inferred_paralog"], "gene-ASPN")
        self.assertFalse(rec["locus_reattributed"])  # owner == recovered-from
        self.assertIn("paralog_misassignment", {f["type"] for f in flags})
        self.assertEqual(summary["n_paralog_misassignment"], 1)
        self.assertEqual(summary["n_owner_search_fumble"], 0)

    def test_homology_paralog_owner_flags_not_reattributes(self):
        # The decorin/biglycan case: a HIGH 'GOI_DCN' call whose RBH-best is a homology-family
        # paralog (home_gene-BGN, is_goi=False) must be relabeled paralog_not_goi and flagged
        # paralog_misassignment — NOT counted as a locus reattribution to a searched locus.
        panel_meta = {
            "locus_1|gene-DCN": {"locus": "locus_1", "gene": "gene-DCN", "is_goi": True},
            "home_gene-BGN|gene-BGN": {"locus": "home_gene-BGN", "gene": "gene-BGN",
                                       "is_goi": False},
        }
        rec = _dedup_record("cow", "NC_X", 100, 500, provenance=["locus_1"])  # HIGH
        goi_dedup = {"records": [rec]}
        rows = [_ownership_row("cow", "NC_X", 100, 500,
                               best_label="home_gene-BGN|gene-BGN", best_bit=300.0,
                               second_label="locus_1|gene-DCN", second_bit=200.0,
                               locus="locus_1")]
        flags, summary = build_locus_ownership(goi_dedup, rows, panel_meta, {})
        self.assertEqual(rec["goi_class"], "paralog_not_goi")
        self.assertEqual(rec["inferred_paralog"], "gene-BGN")
        self.assertFalse(rec["locus_reattributed"])           # paralog owner != reattribution
        types = {f["type"] for f in flags}
        self.assertIn("paralog_misassignment", types)
        self.assertNotIn("locus_reattributed", types)
        self.assertEqual(summary["n_paralog_misassignment"], 1)
        self.assertEqual(summary["n_reattributed"], 0)

    def test_synteny_tiebreak_flips_near_tie(self):
        # Near-tied panel scores between two loci, neither the GOI; the flanking neighbourhood on
        # the target chrom favours locus_4 -> ownership flips to locus_4.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        goi_dedup = {"records": [rec]}
        rows = [_ownership_row("frog", "NC_q1", 200, 600,
                               best_label="locus_2|gene-TWO", best_bit=200.0,
                               second_label="locus_4|gene-OMD", second_bit=195.0,  # gap 5 < 10
                               locus="locus_1")]
        flanking = {("frog", "locus_2", "NC_q1"): 1, ("frog", "locus_4", "NC_q1"): 6}
        meta = dict(self._panel_meta(), **{"locus_2|gene-TWO": {"locus": "locus_2", "gene": "gene-TWO", "is_goi": False}})
        flags, summary = build_locus_ownership(goi_dedup, rows, meta, flanking, tiebreak_gap=10.0)
        self.assertEqual(rec["owning_locus"], "locus_4")
        self.assertTrue(rec["ownership_tiebroken"])
        self.assertEqual(summary["n_tiebreak_applied"], 1)

    def test_synteny_does_not_take_a_near_tie_from_the_goi(self):
        # The same near tie with the GOI on top: only a non-GOI gene that beats the GOI by the
        # gap owns a call, so the flanking of locus_4 changes nothing (fix-phase item 10).
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        rows = [_ownership_row("frog", "NC_q1", 200, 600,
                               best_label="locus_1|gene-DCN", best_bit=200.0,
                               second_label="locus_4|gene-OMD", second_bit=195.0,
                               locus="locus_1")]
        flanking = {("frog", "locus_1", "NC_q1"): 1, ("frog", "locus_4", "NC_q1"): 6}
        flags, summary = build_locus_ownership({"records": [rec]}, rows, self._panel_meta(), flanking, tiebreak_gap=10.0)
        self.assertEqual((rec["owning_locus"], rec["owning_gene"], rec["goi_class"]), ("locus_1", "gene-DCN", "confident_goi"))
        self.assertNotIn("ownership_near_tie", rec)       # the GOI is on top: nothing to tie with
        self.assertFalse(rec["ownership_tiebroken"])
        self.assertEqual((summary["n_near_tie"], summary["n_tiebreak_applied"]), (0, 0))
        self.assertEqual(flags, [])

    def test_a_paralog_must_beat_the_goi_by_the_gap(self):
        def run(paralog_bit):
            rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
            rows = [_ownership_row("frog", "NC_q1", 200, 600,
                                   best_label="locus_4|gene-OMD", best_bit=paralog_bit,
                                   second_label="locus_1|gene-DCN", second_bit=200.0,
                                   locus="locus_1")]
            flags, _s = build_locus_ownership({"records": [rec]}, rows, self._panel_meta(), {}, tiebreak_gap=10.0)
            return rec, flags
        rec, flags = run(209.9)
        self.assertEqual((rec["goi_class"], rec["owning_is_goi_gene"], rec["owning_gene"]), ("confident_goi", True, "gene-DCN"))
        self.assertTrue(rec["ownership_near_tie"])
        self.assertEqual(flags, [])
        rec, flags = run(210.0)
        self.assertEqual((rec["goi_class"], rec["owning_is_goi_gene"], rec["inferred_paralog"]), ("paralog_not_goi", False, "gene-OMD"))
        self.assertEqual([f["type"] for f in flags], ["paralog_misassignment"])

    def test_the_goi_entry_is_read_from_all_panel_scores(self):
        # Two paralogs on top and the GOI third: the table's two best do not name it.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        row = _ownership_row("frog", "NC_q1", 200, 600, best_label="locus_4|gene-OMD", best_bit=205.0,
                             second_label="locus_4|gene-ASPN", second_bit=203.0, locus="locus_1")
        row["panel_scores_d"] = {"locus_4|gene-OMD": 205.0, "locus_4|gene-ASPN": 203.0, "locus_1|gene-DCN": 200.0}
        build_locus_ownership({"records": [rec]}, [row], self._panel_meta(), {}, tiebreak_gap=10.0)
        self.assertEqual((rec["goi_class"], rec["owning_gene"]), ("confident_goi", "gene-DCN"))
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        row["panel_scores_d"]["locus_1|gene-DCN"] = 190.0
        build_locus_ownership({"records": [rec]}, [row], self._panel_meta(), {}, tiebreak_gap=10.0)
        self.assertEqual((rec["goi_class"], rec["inferred_paralog"]), ("paralog_not_goi", "gene-OMD"))

    def test_a_table_that_names_no_goi_entry_relabels_as_before(self):
        # An older table with the GOI third: its score is not known, so the relabel stands.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        rows = [_ownership_row("frog", "NC_q1", 200, 600, best_label="locus_4|gene-OMD", best_bit=205.0,
                               second_label="locus_4|gene-ASPN", second_bit=203.0, locus="locus_1")]
        build_locus_ownership({"records": [rec]}, rows, self._panel_meta(), {}, tiebreak_gap=10.0)
        self.assertEqual((rec["goi_class"], rec["inferred_paralog"]), ("paralog_not_goi", "gene-OMD"))

    def _copy_meta(self):
        return dict(self._panel_meta(), **{"home_gene-C|gene-C": {"locus": "home_gene-C", "gene": "gene-C", "is_goi": True,
                                                                  "goi_found_by": "own_copy"},
                                           "locus_1|gene-DCN": {"locus": "locus_1", "gene": "gene-DCN", "is_goi": True,
                                                                "goi_found_by": "goi_info"}})

    def _copy_case(self, copy_bit, provenance=("locus_1",)):
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=list(provenance))
        row = _ownership_row("frog", "NC_q1", 200, 600, best_label="home_gene-C|gene-C", best_bit=copy_bit,
                             second_label="locus_1|gene-DCN", second_bit=300.0, locus="locus_1")
        row["panel_scores_d"] = {"home_gene-C|gene-C": copy_bit, "locus_1|gene-DCN": 300.0, "locus_4|gene-OMD": 100.0}
        flags, summary = build_locus_ownership({"records": [rec]}, [row], self._copy_meta(), {}, tiebreak_gap=10.0)
        return rec, flags, summary

    def test_a_home_copy_the_call_does_not_favour_is_the_goi(self):
        # A copy on top, by less than the bar of the search's p*: the call stays the GOI's, at a
        # searched locus, with no flag: no search started from the pseudo-locus home_gene-C.
        for bit in (300.0, 313.0, 325.0):
            rec, flags, _s = self._copy_case(bit)
            self.assertEqual((rec["goi_class"], rec["owning_is_goi_gene"], rec["owning_gene"], rec["owning_locus"]),
                             ("confident_goi", True, "gene-C", "locus_1"))
            self.assertFalse(rec["locus_reattributed"])
            self.assertEqual(flags, [])

    def test_a_home_copy_the_call_favours_is_a_paralog(self):
        # 26 raw = 10 bits: the call is that copy's, as in the search's specificity.
        rec, flags, _s = self._copy_case(326.0)
        self.assertEqual((rec["goi_class"], rec["owning_is_goi_gene"], rec["inferred_paralog"]),
                         ("paralog_not_goi", False, "gene-C"))
        self.assertEqual([f["type"] for f in flags], ["paralog_misassignment"])

    def test_a_near_tie_held_by_an_own_copy_stays_at_a_searched_locus(self):
        # A paralog on top, within the gap of an own copy that is itself above the GOI's own entry: the
        # call stays the GOI's, and not at the copy's pseudo-locus home_gene-C.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_1"])
        row = _ownership_row("frog", "NC_q1", 200, 600, best_label="locus_4|gene-OMD", best_bit=1565.0,
                             second_label="home_gene-C|gene-C", second_bit=1559.0, locus="locus_1")
        row["panel_scores_d"] = {"locus_4|gene-OMD": 1565.0, "home_gene-C|gene-C": 1559.0, "locus_1|gene-DCN": 1548.0}
        flags, summary = build_locus_ownership({"records": [rec]}, [row], self._copy_meta(), {}, tiebreak_gap=10.0)
        self.assertEqual((rec["goi_class"], rec["owning_is_goi_gene"], rec["owning_gene"], rec["owning_locus"]),
                         ("confident_goi", True, "gene-C", "locus_1"))
        self.assertTrue(rec["ownership_near_tie"])
        self.assertEqual((flags, summary["n_near_tie"]), ([], 1))

    def test_no_rows_gives_the_same_summary_keys(self):
        _f, summary = build_locus_ownership({"records": []}, [], {}, {})
        self.assertEqual(summary["n_near_tie"], 0)

    def test_a_near_tie_raises_no_locus_flag(self):
        # Found only from a paralog's locus, the paralog 5 above the GOI: the call stays the GOI's,
        # and that says nothing about the GOI's locus having failed.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_4"])
        rows = [_ownership_row("frog", "NC_q1", 200, 600, best_label="locus_4|gene-OMD", best_bit=200.0,
                               second_label="locus_1|gene-DCN", second_bit=195.0, locus="locus_4")]
        flanking = {("frog", "locus_1", "NC_q1"): 1, ("frog", "locus_4", "NC_q1"): 6}
        flags, summary = build_locus_ownership({"records": [rec]}, rows, self._panel_meta(), flanking, tiebreak_gap=10.0)
        self.assertEqual((rec["goi_class"], rec["owning_locus"], rec["locus_reattributed"]), ("confident_goi", "locus_1", False))
        self.assertTrue(rec["ownership_near_tie"])
        self.assertEqual((flags, summary["n_near_tie"], summary["n_tiebreak_applied"], summary["n_reattributed"]), ([], 1, 0, 0))

    def test_decisive_gap_ignores_tiebreak(self):
        # Large gap -> RBH wins outright even if flanking points elsewhere.
        rec = _dedup_record("frog", "NC_q1", 200, 600, provenance=["locus_4"])
        goi_dedup = {"records": [rec]}
        rows = [_ownership_row("frog", "NC_q1", 200, 600,
                               best_label="locus_1|gene-DCN", best_bit=375.0,
                               second_label="locus_4|gene-OMD", second_bit=104.0,
                               locus="locus_4")]
        flanking = {("frog", "locus_1", "NC_q1"): 0, ("frog", "locus_4", "NC_q1"): 9}
        flags, summary = build_locus_ownership(goi_dedup, rows, self._panel_meta(), flanking,
                                               tiebreak_gap=10.0)
        self.assertEqual(rec["owning_locus"], "locus_1")
        self.assertFalse(rec["ownership_tiebroken"])

    def test_rescue_tagged_genome_still_matches(self):
        """A rescue-pass ownership row must reach the record it belongs to.

        RESCUE_GOI_HULL's ownership task is tagged "<genome>.hull_rescue" so its
        output filename cannot collide with the seeded task's, and that tag is
        written verbatim into the TSV's `genome` column. The dedup records are
        already canonicalized, so build_locus_ownership has to canonicalize the
        row's genome too — otherwise the key never matches and every rescue-derived
        call silently skips the RBH check, which is the guard against exactly the
        mislabelled calls the rescue path is most likely to produce.
        """
        for tag in ("cow", "cow.hull_rescue", "cow.rescue", "cow.fna"):
            with self.subTest(genome_column=tag):
                rec = _dedup_record("cow", "NC_t5", 100, 500, provenance=["locus_1"])
                rows = [_ownership_row(tag, "NC_t5", 100, 500,
                                       best_label="locus_9|home_gene-BGN", best_bit=950.0,
                                       second_label="locus_1|gene-DCN", second_bit=507.0,
                                       locus="locus_1")]
                flags, summary = build_locus_ownership({"records": [rec]}, rows,
                                                       self._panel_meta(), {})
                self.assertEqual(summary["evaluated"], 1)
                self.assertEqual(rec["goi_class"], "paralog_not_goi")

    def test_no_ownership_rows_is_noop(self):
        rec = _dedup_record("cow", "NC_t5", 100, 500, provenance=["locus_4"])
        goi_dedup = {"records": [rec]}
        flags, summary = build_locus_ownership(goi_dedup, [], self._panel_meta(), {})
        self.assertEqual(flags, [])
        self.assertEqual(summary["evaluated"], 0)
        self.assertNotIn("owning_locus", rec)


if __name__ == "__main__":
    unittest.main()
