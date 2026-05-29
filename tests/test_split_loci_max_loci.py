#!/usr/bin/env python3
"""Regression tests for bin/split_loci.py locus cap + bit-score filter (docs/TODO.md §1d).

Background. Xaver ran SynVoy with Photinus pyralis luciferase (XP_031329057.1), an
AMP-binding / acyl-CoA synthetase family member with dozens of paralogs. split_loci.py's
old filter OR'd "has >=2 HSPs" with an e-value ratio test whose `primary_best == 0`
short-circuit degenerated to "keep anything <= 1e-10" once the primary hit's BLAST e-value
saturated at 0.0. All 49 paralog loci leaked through and the run took ~55 h.

The fix ranks secondary loci by best_bit / primary_bit (bit scores do not saturate) and
caps the result at --max_loci. The two real luciferase loci (Luc1 on VVIM01000001.1 and
Luc2 on VVIM01000009.1) rank 1st and 2nd by bit score, so they must always survive.

These tests drive the real LOCATE_GENE chain: the awk transform from locate_gene.nf, then
merge_hits.py, then split_loci.py, using the real hits_blast.txt captured from that run.
"""
import os
import subprocess
import sys
import unittest

BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
MERGE_HITS = os.path.join(BIN_DIR, "merge_hits.py")
SPLIT_LOCI = os.path.join(BIN_DIR, "split_loci.py")
LUCIFERASE_HITS = os.path.join(FIXTURES, "luciferase_hits_blast.txt")

# Coordinates of the two genuine luciferase loci (GFF-confirmed; see TODO §1d).
LUC1_CHROM, LUC1_COORD = "VVIM01000001.1", 28341575   # PPYR_00001, product="luciferase"
LUC2_CHROM, LUC2_COORD = "VVIM01000009.1", 30175000   # PPYR_00002, "diffuse expressed luciferase 2"
# A representative weak AMP-binding paralog that must be dropped (bit 183 = 0.44x primary).
PARALOG_CHROM, PARALOG_COORD = "VVIM01000005.1", 33530000


def _awk_to_bed(hits_blast_path, out_path, with_bit=True):
    """Replicate the locate_gene.nf awk that turns BLAST outfmt 6 into bed-like rows.

    outfmt 6 cols (1-based): 2=sseqid 9=sstart 10=send 11=evalue 12=bitscore.
    Emits: chrom, start(0-based), end, "blast", evalue, strand[, bitscore].
    """
    with open(hits_blast_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            chrom, s, e = p[1], int(p[8]), int(p[9])
            strand = "+"
            if s > e:
                s, e, strand = e, s, "-"
            cols = [chrom, str(s - 1), str(e), "blast", p[10], strand]
            if with_bit:
                cols.append(p[11])
            fout.write("\t".join(cols) + "\n")


def _run_split(tmpdir, hits_blast_path, with_bit=True, extra_args=None):
    """Run awk -> merge_hits.py -> split_loci.py; return (loci, stdout, stderr)."""
    blast_bed = os.path.join(tmpdir, "blast.bed")
    mmseqs_bed = os.path.join(tmpdir, "mmseqs.bed")
    home_bed = os.path.join(tmpdir, "home_gene_location.bed")
    open(mmseqs_bed, "w").close()
    _awk_to_bed(hits_blast_path, blast_bed, with_bit=with_bit)

    subprocess.run(
        [sys.executable, MERGE_HITS, "--mmseqs", mmseqs_bed, "--blast", blast_bed,
         "--max_evalue", "0.01", "--output", home_bed],
        check=True, capture_output=True, text=True,
    )
    proc = subprocess.run(
        [sys.executable, SPLIT_LOCI, "--bed", home_bed, "--output_prefix",
         os.path.join(tmpdir, "locus")] + (extra_args or []),
        check=True, capture_output=True, text=True,
    )
    loci = []
    for fn in sorted(os.listdir(tmpdir)):
        if not (fn.startswith("locus_") and fn.endswith(".bed")):
            continue
        chrom = None
        starts, ends = [], []
        with open(os.path.join(tmpdir, fn)) as f:
            for line in f:
                parts = line.split("\t")
                chrom = parts[0]
                starts.append(int(parts[1]))
                ends.append(int(parts[2]))
        loci.append({"file": fn, "chrom": chrom, "min": min(starts), "max": max(ends)})
    return loci, proc.stdout, proc.stderr


def _covers(loci, chrom, coord):
    return any(l["chrom"] == chrom and l["min"] <= coord <= l["max"] for l in loci)


class TestSplitLociLuciferase(unittest.TestCase):
    def test_family_capped_real_run(self):
        """Defaults: the 49+ luciferase paralog loci collapse to <=5, keeping Luc1 + Luc2."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            loci, out, err = _run_split(d, LUCIFERASE_HITS)
            self.assertLessEqual(len(loci), 5, "max_loci cap (5) not enforced")
            self.assertTrue(_covers(loci, LUC1_CHROM, LUC1_COORD),
                            "Luc1 locus (VVIM01000001.1) was dropped")
            self.assertTrue(_covers(loci, LUC2_CHROM, LUC2_COORD),
                            "Luc2 locus (VVIM01000009.1) was dropped")
            self.assertFalse(_covers(loci, PARALOG_CHROM, PARALOG_COORD),
                             "weak AMP-binding paralog (ratio 0.44) leaked through")
            # The large-paralog-family warning must fire on stderr (>10 pre-filter loci).
            self.assertIn("WARNING", err)
            self.assertIn("large paralog family", err)
            # The vast majority of candidate loci are dropped, each logged to stdout.
            self.assertGreaterEqual(out.count("Filtered out locus"), 40)

    def test_max_loci_two_keeps_both_real_loci(self):
        """A tighter cap keeps exactly the two genuine luciferase loci."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            loci, _out, _err = _run_split(d, LUCIFERASE_HITS, extra_args=["--max_loci", "2"])
            self.assertEqual(len(loci), 2)
            self.assertTrue(_covers(loci, LUC1_CHROM, LUC1_COORD))
            self.assertTrue(_covers(loci, LUC2_CHROM, LUC2_COORD))

    def test_no_bit_column_still_capped(self):
        """6-column BEDs (no bit score) fall back to the e-value-ranked max_loci cap."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            loci, _out, _err = _run_split(d, LUCIFERASE_HITS, with_bit=False)
            self.assertLessEqual(len(loci), 5)
            # Luc1/Luc2 also have the two best e-values, so the cap still keeps them.
            self.assertTrue(_covers(loci, LUC1_CHROM, LUC1_COORD))
            self.assertTrue(_covers(loci, LUC2_CHROM, LUC2_COORD))


class TestSplitLociZeroEvalueShortCircuit(unittest.TestCase):
    """The specific bug: a 0.0 primary e-value must not let weak distant paralogs leak."""

    def _write_bed(self, path, rows):
        with open(path, "w") as f:
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")

    def test_zero_evalue_primary_does_not_leak_weak_paralog(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "home.bed")
            # Primary: e=0.0, bit 400. Distant paralog on another chrom: e=1e-15 (passes the
            # old <=1e-10 gate) but bit 60 -> ratio 0.15, must be dropped now.
            self._write_bed(bed, [
                ["chr1", 1000, 2000, "gene_loc", "0.0", "+", "400"],
                ["chr2", 5000, 6000, "gene_loc", "1e-15", "+", "60"],
            ])
            proc = subprocess.run(
                [sys.executable, SPLIT_LOCI, "--bed", bed,
                 "--output_prefix", os.path.join(d, "locus")],
                check=True, capture_output=True, text=True,
            )
            kept = [f for f in os.listdir(d) if f.startswith("locus_") and f.endswith(".bed")]
            self.assertEqual(len(kept), 1, "weak paralog leaked past a 0.0-evalue primary")
            with open(os.path.join(d, kept[0])) as f:
                self.assertIn("chr1", f.read())


if __name__ == "__main__":
    unittest.main()
