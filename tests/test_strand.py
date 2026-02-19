
import unittest
import sys
import os
import tempfile
import shutil
import subprocess

class TestClusterGRS(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.script = os.path.abspath("bin/cluster_grs.py")
        
    def tearDown(self):
        shutil.rmtree(self.test_dir)
        
    def create_inputs(self, hits_data, bed_data):
        hits_file = os.path.join(self.test_dir, "hits.m8")
        bed_file = os.path.join(self.test_dir, "synteny.bed")
        out_file = os.path.join(self.test_dir, "out.bed")
        
        with open(hits_file, 'w') as f:
            for h in hits_data:
                # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
                f.write("\t".join(map(str, h)) + "\n")
                
        with open(bed_file, 'w') as f:
            for b in bed_data:
                # chr, start, end, name, score, strand
                f.write("\t".join(map(str, b)) + "\n")
                
        return hits_file, bed_file, out_file

    def run_script(self, hits, bed, out):
        cmd = [
            sys.executable, self.script,
            "--hits", hits,
            "--synteny_bed", bed,
            "--flanking_count", "3",
            "--output", out,
            "--min_score", "0.1" 
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def test_perfect_synteny(self):
        # Home: A(+) B(+) C(+)
        bed = [
            ["c1", 100, 200, "A", 0, "+"],
            ["c1", 300, 400, "B", 0, "+"],
            ["c1", 500, 600, "C", 0, "+"]
        ]
        # Target: A(+) B(+) C(+)
        # Hits: tstart < tend for + strand
        hits = [
            ["A", "t1", 100, 100, 0, 0, 0, 0, 1000, 1100, 0, 100],
            ["B", "t1", 100, 100, 0, 0, 0, 0, 1200, 1300, 0, 100],
            ["C", "t1", 100, 100, 0, 0, 0, 0, 1400, 1500, 0, 100]
        ]
        
        h, b, o = self.create_inputs(hits, bed)
        res = self.run_script(h, b, o)
        
        print("Perfect Output:", res.stdout)
        self.assertIn("score: 1.00", res.stdout.lower())
        self.assertIn("consistency: 1.00", res.stdout.lower())
        self.assertIn("strand: 1.00", res.stdout.lower())

    def test_inverted_synteny(self):
        # Home: A(+) B(+) C(+)
        bed = [
            ["c1", 100, 200, "A", 0, "+"],
            ["c1", 300, 400, "B", 0, "+"],
            ["c1", 500, 600, "C", 0, "+"]
        ]
        # Target: C(-) B(-) A(-)
        # Hits: tstart > tend for - strand
        # Order on genome: A(3000) B(2000) C(1000)? No.
        # Genome: ... C ... B ... A ... 
        # C is at 1000, B at 1200, A at 1400 ?? 
        # If Inverted: A should be at HIGH coord, C at LOW coord?
        # A, B, C order in query.
        # Inverted Target: C, B, A order.
        # C at 1000, B at 1200, A at 1400.
        # But strands are minus.
        
        hits = [
            ["C", "t1", 100, 100, 0, 0, 0, 0, 1100, 1000, 0, 100], # - strand
            ["B", "t1", 100, 100, 0, 0, 0, 0, 1300, 1200, 0, 100], # - strand
            ["A", "t1", 100, 100, 0, 0, 0, 0, 1500, 1400, 0, 100]  # - strand
        ]
        
        h, b, o = self.create_inputs(hits, bed)
        res = self.run_script(h, b, o)
        
        print("Inverted Output:", res.stdout)
        # Unique: 3/3 (1.0)
        # Consistency: 
        # Order indices: C(2) -> B(1) -> A(0).
        # Diffs: |2-1|=1, |1-0|=1. Consistent!
        # Strand: All (-). Expected (+). Ratio 0/3 Same?
        # But Inverted logic: We take MAX(Same, Diff).
        # Diff = 3/3 = 1.0.
        # Should be perfect score.
        
        self.assertIn("score: 1.00", res.stdout.lower())
        
    def test_bad_strand(self):
        # Home: A(+) B(+) C(+)
        bed = [
            ["c1", 100, 200, "A", 0, "+"],
            ["c1", 300, 400, "B", 0, "+"],
            ["c1", 500, 600, "C", 0, "+"]
        ]
        # Target: Mixed
        hits = [
            ["A", "t1", 100, 100, 0, 0, 0, 0, 1000, 1100, 0, 100], # +
            ["B", "t1", 100, 100, 0, 0, 0, 0, 1300, 1200, 0, 100], # - (BAD)
            ["C", "t1", 100, 100, 0, 0, 0, 0, 1400, 1500, 0, 100]  # +
        ]
        
        h, b, o = self.create_inputs(hits, bed)
        res = self.run_script(h, b, o)
        
        print("Mixed Output:", res.stdout)
        # Strand cons: 2 same, 1 diff.
        # Score should be lower than 1.0.
        
        # Verify score is < 1.0
        # Expecting something like 0.8 or 0.7
        self.assertFalse("score: 1.00" in res.stdout.lower())

if __name__ == '__main__':
    unittest.main()
