
import unittest
import sys
import os
import shutil
import tempfile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Add bin to path
sys.path.append(os.path.abspath("bin"))
from iterative_search_runner import batch_rbh_check


@unittest.skipUnless(
    shutil.which("mmseqs"),
    "mmseqs2 binary not on PATH (CI installs only Python deps; full env via conda)",
)
class TestRBH(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.home_fa = os.path.join(self.test_dir, "home.faa")
        self.home_db = os.path.join(self.test_dir, "home_db")
        
        # Create Home DB with GeneA and GeneX
        with open(self.home_fa, 'w') as f:
            f.write(">GeneA\nMKLLLL\n")
            f.write(">GeneX\nMKVVVV\n")
            
        os.system(f"mmseqs createdb {self.home_fa} {self.home_db}")
        
    def tearDown(self):
        shutil.rmtree(self.test_dir)
        
    def test_rbh_filtering(self):
        # Case 1: GeneA -> hits GeneA (Valid)
        rec1 = SeqRecord(Seq("MKLLLL"), id="GeneA|Variant1", description="")
        
        # Case 2: GeneB -> hits GeneX (Invalid, mismatch)
        # Note: GeneB expects to find GeneB. But best hit is GeneX (seq match).
        rec2 = SeqRecord(Seq("MKVVVV"), id="GeneB|Variant1", description="")
        
        candidates = [rec1, rec2]
        unique_map = {
            "GeneA|Variant1": "GeneA",
            "GeneB|Variant1": "GeneB"
        }
        
        valid = batch_rbh_check(candidates, self.home_db, unique_map)
        
        print(f"Valid IDs: {valid}")
        self.assertIn("GeneA|Variant1", valid)
        self.assertNotIn("GeneB|Variant1", valid)

if __name__ == '__main__':
    unittest.main()
