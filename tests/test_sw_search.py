
import unittest
import os
import shutil
import tempfile
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add bin directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

# Import the function to test
from smith_waterman_search import smith_waterman_ssearch

class TestSWSearch(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.query_file = os.path.join(self.test_dir, "query.faa")
        self.target_file = os.path.join(self.test_dir, "target.fna")
        self.output_file = os.path.join(self.test_dir, "output.tsv")

        # Create a simple protein query
        self.query_seq = "MQIFVKTLTGKTITLEVEPSDTIE" * 2 # 48 AA
        SeqIO.write(SeqRecord(Seq(self.query_seq), id="query1", description="test query"), 
                   self.query_file, "fasta")

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_tandem_duplication(self):
        """
        Verify that SW finds two hits for a tandem duplication.
        Target DNA contains TWO copies of the query sequence separated by a spacer.
        """
        # Create target DNA with 2 copies
        # Translate query to DNA for embedding
        # Only works if we map back roughly. Let's just create DNA that translates to the query.
        # Simple translation: M -> ATG, Q -> CAA, etc.
        # Helper to reverse translate naively
        aa_to_dna = {
            'M': 'ATG', 'Q': 'CAA', 'I': 'ATT', 'F': 'TTT', 'V': 'GTT', 'K': 'AAA',
            'T': 'ACA', 'L': 'TTA', 'G': 'GGA', 'E': 'GAA', 'P': 'CCA', 'S': 'TCA',
            'D': 'GAT', 'A': 'GCA', 'R': 'AGA', 'N': 'AAC', 'H': 'CAT', 'C': 'TGT',
            'W': 'TGG', 'Y': 'TAT'
        }
        
        dna_seq = "".join([aa_to_dna.get(aa, 'NNN') for aa in self.query_seq])
        
        # Construct target: [spacer] [copy1] [spacer] [copy2] [spacer]
        spacer = "N" * 100
        full_target = spacer + dna_seq + spacer + dna_seq + spacer
        
        SeqIO.write(SeqRecord(Seq(full_target), id="chr1", description="tandem target"), 
                   self.target_file, "fasta")

        # convert to string path for subprocess
        smith_waterman_ssearch(self.query_file, self.target_file, self.output_file, threads=1)

        # Parse output
        hits = []
        if os.path.exists(self.output_file):
            with open(self.output_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 12:
                            hits.append(parts)
        
        # Should find at least 2 hits
        print(f"Found {len(hits)} hits")
        self.assertTrue(len(hits) >= 2, f"Expected >= 2 hits for tandem duplication, found {len(hits)}")
        
        # Verify coordinates are distinct
        # ssearch36 m8 format: q, t, pid, aln, mis, gap, qs, qe, ts, te, eval, bit
        # ts=8, te=9
        starts = sorted([int(h[8]) for h in hits])
        print(f"Hit start positions: {starts}")
        
        # Copy 1 start approx 100
        # Copy 2 start approx 100 + len(dna) + 100
        # Check separation
        self.assertTrue(starts[1] - starts[0] > 50, "Hits should be separated by spacer")

if __name__ == '__main__':
    unittest.main()
