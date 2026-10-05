"""Model proteins built from miniprot CDS must be in frame across every intron phase.

Found 2026-09-11 (family benchmark): ``annotate_using_miniprot`` translated each CDS
from its first base and trimmed the GFF3 phase only on the FIRST CDS, and the runner
built the model protein by concatenating those per-exon translations. Every exon behind
a phase-1/2 intron was read in the wrong frame, so the protein SynVoy wrote for a
multi-exon GOI (regions/*.faa -> seeds, tree, paralog RBH, phylo placement) was largely
garbage while the GFF coordinates and Identity were right. A 64 %-identity Harpegnathos
hyaluronidase model came out at 47.6 %; all 113 multi-exon HIGH GOI models across the
APYR/DPP4/HYAL/KAZA/PLA2 runs had at least one exon out of frame.

The model protein is now miniprot's own translation (--trans ##STA), and each exon's
``seq`` honours its own phase.
"""
from __future__ import annotations

import random
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as ag  # noqa: E402
from sequence_utils import reverse_complement, translate  # noqa: E402

CODON = {  # one fixed codon per residue, none of them a stop
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
    'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
    'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}


def _protein(n, seed=5):
    rng = random.Random(seed)
    return 'M' + ''.join(rng.choice('ARNDCQEGHILKFPSTWYV') for _ in range(n - 1))


def _random_dna(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _gene(protein, cuts, seed=11):
    """Genomic + strand layout: flank | exon1 | GT..AG intron | exon2 | ... | flank.
    ``cuts`` are nucleotide offsets into the CDS where introns sit; offsets that are
    not multiples of 3 make phase-1/2 introns."""
    cds = ''.join(CODON[a] for a in protein) + 'TAA'
    pieces, prev = [], 0
    for c in cuts + [len(cds)]:
        pieces.append(cds[prev:c])
        prev = c
    left = _random_dna(300, seed)
    genome, coords, pos = left, [], len(left)
    for i, ex in enumerate(pieces):
        coords.append((pos, pos + len(ex)))           # 0-based half-open
        genome += ex
        pos += len(ex)
        if i < len(pieces) - 1:
            intron = 'GTAAGT' + _random_dna(388, seed + i + 1) + 'TTTCAG'
            genome += intron
            pos += len(intron)
    genome += _random_dna(300, seed + 99)
    return genome, coords, pieces


def _miniprot_gff(genome_len, coords, pieces, protein, strand, sta):
    """miniprot-style --gff --trans output for one model on ``strand``."""
    lines = ["##gff-version 3", "##PAF\tq\t0\t0\t0\t" + strand, f"##STA\t{sta}"]
    lo = min(s for s, _ in coords)
    hi = max(e for _, e in coords)
    if strand == '-':
        # the genome handed to miniprot is the reverse complement of the + layout
        coords = [(genome_len - e, genome_len - s) for s, e in coords]
        lo, hi = genome_len - hi, genome_len - lo
    lines.append(f"chr\tminiprot\tmRNA\t{lo + 1}\t{hi}\t100\t{strand}\t.\t"
                 f"ID=MP000001;Rank=1;Identity=1.0000;Target=q 1 {len(protein)}")
    # coding order = the + layout order; phase from cumulative length
    cum, q = 0, 1
    for (s, e), piece in zip(coords, pieces):
        phase = (3 - cum % 3) % 3
        aa_len = max(1, (len(piece) - phase) // 3)
        lines.append(f"chr\tminiprot\tCDS\t{s + 1}\t{e}\t100\t{strand}\t{phase}\t"
                     f"Parent=MP000001;Rank=1;Identity=1.0000;Target=q {q} {q + aa_len - 1}")
        cum += len(piece)
        q += aa_len
    return "\n".join(lines) + "\n"


@pytest.fixture
def fake_miniprot(monkeypatch):
    def install(stdout):
        monkeypatch.setattr(ag, "MINIPROT_AVAILABLE", True)
        monkeypatch.setattr(ag.subprocess, "run",
                            lambda *a, **k: SimpleNamespace(returncode=0, stdout=stdout, stderr=""))
    return install


CUTS = [100, 203]   # 100 = 33 codons + 1 -> phase-2 exon 2; 203 -> phase-1 exon 3


class TestParse:
    def test_sta_protein_attached_to_its_model(self):
        gff = ("##gff-version 3\n##PAF\tq\n##STA\tMKVL*\n"
               "c\tminiprot\tmRNA\t1\t30\t50\t+\t.\tID=MP000001;Rank=1;Identity=0.9\n"
               "c\tminiprot\tCDS\t1\t30\t50\t+\t0\tParent=MP000001;Target=q 1 10\n"
               "##PAF\tq\n##STA\tMAAA\n"
               "c\tminiprot\tmRNA\t100\t130\t40\t+\t.\tID=MP000002;Rank=2;Identity=0.8\n"
               "c\tminiprot\tCDS\t100\t130\t40\t+\t0\tParent=MP000002;Target=q 1 10\n")
        models = ag._parse_miniprot_gff(gff)
        assert [m['protein'] for m in models] == ["MKVL*", "MAAA"]

    def test_no_trans_output_gives_empty_protein(self):
        gff = ("c\tminiprot\tmRNA\t1\t30\t50\t+\t.\tID=MP000001;Rank=1;Identity=0.9\n"
               "c\tminiprot\tCDS\t1\t30\t50\t+\t0\tParent=MP000001;Target=q 1 10\n")
        assert ag._parse_miniprot_gff(gff)[0]['protein'] == ""


@pytest.mark.parametrize("strand", ["+", "-"])
class TestPhaseAwareTranslation:
    def test_every_exon_in_frame_and_protein_from_sta(self, fake_miniprot, strand):
        protein = _protein(90)
        genome, coords, pieces = _gene(protein, CUTS)
        # the layout really has phase-2 and phase-1 introns (else the test proves nothing)
        assert len(pieces[0]) % 3 == 1 and (len(pieces[0]) + len(pieces[1])) % 3 == 2
        target = genome if strand == '+' else reverse_complement(genome)
        fake_miniprot(_miniprot_gff(len(genome), coords, pieces, protein, strand, protein + "*"))

        exons, full = ag.annotate_using_miniprot(protein, target, "chr", sensitive=True)

        assert full == protein
        assert [e['phase'] for e in exons] == [0, 2, 1]
        for e in exons:
            assert e['seq'] and e['seq'] in protein, (e['exon_num'], e['seq'])

    def test_without_sta_falls_back_to_in_frame_exons(self, fake_miniprot, strand):
        protein = _protein(90)
        genome, coords, pieces = _gene(protein, CUTS)
        target = genome if strand == '+' else reverse_complement(genome)
        gff = _miniprot_gff(len(genome), coords, pieces, protein, strand, "")
        fake_miniprot("\n".join(l for l in gff.splitlines() if not l.startswith("##STA")))

        exons, full = ag.annotate_using_miniprot(protein, target, "chr", sensitive=True)

        # the two intron-split codons are the only residues lost
        assert len(full) == len(protein) - 2
        for e in exons:
            assert e['seq'] in protein

    def test_first_cds_only_trimming_would_break_frame(self, strand):
        """Pin why the fix matters: the old per-exon construction is NOT the protein."""
        protein = _protein(90)
        genome, coords, pieces = _gene(protein, CUTS)
        old = ''.join(translate(p).replace('*', '') for p in pieces)
        assert protein[:33] in old            # exon 1 is fine...
        assert protein[34:67] not in old      # ...exon 2 behind the phase-2 intron is not


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
def test_real_miniprot_phase1_and_phase2_introns():
    protein = _protein(180, seed=21)
    genome, coords, pieces = _gene(protein, [100, 302], seed=31)
    exons, full = ag.annotate_using_miniprot(protein, genome, "chr", sensitive=True)
    assert len(exons) == 3
    assert full == protein
    assert all(e['seq'] in protein for e in exons)
