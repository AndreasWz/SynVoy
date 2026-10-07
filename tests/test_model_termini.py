"""Gene models must carry the gene's own start and stop codon.

Measured 2026-09-15 against the curated annotations of the family benchmark (216 recovered
genes, 1,439 curated exons): internal splice junctions are reproduced exactly wherever the
model spans them, but the termini are systematically short, because miniprot reports only the
part of the query it aligned. Every apyrase model was short by exactly 3 codons at the N
terminus (so it had no initiator Met) and by 2 codons plus the stop codon at the C terminus,
while all ten internal junctions were byte-identical to the curated gene. Across seven arms
and ~200 models, ``StartCodon`` was set 11 times and ``StopCodon`` 121 times.

``refine_model_termini`` walks in frame outwards from the terminal CDS to close those ends,
and refuses to do so wherever the extension would be invention: at a splice site, across an
in-frame stop, past the unaligned query residues, or on a terminus that is already complete.
"""
from __future__ import annotations

import random
import shutil
import subprocess
import sys
from pathlib import Path

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
PROT = "MANDCQEGHIKF"          # prot[2] = N -> 'AAT', so no AG in front of the model start
FLANK = "CCT" * 40             # proline codons: no stop codon in any frame of the flank


def _plus_layout():
    """(sequence, cds_start, cds_end) of a single-exon gene ATG..TAA on the + strand."""
    cds = "".join(CODON[a] for a in PROT) + "TAA"
    seq = FLANK + cds + FLANK
    return seq, len(FLANK), len(FLANK) + len(cds)


def _cds(gstart, gend, phase=0, qstart=4, qend=10):
    return [{"gstart": gstart, "gend": gend, "phase": phase, "qstart": qstart, "qend": qend}]


class TestPlusStrand:
    def test_recovers_start_and_stop(self):
        seq, s, e = _plus_layout()
        cds = _cds(s + 9, e - 9)           # 3 codons short at the N end, 2 + stop at the C end
        prefix, suffix = ag.refine_model_termini(cds, seq, '+')
        assert (cds[0]["gstart"], cds[0]["gend"]) == (s, e)
        assert prefix == PROT[:3]
        assert suffix == PROT[-2:]
        assert seq[cds[0]["gstart"]:cds[0]["gstart"] + 3] == "ATG"
        assert seq[cds[0]["gend"] - 3:cds[0]["gend"]] == "TAA"

    def test_complete_model_is_left_alone(self):
        seq, s, e = _plus_layout()
        cds = _cds(s, e)
        assert ag.refine_model_termini(cds, seq, '+') == ('', '')
        assert (cds[0]["gstart"], cds[0]["gend"]) == (s, e)

    def test_splice_acceptor_blocks_the_start_walk(self):
        seq, s, e = _plus_layout()
        seq = seq[:s + 7] + "AG" + seq[s + 9:]      # AG immediately 5' of the model start
        cds = _cds(s + 9, e - 9)
        prefix, suffix = ag.refine_model_termini(cds, seq, '+')
        assert prefix == ''                          # the missing part is an upstream EXON
        assert cds[0]["gstart"] == s + 9
        assert suffix == PROT[-2:]                   # the C terminus is unaffected

    def test_splice_donor_blocks_the_stop_walk(self):
        seq, s, e = _plus_layout()
        seq = seq[:e - 9] + "GT" + seq[e - 7:]       # GT immediately 3' of the model end
        cds = _cds(s + 9, e - 9)
        prefix, suffix = ag.refine_model_termini(cds, seq, '+')
        assert suffix == ''
        assert cds[0]["gend"] == e - 9
        assert prefix == PROT[:3]

    def test_walk_refuses_to_cross_an_in_frame_stop(self):
        seq, s, e = _plus_layout()
        seq = seq[:s + 3] + "TAA" + seq[s + 6:]      # stop between the true ATG and the model
        cds = _cds(s + 9, e - 9)
        prefix, _ = ag.refine_model_termini(cds, seq, '+')
        assert prefix == ''
        assert cds[0]["gstart"] == s + 9

    def test_non_zero_phase_start_is_left_alone(self):
        seq, s, e = _plus_layout()
        cds = _cds(s + 9, e - 9, phase=1)
        prefix, _ = ag.refine_model_termini(cds, seq, '+')
        assert prefix == ''
        assert cds[0]["gstart"] == s + 9

    def test_bound_stops_a_distant_start(self):
        seq, s, e = _plus_layout()
        cds = _cds(s + 9, e - 9)
        prefix, _ = ag.refine_model_termini(cds, seq, '+', max_start_codons=2)
        assert prefix == ''                          # the ATG is 3 codons away
        assert cds[0]["gstart"] == s + 9

    def test_partial_codon_at_the_end_is_left_alone(self):
        seq, s, e = _plus_layout()
        cds = _cds(s + 9, e - 10)                    # coding length not a whole codon count
        _, suffix = ag.refine_model_termini(cds, seq, '+')
        assert suffix == ''
        assert cds[0]["gend"] == e - 10


class TestMinusStrand:
    def test_recovers_start_and_stop(self):
        fwd, s, e = _plus_layout()
        seq = reverse_complement(fwd)
        n = len(seq)
        lo, hi = n - e, n - s                        # the gene's span in the flipped sequence
        cds = _cds(lo + 9, hi - 9)                   # short at both ends, coding order is 3'->5'
        prefix, suffix = ag.refine_model_termini(cds, seq, '-')
        assert (cds[0]["gstart"], cds[0]["gend"]) == (lo, hi)
        assert prefix == PROT[:3]
        assert suffix == PROT[-2:]
        assert reverse_complement(seq[hi - 3:hi]) == "ATG"
        assert reverse_complement(seq[lo:lo + 3]) == "TAA"

    def test_splice_acceptor_blocks_the_start_walk(self):
        fwd, s, e = _plus_layout()
        seq = reverse_complement(fwd)
        n = len(seq)
        lo, hi = n - e, n - s
        edge = hi - 9
        seq = seq[:edge] + reverse_complement("AG") + seq[edge + 2:]
        cds = _cds(lo + 9, edge)
        prefix, _ = ag.refine_model_termini(cds, seq, '-')
        assert prefix == ''
        assert cds[0]["gend"] == edge


def _rand_protein(n, seed=5):
    """A non-repetitive protein: repeats make miniprot pick a spurious placement."""
    rng = random.Random(seed)
    return 'M' + ''.join(rng.choice('ARNDCQEGHILKFPSTWYV') for _ in range(n - 1))


def _multi_exon_gene(protein, cuts, seed=11):
    """flank | exon1 | GT..AG intron | exon2 | ... | flank on the + strand."""
    rng = random.Random(seed)
    rand = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    cds = "".join(CODON[a] for a in protein) + "TAA"
    pieces, prev = [], 0
    for c in list(cuts) + [len(cds)]:
        pieces.append(cds[prev:c])
        prev = c
    genome, coords = rand(300), []
    pos = len(genome)
    for i, ex in enumerate(pieces):
        coords.append((pos, pos + len(ex)))
        genome += ex
        pos += len(ex)
        if i < len(pieces) - 1:
            intron = "GTAAGT" + rand(388) + "TTTCAG"
            genome += intron
            pos += len(intron)
    return genome + rand(300), coords, pieces


@pytest.mark.skipif(shutil.which("miniprot") is None, reason="miniprot not on PATH")
class TestRealMiniprot:
    def test_model_reaches_the_start_and_stop_codon(self):
        """A query trimmed at both termini — what a diverged ortholog looks like to
        miniprot — must still yield a model with the target gene's own ATG and stop."""
        protein = _rand_protein(140)
        genome, coords, _ = _multi_exon_gene(protein, [100, 203])  # phase-2 + phase-1 introns
        query = protein[4:-3]                                      # termini unaligned
        exons, model_protein = ag.annotate_using_miniprot(
            query, genome, "chr", sensitive=True)
        assert exons, "miniprot built no model"
        assert exons[0]["gstart"] == coords[0][0], "model does not start at the gene's ATG"
        assert exons[-1]["gend"] == coords[-1][1], "model does not end at the stop codon"
        assert exons[0]["has_start_codon"] is True
        assert exons[-1]["has_stop_codon"] is True
        assert model_protein == protein, "refined model protein is not the gene's protein"

    def test_untrimmed_query_is_unchanged_by_refinement(self):
        """The same gene with its own protein as the query: the model already reaches both
        ends, so refinement must not move anything."""
        protein = _rand_protein(140, seed=7)
        genome, coords, _ = _multi_exon_gene(protein, [100, 203])
        exons, model_protein = ag.annotate_using_miniprot(
            protein, genome, "chr", sensitive=True)
        assert exons
        assert exons[0]["gstart"] == coords[0][0]
        assert exons[-1]["gend"] == coords[-1][1]
        assert model_protein == protein


class TestGffAttributes:
    """The model-level GFF attributes must say whether the model has its own termini.
    The per-exon flags existed before but only reached CDS lines, so a model that began
    three codons inside the gene was reported `complete` with nothing to show otherwise."""

    def _attrs(self, start, stop):
        sys.path.insert(0, str(ROOT / "bin"))
        import iterative_search_runner as isr
        exons = [{"has_start_codon": start, "has_stop_codon": False},
                 {"has_start_codon": False, "has_stop_codon": stop}]
        return isr._goi_feature_attrs(
            {"ID": "m1"}, evidence_type="exon_annotation", identity=70.0,
            exon_count=2, query_cov=0.9, flanking_support=8, exons=exons)

    def test_complete_termini_reported(self):
        a = self._attrs(True, True)
        assert a["StartCodon"] == "yes" and a["StopCodon"] == "yes"

    def test_missing_termini_reported(self):
        a = self._attrs(False, False)
        assert a["StartCodon"] == "no" and a["StopCodon"] == "no"

    def test_absent_exon_list_omits_the_attributes(self):
        sys.path.insert(0, str(ROOT / "bin"))
        import iterative_search_runner as isr
        a = isr._goi_feature_attrs({"ID": "m1"}, evidence_type="tandem_copy",
                                   identity=70.0, exon_count=1, query_cov=0.9,
                                   flanking_support=8)
        assert "StartCodon" not in a and "StopCodon" not in a


class TestFusedModelSupersession:
    """A model that spans several times more genome than a model of the same gene at the
    same place, while aligning barely more of the query, has chained in sequence that is
    not part of the gene. Family benchmark 2026-09-15: an 11.4-kb hull rescue (query
    coverage 0.988) outranked the 2.6-kb in-block model (0.904) that matched the curated
    gene exon for exon, because it was HIGH and the compact one MEDIUM."""

    @staticmethod
    def _ann(start, end, conf, ident, qcov, goi_class="probable_goi", source="locus_1"):
        return {"genome": "G", "chrom": "c1", "start": start, "end": end,
                "confidence": conf, "identity": ident, "goi_class": goi_class,
                "original_goi_class": goi_class, "target_gene": "gm", "source": source,
                "model_status": "complete", "query_cov": qcov, "mrna_id": f"m{start}"}

    def _dedupe(self, anns):
        sys.path.insert(0, str(ROOT / "bin"))
        import generate_report as gr
        return gr.dedupe_goi_annotations(anns)

    def test_compact_model_wins_when_the_wide_one_adds_no_coverage(self):
        wide = self._ann(245114, 256497, "HIGH", "51.7", "0.988", "synteny_hull_rescue")
        compact = self._ann(253921, 256527, "MEDIUM", "46.7", "0.904")
        out = self._dedupe([wide, compact])
        assert len(out["records"]) == 1
        rec = out["records"][0]
        assert (rec["start"], rec["end"]) == (253921, 256527)
        assert rec["superseded_wide_model"] is True
        assert out["wide_models_superseded"] == 1

    def test_wide_model_kept_when_it_really_adds_coverage(self):
        wide = self._ann(245114, 256497, "HIGH", "51.7", "0.950", "synteny_hull_rescue")
        fragment = self._ann(253921, 256527, "MEDIUM", "46.7", "0.300")
        out = self._dedupe([wide, fragment])
        assert len(out["records"]) == 1
        assert out["records"][0]["start"] == 245114
        assert out["wide_models_superseded"] == 0

    def test_unrecorded_coverage_does_not_trigger_supersession(self):
        wide = self._ann(245114, 256497, "HIGH", "51.7", "", "synteny_hull_rescue")
        compact = self._ann(253921, 256527, "MEDIUM", "46.7", "0.904")
        out = self._dedupe([wide, compact])
        assert out["records"][0]["start"] == 245114
        assert out["wide_models_superseded"] == 0

    def test_similar_spans_still_use_the_confidence_rule(self):
        a = self._ann(1000, 3000, "HIGH", "60.0", "0.90")
        b = self._ann(1010, 3010, "MEDIUM", "58.0", "0.92")
        out = self._dedupe([a, b])
        assert len(out["records"]) == 1
        assert out["records"][0]["confidence"] == "HIGH"
        assert out["wide_models_superseded"] == 0

    def test_disjoint_models_are_not_merged(self):
        a = self._ann(1000, 3000, "HIGH", "60.0", "0.90")
        b = self._ann(40000, 42000, "HIGH", "58.0", "0.88")
        out = self._dedupe([a, b])
        assert len(out["records"]) == 2


class TestRescuePathRefinement:
    """The rescue passes call miniprot directly, so they bypassed the refinement the main
    path applies. Boundary audit 2026-09-15: for 8 of 31 recovered apyrase genes the
    coordinates the report showed came from an unrefined hull-rescue model, 9-24 bp short
    of the curated gene, while the refined in-block model matched exactly."""

    @staticmethod
    def _rows(seq_len, gstart_1based, gend_1based, strand):
        mrna = ["chr", "miniprot", "mRNA", str(gstart_1based), str(gend_1based), "100",
                strand, ".", "ID=MP1;Identity=0.9"]
        cds = [["chr", "miniprot", "CDS", str(gstart_1based), str(gend_1based), "100",
                strand, "0", "Parent=MP1"]]
        return mrna, cds

    def test_plus_strand_rows_are_extended(self):
        sys.path.insert(0, str(ROOT / "bin"))
        import rescue_strong_synteny as rss
        seq, s, e = _plus_layout()
        mrna, cds = self._rows(len(seq), s + 10, e - 9, '+')   # 1-based: 3 codons in
        rss._refine_rescue_termini(seq, mrna, cds)
        assert (int(cds[0][3]), int(cds[0][4])) == (s + 1, e)  # back to the gene's own ends
        assert (int(mrna[3]), int(mrna[4])) == (s + 1, e)      # mRNA span follows
        assert seq[int(cds[0][3]) - 1:int(cds[0][3]) + 2] == "ATG"
        assert seq[int(cds[0][4]) - 3:int(cds[0][4])] == "TAA"

    def test_minus_strand_rows_are_extended(self):
        sys.path.insert(0, str(ROOT / "bin"))
        import rescue_strong_synteny as rss
        fwd, s, e = _plus_layout()
        seq = reverse_complement(fwd)
        n = len(seq)
        lo, hi = n - e, n - s
        mrna, cds = self._rows(n, lo + 10, hi - 9, '-')
        rss._refine_rescue_termini(seq, mrna, cds)
        assert (int(cds[0][3]), int(cds[0][4])) == (lo + 1, hi)
        assert reverse_complement(seq[hi - 3:hi]) == "ATG"

    def test_guards_still_apply(self):
        """A splice acceptor in front of the model blocks the walk here too."""
        sys.path.insert(0, str(ROOT / "bin"))
        import rescue_strong_synteny as rss
        seq, s, e = _plus_layout()
        seq = seq[:s + 7] + "AG" + seq[s + 9:]
        mrna, cds = self._rows(len(seq), s + 10, e - 9, '+')
        rss._refine_rescue_termini(seq, mrna, cds)
        assert int(cds[0][3]) == s + 10                        # start unchanged
        assert int(cds[0][4]) == e                             # stop still recovered

    def test_empty_input_is_a_noop(self):
        sys.path.insert(0, str(ROOT / "bin"))
        import rescue_strong_synteny as rss
        rss._refine_rescue_termini("", None, [])                # must not raise
