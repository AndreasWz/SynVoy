"""Home hits chained into one row per gene copy (batch of 2026-10-15, HG-1A; bin/merge_hits.py).

The 1-kb merge joined hits on one strand that began less than 1 kb apart, whatever their query
coordinates: two copies 500 bp apart became one interval, and the exons of a gene with introns
over 1 kb stayed apart. Now a hit continues a chain when it continues the query
(sequence_utils.continues_query), lies within the intron bound, and the two are not in two
annotated genes (nor on both sides of one). A row covers its hits and the duplicates they
replaced (the other program's hits of the same exons), with their best E-value and bit score.
All synthetic."""
import os
import subprocess
import sys

import pytest

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)
import merge_hits as mh  # noqa: E402


def _table(path, rows):
    """rows: (chrom, qstart, qend, sstart, send, evalue, bits); 12-column outfmt 6."""
    with open(path, "w") as fh:
        for chrom, qs, qe, ss, se, ev, bits in rows:
            fh.write(f"q\t{chrom}\t90.0\t100\t5\t0\t{qs}\t{qe}\t{ss}\t{se}\t{ev}\t{bits}\n")
    return str(path)


def _gff(path, genes):
    """genes: (gene id, chrom, strand, [(cds start, cds end), ...]) 1-based."""
    lines = ["##gff-version 3"]
    for gid, chrom, strand, cds in genes:
        lo, hi = min(s for s, _e in cds), max(e for _s, e in cds)
        lines.append(f"{chrom}\tt\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID=gene-{gid};gene_biotype=protein_coding")
        lines.append(f"{chrom}\tt\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID=rna-{gid};Parent=gene-{gid}")
        for s, e in cds:
            lines.append(f"{chrom}\tt\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=cds-{gid};Parent=rna-{gid}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _run(tmp_path, tables, gff=None, extra=()):
    out = tmp_path / "home_gene_location.bed"
    cmd = [sys.executable, os.path.join(BIN, "merge_hits.py"), "--output", str(out)]
    for t in tables:
        cmd += ["--hits", t]
    if gff:
        cmd += ["--gff", gff]
    p = subprocess.run(cmd + list(extra), capture_output=True, text=True, check=True)
    rows = [ln.rstrip("\n").split("\t") for ln in open(out) if ln.strip()]
    return rows, p.stdout


def _spans(rows):
    return [(r[0], int(r[1]), int(r[2]), r[5]) for r in rows]


def test_two_exons_that_continue_the_query_are_one_copy(tmp_path):
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 101, 200, 6001, 6300, "1e-35", 140)])
    rows, out = _run(tmp_path, [t])
    assert _spans(rows) == [("c", 1000, 6300, "+")]           # 4.7 kb apart: two rows before
    assert rows[0][3] == "gene_loc" and rows[0][4] == "1e-40" and rows[0][6] == "150"
    assert "2 hit(s)" in out and "1 gene cop" in out


def test_two_full_length_copies_close_together_stay_apart(tmp_path):
    t = _table(tmp_path / "b.tsv", [("c", 1, 200, 1001, 1600, "1e-80", 300),
                                    ("c", 1, 200, 2101, 2700, "1e-70", 280)])
    rows, _ = _run(tmp_path, [t])
    assert _spans(rows) == [("c", 1000, 1600, "+"), ("c", 2100, 2700, "+")]   # one row before


def test_the_minus_strand_is_read_in_coding_order(tmp_path):
    # Exon 1 (query 1-100) lies downstream on the minus strand.
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 6300, 6001, "1e-40", 150),
                                    ("c", 101, 200, 1300, 1001, "1e-35", 140)])
    rows, _ = _run(tmp_path, [t])
    assert _spans(rows) == [("c", 1000, 6300, "-")]
    # The same two hits in genomic order would read the query backwards: no chain.
    t2 = _table(tmp_path / "b2.tsv", [("c", 101, 200, 6300, 6001, "1e-40", 150),
                                      ("c", 1, 100, 1300, 1001, "1e-35", 140)])
    rows2, _ = _run(tmp_path, [t2])
    assert len(rows2) == 2


def test_both_programs_on_one_exon_give_one_hit_with_their_span_and_best_e(tmp_path):
    mm = _table(tmp_path / "m.m8", [("c", 1, 100, 1001, 1300, "1e-50", 200)])
    bl = _table(tmp_path / "b.tsv", [("c", 1, 98, 995, 1300, "1e-60", 190)])
    rows, out = _run(tmp_path, [bl, mm])
    assert _spans(rows) == [("c", 994, 1300, "+")]           # the union of the two
    assert rows[0][4] == "1e-60" and rows[0][6] == "200"     # best E and best bits
    assert "1 hit(s)" in out


def test_a_chain_never_joins_two_annotated_genes(tmp_path):
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 101, 200, 5001, 5300, "1e-35", 140)])
    gff = _gff(tmp_path / "h.gff", [("A", "c", "+", [(1001, 1300)]), ("B", "c", "+", [(5001, 5300)])])
    assert len(_run(tmp_path, [t], gff)[0]) == 2
    assert len(_run(tmp_path, [t])[0]) == 1                   # without the annotation: one copy


def test_a_gene_between_two_hits_stops_the_chain_but_a_nested_one_does_not(tmp_path):
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 101, 200, 6001, 6300, "1e-35", 140)])
    between = _gff(tmp_path / "g1.gff", [("C", "c", "+", [(3001, 3500)])])
    assert len(_run(tmp_path, [t], between)[0]) == 2
    # Gene A holds both hits; gene N lies in A's intron.
    nested = _gff(tmp_path / "g2.gff", [("A", "c", "+", [(1001, 1300), (6001, 6300)]),
                                        ("N", "c", "-", [(3001, 3500)])])
    assert _spans(_run(tmp_path, [t], nested)[0]) == [("c", 1000, 6300, "+")]


def test_the_intron_bound(tmp_path):
    hits = [{"chrom": "c", "start": 1000, "end": 1300}]
    models = {"gene-A": {"chrom": "c", "strand": "+", "exons": [(1001, 1300, 0), (4301, 4600, 0)]}}
    assert mh.resolve_max_intron("auto", hits, models, 5.0, 10000)[0] == 15000     # 3 kb x 5
    small = {"gene-A": {"chrom": "c", "strand": "+", "exons": [(1001, 1300, 0), (1501, 1600, 0)]}}
    assert mh.resolve_max_intron("auto", hits, small, 5.0, 10000)[0] == 10000      # the floor
    huge = {"gene-A": {"chrom": "c", "strand": "+", "exons": [(1001, 1300, 0), (301_301, 301_600, 0)]}}
    assert mh.resolve_max_intron("auto", hits, huge, 5.0, 10000)[0] == mh.MAX_INTRON_CAP == 1_100_000
    elsewhere = {"gene-A": {"chrom": "d", "strand": "+", "exons": [(1001, 1300, 0), (4301, 4600, 0)]}}
    assert mh.resolve_max_intron("auto", hits, elsewhere, 5.0, 10000)[0] == 10000  # not overlapped
    assert mh.resolve_max_intron("auto", hits, {}, 5.0, 10000)[0] == 10000         # no annotation
    assert mh.resolve_max_intron("2500", hits, models, 5.0, 10000)[0] == 2500      # the user's number
    # A gap beyond the bound starts a new copy.
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 101, 200, 6001, 6300, "1e-35", 140)])
    assert len(_run(tmp_path, [t], extra=["--max_intron", "4000"])[0]) == 2
    assert len(_run(tmp_path, [t], extra=["--max_intron", "auto", "--intron_floor", "4000"])[0]) == 2


def test_an_exon_alignment_that_runs_into_the_intron_still_continues(tmp_path):
    # The next exon's hit starts 19 residues before the previous one ends (its alignment ran
    # into the intron): one copy. A hit that restarts the query 43 residues back: a new copy.
    assert mh.HOME_CHAIN_QUERY_TOL == 20
    t = _table(tmp_path / "b.tsv", [("c", 1, 120, 1001, 1360, "1e-40", 150),
                                    ("c", 101, 300, 3001, 3600, "1e-35", 140)])
    assert len(_run(tmp_path, [t])[0]) == 1
    t2 = _table(tmp_path / "b2.tsv", [("c", 1, 120, 1001, 1360, "1e-40", 150),
                                      ("c", 99, 300, 3001, 3606, "1e-35", 140)])
    assert len(_run(tmp_path, [t2])[0]) == 2                  # 21 residues: past the tolerance
    t3 = _table(tmp_path / "b3.tsv", [("c", 57, 103, 1001, 1141, "1e-10", 60),
                                      ("c", 60, 109, 3001, 3150, "1e-10", 60)])
    assert len(_run(tmp_path, [t3])[0]) == 2                  # a repeat: 43 residues


def test_a_hit_goes_to_the_chain_whose_query_end_is_closest(tmp_path):
    # Two partial copies (query 1-100 and 1-150, both on +), then a hit of query 151-250:
    # it continues the second (ends at 150), not the first (ends at 100).
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 1, 150, 2001, 2450, "1e-45", 170),
                                    ("c", 151, 250, 5001, 5300, "1e-35", 140)])
    rows, _ = _run(tmp_path, [t])
    assert _spans(rows) == [("c", 1000, 1300, "+"), ("c", 2000, 5300, "+")]


def test_no_hit_writes_an_empty_file(tmp_path):
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "0.5", 30)])
    rows, out = _run(tmp_path, [t, str(tmp_path / "missing.m8")])
    assert rows == [] and "No hit with E <= 0.01" in out


@pytest.mark.parametrize("gff_name", ["NO_GFF", "missing.gff"])
def test_without_a_usable_annotation_the_floor_applies(tmp_path, gff_name):
    t = _table(tmp_path / "b.tsv", [("c", 1, 100, 1001, 1300, "1e-40", 150),
                                    ("c", 101, 200, 9001, 9300, "1e-35", 140)])
    rows, out = _run(tmp_path, [t], str(tmp_path / gff_name))
    assert len(rows) == 1 and "no annotation: the floor of 10000 bp" in out
