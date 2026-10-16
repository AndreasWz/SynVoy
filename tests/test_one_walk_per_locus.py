"""One flanking walk per home locus (batch of 2026-10-15, HG-1B).

split_loci.py snaps the hit intervals of a kept locus to the annotated genes that cover at least
half of them (the coding span of the longest model, the hit span kept in columns 8-9), merges
intervals on one gene and writes the interval with the best copy of the gene first, named
best_interval. extract_flanking_genes.py walks once per locus from that row (else the first
row); genes at the locus's other intervals are not anchors and are written as GOI_ rows without
a protein; the anchors' goi= field still spans the hits. Until 2026-10-16 every row had a walk
of its own and a GOI pseudo-gene. All synthetic."""
import os
import random
import subprocess
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)
import split_loci as sl  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT'}


# ------------------------------------------------------------------ split_loci
def _entry(start, end, evalue, bit, chrom="c"):
    return {"chrom": chrom, "start": start, "end": end, "evalue": evalue, "bit": bit,
            "line": f"{chrom}\t{start}\t{end}\tgene_loc\t{evalue}\t+\t{bit}"}


def _hit(start, end, qs, qe, bits, chrom="c"):
    return {"chrom": chrom, "start": start, "end": end, "qstart": qs, "qend": qe, "bits": bits,
            "strand": "+", "evalue": 1e-10}


SPANS = {"c": [(1000, 1460, "gene-A"), (3000, 3460, "gene-B"), (6000, 6460, "gene-C")]}


def test_rows_snap_to_genes_and_the_best_copy_comes_first():
    cl = [_entry(3000, 3454, 1e-55, 150), _entry(1000, 1454, 1e-50, 180),
          _entry(6000, 6200, 1e-20, 80), _entry(8000, 8300, 1e-10, 50)]
    hits = [_hit(3000, 3454, 1, 150, 150), _hit(1000, 1454, 1, 150, 180),
            _hit(6000, 6200, 1, 60, 80), _hit(8000, 8300, 1, 100, 50)]
    lines = sl.locus_lines(cl, hits, SPANS)
    # gene-A's interval holds the best copy (180 bits; the best E-value is gene-B's).
    assert lines[0] == "c\t1000\t1460\tbest_interval\t1e-50\t+\t180\t1000\t1454"
    assert lines[1] == "c\t3000\t3460\tgene_loc\t1e-55\t+\t150\t3000\t3454"
    assert lines[2] == "c\t6000\t6460\tgene_loc\t1e-20\t+\t80\t6000\t6200"   # half covered: snapped
    assert lines[3] == "c\t8000\t8300\tgene_loc\t1e-10\t+\t50"             # no gene: as it was
    assert sum(1 for ln in lines if ln.split("\t")[3] == sl.BEST_INTERVAL) == 1


def test_two_intervals_on_one_gene_become_one_row():
    cl = [_entry(1000, 1200, 1e-30, 90), _entry(1250, 1454, 1e-40, 100), _entry(3000, 3454, 1e-35, 120)]
    hits = [_hit(1000, 1200, 1, 60, 90), _hit(1250, 1454, 61, 150, 100), _hit(3000, 3454, 1, 150, 120)]
    lines = sl.locus_lines(cl, hits, SPANS)
    # gene-A's two pieces are one copy of 190 bits: the best; merged with the best E and bits.
    assert lines == ["c\t1000\t1460\tbest_interval\t1e-40\t+\t100\t1000\t1454",
                     "c\t3000\t3460\tgene_loc\t1e-35\t+\t120\t3000\t3454"]


def test_a_gene_covering_less_than_half_does_not_take_the_row():
    lines = sl.locus_lines([_entry(1300, 1900, 1e-20, 60), _entry(3000, 3454, 1e-30, 90)],
                           [_hit(1300, 1900, 1, 200, 60), _hit(3000, 3454, 1, 150, 90)], SPANS)
    assert lines[1] == "c\t1300\t1900\tgene_loc\t1e-20\t+\t60"              # 160 of 600 bp


def test_one_row_has_no_marker_and_no_gff_no_snapping():
    assert sl.locus_lines([_entry(3000, 3454, 1e-55, 150)], [_hit(3000, 3454, 1, 150, 150)], SPANS) == \
        ["c\t3000\t3460\tgene_loc\t1e-55\t+\t150\t3000\t3454"]
    two = [_entry(3000, 3454, 1e-55, 150), _entry(1000, 1454, 1e-50, 180)]
    assert sl.locus_lines(two, [], {}) == ["c\t3000\t3454\tbest_interval\t1e-55\t+\t150",
                                          "c\t1000\t1454\tgene_loc\t1e-50\t+\t180"]  # by E without hits
    assert sl.coding_spans("NO_GFF", {"c"}) == {} and sl.coding_spans(None, {"c"}) == {}


# ------------------------------------------------------------------ extract_flanking
def _protein(rng, n=100):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _build(tmp_path, rng):
    """One scaffold: three flanking genes, a tandem array of three genes (T1 is the GOI, T2 and
    T3 are unrelated to it, so only the interval rule keeps them out), three more flanking
    genes. Each gene a single 303-bp CDS (100 codons + stop)."""
    names = ["L3", "L2", "L1", "T1", "T2", "T3", "R1", "R2", "R3"]
    starts = [1000, 4000, 7000, 10000, 13000, 16000, 19000, 22000, 25000]
    prots = {n: _protein(rng) for n in names}
    seq = list("".join(rng.choice("ACGT") for _ in range(28000)))
    gff = ["##gff-version 3"]
    for n, s in zip(names, starts):
        cds = "".join(CODON[a] for a in prots[n]) + "TAA"
        seq[s - 1:s - 1 + len(cds)] = cds
        e = s + len(cds) - 1
        gff += [f"c\tt\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{n};Name={n};gene_biotype=protein_coding",
                f"c\tt\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{n};Parent=gene-{n}",
                f"c\tt\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{n};Parent=rna-{n}"]
    (tmp_path / "g.fna").write_text(">c\n" + "".join(seq) + "\n")
    (tmp_path / "g.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "goi.faa").write_text(f">GOI_q\n{prots['T1']}\n")


def _flank(tmp_path, bed_rows):
    (tmp_path / "locus_2.bed").write_text("".join(r + "\n" for r in bed_rows))
    p = subprocess.run([sys.executable, os.path.join(BIN, "extract_flanking_genes.py"),
                        "--bed", str(tmp_path / "locus_2.bed"), "--gff", str(tmp_path / "g.gff"),
                        "--genome", str(tmp_path / "g.fna"), "--n_flank", "3", "--min_size", "0",
                        "--exon_mode", "true", "--goi_faa", str(tmp_path / "goi.faa"),
                        "--goi_relative_evalue", "0.00001", "--goi_relative_statistic", "sw",
                        "--out_bed", str(tmp_path / "out.bed"), "--out_faa", str(tmp_path / "out.faa")],
                       capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    bed = [ln.rstrip("\n").split("\t") for ln in open(tmp_path / "out.bed") if ln.strip()]
    faa = [ln.rstrip("\n") for ln in open(tmp_path / "out.faa") if ln.startswith(">")]
    return p.stdout, bed, faa


def test_one_walk_from_the_best_interval(tmp_path):
    _build(tmp_path, random.Random(7))
    out, bed, faa = _flank(tmp_path, [
        "c\t9999\t10302\tbest_interval\t1e-50\t+\t190\t9999\t10296",
        "c\t12999\t13302\tgene_loc\t1e-45\t+\t180\t12999\t13296",
        "c\t15999\t16302\tgene_loc\t1e-40\t+\t170\t15999\t16296"])
    assert out.count("[flank-filter] Selected") == 1
    assert "Skipping downstream gene at another hit interval of this locus: T2" in out
    anchors = [r[3] for r in bed if not r[3].startswith("GOI_")]
    assert anchors == ["gene-L1", "gene-L2", "gene-L3", "gene-R1", "gene-R2", "gene-R3"] or \
        sorted(anchors) == sorted(["gene-L1", "gene-L2", "gene-L3", "gene-R1", "gene-R2", "gene-R3"])
    goi = {r[3]: (r[1], r[2], r[6]) for r in bed if r[3].startswith("GOI_")}
    assert goi == {"GOI_c_9999": ("9999", "10302", "GOI_c_9999"),       # the walk's own pseudo-gene
                   "GOI_c_12999": ("12999", "13302", "T2"), "GOI_c_15999": ("15999", "16302", "T3")}
    assert not any(h[1:].startswith("GOI_") for h in faa)                # no protein for GOI rows
    assert not any("gene-T2" in h or "gene-T3" in h for h in faa)
    full = [h for h in faa if "full_length_protein" in h]
    assert full and all(h.endswith("goi=9999-16296") for h in full)    # the hit span


def test_a_bed_without_a_marker_walks_from_its_first_row(tmp_path):
    _build(tmp_path, random.Random(8))
    out, bed, faa = _flank(tmp_path, ["c\t9999\t10302\tgene_loc\t1e-50\t+\t190",
                                      "c\t12999\t13302\tgene_loc\t1e-45\t+\t180"])
    assert out.count("[flank-filter] Selected") == 1
    assert "gene-T2" not in [r[3] for r in bed]
    assert any(r[3] == "GOI_c_12999" and r[6] == "T2" for r in bed)
    full = [h for h in faa if "full_length_protein" in h]
    assert full and all(h.endswith("goi=9999-13302") for h in full)    # the rows' span
