"""The flanking walk's limit per side (batch of 2026-10-15, HG-4; bin/extract_flanking_genes.py).

Before, a side stopped after 5 x n_flank genes, whatever they were: a family larger than that,
or a run of tRNA genes, could end it short of anchors. Now a side stops only when it has
visited more than 5 x n_flank genes AND n_flank + R protein-coding genes, R being the coding
genes on that side that overlap a home hit of the query (relatives the walk will skip). So a
walk can only get longer; non-coding genes and the locus's own genes do not count. All
synthetic."""
import os
import random
import subprocess
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT'}
PITCH = 1500


def _protein(rng, n=100):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _mutate(rng, seq, frac):
    return "M" + "".join(rng.choice(AA) if rng.random() < frac else c for c in seq[1:])


def _build(tmp_path, layout, seed=1, with_hits=True):
    """layout: kinds left to right, one of 'goi', 'rel' (a relative: 80 % identity to the GOI),
    'gene' (unrelated), 'trna' (a gene without CDS). Every gene 303 bp at PITCH spacing on
    scaffold c; with_hits writes one home hit per 'goi' and 'rel' gene."""
    rng = random.Random(seed)
    goi = _protein(rng)
    seq = list("".join(rng.choice("ACGT") for _ in range(PITCH * (len(layout) + 2))))
    gff, hits, names = ["##gff-version 3"], [], []
    for i, kind in enumerate(layout):
        s = PITCH * (i + 1)
        prot = goi if kind == "goi" else _mutate(rng, goi, 0.2) if kind == "rel" else _protein(rng)
        cds = "".join(CODON[a] for a in prot) + "TAA"
        seq[s - 1:s - 1 + len(cds)] = cds
        e = s + len(cds) - 1
        name = f"{kind}{i}"
        names.append(name)
        bio = "tRNA" if kind == "trna" else "protein_coding"
        gff.append(f"c\tt\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{name};Name={name};gene_biotype={bio}")
        if kind == "trna":
            gff.append(f"c\tt\ttRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{name};Parent=gene-{name}")
            continue
        gff += [f"c\tt\tmRNA\t{s}\t{e}\t.\t+\t.\tID=rna-{name};Parent=gene-{name}",
                f"c\tt\tCDS\t{s}\t{e}\t.\t+\t0\tID=cds-{name};Parent=rna-{name}"]
        if kind in ("goi", "rel"):
            hits.append(f"q\tc\t80.0\t100\t20\t0\t1\t100\t{s}\t{e - 3}\t1e-30\t150")
    (tmp_path / "g.fna").write_text(">c\n" + "".join(seq) + "\n")
    (tmp_path / "g.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "goi.faa").write_text(f">GOI_q\n{goi}\n")
    (tmp_path / "hits.tsv").write_text("".join(h + "\n" for h in hits) if with_hits else "")
    i = layout.index("goi")
    s = PITCH * (i + 1)
    (tmp_path / "locus_1.bed").write_text(f"c\t{s - 1}\t{s + 302}\tgene_loc\t1e-30\t+\t150\n")
    return names


def _walk(tmp_path, n_flank):
    p = subprocess.run([sys.executable, os.path.join(BIN, "extract_flanking_genes.py"),
                        "--bed", str(tmp_path / "locus_1.bed"), "--gff", str(tmp_path / "g.gff"),
                        "--genome", str(tmp_path / "g.fna"), "--n_flank", str(n_flank),
                        "--min_size", "0", "--exon_mode", "true", "--goi_faa", str(tmp_path / "goi.faa"),
                        "--goi_relative_evalue", "0.00001", "--home_hits", str(tmp_path / "hits.tsv"),
                        "--max_evalue", "0.01", "--out_bed", str(tmp_path / "out.bed"),
                        "--out_faa", str(tmp_path / "out.faa")],
                       capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    anchors = [ln.split("\t")[3] for ln in open(tmp_path / "out.bed")
               if ln.strip() and not ln.split("\t")[3].startswith("GOI_")]
    return p.stdout, p.stderr, anchors


def test_a_family_larger_than_the_old_bound_no_longer_ends_the_side(tmp_path):
    # n_flank 3: the old bound is 15 genes; 20 relatives with home hits come first downstream.
    names = _build(tmp_path, ["gene"] * 3 + ["goi"] + ["rel"] * 20 + ["gene"] * 3)
    out, err, anchors = _walk(tmp_path, 3)
    downstream = [f"gene-{n}" for n in names[24:27]]
    assert all(a in anchors for a in downstream)
    assert "downstream: 3 anchor(s); visited 23 protein-coding gene(s) (20 relative(s) skipped)" in out
    assert "ended: anchors found (limit: more than 15 genes and 23 protein-coding genes, R = 20)" in out


def test_without_hits_the_old_bound_still_holds_and_can_end_a_side(tmp_path):
    names = _build(tmp_path, ["gene"] * 3 + ["goi"] + ["rel"] * 20 + ["gene"] * 3, with_hits=False)
    out, err, anchors = _walk(tmp_path, 3)
    assert not any(f"gene-{n}" in anchors for n in names[24:27])     # as before: 15 genes, then stop
    assert "downstream: 0 anchor(s)" in out and "ended: limit" in out
    assert "WARNING: the downstream walk reached its limit with 0 of 3 anchors" in err


def test_twelve_relatives_without_hits_still_reach_n(tmp_path):
    names = _build(tmp_path, ["gene"] * 10 + ["goi"] + ["rel"] * 12 + ["gene"] * 10, with_hits=False)
    out, _err, anchors = _walk(tmp_path, 10)
    assert all(f"gene-{n}" in anchors for n in names[23:33])          # R = 0, the old 50 suffice
    assert "downstream: 10 anchor(s)" in out


def test_non_coding_genes_no_longer_use_up_the_limit(tmp_path):
    names = _build(tmp_path, ["gene"] * 3 + ["goi"] + ["trna"] * 25 + ["gene"] * 3)
    out, _err, anchors = _walk(tmp_path, 3)
    assert all(f"gene-{n}" in anchors for n in names[29:32])
    assert "downstream: 3 anchor(s); visited 3 protein-coding gene(s) (0 relative(s) skipped) and 25 non-coding" in out


def test_a_short_walk_is_unchanged(tmp_path):
    names = _build(tmp_path, ["gene"] * 4 + ["goi"] + ["rel", "gene", "gene", "gene", "gene"])
    out, _err, anchors = _walk(tmp_path, 3)
    assert anchors == [f"gene-{n}" for n in names[1:4]] + [f"gene-{n}" for n in names[6:9]]
    assert "upstream: 3 anchor(s); visited 3 protein-coding gene(s) (0 relative(s) skipped)" in out
    assert out.count("[flank-filter] Selected") == 1
