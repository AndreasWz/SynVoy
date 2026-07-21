"""Tests for bin/stage_target_gffs.py — Pro-Mode target-annotation staging.

The load-bearing property is the DESTINATION FILENAME. Two independent
consumers probe for an annotation next to the genome:

  * iterative_search_runner.find_native_annotation_path() — its first candidate
    is Path(genome).with_suffix('.gff'); a hit stamps TargetGene/TargetProduct
    onto every model and feeds the family-consistency gate.
  * borrow_annotations.find_annotated_target() — probes fna.with_suffix('.gff').

...while cluster_regions.nf resolves the genome with a `find -name
"<genome_name>*"` glob where <genome_name> INCLUDES the extension. So the name
must simultaneously be found by the first two and invisible to the third.
test_staged_name_is_invisible_to_cluster_regions_glob pins that.
"""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import stage_target_gffs as sts  # noqa: E402


# Mirrors RefSeq layout: the gene= attribute is repeated on the mRNA/CDS rows,
# which is what load_native_annotation_index() reads for TargetGene.
GFF_BODY = """##gff-version 3
chr1\tRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-DCN;Name=DCN;gene=DCN
chr1\tRefSeq\tmRNA\t1000\t2000\t.\t+\t.\tID=rna-DCN;Parent=gene-DCN;gene=DCN;product=decorin
chr1\tRefSeq\tCDS\t1000\t2000\t.\t+\t0\tID=cds-DCN;Parent=rna-DCN;gene=DCN;product=decorin
"""


def _mkgenome(d: Path, name: str) -> Path:
    p = d / name
    p.write_text(">chr1\nACGTACGTAC\n")
    return p


def _mkgff(d: Path, name: str, body: str = GFF_BODY) -> Path:
    p = d / name
    p.write_text(body)
    return p


def _run(genomes_dir: Path, gffs, extra=None):
    cmd = [sys.executable, str(BIN / "stage_target_gffs.py"),
           "--genomes_dir", str(genomes_dir)]
    if extra:
        cmd += extra
    cmd += ["--gffs"] + [str(g) for g in gffs]
    return subprocess.run(cmd, capture_output=True, text=True)


# --- naming contract -------------------------------------------------------

@pytest.mark.parametrize("genome,expected", [
    ("cow.fna", "cow.gff"),
    ("cow.fa", "cow.gff"),
    ("cow.fasta", "cow.gff"),
    ("cow.fna.gz", "cow.fna.gff"),
    ("GCF_002263795.1_ARS-UCD1.2_genomic.fna",
     "GCF_002263795.1_ARS-UCD1.2_genomic.gff"),
])
def test_native_annotation_name(genome, expected):
    assert sts.native_annotation_name(genome) == expected


def test_staged_name_matches_search_layer_probe():
    """The name we write must be what find_native_annotation_path() looks for."""
    sys.path.insert(0, str(BIN))
    from iterative_search_runner import find_native_annotation_path

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        genome = _mkgenome(d, "cow.fna")
        (d / sts.native_annotation_name(genome.name)).write_text(GFF_BODY)
        assert find_native_annotation_path(str(genome)) == str(d / "cow.gff")


def test_staged_name_is_invisible_to_cluster_regions_glob(tmp_path):
    """cluster_regions.nf globs '<genome_name>*' where genome_name='cow.fna'.

    Naming the annotation cow.fna.gff would make that glob able to return the
    GFF instead of the genome. cow.gff must not match.
    """
    genome_name = "cow.fna"
    staged = sts.native_annotation_name(genome_name)
    assert not staged.startswith(genome_name), (
        f"{staged} would be matched by find -name '{genome_name}*'"
    )


# --- matching --------------------------------------------------------------

def test_exact_stem_match(tmp_path):
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    _mkgenome(gdir, "mouse.fna")
    gffs = [_mkgff(adir, "cow.gff"), _mkgff(adir, "mouse.gff")]

    r = _run(gdir, gffs)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (gdir / "cow.gff").exists()
    assert (gdir / "mouse.gff").exists()


def test_accession_match_across_differing_suffixes(tmp_path):
    """Genome GCF_x.fna vs annotation GCF_x_ARS-UCD1.2_genomic.gff."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "GCF_002263795.1.fna")
    gff = _mkgff(adir, "GCF_002263795.1_ARS-UCD1.2_genomic.gff")

    r = _run(gdir, [gff])
    assert r.returncode == 0, r.stdout + r.stderr
    assert (gdir / "GCF_002263795.1.gff").exists()


def test_gzipped_gff_is_decompressed(tmp_path):
    """borrow_annotations only probes plain .gff, so .gz must be expanded."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    gz = adir / "cow.gff.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(GFF_BODY)

    r = _run(gdir, [gz])
    assert r.returncode == 0, r.stdout + r.stderr
    dest = gdir / "cow.gff"
    assert dest.exists() and not dest.is_symlink()
    assert "gene-DCN" in dest.read_text()


def test_unmatched_gff_is_fatal_by_default(tmp_path):
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    gff = _mkgff(adir, "chicken.gff")

    r = _run(gdir, [gff])
    assert r.returncode == 1
    assert "no genome matches" in r.stdout
    assert not (gdir / "chicken.gff").exists()


def test_unmatched_gff_tolerated_with_flag(tmp_path):
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    _mkgenome(gdir, "mouse.fna")
    gffs = [_mkgff(adir, "cow.gff"), _mkgff(adir, "chicken.gff")]

    r = _run(gdir, gffs, extra=["--allow_unmatched"])
    assert r.returncode == 0, r.stdout + r.stderr
    assert (gdir / "cow.gff").exists()


def test_partial_annotation_is_allowed(tmp_path):
    """Annotating only some targets is legitimate and must not fail."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    _mkgenome(gdir, "mouse.fna")
    gff = _mkgff(adir, "cow.gff")

    r = _run(gdir, [gff])
    assert r.returncode == 0, r.stdout + r.stderr
    assert (gdir / "cow.gff").exists()
    assert not (gdir / "mouse.gff").exists()
    assert "without annotation" in r.stdout


def test_sentinel_is_a_noop(tmp_path):
    """The workflow passes assets/sentinels/NO_GFFS when --target_gffs is unset."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    _mkgenome(gdir, "cow.fna")
    sentinel = ROOT / "assets" / "sentinels" / "NO_GFFS"

    r = _run(gdir, [sentinel])
    assert r.returncode == 0, r.stdout + r.stderr
    assert not (gdir / "cow.gff").exists()
    assert "No target GFFs supplied" in r.stdout


def test_no_gffs_at_all_is_a_noop(tmp_path):
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    _mkgenome(gdir, "cow.fna")
    r = _run(gdir, [])
    assert r.returncode == 0, r.stdout + r.stderr


def test_double_claim_is_reported(tmp_path):
    """Two annotations resolving to one genome must not silently overwrite."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    a1 = tmp_path / "a"
    a1.mkdir()
    a2 = tmp_path / "b"
    a2.mkdir()
    _mkgenome(gdir, "cow.fna")
    gffs = [_mkgff(a1, "cow.gff"), _mkgff(a2, "cow.gff3")]

    r = _run(gdir, gffs)
    assert r.returncode == 1
    assert "already annotated" in r.stdout


def test_featureless_gff_is_flagged(tmp_path):
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    _mkgenome(gdir, "cow.fna")
    gff = _mkgff(adir, "cow.gff", body="##gff-version 3\n")

    r = _run(gdir, [gff])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "no usable features" in r.stdout


def test_species_map_is_unaffected_by_staged_gffs(tmp_path):
    """GFFs land in the genomes dir; list_genomes must still see only FASTAs."""
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    _mkgenome(gdir, "cow.fna")
    _mkgff(gdir, "cow.gff")
    assert [p.name for p in sts.list_genomes(gdir)] == ["cow.fna"]


# --- annotation actually reaches the search layer --------------------------

def test_staged_annotation_is_indexable(tmp_path):
    """End-to-end: staged file parses into the native-annotation index that
    supplies TargetGene/TargetProduct."""
    from iterative_search_runner import (find_native_annotation_path,
                                         load_native_annotation_index,
                                         lookup_native_annotation)
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    adir = tmp_path / "gffs"
    adir.mkdir()
    genome = _mkgenome(gdir, "cow.fna")
    gff = _mkgff(adir, "cow.gff")

    r = _run(gdir, [gff])
    assert r.returncode == 0, r.stdout + r.stderr

    found = find_native_annotation_path(str(genome))
    assert found is not None
    index = load_native_annotation_index(found)
    hit = lookup_native_annotation(index, "chr1", 1200, 1500, "+")
    assert hit is not None
    assert hit["label"] == "DCN"
    assert hit["product"] == "decorin"
