"""A home genome without its own GFF (batch of 2026-10-15, HG-5a): the effective home GFF is
the gene predictions merged with annotations borrowed from target genomes.

* One borrowed ``gene`` row switched the CDS-only fallback off, so every prediction was lost:
  CDS that belong to no feature of the file are now kept as pseudo genes in a mixed file too.
* Predictions are written with gene and mRNA parents: readers that group CDS by Parent
  (gff_to_faa.py, sequence_utils.cds_models_by_gene) dropped the bare CDS rows.
* AUGUSTUS proteins get their region and global coordinates (Prodigal ids name the region,
  AUGUSTUS ids do not).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import extract_flanking_genes as efg  # noqa: E402
import prodigal_on_regions as por  # noqa: E402
from sequence_utils import cds_models_by_gene, parse_fasta  # noqa: E402

PRODIGAL = [  # Prodigal style: one CDS row per gene, no Parent
    "chr1\tProdigal_v2.6.3\tCDS\t1001\t1600\t50.1\t+\t0\tID=1_1;partial=00;start_type=ATG",
    "chr1\tProdigal_v2.6.3\tCDS\t3001\t3900\t60.2\t-\t0\tID=1_2;partial=00;start_type=ATG",
    "chr1\tProdigal_v2.6.3\tCDS\t6001\t6300\t40.3\t+\t0\tID=1_3;partial=00;start_type=GTG",
]


def _borrowed(n):
    rows = []
    for i in range(n):
        s, e = 10001 + 2000 * i, 10900 + 2000 * i
        rows += [f"chr1\tborrowed\tgene\t{s}\t{e}\t.\t+\t.\tID=b{i};Name=borrowed{i}",
                 f"chr1\tborrowed\tmRNA\t{s}\t{e}\t.\t+\t.\tID=b{i}.t1;Parent=b{i}",
                 f"chr1\tborrowed\tCDS\t{s}\t{e}\t.\t+\t0\tID=b{i}.cds;Parent=b{i}.t1"]
    return rows


def _genes(tmp_path, rows):
    gff = tmp_path / "home.gff"
    gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return efg.parse_gff_for_genes(str(gff))


def test_predictions_survive_borrowed_gene_rows(tmp_path):
    assert len(_genes(tmp_path, PRODIGAL)) == 3                        # predictions alone
    for n_borrowed in (1, 2):
        genes = _genes(tmp_path, PRODIGAL + _borrowed(n_borrowed))
        assert len(genes) == 3 + n_borrowed, [g["id"] for g in genes]
        assert {g["id"] for g in genes} >= {"1_1", "1_2", "1_3"}
        assert all(g["cds_parts"] for g in genes)


def test_a_cds_whose_parent_no_row_defines_becomes_a_gene(tmp_path):
    rows = _borrowed(1) + ["chr1\tx\tCDS\t20001\t20300\t.\t+\t0\tID=c9;Parent=gone"]   # unattached
    genes = _genes(tmp_path, rows)
    assert sorted(g["id"] for g in genes) == ["b0", "gone"]


def test_predictions_are_written_with_gene_and_mrna_parents(tmp_path):
    raw = tmp_path / "pred.gff"
    raw.write_text("\n".join(r.replace("chr1", "chr1__1000__9000") for r in PRODIGAL) + "\n")
    out = tmp_path / "home_predicted.gff"
    por.adjust_gff(str(raw), str(out), {"chr1__1000__9000": ("chr1", 1000, 9000)})
    rows = [r.split("\t") for r in out.read_text().splitlines() if r and not r.startswith("#")]
    assert [r[2] for r in rows] == ["gene", "mRNA", "CDS"] * 3
    assert rows[0][3:5] == ["2001", "2600"] and rows[0][8] == "ID=chr1__1000__9000_1_1"
    assert rows[2][8].startswith("ID=chr1__1000__9000_1_1.cds;Parent=chr1__1000__9000_1_1.t1;partial=00")
    models = cds_models_by_gene(str(out))
    assert sorted(models) == ["chr1__1000__9000_1_1", "chr1__1000__9000_1_2", "chr1__1000__9000_1_3"]
    assert efg.parse_gff_for_genes(str(out))[0]["id"] == "chr1__1000__9000_1_1"


def test_augustus_proteins_get_their_region_and_global_coordinates(tmp_path):
    region = "chr2__5000__20000"
    gff = tmp_path / "aug.gff"
    gff.write_text(f"{region}\tAUGUSTUS\tgene\t101\t400\t0.9\t+\t.\tID=g1\n"
                   f"{region}\tAUGUSTUS\ttranscript\t101\t400\t0.9\t+\t.\tID=g1.t1;Parent=g1\n")
    faa = tmp_path / "aug.faa"
    faa.write_text(">g1 # 101 # 400 # 1 # ID=g1\nMKVLAAGIVLLW\n"
                   ">chr3__0__900_1 # 10 # 99 # -1 # ID=1_1;partial=00\nMKKLLV\n")
    out = tmp_path / "home_proteome.faa"
    id_map = {region: ("chr2", 5000, 20000), "chr3__0__900": ("chr3", 0, 900)}
    por.adjust_proteins(str(faa), str(out), id_map, por.gene_regions(str(gff)))
    headers = [h for h, _i, _s in parse_fasta(str(out))]
    assert headers[0] == f"{region}_g1 # 5101 # 5400 # 1 # ID={region}_g1"
    assert headers[1] == "chr3__0__900_1 # 10 # 99 # -1 # ID=chr3__0__900_1_1;partial=00"   # Prodigal as before
