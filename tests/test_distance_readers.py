"""Steps after the search read the distances the search used for a genome from the
`##synvoy-distances` line of its result GFF (region scorer, hull rescue), and the steps
that read model ids understand the two id forms added on 2026-10-09
(`..._l1m2_...` = second model of a locus, `..._r1_...` = found with a relative's model).
"""
import os
import subprocess
import sys

import pytest

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, BIN)

import phylo_placement_check as ppc  # noqa: E402
import rescue_goi_hull as hull  # noqa: E402
import synvoy_scale  # noqa: E402
import synvoy_tree  # noqa: E402

CLUSTER = os.path.join(BIN, "cluster_grs.py")


def _header(**kw):
    d = {"cluster_distance": 20_000, "region_padding": 20_000, "goi_margin": 5_000,
         "synteny_bridge_max_gap": 48_000, "dispersed_envelope_margin": 22_000,
         "expansion": 1.5, "home_pitch_bp": 2_400, "home_span_bp": 46_000, "home_anchors": 20,
         "source": {}}
    d.update(kw)
    return synvoy_scale.format_gff_header(d)


def _flank(chrom, start, name, conf="HIGH"):
    return (f"{chrom}\tSynVoy\tmRNA\t{start}\t{start + 1500}\t80.0\t+\t.\t"
            f"ID={name}|t_b0_fl1_flank_ann;Name={name};SynVoyRole=flanking;Confidence={conf};Identity=80.0")


# ------------------------------------------------------------------ region scorer

def _cluster_inputs(tmp_path, header):
    """Six anchors, 12 kb apart: one block at a 20 kb gap, six pieces at 5 kb."""
    m8, gff, bed = [], ["##gff-version 3"] + ([header] if header else []), []
    for j in range(6):
        pos = 200_000 + j * 12_000
        m8.append("\t".join([f"g{j}", "chrT", "70.0", "400", "1", "0", "1", "400",
                             str(pos), str(pos + 1500), "1e-30", "200"]))
        gff.append(_flank("chrT", pos, f"g{j}"))
        bed.append(f"chrHome\t{j * 3000}\t{j * 3000 + 900}\tg{j}\t.\t+")
    (tmp_path / "hits.m8").write_text("\n".join(m8) + "\n")
    (tmp_path / "target.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "synteny.bed").write_text("\n".join(bed) + "\n")


def _run_cluster(tmp_path, *extra):
    out = tmp_path / "regions.bed"
    r = subprocess.run([sys.executable, CLUSTER, "--hits", str(tmp_path / "hits.m8"),
                        "--target_gff", str(tmp_path / "target.gff"),
                        "--synteny_bed", str(tmp_path / "synteny.bed"),
                        "--output", str(out), "--scores_output", str(tmp_path / "scores.tsv"), *extra],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1500:]
    rows = [line.split("\t") for line in out.read_text().splitlines() if line.strip()]
    return rows, r.stderr


def test_region_scorer_uses_the_block_gap_of_the_search(tmp_path):
    _cluster_inputs(tmp_path, _header(cluster_distance=20_000))
    rows, err = _run_cluster(tmp_path)
    assert "cluster distance 20000 bp (from the search)" in err
    assert len(rows) == 1
    assert int(rows[0][1]) <= 200_000 and int(rows[0][2]) >= 260_000 + 1500


def test_region_scorer_follows_a_smaller_gap_in_the_header(tmp_path):
    _cluster_inputs(tmp_path, _header(cluster_distance=5_000))
    _rows, err = _run_cluster(tmp_path)
    assert "cluster distance 5000 bp (from the search)" in err


@pytest.mark.parametrize("value", ["auto", "-1", "0"])
def test_region_scorer_without_a_header_uses_the_former_150_kb(tmp_path, value):
    _cluster_inputs(tmp_path, None)
    rows, err = _run_cluster(tmp_path, "--cluster_distance", value)
    assert "cluster distance 150000 bp (fixed: no distances in the target GFF)" in err
    assert len(rows) == 1


def test_region_scorer_takes_a_number_as_it_is(tmp_path):
    _cluster_inputs(tmp_path, _header(cluster_distance=20_000))
    _rows, err = _run_cluster(tmp_path, "--cluster_distance", "5000")
    assert "cluster distance" not in err            # nothing derived, nothing to report


# -------------------------------------------------------------------- hull rescue

class _Args:
    def __init__(self, gff, **kw):
        self.region_gff = gff
        self.cluster_max_gap = self.window_pad = self.max_window = "auto"
        self.__dict__.update(kw)


def test_hull_distances_come_from_the_search(tmp_path):
    gff = tmp_path / "r.gff"
    gff.write_text("##gff-version 3\n" + _header() + "\n")
    args = _Args(str(gff))
    hull.resolve_hull_distances(args)
    assert args.cluster_max_gap == 48_000                 # the gap a block may bridge
    assert args.window_pad == 5_000                       # the GOI margin
    assert args.max_window == 2 * 46_000 * 1.5            # two home neighbourhoods, expanded


def test_hull_distances_without_a_header_are_the_former_fixed_values(tmp_path):
    gff = tmp_path / "r.gff"
    gff.write_text("##gff-version 3\n")
    args = _Args(str(gff))
    hull.resolve_hull_distances(args)
    assert (args.cluster_max_gap, args.window_pad, args.max_window) == (3_000_000, 100_000, 20_000_000)


def test_hull_distance_set_by_the_user_wins(tmp_path):
    gff = tmp_path / "r.gff"
    gff.write_text("##gff-version 3\n" + _header() + "\n")
    args = _Args(str(gff), window_pad="100000", max_window=2e7)
    hull.resolve_hull_distances(args)
    assert (args.cluster_max_gap, args.window_pad, args.max_window) == (48_000, 100_000, 20_000_000)


def test_hull_with_bridging_switched_off_keeps_a_cluster_gap(tmp_path):
    gff = tmp_path / "r.gff"
    gff.write_text("##gff-version 3\n" + _header(synteny_bridge_max_gap=0) + "\n")
    args = _Args(str(gff))
    hull.resolve_hull_distances(args)
    assert args.cluster_max_gap == 3_000_000


def _hull_gff(tmp_path, goi_start):
    rows = ["##gff-version 3"]
    for i in range(4):
        rows.append(f"chr1\tx\tmRNA\t{1_000_000 + i * 200_000}\t{1_010_000 + i * 200_000}\t.\t+\t.\t"
                    f"ID=f{i};SynVoyRole=flanking;Confidence=HIGH")
    rows.append(f"chr1\tx\tgene\t{goi_start}\t{goi_start + 5000}\t.\t+\t.\tID=GOI_good;SynVoyRole=goi;Confidence=HIGH")
    p = tmp_path / "region.gff"
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def _run_hull(tmp_path, monkeypatch, goi_start, pad):
    (tmp_path / "q.faa").write_text(">GOI_X\n" + "MKTAYIAKQR" * 10 + "\n")
    searched = []
    monkeypatch.setattr(hull, "_extract_window", lambda genome, chrom, ws, we, samtools: searched.append((ws, we)) or "")
    monkeypatch.setattr(sys, "argv", [
        "rescue_goi_hull.py", "--region_gff", _hull_gff(tmp_path, goi_start),
        "--target_genome", str(tmp_path / "absent.fna"), "--query", str(tmp_path / "q.faa"),
        "--genome_name", "t", "--output", str(tmp_path / "out.gff"),
        "--cluster_max_gap", "3000000", "--window_pad", str(pad), "--max_window", "20000000"])
    assert hull.main() in (0, None)
    return searched


def test_goi_just_outside_the_hull_but_inside_the_searched_window_blocks_the_rescue(tmp_path, monkeypatch):
    """Hull = 1.00-1.61 Mb; a HIGH model 40 kb beyond its end lies in the DNA the rescue
    would search (hull + 100 kb). It used to be modelled a second time."""
    assert _run_hull(tmp_path, monkeypatch, goi_start=1_650_000, pad=100_000) == []


def test_goi_beyond_the_searched_window_does_not_block_the_rescue(tmp_path, monkeypatch):
    assert _run_hull(tmp_path, monkeypatch, goi_start=1_650_000, pad=10_000) == [(990_000, 1_620_000)]


# ------------------------------------------------------------------------ model ids

IDS = [
    ("GOI_Melt|Apis_florea_fna_b0_l1_exon_ann", "Apis_florea"),
    ("GOI_Melt|Apis_florea_fna_b0_l1m2_exon_ann", "Apis_florea"),
    ("GOI_Melt|Apis_florea_fna_b12_l3m10_exon_ann", "Apis_florea"),
    ("GOI_STE2|Lachancea_thermotolerans_fna_b15_r1_fallback", "Lachancea_thermotolerans"),
    ("GOI_X|Apis_florea_fa_b3_r2m3_exon_ann", "Apis_florea"),
    ("GOI_copy_1|Apis_cerana_fna_b0_l1", "Apis_cerana"),
    ("GOI_copy_1|Apis_cerana_fna_b0_l1_by_GOI_Melt", "Apis_cerana"),
    ("GOI_copy_2|Apis_cerana_fna_b1_r2_by_GOI_orf19.3940.1", "Apis_cerana"),
    ("GOI_Melt|Apis_florea_fna_b0_l1_tc2", "Apis_florea"),
    ("GOI_Melt|Apis_florea_fna_b0_l1_fallback_3", "Apis_florea"),
    ("GOI_Melt|Apis_florea_fna_b0_l1_extra3", "Apis_florea"),
    ("GOI_A|Apis_florea_fna_b2_fl1_flank_ann", "Apis_florea"),
]


@pytest.mark.parametrize("model_id, species", IDS)
def test_tree_helper_reads_the_species_from_every_id_form(model_id, species):
    assert synvoy_tree.species_from_leaf(model_id) == species


@pytest.mark.parametrize("model_id, species", IDS)
def test_placement_check_reads_the_species_from_every_id_form(model_id, species):
    assert ppc.species_from_model_id(model_id) == species
    stems = ["Apis_florea.fna", "Apis_florea.fa", "Apis_cerana.fna", "Lachancea_thermotolerans.fna"]
    resolved = ppc.species_from_model_id(model_id, known_stems=stems)
    assert resolved.replace(".", "_").startswith(species)


def test_a_strain_tag_that_looks_like_a_locus_tag_is_not_cut_from_the_species():
    """`_r1` inside the genome name stays: only the tags after the block id are cut."""
    assert synvoy_tree.species_from_leaf("GOI_X|Yeast_r1_fna_b0_l1_exon_ann") == "Yeast_r1"
