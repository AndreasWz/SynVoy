"""Fix-phase item 4: the melittin scorer counts decoy models as calls (score_coordinates.py).

A decoy-mode run writes the models of its reversed and shuffled queries to decoys/. The scorer took
every *.gff* under the run, so they read as calls: Round 2's arm A hit 6 truth models with them and 3
without them. Decoy models are the search's null test, wherever they stand, and never calls.
"""
import importlib.util
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    "score_coordinates_decoys", os.path.join(ROOT, "scripts", "benchmark", "score_coordinates.py"))
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)


def _row(chrom, start, end, attrs, kind="gene"):
    return "\t".join([chrom, "tandem_copy", kind, str(start), str(end), "50", "+", ".", attrs]) + "\n"


def _goi(chrom, start, end, ident, parent):
    return _row(chrom, start, end, f"ID={ident};SynVoy_Parent={parent};SynVoyRole=goi;Confidence=LOW")


def _write(path, *rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.writelines(rows)


def _truth(species="Sp_one", chrom="c1", start=1000, end=2000):
    return [{"species": species, "chrom": chrom, "start": start, "end": end, "model_id": "m1",
             "truth_kind": "gene", "confidence": "curated", "chrom_keys": {chrom}}]


@pytest.fixture
def run(tmp_path):
    """A decoy-mode run: one real model away from the truth, one decoy model on it."""
    r = tmp_path / "run"
    _write(str(r / "plot_inputs_synteny_block_locus_1" / "Sp_one.fna.gff"),
           _goi("c1", 5000, 5100, "GOI_x|Sp_one_fna_b0_l1_tc1", "GOI_x"))
    _write(str(r / "decoys" / "locus_1" / "locus_1_Sp_one.fna.decoy.gff"),
           _goi("c1", 1200, 1400, "GOI_x~rev|Sp_one_fna_b0_l1_tc1", "GOI_x~rev"))
    return str(r)


def test_decoy_folder_is_not_gathered(run):
    paths = sc.gather_gffs(run)
    assert len(paths) == 1 and os.sep + "decoys" + os.sep not in paths[0]
    assert len(sc.gather_gffs(run, decoys=True)) == 2


def test_only_a_folder_inside_the_run_counts(tmp_path):
    """A run whose own path holds 'decoys' is scored whole."""
    r = tmp_path / "decoys" / "b27_decoys"
    _write(str(r / "regions" / "Sp_one.fna.gff"), _goi("c1", 1, 2, "GOI_x|a", "GOI_x"))
    assert len(sc.gather_gffs(str(r))) == 1


def test_a_run_scores_without_its_decoy_models(run):
    calls = sc.load_calls_from_gff(sc.gather_gffs(run))
    assert [c["id"] for c in calls] == ["GOI_x|Sp_one_fna_b0_l1_tc1"]
    _, _, summary = sc.score(_truth(), calls, searched_chroms={"c1"})
    assert summary["truth_models_hit"] == 0 and summary["calls_on_truth"] == 0


def test_the_old_gather_would_have_hit_the_truth(run):
    """What the fix removes: the same run with the decoy file read as calls hits the truth model."""
    decoy = [p for p in sc.gather_gffs(run, decoys=True) if "decoys" in p]
    rows = open(decoy[0]).read().replace("~rev", "")
    real = os.path.join(os.path.dirname(os.path.dirname(decoy[0])), "copy.gff")
    with open(real, "w") as fh:
        fh.write(rows)
    calls = sc.load_calls_from_gff([real])
    _, _, summary = sc.score(_truth(), calls, searched_chroms={"c1"})
    assert summary["truth_models_hit"] == 1


@pytest.mark.parametrize("ident,parent,decoy", [
    ("GOI_x~rev|g_b0_l1", "GOI_x~rev", True),
    ("GOI_x~shuf1|g_b0_l1", "GOI_x~shuf1", True),
    ("GOI_x~shuf2|g_b0_r1_tc2", "GOI_x~shuf2", True),
    ("GOI_x~rev~rel|g_b0_l1", "GOI_x~rev~rel", True),
    ("GOI_copy_1|g_b0_l1", "GOI_x~rev", True),
    ("GOI_x|g_b0_l1", "GOI_x", False),
    ("GOI_x~rel|g_b0_r2", "GOI_x~rel", False),
    ("GOI_review|g_b0_l1", "GOI_review", False),
    ("GOI_shuffle-1|g_b0_l1", "GOI_shuffle-1", False),
])
def test_decoy_ids(ident, parent, decoy):
    assert sc.is_decoy_model({"ID": ident, "SynVoy_Parent": parent}) is decoy


def test_a_decoy_model_outside_the_folder_is_no_call(tmp_path):
    """A single GFF, or a leak into the result GFFs: the model id decides, not only the folder."""
    p = str(tmp_path / "one.gff")
    _write(p, _goi("c1", 1200, 1400, "GOI_x~shuf2|a", "GOI_x~shuf2"), _goi("c1", 1200, 1400, "GOI_x~rel|a", "GOI_x~rel"))
    assert [c["id"] for c in sc.load_calls_from_gff(sc.gather_gffs(p))] == ["GOI_x~rel|a"]


def test_a_scaffold_seen_only_in_decoys_is_a_real_miss(tmp_path):
    """A decoy model's scaffold proves the run used the curated assembly: a truth model there is a
    miss. Reading the folder for the calls only would call it an assembly mismatch and drop it
    from the denominator."""
    r = tmp_path / "run"
    (r / "regions").mkdir(parents=True)
    (r / "regions" / "Sp_one.fna.regions.bed").write_text("")
    _write(str(r / "regions" / "Sp_one.fna.gff"), _goi("other", 5000, 5100, "GOI_x|a", "GOI_x"))
    _write(str(r / "decoys" / "locus_1" / "locus_1_Sp_one.fna.decoy.gff"),
           _goi("c1", 5000, 5100, "GOI_x~rev|a", "GOI_x~rev"))
    calls = sc.load_calls_from_gff(sc.gather_gffs(str(r)))
    searched = sc.load_searched_chroms(sc.gather_gffs(str(r), decoys=True))
    _, _, summary = sc.score(_truth(), calls, searched_chroms=searched,
                             searched_species=sc.load_searched_species(str(r)))
    assert summary["truth_models_real"] == 1 and summary["truth_models_assembly_mismatch"] == 0
    assert summary["truth_models_hit"] == 0
    # The narrower set is what main() no longer passes:
    narrow = sc.load_searched_chroms(sc.gather_gffs(str(r)))
    _, _, summary = sc.score(_truth(), calls, searched_chroms=narrow,
                             searched_species=sc.load_searched_species(str(r)))
    assert summary["truth_models_assembly_mismatch"] == 1
