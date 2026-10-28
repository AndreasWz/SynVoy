"""An Augustus model that is not installed stops the gene prediction (fix-phase item 16, 2026-10-28).

gene_predictor.AUGUSTUS_SPECIES_MAP mapped three hints to models Augustus does not have
(bombus -> bombus_terrestris1, silkworm -> bombyx_mori, mouse -> mouse). Augustus then exited 1,
predict_genes fell back to Prodigal with one warning, and the task went on with Prodigal's
intron-blind models (measured on a 40 kb piece of the Candida home genome: 60 Prodigal genes
instead of 41 Augustus ones). Now resolve_augustus_species checks the model's parameter file:
a model Augustus cannot load raises AugustusModelNotInstalled with the hint, the model and the
nearest installed models, and both prediction scripts stop with that message before they predict;
the runner stops when a step that runs Augustus on target windows is on. `bombus` maps to the
installed bombus_terrestris2; aliases and models the user installed pass.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import gene_predictor as gp  # noqa: E402

AUGUSTUS = shutil.which("augustus")
real_augustus = pytest.mark.skipif(not AUGUSTUS or not gp.augustus_config_dir(),
                                   reason="Augustus not installed")


def _config(tmp_path: Path, models=(), loose=()) -> str:
    """A fake Augustus config folder: species/<m>/<m>_parameters.cfg for each of ``models``,
    <m>_parameters.cfg in the folder itself for each of ``loose``."""
    cfg = tmp_path / "config"
    (cfg / "species").mkdir(parents=True)
    for m in models:
        (cfg / "species" / m).mkdir()
        (cfg / "species" / m / f"{m}_parameters.cfg").write_text("#\n")
    for m in loose:
        (cfg / f"{m}_parameters.cfg").write_text("#\n")
    return str(cfg)


@pytest.fixture
def fake_augustus(tmp_path, monkeypatch):
    """Augustus 'installed' with a fake config folder holding fly, bombus_terrestris2,
    bombus_impatiens1, heliconius_melpomene1 and a user's model my_ant_model; no taxonomy."""
    cfg = _config(tmp_path, models=("fly", "bombus_terrestris2", "bombus_impatiens1",
                                    "heliconius_melpomene1", "my_ant_model"))
    (Path(cfg) / "species" / "half_model").mkdir()           # a folder without its parameters
    monkeypatch.setattr(gp, "check_augustus_available", lambda: True)
    monkeypatch.setattr(gp, "augustus_config_dir", lambda: cfg)
    monkeypatch.setattr(gp, "installed_augustus_models", lambda help_text=None, config_dir=None: {
        "fly": "Drosophila melanogaster", "bombus_terrestris2": "Bombus terrestris",
        "bombus_impatiens1": "Bombus impatiens", "heliconius_melpomene1": "Heliconius melpomene"})
    monkeypatch.setattr(gp, "offline_ncbi_taxa", lambda: None)
    return cfg


def test_the_parameter_file_decides(tmp_path):
    cfg = _config(tmp_path, models=("fly",), loose=("loose_model",))
    (Path(cfg) / "species" / "empty").mkdir()
    assert gp.augustus_model_installed("fly", cfg)
    assert gp.augustus_model_installed("loose_model", cfg)     # where Augustus looks next
    assert not gp.augustus_model_installed("empty", cfg)       # a folder is not a model
    assert not gp.augustus_model_installed("bombus_terrestris1", cfg)
    assert not gp.augustus_model_installed("../config/species/fly", cfg)
    assert not gp.augustus_model_installed("", cfg)


def test_installed_names_and_hints_resolve(fake_augustus):
    assert gp.resolve_augustus_species("bombus") == "bombus_terrestris2"   # the corrected entry
    assert gp.resolve_augustus_species("fly") == "fly"
    assert gp.resolve_augustus_species("Drosophila") == "fly"
    assert gp.resolve_augustus_species("") == "fly"
    assert gp.resolve_augustus_species("my_ant_model") == "my_ant_model"   # a model the user installed
    assert gp.resolve_augustus_species(" fly ") == "fly"


@pytest.mark.parametrize("hint, model", [("silkworm", "bombyx_mori"), ("mouse", "mouse"),
                                         ("half_model", "half_model"), ("no_such_model", "no_such_model")])
def test_a_model_that_is_not_installed_stops(fake_augustus, hint, model):
    with pytest.raises(gp.AugustusModelNotInstalled) as err:
        gp.resolve_augustus_species(hint)
    msg = str(err.value)
    assert f"augustus_species '{hint}'" in msg
    assert f"{model}_parameters.cfg" in msg
    if model != hint:
        assert f"maps to the Augustus model '{model}'" in msg
    assert "augustus --species=help" in msg


def test_the_message_names_the_nearest_models_by_spelling(fake_augustus):
    with pytest.raises(gp.AugustusModelNotInstalled) as err:
        gp.resolve_augustus_species("bombus_terrestris1")
    assert "Nearest installed models: bombus_terrestris2 (Bombus terrestris)" in str(err.value)
    assert ("4 models are installed: the identifiers `augustus --species=help` lists (not its "
            "aliases in parentheses) and the Genus_species folders") in str(err.value)


class FakeTaxonomy:
    LINEAGES = {"Bombyx mori": [1, 2, 50, 70, 7091], "Heliconius melpomene": [1, 2, 50, 70, 34740],
                "Bombus terrestris": [1, 2, 50, 60, 30195], "Bombus impatiens": [1, 2, 50, 60, 132113],
                "Drosophila melanogaster": [1, 2, 50, 80, 7227]}

    def lineage(self, name):
        key = " ".join(str(name).replace("_", " ").split())
        key = next((k for k in self.LINEAGES if k.lower() == key.lower()), None)
        return {"lineage": self.LINEAGES[key], "matched": key} if key else None

    def describe(self, taxid):
        return {70: ("Obtectomera", "no rank"), 60: ("Bombus", "genus")}.get(taxid, ("", ""))


def test_the_nearest_model_by_taxon():
    models = {"fly": "Drosophila melanogaster", "heliconius_melpomene1": "Heliconius melpomene",
              "bombus_terrestris2": "Bombus terrestris", "bombus_impatiens1": "Bombus impatiens"}
    near = gp.nearest_installed_models("silkworm", "bombyx_mori", models, FakeTaxonomy())
    assert near[0].startswith("heliconius_melpomene1 (Heliconius melpomene; nearest by taxon, "
                              "shared taxon Obtectomera)")
    assert gp.nearest_installed_models("mouse", "mouse", models, None) == []
    # a model id's version digits are not part of its species: the same species comes first
    near = gp.nearest_installed_models("bombus_terrestris1", "bombus_terrestris1", models, FakeTaxonomy())
    assert near[0].startswith("bombus_terrestris2 (Bombus terrestris; nearest by taxon")


def test_without_augustus_nothing_is_checked(monkeypatch):
    monkeypatch.setattr(gp, "check_augustus_available", lambda: False)
    assert gp.resolve_augustus_species("silkworm") == "bombyx_mori"


def test_the_config_folder_is_the_one_augustus_reads(tmp_path, monkeypatch):
    """With $AUGUSTUS_CONFIG_PATH set, Augustus reads that folder only, also when it holds no
    species folder; the folder next to the binary counts only when the variable is unset."""
    prefix = tmp_path / "prefix"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "config" / "species" / "fly").mkdir(parents=True)
    (prefix / "config" / "species" / "fly" / "fly_parameters.cfg").write_text("#\n")
    monkeypatch.setattr(gp.shutil, "which", lambda name: str(prefix / "bin" / name))
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    assert gp.augustus_config_dir() == str(prefix / "config")
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(empty))
    assert gp.augustus_config_dir() == str(empty)               # not the folder next to the binary
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", "")
    assert gp.augustus_config_dir() == "/"                      # Augustus then looks in /species
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(empty))
    monkeypatch.setattr(gp, "check_augustus_available", lambda: True)
    monkeypatch.setattr(gp, "installed_augustus_models", lambda help_text=None, config_dir=None: {})
    monkeypatch.setattr(gp, "offline_ncbi_taxa", lambda: None)
    with pytest.raises(gp.AugustusModelNotInstalled) as err:
        gp.resolve_augustus_species("fly")
    assert f"only from $AUGUSTUS_CONFIG_PATH ({empty})" in str(err.value)


def test_the_runner_stops_only_when_a_step_runs_augustus(fake_augustus, monkeypatch):
    """Augustus runs on target windows only in the structural discovery, which runs with the PLM
    or the structural search on; either search alone, or discovery alone, runs no Augustus."""
    import iterative_search_runner as runner
    off = argparse.Namespace(enable_plm_search=False, enable_structural_search=False,
                             enable_structural_discovery=False, augustus_species="silkworm")

    def args(**kw):
        return argparse.Namespace(**{**vars(off), **kw})

    assert runner.augustus_steps_model_error(off) is None        # the search never runs Augustus
    for step in ("enable_plm_search", "enable_structural_search", "enable_structural_discovery"):
        assert runner.augustus_steps_model_error(args(**{step: True})) is None
    for search in ("enable_plm_search", "enable_structural_search"):
        on = args(**{search: True}, enable_structural_discovery=True)
        assert "bombyx_mori" in runner.augustus_steps_model_error(on)
        assert runner.augustus_steps_model_error(
            argparse.Namespace(**{**vars(on), "augustus_species": "bombus"})) is None
    monkeypatch.setattr(runner, "GPU_IN_WORKER", True)           # the window searches in the worker
    assert "bombyx_mori" in runner.augustus_steps_model_error(args(enable_plm_search=True))
    main = (BIN / "iterative_search_runner.py").read_text().split("def main(", 1)[1]
    assert "augustus_steps_model_error(args)" in main.split("PLM search validation", 1)[0]


@real_augustus
def test_the_map_names_installed_models_but_silkworm_and_mouse():
    cfg = gp.augustus_config_dir()
    missing = sorted(h for h, m in gp.AUGUSTUS_SPECIES_MAP.items() if not gp.augustus_model_installed(m, cfg))
    assert missing == ["mouse", "silkworm"]
    assert gp.resolve_augustus_species("bombus") == "bombus_terrestris2"
    assert gp.resolve_augustus_species("yeast") == "saccharomyces"    # an alias of the help list
    with pytest.raises(gp.AugustusModelNotInstalled):
        gp.resolve_augustus_species("silkworm")


@real_augustus
@pytest.mark.parametrize("script", ["prodigal_on_regions.py", "extract_flanking_genes.py"])
def test_both_prediction_scripts_stop_before_they_predict(tmp_path, script):
    """Both scripts stop with the message (exit 1) before they predict (prodigal_on_regions.py
    before it reads the genome); with --gene_predictor prodigal the hint is not used."""
    import random
    rng = random.Random(16)
    genome = tmp_path / "g.fna"
    genome.write_text(">c1\n" + "".join(rng.choice("ACGT") for _ in range(4000)) + "\n")
    bed = tmp_path / "goi.bed"
    bed.write_text("c1\t1500\t1800\tbest_interval\t1e-30\t+\n")
    if script == "prodigal_on_regions.py":
        base = [sys.executable, str(BIN / script), "--genome", str(genome), "--goi_bed", str(bed),
                "--window", "1000", "--output_faa", str(tmp_path / "o.faa"),
                "--output_gff", str(tmp_path / "o.gff")]
    else:
        base = [sys.executable, str(BIN / script), "--bed", str(bed), "--gff", "NO_GFF",
                "--genome", str(genome), "--pred_flank_window", "1000",
                "--out_bed", str(tmp_path / "o.bed"), "--out_faa", str(tmp_path / "o.faa")]
    stop = subprocess.run(base + ["--augustus_species", "silkworm"], capture_output=True, text=True,
                          cwd=tmp_path)
    assert stop.returncode == 1
    assert "ERROR: augustus_species 'silkworm' maps to the Augustus model 'bombyx_mori'" in stop.stderr
    if gp.offline_ncbi_taxa() is not None:     # the taxon route needs ~/.etetoolkit/taxa.sqlite
        assert "Nearest installed models:" in stop.stderr
    else:                                      # no installed model is within difflib's 0.6 of bombyx_mori
        assert "No installed model is near it by taxon or by spelling." in stop.stderr
    if shutil.which("prodigal"):
        ok = subprocess.run(base + ["--augustus_species", "silkworm", "--gene_predictor", "prodigal"],
                            capture_output=True, text=True, cwd=tmp_path)
        assert ok.returncode == 0, ok.stderr
