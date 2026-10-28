"""The Augustus model of a home genome without a GFF (profile stage 5c).

Without a GFF the home genome's genes are predicted (PREPARE_HOME_PROTEOME, and
extract_flanking_genes.py when the borrowed annotation holds none). Until 2026-10-27 that used
Augustus's fruit-fly model for every home species. With any --modules other than legacy, the
resolver now takes the installed model nearest the home species in the NCBI taxonomy on disk,
records it in effective_params.json (settings, source, derivation), and stops the run when a
eukaryotic home would fall back to Prodigal because Augustus is missing. A model the user names
wins; without taxonomy the model stays 'fly' with a warning; under legacy nothing changes.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import gene_predictor as gp  # noqa: E402
import resolve_effective_params as rep  # noqa: E402
from sequence_utils import offline_ncbi_taxa  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
NCBI = offline_ncbi_taxa()
MODELS = gp.installed_augustus_models()
real_taxonomy = pytest.mark.skipif(NCBI is None, reason="no NCBI taxonomy database on disk")
real_models = pytest.mark.skipif(not MODELS, reason="Augustus not installed")


class FakeTaxonomy:
    """root 1 > cellular 2 > Eukaryota 2759 > ... ; Bacteria 3 under cellular."""
    LINEAGES = {
        "Drosophila melanogaster": [1, 2, 2759, 50, 51, 7227],
        "Drosophila simulans": [1, 2, 2759, 50, 51, 7240],
        "Apis mellifera": [1, 2, 2759, 50, 60, 7460],
        "Apis dorsata": [1, 2, 2759, 50, 60, 7462],
        "Homo sapiens": [1, 2, 2759, 70, 9606],
        "Escherichia coli": [1, 2, 3, 562],
    }
    NAMES = {2759: ("Eukaryota", "superkingdom"), 50: ("Insecta", "class"), 60: ("Apis", "genus"),
             51: ("Drosophila", "genus"), 7227: ("Drosophila melanogaster", "species"),
             2: ("cellular organisms", "no rank")}

    def lineage(self, name):
        key = " ".join(str(name).replace("_", " ").split())
        lin = self.LINEAGES.get(key)
        return {"lineage": lin, "matched": key} if lin else None

    def describe(self, taxid):
        return self.NAMES.get(taxid, ("", ""))


FAKE_MODELS = {"fly": "Drosophila melanogaster", "honeybee1": "Apis mellifera",
               "adorsata": "Apis dorsata", "human": "Homo sapiens", "E_coli_K12": "Escherichia coli",
               "nowhere": "Nonexistus nullus"}


def test_the_model_of_the_species_itself_and_the_deepest_shared_taxon():
    tax = FakeTaxonomy()
    model, info = gp.nearest_augustus_model("Drosophila_melanogaster", FAKE_MODELS, tax)
    assert model == "fly" and info["same_species"] and info["eukaryote"] is True
    model, info = gp.nearest_augustus_model("Drosophila simulans", FAKE_MODELS, tax)
    assert model == "fly" and info["shared_taxon"] == "Drosophila" and info["shared_rank"] == "genus"
    assert not info["same_species"]
    # a tie (two Apis models for a third Apis species) goes to the first model id
    tax.LINEAGES["Apis cerana"] = [1, 2, 2759, 50, 60, 7461]
    model, _ = gp.nearest_augustus_model("Apis cerana", FAKE_MODELS, tax)
    assert model == "adorsata"
    model, info = gp.nearest_augustus_model("Escherichia coli", FAKE_MODELS, tax)
    assert model == "E_coli_K12" and info["eukaryote"] is False


def test_no_lineage_gives_no_model():
    model, info = gp.nearest_augustus_model("Notarealus specius", FAKE_MODELS, FakeTaxonomy())
    assert model is None and info["eukaryote"] is None and "not in the taxonomy" in info["reason"]
    model, info = gp.nearest_augustus_model("Homo sapiens", FAKE_MODELS, None)
    assert model is None


def test_installed_models_from_the_help_list_and_the_config_folder(tmp_path):
    help_text = ("identifier | species\n-----|-----\n"
                 "fly                    | Drosophila melanogaster\n"
                 "(fly_exp)              | Drosophila melanogaster\n"
                 "honeybee1              | Apis mellifera\n"
                 "chlorella              | Chlorella sp.\n")
    species = tmp_path / "species"
    for name in ("fly", "honeybee1", "Anopheles_gambiae", "generic", "fly_exp"):
        (species / name).mkdir(parents=True)
        (species / name / f"{name}_parameters.cfg").write_text("x")
    models = gp.installed_augustus_models(help_text, str(tmp_path))
    # chlorella is listed but not installed; Anopheles_gambiae is installed but not listed
    assert models == {"fly": "Drosophila melanogaster", "honeybee1": "Apis mellifera",
                      "Anopheles_gambiae": "Anopheles gambiae"}


@real_taxonomy
@real_models
@pytest.mark.parametrize("species,model", [
    ("Drosophila_melanogaster", "fly"),
    ("Apis mellifera", "honeybee1"),
    ("Homo sapiens", "human"),
    ("Candida albicans", "candida_albicans"),
    ("Candida albicans SC5314", "candida_albicans"),
    ("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)", "saccharomyces_cerevisiae_S288C"),
])
def test_real_lineages(species, model):
    if model not in MODELS:
        pytest.skip(f"{model} not installed")
    got, info = gp.nearest_augustus_model(species, MODELS, gp.TaxonomyLookup(NCBI))
    assert got == model and info["eukaryote"] is True


@real_taxonomy
def test_an_unknown_species_keeps_fly_with_a_warning():
    out = rep.home_gene_model("Notarealus specius", "auto", "auto", models=MODELS or FAKE_MODELS,
                              taxonomy=gp.TaxonomyLookup(NCBI), augustus_present=True)
    assert (out["value"], out["source"]) == ("fly", "default")
    assert out["warning"] and "Notarealus specius" in out["warning"] and out["stop"] is None


def test_the_rule_of_the_resolver():
    tax = FakeTaxonomy()
    out = rep.home_gene_model("Apis_mellifera", "auto", "auto", models=FAKE_MODELS, taxonomy=tax,
                              augustus_present=True)
    assert (out["value"], out["source"], out["warning"], out["stop"]) == (
        "honeybee1", "home_lineage", None, None)
    assert out["derivation"]["rule"] == "home_lineage" and out["derivation"]["default"] == "fly"
    assert "nearest the home species Apis mellifera: the same species" in out["derivation"]["line"]
    # a model the user named wins
    out = rep.home_gene_model("Apis_mellifera", "fly", "auto", models=FAKE_MODELS, taxonomy=tax,
                              augustus_present=True)
    assert (out["value"], out["source"]) == ("fly", "user_set")
    # without the database: 'fly' and a warning
    out = rep.home_gene_model("Apis_mellifera", "auto", "auto", models=FAKE_MODELS, taxonomy=None,
                              augustus_present=True)
    assert (out["value"], out["source"]) == ("fly", "default") and "taxonomy" in out["warning"]


def test_a_eukaryote_without_augustus_stops_auto():
    tax = FakeTaxonomy()
    out = rep.home_gene_model("Homo sapiens", "auto", "auto", models={}, taxonomy=tax,
                              augustus_present=False)
    assert out["stop"] and "Prodigal cannot predict genes with introns" in out["stop"]
    for species, predictor in [("Homo sapiens", "prodigal"), ("Escherichia coli", "auto"),
                               ("Notarealus specius", "auto")]:
        out = rep.home_gene_model(species, "auto", predictor, models={}, taxonomy=tax,
                                  augustus_present=False)
        assert out["stop"] is None, (species, predictor)
    out = rep.home_gene_model("Homo sapiens", "auto", "auto", models=FAKE_MODELS, taxonomy=tax,
                              augustus_present=True)
    assert out["stop"] is None


def _resolver(tmp_path, *extra, env=None):
    defaults = tmp_path / "defaults.json"
    defaults.write_text(json.dumps({"goi_margin_genes": 2}))
    out = tmp_path / f"ep_{len(list(tmp_path.glob('ep_*')))}.json"
    proc = subprocess.run([sys.executable, str(ROOT / "bin" / "resolve_effective_params.py"),
                           "--defaults", str(defaults), "--output", str(out), *extra],
                          capture_output=True, text=True, env=env)
    return proc, (json.loads(out.read_text()) if out.exists() else None)


def test_legacy_and_a_home_gff_record_nothing(tmp_path):
    for extra in (("--modules", "legacy", "--home_gff_present", "false", "--home_species", "Homo sapiens"),
                  ("--modules", "legacy"),
                  ("--modules", "auto", "--home_gff_present", "true", "--home_species", "Homo sapiens"),
                  ("--modules", "auto")):
        proc, ep = _resolver(tmp_path, *extra)
        assert proc.returncode == 0, proc.stderr
        assert "augustus_species" not in ep["settings"] and "augustus_species" not in ep["derivation"]
        assert "HOME-GENES" not in proc.stdout


@real_taxonomy
@real_models
def test_the_resolver_records_the_choice(tmp_path):
    proc, ep = _resolver(tmp_path, "--modules", "auto", "--home_gff_present", "false",
                         "--home_species", "Drosophila_melanogaster")
    assert proc.returncode == 0, proc.stderr
    assert ep["settings"]["augustus_species"] == "fly"
    assert ep["source"]["augustus_species"] == "home_lineage"
    assert ep["derivation"]["augustus_species"]["shared_rank"] == "species"
    assert "HOME-GENES: augustus_species fly" in proc.stdout


@real_taxonomy
def test_the_resolver_stops_a_eukaryote_without_augustus(tmp_path):
    env = dict(os.environ, PATH="/usr/bin:/bin")
    if shutil.which("augustus", path=env["PATH"]):
        pytest.skip("augustus on the system PATH")
    proc, ep = _resolver(tmp_path, "--modules", "auto", "--home_gff_present", "false",
                         "--home_species", "Homo sapiens", env=env)
    assert proc.returncode == 2 and "needs Augustus" in proc.stderr
    assert ep["settings"]["augustus_species"] == "fly"            # the JSON is written first


def test_the_modules_read_the_model_through_the_helper():
    for nf in ("prepare_home.nf", "extract_flanking.nf", "iterative_search.nf"):
        text = (ROOT / "modules" / nf).read_text()
        assert "params.augustus_species" not in text, nf
        assert "HomeGenes.augustusSpecies(params, settings)" in text, nf
    assert "HomeGenes.resolverArgs(params, home_species, home_gff_present)" in (
        ROOT / "modules" / "resolve_effective_params.nf").read_text()
    config = (ROOT / "nextflow.config").read_text()
    assert "    augustus_species = 'auto'\n" in config


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_the_groovy_helper(tmp_path):
    """Under legacy 'auto' is 'fly' and the resolver gets no argument, as before 5c."""
    (tmp_path / "lib").mkdir()
    for f in ("GenomeAxis.groovy", "HomeGenes.groovy"):
        shutil.copy(ROOT / "lib" / f, tmp_path / "lib" / f)
    (tmp_path / "t.nf").write_text(
        "workflow {\n"
        "    [['legacy', 'auto', [:]], ['legacy', 'human', [:]], ['auto', 'auto', [:]],\n"
        "     ['auto', 'auto', [augustus_species: 'candida_albicans']], ['auto', 'human', [augustus_species: 'human']]].each { c ->\n"
        "        def p = [modules: c[0], augustus_species: c[1], gene_predictor: 'auto']\n"
        "        println \"RESULT ${HomeGenes.augustusSpecies(p, c[2])}|${HomeGenes.resolverArgs(p, \"Candida albicans (strain O'Hara)\", false)}|\"\n"
        "    }\n}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    args = (" --home_gff_present false --home_species 'Candida albicans (strain O'\"'\"'Hara)'"
            " --augustus_species '{}' --gene_predictor auto")
    assert rows == ["fly||", "human||", "fly|" + args.format("auto") + "|",
                    "candida_albicans|" + args.format("auto") + "|",
                    "human|" + args.format("human") + "|"], proc.stdout + proc.stderr
