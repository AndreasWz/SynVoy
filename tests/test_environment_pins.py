"""Fix e6 (9 Nov 2026): environment.yml pins the tools to the versions the benchmark's container ran.

Every run records the versions it used (effective_params.json: tool_versions, run_resources.nextflow_version; fix-phase
d8). Until 2026-11-09 environment.yml gave lower bounds only, so two machines could run other versions of the same
code. r3 (40 LRZ jobs on fce4ce3, 256 effective_params.json) recorded one set, R3 below; each recorded tool is pinned
exactly. The file is read line by line (CI installs no PyYAML).
"""
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

R3 = {"augustus": "3.5.0", "blastp": "2.17.0+", "datasets": "18.30.1", "iqtree": "3.1.2", "mafft": "7.526",
      "miniprot": "0.18-r281", "mmseqs": "18.8cc5c", "parasail": "1.3.4", "prodigal": "2.6.3", "python": "3.12.13",
      "samtools": "1.23", "tblastn": "2.17.0+"}
R3_NEXTFLOW = "26.04.3"
PACKAGE = {"augustus": "augustus", "blastp": "blast", "tblastn": "blast", "datasets": "ncbi-datasets-cli",
           "iqtree": "iqtree", "mafft": "mafft", "miniprot": "miniprot", "mmseqs": "mmseqs2",
           "parasail": "parasail-python", "prodigal": "prodigal", "python": "python", "samtools": "samtools"}


def _specs():
    out = {}
    for line in open(os.path.join(ROOT, "environment.yml")):
        m = re.match(r"^\s*-\s+([A-Za-z0-9_.\-]+)\s*(\S*)", line.split("#")[0])
        if m:
            out[m.group(1)] = m.group(2)
    return out


def test_every_recorded_tool_is_pinned_to_r3s_version():
    specs = _specs()
    for tool, version in R3.items():
        want = re.sub(r"(\+|-r\d+)$", "", version)          # blastp's "2.17.0+", miniprot's "0.18-r281"
        assert specs.get(PACKAGE[tool]) == f"=={want}", (tool, specs.get(PACKAGE[tool]))
    assert specs.get("nextflow") == f"=={R3_NEXTFLOW}"


def test_the_record_and_the_pins_name_the_same_tools():
    """A tool the run records but environment.yml does not pin would let its version drift."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("rep", os.path.join(ROOT, "bin", "resolve_effective_params.py"))
    rep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rep)
    recorded = {name for name, *_rest in rep.TOOL_VERSION_COMMANDS} | {"python", "parasail"}
    assert recorded - {"iqtree2"} == set(R3) - {"iqtree2"}
