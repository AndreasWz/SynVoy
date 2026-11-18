"""The conda environment of the steps (items 24 and 26 of the bug list of 16 Nov 2026).

Item 24: the packages of the language-model and structure search steps (environment-ml.yml) did not
reach the steps under the conda profiles. Every conda profile set process.conda to environment.yml,
and the header of environment-ml.yml said to install it into synvoy_env, an environment no step runs
in. Now environment-ml.yml holds environment.yml's packages (same versions) and the ML packages, and
one closure in nextflow.config gives the search (ITERATIVE_SEARCH, the only step that runs these
layers) environment-ml.yml when --enable_plm_search or --enable_structural_search is on; a
process.conda of the user's own (-c) still wins, for every step.

Review of the series (fix round 1): the closure first gave environment-ml.yml to every step, so
switching a layer on or off under -resume rebuilt every step's environment and reran every step; and
a comment edit of environment.yml changed its text, which names the environment (Nextflow hashes the
file's whole text), so every conda run would have built a new environment and rerun every step.
Every step but the search keeps environment.yml, and environment.yml keeps its text.

Item 26: conda.useMamba was fixed to true, so a machine without mamba stopped the first run with
"mamba: command not found". It is true now only when an executable mamba is on the PATH of the
machine that runs Nextflow; otherwise conda builds the environment.

The file checks read the YAML line by line (CI installs no PyYAML); the Nextflow checks are skipped
without Nextflow on PATH.
"""
import hashlib
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEXTFLOW = shutil.which("nextflow")
needs_nextflow = pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
ML_PACKAGES = {"foldseek", "pytorch", "transformers", "sentencepiece", "protobuf"}


def _deps(name):
    out = {}
    inside = False
    for line in open(os.path.join(ROOT, name)):
        code = line.split("#")[0].rstrip()
        if re.match(r"^dependencies:\s*$", code):
            inside = True
            continue
        if inside and re.match(r"^\S", code):
            break
        m = re.match(r"^\s*-\s+([A-Za-z0-9_.\-]+)\s*(.*)$", code)
        if inside and m:
            out[m.group(1)] = m.group(2).strip()
    return out


def _config_code():
    return "\n".join(re.sub(r"(?<![:\"'])//.*$", "", l) for l in open(os.path.join(ROOT, "nextflow.config")))


def test_the_ml_environment_holds_every_package_of_the_main_one_with_its_version():
    main, ml = _deps("environment.yml"), _deps("environment-ml.yml")
    assert len(main) >= 15 and "parasail-python" in main
    assert {k: ml.get(k) for k in main} == main
    assert ML_PACKAGES <= set(ml)
    assert not ML_PACKAGES & set(main)


def test_one_closure_chooses_the_environment():
    code = _config_code()
    profiles = code[code.index("\nprofiles {"):]
    assert not re.search(r"^\s*process\.conda\s*=", code, re.M)       # no profile fixes it
    m = re.search(r"^\s*conda\s*=\s*\{(.*)\}\s*$", profiles, re.M)
    assert m, "the process block's conda closure is missing"
    closure = m.group(1)
    for text in ("environment-ml.yml", "environment.yml", "params.enable_plm_search.toString().toBoolean()",
                 "params.enable_structural_search.toString().toBoolean()", "task.process == 'ITERATIVE_SEARCH'"):
        assert text in closure, text
    # A withName selector would outrank a process.conda of the user's own (-c) for the search.
    for body in _selector_bodies(code, "ITERATIVE_SEARCH"):
        assert not re.search(r"^\s*conda\s*=", body, re.M), \
            "the search's environment is chosen in the generic closure, not in a withName selector"


def _selector_bodies(code, name):
    """The text inside each withName: '<name>' { ... } block (braces counted)."""
    bodies = []
    for m in re.finditer(r"withName:\s*'%s'\s*\{" % re.escape(name), code):
        depth, i = 1, m.end()
        while depth and i < len(code):
            depth += {"{": 1, "}": -1}.get(code[i], 0)
            i += 1
        bodies.append(code[m.end():i - 1])
    assert bodies, name
    return bodies


# sha256 of environment.yml as of cb1425f (the pins of fix e6). Nextflow names a conda environment by a hash
# of the file's whole text, comments included, and that name is part of every task's hash: any edit of the
# file, even of a comment, builds a new environment on the next run and reruns every step under -resume.
# When the packages change on purpose, update this value in the same change and say so in the CHANGELOG.
ENVIRONMENT_YML_SHA256 = "72a623886e26cb7b0f94e59b5b9e1275c78294f12569beb614dde88d145f05c9"


def test_environment_yml_keeps_its_text():
    with open(os.path.join(ROOT, "environment.yml"), "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    assert digest == ENVIRONMENT_YML_SHA256, (
        "environment.yml changed: every conda run builds a new environment and reruns every step; "
        "if that is intended, update ENVIRONMENT_YML_SHA256")


def test_no_profile_fixes_use_mamba():
    code = _config_code()
    assert not re.search(r"conda\.useMamba\s*=\s*true", code)
    assert re.search(r"^\s*useMamba\s*=\s*\(System\.getenv\('PATH'\)", code, re.M)


def _stub_project(tmp_path):
    shutil.copy(os.path.join(ROOT, "nextflow.config"), tmp_path / "nextflow.config")
    shutil.copytree(os.path.join(ROOT, "conf"), tmp_path / "conf")     # the config includes its files
    # The search and one other step (any name but ITERATIVE_SEARCH), each printing its environment.
    (tmp_path / "main.nf").write_text(
        'process ITERATIVE_SEARCH {\n    output: stdout\n    script:\n'
        '    "echo TASKCONDA=${task.process}=${task.conda}"\n}\n'
        'process LOCATE_GENE {\n    output: stdout\n    script:\n'
        '    "echo TASKCONDA=${task.process}=${task.conda}"\n}\n'
        'workflow {\n    ITERATIVE_SEARCH() | view\n    LOCATE_GENE() | view\n}\n')
    (tmp_path / "off.config").write_text("conda.enabled = false\n")   # no environment is built


def _task_conda(tmp_path, *args):
    proc = subprocess.run([NEXTFLOW, "run", "main.nf", "-ansi-log", "false", "-profile", "standard",
                           "-c", "off.config", *args], cwd=tmp_path, capture_output=True, text=True,
                          timeout=300)
    found = dict(re.findall(r"TASKCONDA=(\w+)=(\S+)", proc.stdout))
    assert proc.returncode == 0 and set(found) == {"ITERATIVE_SEARCH", "LOCATE_GENE"}, \
        proc.stdout[-2000:] + proc.stderr[-2000:]
    return {step: os.path.basename(env) if env.endswith(".yml") else env for step, env in found.items()}


@needs_nextflow
def test_only_the_search_gets_the_ml_environment_when_an_ml_step_is_on(tmp_path):
    _stub_project(tmp_path)
    plain = {"ITERATIVE_SEARCH": "environment.yml", "LOCATE_GENE": "environment.yml"}
    ml = {"ITERATIVE_SEARCH": "environment-ml.yml", "LOCATE_GENE": "environment.yml"}
    assert _task_conda(tmp_path) == plain
    assert _task_conda(tmp_path, "--enable_plm_search", "true") == ml
    assert _task_conda(tmp_path, "--enable_structural_search", "true") == ml
    assert _task_conda(tmp_path, "--enable_structural_search", "false") == plain
    (tmp_path / "own.config").write_text('process.conda = "/opt/own_env"\n')
    assert _task_conda(tmp_path, "-c", "own.config", "--enable_plm_search", "true") == \
        {"ITERATIVE_SEARCH": "/opt/own_env", "LOCATE_GENE": "/opt/own_env"}


@needs_nextflow
def test_mamba_is_used_only_when_it_is_on_the_path(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    (tools / "mamba").write_text("#!/bin/sh\nexit 0\n")
    (tools / "mamba").chmod(0o755)
    base = os.pathsep.join([os.path.dirname(NEXTFLOW), os.path.dirname(shutil.which("java") or NEXTFLOW),
                            "/usr/bin", "/bin"])
    assert not any(os.access(os.path.join(d, "mamba"), os.X_OK) for d in base.split(os.pathsep)), \
        "a mamba next to nextflow or java: this check cannot tell the two cases apart here"

    def use_mamba(path):
        proc = subprocess.run([NEXTFLOW, "config", ROOT, "-profile", "hpc_conda", "-flat"], cwd=tmp_path,
                              capture_output=True, text=True, timeout=300, env=dict(os.environ, PATH=path))
        assert proc.returncode == 0, proc.stderr[-2000:]
        return re.search(r"^conda\.useMamba = (\w+)$", proc.stdout, re.M).group(1)
    assert use_mamba(base) == "false"
    assert use_mamba(str(tools) + os.pathsep + base) == "true"
