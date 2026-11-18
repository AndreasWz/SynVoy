"""The conda environment of the steps (items 24 and 26 of the bug list of 16 Nov 2026).

Item 24: the packages of the language-model and structure search steps (environment-ml.yml) did not
reach the steps under the conda profiles. Every conda profile set process.conda to environment.yml,
and the header of environment-ml.yml said to install it into synvoy_env, an environment no step runs
in. Now environment-ml.yml holds environment.yml's packages (same versions) and the ML packages, and
one closure in nextflow.config gives the steps environment-ml.yml when --enable_plm_search or
--enable_structural_search is on; a process.conda of the user's own (-c) still wins.

Item 26: conda.useMamba was fixed to true, so a machine without mamba stopped the first run with
"mamba: command not found". It is true now only when an executable mamba is on the PATH of the
machine that runs Nextflow; otherwise conda builds the environment.

The file checks read the YAML line by line (CI installs no PyYAML); the Nextflow checks are skipped
without Nextflow on PATH.
"""
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
                 "params.enable_structural_search.toString().toBoolean()"):
        assert text in closure, text


def test_no_profile_fixes_use_mamba():
    code = _config_code()
    assert not re.search(r"conda\.useMamba\s*=\s*true", code)
    assert re.search(r"^\s*useMamba\s*=\s*\(System\.getenv\('PATH'\)", code, re.M)


def _stub_project(tmp_path):
    shutil.copy(os.path.join(ROOT, "nextflow.config"), tmp_path / "nextflow.config")
    shutil.copytree(os.path.join(ROOT, "conf"), tmp_path / "conf")     # the config includes its files
    (tmp_path / "main.nf").write_text('process A {\n    output: stdout\n    script:\n'
                                      '    "echo TASKCONDA=${task.conda}"\n}\nworkflow {\n    A() | view\n}\n')
    (tmp_path / "off.config").write_text("conda.enabled = false\n")   # no environment is built


def _task_conda(tmp_path, *args):
    proc = subprocess.run([NEXTFLOW, "run", "main.nf", "-ansi-log", "false", "-profile", "standard",
                           "-c", "off.config", *args], cwd=tmp_path, capture_output=True, text=True,
                          timeout=300)
    found = re.findall(r"TASKCONDA=(\S+)", proc.stdout)
    assert proc.returncode == 0 and found, proc.stdout[-2000:] + proc.stderr[-2000:]
    return os.path.basename(found[0]) if found[0].endswith(".yml") else found[0]


@needs_nextflow
def test_the_steps_get_the_ml_environment_when_an_ml_step_is_on(tmp_path):
    _stub_project(tmp_path)
    assert _task_conda(tmp_path) == "environment.yml"
    assert _task_conda(tmp_path, "--enable_plm_search", "true") == "environment-ml.yml"
    assert _task_conda(tmp_path, "--enable_structural_search", "true") == "environment-ml.yml"
    assert _task_conda(tmp_path, "--enable_structural_search", "false") == "environment.yml"
    (tmp_path / "own.config").write_text('process.conda = "/opt/own_env"\n')
    assert _task_conda(tmp_path, "-c", "own.config", "--enable_plm_search", "true") == "/opt/own_env"


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
