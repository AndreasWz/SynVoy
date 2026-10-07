"""Static wiring checks between nextflow.config, the .nf files and the bin/ CLIs.

Nothing here runs Nextflow. Each test pins one way the three layers have drifted
apart before:

* a module passing a flag its script does not accept (the task dies with
  "unrecognized arguments" -- only at run time, only on that code path);
* a ``params.x`` read that no config declares (it is ``null`` at run time);
* a boolean param tested by Groovy truthiness. Nextflow 26 hands every
  ``--name value`` over as a *String*, so ``--disable_x false`` is the truthy
  string "false": ``when: !params.disable_x`` then switches the step OFF, and the
  documented ``--require_paralog_panel false`` did nothing;
* a ``settings.<key>`` read that main.nf never puts into the settings map;
* a kill switch that exists in a script but cannot be reached from the pipeline.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NF_FILES = [ROOT / "main.nf"] + sorted((ROOT / "modules").glob("*.nf")) \
    + sorted((ROOT / "subworkflows").glob("*.nf"))
CONFIG = ROOT / "nextflow.config"


# ───────────────────────────── helpers ──────────────────────────────

def _strip_comment(line: str) -> str:
    """Drop a trailing // comment that is not part of a URL or a string."""
    m = re.search(r"(?<![:\"'])//", line)
    return line[:m.start()] if m else line


def _code(path: Path) -> list[str]:
    return [_strip_comment(l) for l in path.read_text().splitlines()]


def _declared_params() -> dict[str, str]:
    """name -> literal default, from the top-level ``params {}`` block."""
    out: dict[str, str] = {}
    depth, inside = 0, False
    for line in _code(CONFIG):
        if not inside:
            if re.match(r"^params\s*\{", line):
                inside, depth = True, line.count("{") - line.count("}")
            continue
        if depth == 1:
            m = re.match(r"^\s*([A-Za-z_]\w*)\s*=\s*(.*?)\s*$", line)
            if m:
                out.setdefault(m.group(1), m.group(2))
        depth += line.count("{") + line.count("[") - line.count("}") - line.count("]")
        if depth <= 0:
            break
    return out


def _params_set_anywhere() -> set[str]:
    """Params assigned in any config (profiles included), e.g. internal ``_host_cpus``."""
    names = set(_declared_params())
    configs = [CONFIG] + sorted((ROOT / "conf").rglob("*.config"))
    for cfg in configs:
        text = "\n".join(_code(cfg))
        names |= set(re.findall(r"\bparams\.([A-Za-z_]\w*)\s*=(?!=)", text))
        for m in re.finditer(r"\bparams\s*\{", text):
            depth, j = 0, m.end() - 1
            while j < len(text):
                depth += {"{": 1, "}": -1}.get(text[j], 0)
                if depth == 0:
                    break
                j += 1
            names |= set(re.findall(r"^\s*([A-Za-z_]\w*)\s*=(?!=)", text[m.end():j], re.M))
    return names


def _argparse_flags(script: Path) -> set[str]:
    flags = set()
    for node in ast.walk(ast.parse(script.read_text())):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "add_argument":
            flags |= {a.value for a in node.args
                      if isinstance(a, ast.Constant) and isinstance(a.value, str)
                      and a.value.startswith("--")}
    return flags


SCRIPTS = {p.name: _argparse_flags(p) for p in sorted((ROOT / "bin").glob("*.py"))}
FLAG = re.compile(r"(?<![\w-])(--[A-Za-z][A-Za-z0-9_-]*)")
VAR = re.compile(r"\$\{?([A-Za-z_]\w*)\}?")


def _invocations():
    """Yield (nf file, process, script name, flags) for every bin/ script call.

    A call is the line naming ``<script>.py`` plus its backslash-continued lines.
    Flags held in Groovy variables (``def x = cond ? "--flag v" : ""``) are
    resolved through the ``def`` lines of the same process.
    """
    for nf in NF_FILES:
        text = nf.read_text()
        for pm in re.finditer(r"^process\s+(\w+)\s*\{", text, re.M):
            depth, j = 0, pm.end() - 1
            while j < len(text):
                depth += {"{": 1, "}": -1}.get(text[j], 0)
                if depth == 0:
                    break
                j += 1
            body = text[pm.end():j].splitlines()
            defs: dict[str, set[str]] = {}
            for line in body:
                dm = re.match(r"\s*def\s+([A-Za-z_]\w*)\s*=\s*(.*)$", line)
                if dm:
                    defs[dm.group(1)] = set(FLAG.findall(dm.group(2)))
            i = 0
            while i < len(body):
                line = body[i]
                stripped = line.strip()
                sm = re.search(r"(?:^|[\s/])([a-z_0-9]+\.py)\b", line)
                if not sm or sm.group(1) not in SCRIPTS or stripped.startswith(("#", "//", "def ")):
                    i += 1
                    continue
                chunk = [line]
                while chunk[-1].rstrip().endswith("\\") and i + 1 < len(body):
                    i += 1
                    chunk.append(body[i])
                joined = "\n".join(l for l in chunk if not l.strip().startswith("#"))
                flags = set(FLAG.findall(joined))
                for var in VAR.findall(joined):
                    flags |= defs.get(var, set())
                yield nf.name, pm.group(1), sm.group(1), flags
                i += 1


INVOCATIONS = list(_invocations())


# ───────────────────────────── tests ──────────────────────────────

def test_the_parser_sees_the_pipeline():
    """Guard the guard: if the .nf layout changes and nothing is parsed any more,
    every test below would pass vacuously."""
    scripts_called = {s for _f, _p, s, _fl in INVOCATIONS}
    assert len(INVOCATIONS) >= 30, len(INVOCATIONS)
    assert {"iterative_search_runner.py", "cluster_grs.py", "plot_synteny.py",
            "generate_report.py", "resolve_effective_params.py"} <= scripts_called
    by_script = {s: fl for _f, _p, s, fl in INVOCATIONS}
    assert len(by_script["iterative_search_runner.py"]) > 80
    assert len(_declared_params()) > 150


@pytest.mark.parametrize("nf,process,script,flags", INVOCATIONS,
                         ids=[f"{p}:{s}" for _f, p, s, _fl in INVOCATIONS])
def test_modules_only_pass_flags_the_script_accepts(nf, process, script, flags):
    unknown = sorted(flags - SCRIPTS[script])
    assert not unknown, (
        f"{nf}:{process} passes {unknown} to {script}, which does not define them "
        f"(argparse would exit with 'unrecognized arguments')")


def test_every_param_read_is_declared():
    known = _params_set_anywhere()
    undeclared = {}
    for nf in NF_FILES + [CONFIG]:
        for n, line in enumerate(_code(nf), 1):
            for name in re.findall(r"\bparams\.([A-Za-z_]\w*)", line):
                if name not in known and name not in ("get", "put", "containsKey"):
                    undeclared.setdefault(name, f"{nf.name}:{n}")
    assert not undeclared, f"params read but declared nowhere (null at run time): {undeclared}"


def _boolean_params() -> set[str]:
    return {k for k, v in _declared_params().items() if v in ("true", "false")}


def test_boolean_params_are_never_tested_by_raw_truthiness():
    """`params.flag ? a : b`, `!params.flag`, `params.a || params.b` and a bare
    `params.flag` process input are all wrong for a CLI-supplied "false"."""
    bools = _boolean_params()
    assert len(bools) > 30
    offenders = []
    files = NF_FILES + [CONFIG] + sorted((ROOT / "conf").rglob("*.config"))
    for path in files:
        for n, line in enumerate(_code(path), 1):
            for m in re.finditer(r"\bparams\.(\w+)", line):
                name = m.group(1)
                if name not in bools:
                    continue
                pre, post = line[:m.start()], line[m.end():]
                if re.match(r"\s*=(?!=)", post):                      # config assignment
                    continue
                if pre.rstrip().endswith("paramBool(") \
                        or post.startswith(".toString().toBoolean()"):
                    continue
                if pre.rstrip().endswith("${") and post.lstrip().startswith("}"):
                    continue                                          # text for a str2bool CLI
                if path.name == "main.nf" and name == "prefer_large_genes":
                    continue      # handed to EXTRACT_FLANKING, which passes it as text
                offenders.append(f"{path.name}:{n}: {line.strip()[:110]}")
    assert not offenders, (
        "boolean params tested without coercion (wrap in paramBool(...) in main.nf or "
        ".toString().toBoolean() elsewhere):\n  " + "\n  ".join(offenders))


def test_settings_keys_read_by_modules_are_provided_by_main_nf():
    main = (ROOT / "main.nf").read_text()
    m = re.search(r"def preset_default_keys = \[(.*?)\]", main, re.S)
    assert m
    provided = set(re.findall(r"'([a-z_]+)'", m.group(1)))
    read = {}
    for nf in NF_FILES:
        for n, line in enumerate(_code(nf), 1):
            for key in re.findall(r"\bsettings\.([A-Za-z_]\w*)", line):
                read.setdefault(key, f"{nf.name}:{n}")
    missing = {k: v for k, v in read.items() if k not in provided}
    assert read, "no settings.<key> reads found -- parser out of date?"
    assert not missing, f"settings keys read but never provided: {missing}"


@pytest.mark.parametrize("param,module,flag,script", [
    ("disable_ambiguous_tier", "iterative_search.nf", "--disable_ambiguous_tier",
     "iterative_search_runner.py"),
    ("disable_distant_synteny_rescue", "cluster_regions.nf",
     "--disable_distant_synteny_rescue", "cluster_grs.py"),
    ("legacy_strand_score", "cluster_regions.nf", "--legacy_strand_score", "cluster_grs.py"),
    ("synteny_bridge_two_sided", "iterative_search.nf", "--synteny_bridge_two_sided",
     "iterative_search_runner.py"),
    ("phylo_placement_promote", "generate_report.nf", "--phylo_placement_promote",
     "generate_report.py"),
])
def test_documented_switches_reach_their_script(param, module, flag, script):
    """Each of these was once a switch a user could set that changed nothing,
    because the script had the flag and the module never passed it."""
    assert param in _declared_params(), f"{param} is not declared in nextflow.config"
    assert flag in SCRIPTS[script]
    text = (ROOT / "modules" / module).read_text()
    assert f"params.{param}" in text and flag in text, f"{module} does not wire {param}"


def test_plot_extra_args_passthrough_is_wired():
    assert _declared_params().get("plot_extra_args") in ("''", '""')
    text = (ROOT / "modules" / "plot_synteny.nf").read_text()
    assert "params.plot_extra_args" in text and "${plot_extra_args}" in text


def test_llm_key_is_not_exported_as_an_empty_env_var():
    """`env { LLM_API_KEY = ... ?: '' }` made Nextflow warn about an empty variable
    for every task (12,929 of 16,528 console lines in one cluster run)."""
    text = "\n".join(_code(CONFIG))
    m = re.search(r"^env\s*\{(.*?)^\}", text, re.S | re.M)
    assert m, "env {} block not found"
    assert "LLM_API_KEY" not in m.group(1)
    assert "PYTHONHASHSEED" in m.group(1)
    for scope in ("docker", "singularity", "apptainer"):
        assert re.search(rf"^{scope}\.envWhitelist\s*=.*LLM_API_KEY", text, re.M), scope
