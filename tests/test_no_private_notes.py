"""No user-visible text cites a private working note (item 29 of the bug list of 16 Nov 2026).

The repository is public; the notes it was developed with are not: docs/TODO.md,
docs/STATE_OF_THE_PROJECT.md, docs/HOW_IT_WORKS.md, CLAUDE.md, docs/notes/, the git-ignored tmp/,
a developer's home folder, the section marks of those notes (a paragraph sign followed by 1e or
A1) and their working-note ids (fix-phase items, review fixes, batches). A log line, a --help
text, a report note or a config comment that cites one points the user at nothing; a pointer names
a section of docs/ALGORITHM.md, docs/OUTPUT.md or docs/PARAMETERS.md instead.

Static checks over what a user reads or runs:
* every runtime string of bin/*.py (log, print and error messages, --help texts, report notes),
  and the module docstring a script shows as its --help description; other docstrings are not
  checked;
* the log, warning and error lines of main.nf, subworkflows/, modules/ and lib/, their //
  comments cut (code comments there are not checked; a bash comment inside a module's script is
  part of the task's hash, so a change there reruns the task under -resume);
* the whole text of nextflow.config, conf/, the demo configs and scripts, the launcher,
  install.sh, scripts/*.sh, the sentinel files, the CI workflow and the shipped docs;
* no default path into a developer's home folder in scripts/;
* every parameter-like name in a nextflow.config comment is a parameter or a removed one.

environment.yml is not checked: Nextflow names the conda environment by the file's whole text
(tests/test_conda_environments.py pins it), so its comments change only with its packages.
"""
import ast
import glob
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PRIVATE = re.compile(
    r"§ ?1[a-z]\b|§ ?1\.\d|§A\d|§1[5-9]\b|§2\d\b"            # section marks of the private notes
    r"|TODO\.md|TODO_JUN|STATE_OF_THE_PROJECT|HOW_IT_WORKS|CLAUDE\.md|docs/notes/"
    r"|CAN_BRIDGE_ANALYSIS|NEXT_SESSION_FAMILY|docs/analysis/|(?<![\w.])_local/"
    r"|(?<![\w./$~{}-])tmp/[\w-]"                             # the repository's git-ignored tmp/
    r"|(?<![\w.])/home/[A-Za-z_][\w.-]*/|LRZ job|\bjobs? \d{6,}"   # a home folder
    r"|fix[- ]phase|fix item \d|bug hunt|review fix|\bRF-\d|\bHG-\d|FINAL step \d"
    r"|batch [A-Z]+-?\d|batch of 20\d\d|freeze fix|label rework|preset plan|preset session"
    r"|control arm|\bQW\d")


def _rel(path):
    return os.path.relpath(path, ROOT)


def _files(*patterns):
    out = []
    for pattern in patterns:
        out += sorted(glob.glob(os.path.join(ROOT, pattern), recursive=True))
    return [p for p in out if os.path.isfile(p) and "__pycache__" not in p]


def _hits(text, where):
    return [f"{where}: {m.group(0)!r} in {text[max(0, m.start() - 50):m.end() + 50]!r}"
            for m in PRIVATE.finditer(text)]


def test_the_pattern_finds_what_it_is_for():
    for text in ("see docs/TODO.md", "§1x wrote 3 rows", "§A4 rank binning", "(fix-phase item 18)",
                 "(review fix RF-2)", "(batch G3A; off by default)", "FINAL step 5",
                 "STATE_OF_THE_PROJECT.md F9", "tmp/session_notes/x", "/home/someone/x",
                 "docs/HOW_IT_WORKS.md", "CLAUDE.md §17"):
        assert PRIVATE.search(text), text
    for text in ("docs/ALGORITHM.md 7.11", "docs/PARAMETERS.md § 20", "§ 7.8", "/tmp/x", "work/tmp/x",
                 "$TMPDIR/tmp/x", "the A1 module", "lrz_ai profile", "as in runs made before 2026-10-08",
                 "results/my_run", "~/.etetoolkit/taxa.sqlite"):
        assert not PRIVATE.search(text), text


# ----------------------------------------------------------------------------------- bin/*.py
def _docstring_nodes(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                out.add(body[0].value)
    return out


def _help_description(tree, text):
    """The part of the module docstring that the script prints as its --help description."""
    doc = ast.get_docstring(tree, clean=False)
    if not doc:
        return ""
    if re.search(r"description\s*=\s*__doc__\s*\.\s*split", text):
        return doc.split("\n\n")[0]
    if re.search(r"description\s*=\s*__doc__\b", text):
        return doc
    return ""


def test_bin_strings_cite_no_private_note():
    bad = []
    for path in _files("bin/*.py"):
        text = open(path, encoding="utf-8").read()
        tree = ast.parse(text)
        docstrings = _docstring_nodes(tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) and node not in docstrings:
                bad += _hits(node.value, f"{_rel(path)}:{node.lineno}")
        bad += _hits(_help_description(tree, text), f"{_rel(path)} (--help description)")
    assert not bad, "\n".join(bad)


def test_a_help_description_from_the_docstring_is_found():
    text = '"""First paragraph.\n\nSecond, docs/TODO.md."""\nap = ArgumentParser(description=__doc__)\n'
    assert "docs/TODO.md" in _help_description(ast.parse(text), text)
    text = text.replace("description=__doc__", "description=__doc__.split('\\n\\n')[0]")
    assert _help_description(ast.parse(text), text) == "First paragraph."


# ------------------------------------------------------------------------ Nextflow and Groovy
SINK = re.compile(r"\blog\.(?:info|warn|error|debug)\b|\berror\s*[(\"]|\bprintln\b|System\.err"
                  r"|\b(?:errors|warnings|validationErrors|validationWarnings)\s*<<|\bthrow\b")


def _code(line):
    return re.sub(r"(?<![:\"'])//.*$", "", line.rstrip("\n"))


def _message_lines(path):
    """(line number, code) of every log, warning or error line, with the lines that continue its
    message (a line that ends in '+')."""
    lines = [_code(line) for line in open(path, encoding="utf-8")]
    i = 0
    while i < len(lines):
        if not SINK.search(lines[i]):
            i += 1
            continue
        j = i
        while lines[j].rstrip().endswith("+") and j + 1 < len(lines):
            j += 1
        for k in range(i, j + 1):
            yield k + 1, lines[k]
        i = j + 1


def test_nextflow_messages_cite_no_private_note():
    bad = []
    for path in _files("main.nf", "subworkflows/*.nf", "modules/*.nf", "lib/*.groovy"):
        for number, code in _message_lines(path):
            bad += _hits(code, f"{_rel(path)}:{number}")
    assert not bad, "\n".join(bad)


def test_a_continued_message_is_read_to_its_end(tmp_path):
    nf = tmp_path / "x.nf"
    nf.write_text('    log.warn "first part " +\n        "second part, docs/TODO.md"\n'
                  '    // a comment, docs/TODO.md\n    def y = 1\n')
    assert [n for n, code in _message_lines(str(nf)) if PRIVATE.search(code)] == [2]


# ------------------------------------------------------------- configs, scripts, assets, docs
WHOLE_TEXT = ("nextflow.config", "conf/**/*.config", "demo/*.config", "demo/*.sh", "demo/README.md",
              "install.sh", "run_synvoy.sh", "scripts/*.sh", "assets/sentinels/*",
              ".github/workflows/*.yml", "README.md", "CHANGELOG.md", "docs/*.md",
              "scripts/benchmark/README.md")


def test_configs_scripts_assets_and_docs_cite_no_private_note():
    files = _files(*WHOLE_TEXT)
    assert len(files) > 20
    bad = []
    for path in files:
        for number, line in enumerate(open(path, encoding="utf-8"), 1):
            bad += _hits(line, f"{_rel(path)}:{number}")
    assert not bad, "\n".join(bad)


HOME_FOLDER = re.compile(r"(?:^|[\s\"'=:(])/home/[A-Za-z_][\w.-]*/")


def test_scripts_default_to_no_home_folder():
    bad = []
    for path in _files("scripts/**/*.py", "scripts/**/*.sh", "demo/*.sh", "install.sh", "run_synvoy.sh"):
        for number, line in enumerate(open(path, encoding="utf-8"), 1):
            if HOME_FOLDER.search(line):
                bad.append(f"{_rel(path)}:{number}: {line.strip()[:160]}")
    assert not bad, "\n".join(bad)


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_the_repository_default_of_the_benchmark_scripts_is_their_own_checkout(tmp_path):
    """REPO comes from the script's own place, so the scripts work in any checkout."""
    for name in ("run_mcscanx.sh", "run_toga_tier3.sh", "run_toga2_tier3.sh"):
        path = os.path.join(ROOT, "scripts", "benchmark", name)
        line = next(x for x in open(path, encoding="utf-8") if x.startswith("REPO="))
        probe = tmp_path / "scripts" / "benchmark" / name
        probe.parent.mkdir(parents=True, exist_ok=True)
        probe.write_text("#!/bin/bash\n" + line + 'echo "$REPO"\n')
        out = subprocess.run(["bash", str(probe)], capture_output=True, text=True, cwd="/",
                             env={k: v for k, v in os.environ.items() if k != "REPO"})
        assert out.stdout.strip() == str(tmp_path), (name, out.stdout, out.stderr)


# ------------------------------------------------------ parameter names in nextflow.config comments
PARAM_LIKE = re.compile(r"(?<![-\w./])([a-z][a-z0-9]*(?:_[a-z0-9]+)+)\b")
PARAM_WORD = re.compile(r"(?:^|_)(?:min|max|threshold|timeout|limit|margin|floor|frac|cap|gap)(?:_|$)"
                        r"|^(?:disable|enable|allow|require)_")


def _declared_and_removed():
    config = open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8").read()
    code = "\n".join(_code(line) for line in config.splitlines())
    names = set(re.findall(r"^\s*([A-Za-z_]\w*)\s*=(?!=)", code.split("profiles {")[0], re.M))
    names |= set(re.findall(r"\bparams\.([A-Za-z_]\w*)\s*=(?!=)", code))
    main = open(os.path.join(ROOT, "main.nf"), encoding="utf-8").read()
    block = re.search(r"^def removedParams\(\) \{(.*?)^\}", main, re.S | re.M).group(1)
    names |= set(re.findall(r"^\s+([a-z_]+)\s*:\s*\"", block, re.M))
    names |= set(re.findall(r"'([a-z_]+)'", block))
    return names


def test_parameters_named_in_config_comments_exist():
    known = _declared_and_removed()
    assert {"max_loci", "goi_hull_min_flanking", "pub_width_mm"} <= known
    bad = []
    text = open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8").read()
    for number, line in enumerate(text.splitlines(), 1):
        m = re.search(r"(?<![:\"'])//(.*)$", line)
        if not m:
            continue
        for name in PARAM_LIKE.findall(m.group(1)):
            if PARAM_WORD.search(name) and name not in known:
                bad.append(f"nextflow.config:{number}: {name}")
    assert not bad, "\n".join(bad)
