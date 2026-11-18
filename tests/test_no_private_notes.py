"""No user-visible text cites a private working note.

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

Only files of the repository are read: in a git checkout the files git tracks (git ls-files), and
in a plain export, or a copy that lies inside another checkout's ignored folder, the files no
.gitignore ignores. The local notes a developer keeps git-ignored in the checkout (docs/TODO.md and
the like) are what the checks are about, so they never fail them.

environment.yml is not checked: Nextflow names the conda environment by the file's whole text
(tests/test_conda_environments.py pins it), so its comments change only with its packages.
"""
import ast
import functools
import glob
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# A path into a home folder: /home/<name>/ at the start, or after a space, a quote, '=', ':' or '('
# (not $DEMO/home/... or local_data/demo/home/..., relative paths with a folder named home).
HOME_PATH = r"(?:^|[\s\"'=:(])/home/[A-Za-z_][\w.-]*/"

PRIVATE = re.compile(
    r"§ ?1[a-z]\b|§ ?1\.\d|§A\d|§1[5-9]\b|§2\d\b"            # section marks of the private notes
    r"|TODO\.md|TODO_JUN|STATE_OF_THE_PROJECT|HOW_IT_WORKS|CLAUDE\.md|docs/notes/"
    r"|CAN_BRIDGE_ANALYSIS|NEXT_SESSION_FAMILY|docs/analysis/|(?<![\w.])_local/"
    r"|(?<![\w./$~{}-])tmp/[\w-]"                             # the repository's git-ignored tmp/
    r"|" + HOME_PATH + r"|LRZ job|\bjobs? \d{6,}"
    r"|fix[- ]phase|fix item \d|bug hunt|review fix|\bRF-\d|\bHG-\d|FINAL step \d"
    r"|batch [A-Z]+-?\d|batch of 20\d\d|freeze fix|label rework|preset plan|preset session"
    r"|control arm|\bQW\d")


def _rel(path, root=ROOT):
    return os.path.relpath(path, root).replace(os.sep, "/")


def _git_env():
    """The environment of a git call without GIT_* variables: a git hook that runs the tests sets
    GIT_DIR or GIT_INDEX_FILE, which would point the call at another repository."""
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


@functools.lru_cache(maxsize=None)
def _tracked_files(root):
    """The files git tracks in ``root`` ('/'-separated, relative to it), or None when ``root`` is
    not the top of a git checkout: a plain export, or a copy inside another checkout's ignored
    folder, where git would answer for the files of that other checkout. Read-only git only."""
    git = shutil.which("git")
    if git is None:
        return None
    try:
        top = subprocess.run([git, "-C", root, "rev-parse", "--show-toplevel"], capture_output=True,
                             text=True, env=_git_env(), timeout=60)
        if top.returncode != 0 or os.path.realpath(top.stdout.strip()) != os.path.realpath(root):
            return None
        listed = subprocess.run([git, "-C", root, "ls-files", "-z"], capture_output=True,
                                env=_git_env(), timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    if listed.returncode != 0:
        return None
    return frozenset(p for p in listed.stdout.decode("utf-8", "surrogateescape").split("\0") if p)


def _glob_regex(pattern):
    """A .gitignore pattern as a regular expression over a '/'-separated path: '*' and '?' stop
    at '/', '**/' is any number of folders, a trailing '/**' everything inside, '[...]' a class
    ('!' negates it), a backslash quotes the next character."""
    out, i = [], 0
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        elif pattern[i] == "[" and "]" in pattern[i + 2:]:
            j = pattern.index("]", i + 2)
            body = pattern[i + 1:j]
            out.append("[" + ("^" + body[1:] if body.startswith("!") else body) + "]")
            i = j + 1
        elif pattern[i] == "\\" and i + 1 < len(pattern):
            out.append(re.escape(pattern[i + 1]))
            i += 2
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out))


class _GitIgnore:
    """Whether a path ('/'-separated, relative to ``root``) is ignored by the .gitignore files of
    ``root`` and its folders, read as git reads them: blank and '#' lines skipped; '!' includes
    again; a trailing '/' matches folders only; a pattern with a '/' before its end is relative to
    the folder of its .gitignore, one without matches a name at any depth; a later rule wins, a
    deeper file's rules come after; a path inside an ignored folder is ignored."""

    def __init__(self, root):
        self.root = root
        self._rules = {}

    def rules(self, folder):
        if folder not in self._rules:
            rules = []
            path = os.path.join(self.root, folder, ".gitignore")
            if os.path.isfile(path):
                with open(path, encoding="utf-8") as fh:
                    for line in fh:
                        line = line.rstrip()
                        if not line or line.startswith("#"):
                            continue
                        negate = line.startswith("!")
                        line = line[1:] if negate else line
                        folder_only = line.endswith("/")
                        line = line.rstrip("/")
                        if line:
                            rules.append((_glob_regex(line.lstrip("/")), negate, folder_only, "/" in line))
            self._rules[folder] = rules
        return self._rules[folder]

    def _excluded(self, parts, is_folder):
        verdict = False
        for depth in range(len(parts)):
            below = "/".join(parts[depth:])
            for regex, negate, folder_only, anchored in self.rules("/".join(parts[:depth])):
                if folder_only and not is_folder:
                    continue
                if regex.fullmatch(below if anchored else parts[-1]):
                    verdict = not negate
        return verdict

    def __call__(self, rel):
        parts = rel.split("/")
        return any(self._excluded(parts[:i], i < len(parts)) for i in range(1, len(parts) + 1))


def _files(*patterns, root=ROOT):
    """The files of the repository the glob patterns match under ``root``: those git tracks when
    ``root`` is the top of a git checkout, else those no .gitignore ignores."""
    out = []
    for pattern in patterns:
        out += sorted(glob.glob(os.path.join(root, pattern), recursive=True))
    out = [p for p in out if os.path.isfile(p) and "__pycache__" not in p]
    tracked = _tracked_files(root)
    if tracked is not None:
        return [p for p in out if _rel(p, root) in tracked]
    ignored = _GitIgnore(root)
    return [p for p in out if not ignored(_rel(p, root))]


def _hits(text, where):
    return [f"{where}: {m.group(0)!r} in {text[max(0, m.start() - 50):m.end() + 50]!r}"
            for m in PRIVATE.finditer(text)]


def test_the_pattern_finds_what_it_is_for():
    for text in ("see docs/TODO.md", "§1x wrote 3 rows", "§A4 rank binning", "(fix-phase item 18)",
                 "(review fix RF-2)", "(batch G3A; off by default)", "FINAL step 5",
                 "STATE_OF_THE_PROJECT.md F9", "tmp/session_x/y", "/home/someone/x",
                 "PATH=/home/someone/bin", "docs/HOW_IT_WORKS.md", "CLAUDE.md §17"):
        assert PRIVATE.search(text), text
    for text in ("docs/ALGORITHM.md 7.11", "docs/PARAMETERS.md § 20", "§ 7.8", "/tmp/x", "work/tmp/x",
                 "$TMPDIR/tmp/x", "the A1 module", "lrz_ai profile", "as in runs made before 2026-10-08",
                 "results/my_run", "~/.etetoolkit/taxa.sqlite", '"$DEMO/home/genome.fna"',
                 "local_data/demo/home/genome.gff", "${GENOMES_DIR}/home/genome.gff"):
        assert not PRIVATE.search(text), text


# ------------------------------------------------------------------- which files are read
def _write(path, text="x\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _git(repo, *args):
    """git init or git add in a test repository under tmp_path, without the user's or the
    system's git configuration."""
    env = dict(_git_env(), HOME=str(repo), GIT_CONFIG_NOSYSTEM="1")
    subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, env=env)


needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="needs git")


@needs_git
def test_in_a_git_checkout_only_the_tracked_files_are_read(tmp_path):
    """git ls-files lists the index, so no commit is needed: an ignored local note, and a new file
    nobody added, are not read; a file added in spite of .gitignore is."""
    repo = tmp_path / "repo"
    for rel in ("docs/A.md", "docs/TODO.md", "docs/new.md", "docs/forced.md"):
        _write(repo / rel, "see docs/TODO.md\n")
    _write(repo / ".gitignore", "docs/TODO.md\ndocs/forced.md\n")
    _git(repo, "init", "-q")
    _git(repo, "add", ".gitignore", "docs/A.md")
    _git(repo, "add", "-f", "docs/forced.md")
    assert _tracked_files(str(repo)) == {".gitignore", "docs/A.md", "docs/forced.md"}
    assert [_rel(p, str(repo)) for p in _files("docs/*.md", root=str(repo))] == ["docs/A.md",
                                                                                 "docs/forced.md"]


@needs_git
def test_a_copy_inside_another_checkouts_ignored_folder_reads_its_own_gitignore(tmp_path):
    """git answers there for the outer checkout, whose index holds none of the copy's files, so
    the copy is read as a plain export: every file but those its .gitignore ignores."""
    outer = tmp_path / "outer"
    _write(outer / ".gitignore", "tmp/\n")
    _git(outer, "init", "-q")
    _git(outer, "add", ".gitignore")
    copy = outer / "tmp" / "copy"
    _write(copy / ".gitignore", "docs/TODO.md\n")
    for rel in ("docs/A.md", "docs/B.md", "docs/TODO.md"):
        _write(copy / rel)
    assert _tracked_files(str(copy)) is None
    assert [_rel(p, str(copy)) for p in _files("docs/*.md", root=str(copy))] == ["docs/A.md",
                                                                                 "docs/B.md"]


def test_a_plain_export_skips_what_its_gitignore_files_ignore(tmp_path):
    export = tmp_path / "export"
    _write(export / ".gitignore", "# local material\nwork/\n*.log\n!keep.log\ndocs/TODO.md\n"
                                  "docs/notes/\nresults_*/\n_local/\nscripts/lrz/\na/**/deep.md\n")
    _write(export / "web" / ".gitignore", "var/\n")
    ignored = _GitIgnore(str(export))
    for rel in ("docs/TODO.md", "docs/notes/x.md", "work/x.md", "a/work/x.md", "x.log", "b/x.log",
                "work/keep.log", "results_2/x.md", "_local/x.md", "scripts/lrz/x.sh", "a/deep.md",
                "a/b/c/deep.md", "web/var/x.json", "web/a/var/x.json"):
        assert ignored(rel), rel
    for rel in ("docs/ALGORITHM.md", "docs/x/TODO.md", "keep.log", "work", "results/x.md",
                "scripts/x.sh", "b/deep.md", "var/x.json", "web/x.json"):
        assert not ignored(rel), rel
    for rel in ("docs/A.md", "docs/TODO.md", "docs/notes/n.md", "scripts/lrz/s.sh", "scripts/r.sh"):
        _write(export / rel)
    assert [_rel(p, str(export)) for p in _files("docs/**/*.md", "scripts/**/*.sh", root=str(export))] == [
        "docs/A.md", "scripts/r.sh"]


def test_the_repositorys_gitignore_keeps_the_local_notes_out():
    """What the fallback reads in an export of this repository: the local notes are ignored, the
    shipped docs are not."""
    ignored = _GitIgnore(ROOT)
    for rel in ("docs/TODO.md", "docs/HOW_IT_WORKS.md", "docs/STATE_OF_THE_PROJECT.md", "CLAUDE.md",
                "docs/notes/x.md", "_local/x.md", "tmp/x.md", "scripts/lrz/x.sh"):
        assert ignored(rel), rel
    for rel in ("docs/ALGORITHM.md", "docs/OUTPUT.md", "README.md", "scripts/replot_run.sh",
                "nextflow.config"):
        assert not ignored(rel), rel


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


HOME_FOLDER = re.compile(HOME_PATH)


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
