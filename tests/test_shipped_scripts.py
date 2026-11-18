"""The run hints install.sh prints at its end, and the launcher's own usage text.

install.sh ended with `nextflow run main.nf ... -profile standard` for the human gene LY6E (a
vertebrate gene, which needs about 32 GB of memory), past the launcher's checks and its
memory-safe default profile; its pro-mode command had no --home_species, and the placeholder file
name genome.fna names no species, so that run stopped at its start; and it pointed to a USAGE.md
that is not in the repository root. The hints now go through ./run_synvoy.sh: the check of the
installation, the quickstart's Easy Mode run and a Pro Mode command with --home_species, then the
docs to read. run_synvoy.sh's header showed the same two commands; it shows the README's now.
"""
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
needs_bash = pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")


def _read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return fh.read()


def _hints_source():
    """install.sh after its final check: the PASS/FAIL summary and what it prints after it."""
    text = _read("install.sh")
    return text[text.index('if [ "$FAIL" -eq 0 ]'):]


def _printed_hints():
    src = _hints_source()
    tail = src[src.index("\nfi\n") + len("\nfi\n"):]
    out = subprocess.run(["bash", "-c", "ENV_NAME=synvoy_env\n" + tail], capture_output=True,
                         text=True, check=True)
    return out.stdout


def _commands(text):
    """The launcher commands in a text, a line that ends in a backslash joined with the next."""
    joined = re.sub(r"\\\n[#\s]*", " ", text)
    return [re.sub(r"\s+", " ", line.strip(" #")) for line in joined.splitlines()
            if line.strip(" #").startswith("./run_synvoy.sh")]


def _options(command):
    return re.findall(r"(?<!\S)--([a-z][\w-]*)", command)


def _declared_params():
    config = _read("nextflow.config")
    code = "\n".join(re.sub(r"(?<![:\"'])//.*$", "", line) for line in config.splitlines())
    return set(re.findall(r"^\s*([A-Za-z_]\w*)\s*=(?!=)", code.split("profiles {")[0], re.M))


def test_install_hints_go_through_the_launcher():
    tail = _hints_source()
    assert "./run_synvoy.sh" in tail
    for old in ("nextflow run main.nf", "Q16553", "-profile standard", "conda activate"):
        assert old not in tail, old


@needs_bash
def test_install_sh_parses():
    subprocess.run(["bash", "-n", os.path.join(ROOT, "install.sh")], check=True)


@needs_bash
def test_every_doc_the_install_hints_name_exists():
    printed = _printed_hints()
    docs = re.findall(r"docs/[\w.-]+\.md", printed)
    assert {"docs/QUICKSTART.md", "docs/USAGE.md", "docs/PARAMETERS.md"} <= set(docs)
    for doc in docs:
        assert os.path.isfile(os.path.join(ROOT, doc)), doc
    assert "USAGE.md for the full parameter reference" not in printed


@needs_bash
def test_the_install_hints_are_commands_the_launcher_runs():
    commands = _commands(_printed_hints())
    assert commands[0] == "./run_synvoy.sh --check-only"
    assert len(commands) == 3, commands
    declared = _declared_params()
    for command in commands[1:]:
        unknown = [o for o in _options(command) if o not in declared]
        assert not unknown, (command, unknown)
    easy, pro = commands[1], commands[2]
    assert "--mode easy" in easy and "--query_id" in easy
    # the quickstart's own command, word for word
    assert easy in _read("docs/QUICKSTART.md")
    # a pro-mode run stops at its start without a home species that the file name cannot give
    assert "--mode pro" in pro
    for option in ("query", "home_species", "home_genome", "target_genomes", "outdir"):
        assert option in _options(pro), option


def _launcher_usage():
    text = _read("run_synvoy.sh")
    head = text.split("set -o pipefail")[0]
    return head[head.index("Usage"):]


def test_the_launcher_usage_shows_the_readme_commands():
    usage = _launcher_usage()
    assert "Q16553" not in usage
    commands = _commands(usage)
    assert "./run_synvoy.sh --check-only" in commands
    pros = [c for c in commands if "--mode pro" in c]
    assert pros and all("--home_species" in c for c in pros)
    declared = _declared_params()
    for command in commands:
        unknown = [o for o in _options(command) if o not in declared and o != "check-only"]
        assert not unknown, (command, unknown)
    readme = _read("README.md")
    easy = [c for c in commands if "--mode easy" in c and "-profile" not in c]
    assert easy and easy[0] in re.sub(r"\s+", " ", readme)
