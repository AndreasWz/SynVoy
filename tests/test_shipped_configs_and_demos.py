"""The configs, demos, scripts and the image the repository ships (bug list of 16 Nov 2026).

* item 39a: demo/run_uniprot.sh searched a genome that nothing downloaded (Anopheles gambiae), and
  demo/run_demo_insect.sh took every genome in the folder. The insect fetch script downloads every
  target of the two demos, the insect demo names its three targets, and the uniprot demo stops
  with a hint when one is missing.
* item 39c: the profiles test, test_melettin and test_tetramorium needed a test_data/ the
  repository does not hold and set a parameter SynVoy does not have (gene); they are gone.
* item 39d: scripts/slurm_submit.sh defaulted to hpc_singularity without a slot for the image.
* item 25: the image put its tools on the PATH only through its Docker entrypoint, which
  Singularity and Apptainer do not run.
"""
import os
import re
import shutil
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(rel):
    return open(os.path.join(ROOT, rel)).read()


def _code(text):
    return "\n".join(re.sub(r"(?<![:\"'])//.*$", "", line) for line in text.splitlines())


def _declared():
    text = _code(_read("nextflow.config"))
    start = text.index("params {")
    depth, j = 0, text.index("{", start)
    while True:
        depth += {"{": 1, "}": -1}.get(text[j], 0)
        if depth == 0:
            break
        j += 1
    return set(re.findall(r"^\s{4}([A-Za-z_]\w*)\s*=", text[start:j], re.M))


# ----------------------------------------------------------------------- 39c
def test_every_param_a_shipped_config_sets_exists():
    declared = _declared()
    configs = [os.path.join(dp, f) for dp, _d, fs in os.walk(os.path.join(ROOT, "conf")) for f in fs
               if f.endswith(".config")] + [os.path.join(ROOT, "demo", f) for f in os.listdir(os.path.join(ROOT, "demo"))
                                            if f.endswith(".config")]
    assert configs
    unknown = {}
    for path in sorted(configs):
        text = _code(open(path).read())
        names = set(re.findall(r"\bparams\.([A-Za-z_]\w*)\s*=(?!=)", text))
        for m in re.finditer(r"\bparams\s*\{(.*?)\n\}", text, re.S):
            names |= set(re.findall(r"^\s*([A-Za-z_]\w*)\s*=(?!=)", m.group(1), re.M))
        bad = sorted(names - declared)
        if bad:
            unknown[os.path.relpath(path, ROOT)] = bad
    assert not unknown, unknown


def test_every_included_config_exists_and_the_test_profiles_are_gone():
    text = _code(_read("nextflow.config"))
    for inc in re.findall(r"includeConfig\s+'([^']+)'", text):
        assert os.path.isfile(os.path.join(ROOT, inc)), inc
    for name in ("test", "test_melettin", "test_tetramorium"):
        assert not re.search(rf"^\s{{4}}{name}\s*\{{", text, re.M), name
        assert not os.path.exists(os.path.join(ROOT, "conf", f"{name}.config"))


# ----------------------------------------------------------------------- 39a
def _fetched_insect_names():
    return set(re.findall(r"\[GC[AF]_\d+\.\d+\]=(\w+)", _read("demo/fetch_data_insect.sh")))


def test_the_insect_data_holds_every_target_of_both_demos():
    fetched = _fetched_insect_names()
    uniprot = _read("demo/run_uniprot.sh")
    named = set(re.findall(r"\$DEMO/(?:genomes|target_gffs)/(\w+)\.(?:fna|gff)", uniprot))
    assert named == {"Drosophila_pseudoobscura", "Anopheles_gambiae"}
    assert named <= fetched
    insect = _read("demo/run_demo_insect.sh")
    assert "genomes/*.fna" not in insect and "target_gffs/*.gff" not in insect
    targets = set(re.search(r"for n in ([\w ]+); do", insect).group(1).split())
    assert targets == {"Drosophila_simulans", "Drosophila_yakuba", "Drosophila_pseudoobscura"}
    assert targets <= fetched


def test_the_demos_stop_when_a_target_is_missing(tmp_path):
    for script in ("demo/run_uniprot.sh", "demo/run_demo_insect.sh"):
        text = _read(script)
        assert 'missing. Run demo/fetch_data_insect.sh first."; exit 1; }' in text, script
        assert subprocess.run(["bash", "-n", os.path.join(ROOT, script)]).returncode == 0
    # the insect demo with its query in place but no target genomes: it stops before the launcher
    repo = tmp_path / "repo"
    (repo / "demo").mkdir(parents=True)
    shutil.copy(os.path.join(ROOT, "demo", "run_demo_insect.sh"), repo / "demo")
    (repo / "local_data" / "demo_insect" / "query").mkdir(parents=True)
    (repo / "local_data" / "demo_insect" / "query" / "Adh.faa").write_text(">q\nMKT\n")
    (repo / "run_synvoy.sh").write_text("#!/bin/sh\necho LAUNCHED\n")
    (repo / "run_synvoy.sh").chmod(0o755)
    proc = subprocess.run(["bash", str(repo / "demo" / "run_demo_insect.sh"), "Adh"], capture_output=True,
                          text=True, timeout=60, env={"PATH": os.environ["PATH"], "HOME": str(tmp_path)})
    assert proc.returncode == 1 and "LAUNCHED" not in proc.stdout
    assert "Drosophila_simulans.fna missing. Run demo/fetch_data_insect.sh first." in proc.stdout


# ----------------------------------------------------------------------- 39d
def test_slurm_submit_has_a_slot_for_the_image(tmp_path):
    text = _read("scripts/slurm_submit.sh")
    assert re.search(r'^CONTAINER=""', text, re.M)
    assert 'EXTRA_ARGS+=(--docker_container "${CONTAINER}")' in text
    assert '"${EXTRA_ARGS[@]}"' in text and "--mode easy" in text
    assert subprocess.run(["bash", "-n", os.path.join(ROOT, "scripts", "slurm_submit.sh")]).returncode == 0
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "nextflow").write_text("#!/bin/sh\necho NEXTFLOW_RAN\n")
    (stub / "nextflow").chmod(0o755)
    env = {"PATH": f"{stub}:/usr/bin:/bin", "HOME": str(tmp_path), "SCRATCH": str(tmp_path)}
    proc = subprocess.run(["bash", os.path.join(ROOT, "scripts", "slurm_submit.sh")], cwd=tmp_path,
                          capture_output=True, text=True, timeout=60, env=env)
    assert proc.returncode == 2 and "NEXTFLOW_RAN" not in proc.stdout
    assert "PROFILE=hpc_singularity needs CONTAINER=/path/to/synvoy.sif" in proc.stderr


# ----------------------------------------------------------------------- 25
def test_the_image_sets_its_tools_without_its_entrypoint():
    text = _read("Dockerfile")
    env = " ".join(l.strip() for l in text[text.index("ENV PATH="):].split("\n\n")[0].splitlines())
    path = re.search(r'PATH="([^"]+)"', env).group(1).split(":")
    assert path[:2] == ["/opt/conda/bin", "/opt/synvoy/bin"] and path[-1] == "$PATH"
    assert 'AUGUSTUS_CONFIG_PATH="/opt/conda/config/"' in env and 'CONDA_PREFIX="/opt/conda"' in env
    assert "micromamba install -y -n base" in text           # the base environment is /opt/conda
    install = _read("docs/INSTALL.md")
    assert "%environment" not in install and "docker-archive://synvoy.tar" in install
