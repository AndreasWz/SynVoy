"""The gene tree is reproducible by default.

IQ-TREE 3.1.2 gives the same tree on a rerun of the same alignment only in one thread; with more,
branch lengths and supports can differ. compute_tree_cpus, which sets the tree step's threads, was 2
by default (and every CPU under docker_max), so a default run that inferred a tree could give another
one on a rerun. Now the default is 1, docker_max follows the parameter, the tree's info records
whether a rerun gives the same tree (reproducible), compute_tree.py warns above one thread, and the
report warns for such a tree.
"""
import json
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)


def _config():
    return open(os.path.join(ROOT, "nextflow.config")).read()


def _profile(text, name):
    m = re.search(rf"^    {name} \{{\n(.*?)^    \}}\n", text, re.S | re.M)
    assert m, name
    return m.group(1)


def test_the_default_is_one_thread():
    assert re.search(r"^    compute_tree_cpus = 1\s*$", _config(), re.M)


def test_docker_max_follows_the_parameter():
    body = _profile(_config(), "docker_max")
    assert "params.compute_tree_cpus =" not in body
    block = re.search(r"withName: 'COMPUTE_TREE' \{(.*?)\}", body, re.S).group(1)
    assert "cpus = params.compute_tree_cpus as Integer" in block
    top = re.search(r"withName: 'COMPUTE_TREE' \{(.*?)\}", _config().split("\nprofiles {")[1].split("\n}\n", 1)[1],
                    re.S).group(1)
    assert "cpus = params.compute_tree_cpus as Integer" in top


def _tree_main(tmp_path, monkeypatch, threads):
    import compute_tree as ct
    faa = tmp_path / "goi.faa"
    faa.write_text("".join(f">GOI_x|s{i}\nMKTAYIAKQR{'A' * i}\n" for i in range(5)))
    seen = []

    def run(cmd, *a, **k):
        seen.append(cmd)
        if "--prefix" in cmd:
            open(cmd[cmd.index("--prefix") + 1] + ".treefile", "w").write("(a,b,c);\n")
        else:
            k["stdout"].write(">a\nM\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(ct.subprocess, "run", run)
    monkeypatch.setattr(sys, "argv", ["compute_tree.py", "--input", str(faa), "--output",
                                      str(tmp_path / "locus_1_tree.nwk"), "--threads", str(threads)])
    ct.main()
    iq = next(c for c in seen if "--prefix" in c)
    return json.loads((tmp_path / "locus_1_tree.info.json").read_text()), iq[iq.index("-T") + 1]


def test_one_thread_is_recorded_as_reproducible(tmp_path, monkeypatch, capsys):
    info, t = _tree_main(tmp_path, monkeypatch, 1)
    assert info["threads"] == 1 and info["reproducible"] is True and t == "1"
    assert "WARNING" not in capsys.readouterr().err


def test_two_threads_are_recorded_and_warned(tmp_path, monkeypatch, capsys):
    info, t = _tree_main(tmp_path, monkeypatch, 2)
    assert info["threads"] == 2 and info["reproducible"] is False and t == "2"
    assert "IQ-TREE runs with 2 threads: a rerun on the same alignment can give" in capsys.readouterr().err


def _report(tmp_path, trees):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    (tmp_path / "ti").mkdir()
    for locus, info in trees.items():
        (tmp_path / "ti" / f"{locus}_tree.info.json").write_text(json.dumps(info))
    p = subprocess.run([sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
                        "--output", str(tmp_path / "r.json"), "--allow-empty",
                        "--tree_info_dir", str(tmp_path / "ti")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return json.loads((tmp_path / "r.json").read_text())


def test_the_report_warns_for_a_tree_of_several_threads(tmp_path):
    boot = {"status": "ultrafast", "replicates": 1000, "reason": ""}
    rep = _report(tmp_path, {"locus_1": {"locus": "locus_1", "seed": 1, "threads": 2, "bootstrap": boot},
                             "locus_2": {"locus": "locus_2", "seed": 1, "threads": 1, "bootstrap": boot},
                             "locus_3": {"locus": "locus_3", "seed": 1, "threads": 4,
                                         "bootstrap": dict(boot, status="placeholder_tree")}})
    warned = [w for w in rep["warnings"] if "IQ-TREE threads" in w]
    assert warned == ["The locus_1 tree was inferred with 2 IQ-TREE threads: a rerun can give other branch "
                      "lengths and supports (--compute_tree_cpus 1 gives the same tree) (trees)."]
