"""Every task of the launcher's default profile fits a machine sold with 8 GB.

Nextflow's local executor does not start a task that asks for more memory than the machine has, and
a laptop sold with 8 GB reports about 7.7 GB. PHYLO_SORT asked for 8 GB in every profile, so a
default run (auto,low_mem) stopped there. Its largest need is the first run on a machine without
ete3's taxonomy database, which ete3 then builds with the whole NCBI taxonomy in memory (about
2.6 GB, measured on 10 % to 35 % of the 2.8 million rows and scaled); with the database on disk it
needs about 130 MB. It asks for 4 GB now. These tests read the requests statically (no Nextflow): the
fixed ones of nextflow.config and the modules, with the low_mem profile's values for the params they
read, and the search's request, which the module caps at the local executor's memory.
"""
import os
import re
import shutil
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(ROOT, "nextflow.config")
FITS_GB = 4.0          # the most any task of auto,low_mem asks for at its first attempt
MACHINE_GB = 7.5       # a laptop sold with 8 GB reports about 7.7 GB


def _code():
    return "\n".join(re.sub(r"(?<![:\"'])//.*$", "", line) for line in open(CONFIG).read().splitlines())


def _gb(text):
    m = re.match(r"^\s*'?([0-9.]+)\s*([KMGT])B?'?\s*$", text.strip().strip('"'))
    assert m, text
    return float(m.group(1)) * {"K": 1 / 2 ** 20, "M": 1 / 1024, "G": 1, "T": 1024}[m.group(2)]


def _block(text, header):
    start = text.index(header)
    depth, j = 0, text.index("{", start)
    while True:
        depth += {"{": 1, "}": -1}.get(text[j], 0)
        if depth == 0:
            return text[start:j + 1]
        j += 1


def _params():
    """The param values under auto,low_mem: the params block, then the two profiles."""
    text = _code()
    values = dict(re.findall(r"^\s{4}([a-z_]+)\s*=\s*([^\n]+?)\s*$", _block(text, "params {"), re.M))
    for profile in ("auto {", "low_mem {"):
        values.update(re.findall(r"params\.([a-z_]+)\s*=\s*([^\n]+?)\s*$", _block(text, profile), re.M))
    return values


def _top_process_block():
    text = _code()
    after_profiles = text[text.index(_block(text, "profiles {")) + len(_block(text, "profiles {")):]
    return _block(after_profiles, "process {")


def _requests():
    """{process: (memory GB, cpus)} of the fixed requests under auto,low_mem."""
    params = _params()
    block = _top_process_block()

    def value(expr):
        expr = expr.strip()
        m = re.match(r"^params\.([a-z_]+)(?:\s+as\s+Integer)?$", expr)
        return params[m.group(1)].strip("'\"") if m else expr.strip("'\"")

    default_mem = _gb(value(re.search(r"^\s{4}memory\s*=\s*(.+)$", block, re.M).group(1)))
    default_cpus = int(value(re.search(r"^\s{4}cpus\s*=\s*(.+)$", block, re.M).group(1)))
    out = {"(default)": (default_mem, default_cpus)}
    for nf in sorted(os.listdir(os.path.join(ROOT, "modules"))):     # a module's own directive ...
        text = open(os.path.join(ROOT, "modules", nf)).read()
        for proc, body in re.findall(r"^process (\w+) \{(.*?)^\}", text, re.S | re.M):
            m = re.search(r"^\s+memory\s+'([^']+)'", body, re.M)
            if m:
                out[proc] = (_gb(m.group(1)), default_cpus)
    for name, body in re.findall(r"withName: '(\w+)' \{(.*?)\n    \}", block, re.S):   # ... the config's
        mem = re.search(r"^\s*memory\s*=\s*(.+)$", body, re.M)
        cpus = re.search(r"^\s*cpus\s*=\s*(.+)$", body, re.M)
        own = out.get(name, (default_mem, default_cpus))
        out[name] = (_gb(value(mem.group(1))) if mem else own[0],
                     int(value(cpus.group(1))) if cpus else own[1])
    # The search's memory is the module's request (lib/SearchMemory.groovy), capped by the local
    # executor's memory: test_the_search_request_stays_within_the_local_executors_memory.
    out.pop("ITERATIVE_SEARCH", None)
    return out


def test_every_fixed_request_of_the_default_profile_fits_an_8_gb_laptop():
    requests = _requests()
    assert {"(default)", "PHYLO_SORT", "LOCATE_GENE", "COMPUTE_TREE", "FETCH_RELATED_GENOMES"} <= set(requests)
    too_big = {p: r for p, r in requests.items() if r[0] > FITS_GB or r[1] > 2}
    assert not too_big, too_big
    assert max(r[0] for r in requests.values()) <= MACHINE_GB


def test_phylo_sort_asks_for_its_measured_need():
    mem, _cpus = _requests()["PHYLO_SORT"]
    assert 3.0 <= mem <= FITS_GB          # the taxonomy build: about 2.6 GB


def test_the_search_request_stays_within_the_local_executors_memory():
    nf = open(os.path.join(ROOT, "modules", "iterative_search.nf")).read()
    assert "memory { def m = SearchMemory.request(params, genomes, task)" in nf
    lib = open(os.path.join(ROOT, "lib", "SearchMemory.groovy")).read()
    assert "if ((task.executor ?: 'local').toString() == 'local')" in lib and "mem = local" in lib
    assert "withName: 'ITERATIVE_SEARCH'" in _top_process_block()
    body = re.search(r"withName: 'ITERATIVE_SEARCH' \{(.*?)\n    \}", _top_process_block(), re.S).group(1)
    assert not re.search(r"^\s*memory\s*=", body, re.M)    # the module's request is not overridden


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")
def test_nextflow_resolves_the_same_phylo_sort_request(tmp_path):
    proc = subprocess.run(["nextflow", "config", ROOT, "-profile", "auto,low_mem", "-flat"],
                          capture_output=True, text=True, timeout=300, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "process.'withName:PHYLO_SORT'.memory = '4 GB'" in proc.stdout or \
        "process.withName:PHYLO_SORT.memory = '4 GB'" in proc.stdout, \
        [l for l in proc.stdout.splitlines() if "PHYLO_SORT" in l]
