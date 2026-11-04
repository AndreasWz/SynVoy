"""Fix-phase d6 (2 Nov; item 22): no resource limit changes a result without the result saying so.

- The MMseqs2 low-memory ladders (the genome search and each window) record a lower rung: limit counters
  (mmseqs_genome_lower_rung, mmseqs_window_lower_rung), a ##synvoy-reduced-search line in the genome's GFF, a report
  warning. A window whose ladder is spent stops the search (it kept its tblastn and Smith-Waterman hits only).
- effective_params.json records the resources the run used (run_resources: the Nextflow profile, MMseqs2 sensitivity and
  split limit, CPUs, memory, skip_tree), and the report's run_settings carries them.
- IQ-TREE runs with a fixed seed; the seed, threads and the bootstrap, also one dropped after a failure, are recorded
  (<locus>_tree.info.json), named in the report (trees, a warning) and on the tree page.
- A memory kill exits 137 (retried): a killed miniprot in the rescues and in the search, a killed Augustus, a search
  worker killed under the pool. docker_max's search memory grows with the attempt.
"""
import json
import os
import subprocess
import sys
import types

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import iterative_search_runner as isr  # noqa: E402
import sequence_utils as su  # noqa: E402


def _proc(rc, stderr=""):
    return types.SimpleNamespace(returncode=rc, stdout="", stderr=stderr)


OOM = "Error: Prefilter died\nCannot fit databases into 2G. Please use a computer with more main memory.\n"


def _args(tmp_path):
    return types.SimpleNamespace(output_dir=str(tmp_path), mmseqs_sens=7.5, mmseqs_split_memory_limit="2G",
                                 evalue=0.01, mmseqs_verbosity=0, aug_relaxed_evalue_cap=10.0,
                                 aug_relaxed_evalue_mult=1000.0, window_mmseqs_min_length=30,
                                 genome_ctx={})


def _fake_mmseqs(fail_first, record):
    calls = {"search": 0}

    def run(cmd, *a, **k):
        if cmd[1] == "search":
            calls["search"] += 1
            record.append(cmd)
            if calls["search"] <= fail_first:
                return _proc(1, OOM)
        if cmd[1] == "convertalis":
            open(cmd[5], "w").close()
        return _proc(0)
    return run


# --------------------------------------------------------------------------- the ladders

def test_the_genome_ladder_hands_back_the_rung_that_ran(tmp_path, monkeypatch):
    cmds = []
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(1, cmds))
    used = {}
    ok, _d, _r = isr.run_mmseqs_easy_search_with_retries("q.faa", "g.fna", str(tmp_path / "h.m8"), str(tmp_path / "t"),
                                                         _args(tmp_path), 4, "G", used=used)
    assert ok and used == {"rung": "lowmem_retry_1", "sensitivity": 7.0, "threads": 1, "split": "2G", "lower": True}
    assert cmds[1][cmds[1].index("-s") + 1] == "7.0"


def test_the_first_rung_is_no_lower_rung(tmp_path, monkeypatch):
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(0, []))
    used = {}
    isr.run_mmseqs_easy_search_with_retries("q.faa", "g.fna", str(tmp_path / "h.m8"), str(tmp_path / "t"),
                                            _args(tmp_path), 4, "G", used=used)
    assert used["lower"] is False and used["sensitivity"] == 7.5


def test_a_window_on_a_lower_rung_is_counted_and_named(tmp_path, monkeypatch):
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(1, []))
    args = _args(tmp_path)
    ok, _d = isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"), str(tmp_path / "t"),
                                                   args, 1, "G_b3")
    assert ok
    assert args.genome_ctx["limits"]["mmseqs_window_lower_rung"] == 1
    assert args.genome_ctx["reduced_search"]["windows"] == [
        {"window": "G_b3", "rung": "aug_lowmem_auto", "sensitivity": 7.0, "threads": 1, "split": "2G"}]


def test_a_spent_window_ladder_stops_the_search(tmp_path, monkeypatch):
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", lambda **k: (False, "Prefilter died"))
    monkeypatch.setattr(isr.subprocess, "run", lambda *a, **k: _proc(0))
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6,
             min_identity=10.0, aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=False, deterministic_goi_search=True, genome_ctx={})
    with pytest.raises(su.StepFailed) as stop:
        isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": "M" + "A" * 119}], "G_b3",
                                 types.SimpleNamespace(**a), "t1", 1, region_len=2000)
    assert "MMseqs2 search of a window failed: Prefilter died" in str(stop.value)
    assert not [p for p in os.listdir(tmp_path) if p.startswith("tmp_")]


def test_the_counters_are_limit_columns():
    assert {"mmseqs_genome_lower_rung", "mmseqs_window_lower_rung"} <= set(isr.LIMIT_KEYS)
    assert isr.REDUCED_SEARCH_HEADER == "##synvoy-reduced-search"


def _report(tmp_path, extra_cols=None, effective=None, trees=None):
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "hits").mkdir()
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    cols = {"genome": "Gx", "median_flank_identity": "80.0", "n_flanking_genes": "20", "blocks_searched": "3",
            "significant_blocks": "1", "best_block_esyn": "1e-9", "best_block_context": "complete",
            "truncated_blocks": "0", "genetic_code": "1"}
    cols.update({f"limit_{k}": "0" for k in isr.LIMIT_KEYS})
    cols.update(extra_cols or {})
    (stage / "hits" / "locus_1_Gx.fna.target_policy.tsv").write_text(
        "\t".join(cols) + "\n" + "\t".join(cols.values()) + "\n")
    argv = [sys.executable, os.path.join(BIN, "generate_report.py"), "--results_dir", str(stage),
            "--output", str(tmp_path / "r.json"), "--allow-empty"]
    if effective:
        (tmp_path / "e.json").write_text(json.dumps(effective))
        argv += ["--effective_params", str(tmp_path / "e.json")]
    if trees:
        (tmp_path / "ti").mkdir()
        for locus, info in trees.items():
            (tmp_path / "ti" / f"{locus}_tree.info.json").write_text(json.dumps(info))
        argv += ["--tree_info_dir", str(tmp_path / "ti")]
    p = subprocess.run(argv, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return json.loads((tmp_path / "r.json").read_text())


def test_the_report_warns_on_a_lower_rung(tmp_path):
    rep = _report(tmp_path, {"limit_mmseqs_window_lower_rung": "2"})
    assert rep["limits"]["counts"]["mmseqs_window_lower_rung"] == {"total": 2, "genomes": ["Gx"]}
    assert any("lower rung of its low-memory ladder" in w and "2 window search(es) (Gx)" in w for w in rep["warnings"])


def test_no_lower_rung_no_warning(tmp_path):
    assert not any("lower rung" in w for w in _report(tmp_path)["warnings"])


# --------------------------------------------------------------------------- run_resources

def test_the_resolver_writes_the_run_resources(tmp_path):
    out = tmp_path / "e.json"
    p = subprocess.run([sys.executable, os.path.join(BIN, "resolve_effective_params.py"), "--defaults", "{}",
                        "--run_resource", "nextflow_profile=auto,low_mem", "--run_resource", "mmseqs_sensitivity=7.0",
                        "--run_resource", "iterative_search_memory=4 GB", "--output", str(out)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert json.loads(out.read_text())["run_resources"] == {
        "nextflow_profile": "auto,low_mem", "mmseqs_sensitivity": "7.0", "iterative_search_memory": "4 GB"}


def test_the_report_carries_the_run_resources(tmp_path):
    rep = _report(tmp_path, effective={"settings": {}, "run_resources": {"nextflow_profile": "auto,low_mem"}})
    assert rep["run_settings"]["run_resources"] == {"nextflow_profile": "auto,low_mem"}


def test_the_module_passes_the_resources():
    nf = open(os.path.join(ROOT, "modules", "resolve_effective_params.nf")).read()
    for key in ("nextflow_profile: workflow.profile", "mmseqs_sensitivity: params.mmseqs_sensitivity",
                "mmseqs_split_memory_limit:", "iterative_search_cpus:", "skip_tree:", "${resources_arg}"):
        assert key in nf


# --------------------------------------------------------------------------- IQ-TREE

def test_iqtree_runs_with_a_fixed_seed(tmp_path, monkeypatch):
    import compute_tree as ct
    seen = []

    def run(cmd, *a, **k):
        seen.append(cmd)
        open(cmd[cmd.index("--prefix") + 1] + ".treefile", "w").write("(a,b,c);\n")
        return _proc(0)
    monkeypatch.setattr(ct.subprocess, "run", run)
    ct.run_iqtree("x.aln", str(tmp_path / "locus_1_tree.nwk"), threads=2, bootstrap_reps=1000)
    assert seen[0][seen[0].index("--seed") + 1] == str(ct.IQTREE_SEED) and "-B" in seen[0]


def _tree_main(tmp_path, monkeypatch, fail_with_bootstrap):
    import compute_tree as ct
    faa = tmp_path / "goi.faa"
    faa.write_text("".join(f">GOI_x|s{i}\nMKTAYIAKQR{'A' * i}\n" for i in range(5)))
    monkeypatch.setattr(ct, "run_mafft", lambda *a, **k: None)

    def iq(aln, out, threads=1, bootstrap_reps=1000, seed=1):
        if bootstrap_reps and fail_with_bootstrap:
            raise subprocess.CalledProcessError(2, ["iqtree2"], output="", stderr="ERROR: bootstrap failed here")
        open(out, "w").write("(a,b,c);\n")
    monkeypatch.setattr(ct, "run_iqtree", iq)
    monkeypatch.setattr(sys, "argv", ["compute_tree.py", "--input", str(faa), "--output",
                                      str(tmp_path / "locus_2_tree.nwk"), "--threads", "3"])
    ct.main()
    return json.loads((tmp_path / "locus_2_tree.info.json").read_text())


def test_the_tree_records_how_it_was_made(tmp_path, monkeypatch):
    info = _tree_main(tmp_path, monkeypatch, False)
    assert info == {"locus": "locus_2", "seed": 1, "threads": 3, "sequences": 5,
                    "bootstrap": {"status": "ultrafast", "replicates": 1000, "reason": ""}}


def test_a_dropped_bootstrap_is_recorded(tmp_path, monkeypatch):
    info = _tree_main(tmp_path, monkeypatch, True)
    assert info["bootstrap"]["status"] == "dropped_after_failure" and info["bootstrap"]["replicates"] == 0
    assert "1000 ultrafast replicates failed (exit 2: ERROR: bootstrap failed here)" in info["bootstrap"]["reason"]


def test_the_report_names_a_dropped_bootstrap(tmp_path):
    boot = {"status": "dropped_after_failure", "replicates": 0, "reason": "IQ-TREE with 1000 ... failed"}
    rep = _report(tmp_path, trees={"locus_1": {"locus": "locus_1", "seed": 1, "threads": 2, "bootstrap": boot}})
    assert rep["trees"]["locus_1"]["bootstrap"]["status"] == "dropped_after_failure"
    assert any("The locus_1 tree has no branch support" in w for w in rep["warnings"])


def test_the_tree_page_names_the_support(tmp_path):
    import plot_synteny as ps
    rep = tmp_path / "r.json"
    rep.write_text(json.dumps({"trees": {"locus_1": {"seed": 1, "threads": 2, "bootstrap": {
        "status": "dropped_after_failure", "replicates": 0, "reason": "IQ-TREE ... failed"}}}}))
    assert ps._tree_support_note(str(rep), "synteny_block_locus_1_tree.html") == "none: IQ-TREE ... failed"
    assert ps._tree_support_note(str(rep), "synteny_block_locus_2_tree.html") is None
    assert ps._tree_support_note(None, "synteny_block_locus_1_tree.html") is None


# --------------------------------------------------------------------------- memory kills

@pytest.mark.parametrize("rc,stderr,killed", [(-9, "", True), (137, "", True), (1, "MemoryError", True),
                                              (1, "terminate called after throwing an instance of 'std::bad_alloc'", True),
                                              (1, "some other error", False), (0, "", False)])
def test_what_counts_as_a_memory_kill(rc, stderr, killed):
    assert su.memory_killed(rc, stderr) is killed


def test_a_memory_kill_exits_137_and_a_guard_3():
    def oom():
        raise su.MemoryKilled("miniprot ran out of memory")
    for wrapper in (su.exit_step_failed, su.exit_on_memory_kill):
        with pytest.raises(SystemExit) as exc:
            wrapper(oom)
        assert exc.value.code == 137
    with pytest.raises(SystemExit) as exc:
        su.exit_step_failed(lambda: (_ for _ in ()).throw(su.StepFailed("x")))
    assert exc.value.code == 3


def test_a_killed_rescue_miniprot_is_no_no_model(tmp_path, monkeypatch):
    import rescue_strong_synteny as rss
    monkeypatch.setattr(rss.subprocess, "run", lambda *a, **k: _proc(-9))
    with pytest.raises(su.MemoryKilled, match="ran out of memory"):
        rss._run_miniprot_relaxed("w.fna", "q.faa", 0.05, 60, place="G chr1:1-100")
    monkeypatch.setattr(rss.subprocess, "run", lambda *a, **k: _proc(1, "bad index"))
    with pytest.raises(RuntimeError, match=r"failed \(exit 1\): bad index"):
        rss._run_miniprot_relaxed("w.fna", "q.faa", 0.05, 60, place="G chr1:1-100")


def test_a_killed_search_miniprot_raises(monkeypatch):
    import annotate_goi_exons as ag
    monkeypatch.setattr(ag, "MINIPROT_AVAILABLE", True)
    monkeypatch.setattr(ag.subprocess, "run", lambda *a, **k: _proc(-9))
    with pytest.raises(su.MemoryKilled, match="ran out of memory"):
        ag._run_miniprot("M" + "A" * 50, "ACGT" * 100, "c1", sensitive=True)


def test_a_killed_augustus_raises(monkeypatch, tmp_path):
    import gene_predictor as gp
    monkeypatch.setattr(gp, "resolve_augustus_species", lambda s: s)
    monkeypatch.setattr(gp.subprocess, "run", lambda *a, **k: _proc(-9))
    with pytest.raises(su.MemoryKilled, match="Augustus .* ran out of memory"):
        gp.run_augustus(str(tmp_path / "x.fa"), str(tmp_path / "x.faa"), str(tmp_path / "x.gff"), species="fly")


def test_the_wave_exits_137_on_any_memory_kill_and_a_killed_worker_is_one():
    src = open(isr.__file__).read()
    assert "if all('ran out of memory' in msg for _, msg in wave_errors):\n                sys.exit(137)" in src
    assert "isinstance(exc, concurrent.futures.BrokenExecutor)" in src
    assert "often because it ran out of memory (retried with more" in src


def test_docker_max_lets_the_search_memory_grow():
    cfg = open(os.path.join(ROOT, "nextflow.config")).read()
    block = cfg[cfg.index("    docker_max {"):cfg.index("    singularity {")]
    search = block[block.index("withName: 'ITERATIVE_SEARCH'"):]
    search = search[:search.index("}")]
    assert "memory =" not in search


def test_a_rung_that_changes_only_threads_is_no_lower_rung(tmp_path, monkeypatch):
    """First review: with sensitivity 7.0 the first retry differs only in its thread count, which changes only the row
    order; it is not counted as a lower rung."""
    monkeypatch.setattr(isr.subprocess, "run", _fake_mmseqs(1, []))
    args = _args(tmp_path)
    args.mmseqs_sens = 7.0
    used = {}
    isr.run_mmseqs_easy_search_with_retries("q.faa", "g.fna", str(tmp_path / "h.m8"), str(tmp_path / "t"), args, 4, "G",
                                            used=used)
    assert used["rung"] == "lowmem_retry_1" and used["lower"] is False
    isr.run_augmented_mmseqs_with_retries("v.faa", "r.fna", str(tmp_path / "a.m8"), str(tmp_path / "t2"), args, 1, "G_b3")
    assert "limits" not in args.genome_ctx and "reduced_search" not in args.genome_ctx


def test_a_lower_rung_without_hits_is_still_counted():
    src = open(isr.__file__).read()
    assert "{'limits': {'mmseqs_genome_lower_rung': 1}} if genome_rung.get('lower') else None" in src


def test_borrow_annotations_retries_a_memory_kill():
    assert "exit_on_memory_kill(main)" in open(os.path.join(BIN, "borrow_annotations.py")).read()
