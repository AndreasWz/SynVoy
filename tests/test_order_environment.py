"""Fix-phase d8 (2 Nov; items 24 and 33): no result depends on the order things arrive in, the string hash, the locale or a
target's file size.

- Tied home anchors: synvoy_scale built its gene table from a set (string-hash order), and synvoy_evalue.ranks_around_goi
  sorted on one coordinate, so genes tied on it kept that order. PYTHONHASHSEED = 0 fixed it for one Python version only
  (3.10 hashes strings with SipHash-2-4, 3.11 and later with SipHash-1-3).
- The locale: STAGE_GENOMES's globs and sorts, PLOT_SYNTENY's sort of the home scaffolds, ASSESS_GENOME_QUALITY's find.
- File-system order: the GFF globs of the two downloaders. Completion order: the target genomes' and the locus BEDs'
  collect().
- The tools' versions are recorded (tool_versions in effective_params.json and the report's run_settings).
- Item 33: a target's size in the distances is its bases, not its file size (test_search_core_round2.py).
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
sys.path.insert(0, BIN)
import synvoy_evalue as se  # noqa: E402
import synvoy_scale as ss  # noqa: E402

# Six anchors: two pairs tied on their start and end, one on each side of the GOI.
HEADERS = [f"{g} home={s}-{e} genome=1000000 goi=5000-6000" for g, s, e in (
    ("gene-zeta", 1000, 2000), ("gene-alpha", 1000, 2000), ("gene-mid", 3000, 4000),
    ("gene-yak", 8000, 9000), ("gene-beta", 8000, 9000), ("gene-far", 10000, 11000))]


def test_the_gene_table_is_in_name_order():
    g = ss.home_geometry_from_headers(HEADERS)
    assert list(g["genes"]) == sorted(g["genes"])


def test_tied_anchors_rank_by_name():
    g = ss.home_geometry_from_headers(HEADERS)
    ranks = se.ranks_around_goi(g["genes"], g["goi"])
    assert ranks == {"gene-mid": -1, "gene-alpha": -2, "gene-zeta": -3,
                     "gene-beta": 1, "gene-yak": 2, "gene-far": 3}


def test_a_closer_start_ranks_first_on_a_tied_end():
    genes = {"a": {"start": 100, "end": 900}, "b": {"start": 500, "end": 900}}
    assert se.ranks_around_goi(genes, (1000, 2000)) == {"b": -1, "a": -2}


def test_the_ranks_do_not_follow_the_hash_seed():
    code = ("import sys, json; sys.path.insert(0, %r); import synvoy_scale as ss, synvoy_evalue as se; "
            "g = ss.home_geometry_from_headers(%r); "
            "print(json.dumps([list(g['genes']), g['order'], se.ranks_around_goi(g['genes'], g['goi'])]))" % (BIN, HEADERS))
    outs = set()
    for seed in ("0", "1", "2", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        outs.add(subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True).stdout)
    assert len(outs) == 1


def _nf(name):
    return open(os.path.join(ROOT, "modules", name)).read()


def test_the_shell_sorts_and_globs_are_byte_ordered():
    assert "export LC_ALL=C\n    mkdir staged_genomes" in _nf("stage_genomes.nf")
    plot = _nf("plot_synteny.nf")
    assert "| LC_ALL=C sort -u | head -1" in plot and plot.count("LC_ALL=C sort -n") == 2
    assert "| LC_ALL=C sort > genome_files.txt" in _nf("assess_quality.nf")


def test_the_collects_are_sorted():
    main = open(os.path.join(ROOT, "main.nf")).read()
    assert "resolveTargetGenomeFiles(params.target_genomes)).collect(sort: true)" in main
    sub = open(os.path.join(ROOT, "subworkflows", "adjudicate_and_report.nf")).read()
    assert "locus_beds.flatten().collect(sort: true)" in sub


def test_the_plots_region_beds_are_in_name_order():
    """The loudness session's list item 1: PLOT_SYNTENY's region BEDs came in task completion order
    (groupTuple of one CLUSTER_REGIONS task per genome and locus), so the plot task's command line
    and its -resume key differed from run to run (39 of 40 plot commands between two arms)."""
    main = open(os.path.join(ROOT, "main.nf")).read()
    assert "def beds = entry[4].toSorted { b -> b.name }" in main


def test_the_downloaders_take_the_gff_in_name_order(tmp_path):
    import fetch_related_genomes as fr
    d = tmp_path / "x"
    (d / "b").mkdir(parents=True)
    (d / "a").mkdir()
    for sub in ("b", "a"):
        (d / sub / "genomic.gff").write_text(f"##gff-version 3\nc\t{sub}\tgene\t1\t9\t.\t+\t.\tID=g\n")
    target = tmp_path / "t.gff"
    assert fr._extract_best_gff(d, target) == "ok"
    assert "\ta\tgene" in target.read_text()
    src = open(os.path.join(BIN, "download_genome.py")).read()
    assert 'for f in sorted(extract_dir.rglob(pat))]' in src


def test_the_tool_versions_are_recorded(tmp_path):
    out = tmp_path / "e.json"
    subprocess.run([sys.executable, os.path.join(BIN, "resolve_effective_params.py"), "--defaults", "{}",
                    "--tool_versions", "--output", str(out)], check=True, capture_output=True)
    versions = json.loads(out.read_text())["tool_versions"]
    assert {"python", "parasail", "mmseqs", "tblastn", "miniprot", "samtools", "mafft", "prodigal"} <= set(versions)
    assert versions["python"] == sys.version.split()[0]
    assert "--tool_versions" in _nf("resolve_effective_params.nf")
    assert "nextflow_version: workflow.nextflow.version" in _nf("resolve_effective_params.nf")
