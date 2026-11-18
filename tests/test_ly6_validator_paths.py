"""The LY6 ground-truth check takes its paths as arguments (freeze item 44;
scripts/benchmark/run_ly6_ground_truth_check.py).

It read the curated models from ground_truth/ly6_3ftx/Ly6, relative to the current directory
(a folder that exists nowhere in the repository), wrote its reports into that truth folder, and
ran the validator as "python scripts/benchmark/validate_ly6_results.py": whatever python is on
PATH, from the repository root only. Without the folder it found no file and ended in a
traceback. Every run here is started from tmp_path, not the repository root. All synthetic."""
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SCRIPT = os.path.join(ROOT, "scripts", "benchmark", "run_ly6_ground_truth_check.py")


def _run(tmp_path, *args, env=None):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True,
                          cwd=tmp_path, env=env)


def _run_folder(tmp_path):
    pred = tmp_path / "run" / "plot_inputs_synteny_block_locus_1"
    pred.mkdir(parents=True)
    (pred / "x.fna.gff").write_text(
        "##gff-version 3\n"
        "NC_1\tSynVoy\tmRNA\t1000\t1600\t.\t+\t.\tID=GOI_q|x_b0_l1;SynVoyRole=goi;Confidence=HIGH\n"
        "NC_1\tSynVoy\tCDS\t1000\t1600\t.\t+\t0\tID=GOI_q|x_b0_l1.cds;Parent=GOI_q|x_b0_l1;SynVoyRole=goi\n")
    return tmp_path / "run"


def _truth_folder(tmp_path):
    ly6 = tmp_path / "truth" / "Ly6"
    ly6.mkdir(parents=True)
    (ly6 / "Hosa_NC_1.gff").write_text(
        "##gff-version 3\n"
        "NC_1\tcurated\tgene\t1100\t1500\t.\t+\t.\tID=g1;Name=LY6E\n"
        "NC_1\tcurated\tCDS\t1100\t1500\t.\t+\t0\tID=c1;Parent=g1;Name=LY6E\n")
    (ly6 / "Hosa_NC_1_extraction.gff").write_text("##gff-version 3\n")
    return tmp_path / "truth"


def test_without_gt_dir_the_script_stops_with_usage(tmp_path):
    p = _run(tmp_path, "--outdir", str(_run_folder(tmp_path)))
    assert p.returncode == 2 and "--gt-dir" in p.stderr and "Traceback" not in p.stderr


def test_a_gt_dir_without_curated_gffs_stops(tmp_path):
    empty = tmp_path / "empty" / "Ly6"
    empty.mkdir(parents=True)
    (empty / "only_extraction.gff").write_text("##gff-version 3\n")    # *_extraction.gff do not count
    (empty / "x_extraction.gff").write_text("##gff-version 3\n")
    p = _run(tmp_path, "--gt-dir", str(tmp_path / "empty"), "--outdir", str(_run_folder(tmp_path)))
    assert p.returncode == 2 and "no curated GFF" in p.stderr and "Traceback" not in p.stderr
    p = _run(tmp_path, "--gt-dir", str(tmp_path / "missing"), "--outdir", str(tmp_path / "run"))
    assert p.returncode == 2 and "Traceback" not in p.stderr


def test_a_run_folder_without_plot_inputs_stops(tmp_path):
    p = _run(tmp_path, "--gt-dir", str(_truth_folder(tmp_path)), "--outdir", str(tmp_path / "nothing"))
    assert p.returncode == 2 and "plot_inputs_synteny_block_locus_1" in p.stderr
    assert "Traceback" not in p.stderr


def test_the_report_goes_to_the_report_dir_and_the_truth_is_untouched(tmp_path):
    truth = _truth_folder(tmp_path)
    before = sorted((p.name, p.read_text()) for p in (truth / "Ly6").iterdir())
    run = _run_folder(tmp_path)
    # No python on PATH: the validator has to run under this interpreter, from its own folder.
    env = {k: v for k, v in os.environ.items() if k != "PATH"}
    env["PATH"] = str(tmp_path / "no_bin_here")
    p = _run(tmp_path, "--gt-dir", str(truth), "--outdir", str(run),
             "--report-dir", str(tmp_path / "reports"), env=env)
    assert p.returncode == 0, p.stderr
    report = json.loads((tmp_path / "reports" / "ly6_validation_report.json").read_text())
    assert (tmp_path / "reports" / "ly6_validation_report.md").exists()
    assert len(report["results"]) == 1 and report["results"][0]["cds_tracking_pct"] > 0
    assert report["results"][0]["match"]["gt_file"] == "Hosa_NC_1.gff"
    assert report["skipped_extractions"] == ["Hosa_NC_1_extraction.gff"]
    assert sorted((p.name, p.read_text()) for p in (truth / "Ly6").iterdir()) == before
    assert not (run / "ly6_validation_report.json").exists()


def test_the_ly6_folder_itself_and_the_default_report_dir(tmp_path):
    truth = _truth_folder(tmp_path) / "Ly6"
    run = _run_folder(tmp_path)
    p = _run(tmp_path, "--gt-dir", str(truth), "--outdir", str(run))
    assert p.returncode == 0, p.stderr
    assert json.loads((run / "ly6_validation_report.json").read_text())["results"]
