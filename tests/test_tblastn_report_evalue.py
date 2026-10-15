"""tblastn in a search window reports to E 1000 (iterative_search_runner.TBLASTN_REPORT_EVALUE).

The cutoff is not only an output filter: tblastn links HSPs by sum statistics only among
the ones it keeps, so the old cutoff of 10 (min(aug_relaxed_evalue_cap, search_evalue x
aug_relaxed_evalue_mult)) also changed the E-value of what it reported. The kept threshold
of the window search (E <= search_evalue x aug_relaxed_parse_evalue_mult = 0.1, identity
>= 15 %) is applied when the output is read, as before.
"""
import os
import sys
import types

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "bin"))
import iterative_search_runner as isr  # noqa: E402


def _args(tmp_path, **over):
    a = dict(output_dir=str(tmp_path), evalue=0.01, aug_relaxed_evalue_mult=1000.0,
             aug_relaxed_evalue_cap=10.0, aug_relaxed_parse_evalue_mult=10.0,
             aug_relaxed_identity_min=15.0, aug_relaxed_identity_factor=0.6, min_identity=10.0,
             aug_relaxed_length_min=15, aug_relaxed_length_div=2, min_length=10,
             enable_smith_waterman=False, deterministic_goi_search=True)
    a.update(over)
    return types.SimpleNamespace(**a)


def _run(tmp_path, monkeypatch, args):
    region = tmp_path / "region.fa"
    region.write_text(">region_seq\n" + "ACGT" * 500 + "\n")
    seen = {}

    def fake_run(cmd, *a, **k):
        if cmd[0] == "tblastn":
            seen["tblastn"] = list(cmd)
            out = cmd[cmd.index("-out") + 1]
            with open(out, "w") as fh:
                # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                fh.write("GOI_q\tregion_seq\t30.0\t40\t28\t0\t1\t40\t100\t219\t0.05\t35.0\n")
                fh.write("GOI_q\tregion_seq\t30.0\t40\t28\t0\t41\t80\t400\t519\t0.5\t28.0\n")
                fh.write("GOI_q\tregion_seq\t12.0\t40\t35\t0\t1\t40\t800\t919\t0.01\t30.0\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(isr.subprocess, "run", fake_run)
    monkeypatch.setattr(isr, "run_augmented_mmseqs_with_retries", lambda **k: (False, "skipped in test"))
    hits = isr.run_augmented_search(str(region), [{"id": "GOI_q", "seq": "M" + "A" * 79}],
                                    "g1", args, "t1", 1)
    return seen, hits


def test_tblastn_reports_to_1000(tmp_path, monkeypatch):
    seen, _hits = _run(tmp_path, monkeypatch, _args(tmp_path))
    cmd = seen["tblastn"]
    assert float(cmd[cmd.index("-evalue") + 1]) == 1000.0


def test_cutoff_does_not_follow_the_window_settings(tmp_path, monkeypatch):
    """The MMseqs2 cutoff settings no longer reach tblastn (they gave 10 at the defaults)."""
    seen, _hits = _run(tmp_path, monkeypatch, _args(tmp_path, aug_relaxed_evalue_cap=0.5,
                                                    aug_relaxed_evalue_mult=10.0))
    cmd = seen["tblastn"]
    assert float(cmd[cmd.index("-evalue") + 1]) == isr.TBLASTN_REPORT_EVALUE == 1000.0


def test_kept_threshold_is_unchanged(tmp_path, monkeypatch):
    """Of the reported rows only E <= 0.1 and identity >= 15 % are kept."""
    _seen, hits = _run(tmp_path, monkeypatch, _args(tmp_path))
    kept = sorted((h["start"], h["evalue"]) for h in hits)
    assert kept == [(99, 0.05)]
    assert all(h["method"] == "tblastn" for h in hits)
