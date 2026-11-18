"""The GOI loci a window's cap cuts are logged with their positions.

cap_goi_parent_loci keeps the goi_candidate_regions_per_query strongest GOI loci of one query in a
window (2; the local_family module: home copies + 1). The cut loci were counted (regions_cap,
regions_cap_rel, regions_cap_decoy) but their positions were logged nowhere. Now one INFO line per
window and query names them, the strongest first, with the weakest locus kept; decoys are left out.
The loci kept, and their order, do not change.
"""
import logging
import os
import sys
import types

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402


def _args(cap=2):
    return types.SimpleNamespace(genome_ctx={}, goi_candidate_regions_per_query=cap)


def _locus(gstart, bits, n=1, strand="+", qcov=0.5):
    return [{"gstart": gstart + 100 * i, "gend": gstart + 100 * i + 90, "bits": bits, "strand": strand,
             "qcov": qcov} for i in range(n)]


def _rank(hits):
    return (hits[0]["qcov"], max(h["bits"] for h in hits), len(hits))


WHERE = ("G.fna", 3, "chr2", 1_000_000)


def _lines(caplog):
    return [r.getMessage() for r in caplog.records if "GOI loci" in r.getMessage()]


def test_the_cut_loci_are_named_and_the_kept_ones_do_not_change(caplog):
    loci = [_locus(0, 30.0), _locus(5_000, 80.0, n=2, strand="-"), _locus(9_000, 50.0), _locus(20_000, 20.0)]
    old = sorted(loci, key=_rank, reverse=True)[:2]
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        kept = isr.cap_goi_parent_loci(_args(), "GOI_q", loci, _rank, where=WHERE)
    assert kept == old
    assert _lines(caplog) == [
        "[G.fna] block 3: GOI_q: kept 2 of 4 GOI loci (goi_candidate_regions_per_query 2); cut: "
        "chr2:1000001-1000090 (+, 1 hits, qcov 0.50, 30.0 bits); "
        "chr2:1020001-1020090 (+, 1 hits, qcov 0.50, 20.0 bits); "
        "weakest kept: chr2:1009001-1009090 (+, 1 hits, qcov 0.50, 50.0 bits)"]
    assert args_counts(loci) == {"regions_cap": 2}


def args_counts(loci, parent="GOI_q"):
    a = _args()
    isr.cap_goi_parent_loci(a, parent, loci, _rank, where=WHERE)
    return a.genome_ctx.get("limits", {})


def test_nothing_cut_logs_nothing(caplog):
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.cap_goi_parent_loci(_args(cap=4), "GOI_q", [_locus(0, 30.0), _locus(9_000, 50.0)], _rank, where=WHERE)
    assert _lines(caplog) == []


def test_a_decoy_is_counted_and_not_logged(caplog):
    loci = [_locus(0, 30.0), _locus(5_000, 80.0), _locus(9_000, 50.0)]
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        assert args_counts(loci, parent="GOI_q~rev") == {"regions_cap_decoy": 1}
    assert _lines(caplog) == []


def test_a_long_cut_list_names_the_first_ten(caplog):
    loci = [_locus(10_000 * i, float(100 - i)) for i in range(14)]
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.cap_goi_parent_loci(_args(), "GOI_q~rel", loci, _rank, where=WHERE)
    (line,) = _lines(caplog)
    assert "kept 2 of 14 GOI loci" in line and line.count("bits)") == 11 and "; and 2 more;" in line


def test_without_a_place_the_cap_logs_nothing(caplog):
    with caplog.at_level(logging.INFO, logger=isr.logger.name):
        isr.cap_goi_parent_loci(_args(), "GOI_q", [_locus(0, 30.0)] * 3, _rank)
    assert _lines(caplog) == []


def test_the_window_passes_its_place():
    src = open(os.path.join(ROOT, "bin", "iterative_search_runner.py")).read()
    body = src.split("def process_region_block(", 1)[1].split("\ndef ", 1)[0]
    assert "where=(genome_name, block_idx, chrom, w_start))" in body
