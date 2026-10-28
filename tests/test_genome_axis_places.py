"""The genome axis, anchor places (profile stage 5a).

A setting that should differ per target genome is set inside the search, where each genome's
distances are derived; the preset modules see only the home side. With --genome_axis (the
pipeline passes it with any --modules other than legacy) and anchor_max_places left to 'auto',
a genome in which the block cap (max_blocks_per_genome) binds gets its blocks built again
once, from the 5 best places of each flanking gene; the cap then cuts as before if it still
binds. Without the axis, or with a value the user set, the blocks are built as before.
Each genome's value is logged, written to target_policy.tsv with its source and to the result
GFF. The synthetic genome below has two flanking genes: their true neighbourhood on chr1 and
seven chance places each on chr2-chr8, ranked so that the 5 best places of the two genes meet
only on chr1 and chr5.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402
import select_modules as sm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


def _hit(query, chrom, start, bits):
    return {"query": query, "chrom": chrom, "start": start, "end": start + 900, "strand": "+",
            "pident": 60.0, "alnlen": 300, "evalue": 1e-20, "bits": float(bits),
            "qstart": 1, "qend": 300}


def _genome():
    """chr1: the neighbourhood (A and B, 500 bits each); chr2-chr8: chance places, A strongest
    on chr2, B strongest on chr8."""
    hits = [_hit("gene-A", "chr1", 100_000, 500), _hit("gene-B", "chr1", 110_000, 500)]
    for i, c in enumerate(range(2, 9)):
        hits.append(_hit("gene-A", f"chr{c}", 100_000, 50 - i))
        hits.append(_hit("gene-B", f"chr{c}", 110_000, 44 + i))
    return hits


def _args(**kw):
    a = argparse.Namespace(max_intron=20_000, cluster_distance=50_000, region_padding=30_000,
                           min_block_genes=2, max_blocks_per_genome=3, genome_axis=False,
                           anchor_max_places=None, genome_ctx={})
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def _opts(proxies=(), goi_seeded=None):
    proxies = set(proxies)

    def seeded(b):
        genes = set(b.get("genes", []))
        return bool(genes & proxies) or any(isr.is_goi_query_id(g) for g in genes)

    return dict(home_rank=None, bridge={}, gene_gap=None, goi_seeded=goi_seeded or seeded,
                proxy_parents=proxies, gw_any_genome=False)


def _chroms(blocks):
    return [b["chrom"] for b in blocks]


def test_without_the_axis_the_cap_cuts_as_before():
    hits = _genome()
    blocks, values = isr.genome_axis_blocks(hits, hits, [], _args(), "g", **_opts())
    assert values is None
    assert _chroms(blocks) == ["chr1", "chr2", "chr3"]       # merge order, cut at 3
    same, n = isr.build_genome_blocks(hits, _args(), "g", **_opts())
    assert _chroms(same) == _chroms(blocks) and n == 8


def test_cap_binding_rebuilds_from_the_best_places(caplog):
    hits = _genome()
    a = _args(genome_axis=True)
    with caplog.at_level("INFO"):
        blocks, values = isr.genome_axis_blocks(hits, hits, [], a, "g", **_opts())
    assert _chroms(blocks) == ["chr1", "chr5"]
    assert values == {"anchor_max_places": {"value": 5, "source": "genome_axis",
                                            "blocks_before_cap": 8, "blocks_rebuilt": 2}}
    text = caplog.text
    assert "the block cap (3) bound at 8 blocks" in text
    assert "anchor_max_places 5 (genome_axis; the cap of 3 bound at 8 blocks; rebuilt" in text
    assert "Capping blocks" not in text                      # the first build is not cut
    assert not (a.genome_ctx.get("limits") or {}).get("blocks_capped")


def test_rebuild_is_capped_when_the_cap_still_binds(caplog):
    hits = _genome()
    a = _args(genome_axis=True, max_blocks_per_genome=1)
    with caplog.at_level("INFO"):
        blocks, values = isr.genome_axis_blocks(hits, hits, [], a, "g", **_opts())
    assert _chroms(blocks) == ["chr1"]
    assert values["anchor_max_places"] == {"value": 5, "source": "genome_axis",
                                           "blocks_before_cap": 8, "blocks_rebuilt": 2}
    assert "Capping blocks: 2 -> 1" in caplog.text          # the rebuild's cap, as before
    assert "rebuilt from the 5 best places: 2 blocks, capped to 1" in caplog.text
    assert a.genome_ctx["limits"]["blocks_capped"] == 1      # counted once, for the rebuild


def test_cap_not_reached_changes_nothing():
    hits = _genome()
    off, _ = isr.genome_axis_blocks(hits, hits, [], _args(max_blocks_per_genome=10), "g", **_opts())
    on, values = isr.genome_axis_blocks(hits, hits, [], _args(genome_axis=True, max_blocks_per_genome=10),
                                        "g", **_opts())
    assert _chroms(on) == _chroms(off) == [f"chr{c}" for c in range(1, 9)]
    assert values == {"anchor_max_places": {"value": 0, "source": "default", "blocks_before_cap": 8}}


def test_no_cap_and_no_flanking_anchors_do_not_rebuild():
    hits = _genome()
    blocks, values = isr.genome_axis_blocks(hits, hits, [], _args(genome_axis=True, max_blocks_per_genome=0),
                                            "g", **_opts())
    assert len(blocks) == 8 and values["anchor_max_places"]["source"] == "default"
    # seeds without flanking anchors (the GOI's own hits): nothing to limit, the cap cuts
    goi = [dict(h, query="GOI_X") for h in hits]
    blocks, values = isr.genome_axis_blocks(goi, [], [], _args(genome_axis=True, min_block_genes=1),
                                            "g", **_opts())
    assert len(blocks) == 3 and values["anchor_max_places"] == {
        "value": 0, "source": "default", "blocks_before_cap": 8}


def test_a_value_the_user_set_wins():
    hits = _genome()
    # 0: every place for every genome; the cap cuts as before
    blocks, values = isr.genome_axis_blocks(hits, hits, [], _args(genome_axis=True, anchor_max_places=0),
                                            "g", **_opts())
    assert _chroms(blocks) == ["chr1", "chr2", "chr3"]
    assert values["anchor_max_places"] == {"value": 0, "source": "user_set", "blocks_before_cap": 8}
    # 2: the search limits the seeds itself (process_single_genome), there is no rebuild
    seeds = isr.limit_anchor_places(hits, 2)
    blocks, values = isr.genome_axis_blocks(seeds, hits, [], _args(genome_axis=True, anchor_max_places=2),
                                            "g", **_opts())
    assert _chroms(blocks) == ["chr1"]
    assert values["anchor_max_places"] == {"value": 2, "source": "user_set", "blocks_before_cap": 1}


def test_goi_proxies_keep_every_place_in_the_rebuild():
    """gene-P is the GOI's own home gene among the anchors: seven single-gene places, each a
    GOI-seeded block. The rebuild limits A and B, never P."""
    hits = _genome() + [_hit("gene-P", f"chr{c}", 50_000, 60 - c) for c in range(9, 16)]
    a = _args(genome_axis=True, max_blocks_per_genome=10)
    blocks, values = isr.genome_axis_blocks(hits, hits, [], a, "g", **_opts(proxies={"gene-P"}))
    assert values["anchor_max_places"]["source"] == "genome_axis"
    assert values["anchor_max_places"]["blocks_before_cap"] == 15
    assert sorted(_chroms(blocks), key=lambda c: int(c[3:])) == ["chr1", "chr5"] + [f"chr{c}" for c in range(9, 16)]
    # without the exemption the limit would cut P to its 5 best places
    assert len([h for h in isr.limit_anchor_places(hits, 5) if h["query"] == "gene-P"]) == 5
    kept = isr.limit_anchor_places(hits, 5, keep={"gene-P"})
    assert [h for h in kept if h["query"] == "gene-P"] == [h for h in hits if h["query"] == "gene-P"]
    assert [id(h) for h in kept] == [id(h) for h in hits if h in kept]      # input order kept


def test_rescued_goi_hits_stay_seeds_of_the_rebuild():
    hits = _genome()
    rescued = [_hit("GOI_X", "chr20", 10_000, 300)]
    a = _args(genome_axis=True)
    blocks, values = isr.genome_axis_blocks(hits + rescued, hits, rescued, a, "g", **_opts())
    assert values["anchor_max_places"]["blocks_before_cap"] == 9
    assert "chr20" in _chroms(blocks)                         # a GOI-seeded block, kept


def test_the_e_value_places_are_not_touched():
    """genome_ctx['places'] (the places of every anchor, for E_syn) is counted before the
    blocks; the rebuild changes the blocks only."""
    hits = _genome()
    ctx = {"places": {"gene-A": 8, "gene-B": 8}}
    a = _args(genome_axis=True, genome_ctx=ctx)
    isr.genome_axis_blocks(hits, hits, [], a, "g", **_opts())
    assert a.genome_ctx["places"] == {"gene-A": 8, "gene-B": 8}


def test_the_value_is_recorded(tmp_path):
    def table(**kw):
        a = argparse.Namespace(output_dir=str(tmp_path), decoy_mode=False, **kw)
        isr._write_target_policy(a, "locus_1_", "g.fna", 61.0, 9, getattr(a, "genome_ctx", None))
        with open(tmp_path / "hits" / "locus_1_g.fna.target_policy.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows[0]

    legacy = table(genome_axis=False, anchor_max_places=None)
    assert list(legacy) == isr.TARGET_POLICY_COLUMNS + isr.LIMIT_COLUMNS      # as before
    values = {"anchor_max_places": {"value": 5, "source": "genome_axis",
                                    "blocks_before_cap": 112, "blocks_rebuilt": 64}}
    row = table(genome_axis=True, anchor_max_places=None, genome_ctx={"axis": values})
    assert list(row) == isr.TARGET_POLICY_COLUMNS + isr.GENOME_AXIS_COLUMNS + isr.LIMIT_COLUMNS
    assert (row["anchor_max_places"], row["anchor_max_places_source"], row["blocks_before_cap"],
            row["blocks_rebuilt"]) == ("5", "genome_axis", "112", "64")
    # before the blocks (or a genome that ends early): the value it starts from
    row = table(genome_axis=True, anchor_max_places=None)
    assert (row["anchor_max_places"], row["anchor_max_places_source"], row["blocks_before_cap"]) == (
        "0", "default", "")
    row = table(genome_axis=True, anchor_max_places=7)
    assert (row["anchor_max_places"], row["anchor_max_places_source"]) == ("7", "user_set")
    header = isr.format_genome_axis_header(values)
    assert header.startswith(isr.GENOME_AXIS_HEADER + " ")
    assert json.loads(header[len(isr.GENOME_AXIS_HEADER):]) == values


def test_auto_is_parsed_and_no_module_may_set_the_key():
    assert isr.count_or_auto("auto") is None and isr.count_or_auto(" AUTO ") is None
    assert isr.count_or_auto("0") == 0 and isr.count_or_auto("5") == 5
    with pytest.raises(ValueError):
        isr.count_or_auto("five")
    assert "anchor_max_places" in sm.GENOME_AXIS_KEYS <= sm.FORBIDDEN_FOR_MODULES
    help_text = subprocess.run([sys.executable, str(ROOT / "bin" / "iterative_search_runner.py"), "--help"],
                               capture_output=True, text=True).stdout
    assert "--genome_axis" in help_text and "--anchor_max_places" in help_text


def test_the_module_passes_the_axis_through_the_helper():
    nf = (ROOT / "modules" / "iterative_search.nf").read_text()
    assert "--anchor_max_places ${GenomeAxis.anchorMaxPlaces(params)}${axis_args} \\\\" in nf
    assert "params.anchor_max_places" not in nf
    config = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*anchor_max_places = 'auto'$", config, re.M)


@pytest.mark.skipif(NEXTFLOW is None, reason="nextflow not on PATH")
def test_the_switch_and_the_legacy_value(tmp_path):
    """Under legacy the search gets '--anchor_max_places 0' and no switch, as before the axis."""
    (tmp_path / "lib").mkdir()
    shutil.copy(ROOT / "lib" / "GenomeAxis.groovy", tmp_path / "lib" / "GenomeAxis.groovy")
    cases = [[None, "auto"], ["", "auto"], ["legacy", "auto"], [" Legacy ", "auto"], ["auto", "auto"],
             ["none", "auto"], ["auto,-dispersed_family", "auto"], ["short_query,B2", "auto"],
             ["legacy", 5], ["auto", 0], ["auto", 5], ["legacy", "0"]]
    lines = ",\n".join(f"        [{'null' if m is None else repr(m)}, {v!r}]" for m, v in cases)
    (tmp_path / "t.nf").write_text(
        "workflow {\n    [\n" + lines + "\n    ].each { c ->\n"
        "        def p = [modules: c[0], anchor_max_places: c[1]]\n"
        "        println \"RESULT ${GenomeAxis.on(p)}|${GenomeAxis.anchorMaxPlaces(p)}|${GenomeAxis.args(p)}|\"\n"
        "    }\n}\n")
    proc = subprocess.run([NEXTFLOW, "run", "t.nf", "-ansi-log", "false"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=300)
    rows = [line[len("RESULT "):] for line in proc.stdout.splitlines() if line.startswith("RESULT ")]
    on = "true|auto| --genome_axis true|"
    assert rows == ["false|0||", "false|0||", "false|0||", "false|0||", on, on, on, on,
                    "false|5||", "true|0| --genome_axis true|", "true|5| --genome_axis true|",
                    "false|0||"], proc.stdout + proc.stderr
