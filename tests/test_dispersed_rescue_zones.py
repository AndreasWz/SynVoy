"""The dispersed-GOI rescue is chosen against the GOI zones of the blocks built without it
(iterative_search_runner.blocks_with_rescue, rescue_selector, select_dispersed_goi_seeds).

Before, the rescue ran before any block existed and skipped every strong GOI hit within the GOI
margin of an anchor, "already covered by a flanking block". A lone anchor forms no block (a block
needs min_block_genes = 2 genes), so a strong hit beside it was searched by nobody: on a sequence
with fewer than 3 anchors it was refused for "too few neighbours", on one with 3 or more it fell
into the dead zone. Now the blocks are built first; a passing hit inside one of their GOI zones is
left to that block, and every other passing hit is a seed inside the envelope of a sequence with
>= 3 anchors, or within the GOI margin of an anchor on a sequence with fewer. Far hits keep the
3-anchor envelope, the 9-pitch reach and min_block_genes = 2.

Distances as genome_scale_args gives them for a gene pitch of 3 kb: block gap and window margin
21 kb, GOI margin 6 kb, envelope margin 27 kb. Synthetic ids only.
"""
import argparse
import logging
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402

PITCH = 3000
BLOCK_GAP, PADDING, GOI_MARGIN, ENVELOPE = 7 * PITCH, 7 * PITCH, 2 * PITCH, 9 * PITCH


def flank(gene, chrom, start, length=1000, bits=500.0):
    return {"query": gene, "chrom": chrom, "start": start, "end": start + length, "strand": "+",
            "pident": 80.0, "evalue": 1e-50, "alnlen": 900, "bits": bits, "qstart": 1, "qend": 300}


def goi(chrom, start, length=600, pident=60.0, evalue=1e-120, alnlen=834, bits=450.0, query="GOI_q"):
    return {"query": query, "chrom": chrom, "start": start, "end": start + length, "strand": "+",
            "pident": pident, "evalue": evalue, "alnlen": alnlen, "bits": bits, "qstart": 1, "qend": 278}


def make_args(**kw):
    a = argparse.Namespace(max_intron=20_000, cluster_distance=BLOCK_GAP, region_padding=PADDING,
                           min_block_genes=2, max_blocks_per_genome=80, genome_axis=False,
                           anchor_max_places=None, genome_ctx={}, goi_margin=GOI_MARGIN,
                           gap_search_window=50_000, dispersed_envelope_margin=ENVELOPE)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def goi_seeded(b):
    return any(isr.is_goi_query_id(g) for g in b.get("genes", []))


def opts(home_rank=None):
    return dict(home_rank=home_rank, bridge={}, gene_gap=None, goi_seeded=goi_seeded,
                proxy_parents=set(), gw_any_genome=False)


def canonical(hits):
    """parse_hits' total order, then one anchor per place, as the search reads its hits."""
    hits = sorted(hits, key=lambda h: (h["query"], h["chrom"], h["start"], h["end"], h["strand"], -h["bits"]))
    return isr.resolve_anchor_hits(hits)


def search(hits, home_rank=None, **kw):
    """The genome's blocks as process_single_genome builds them, the rescue recorded
    (apply_rescue), its searched zones taken as the GOI cores and the genome's end counted."""
    logging.disable(logging.CRITICAL)
    try:
        a = make_args(**kw)
        hits = canonical(hits)
        flanking = [h for h in hits if not isr.is_goi_query_id(h["query"])]
        select = isr.rescue_selector(hits, flanking, a, "g.fna")
        blocks, _values, rescue = isr.genome_axis_blocks(flanking, flanking, a, "g.fna",
                                                        select_rescue=select, **opts(home_rank))
        isr.apply_rescue(a, "g.fna", rescue, blocks)
        a.genome_ctx["goi_cores"] = isr.goi_zones(blocks, GOI_MARGIN)
        isr.log_genome_limits(a, "g.fna")
        return a, blocks, rescue
    finally:
        logging.disable(logging.NOTSET)


def searched(a, h):
    return isr.overlaps_any(h["chrom"], h["start"], h["end"], a.genome_ctx["goi_cores"])


def spans(blocks):
    return [(b["chrom"], b["start"], b["end"], b["genes"]) for b in blocks]


ENVELOPE_ELSEWHERE = [flank(f"gene-{g}", "chr1", 100_000 + i * 10_000) for i, g in enumerate("ABCD")]
LONE = flank("gene-E", "ctg7", 50_000)


def test_a_strong_hit_beside_a_lone_neighbour_is_a_seed():
    hit = goi("ctg7", 51_200)                  # 200 bp after gene-E
    a, blocks, rescue = search(ENVELOPE_ELSEWHERE + [LONE, hit])
    assert [(s["rescue_rule"], s["rescue_anchors"], s["rescue_anchor_bp"]) for s in rescue.seeds] == [
        ("near_neighbour", 1, 200)]
    assert ("ctg7", 50_000, 51_800, ["GOI_q", "gene-E"]) in spans(blocks)
    assert searched(a, hit)
    assert rescue.rejected == [] and isr.rescue_table_rows(a.genome_ctx, a.genome_ctx["goi_cores"])[0] == []
    assert [b["dispersed_rescue"] for b in blocks if b["chrom"] == "ctg7"] == [
        {"rule": "near_neighbour", "anchors": 1, "anchor_bp": 200}]


def test_a_hit_beside_a_lone_neighbour_of_an_envelope_sequence_leaves_the_dead_zone():
    block = [flank(f"gene-{g}", "chrB", 100_000 + i * 5_000) for i, g in enumerate("ABC")]
    lone = flank("gene-D", "chrB", 400_000)    # beyond twice the window margin: no block of its own
    hit = goi("chrB", 401_200)
    a, blocks, rescue = search(block + [lone, hit])
    assert [(s["rescue_rule"], s["rescue_anchors"]) for s in rescue.seeds] == [("envelope", 4)]
    assert ("chrB", 400_000, 401_800, ["GOI_q", "gene-D"]) in spans(blocks)
    assert searched(a, hit) and not a.genome_ctx.get("limits", {}).get("dead_zone")


def test_without_an_envelope_anywhere_the_hit_beside_a_fallback_block_is_left_to_it():
    """Every block is single-gene, so the fallback keeps them all: the hit beside gene-E lies in
    its GOI zone (no refusal record any more; before, it was listed as refused and modelled), the
    far hit is refused, and the quality drops are counted in this genome too."""
    beside, far = goi("ctg7", 51_200), goi("ctg7", 500_000)
    weak = goi("ctg7", 52_000, pident=30.0)
    a, blocks, rescue = search([LONE, beside, far, weak])
    assert spans(blocks) == [("ctg7", 50_000, 51_000, ["gene-E"])]
    assert rescue.seeds == [] and searched(a, beside)
    assert rescue.drops == {"dispersed_identity": 1, "dispersed_in_goi_zone": 1}
    assert [(r["start"], r["rejected_by"], r["envelope"]) for r in rescue.rejected] == [
        (500_000, "no_flanking_envelope_in_genome", "none_in_genome")]


def test_the_refusals_of_far_hits_are_unchanged():
    two = [flank("gene-F", "ctg9", 10_000), flank("gene-G", "ctg9", 300_000)]
    far_on_two = goi("ctg9", 150_000)          # far from both anchors of a 2-anchor contig
    outside = goi("chr1", 900_000)              # outside chr1's envelope
    _a, _blocks, rescue = search(ENVELOPE_ELSEWHERE + two + [far_on_two, outside])
    assert rescue.seeds == []
    assert sorted((r["chrom"], r["rejected_by"]) for r in rescue.rejected) == [
        ("chr1", "outside_flanking_envelope"), ("ctg9", "chromosome_lacks_min_anchors")]
    # without a near radius no hit is seeded beside a lone anchor
    hits = canonical(ENVELOPE_ELSEWHERE + [LONE, goi("ctg7", 51_200)])
    flanking = [h for h in hits if not isr.is_goi_query_id(h["query"])]
    rejected = []
    assert isr.select_dispersed_goi_seeds(hits, flanking, near_radius=0, region_padding=PADDING,
                                          envelope_margin=ENVELOPE, rejected_out=rejected) == []
    assert [r["rejected_by"] for r in rejected] == ["chromosome_lacks_min_anchors"]


def test_a_hit_inside_a_goi_zone_is_no_seed():
    """A hit more than the GOI margin from every anchor but inside the block +- GOI margin (it
    used to be a seed, which made the block GOI-seeded and extended it)."""
    wide = [flank("gene-A", "chrW", 100_000), flank("gene-B", "chrW", 115_000), flank("gene-C", "chrW", 120_000)]
    hit = goi("chrW", 108_000)
    a, blocks, rescue = search(wide + [hit])
    assert rescue.seeds == [] and rescue.drops == {"dispersed_in_goi_zone": 1}
    assert spans(blocks) == [("chrW", 100_000, 121_000, ["gene-A", "gene-B", "gene-C"])]
    assert not goi_seeded(blocks[0]) and "dispersed_rescue" not in blocks[0]
    assert searched(a, hit)


def test_a_far_seed_keeps_the_first_pass_fallback():
    """Three single-gene blocks (the fallback keeps them) and a far hit inside their envelope: the
    seed's block passes the filter, and the first pass's blocks are kept with it, so the hit beside
    gene-A stays in a searched zone (before, the fallback did not run and it fell into the dead
    zone)."""
    single = [flank(f"gene-{g}", "chrZ", 100_000 + i * 100_000) for i, g in enumerate("ABC")]
    beside, far = goi("chrZ", 101_200), goi("chrZ", 250_000)
    a, blocks, rescue = search(single + [beside, far])
    assert [s["start"] for s in rescue.seeds] == [250_000]
    assert spans(blocks) == [("chrZ", 100_000, 101_000, ["gene-A"]), ("chrZ", 200_000, 201_000, ["gene-B"]),
                             ("chrZ", 300_000, 301_000, ["gene-C"]), ("chrZ", 250_000, 250_600, ["GOI_q"])]
    assert searched(a, beside) and searched(a, far)
    assert not a.genome_ctx.get("limits", {}).get("dead_zone")


def test_nothing_to_rescue_gives_the_blocks_of_the_flanking_hits():
    conserved = [flank(f"gene-{g}", "chrK", 200_000 + i * 4_000) for i, g in enumerate("ABCDEF")]
    hits = canonical(conserved + [goi("chrK", 210_500)])
    flanking = [h for h in hits if not isr.is_goi_query_id(h["query"])]
    a = make_args()
    blocks, _values, rescue = isr.genome_axis_blocks(
        flanking, flanking, a, "g", select_rescue=isr.rescue_selector(hits, flanking, a, "g"), **opts())
    plain, _n = isr.build_genome_blocks(flanking, make_args(), "g", **opts())
    assert spans(blocks) == spans(plain) and rescue.seeds == []


def _many_blocks():
    """81 two-gene blocks on contigs c000-c080, then a block (chrB, A-C) whose GOI zone holds a
    passing GOI hit, and a lone neighbour (gene-D) with a strong hit beside it."""
    hits = []
    for i in range(81):
        hits += [flank(f"gene-M{i}", f"c{i:03d}", 10_000), flank(f"gene-N{i}", f"c{i:03d}", 14_000)]
    hits += [flank(f"gene-{g}", "chrB", 100_000 + i * 5_000) for i, g in enumerate("ABC")]
    hits += [flank("gene-D", "chrB", 400_000), goi("chrB", 401_200), goi("chrB", 104_000)]
    return hits


def test_the_cap_keeps_the_seeded_block_and_the_block_holding_an_in_zone_hit():
    for home_rank in (None, {f"gene-{g}": i for i, g in enumerate("ABCD")}):
        a, blocks, rescue = search(_many_blocks(), home_rank=home_rank)
        assert len(blocks) == 80
        assert ("chrB", 400_000, 401_800, ["GOI_q", "gene-D"]) in spans(blocks)
        assert any(b["chrom"] == "chrB" and "gene-A" in b["genes"] for b in blocks)
        assert a.genome_ctx["limits"]["blocks_capped"] == 3
        assert not a.genome_ctx["limits"].get("dead_zone")
        assert [h["start"] for h in rescue.in_zone] == [104_000]


def test_the_cap_keeps_the_build_order_and_does_nothing_below_the_cap():
    blocks = [{"chrom": f"c{i}", "start": 0, "end": 10, "genes": ["g1", "g2"], "genes_count": 2}
              for i in range(5)]
    blocks[3]["genes"] = ["GOI_q", "g1"]
    a = make_args(max_blocks_per_genome=3)
    kept = isr.cap_genome_blocks(list(blocks), a, "g", goi_seeded, set(), None)
    assert [b["chrom"] for b in kept] == ["c0", "c1", "c3"]       # c3 survives, order kept
    a = make_args(max_blocks_per_genome=5)
    assert isr.cap_genome_blocks(list(blocks), a, "g", goi_seeded, set(), None) == blocks


def test_a_genome_without_flanking_hits_keeps_the_blocks_of_its_strongest_hits():
    """In a genome without a flanking hit the GOI's own hits seed every block; when the cap binds,
    the blocks holding the strongest hits survive (not the first ones by sequence name), in the
    build's order."""
    hits = [goi(f"chr{i:03d}", 10_000, bits=100.0 + i) for i in range(82)]
    a = make_args(min_block_genes=2)
    blocks, _values, rescue = isr.genome_axis_blocks(hits, [], a, "g", select_rescue=None,
                                                    survivor_key=isr.strongest_seed_key(hits), **opts())
    chroms = [b["chrom"] for b in blocks]
    assert rescue is None and len(blocks) == 80
    assert "chr000" not in chroms and "chr001" not in chroms
    assert chroms == sorted(chroms) and chroms[-1] == "chr081"
    key = isr.strongest_seed_key(hits)
    assert key({"chrom": "chr005", "start": 9_000, "end": 11_000}) == -105.0
    assert key({"chrom": "chr005", "start": 20_000, "end": 21_000}) == 0.0


def test_a_neighbour_counted_on_each_sequence_is_measured_only():
    """gene-C's best hit is a paralog on chrY; its true copy on chrX (150 bits) does not
    count, so chrX carries two anchors and the strong hit between them is refused. The record says
    so (2 anchors, 3 genes with a hit on chrX), and the counter counts the refusal that counting
    each gene on every sequence would turn into a seed; no call changes."""
    flanking = [flank("gene-A", "chrX", 100_000), flank("gene-B", "chrX", 900_000),
                flank("gene-C", "chrX", 500_000, bits=150.0), flank("gene-C", "chrY", 50_000, bits=160.0)]
    hit = goi("chrX", 300_000)
    hits = canonical(flanking + [hit])
    drops, rejected = {}, []
    seeds = isr.select_dispersed_goi_seeds(hits, [h for h in hits if h["query"] != "GOI_q"],
                                           near_radius=GOI_MARGIN, region_padding=PADDING,
                                           envelope_margin=200_000, rejected_out=rejected, drops=drops)
    assert seeds == []
    assert [(r["rejected_by"], r["anchors_on_sequence"], r["genes_on_sequence"]) for r in rejected] == [
        ("no_flanking_envelope_in_genome", 2, 3)]
    assert drops == {"dispersed_refused_per_piece_pass": 1}
    assert "dispersed_refused_per_piece_pass" in isr.LIMIT_KEYS


def test_the_blocks_and_seeds_do_not_depend_on_the_input_order():
    base = (ENVELOPE_ELSEWHERE + [LONE] + [flank(f"gene-{g}", "chrB", 100_000 + i * 5_000) for i, g in enumerate("ABC")]
            + [flank("gene-D", "chrB", 400_000)]
            + [flank(f"gene-{g}2", "chrZ", 100_000 + i * 100_000) for i, g in enumerate("ABC")]
            + [goi("ctg7", 51_200), goi("chrB", 401_200), goi("chrZ", 250_000), goi("ctg7", 51_300, query="GOI_q|exon_1")])

    def outcome(hits):
        _a, blocks, rescue = search(hits)
        return spans(blocks), [(s["chrom"], s["start"], s["query"]) for s in rescue.seeds]

    ref = outcome(list(base))
    for seed in range(20):
        shuffled = list(base)
        random.Random(seed).shuffle(shuffled)
        assert outcome(shuffled) == ref


# ---------------------------------------------------------------- the genome's search, end to end
def _genome_search(tmp_path, monkeypatch, hits):
    """process_single_genome on a small genome, with the genome-wide search and the window search
    replaced: the window stub records each block it is given, registers the block's GOI core as
    the real one does, and returns a LOW model on ctg7 (and one over the refused hit on ctg9)."""
    from sequence_utils import write_fasta
    for d in ("hits", "regions"):
        (tmp_path / d).mkdir(exist_ok=True)
    genome = str(tmp_path / "Target.fna")
    write_fasta([("chr1", "A" * 400_000), ("ctg7", "C" * 80_000), ("ctg9", "G" * 900_000)], genome)
    db = str(tmp_path / "db.faa")
    write_fasta([(f"gene-{g}", "MKV" * 50) for g in "ABCDE"] + [("GOI_q", "MKT" * 100)], db)
    seen = []

    def window(block_idx, block, _hits, _seqs, _db, genome_name, args, *_a, **_k):
        m = isr.goi_core_margin(args)
        args.genome_ctx.setdefault("goi_cores", []).append((block["chrom"], max(0, block["start"] - m), block["end"] + m))
        args.genome_ctx.setdefault("block_context", {})[block_idx] = ""
        isr.record_block_row(args, block_idx, block, block["start"], block["end"], 1, None)
        seen.append((block["chrom"], block.get("dispersed_rescue")))
        if block["chrom"] != "ctg7":
            return [], []
        attrs = "SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=LOW;EvidenceType=exon_annotation"
        return [], [f"ctg7\tx\tmRNA\t51201\t51800\t.\t+\t.\tID=GOI_q|Target_fna_b{block_idx}_l1;{attrs}",
                    f"ctg9\tx\tmRNA\t450001\t450600\t.\t+\t.\tID=GOI_q|Target_fna_b{block_idx}_l2;{attrs}"]

    ordered = sorted(hits, key=lambda h: (h["query"], h["chrom"], h["start"], h["end"], h["strand"], -h["bits"]))
    monkeypatch.setattr(isr, "run_mmseqs_easy_search_with_retries", lambda **k: (True, "", False))
    monkeypatch.setattr(isr, "parse_hits", lambda *a, **k: [dict(h) for h in ordered])
    monkeypatch.setattr(isr, "process_region_block", window)
    args = argparse.Namespace(output_dir=str(tmp_path), min_identity=10.0, min_length=10, evalue=0.01,
                              cluster_distance=BLOCK_GAP, max_intron=20_000, gap_search_window=50_000,
                              region_padding=PADDING, goi_margin=GOI_MARGIN, synteny_bridge_max_gap=0,
                              quiet_subtools=True, min_block_genes=2, max_blocks_per_genome=80,
                              dispersed_envelope_margin=ENVELOPE, prefix="locus_1")
    result = isr.process_single_genome(genome, db, args, None, "locus_1_", 1)
    return result, seen


def _table(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_the_search_seeds_the_rescue_and_records_what_it_did_not_take(tmp_path, monkeypatch):
    hits = ENVELOPE_ELSEWHERE + [LONE, goi("ctg7", 51_200), goi("ctg9", 450_000),
                                 goi("ctg9", 700_000, pident=25.0, bits=60.0)]
    (name, _genes, _extra, error), seen = _genome_search(tmp_path, monkeypatch, hits)
    assert error is None and name == "Target.fna"
    assert seen == [("chr1", None), ("ctg7", {"rule": "near_neighbour", "anchors": 1, "anchor_bp": 200})]
    rows = _table(tmp_path / "hits" / "locus_1_Target.fna.nonsyntenic.tsv")
    assert [(r["outcome"], r["rejected_by"], r["chrom"], r["start"], r["locus"], r["modelled_here"]) for r in rows] == [
        ("refused", "chromosome_lacks_min_anchors", "ctg9", "450000", "locus_1", "LOW"),
        ("dropped", "dispersed_identity", "ctg9", "700000", "locus_1", "")]
    policy = _table(tmp_path / "hits" / "locus_1_Target.fna.target_policy.tsv")[0]
    assert (policy["rescue_dropped_outside_zones"], policy["limit_dispersed_identity"], policy["limit_dead_zone"]) == (
        "1", "1", "0")
    blocks = _table(tmp_path / "blocks" / "locus_1_Target.fna.blocks.tsv")
    assert [(r["chrom"], r["goi_seeded"], r["dispersed_rescue"]) for r in blocks] == [
        ("chr1", "0", ""), ("ctg7", "1", "near_neighbour")]


def test_a_genome_without_flanking_hits_is_flagged_and_warned_about(tmp_path, monkeypatch, caplog):
    hits = [goi("ctg7", 51_200), goi("ctg9", 450_000)]
    with caplog.at_level(logging.WARNING, logger=isr.logger.name):
        (_name, _genes, _extra, error), seen = _genome_search(tmp_path, monkeypatch, hits)
    assert error is None
    assert seen == [("ctg7", isr.NO_FLANKING_RESCUE), ("ctg9", isr.NO_FLANKING_RESCUE)]
    assert any("No flanking gene was found in this genome" in r.getMessage() for r in caplog.records)
    assert not (tmp_path / "hits" / "locus_1_Target.fna.nonsyntenic.tsv").exists()
    policy = _table(tmp_path / "hits" / "locus_1_Target.fna.target_policy.tsv")[0]
    assert (policy["n_flanking_genes"], policy["blocks_searched"]) == ("0", "2")
