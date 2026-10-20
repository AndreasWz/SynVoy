"""blocks/<genome>.blocks.tsv (choosing the figures' regions with the label's scoring, step 1,
2026-10-19): one row per block the search cut a window for, with its neighbourhood, its
E_syn at full precision and the final labels of its GOI models. It changes no call.

The fixture runs the genome's end as the search does: E_syn from block_place_terms,
BlockEsyn written by annotate_place_evalues (blastp stubbed), the multiplier
(finalize_place_evalues) and the final labels (relabel_goi_lines).
"""
import csv
import os
import re
import sys
from types import SimpleNamespace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import iterative_search_runner as isr  # noqa: E402
import synvoy_evalue  # noqa: E402

GENOME = "Candida_b2.fna"     # the genome tag itself holds "_b2_"
TAG = "Candida_b2_fna"
E_OF = {"A": 1e-20, "B": 1e-4, "C": 0.5, "W": 1e-8}   # blastp E per candidate (stub)


def _args(tmp_path):
    ranks = {f"g{i}": i for i in range(1, 11)}
    return SimpleNamespace(output_dir=str(tmp_path),
                           genome_ctx={"genome_len": 20_000_000, "places": {g: 2 for g in ranks},
                                       "block_esyn": {}},
                           home_rank=ranks, goi_ranks={}, home_geometry={"pitch": 20_000},
                           expansion=1.0)


def _cand(model_id, chrom, start, end, subject, role="goi", qcov=0.9, bits=80.0):
    attrs = (f"ID={model_id};Name=GOI_q;SynVoy_Parent=GOI_q;SynVoyRole={role};"
             f"EvidenceType=exon_annotation;Confidence=LOW;QueryCoverage={qcov};Bits={bits}")
    return {"start": start, "end": end, "label_subjects": [subject],
            "gff": [f"{chrom}\texon_annotation\tmRNA\t{start}\t{end}\t.\t+\t.\t{attrs}"]}


def _genome_end(tmp_path, monkeypatch):
    """Two blocks searched as the search does (block 0 significant, block 1 not), block 2
    cut out but ended before its candidates were judged; GOI calls judged to the end."""
    monkeypatch.setattr(isr.synvoy_evalue, "blastp_evalue",
                        lambda q, subjects, space: (E_OF[subjects[0]], 50.0))
    args = _args(tmp_path)
    blocks = {0: ({"chrom": "c1", "start": 100_000, "end": 300_000, "genes": ["g2", "g3"]},
                  (50_000, 350_000), 7),
              1: ({"chrom": "c2", "start": 1_000_000, "end": 1_050_000, "genes": ["GOI_q"]},
                  (950_000, 1_100_000), 2),
              2: ({"chrom": "c3", "start": 10_000, "end": 20_000, "genes": ["g5", "g6"]},
                  (0, 70_000), 2)}
    terms = {}
    for b, (block, (h0, h1), k) in blocks.items():
        isr.record_block_row(args, b, block, h0, h1, k + 1, k)
        if b < 2:
            terms[b] = isr.block_place_terms(args, b, h0, h1, h0, h1, k, [], hood_len=h1 - h0)
    terms["gw"] = isr.block_place_terms(args, 5000, 5_000_000, 5_010_000, 5_000_000, 5_010_000,
                                        0, [], genome_wide=True)
    args.genome_ctx["block_context"] = {0: "", 1: "right", 2: "left"}
    cands = {0: [_cand(f"GOI_q|{TAG}_b0_l1_exon_ann", "c1", 200_001, 201_500, "A"),
                 _cand(f"GOI_q|{TAG}_b0_l2_exon_ann", "c1", 250_001, 250_900, "B", bits=40.0)],
             1: [_cand(f"GOI_q|{TAG}_b1_l1_exon_ann", "c2", 1_020_001, 1_020_600, "C")],
             "gw": [_cand(f"GOI_q|{TAG}_b5000_l1_exon_ann", "c9", 5_005_001, 5_005_900, "W")]}
    lines = []
    for key, cs in cands.items():
        isr.annotate_place_evalues(cs, {"GOI_q": {"seq": "M" * 120}}, {}, terms[key], "", 0)
        for c in cs:
            mid = isr._parse_gff_attr_string(c["gff"][0].split("\t")[8])["ID"]
            args.genome_ctx.setdefault("judge_numbers", {})[mid] = c["judge_numbers"]
            lines.extend(c["gff"])
    lines.append(f"c1\tflanking_miniprot\tmRNA\t120001\t121000\t.\t+\t.\tID=g2|{TAG}_b0_fl1_flank_ann;"
                 "SynVoyRole=flanking;Confidence=HIGH")
    lines = isr.finalize_place_evalues(lines, args.genome_ctx)
    lines, _counts = isr.relabel_goi_lines(lines, args.genome_ctx["judge_numbers"])
    isr._write_block_table(args, "locus_1_", GENOME, lines, args.genome_ctx,
                           lambda b: any(g.startswith("GOI_") for g in b["genes"]))
    path = tmp_path / "blocks" / f"locus_1_{GENOME}.blocks.tsv"
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {int(r["block"]): r for r in reader}
        header = reader.fieldnames
    return terms, lines, rows, header


def test_one_row_per_searched_block_with_its_esyn(tmp_path, monkeypatch):
    terms, _lines, rows, header = _genome_end(tmp_path, monkeypatch)
    assert header == isr.BLOCK_TABLE_COLUMNS
    assert sorted(rows) == [0, 1, 2]                      # no row for the genome-wide window
    sig, weak, ended = rows[0], rows[1], rows[2]
    assert terms[0]["significant"] and not terms[1]["significant"]
    assert float(sig["esyn"]) == terms[0]["esyn"] and sig["esyn"] == repr(terms[0]["esyn"])
    assert float(weak["esyn"]) == terms[1]["esyn"]
    assert (sig["significant"], weak["significant"]) == ("1", "0")
    assert (ended["esyn"], ended["significant"], ended["best_call"]) == ("inf", "0", "")
    assert (sig["chrom"], sig["start"], sig["end"], sig["hood_start"], sig["hood_end"]) == (
        "c1", "100000", "300000", "50000", "350000")
    assert (sig["flanking_support"], sig["collinear_support"]) == ("8", "7")
    assert (sig["goi_seeded"], weak["goi_seeded"]) == ("0", "1")
    assert (sig["context"], weak["context"], ended["context"]) == ("", "right", "left")


def test_counts_and_best_call_are_the_final_labels(tmp_path, monkeypatch):
    _terms, _lines, rows, _header = _genome_end(tmp_path, monkeypatch)
    sig, weak = rows[0], rows[1]
    assert (sig["best_call"], sig["best_call_id"]) == ("HIGH", f"GOI_q|{TAG}_b0_l1_exon_ann")
    assert [sig[c] for c in ("n_high", "n_medium", "n_ambiguous", "n_low")] == ["1", "1", "0", "0"]
    assert (weak["best_call"], weak["n_low"], weak["n_high"]) == ("LOW", "1", "0")


def test_block_esyn_of_every_goi_row_matches_its_block_row(tmp_path, monkeypatch):
    _terms, lines, rows, _header = _genome_end(tmp_path, monkeypatch)
    checked = 0
    for line in lines:
        cols = line.split("\t")
        a = isr._parse_gff_attr_string(cols[8])
        if a.get("SynVoyRole") != "goi":
            continue
        b = isr.block_index_of(a["ID"], TAG)
        if b not in rows:
            assert b == 5000 and a["BlockEsyn"] == "inf"   # the genome-wide window
            continue
        assert a["BlockEsyn"] == f"{float(rows[b]['esyn']):.3g}"
        checked += 1
    assert checked == 3


def test_header_only_where_the_search_reached_a_genome_without_a_block(tmp_path):
    isr._write_block_table(SimpleNamespace(output_dir=str(tmp_path)), "", "G.fna")
    with open(tmp_path / "blocks" / "G.fna.blocks.tsv") as fh:
        assert fh.read().rstrip("\n").split("\t") == isr.BLOCK_TABLE_COLUMNS


def test_block_index_of_reads_the_block_after_the_genome_tag():
    assert isr.block_index_of(f"GOI_q|{TAG}_b13_l1_exon_ann", TAG) == 13
    assert isr.block_index_of(f"GOI_copy_3|{TAG}_b4_l1", TAG) == 4
    assert isr.block_index_of(f"GOI_q|{TAG}_b1_r1_fallback", TAG) == 1
    assert isr.block_index_of("GOI_q|Other_fna_b1_l1_fallback", TAG) is None
    assert isr.block_index_of(f"gene-X|{TAG}_rearranged", TAG) is None
    assert isr.block_index_of("", TAG) is None


def test_the_search_records_and_writes_the_table():
    src = open(os.path.join(os.path.dirname(__file__), "..", "bin", "iterative_search_runner.py")).read()
    assert "record_block_row(args, block_idx, block, hood_start, hood_end," in src
    assert "_write_block_table(args, prefix, genome_name, all_gff_lines, args.genome_ctx, _goi_seeded)" in src
    module = open(os.path.join(os.path.dirname(__file__), "..", "modules", "iterative_search.nf")).read()
    assert "pattern: 'iterative_results/blocks/*'" in module
    # Nextflow publishes a file inside a declared output folder only with that folder (the
    # table sat in hits/ first and was never published): no output folder may hold it.
    outputs = re.findall(r'path\("(iterative_results/[^"*]+)"\)', module)
    assert outputs and not any("iterative_results/blocks/".startswith(o.rstrip('/') + '/')
                               for o in outputs)
    assert synvoy_evalue.SIGNIFICANT_BLOCK_ESYN == 0.1
