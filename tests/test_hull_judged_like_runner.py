"""Review fix RF-1d (batch of 2026-10-15): the hull rescue is judged like the runner.

* the label protein goes through synvoy_label.frame_guard (FrameFixed / FrameBroken; a
  broken frame is judged on the translation's pieces between stops);
* completeness is the query span of the miniprot alignment (Target=), not CDS aa / query
  length, which reaches 1.0 for any long model;
* the specificity gets the GOI's home neighbourhood (it got None);
* the 40 % identity gate stays as a cost filter in front of the judge (decision, 15 Oct).
* the place is the main path's (older item 51): the core or window of a significant block of the
  genome's block table that holds the model, its E-value there times SigBlocks; else the genome.
  It was always the genome, ~G/W times the E of the same model in a block's core, without the
  positional tie-break: a hull twin read AMBIGUOUS where the main path held MEDIUM (item 41).
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

import rescue_goi_hull as hull  # noqa: E402
import synvoy_scale  # noqa: E402

CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'}


def _dna(p):
    return "".join(CODON[a] for a in p)


def _rows(spans, strand="+"):
    return [["w", "miniprot", "CDS", str(s), str(e), ".", strand, "0", "Parent=m"] for s, e in spans]


def test_coverage_is_the_query_span_of_the_alignment():
    mrna = ["w", "miniprot", "mRNA", "1", "3000", "100", "+", ".", "ID=m;Identity=0.9;Target=q 11 59"]
    rows = _rows([(1, 1500), (2001, 3000)])                   # 833 aa of CDS
    assert hull._query_span_coverage(mrna, rows, 100) == 0.49
    assert hull._cds_coverage(rows, 100) == 1.0              # the old measure
    no_target = mrna[:8] + ["ID=m;Identity=0.9"]
    assert hull._query_span_coverage(no_target, rows, 100) == 1.0


def test_the_label_protein_goes_through_the_frame_guard():
    good = "AAA" + _dna("MKTAYIAKQRQISF") + "TAA"
    assert hull.hull_label_protein(good, _rows([(4, 3 + 42)])) == ("MKTAYIAKQRQISF", "")
    cds = "ATGTAGTTGAATTACTCTAAATGA"             # a stop in every frame (table 1)
    protein, flag = hull.hull_label_protein("AAA" + cds + "GGG", _rows([(4, 3 + len(cds))]))
    assert flag == "FrameBroken=1"
    assert protein == hull.translate(cds).rstrip("*") and "*" in protein   # judged in pieces


def test_the_flanking_anchors_of_the_genome_give_the_neighbourhood(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text("##gff-version 3\n"
                   "c1\tx\tmRNA\t100\t900\t.\t+\t.\tID=a|g_b0;SynVoy_Parent=gene-A;SynVoyRole=flanking;Confidence=HIGH\n"
                   "c1\tx\tmRNA\t2000\t2900\t.\t+\t.\tID=b|g_b0;SynVoy_Parent=gene-C;SynVoyRole=flanking;Confidence=LOW\n"
                   "c1\tx\tmRNA\t1000\t1900\t.\t+\t.\tID=GOI_q|g_b0;SynVoy_Parent=GOI_q;SynVoyRole=goi;Confidence=HIGH\n")
    assert hull.flanking_parent_ids(str(gff)) == ["gene-A", "gene-C"]


def test_the_judge_gets_the_neighbourhood_and_the_span_coverage(monkeypatch):
    seen = {}
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", lambda q, s, sp: (1e-30, 200.0))
    monkeypatch.setattr(hull, "sw_bits", lambda q, x: 200.0)

    def fake_spec(proteome, qid, q, x, bits, hood, n_flanking=0, own_gene_id=None):
        seen["hood"], seen["n_flanking"] = hood, n_flanking
        return {"margin": 50.0, "real_competitor": False, "best_other": "", "best_other_in_neighbourhood": None}
    monkeypatch.setattr(hull.synvoy_specificity, "specificity", fake_spec)
    monkeypatch.setattr(hull.synvoy_specificity, "gff_attrs", lambda spec: "SpecMargin=50")
    rows = _rows([(1, 300)])
    conf, _key, attrs, _bits = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, 10 ** 7, object(),
                                                     label_protein="M" * 100, coverage=0.49,
                                                     neighbourhood={"gene-A", "gene-B"}, n_flanking=10)
    assert seen["hood"] == {"gene-A", "gene-B"} and seen["n_flanking"] == 10   # RF-2
    assert conf == "MEDIUM"                                  # significant and specific, 0.49 < 0.7
    conf, *_ = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, 10 ** 7, object(),
                                     label_protein="M" * 100, coverage=0.9, neighbourhood=None)
    assert conf == "HIGH"


def test_the_hull_process_hands_over_the_home_database_folder():
    root = os.path.join(os.path.dirname(__file__), "..")
    module = open(os.path.join(root, "modules", "rescue_goi_hull.nf")).read()
    assert "--home_db_dir ${home_db}" in module
    assert "--goi_info ${goi_info}" in module                    # RF-2b: the own home gene
    sub = open(os.path.join(root, "subworkflows", "adjudicate_and_report.nf")).read()
    assert ("RESCUE_GOI_HULL(hull_inputs_ch, genomes_dir, rescue_query_ch, home_db, goi_info, "
            "effective_params)") in sub                   # profile stage 4: the run's settings
    assert "--effective_params ${effective_params}" in module
    main_nf = open(os.path.join(root, "main.nf")).read()
    assert ("RESOLVE_SETTINGS.out.json,\n            ANNOTATE_GOI.out.info.first(),\n"
            "            ANNOTATE_GOI.out.exons.first()") in main_nf   # item 17: the panel's q (inputs added since follow it)
    assert "home_db.map { d -> d.resolve('home_proteome.faa') }" in sub
    # Older item 51: the search's block table of the same (locus, genome) reaches the hull task.
    assert "ITERATIVE_SEARCH.out.hits,\n            ITERATIVE_SEARCH.out.blocks,\n" in main_nf
    assert "iter_blocks" in sub and ".join(hull_blocks_ch, remainder: true)" in sub
    assert "assets/sentinels/NO_BLOCKS" in sub
    assert os.path.exists(os.path.join(root, "assets", "sentinels", "NO_BLOCKS"))
    assert "path(region_gff), path(blocks_tsv)" in module and "--blocks_tsv ${blocks_tsv}" in module


# ---- the place (older item 51) ------------------------------------------------------------

ROWS = [  # blocks/<locus>_<genome>.blocks.tsv rows: span and neighbourhood 0-based, half-open
    {"chrom": "chr1", "start": "100000", "end": "200000", "hood_start": "50000", "hood_end": "250000",
     "significant": "1"},
    {"chrom": "chr1", "start": "300000", "end": "400000", "hood_start": "250000", "hood_end": "450000",
     "significant": "0"},
    {"chrom": "chr2", "start": "100000", "end": "200000", "hood_start": "0", "hood_end": "300000",
     "significant": "1"},
]
DISTANCES = {"goi_margin": 10000, "search_margin": 30000, "region_padding": 50000}
G = 10_000_000


def test_the_place_is_the_core_or_window_of_a_significant_block():
    assert hull.hull_place("chr1", (95001, 96000), ROWS, DISTANCES, G) == ("core", 120000, 2)
    assert hull.hull_place("chr1", (75001, 76000), ROWS, DISTANCES, G) == ("window", 160000, 2)
    # In the neighbourhood but outside the DNA the search cut out (block +- search margin).
    assert hull.hull_place("chr1", (60001, 61000), ROWS, DISTANCES, G) == ("genome", G, 2)
    # Without a search margin the search cut out the whole neighbourhood.
    no_margin = {"goi_margin": 10000}
    assert hull.hull_place("chr1", (60001, 61000), ROWS, no_margin, G) == ("window", 200000, 2)
    # A block that is not significant predicts no place; nor does one on another scaffold.
    assert hull.hull_place("chr1", (350001, 351000), ROWS, DISTANCES, G) == ("genome", G, 2)
    assert hull.hull_place("chr3", (95001, 96000), ROWS, DISTANCES, G) == ("genome", G, 2)
    # The core ends at the end of the scaffold, as in the search.
    assert hull.hull_place("chr1", (195001, 196000), ROWS, DISTANCES, G, chrom_len=205000) == \
        ("core", 115000, 2)


def test_without_a_block_table_or_distances_the_place_is_the_genome(tmp_path):
    assert hull.hull_place("chr1", (95001, 96000), [], DISTANCES, G) == ("genome", G, 0)
    assert hull.hull_place("chr1", (95001, 96000), ROWS, {}, G) == ("genome", G, 2)    # an older GFF
    assert hull.read_block_rows(None) == [] and hull.read_block_rows(str(tmp_path / "NO_BLOCKS")) == []
    (tmp_path / "NO_BLOCKS").write_text("sentinel: no block table\n")
    assert hull.read_block_rows(str(tmp_path / "NO_BLOCKS")) == []
    table = tmp_path / "locus_1_g.fna.blocks.tsv"
    cols = ["genome", "block", "chrom", "start", "end", "hood_start", "hood_end", "significant"]
    table.write_text("\t".join(cols) + "\n" + "g.fna\t0\tchr1\t100000\t200000\t50000\t250000\t1\n")
    rows = hull.read_block_rows(str(table))
    assert hull.hull_place("chr1", (95001, 96000), rows, DISTANCES, G) == ("core", 120000, 1)


def _fake_terms(monkeypatch, margin=50.0, real=False, near=None):
    """blastp E proportional to the search space (1e-15 per residue pair), bits 200, and fixed
    specificity terms."""
    spaces = []

    def blastp(q, subjects, searchsp):
        spaces.append(searchsp)
        return searchsp * 1e-15, 200.0
    monkeypatch.setattr(hull.synvoy_evalue, "blastp_evalue", blastp)
    monkeypatch.setattr(hull, "sw_bits", lambda q, x: 200.0)
    monkeypatch.setattr(hull.synvoy_specificity, "specificity", lambda *a, **k: {
        "margin": margin, "real_competitor": real, "best_other": "gene-P" if real else "",
        "in_neighbourhood": near})
    monkeypatch.setattr(hull.synvoy_specificity, "gff_attrs", lambda spec: f"SpecMargin={spec['margin']}")
    return spaces


def test_the_judge_takes_the_place_and_the_multiplier(monkeypatch):
    spaces = _fake_terms(monkeypatch)
    rows = _rows([(1, 300)])
    conf, key, attrs, _bits = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, G, object(),
                                                    label_protein="M" * 100, coverage=0.9,
                                                    place="core", place_bp=120000, sig_blocks=3)
    e_core = 100 * (120000 // 3) * 1e-15
    assert spaces[-1] == 100 * (120000 // 3)
    assert attrs.startswith("Basis=significant_specific;Place=core;")
    assert f"EPlace={e_core:.3g};SigBlocks=3;Eeff={3 * e_core:.3g}" in attrs
    assert conf == "HIGH" and key[1] == 3 * e_core
    conf, key, attrs, _bits = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, G, object(),
                                                    label_protein="M" * 100, coverage=0.9)
    e_genome = 100 * (G // 3) * 1e-15
    assert spaces[-1] == 100 * (G // 3)
    assert attrs.startswith(f"Basis=significant_specific;Place=genome;Eeff={e_genome:.3g};")  # as before
    assert "EPlace" not in attrs and "SigBlocks" not in attrs


def test_a_twin_now_reads_the_main_paths_label(monkeypatch):
    """Item 41's case: another home protein explains the model within 3.3 bits, and it lies
    outside the GOI's home neighbourhood. In the core of a significant block the position decides
    for the gene (MEDIUM, as on the main path); in the genome place it cannot (AMBIGUOUS)."""
    _fake_terms(monkeypatch, margin=1.0, real=True, near=False)
    rows = _rows([(1, 300)])
    genome, *_ = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, G, object(),
                                       label_protein="M" * 100, coverage=0.9)
    core, *_ = hull.judge_hull_model("q", "M" * 100, "M" * 100, rows, G, object(),
                                     label_protein="M" * 100, coverage=0.9,
                                     place="core", place_bp=120000, sig_blocks=1)
    assert (genome, core) == ("AMBIGUOUS", "MEDIUM")
    main_rows = {"chr1": [(95001, 96000, "+", "MEDIUM", "GOI_q|g_b0_l1")]}
    assert hull.main_path_twin("chr1", 95001, 96000, "+", core, main_rows) == "GOI_q|g_b0_l1"


def _main(tmp_path, monkeypatch, blocks=None):
    """rescue_goi_hull.main on one hull (four HIGH flanking models on chr1 at 1.0-1.61 Mb), with
    the window, miniprot and blastp replaced: one model near chr1:990 kb."""
    header = synvoy_scale.format_gff_header({"goi_margin": 20000, "search_margin": 60000})
    rows = ["##gff-version 3", header] + [
        f"chr1\tx\tmRNA\t{1_000_000 + i * 200_000}\t{1_010_000 + i * 200_000}\t.\t+\t.\t"
        f"ID=f{i};SynVoyRole=flanking;Confidence=HIGH" for i in range(4)]
    (tmp_path / "region.gff").write_text("\n".join(rows) + "\n")
    (tmp_path / "q.faa").write_text(">GOI_X\n" + "M" * 100 + "\n")
    (tmp_path / "t.fna.fai").write_text("chr1\t5000000\t6\t60\t61\n")
    model = ("chr1\tminiprot\tmRNA\t100\t279\t50\t+\t.\tID=MP1;Identity=0.9000;Target=GOI_X 21 80\n"
             "chr1\tminiprot\tCDS\t100\t279\t50\t+\t0\tParent=MP1;Identity=0.9000;Target=GOI_X 21 80\n")
    monkeypatch.setattr(hull, "_extract_window", lambda *a: "ATGGCC" * 1000)
    monkeypatch.setattr(hull, "_run_miniprot_relaxed", lambda *a: model)
    _fake_terms(monkeypatch)
    extra = ["--blocks_tsv", str(blocks)] if blocks else []
    monkeypatch.setattr(sys, "argv", [
        "rescue_goi_hull.py", "--region_gff", str(tmp_path / "region.gff"),
        "--target_genome", str(tmp_path / "t.fna"), "--query", str(tmp_path / "q.faa"),
        "--genome_name", "t", "--output", str(tmp_path / "out.gff"),
        "--cluster_max_gap", "3000000", "--window_pad", "10000", "--max_window", "20000000",
        "--min_coverage", "0.5", *extra])
    assert hull.main() == 0
    return [ln for ln in (tmp_path / "out.gff").read_text().splitlines() if "\tgene\t" in ln]


def test_the_rescue_writes_the_place_of_the_block_table(tmp_path, monkeypatch):
    table = tmp_path / "locus_1_t.fna.blocks.tsv"
    cols = ["genome", "block", "chrom", "start", "end", "hood_start", "hood_end", "significant"]
    table.write_text("\t".join(cols) + "\n"
                     "t.fna\t0\tchr1\t1000000\t1610000\t900000\t1710000\t1\n"
                     "t.fna\t1\tchr1\t3000000\t3100000\t2900000\t3200000\t1\n")
    gene = _main(tmp_path, monkeypatch, table)
    start, end = (int(x) for x in gene[0].split("\t")[3:5])
    assert len(gene) == 1 and 980000 <= start < end <= 1630000      # in block 0's core
    core_bp = 1630000 - 980000
    assert "Place=core;" in gene[0] and "SigBlocks=2;" in gene[0]
    assert f"EPlace={100 * (core_bp // 3) * 1e-15:.3g};" in gene[0]
    gene = _main(tmp_path, monkeypatch)                                      # no table: the genome
    assert "Place=genome;" in gene[0] and "EPlace" not in gene[0] and "SigBlocks" not in gene[0]
    (tmp_path / "NO_BLOCKS").write_text("sentinel\n")
    assert _main(tmp_path, monkeypatch, tmp_path / "NO_BLOCKS") == gene
