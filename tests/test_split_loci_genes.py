"""Home loci in natural units (bin/split_loci.py, 2026-10-08/09).

* Two places where the query aligns are one locus when few annotated GENES lie between
  them; the fixed 50 kb is only the fall-back without an annotation.
* The annotated home gene that the query is becomes the primary locus.
* A locus is scored as ONE copy of the gene (whole gene, not best exon, not a sum over
  copies), and the primary locus cannot be pushed out by the cap.
"""
import json
import os
import subprocess
import sys

import pytest

BIN = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "bin"))
sys.path.insert(0, BIN)

import split_loci as sl  # noqa: E402

SCRIPT = os.path.join(BIN, "split_loci.py")


def _entry(chrom, start, end, evalue=1e-50, bit=200.0):
    return {"chrom": chrom, "start": start, "end": end, "evalue": evalue, "bit": bit,
            "line": f"{chrom}\t{start}\t{end}\thit\t{evalue}\t+\t{bit:g}"}


def _spans(clusters):
    return [(c[0]["chrom"], min(e["start"] for e in c), max(e["end"] for e in c)) for c in clusters]


def _gff(path, genes, biotype=True):
    """genes: (chrom, start_1based, end, name[, biotype])."""
    lines = ["##gff-version 3"]
    for g in genes:
        chrom, start, end, name = g[:4]
        bt = g[4] if len(g) > 4 else "protein_coding"
        attrs = f"ID=gene-{name};Name={name}" + (f";gene_biotype={bt}" if biotype else "")
        lines.append(f"{chrom}\tRefSeq\tgene\t{start}\t{end}\t.\t+\t.\t{attrs}")
        lines.append(f"{chrom}\tRefSeq\tmRNA\t{start}\t{end}\t.\t+\t.\tID=rna-{name};Parent=gene-{name}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# ------------------------------------------------------------------ genes as the unit

def test_gene_spans_are_sorted_zero_based_and_protein_coding(tmp_path):
    gff = _gff(tmp_path / "h.gff", [("chr1", 5001, 6000, "B"), ("chr1", 1001, 2000, "A"),
                                    ("chr1", 3001, 3500, "NC", "lncRNA"), ("chr2", 11, 90, "C")])
    spans = sl.load_gene_spans(gff, {"chr1"})
    assert spans == {"chr1": [(1000, 2000), (5000, 6000)]}


def test_every_gene_counts_when_the_annotation_states_no_biotype(tmp_path):
    gff = _gff(tmp_path / "h.gff", [("chr1", 1001, 2000, "A"), ("chr1", 3001, 3500, "B")], biotype=False)
    assert sl.load_gene_spans(gff, {"chr1"}) == {"chr1": [(1000, 2000), (3000, 3500)]}


@pytest.mark.parametrize("path", [None, "NO_GFF", "/nonexistent/NO_GFF", "/nonexistent/home.gff"])
def test_no_annotation_gives_no_spans(path):
    assert sl.load_gene_spans(path, {"chr1"}) == {}


def test_genes_between_matches_the_shared_definition():
    import synvoy_scale
    spans = sorted([(100, 200), (300, 400), (450, 700), (800, 900), (950, 990)])
    for left, right in [(200, 800), (150, 800), (200, 650), (0, 2000), (900, 950), (990, 5000)]:
        assert sl.genes_between(spans, left, right) == synvoy_scale.genes_between(spans, left, right)


def test_hits_with_few_genes_between_are_one_locus_whatever_the_distance():
    """Gene-poor region: two members of a family 400 kb apart with two genes between."""
    spans = {"chr1": [(0, 900), (150_000, 160_000), (300_000, 310_000), (520_000, 530_000)]}
    entries = [_entry("chr1", 1_000, 2_000), _entry("chr1", 400_000, 401_000)]
    assert len(sl.split_into_loci(entries, spans, gap_genes=5, gap_bp=None)) == 1
    assert len(sl.split_into_loci(entries, spans, gap_genes=1, gap_bp=None)) == 2
    assert len(sl.split_into_loci(entries, spans, gap_genes=5, gap_bp=50_000)) == 2   # bp wins


def test_hits_with_many_genes_between_are_two_loci_even_when_close():
    """Compact genome: 30 kb hold a dozen genes, so two hits 30 kb apart are far apart."""
    spans = {"chr1": [(2_500 + i * 2_400, 4_500 + i * 2_400) for i in range(12)]}
    entries = [_entry("chr1", 1_000, 2_000), _entry("chr1", 32_000, 33_000)]
    assert len(sl.split_into_loci(entries, spans, gap_genes=5, gap_bp=None)) == 2
    assert len(sl.split_into_loci(entries, {}, gap_genes=5, gap_bp=None)) == 1        # 50 kb fall-back


def test_a_sequence_without_annotated_genes_falls_back_to_50_kb():
    spans = {"chr1": [(10, 20)]}
    entries = [_entry("scaffold9", 1_000, 2_000), _entry("scaffold9", 40_000, 41_000),
               _entry("scaffold9", 120_000, 121_000)]
    assert _spans(sl.split_into_loci(entries, spans, 5, None)) == [
        ("scaffold9", 1_000, 41_000), ("scaffold9", 120_000, 121_000)]


def test_the_gap_is_measured_from_the_furthest_end_of_the_locus():
    """A long first interval that contains the second: the third is judged from its end."""
    entries = [_entry("chr1", 1_000, 90_000), _entry("chr1", 2_000, 3_000), _entry("chr1", 100_000, 101_000)]
    assert len(sl.split_into_loci(entries, {}, 5, None)) == 1


def test_loci_never_span_two_sequences():
    entries = [_entry("chr1", 1_000, 2_000), _entry("chr2", 1_500, 2_500)]
    assert len(sl.split_into_loci(entries, {}, 5, 10**9)) == 2


# ----------------------------------------------------------------- the home gene

def _home_files(tmp_path, rows, method="gff_annotation"):
    bed = tmp_path / "goi_annotation.bed"
    bed.write_text("".join(f"{c}\t{s}\t{e}\texon_{i}\t0\t{strand}\n"
                           for i, (c, s, e, strand) in enumerate(rows, 1)))
    info = tmp_path / "goi_info.json"
    info.write_text(json.dumps({"method": method, "num_exons": len(rows)}))
    return str(bed), str(info)


def test_home_gene_is_the_span_of_its_coding_exons(tmp_path):
    bed, info = _home_files(tmp_path, [("chr1", 5_000, 5_300, "-"), ("chr1", 9_000, 9_400, "-")])
    assert sl.load_home_gene(bed, info) == ("chr1", 5_000, 9_400, "-")


def test_a_home_model_built_by_alignment_is_not_an_annotated_gene(tmp_path):
    bed, info = _home_files(tmp_path, [("chr1", 5_000, 5_300, "+")], method="miniprot")
    assert sl.load_home_gene(bed, info) is None
    assert sl.load_home_gene(str(tmp_path / "absent.bed"), info) is None
    (tmp_path / "empty.bed").write_text("")
    assert sl.load_home_gene(str(tmp_path / "empty.bed"), None) is None


def test_home_gene_replaces_the_hit_clusters_it_overlaps():
    gene = ("chr1", 10_000, 20_000, "-")
    clusters = [[_entry("chr1", 10_500, 11_000, 1e-80, 300.0), _entry("chr1", 19_000, 19_500, 1e-20, 90.0),
                 _entry("chr1", 24_000, 25_000, 1e-10, 60.0)],            # exons + a neighbour's hit
                [_entry("chr4", 500, 900, 1e-30, 120.0)]]
    primary, others = sl.apply_home_gene(clusters, gene)
    assert _spans([primary]) == [("chr1", 10_000, 20_000)]
    assert primary[0]["evalue"] == 1e-80 and primary[0]["bit"] == 300.0
    assert primary[0]["line"].split("\t") == ["chr1", "10000", "20000", "gene_loc", "1e-80", "-", "300"]
    assert _spans(others) == [("chr4", 500, 900)]


def test_home_gene_without_any_hit_still_becomes_the_primary_locus():
    primary, others = sl.apply_home_gene([[_entry("chr4", 500, 900)]], ("chr1", 10_000, 20_000, "+"))
    assert _spans([primary]) == [("chr1", 10_000, 20_000)] and primary[0]["bit"] is None
    assert primary[0]["line"].endswith("\t+\t.") and len(others) == 1


# ---------------------------------------------------------------- score of a locus

def _m8(path, rows):
    """rows: (chrom, qstart, qend, sstart, send, evalue, bits)."""
    path.write_text("".join(
        f"q\t{c}\t50.0\t{qe - qs + 1}\t0\t0\t{qs}\t{qe}\t{ss}\t{se}\t{ev}\t{bits}\n"
        for c, qs, qe, ss, se, ev, bits in rows))
    return str(path)


def test_raw_hits_keep_one_hit_per_stretch_of_dna_and_respect_the_evalue(tmp_path):
    path = _m8(tmp_path / "h.m8", [
        ("chr1", 1, 100, 1001, 1300, 1e-30, 150.0),
        ("chr1", 1, 100, 1300, 1010, 1e-28, 140.0),     # the other program, minus strand order
        ("chr1", 101, 200, 5001, 5300, 1e-20, 90.0),
        ("chr1", 20, 60, 9001, 9120, 5.0, 25.0),         # above the E-value limit
    ])
    hits = sl.load_raw_hits([path, str(tmp_path / "absent.m8"), ""], 0.01)
    assert sorted((h["start"], h["end"], h["bits"]) for h in hits) == [(1000, 1300, 150.0), (5000, 5300, 90.0)]
    assert {(h["qstart"], h["qend"]) for h in hits} == {(1, 100), (101, 200)}


def test_locus_score_is_the_whole_gene_not_the_best_exon(tmp_path):
    """The gene itself: three exons of 80 bits. A relative: one shared domain of 150 bits.
    By best exon the relative outranks the gene; as whole genes it does not."""
    hits = sl.load_raw_hits([_m8(tmp_path / "h.m8", [
        ("chr1", 1, 100, 1001, 1300, 1e-20, 80.0), ("chr1", 101, 200, 2001, 2300, 1e-20, 80.0),
        ("chr1", 201, 300, 3001, 3300, 1e-20, 80.0),
        ("chr7", 101, 200, 501, 800, 1e-40, 150.0),
    ])], 0.01)
    gene = [_entry("chr1", 1_000, 3_300, 1e-20, 80.0)]
    relative = [_entry("chr7", 500, 800, 1e-40, 150.0)]
    sl.score_clusters_by_hits([gene, relative], hits)
    assert sl._cluster_best_bit(gene) == 240.0 and sl._cluster_best_bit(relative) == 150.0


def test_a_tandem_array_is_scored_as_one_copy(tmp_path):
    hits = sl.load_raw_hits([_m8(tmp_path / "h.m8", [
        ("chr1", 1, 300, 1001, 1900, 1e-90, 400.0), ("chr1", 1, 300, 3001, 3900, 1e-80, 380.0),
        ("chr1", 1, 300, 5001, 5900, 1e-70, 350.0),
    ])], 0.01)
    array = [_entry("chr1", 1_000, 5_900)]
    sl.score_clusters_by_hits([array], hits)
    assert sl._cluster_best_bit(array) == 400.0


def test_without_raw_hits_the_interval_score_is_used():
    cl = [_entry("chr1", 1, 9, bit=77.0), _entry("chr1", 20, 30, bit=91.0)]
    sl.score_clusters_by_hits([cl], [])
    assert sl._cluster_best_bit(cl) == 91.0


# ---------------------------------------------------------------------- end to end

def _run(tmp_path, bed_rows, *extra):
    bed = tmp_path / "gene_loc.bed"
    bed.write_text("".join(f"{c}\t{s}\t{e}\tgene_loc\t{ev}\t+\t{bit}\n" for c, s, e, ev, bit in bed_rows))
    out_dir = tmp_path / f"run{len(list(tmp_path.glob('run*')))}"      # a fresh folder per call
    out_dir.mkdir()
    prefix = out_dir / "locus"
    r = subprocess.run([sys.executable, SCRIPT, "--bed", str(bed), "--output_prefix", str(prefix),
                        *extra], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    loci = []
    i = 1
    while os.path.exists(f"{prefix}_{i}.bed"):
        with open(f"{prefix}_{i}.bed") as fh:
            loci.append([line.rstrip("\n").split("\t") for line in fh])
        i += 1
    return loci, r.stdout, r.stderr


BED = [("chr1", 10_500, 11_000, 1e-80, 300), ("chr1", 19_000, 19_500, 1e-20, 90),
       ("chr1", 24_000, 25_000, 1e-60, 250), ("chr4", 500, 900, 1e-70, 280)]


def test_cli_default_without_inputs_is_the_former_behaviour(tmp_path):
    loci, out, _err = _run(tmp_path, BED)
    assert "Locus gap: 50000 bp (no annotation" in out
    assert [[r[:3] for r in locus] for locus in loci] == [
        [["chr1", "10500", "11000"], ["chr1", "24000", "25000"], ["chr1", "19000", "19500"]],
        [["chr4", "500", "900"]]]


def test_cli_annotated_home_gene_is_the_primary_locus(tmp_path):
    bed, info = _home_files(tmp_path, [("chr1", 10_400, 11_000, "-"), ("chr1", 19_000, 19_600, "-")])
    loci, out, _err = _run(tmp_path, BED, "--home_gene_bed", bed, "--home_gene_info", info)
    assert "Primary locus = the annotated home gene chr1:10400-19600 (-); 1 further place" in out
    assert loci[0] == [["chr1", "10400", "19600", "gene_loc", "1e-80", "-", "300"]]
    assert [r[0] for r in loci[1]] == ["chr4"]


def test_cli_gene_only_keeps_the_annotated_locus_alone(tmp_path):
    bed, info = _home_files(tmp_path, [("chr1", 10_400, 19_600, "+")])
    loci, _out, _err = _run(tmp_path, BED, "--home_gene_bed", bed, "--home_gene_info", info,
                            "--home_locus_mode", "gene_only")
    assert len(loci) == 1 and loci[0][0][:3] == ["chr1", "10400", "19600"]


def test_cli_gene_only_without_an_annotated_gene_warns_and_uses_the_alignments(tmp_path):
    loci, _out, err = _run(tmp_path, BED, "--home_locus_mode", "gene_only")
    assert "no annotated home gene matched" in err and len(loci) == 2


def test_cli_alignment_mode_ignores_the_home_gene(tmp_path):
    bed, info = _home_files(tmp_path, [("chr1", 10_400, 19_600, "+")])
    loci, out, _err = _run(tmp_path, BED, "--home_gene_bed", bed, "--home_gene_info", info,
                           "--home_locus_mode", "alignment")
    assert "Primary locus" not in out and len(loci[0]) == 3


def test_cli_counts_genes_when_an_annotation_is_given(tmp_path):
    """24 kb to the neighbouring hit, but eight genes lie between: a locus of its own."""
    gff = _gff(tmp_path / "h.gff", [("chr1", 19_601 + i * 500, 19_900 + i * 500, f"G{i}") for i in range(8)])
    loci, out, _err = _run(tmp_path, BED[:3], "--gff", gff)
    assert "more than 5 annotated gene(s)" in out
    assert [len(locus) for locus in loci] == [2, 1]
    loci, out, _err = _run(tmp_path, BED[:3], "--gff", gff, "--locus_gap", "50000")
    assert "50000 bp (set by the user)" in out and [len(locus) for locus in loci] == [3]


def test_cli_relative_with_one_strong_domain_is_dropped_by_the_whole_gene_score(tmp_path):
    hits = _m8(tmp_path / "hits.m8", [
        ("chr1", 1, 100, 10_501, 11_000, 1e-20, 80.0), ("chr1", 101, 200, 14_001, 14_500, 1e-20, 80.0),
        ("chr1", 201, 300, 19_001, 19_500, 1e-20, 80.0), ("chr1", 301, 400, 19_701, 20_000, 1e-20, 80.0),
        ("chr4", 101, 200, 501, 900, 1e-40, 150.0),
    ])
    bed_rows = [("chr1", 10_500, 20_000, 1e-20, 80), ("chr4", 500, 900, 1e-40, 150)]
    loci, out, _err = _run(tmp_path, bed_rows, "--hits", hits)
    assert len(loci) == 1 and loci[0][0][0] == "chr1"
    assert "Filtered out locus chr4:500-900" in out and "0.47x primary 320" in out

    loci_old, _out, _err = _run(tmp_path, bed_rows)                 # ranked by best exon, as before
    assert [locus[0][0] for locus in loci_old] == ["chr4", "chr1"]


def test_cli_cap_never_removes_the_primary_locus(tmp_path):
    """Six relatives that each score higher than the annotated gene: the gene stays."""
    bed, info = _home_files(tmp_path, [("chr1", 10_400, 19_600, "+")])
    rows = [("chr1", 10_500, 11_000, 1e-20, 90)] + [
        (f"chr{i}", 500, 900, f"1e-{40 + i}", 200 + i) for i in range(2, 8)]
    loci, out, _err = _run(tmp_path, rows, "--home_gene_bed", bed, "--home_gene_info", info,
                           "--max_loci", "3")
    assert len(loci) == 3 and loci[0][0][:3] == ["chr1", "10400", "19600"]
    assert sorted(locus[0][0] for locus in loci[1:]) == ["chr6", "chr7"]      # the two strongest
    assert out.count("exceeds --max_loci 3") == 4


# ------------------------------------------- the preset chooser counts loci the same way

import profile_hits as ph  # noqa: E402

PROFILE = os.path.join(BIN, "profile_hits.py")


def _phit(chrom, tstart, tend, qstart=1, qend=100, bits=100.0):
    return {"chrom": chrom, "tstart": tstart, "tend": tend, "qstart": qstart, "qend": qend,
            "bits": bits, "evalue": 1e-30, "pident": 60.0}


def test_profile_and_locus_step_split_the_same_hits_the_same_way():
    """One definition of a home locus: the hit profile (which chooses the preset) and
    split_loci.py must not disagree about how many places the query aligns at."""
    cases = [
        # gene-poor: 400 kb apart, two genes between
        ({"chr1": [(0, 900), (150_000, 160_000), (300_000, 310_000), (520_000, 530_000)]},
         [(1_000, 2_000), (400_000, 401_000)]),
        # compact: 30 kb apart, twelve genes between
        ({"chr1": [(2_500 + i * 2_400, 4_500 + i * 2_400) for i in range(12)]},
         [(1_000, 2_000), (32_000, 33_000)]),
        # no genes annotated on this sequence: 50 kb decides
        ({"chrOther": [(10, 20)]}, [(1_000, 2_000), (40_000, 41_000), (120_000, 121_000)]),
    ]
    for spans, places in cases:
        by_split = sl.split_into_loci([_entry("chr1", s, e) for s, e in places], spans, 5, None)
        by_profile = ph._cluster_loci([_phit("chr1", s, e) for s, e in places],
                                      sl.FALLBACK_GAP_BP, spans, 5)
        assert len(by_profile) == len(by_split)
    assert [len(sl.split_into_loci([_entry("chr1", s, e) for s, e in c[1]], c[0], 5, None))
            for c in cases] == [1, 2, 2]


def test_profile_tandem_pair_in_a_gene_poor_region_is_one_locus_with_two_copies():
    """Two copies 150 kb apart with one gene between: at the former fixed 100 kb these
    were two strong loci ("relatives elsewhere"); counted in genes it is one locus that
    holds the query twice."""
    hits = [_phit("chr1", 1_000, 1_900, 1, 300, 400.0), _phit("chr1", 151_000, 151_900, 1, 300, 380.0)]
    spans = {"chr1": [(60_000, 90_000)]}
    fixed = ph.build_profile(hits, 300, 100_000, 0.5)
    by_genes = ph.build_profile(hits, 300, sl.FALLBACK_GAP_BP, 0.5, gene_spans=spans, gap_genes=5)
    assert (fixed["n_independent_loci"], fixed["n_strong_loci"], fixed["primary_copies"]) == (2, 2, 1.0)
    assert (by_genes["n_independent_loci"], by_genes["n_strong_loci"], by_genes["primary_copies"]) == (1, 1, 2.0)
    assert by_genes["top_bit"] == 400.0                       # still the score of ONE copy


def _profile_cli(tmp_path, *extra):
    hits = tmp_path / "hits.m8"
    hits.write_text("q\tchr1\t60.0\t300\t0\t0\t1\t300\t1001\t1900\t1e-90\t400\n"
                    "q\tchr1\t58.0\t300\t0\t0\t1\t300\t151001\t151900\t1e-80\t380\n")
    (tmp_path / "q.faa").write_text(">q\n" + "M" * 300 + "\n")
    out = tmp_path / "profile.json"
    r = subprocess.run([sys.executable, PROFILE, "--mmseqs", str(hits), "--query", str(tmp_path / "q.faa"),
                        "--output", str(out), *extra], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.loads(out.read_text())


def test_profile_cli_counts_genes_with_an_annotation(tmp_path):
    gff = _gff(tmp_path / "h.gff", [("chr1", 60_001, 90_000, "G1")])
    p = _profile_cli(tmp_path, "--gff", gff)
    assert (p["locus_unit"], p["n_independent_loci"], p["primary_copies"]) == ("genes", 1, 2.0)


def test_profile_cli_without_an_annotation_uses_the_same_50_kb_as_the_locus_step(tmp_path):
    for extra in ((), ("--gff", "NO_GFF"), ("--gff", str(tmp_path / "absent.gff"))):
        p = _profile_cli(tmp_path, *extra)
        assert (p["locus_unit"], p["n_independent_loci"]) == ("bp (no annotation)", 2)


def test_profile_cli_gap_in_bp_wins(tmp_path):
    gff = _gff(tmp_path / "h.gff", [("chr1", 20_001 + i * 1_000, 20_900 + i * 1_000, f"G{i}") for i in range(30)])
    assert _profile_cli(tmp_path, "--gff", gff)["n_independent_loci"] == 2          # 30 genes between
    p = _profile_cli(tmp_path, "--gff", gff, "--locus_gap", "200000")
    assert (p["locus_unit"], p["n_independent_loci"]) == ("bp", 1)
