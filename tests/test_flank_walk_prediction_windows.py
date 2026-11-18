"""The flank walk over a home genome without a GFF (bin/extract_flanking_genes.py).

Such a home genome is annotated only in stretches: genes are predicted around each home hit of
the query (prodigal_on_regions.py, pred_flank_window on each side) and gene models of a target
are transferred into one window around the locus (borrow_annotations.py). The walk stepped by
index through every gene of the scaffold, so once the stretch around the locus was used up it
took the genes predicted around another hit, however far away, as neighbours: their distance
inflated the gene pitch and every window of the search. The GFF's header now names the
stretches (##synvoy-annotated-window), and the walk keeps to the one that holds the locus.
All synthetic."""
import os
import random
import subprocess
import sys

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")
sys.path.insert(0, BIN)
import extract_flanking_genes as efg  # noqa: E402
import prodigal_on_regions as por  # noqa: E402
from sequence_utils import annotated_window_line, read_annotated_windows  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {'A': 'GCT', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT', 'G': 'GGT', 'H': 'CAT',
         'I': 'ATT', 'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA',
         'R': 'CGT', 'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAT'}
SCAFFOLD = 60_000
GOI_AT = 6001
# Window A = [0, 12000): three predicted genes on each side of the GOI. Window B = [40000, 56000):
# ten genes predicted around another home hit, 30 kb further on.
WINDOW_A = [1501, 3001, 4501, 7501, 9001, 10501]
WINDOW_B = [41001 + 1500 * k for k in range(10)]


def _protein(rng, n=100):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def _build(tmp_path, windows, extra=(), borrowed=()):
    """One scaffold with the GOI, window A's and window B's predicted genes (303 bp each), the
    genes at ``extra`` (predicted) and ``borrowed`` (transferred models); ``windows``: header
    directives (start, end, kind), 0-based half-open. Returns the gene id at each start."""
    rng = random.Random(5)
    goi = _protein(rng)
    seq = list("".join(rng.choice("ACGT") for _ in range(SCAFFOLD)))
    gff = ["##gff-version 3"] + [annotated_window_line("c", s, e, k) for s, e, k in windows]
    ids = {}
    starts = sorted([(GOI_AT, "goi")] + [(s, "pred") for s in WINDOW_A + WINDOW_B + list(extra)]
                    + [(s, "borrowed") for s in borrowed])
    for n, (s, kind) in enumerate(starts, 1):
        prot = goi if kind == "goi" else _protein(rng)
        cds = "".join(CODON[a] for a in prot) + "TAA"
        seq[s - 1:s - 1 + len(cds)] = cds
        e = s + len(cds) - 1
        gid = f"borrowed_t{n}" if kind == "borrowed" else f"c__{s // 1000}_g{n}"
        src = "borrowed_annotation" if kind == "borrowed" else "AUGUSTUS"
        gff += [f"c\t{src}\tgene\t{s}\t{e}\t.\t+\t.\tID={gid}",
                f"c\t{src}\tmRNA\t{s}\t{e}\t.\t+\t.\tID={gid}.t1;Parent={gid}",
                f"c\t{src}\tCDS\t{s}\t{e}\t.\t+\t0\tID={gid}.cds;Parent={gid}.t1"]
        ids[s] = gid
    (tmp_path / "g.fna").write_text(">c\n" + "".join(seq) + "\n")
    (tmp_path / "home.gff").write_text("\n".join(gff) + "\n")
    (tmp_path / "goi.faa").write_text(f">GOI_q\n{goi}\n")
    (tmp_path / "hits.tsv").write_text(f"q\tc\t100.0\t100\t0\t0\t1\t100\t{GOI_AT}\t{GOI_AT + 299}\t1e-60\t200\n")
    (tmp_path / "locus_1.bed").write_text(f"c\t{GOI_AT - 1}\t{GOI_AT + 302}\tbest_interval\t1e-60\t+\t200\n")
    return ids


def _walk(tmp_path, n_flank=10):
    p = subprocess.run([sys.executable, os.path.join(BIN, "extract_flanking_genes.py"),
                        "--bed", str(tmp_path / "locus_1.bed"), "--gff", str(tmp_path / "home.gff"),
                        "--genome", str(tmp_path / "g.fna"), "--n_flank", str(n_flank),
                        "--min_size", "0", "--exon_mode", "true", "--goi_faa", str(tmp_path / "goi.faa"),
                        "--goi_relative_evalue", "0.00001", "--home_hits", str(tmp_path / "hits.tsv"),
                        "--max_evalue", "0.01", "--out_bed", str(tmp_path / "out.bed"),
                        "--out_faa", str(tmp_path / "out.faa")],
                       capture_output=True, text=True, cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    rows = [ln.split("\t") for ln in open(tmp_path / "out.bed") if ln.strip()]
    anchors = [(r[3], int(r[1])) for r in rows if not r[3].startswith("GOI_")]
    return p.stdout, p.stderr, anchors


A_AND_B = [(0, 12000, "predicted"), (40000, 56000, "predicted")]


def test_the_walk_stays_in_the_prediction_window_of_the_locus(tmp_path):
    _build(tmp_path, A_AND_B)
    out, err, anchors = _walk(tmp_path)
    assert len(anchors) == 6 and all(start < 12000 for _gid, start in anchors)
    assert "downstream: 3 anchor(s)" in out and "ended: end of the annotated window" in out
    assert "the walk stays inside it" in out and "c:0-12000" in out
    assert "WARNING: the downstream walk reached the end of the stretch SynVoy annotated" in err
    assert "--pred_flank_window" in err
    assert "WARNING: the upstream walk" not in err          # upstream ends at the scaffold's start


def test_a_real_annotation_is_walked_as_before(tmp_path):
    _build(tmp_path, [])                                    # no directives: an annotation of the genome
    out, err, anchors = _walk(tmp_path)
    assert len([a for a in anchors if a[1] > GOI_AT]) == 10
    assert len([a for a in anchors if a[1] >= 40000]) == 7
    assert "end of the annotated window" not in out and "stays inside" not in out
    assert "WARNING" not in err


def test_a_borrowed_window_extends_the_stretch(tmp_path):
    ids = _build(tmp_path, A_AND_B + [(0, 30000, "borrowed")], borrowed=[20001])
    out, err, anchors = _walk(tmp_path)
    downstream = [gid for gid, start in anchors if start > GOI_AT]
    assert downstream == [ids[s] for s in (7501, 9001, 10501, 20001)]    # the 4th is the borrowed model
    assert all(start < 30000 for _gid, start in anchors)
    assert "c:0-30000" in out


def test_a_prediction_window_inside_the_transfer_window_is_not_walked_into(tmp_path):
    """borrow_annotations.py wrote the whole home window miniprot searched
    (it grows to about 800 kb) as 'borrowed', so the prediction window of another home hit
    inside it (B, 40-56 kb) joined the locus's stretch and the walk took B's genes past the
    unannotated DNA between. The merged GFF now names the stretch the transferred models
    cover, from the locus to the outermost model (20 kb), and the walk stops there."""
    import borrow_annotations as ba
    ids = _build(tmp_path, A_AND_B, borrowed=[20001])
    gff = tmp_path / "home.gff"
    pred = tmp_path / "pred.gff"
    pred.write_text("\n".join(ln for ln in gff.read_text().splitlines() if "borrowed_" not in ln) + "\n")
    model = {"source": "t", "name": "T", "chrom": "c", "strand": "+", "identity": 0.8, "coverage": 0.9,
             "protein": "M", "exons": [(20000, 20303, 0)]}
    locus = ("c", GOI_AT - 1, GOI_AT + 302)              # the best interval of the locus BED
    ba.write_merged_gff(str(gff), str(pred), set(), [model], "T", locus)
    assert read_annotated_windows(str(gff))["c"] == [
        (0, 12000, "predicted"), (GOI_AT - 1, 20303, "borrowed"), (40000, 56000, "predicted")]
    out, err, anchors = _walk(tmp_path)
    downstream = [gid for gid, start in anchors if start > GOI_AT]
    assert downstream == [ids[s] for s in (7501, 9001, 10501)] + ["borrowed_t"]
    assert not any(start >= 40000 for _gid, start in anchors)
    assert "c:0-20303" in out and "ended: end of the annotated window" in out
    assert "WARNING: the downstream walk reached the end of the stretch SynVoy annotated" in err


def test_touching_windows_are_one_stretch(tmp_path):
    ids = _build(tmp_path, [(0, 12000, "predicted"), (12000, 24000, "predicted"),
                            (40000, 56000, "predicted")], extra=[15001, 18001])
    out, _err, anchors = _walk(tmp_path)
    downstream = [gid for gid, start in anchors if start > GOI_AT]
    assert downstream == [ids[s] for s in (7501, 9001, 10501, 15001, 18001)]
    assert "c:0-24000" in out


def test_windows_that_miss_the_locus_set_no_limit(tmp_path):
    _build(tmp_path, [(40000, 56000, "predicted")])
    out, err, anchors = _walk(tmp_path)
    assert "none on c holds the locus" in err
    assert len([a for a in anchors if a[1] >= 40000]) == 7                # the walk as before


def test_walk_extent_merges_and_picks_the_stretch_of_the_locus():
    windows = {"c": [(0, 12000, "predicted"), (12000, 24000, "borrowed"), (40000, 56000, "predicted")]}
    assert efg.walk_extent(windows, "c", 5000, 5300) == (0, 24000)
    assert efg.walk_extent(windows, "c", 45000, 45300) == (40000, 56000)
    assert efg.walk_extent(windows, "other", 5000, 5300) is None
    assert efg.walk_extent({}, "c", 5000, 5300) is None
    assert efg.walk_extent(windows, "c", 23000, 40500) == (0, 24000)    # centre in neither: overlap
    assert efg.walk_extent(windows, "c", 23900, 41000) == (40000, 56000)


def test_the_predictions_name_their_windows_sorted(tmp_path):
    id_map = {"c__40000__56000": ("c", 40000, 56000), "b__0__900": ("b", 0, 900),
              "c__0__12000": ("c", 0, 12000)}
    empty = tmp_path / "empty.gff"
    empty.write_text("")
    out = tmp_path / "home_predicted.gff"
    por.adjust_gff(str(empty), str(out), id_map)
    assert out.read_text().splitlines() == [
        "##gff-version 3",
        "##synvoy-annotated-window b 1 900 predicted",
        "##synvoy-annotated-window c 1 12000 predicted",
        "##synvoy-annotated-window c 40001 56000 predicted"]
    raw = tmp_path / "pred.gff"
    raw.write_text("c__0__12000\tProdigal\tCDS\t1001\t1600\t50\t+\t0\tID=1_1\n")
    por.adjust_gff(str(raw), str(out), id_map)
    assert read_annotated_windows(str(out)) == {
        "b": [(0, 900, "predicted")], "c": [(0, 12000, "predicted"), (40000, 56000, "predicted")]}
    assert len([ln for ln in out.read_text().splitlines() if not ln.startswith("#")]) == 3
