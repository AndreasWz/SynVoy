"""Fix-phase d3 (item 9; magic-numbers report #1): the easy-mode home chromosome filter.

fetch_home_genome.filter_chromosomes_only keeps only chromosome-scale records once the home FASTA holds an NC_ or CM record,
or a record of 1 Mb or more that is not NW_ or NZ_. So a scaffold-level RefSeq assembly whose only NC_ record is its
mitochondrion kept only the mitochondrion, and a GOI on an unplaced scaffold of a chromosome-level assembly was deleted with
its scaffold: the run then anchored on whatever else the query hit, with one stdout line. Now an organelle is never evidence of
a chromosome-level assembly, the removed records are kept beside the genome (home_genome.removed.fna, .removed.tsv), and the
step that places the GOI (annotate_goi_exons.py) stops when no kept sequence holds as good a match as a removed one.
Constructed FASTAs only (FETCH_HOME_GENOME downloads from NCBI by species and cannot run on a local file).
"""
import json
import random
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import annotate_goi_exons as ag  # noqa: E402
import fetch_home_genome as fh  # noqa: E402

CODON = {'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT', 'Q': 'CAA', 'E': 'GAA',
         'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT',
         'P': 'CCT', 'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'}
INTRON = 'GT' + 'ACGTTGCA' * 30 + 'AG'
SPACER = 'TAACTAGTTGA' * 40
AA = 'ARNDCQEGHILKFPSTWYV'
THRESH = 1000           # min_chrom_len for the constructed records


def _protein(length, seed):
    rng = random.Random(seed)
    return 'M' + ''.join(rng.choice(AA) for _ in range(length - 1))


def _mutate(protein, identity, seed):
    rng = random.Random(seed)
    pos = rng.sample(range(1, len(protein)), int(round((1 - identity) * (len(protein) - 1))))
    out = list(protein)
    for i in pos:
        out[i] = rng.choice([a for a in AA if a != out[i]])
    return ''.join(out)


Q = _protein(200, 1)
PARALOG = _mutate(Q, 0.72, 2)


def _gene_dna(prot):
    half = len(prot) // 2
    e1 = ''.join(CODON[a] for a in prot[:half])
    e2 = ''.join(CODON[a] for a in prot[half:]) + 'TAA'
    return e1, e2


def _assembly(tmp_path, records):
    """records: [(seqid, header_rest, [(gene_id, name, protein)])]: each gene two exons. Writes the whole assembly
    (home_genome.fna, before the filter), its GFF, and LOCATE-style hits of the query on every gene whose protein is
    not the query itself (what the query hits once its own gene is gone). Returns paths and spans."""
    fasta, gff, hits = [], ["##gff-version 3\n"], []
    for seqid, rest, genes in records:
        seq = SPACER
        for gid, name, prot in genes:
            e1, e2 = _gene_dna(prot)
            s1 = len(seq)
            seq += e1 + INTRON
            s2 = len(seq)
            seq += e2 + SPACER
            spans = [(s1, s1 + len(e1)), (s2, s2 + len(e2))]
            gff.append(f"{seqid}\tsrc\tgene\t{s1 + 1}\t{spans[1][1]}\t.\t+\t.\tID={gid};Name={name}\n")
            gff.append(f"{seqid}\tsrc\tmRNA\t{s1 + 1}\t{spans[1][1]}\t.\t+\t.\tID=rna-{gid};Parent={gid}\n")
            for a, b in spans:
                gff.append(f"{seqid}\tsrc\tCDS\t{a + 1}\t{b}\t.\t+\t0\tID=cds-{gid};Parent=rna-{gid}\n")
            if prot != Q:
                for a, b in spans:
                    hits.append(f"q\t{seqid}\t72.0\t60\t10\t0\t1\t60\t{a + 1}\t{b}\t1e-20\t90\n")
        fasta.append(f">{seqid} {rest}\n{seq}\n")
    (tmp_path / "home_genome.fna").write_text("".join(fasta))
    (tmp_path / "home.gff").write_text("".join(gff))
    (tmp_path / "hits.tsv").write_text("".join(hits))
    (tmp_path / "q.faa").write_text(f">q\n{Q}\n")
    return tmp_path / "home_genome.fna"


def _ids(path):
    return [line[1:].split()[0] for line in open(path) if line.startswith(">")]


def _annotate(tmp_path, monkeypatch, removed=True, gff=True, query_id="GOISYM"):
    monkeypatch.setattr(ag, "fetch_uniprot_names", lambda _q: [])
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    argv = ["annotate_goi_exons.py", "--query", str(tmp_path / "q.faa"), "--genome", str(tmp_path / "home_genome.fna"),
            "--gff", str(tmp_path / "home.gff") if gff else "NO_GFF", "--blast_hits", str(tmp_path / "hits.tsv"),
            "--query_id", query_id, "--max_evalue", "0.01",
            "--output_exons", str(out / "e.faa"), "--output_bed", str(out / "e.bed"),
            "--output_info", str(out / "info.json")]
    if removed:
        argv += ["--removed_sequences", str(tmp_path / "home_genome.removed.fna")]
    monkeypatch.setattr(sys, "argv", argv)
    ag.main()
    return json.loads((out / "info.json").read_text())


@pytest.fixture(autouse=True)
def _need_parasail():
    pytest.importorskip("parasail")


# --------------------------------------------------------------------------- the filter

def test_a_scaffold_assembly_with_a_mitochondrion_keeps_everything(tmp_path):
    """The mitochondrion's NC_ record is no evidence: before, it was the one record kept."""
    p = _assembly(tmp_path, [("NC_045177.1", "Apis sp. mitochondrion, complete genome", []),
                             ("NW_1.1", "scaffold_1", [("gene-G", "GOISYM", Q)]),
                             ("NW_2.1", "scaffold_2", [])])
    s = fh.filter_chromosomes_only(p, min_chrom_len=THRESH)
    assert _ids(p) == ["NC_045177.1", "NW_1.1", "NW_2.1"]
    assert s["removed"] == 0 and s["organelles"] == ["NC_045177.1"]
    removed_fna, removed_tsv = fh.removed_paths(p)
    assert removed_fna.read_text() == "" and removed_tsv.read_text() == "seqid\tlength_bp\n"


@pytest.mark.parametrize("header", ["chloroplast, complete genome", "plastid genome", "apicoplast, complete sequence",
                                    "mitochondrial DNA"])
def test_other_organelles_are_no_evidence_either(tmp_path, header):
    p = _assembly(tmp_path, [("NC_9.1", header, []), ("NW_1.1", "scaffold_1", [])])
    assert fh.filter_chromosomes_only(p, min_chrom_len=THRESH)["removed"] == 0 and len(_ids(p)) == 2


def test_an_organelle_is_kept_in_a_chromosome_assembly(tmp_path):
    p = _assembly(tmp_path, [("NC_1.1", "chromosome 1", []), ("NC_2.1", "mitochondrion, complete genome", []),
                             ("NW_1.1", "unplaced scaffold", [])])
    s = fh.filter_chromosomes_only(p, min_chrom_len=THRESH)
    assert _ids(p) == ["NC_1.1", "NC_2.1"] and s["removed"] == 1


def test_fly_like_NT_arms_stay_and_nothing_is_listed(tmp_path):
    p = _assembly(tmp_path, [("NC_004354.4", "chromosome X", []),
                             ("NT_033777.3", "chromosome 3R", [("gene-G", "GOISYM", Q)] * 3)])
    assert len(open(p).read()) > THRESH
    s = fh.filter_chromosomes_only(p, min_chrom_len=THRESH)
    assert _ids(p) == ["NC_004354.4", "NT_033777.3"] and s["removed"] == 0


def test_human_like_NT_scaffolds_are_removed_and_listed(tmp_path):
    p = _assembly(tmp_path, [("NC_000012.12", "chromosome 12", [("gene-G", "GOISYM", Q)]),
                             ("NT_187361.1", "unlocalized scaffold", []), ("NW_009646201.1", "alt locus", [])])
    s = fh.filter_chromosomes_only(p, min_chrom_len=100_000)
    assert _ids(p) == ["NC_000012.12"]
    removed_fna, removed_tsv = fh.removed_paths(p)
    assert _ids(removed_fna) == ["NT_187361.1", "NW_009646201.1"]
    rows = [r.split("\t") for r in removed_tsv.read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["NT_187361.1", "NW_009646201.1"]
    assert s["removed"] == 2 and s["removed_bp"] == sum(int(r[1]) for r in rows) == len(SPACER) * 2


# --------------------------------------------------------------------------- the GOI step

def _chromosome_level_with_goi_on_a_scaffold(tmp_path):
    p = _assembly(tmp_path, [("NC_1.1", "chromosome 1", [("gene-P", "PARSYM", PARALOG)]),
                             ("NW_7.1", "unplaced genomic scaffold", [("gene-G", "GOISYM", Q)])])
    fh.filter_chromosomes_only(p, min_chrom_len=2000)
    assert _ids(p) == ["NC_1.1"] and _ids(fh.removed_paths(p)[0]) == ["NW_7.1"]
    return p


def test_without_the_removed_records_the_paralog_is_taken_in_silence(tmp_path, monkeypatch):
    """The parent's behaviour (the step had no removed records): the GOI's own gene is skipped, its paralog is the home
    gene, exit 0."""
    _chromosome_level_with_goi_on_a_scaffold(tmp_path)
    info = _annotate(tmp_path, monkeypatch, removed=False)
    assert info["gene_id"] == "gene-P" and info["found_by"] == "position"
    assert "home_filter" not in info


def test_the_goi_found_by_name_on_a_removed_scaffold_stops(tmp_path, monkeypatch, capsys):
    _chromosome_level_with_goi_on_a_scaffold(tmp_path)
    with pytest.raises(SystemExit) as stop:
        _annotate(tmp_path, monkeypatch)
    assert stop.value.code == 1
    err = capsys.readouterr().err
    assert "gene gene-G, found by name" in err and "NW_7.1" in err and "pro mode" in err


def test_a_copy_on_a_removed_alt_locus_does_not_stop(tmp_path, monkeypatch):
    """Human alt loci hold near-identical copies under the same name: the kept chromosome's gene wins first, as before."""
    p = _assembly(tmp_path, [("NC_1.1", "chromosome 1", [("gene-G", "GOISYM", Q)]),
                             ("NW_7.1", "alt locus", [("gene-G-2", "GOISYM", Q)])])
    fh.filter_chromosomes_only(p, min_chrom_len=2000)
    info = _annotate(tmp_path, monkeypatch)
    assert info["gene_id"] == "gene-G" and info["found_by"] == "name"
    assert info["home_filter"]["removed_sequences"] == 1 and info["home_filter"]["removed_bp"] > 0


def test_nothing_removed_records_zero(tmp_path, monkeypatch):
    p = _assembly(tmp_path, [("NC_1.1", "chromosome 1", [("gene-G", "GOISYM", Q)])])
    fh.filter_chromosomes_only(p, min_chrom_len=THRESH)
    info = _annotate(tmp_path, monkeypatch)
    assert info["home_filter"] == {"removed_sequences": 0, "removed_bp": 0}


@pytest.mark.skipif(not shutil.which("miniprot"), reason="miniprot not installed")
def test_without_a_name_a_better_model_on_a_removed_scaffold_stops(tmp_path, monkeypatch, capsys):
    """No GFF: the hits build a model of the paralog on the kept chromosome (72 %); miniprot on the removed records finds
    the query itself, so the step stops and names the scaffold."""
    _chromosome_level_with_goi_on_a_scaffold(tmp_path)
    with pytest.raises(SystemExit) as stop:
        _annotate(tmp_path, monkeypatch, gff=False)
    assert stop.value.code == 1
    err = capsys.readouterr().err
    assert "a miniprot model at NW_7.1:" in err and "score fraction" in err


@pytest.mark.skipif(not shutil.which("miniprot"), reason="miniprot not installed")
def test_without_a_name_a_good_placement_is_not_checked_against_a_worse_copy(tmp_path, monkeypatch):
    """The query's own gene is on the kept chromosome; a removed scaffold holds only the paralog: no stop."""
    p = _assembly(tmp_path, [("NC_1.1", "chromosome 1", [("gene-G", "GOISYM", Q), ("gene-X", "X", _protein(150, 9))]),
                             ("NW_7.1", "unplaced genomic scaffold", [("gene-P", "PARSYM", PARALOG)])])
    fh.filter_chromosomes_only(p, min_chrom_len=2000)
    # hits on the query's own gene (the builder writes hits only for genes other than the query)
    a = [line for line in open(tmp_path / "home.gff") if "\tCDS\t" in line and line.endswith("Parent=rna-gene-G\n")]
    assert len(a) == 2
    (tmp_path / "hits.tsv").write_text("".join(
        f"q\tNC_1.1\t100.0\t100\t0\t0\t1\t100\t{r.split()[3]}\t{r.split()[4]}\t1e-50\t200\n" for r in a))
    info = _annotate(tmp_path, monkeypatch, gff=False)
    assert info["method"] == "hit_annotation" and info["home_filter"]["removed_sequences"] == 1


# --------------------------------------------------------------------------- the report

def _report(tmp_path, info):
    import subprocess
    stage = tmp_path / "staged"
    (stage / "regions").mkdir(parents=True)
    (stage / "regions" / "locus_1__Gx.fna.gff").write_text("##gff-version 3\n")
    gi = tmp_path / "goi_info.json"
    gi.write_text(json.dumps(info))
    out = tmp_path / "report.json"
    proc = subprocess.run([sys.executable, str(ROOT / "bin" / "generate_report.py"), "--results_dir", str(stage),
                           "--goi_info", str(gi), "--output", str(out), "--allow-empty"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return json.loads(out.read_text())


def test_the_report_records_what_the_filter_removed(tmp_path):
    rep = _report(tmp_path, {"method": "gff_annotation", "home_filter": {"removed_sequences": 3, "removed_bp": 1234}})
    assert rep["home_genome_filter"] == {"removed_sequences": 3, "removed_bp": 1234}


def test_a_pro_mode_report_has_no_filter_record(tmp_path):
    assert "home_genome_filter" not in _report(tmp_path, {"method": "gff_annotation"})


@pytest.mark.skipif(not shutil.which("miniprot"), reason="miniprot not installed")
def test_an_unannotated_isoform_on_a_kept_chromosome_does_not_stop(tmp_path, monkeypatch):
    """First review of the session: the kept gene's annotated transcript lacks the query's middle exon and an
    identical copy of the gene lies on a removed scaffold. The placed transcript agrees at a score fraction of 0.59, a
    miniprot model of the copy at about 1; the kept genome around the query's hits holds the same model, so the step goes
    on (it stopped when the copy's model was compared with the annotated transcript)."""
    third = len(Q) // 3
    exons = [Q[:third], Q[third:2 * third], Q[2 * third:]]

    def gene():
        seq, spans = SPACER, []
        for k, part in enumerate(exons):
            dna = "".join(CODON[a] for a in part) + ("TAA" if k == 2 else "")
            spans.append((len(seq), len(seq) + len(dna)))
            seq += dna + (INTRON if k < 2 else SPACER)
        return seq, spans
    s1, sp = gene()
    s2, _ = gene()
    (tmp_path / "home_genome.fna").write_text(f">NC_1.1 chromosome 1\n{s1}\n>NW_7.1 alt locus\n{s2}\n")
    rows = ["##gff-version 3\n", f"NC_1.1\tsrc\tgene\t{sp[0][0] + 1}\t{sp[2][1]}\t.\t+\t.\tID=gene-G;Name=FOO1\n",
            f"NC_1.1\tsrc\tmRNA\t{sp[0][0] + 1}\t{sp[2][1]}\t.\t+\t.\tID=rna-G;Parent=gene-G\n"]
    rows += [f"NC_1.1\tsrc\tCDS\t{a + 1}\t{b}\t.\t+\t0\tID=cds-G;Parent=rna-G\n" for a, b in (sp[0], sp[2])]
    (tmp_path / "home.gff").write_text("".join(rows))
    (tmp_path / "hits.tsv").write_text("".join(f"q\tNC_1.1\t100.0\t100\t0\t0\t1\t100\t{a + 1}\t{b}\t1e-50\t200\n"
                                               for a, b in sp))
    (tmp_path / "q.faa").write_text(f">q\n{Q}\n")
    fh.filter_chromosomes_only(tmp_path / "home_genome.fna", min_chrom_len=10 ** 7)
    assert _ids(fh.removed_paths(tmp_path / "home_genome.fna")[0]) == ["NW_7.1"]
    info = _annotate(tmp_path, monkeypatch, query_id="ZZZ9")
    assert info["gene_id"] == "gene-G" and info["found_by"] == "position"
