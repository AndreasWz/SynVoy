"""The benchmark's SynVoy calls (scripts/benchmark/synvoy_calls.py) come from the report's
records: PRESENT = a counted HIGH or MEDIUM record, AMBIGUOUS = a record that needs curation,
else ABSENT. The parsers it replaces read LOW rows (the class ambiguous_goi_family_member) as
AMBIGUOUS, counted every other class as resolved, and took a region label as a call."""
import importlib.util
import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SCRIPT = os.path.join(ROOT, "scripts", "benchmark", "synvoy_calls.py")
spec = importlib.util.spec_from_file_location("synvoy_calls", SCRIPT)
sc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc)


def _rec(genome, conf, goi_class, start=1001, end=1600, chrom="c1", mid="GOI_q|x_m1", **kw):
    return {"genome": genome, "chrom": chrom, "start": start, "end": end, "confidence": conf,
            "goi_class": goi_class, "mrna_id": mid, **kw}


def _report(records, per_genome=()):
    return {"goi_dedup": {"records": records},
            "annotations": {"per_genome": list(per_genome)},
            "regions": {"per_genome": []}}


def _status(rows):
    return {r["species"]: r["called_status"] for r in rows}


def test_the_probe_cases():
    """GA: a LOW strong-synteny rescue model and a DISP MEDIUM region; GB: LOW rows only; GC: one
    AMBIGUOUS record and a MEDIUM region; GD: a LOW tandem copy and a HIGH region. The old parser
    gave PRESENT, AMBIGUOUS, PRESENT, PRESENT."""
    per_genome = [{"genome": g, "goi_confidence_counts": {"LOW": n}} for g, n in
                  (("GA", 1), ("GB", 3), ("GC", 0), ("GD", 1))]
    report = _report([_rec("GC", "AMBIGUOUS", "syntenic_candidate_unconfirmed", basis="gene_order_only")],
                     per_genome)
    regions = {"GA": "Reg1_DISP_CMEDIUM_S0.45", "GC": "Reg1_G3_CMEDIUM_S0.40", "GD": "Reg1_G9_CHIGH_S0.90"}
    rows = sc.calls(report, regions=regions)
    assert _status(rows) == {"GA": "ABSENT", "GB": "ABSENT", "GC": "AMBIGUOUS", "GD": "ABSENT"}
    gb = next(r for r in rows if r["species"] == "GB")
    assert gb["extra"] == "counted=0;curation=0;low=3"
    assert next(r for r in rows if r["species"] == "GD")["extra"].endswith(";region=Reg1_G9_CHIGH_S0.90")


def test_a_counted_record_is_present_at_its_coordinates():
    report = _report([_rec("Sp.fna", "MEDIUM", "probable_goi", 50_001, 50_900, chrom="NC_1", mid="m9")])
    strands = {("m9", "NC_1", 50_001, 50_900): {"-"}}
    row = sc.calls(report, strands)[0]
    assert (row["species"], row["called_status"], row["locus_chrom"], row["locus_start"], row["locus_end"],
            row["strand"], row["confidence"]) == ("Sp", "PRESENT", "NC_1", "50001", "50900", "-", "MEDIUM")


def test_uncounted_records_are_not_present():
    report = _report([_rec("G1", "HIGH", "paralog_not_goi"),
                      _rec("G2", "HIGH", "identity_coverage_decoupled"),
                      _rec("G3", "MEDIUM", "genome_wide_strong_hit", genome_wide_strong_hit=True),
                      _rec("G4", "AMBIGUOUS", "paralog_not_goi")])
    assert _status(sc.calls(report)) == {"G1": "ABSENT", "G2": "ABSENT", "G3": "ABSENT", "G4": "ABSENT"}


def test_the_first_counted_record_gives_the_locus():
    report = _report([_rec("G", "AMBIGUOUS", "syntenic_candidate_unconfirmed", 1, 100),
                      _rec("G", "HIGH", "confident_goi", 5_001, 5_900),
                      _rec("G", "MEDIUM", "probable_goi", 9_001, 9_900)])
    row = sc.calls(report)[0]
    assert (row["called_status"], row["locus_start"], row["confidence"]) == ("PRESENT", "5001", "HIGH")
    assert row["extra"] == "counted=2;curation=1;low=0"
    assert row["strand"] == "-"                     # no GFF row: the placeholder


def test_the_cli_writes_the_contract(tmp_path):
    results = tmp_path / "run"
    (results / "plot_inputs_locus_1").mkdir(parents=True)
    (results / "regions" / "locus_1").mkdir(parents=True)
    (results / "synvoy_report.json").write_text(json.dumps(
        _report([_rec("G.fna", "HIGH", "confident_goi", 101, 400, mid="GOI_q|G_m1")])))
    (results / "plot_inputs_locus_1" / "G.fna.gff").write_text(
        "c1\tSynVoy\tmRNA\t101\t400\t.\t+\t.\tID=GOI_q|G_m1;SynVoyRole=goi;Confidence=HIGH\n")
    (results / "regions" / "locus_1" / "H.fna.regions.bed").write_text(
        "c2\t0\t1000\tReg1_G4_CMEDIUM_S0.40\t0.40\t+\n")
    out = tmp_path / "bench" / "calls.tsv"
    subprocess.run([sys.executable, SCRIPT, str(results), str(out)], check=True, capture_output=True)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == sc.COLUMNS
    assert lines[1].split("\t")[:8] == ["G", "-", "PRESENT", "c1", "101", "400", "+", "HIGH"]
    assert lines[2].split("\t")[:3] == ["H", "-", "ABSENT"]


def test_the_runners_call_it():
    for name in ("run_synvoy.sh", "run_tp53.sh"):
        text = open(os.path.join(ROOT, "scripts", "benchmark", name)).read()
        assert "scripts/benchmark/synvoy_calls.py" in text
        assert "resolved_goi_annotations" not in text and "ambiguous_goi_annotations" not in text


# --------------------------------------------------------------------------- genome names (review round 2)

def _gr():
    spec = importlib.util.spec_from_file_location("generate_report", os.path.join(ROOT, "bin", "generate_report.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


NAMES = ["Apis_mellifera.fasta", "Apis_mellifera.fna", "Apis_mellifera.fa", "Apis_mellifera",
         "GCA_000001405.29_GRCh38_genomic.fna", "GCF_000002315.6", "locus_2__Sp_x.fna.gff",
         "Sp_x.fna.hull_rescue.gff", "Sp_x.fasta.regions.bed", "Sp_x.fna.gz", "Sp_y_new_genes"]


def test_the_genome_id_is_the_reports():
    gr = _gr()
    for name in NAMES:
        assert sc.report_genome_id(name) == gr.canonical_genome_id(name), name


def test_a_region_bed_names_the_reports_genome(tmp_path):
    """A .fasta target and an NCBI-style long name: one row per genome, the BED's region on it
    (the BED is named after the genome file, the report after its canonical id)."""
    gr = _gr()
    for genome_file in ("Apis_mellifera.fasta", "GCA_000001405.29_GRCh38_genomic.fna", "Apis_mellifera.fna",
                        "Sp_z.fna.gz"):
        results = tmp_path / genome_file / "run"
        (results / "regions" / "locus_1").mkdir(parents=True)
        g = gr.canonical_genome_id(genome_file)
        (results / "synvoy_report.json").write_text(json.dumps(
            _report([_rec(g, "HIGH", "confident_goi", 101, 400, mid="m1")],
                    [{"genome": g, "goi_confidence_counts": {"HIGH": 1}}])))
        (results / "regions" / "locus_1" / f"{genome_file}.regions.bed").write_text(
            "c1\t0\t1000\tReg1_G4_CHIGH_S0.80\t0.80\t+\n")
        out = tmp_path / genome_file / "calls.tsv"
        subprocess.run([sys.executable, SCRIPT, str(results), str(out)], check=True, capture_output=True)
        rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
        assert len(rows) == 1, (genome_file, rows)
        assert rows[0][2] == "PRESENT" and rows[0][-1].endswith(";region=Reg1_G4_CHIGH_S0.80"), rows
