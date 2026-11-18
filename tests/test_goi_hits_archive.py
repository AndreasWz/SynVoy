"""The GOI rows of each genome's genome-wide hit table are published.

hits/<genome>.m8, the genome-wide MMseqs2 table, sits inside the hits/ output folder and so never
reaches the run folder: the hits behind a refusal or a drop of the dispersed-GOI rescue could not
be read after the run. ITERATIVE_SEARCH now writes the rows whose query is a
GOI query to goi_hits/<prefix><genome>.goi.m8.gz (sorted, gzipped with mtime 0) and publishes the
folder. Synthetic ids only.
"""
import gzip
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402

ROWS = [
    "gene-FLANK1\tc1\t80.0\t300\t0\t0\t1\t300\t1000\t1900\t1e-50\t500\n",
    "GOI_q\tc2\t60.0\t200\t0\t0\t5\t204\t9000\t9600\t1e-30\t300\n",
    "GOI_q|exon_2\tc1\t55.0\t80\t0\t0\t1\t80\t4000\t4240\t1e-10\t90\n",
    "gene-FLANK2\tc3\t70.0\t250\t0\t0\t1\t250\t100\t850\t1e-40\t400\n",
    "GOI_q|H_fna_b0_l1_exon_ann\tc2\t58.0\t190\t0\t0\t1\t190\t9010\t9580\t1e-28\t280",   # no newline
]


def _m8(path, rows):
    with open(path, "w") as fh:
        fh.writelines(rows)


def test_the_archive_holds_the_goi_rows_sorted_and_reproducible(tmp_path):
    m8 = tmp_path / "hits" / "G.fna.m8"
    m8.parent.mkdir()
    _m8(m8, ROWS)
    out = tmp_path / "goi_hits" / "locus_1_G.fna.goi.m8.gz"
    assert isr.archive_goi_rows(str(m8), str(out)) == 3
    lines = gzip.decompress(out.read_bytes()).decode().splitlines()
    assert [ln.split("\t")[0] for ln in lines] == ["GOI_q", "GOI_q|H_fna_b0_l1_exon_ann", "GOI_q|exon_2"]
    assert lines == sorted(lines) and all(len(ln.split("\t")) == 12 for ln in lines)
    first = out.read_bytes()
    _m8(m8, list(reversed(ROWS[:-1])) + [ROWS[-1]])        # another row order (MMseqs2 threads)
    isr.archive_goi_rows(str(m8), str(out))
    assert out.read_bytes() == first                       # same bytes: no time stamp, no name
    assert not os.path.exists(str(out) + ".tmp")


def test_no_goi_row_writes_no_file(tmp_path):
    m8 = tmp_path / "G.fna.m8"
    _m8(m8, [ROWS[0], ROWS[3]])
    out = tmp_path / "goi_hits" / "G.fna.goi.m8.gz"
    assert isr.archive_goi_rows(str(m8), str(out)) == 0 and not out.exists()
    assert isr.archive_goi_rows(str(tmp_path / "missing.m8"), str(out)) == 0 and not out.exists()


def test_the_search_archives_before_any_early_return():
    src = open(os.path.join(ROOT, "bin", "iterative_search_runner.py")).read()
    genome = src.split("def process_single_genome(", 1)[1].split("\ndef ", 1)[0]
    call = genome.index("archive_goi_rows(hits_file, os.path.join(args.output_dir, 'goi_hits', "
                        "f\"{prefix}{genome_name}.goi.m8.gz\"))")
    assert genome.index("split_decoy_rows(hits_file") < call       # the real rows only
    assert call < genome.index("No hits found in MMseqs output.")


def test_the_module_publishes_the_folder():
    module = open(os.path.join(ROOT, "modules", "iterative_search.nf")).read()
    assert "pattern: 'iterative_results/goi_hits/*'" in module
    assert '"${params.outdir}/goi_hits/${locus_id}"' in module
    assert 'path("iterative_results/goi_hits/*"), emit: goi_hits, optional: true' in module
    # Nextflow publishes a file inside a declared output folder only with that folder: no
    # output folder may hold goi_hits/.
    outputs = re.findall(r'path\("(iterative_results/[^"*]+)"\)', module)
    assert outputs and not any("iterative_results/goi_hits/".startswith(o.rstrip("/") + "/")
                               for o in outputs)
