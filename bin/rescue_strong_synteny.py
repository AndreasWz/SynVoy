#!/usr/bin/env python3
"""Relaxed-miniprot rescue pass for ``goi_missing_but_strong_synteny`` blocks
(docs/TODO.md §1e follow-up).

cluster_grs.py emits a ``region_class=goi_missing_but_strong_synteny`` row in
``<genome>.scores.tsv`` whenever a syntenic block has many HIGH-confidence
flanking-gene hits but no GOI model. The user sees the block coordinates but
nothing inside them — making it hard to judge whether the orthologue really is
too divergent or whether the strict miniprot pass just missed it.

This script runs a SECOND, much more permissive miniprot pass strictly inside
those blocks against the GOI query and emits any models as standard SynVoy GFF
rows tagged ``EvidenceType=relaxed_miniprot_rescue`` and ``Confidence=LOW``
(so they're visible in the report without polluting the high-confidence picks).
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import (genetic_code_flag, genetic_code_for,  # noqa: E402
                            parse_genetic_code_spec, set_genetic_code,
                            exit_step_failed, guard_failed, StepFailed)
from sequence_utils import parse_gff_attributes  # noqa: E402,F401  (kept for symmetry)


def _read_fasta_subseq(samtools_bin: str, genome_path: str,
                       chrom: str, start: int, end: int) -> str:
    """Extract genome[chrom:start..end] (1-based inclusive) via ``samtools faidx``."""
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        [samtools_bin, "faidx", genome_path, region],
        check=True, capture_output=True, text=True,
    )
    return "".join(line for line in result.stdout.splitlines() if not line.startswith(">"))


def _read_query(query_path: str) -> tuple[str, str]:
    """Return (header_id, sequence) of the first FASTA record in ``query_path``."""
    name = ""
    seq_parts: list[str] = []
    with open(query_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_parts:
                    break
                name = line[1:].split()[0]
            else:
                seq_parts.append(line)
    return name, "".join(seq_parts)


def _run_miniprot_relaxed(window_fa: str, query_fa: str, min_outc: float,
                          timeout: int, param: str = "strong_synteny_rescue_miniprot_timeout",
                          place: str = "") -> str:
    """Run miniprot --gff with very permissive thresholds against a single window.

    Mirrors the relaxed search profile used elsewhere in iterative_search_runner
    (-n 2 / -p 0.2 / -N 50 / --outs 0.2 / --outc <min_outc>). Returns the raw GFF3
    text; the caller parses it. ``timeout`` is the run's guard (pipeline parameter ``param``):
    a run that does not finish raises StepFailed naming ``place`` (genome, block or hull and its
    coordinates), and the script exits STEP_FAILED_EXIT, which the module turns into the end of
    the run (fix-phase item 18). Until 2026-10-28 such a
    window was skipped and its GFF comment read like a real "no model".
    """
    cmd = [
        "miniprot",
        "--gff",
        "-t", "1",
        # Sensitivity: lower seed thresholds, allow weak secondary alignments.
        "-n", "2", "-p", "0.2", "-N", "50",
        "--outs=0.2",
        f"--outc={min_outc}",
        # Short-peptide friendliness — many ant venom precursors are <100 aa.
        "-L", "10",
    ] + genetic_code_flag("miniprot") + [   # the target's genetic code ([] for table 1)
        window_fa, query_fa,
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise guard_failed(f"[rescue] miniprot on {place or 'the window ' + window_fa}", timeout,
                           param) from None
    except FileNotFoundError as exc:   # an installation fault, not a window: the run stops too
        raise StepFailed(f"[rescue] miniprot is not installed here ({exc})") from exc
    if proc.returncode != 0:
        print(f"[rescue] miniprot non-zero exit ({proc.returncode}) on "
              f"{window_fa}: {proc.stderr[:200]}", file=sys.stderr)
        return ""
    return proc.stdout or ""


def _parse_miniprot_gff(gff_text: str):
    """Yield (mRNA_row, [CDS_rows]) tuples from miniprot --gff output."""
    cur_mrna = None
    cur_cds: list[list[str]] = []
    for raw in gff_text.splitlines():
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "mRNA":
            if cur_mrna is not None:
                yield cur_mrna, cur_cds
            cur_mrna = parts
            cur_cds = []
        elif parts[2] == "CDS" and cur_mrna is not None:
            cur_cds.append(parts)
    if cur_mrna is not None:
        yield cur_mrna, cur_cds


def _refine_rescue_termini(window_seq: str, mrna, cds_rows) -> None:
    """Give a rescue model the gene's own start and stop codon, in place.

    The rescue passes call miniprot directly, so they bypass the terminal refinement
    ``annotate_goi_exons`` applies on the main path (CLAUDE.md §21). Measured on the
    2026-09-15 boundary audit: for 8 of 31 recovered apyrase genes the coordinates the
    report shows came from an unrefined hull-rescue model, 9-24 bp short of the curated
    gene at the 5' end, while the refined in-block model of the same gene matched exactly.

    ``mrna`` and ``cds_rows`` are GFF fields in WINDOW coordinates (1-based inclusive,
    ascending); both are updated. Refusal rules and bounds are the shared helper's.
    """
    if not cds_rows or not window_seq:
        return
    try:
        from annotate_goi_exons import refine_model_termini
    except ImportError:
        return
    strand = cds_rows[0][6]
    rows = sorted(cds_rows, key=lambda c: int(c[3]))
    coding = list(reversed(rows)) if strand == "-" else rows
    cds = []
    for r in coding:
        try:
            phase = int(r[7])
        except (TypeError, ValueError):
            phase = 0
        cds.append({"gstart": int(r[3]) - 1, "gend": int(r[4]), "phase": phase,
                    "qstart": 0, "qend": 0})
    refine_model_termini(cds, window_seq, strand)
    for r, c in zip(coding, cds):
        r[3] = str(int(c["gstart"]) + 1)
        r[4] = str(int(c["gend"]))
    mrna[3] = str(min(int(r[3]) for r in cds_rows))
    mrna[4] = str(max(int(r[4]) for r in cds_rows))


def _get_attr(attr_str: str, key: str, default: str = "") -> str:
    for kv in attr_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            if k.strip() == key:
                return v.strip()
    return default


def _identity_from_mrna_attrs(attrs: str) -> float:
    """miniprot writes Identity as a fraction (0.962); every other SynVoy model carries
    percent. Written as is, the rescue models read Identity=1.0 for 96 %."""
    raw = _get_attr(attrs, "Identity", "")
    try:
        value = float(raw)
    except ValueError:
        return 0.0
    return value * 100.0 if value <= 1.0 else value


def _build_rescue_gff_rows(mrna, cds_rows, chrom: str, window_start: int,
                           query_id: str, parent_id: str) -> list[str]:
    """Map miniprot's window-relative coords back to genome coords and produce
    a (gene, mRNA, CDS×N) block of SynVoy-style GFF rows."""
    mp_start = int(mrna[3])
    mp_end = int(mrna[4])
    strand = mrna[6]
    identity = _identity_from_mrna_attrs(mrna[8])
    gstart = window_start + mp_start - 1
    gend = window_start + mp_end - 1
    gene_attrs = (
        f"ID={parent_id};Name={query_id};SynVoyRole=goi;Confidence=LOW;"
        f"EvidenceType=relaxed_miniprot_rescue;GOIClass=relaxed_rescue;"
        f"ModelStatus=rescue;Identity={identity:.1f};TargetGene={query_id}"
    )
    mrna_attrs = (
        f"ID={parent_id}.mRNA;Parent={parent_id};Name={query_id};SynVoyRole=goi;"
        f"Confidence=LOW;EvidenceType=relaxed_miniprot_rescue;GOIClass=relaxed_rescue;"
        f"ModelStatus=rescue;Identity={identity:.1f};TargetGene={query_id}"
    )
    rows = [
        "\t".join([chrom, "SynVoy_rescue", "gene", str(gstart), str(gend),
                   ".", strand, ".", gene_attrs]),
        "\t".join([chrom, "SynVoy_rescue", "mRNA", str(gstart), str(gend),
                   ".", strand, ".", mrna_attrs]),
    ]
    for i, cds in enumerate(cds_rows, start=1):
        cs = window_start + int(cds[3]) - 1
        ce = window_start + int(cds[4]) - 1
        cds_attrs = (
            f"ID={parent_id}.CDS.{i};Parent={parent_id}.mRNA;SynVoyRole=goi;"
            f"Confidence=LOW;EvidenceType=relaxed_miniprot_rescue"
        )
        rows.append("\t".join([chrom, "SynVoy_rescue", "CDS",
                               str(cs), str(ce), ".", strand, cds[7], cds_attrs]))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="cluster_grs <genome>.scores.tsv")
    ap.add_argument("--target_genome", required=True, help="Target genome FASTA (.fa[.gz])")
    ap.add_argument("--query", required=True, help="GOI query protein FASTA")
    ap.add_argument("--genome_name", required=True, help="Logical genome name for ID stems")
    ap.add_argument("--output", required=True, help="Output GFF path")
    ap.add_argument("--samtools_bin", default="samtools")
    ap.add_argument("--miniprot_outc", type=float, default=0.05,
                    help="--outc passed to miniprot (default 0.05 = 5%% query coverage)")
    ap.add_argument("--window_pad", type=int, default=2000,
                    help="bp padding added each side of the block coords")
    ap.add_argument("--miniprot_timeout", type=int, default=3600,
                    help="Guard of the miniprot run of one block window, in seconds (pipeline: "
                         "strong_synteny_rescue_miniprot_timeout): a run that does not finish ends "
                         "the run (exit STEP_FAILED_EXIT)")
    ap.add_argument("--target_genetic_codes", default="",
                    help="Genetic codes of target genomes, 'genome=table,...'; the table of "
                         "--target_genome is used (lib/GeneticCodes.groovy builds this)")
    ap.add_argument("--default_target_genetic_code", type=int, default=1,
                    help="NCBI genetic code of a target genome not named there (default 1)")
    args = ap.parse_args()
    set_genetic_code(genetic_code_for(args.target_genome,
                                      parse_genetic_code_spec(args.target_genetic_codes),
                                      args.default_target_genetic_code))

    # Filter scores.tsv to blocks needing rescue.
    rescue_blocks = []
    with open(args.scores) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("region_class") or "") != "goi_missing_but_strong_synteny":
                continue
            try:
                rescue_blocks.append({
                    "chrom": row["chrom"],
                    "start": int(float(row["start"])),
                    "end": int(float(row["end"])),
                    "name": row.get("region_name", ""),
                })
            except (KeyError, ValueError):
                continue

    if not rescue_blocks:
        # No rescue needed — emit an empty (header-only) GFF so the downstream
        # collect() always has a valid path. Don't fail the channel on empty.
        with open(args.output, "w") as out:
            out.write("##gff-version 3\n")
            out.write(f"# rescue_strong_synteny: 0 blocks to rescue for {args.genome_name}\n")
        return

    query_id, query_seq = _read_query(args.query)
    if not query_seq:
        print(f"[rescue] ERROR: empty query in {args.query}", file=sys.stderr)
        sys.exit(1)

    out_lines = ["##gff-version 3",
                 f"# rescue_strong_synteny: {len(rescue_blocks)} block(s) "
                 f"for {args.genome_name} against query {query_id}"]

    with tempfile.TemporaryDirectory(prefix="synvoy_rescue_") as tmp:
        query_fa = os.path.join(tmp, "query.faa")
        with open(query_fa, "w") as fh:
            fh.write(f">{query_id}\n{query_seq}\n")

        for idx, block in enumerate(rescue_blocks, start=1):
            chrom = block["chrom"]
            ws = max(1, block["start"] - args.window_pad)
            we = block["end"] + args.window_pad
            try:
                subseq = _read_fasta_subseq(args.samtools_bin, args.target_genome,
                                            chrom, ws, we)
            except subprocess.CalledProcessError as e:
                print(f"[rescue] samtools faidx failed for {chrom}:{ws}-{we}: "
                      f"{e.stderr[:200] if e.stderr else e}", file=sys.stderr)
                continue
            if not subseq:
                continue

            window_fa = os.path.join(tmp, f"win_{idx}.fna")
            with open(window_fa, "w") as fh:
                fh.write(f">{chrom}_{ws}_{we}\n{subseq}\n")

            gff_text = _run_miniprot_relaxed(
                window_fa, query_fa, args.miniprot_outc, args.miniprot_timeout,
                "strong_synteny_rescue_miniprot_timeout",
                f"{args.genome_name} block {idx} ({chrom}:{ws}-{we}, {len(subseq):,} bp)",
            )
            if not gff_text.strip():
                out_lines.append(
                    f"# block {idx} ({chrom}:{ws}-{we}) — relaxed miniprot found no model"
                )
                continue

            for mrna, cds_rows in _parse_miniprot_gff(gff_text):
                # Same terminal refinement as the main path (CLAUDE.md §21).
                _refine_rescue_termini(subseq, mrna, cds_rows)
                parent_id = (
                    f"GOI_rescue_{args.genome_name}_locus_{block['name'] or idx}_{idx}"
                )
                out_lines.extend(_build_rescue_gff_rows(
                    mrna, cds_rows, chrom, ws, query_id, parent_id,
                ))

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    # A miniprot run beyond its guard exits STEP_FAILED_EXIT (3), which the module turns into
    # the run's end; it ignores the other failures of the rescue (errorStrategy).
    exit_step_failed(main)
