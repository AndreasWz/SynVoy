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
                          timeout: int) -> str:
    """Run miniprot --gff with very permissive thresholds against a single window.

    Mirrors the relaxed search profile used elsewhere in iterative_search_runner
    (-n 2 / -p 0.2 / -N 50 / --outs 0.2 / --outc <min_outc>). Returns the raw GFF3
    text; the caller parses it.
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
        window_fa, query_fa,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
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


def _get_attr(attr_str: str, key: str, default: str = "") -> str:
    for kv in attr_str.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            if k.strip() == key:
                return v.strip()
    return default


def _identity_from_mrna_attrs(attrs: str) -> float:
    raw = _get_attr(attrs, "Identity", "")
    try:
        return float(raw)
    except ValueError:
        return 0.0


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
    ap.add_argument("--miniprot_timeout", type=int, default=120)
    args = ap.parse_args()

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
            )
            if not gff_text.strip():
                out_lines.append(
                    f"# block {idx} ({chrom}:{ws}-{we}) — relaxed miniprot found no model"
                )
                continue

            for mrna, cds_rows in _parse_miniprot_gff(gff_text):
                parent_id = (
                    f"GOI_rescue_{args.genome_name}_locus_{block['name'] or idx}_{idx}"
                )
                out_lines.extend(_build_rescue_gff_rows(
                    mrna, cds_rows, chrom, ws, query_id, parent_id,
                ))

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
