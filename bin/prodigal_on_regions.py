#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

try:
    from sequence_utils import load_genome, write_fasta, parse_fasta
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import load_genome, write_fasta, parse_fasta


def str2bool(value: str) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    val = value.strip().lower()
    if val in {"true", "1", "yes", "y", "t"}:
        return True
    if val in {"false", "0", "no", "n", "f"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def parse_bed(path):
    regions = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return regions
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                if end > start:
                    regions.append((chrom, start, end))
            except ValueError:
                continue
    return regions


def merge_regions(regions, genome_seqs, window):
    by_chrom = {}
    for chrom, start, end in regions:
        if chrom not in genome_seqs:
            continue
        slen = len(genome_seqs[chrom])
        w_start = max(0, start - window)
        w_end = min(slen, end + window)
        if w_end <= w_start:
            continue
        by_chrom.setdefault(chrom, []).append((w_start, w_end))

    merged = []
    for chrom, spans in by_chrom.items():
        spans.sort(key=lambda x: x[0])
        cur_start, cur_end = spans[0]
        for s, e in spans[1:]:
            if s <= cur_end:
                cur_end = max(cur_end, e)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((chrom, cur_start, cur_end))
    return merged


def write_regions_fasta(regions, genome_seqs, out_fa):
    records = []
    id_map = {}
    for chrom, start, end in regions:
        seq = genome_seqs[chrom][start:end]
        if not seq:
            continue
        seq_id = f"{chrom}__{start}__{end}"
        records.append((seq_id, seq))
        id_map[seq_id] = (chrom, start, end)
    if records:
        write_fasta(records, out_fa)
    else:
        open(out_fa, "w").close()
    return id_map


def prefix_attrs(attr_str, seq_id):
    if not attr_str:
        return attr_str
    out = []
    for attr in attr_str.split(";"):
        if "=" in attr:
            k, v = attr.split("=", 1)
            if k in ("ID", "Parent"):
                v = ",".join(f"{seq_id}_{item}" for item in v.split(",") if item)
            out.append(f"{k}={v}")
        else:
            out.append(attr)
    return ";".join(out)


def adjust_gff(in_gff, out_gff, id_map):
    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        if not os.path.exists(in_gff) or os.path.getsize(in_gff) == 0:
            return
        with open(in_gff, "r") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seq_id = parts[0]
                if seq_id not in id_map:
                    continue
                chrom, offset, _ = id_map[seq_id]
                try:
                    parts[3] = str(int(parts[3]) + offset)
                    parts[4] = str(int(parts[4]) + offset)
                except ValueError:
                    pass
                parts[0] = chrom
                parts[8] = prefix_attrs(parts[8], seq_id)
                out.write("\t".join(parts) + "\n")


def adjust_proteins(in_faa, out_faa, id_map):
    if not os.path.exists(in_faa) or os.path.getsize(in_faa) == 0:
        open(out_faa, "w").close()
        return
    records = []
    for raw_header, _, seq in parse_fasta(in_faa):
        parts = raw_header.split(" # ")
        if len(parts) >= 4:
            token = parts[0]
            base_seqid = token.rsplit("_", 1)[0]
            mapping = id_map.get(base_seqid)
            if mapping:
                _, offset, _ = mapping
                try:
                    parts[1] = str(int(parts[1]) + offset)
                    parts[2] = str(int(parts[2]) + offset)
                except ValueError:
                    pass
                if len(parts) >= 5:
                    parts[4] = prefix_attrs(parts[4], base_seqid)
                raw_header = " # ".join(parts)
        records.append((raw_header, seq))
    write_fasta(records, out_faa)


def run_prodigal(fasta_in, faa_out, gff_out):
    cmd = [
        "prodigal",
        "-i", fasta_in,
        "-a", faa_out,
        "-f", "gff",
        "-o", gff_out,
        "-p", "meta",
        "-q"
    ]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(description="Run Prodigal on GOI-flanking regions only")
    parser.add_argument("--genome", required=True, help="Genome FASTA")
    parser.add_argument("--goi_bed", required=True, help="GOI bed file (hits)")
    parser.add_argument("--window", type=int, default=50000, help="Flanking window around hits")
    parser.add_argument("--output_faa", required=True, help="Output proteins FASTA")
    parser.add_argument("--output_gff", required=True, help="Output GFF")
    parser.add_argument("--fallback_full_genome", type=str2bool, default=False,
                        help="Fallback to full-genome Prodigal if no regions found")
    args = parser.parse_args()

    genome_seqs = load_genome(args.genome)
    regions = parse_bed(args.goi_bed)

    if regions:
        regions = merge_regions(regions, genome_seqs, max(0, args.window))
    elif args.fallback_full_genome:
        regions = [(chrom, 0, len(seq)) for chrom, seq in genome_seqs.items()]
    else:
        # No GOI regions and fallback disabled: produce empty outputs
        open(args.output_faa, "w").close()
        open(args.output_gff, "w").close()
        print("No GOI regions found; skipping Prodigal.", file=sys.stderr)
        return

    if not regions:
        open(args.output_faa, "w").close()
        open(args.output_gff, "w").close()
        print("No valid regions after merging; skipping Prodigal.", file=sys.stderr)
        return

    tmp_dir = tempfile.mkdtemp(prefix="synterra_prodigal_")
    try:
        regions_fa = os.path.join(tmp_dir, "regions.fna")
        tmp_faa = os.path.join(tmp_dir, "prodigal.faa")
        tmp_gff = os.path.join(tmp_dir, "prodigal.gff")

        id_map = write_regions_fasta(regions, genome_seqs, regions_fa)
        if not id_map:
            open(args.output_faa, "w").close()
            open(args.output_gff, "w").close()
            print("No regions produced for Prodigal.", file=sys.stderr)
            return

        run_prodigal(regions_fa, tmp_faa, tmp_gff)

        adjust_proteins(tmp_faa, args.output_faa, id_map)
        adjust_gff(tmp_gff, args.output_gff, id_map)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
