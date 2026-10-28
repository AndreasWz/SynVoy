#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

try:
    from sequence_utils import (
        load_genome, write_fasta, parse_fasta,
        parse_bed as _parse_bed_dicts, str2bool, set_genetic_code,
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        load_genome, write_fasta, parse_fasta,
        parse_bed as _parse_bed_dicts, str2bool, set_genetic_code,
    )

try:
    from gene_predictor import (AugustusModelNotInstalled, check_predictor_available, predict_genes,
                                resolve_augustus_species)
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gene_predictor import (AugustusModelNotInstalled, check_predictor_available, predict_genes,
                                resolve_augustus_species)


def parse_bed(path):
    """Local tuple-returning view (chrom, start, end) over the shared BED parser.

    Drops zero/negative-length rows, which is the contract merge_regions and
    prodigal expect.
    """
    return [
        (r["chrom"], r["start"], r["end"])
        for r in _parse_bed_dicts(path)
        if r["end"] > r["start"]
    ]


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


def _with_parents(parts):
    """A CDS row without a Parent (Prodigal writes one per gene) as gene, mRNA and CDS rows
    (batch HG-5a): readers that group CDS by Parent, gff_to_faa.py among them, dropped it.
    The gene keeps the CDS's id, so a gene is named as before."""
    attrs = [a for a in parts[8].split(";") if a]
    cid = next((a.split("=", 1)[1] for a in attrs if a.startswith("ID=")), "")
    if not cid:
        cid = f"{parts[0]}_{parts[3]}_{parts[4]}"
    rest = [a for a in attrs if not a.startswith("ID=")]
    gene = parts[:2] + ["gene", parts[3], parts[4], ".", parts[6], ".", f"ID={cid}"]
    mrna = parts[:2] + ["mRNA", parts[3], parts[4], ".", parts[6], ".", f"ID={cid}.t1;Parent={cid}"]
    cds = parts[:8] + [";".join([f"ID={cid}.cds", f"Parent={cid}.t1"] + rest)]
    return [gene, mrna, cds]


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
                if parts[2] == "CDS" and not any(a.startswith("Parent=") for a in parts[8].split(";")):
                    for row in _with_parents(parts):
                        out.write("\t".join(row) + "\n")
                    continue
                out.write("\t".join(parts) + "\n")


def gene_regions(in_gff):
    """{gene id: region sequence id} of the gene rows of a predictor's GFF (region-local)."""
    out = {}
    if not os.path.exists(in_gff):
        return out
    with open(in_gff) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                gid = next((a.split("=", 1)[1] for a in parts[8].split(";") if a.startswith("ID=")), "")
                if gid:
                    out[gid] = parts[0]
    return out


def adjust_proteins(in_faa, out_faa, id_map, gene_region=None):
    """The predicted proteins with global coordinates and region-prefixed ids. A Prodigal id
    names its region (``<region>_<n>``); an AUGUSTUS id (``g1``) does not, so its region comes
    from the gene rows of the GFF (``gene_region``, gene_regions()); before (batch HG-5a) such
    proteins kept region-local coordinates and the same ids in every region."""
    if not os.path.exists(in_faa) or os.path.getsize(in_faa) == 0:
        open(out_faa, "w").close()
        return
    records = []
    for raw_header, _, seq in parse_fasta(in_faa):
        parts = raw_header.split(" # ")
        if len(parts) >= 4:
            token = parts[0]
            base_seqid = token.rsplit("_", 1)[0]
            if base_seqid not in id_map and (gene_region or {}).get(token) in id_map:
                base_seqid = gene_region[token]
                parts[0] = f"{base_seqid}_{token}"
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


def run_gene_prediction(fasta_in, faa_out, gff_out, predictor="auto",
                        augustus_species="fly"):
    """Run gene prediction using the unified gene_predictor dispatch."""
    success = predict_genes(
        fasta_in, faa_out, gff_out,
        predictor=predictor,
        augustus_species=augustus_species,
    )
    if not success:
        raise RuntimeError(
            f"Gene prediction failed (predictor={predictor}). "
            f"Check that Augustus or Prodigal is installed."
        )


def run_prodigal(fasta_in, faa_out, gff_out):
    """Legacy wrapper — dispatches through unified gene predictor."""
    run_gene_prediction(fasta_in, faa_out, gff_out, predictor="prodigal")


def main():
    parser = argparse.ArgumentParser(
        description="Run gene prediction (Augustus/Prodigal) on GOI-flanking regions"
    )
    parser.add_argument("--genome", required=True, help="Genome FASTA")
    parser.add_argument("--goi_bed", required=True, help="GOI bed file (hits)")
    parser.add_argument("--window", type=int, default=50000, help="Flanking window around hits")
    parser.add_argument("--output_faa", required=True, help="Output proteins FASTA")
    parser.add_argument("--output_gff", required=True, help="Output GFF")
    parser.add_argument("--fallback_full_genome", type=str2bool, default=False,
                        help="Fallback to full-genome prediction if no regions found")
    parser.add_argument("--gene_predictor", type=str, default="auto",
                        help="Gene predictor: auto, augustus, or prodigal (default: auto)")
    parser.add_argument("--augustus_species", type=str, default="fly",
                        help="Augustus species model (default: fly)")
    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the genome (default 1, the standard code)")
    args = parser.parse_args()
    set_genetic_code(args.genetic_code)
    # An Augustus model that is not installed stops here (fix-phase item 16): Augustus exited 1
    # and the prediction fell back to Prodigal, which cannot predict genes with introns.
    if args.gene_predictor != "prodigal":
        try:
            resolve_augustus_species(args.augustus_species)
        except AugustusModelNotInstalled as exc:
            sys.exit(f"ERROR: {exc}")

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
        print("No GOI regions found; skipping gene prediction.", file=sys.stderr)
        return

    if not regions:
        open(args.output_faa, "w").close()
        open(args.output_gff, "w").close()
        print("No valid regions after merging; skipping gene prediction.", file=sys.stderr)
        return

    tmp_dir = tempfile.mkdtemp(prefix="synvoy_prodigal_")
    try:
        regions_fa = os.path.join(tmp_dir, "regions.fna")
        tmp_faa = os.path.join(tmp_dir, "prodigal.faa")
        tmp_gff = os.path.join(tmp_dir, "prodigal.gff")

        id_map = write_regions_fasta(regions, genome_seqs, regions_fa)
        if not id_map:
            open(args.output_faa, "w").close()
            open(args.output_gff, "w").close()
            print("No regions produced for gene prediction.", file=sys.stderr)
            return

        run_gene_prediction(
            regions_fa, tmp_faa, tmp_gff,
            predictor=args.gene_predictor,
            augustus_species=args.augustus_species,
        )

        adjust_proteins(tmp_faa, args.output_faa, id_map, gene_regions(tmp_gff))
        adjust_gff(tmp_gff, args.output_gff, id_map)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
