#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import os
import subprocess

def parse_gff(gff_file):
    """
    Parse GFF3 file into a list of gene dictionaries.
    Standardizes to BED Coordinates (0-based start, 1-based end, half-open).
        GFF: 1-based start, 1-based end (closed)
        Internal: start - 1, end
    """
    genes = []
    cds_by_parent = {}
    parent_map = {} # Transcript -> Gene
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            feature_type = parts[2]
            
            # Parse attributes
            attributes = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attributes[k] = v
            
            if feature_type == 'gene':
                genes.append({
                    'chrom': parts[0],
                    'start': int(parts[3]) - 1, # Convert 1-based GFF to 0-based BED
                    'end': int(parts[4]),       # 1-based closed to 1-based half-open (unchanged)
                    'strand': parts[6],
                    'type': feature_type,
                    'attrs': attributes,
                    'id': attributes.get('ID')
                })
            
            elif feature_type in ['mRNA', 'transcript', 'mrna']:
                tid = attributes.get('ID')
                gid = attributes.get('Parent')
                if tid:
                     # Check if Parent is present, sometimes might be gene ID itself or absent
                     if gid:
                         parent_map[tid] = gid
                     else:
                         # Orphan transcript?
                         pass
            
            elif feature_type == 'CDS':
                pid = attributes.get('Parent')
                if pid:
                    if pid not in cds_by_parent:
                        cds_by_parent[pid] = []
                    cds_by_parent[pid].append({
                        'chrom': parts[0],
                        'start': int(parts[3]) - 1, # 0-based
                        'end': int(parts[4]),
                        'strand': parts[6],
                        'phase': parts[7]
                    })
                    
    # Map Gene -> List of CDS
    processed_genes = []
    for gene in genes:
        gid = gene['id']
        gene['cds_parts'] = []
        
        # Find transcripts for this gene
        transcripts = [t for t, g in parent_map.items() if g == gid]
        
        # If no transcripts (maybe direct CDS parent?), check CDS directly
        if not transcripts:
            if gid in cds_by_parent:
                gene['cds_parts'] = cds_by_parent[gid]
        else:
            # Pick first transcript (or longest)
            # Just take the first one for now
            # TODO: Improve isoform selection
            best_t = transcripts[0]
            if best_t in cds_by_parent:
                gene['cds_parts'] = cds_by_parent[best_t]
        
        processed_genes.append(gene)

    return sorted(processed_genes, key=lambda x: (x['chrom'], x['start']))

def load_genome(fasta_file):
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

def run_prediction(genome_file, chrom, start, end, output_faa):
    """
    Run Prodigal (or Augustus) on the specific region.
    1. Extract region to temp fasta.
    2. Run predictor.
    3. Return list of genes/proteins.
    """
    print(f"Running gene prediction on {chrom}:{start}-{end}...")
    
    # 1. Extract sequence
    # Use samtools faidx if available, or just python
    # We loaded genome in memory anyway for main script, but here we might want just a chunk.
    # Let's simple use the in-memory genome if mapped, or re-read.
    # The main function loads genome, so pass seq record?
    pass # logic moved to main
    
def main():
    parser = argparse.ArgumentParser(description="Extract flanking genes")
    parser.add_argument("--bed", required=True, help="Input BED file with gene location")
    parser.add_argument("--gff", required=True, help="Genome GFF file or 'NO_GFF'")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--n_flank", type=int, default=10, help="Number of flanking genes")
    parser.add_argument("--min_size", type=int, default=500, help="Min gene size")
    parser.add_argument("--prefer_large", type=str, default="true", help="Prefer large genes")
    parser.add_argument("--out_bed", required=True, help="Output BED")
    parser.add_argument("--out_faa", required=True, help="Output FASTA")
    
    args = parser.parse_args()
    prefer_large = args.prefer_large.lower() == 'true'

    # Load INPUT Genes
    target_regions = []
    with open(args.bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                target_regions.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })

    if not target_regions:
        print("No target regions found.")
        open(args.out_bed, 'w').close()
        open(args.out_faa, 'w').close()
        return

    genome_seqs = load_genome(args.genome)
    
    extracted_genes = []
    
    # MODE SWITCH
    if args.gff == "NO_GFF" or not os.path.exists(args.gff):
        print("No GFF provided. Running gene prediction on flanking regions...")
        
        # For each target region, extract a window (e.g. +/- 50kb)
        FLANK_WINDOW = 50000 
        
        for region in target_regions:
            chrom = region['chrom']
            if chrom not in genome_seqs:
                continue
                
            slen = len(genome_seqs[chrom])
            
            # Define window
            center = (region['start'] + region['end']) // 2
            w_start = max(0, center - FLANK_WINDOW)
            w_end = min(slen, center + FLANK_WINDOW)
            
            # Extract sequence
            subseq = genome_seqs[chrom].seq[w_start:w_end]
            sub_id = f"{chrom}_{w_start}_{w_end}"
            
            # Write temp fasta
            temp_fa = f"temp_{sub_id}.fasta"
            with open(temp_fa, 'w') as tf:
                SeqIO.write(SeqRecord(subseq, id=sub_id, description=""), tf, "fasta")
                
            # Run Prodigal
            # prodigal -i inputs.fna -a proteins.faa -o coords.gff -p meta
            temp_out_faa = f"temp_{sub_id}.faa"
            cmd = ["prodigal", "-i", temp_fa, "-a", temp_out_faa, "-p", "meta", "-q"]
            
            try:
                subprocess.run(cmd, check=True)
                
                # Parse output FAA
                for record in SeqIO.parse(temp_out_faa, "fasta"):
                    # Prodigal header: >id_1 # start # end # ...
                    # We need to map back to genomic coordinates
                    # Header: k12_1 # 34 # 456 # ...
                    parts = record.description.split(" # ")
                    if len(parts) >= 3:
                        local_start = int(parts[1]) # 1-based
                        local_end = int(parts[2])   # 1-based inclusive
                        strand_code = parts[3] # 1 or -1
                        strand = "+" if strand_code == "1" else "-"
                        
                        # Map to global and convert to 0-based BED
                        # global_start (0-based) = w_start (0-based) + (local_start - 1)
                        global_start = w_start + (local_start - 1)
                        # global_end (1-based half-open) = w_start (0-based) + local_end
                        # e.g. w=0. local=1..3. glob=0..3. len=3. Correct.
                        global_end = w_start + local_end
                        
                        extracted_genes.append({
                            'chrom': chrom,
                            'start': global_start,
                            'end': global_end,
                            'strand': strand,
                            'attrs': {'ID': f"pred_{chrom}_{global_start}"},
                            'seq': record.seq # Store seq directly
                        })
                        
            except Exception as e:
                print(f"Gene prediction failed: {e}")
            finally:
                if os.path.exists(temp_fa): os.remove(temp_fa)
                if os.path.exists(temp_out_faa): os.remove(temp_out_faa)

    else:
        # EXISTING LOGIC FOR GFF
        all_genes = parse_gff(args.gff)
        
        for region in target_regions:
            # Filter genes on same chrom
            chrom_genes = [g for g in all_genes if g['chrom'] == region['chrom']]
            
            # Find closest center
            center_idx = -1
            min_dist = float('inf')
            
            reg_center = (region['start'] + region['end']) / 2
            
            for i, gene in enumerate(chrom_genes):
                gene_center = (gene['start'] + gene['end']) / 2
                dist = abs(reg_center - gene_center)
                if dist < min_dist:
                    min_dist = dist
                    center_idx = i
            
            if center_idx == -1: continue
                
            # Window selection logic
            start_idx = max(0, center_idx - args.n_flank)
            end_idx = min(len(chrom_genes), center_idx + args.n_flank + 1)
            
            # Add genes
            for i in range(start_idx, end_idx):
                extracted_genes.append(chrom_genes[i])

    # Write Outputs
    with open(args.out_bed, 'w') as bed_out, open(args.out_faa, 'w') as faa_out:
        seen = set()
        for gene in extracted_genes:
            gid = gene['attrs'].get('ID', f"{gene['chrom']}_{gene['start']}")
            if gid in seen: continue
            seen.add(gid)
            
            # Write BED
            bed_out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{gid}\t.\t{gene['strand']}\n")
            
            # Write FASTA
            # If we already have seq (from prediction), use it. Else extract.
            if 'seq' in gene:
                prot_seq = gene['seq']
            else:
                if gene['chrom'] in genome_seqs:
                    seq_record = genome_seqs[gene['chrom']]
                    
                    if gene.get('cds_parts'):
                        # SPLICED EXTRACTION
                        cds_parts = sorted(gene['cds_parts'], key=lambda x: x['start'])
                        
                        # Concatenate sequence
                        dna_seq = Seq("")
                        for part in cds_parts:
                            # 0-based slicing: [start:end]
                            part_seq = seq_record.seq[part['start']:part['end']]
                            dna_seq += part_seq
                        
                        if gene['strand'] == '-':
                            dna_seq = dna_seq.reverse_complement()
                            
                        # Translate
                        # Pad if not multiple of 3?
                        remainder = len(dna_seq) % 3
                        if remainder:
                            dna_seq = dna_seq[:-remainder] # Truncate for now
                            
                        prot_seq = dna_seq.translate(to_stop=True)
                        
                    else:
                        # NAIVE FALLBACK (Genomic extraction)
                        # 0-based slicing
                        feature_seq = seq_record.seq[gene['start']:gene['end']]
                        if gene['strand'] == '-':
                            feature_seq = feature_seq.reverse_complement()
                        prot_seq = feature_seq.translate(to_stop=True)
                else:
                    continue

            if len(prot_seq) * 3 >= args.min_size: # Approximate check
                SeqIO.write(SeqRecord(prot_seq, id=gid, description=""), faa_out, "fasta")

if __name__ == "__main__":
    main()
