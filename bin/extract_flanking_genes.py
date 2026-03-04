#!/usr/bin/env python3

import argparse
import sys
import os
import subprocess
import re

# Use our own sequence utilities (no BioPython dependency)
try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, reverse_complement, translate, parse_gff as parse_gff_base
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, reverse_complement, translate, parse_gff as parse_gff_base
    )

def parse_gff_for_genes(gff_file):
    """
    Parse GFF3 file into a list of gene dictionaries.
    Standardizes to BED Coordinates (0-based start, 1-based end, half-open).
    
    Detailed Coordinate Handling:
    - GFF3 Standard: 1-based, closed interval [start, end].
    - Python/BED Standard: 0-based, half-open interval [start, end).
    
    Conversion:
    - BED Start = GFF Start - 1
    - BED End   = GFF End (unchanged, as Python slice includes element at (End-1))
    """
    genes = []
    cds_by_parent = {}
    cds_orphans = []
    parent_map = {} # Transcript -> Gene
    transcript_attrs = {} # Transcript ID -> attributes
    
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
                     transcript_attrs[tid] = attributes
                     # Check if Parent is present, sometimes might be gene ID itself or absent
                     if gid:
                         parent_map[tid] = gid
                     else:
                         # Orphan transcript?
                         pass
            
            elif feature_type == 'CDS':
                cds_entry = {
                    'chrom': parts[0],
                    'start': int(parts[3]) - 1, # 0-based
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'phase': parts[7]
                }
                pid = attributes.get('Parent')
                if pid:
                    if pid not in cds_by_parent:
                        cds_by_parent[pid] = []
                    cds_by_parent[pid].append(cds_entry)
                else:
                    # Prodigal-style CDS-only GFFs may have no Parent/gene features.
                    # Keep these so we can build pseudo genes below.
                    cds_orphans.append((attributes.get('ID', ''), cds_entry))
                    
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
            gene['transcript_attrs'] = {}
        else:
            # Pick first transcript (or longest)
            # Just take the first one for now
            # TODO: Improve isoform selection
            best_t = transcripts[0]
            gene['transcript_attrs'] = transcript_attrs.get(best_t, {})
            if best_t in cds_by_parent:
                gene['cds_parts'] = cds_by_parent[best_t]
        
        processed_genes.append(gene)

    # Fallback for CDS-only annotation files (e.g., Prodigal output):
    # synthesize gene-like records from CDS spans so downstream flanking logic
    # can still select neighboring models around the GOI region.
    if not processed_genes and (cds_by_parent or cds_orphans):
        print("GFF contains no gene features; falling back to CDS-derived pseudo genes")

        pseudo_genes = []

        # Parent-grouped CDS blocks (e.g., transcript IDs).
        for pid, parts_list in cds_by_parent.items():
            parts_sorted = sorted(parts_list, key=lambda x: x['start'])
            if not parts_sorted:
                continue
            chrom = parts_sorted[0]['chrom']
            strand = parts_sorted[0]['strand']
            start = min(p['start'] for p in parts_sorted)
            end = max(p['end'] for p in parts_sorted)
            gid = pid if pid else f"cds_group_{len(pseudo_genes)+1}"
            pseudo_genes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'type': 'CDS',
                'attrs': {'ID': gid},
                'id': gid,
                'cds_parts': parts_sorted
            })

        # Orphan CDS entries (no Parent).
        for idx, (oid, cds_entry) in enumerate(cds_orphans, start=1):
            gid = oid if oid else f"cds_orphan_{idx}"
            pseudo_genes.append({
                'chrom': cds_entry['chrom'],
                'start': cds_entry['start'],
                'end': cds_entry['end'],
                'strand': cds_entry['strand'],
                'type': 'CDS',
                'attrs': {'ID': gid},
                'id': gid,
                'cds_parts': [cds_entry]
            })

        processed_genes = pseudo_genes

    return sorted(processed_genes, key=lambda x: (x['chrom'], x['start']))

def _is_generic_id_label(label):
    if not label:
        return True
    txt = str(label).strip()
    if not txt:
        return True
    if re.match(r'^(gene-)?[A-Za-z]{1,8}\d*_\d+$', txt):
        return True
    if re.match(r'^LOC\d+$', txt, re.IGNORECASE):
        return True
    return False

def _is_noninformative_product(product):
    if not product:
        return True
    txt = str(product).strip().lower()
    if not txt:
        return True
    return txt in {
        'hypothetical protein',
        'uncharacterized protein',
        'predicted protein',
    }

def _preferred_gene_label(gene):
    """
    Choose a human-readable label while keeping stable IDs as primary keys.
    Priority: gene symbol/name -> informative product -> locus_tag/ID.
    """
    attrs = gene.get('attrs', {}) or {}
    tattrs = gene.get('transcript_attrs', {}) or {}

    for source in (tattrs, attrs):
        for key in ('gene', 'Name'):
            cand = str(source.get(key, '')).strip()
            if cand and not _is_generic_id_label(cand):
                return cand

    for source in (tattrs, attrs):
        product = str(source.get('product', '')).strip()
        if product and not _is_noninformative_product(product):
            return product

    for source in (tattrs, attrs):
        for key in ('locus_tag', 'Name', 'ID'):
            cand = str(source.get(key, '')).strip()
            if cand:
                return cand

    return str(gene.get('id', '')).strip()

# load_genome now imported from sequence_utils

def main():
    parser = argparse.ArgumentParser(description="Extract flanking genes")
    parser.add_argument("--bed", required=True, help="Input BED file with gene location")
    parser.add_argument("--gff", required=True, help="Genome GFF file or 'NO_GFF'")
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--n_flank", type=int, default=10, help="Number of flanking genes")
    parser.add_argument("--min_size", type=int, default=500, help="Min gene size")
    parser.add_argument("--prefer_large", type=str, default="true", help="Prefer large genes")
    parser.add_argument("--exon_mode", type=str, default="false", 
                        help="If true, output individual exon CDS sequences instead of full protein")
    parser.add_argument("--pred_flank_window", type=int, default=50000,
                        help="Flanking window around GOI hits for Prodigal prediction")
    parser.add_argument("--pred_keep_pct", type=float, default=0.10,
                        help="Fraction of longest Prodigal predictions to keep")
    parser.add_argument("--out_bed", required=True, help="Output BED")
    parser.add_argument("--out_faa", required=True, help="Output FASTA")
    
    args = parser.parse_args()
    prefer_large = args.prefer_large.lower() == 'true'
    exon_mode = args.exon_mode.lower() == 'true'
    
    if exon_mode:
        print("Exon mode enabled: extracting individual CDS exon sequences")

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
        FLANK_WINDOW = max(0, int(args.pred_flank_window))
        
        for region in target_regions:
            chrom = region['chrom']
            if chrom not in genome_seqs:
                continue
                
            slen = len(genome_seqs[chrom])
            
            # Define window
            center = (region['start'] + region['end']) // 2
            w_start = max(0, center - FLANK_WINDOW)
            w_end = min(slen, center + FLANK_WINDOW)
            
            # Extract sequence - genome_seqs is now dict of strings
            subseq = genome_seqs[chrom][w_start:w_end]
            sub_id = f"{chrom}_{w_start}_{w_end}"
            
            # Write temp fasta using our utility
            temp_fa = f"temp_{sub_id}.fasta"
            write_fasta([(sub_id, subseq)], temp_fa)
                
            # Run Prodigal
            # prodigal -i inputs.fna -a proteins.faa -o coords.gff -p meta
            temp_out_faa = f"temp_{sub_id}.faa"
            cmd = ["prodigal", "-i", temp_fa, "-a", temp_out_faa, "-p", "meta", "-q"]
            
            try:
                subprocess.run(cmd, check=True)
                
                # Parse output FAA using our utility
                for header, clean_id, seq in parse_fasta(temp_out_faa):
                    # Prodigal header: id_1 # start # end # strand # ...
                    # We need to map back to genomic coordinates
                    parts = header.split(" # ")
                    if len(parts) >= 4:
                        local_start = int(parts[1]) # 1-based
                        local_end = int(parts[2])   # 1-based inclusive
                        strand_code = parts[3] # 1 or -1
                        strand = "+" if strand_code == "1" else "-"
                        
                        # Map to global and convert to 0-based BED
                        # global_start (0-based) = w_start (0-based) + (local_start - 1)
                        global_start = w_start + (local_start - 1)
                        # global_end (1-based half-open) = w_start (0-based) + local_end
                        global_end = w_start + local_end
                        
                        extracted_genes.append({
                            'chrom': chrom,
                            'start': global_start,
                            'end': global_end,
                            'strand': strand,
                            'attrs': {'ID': f"pred_{chrom}_{global_start}"},
                            'seq': seq  # Store seq directly (string)
                        })
                        
            except Exception as e:
                print(f"Gene prediction failed: {e}")
            finally:
                if os.path.exists(temp_fa): os.remove(temp_fa)
                if os.path.exists(temp_out_faa): os.remove(temp_out_faa)

        # FILTER: Keep only top 10% longest Prodigal predictions
        # Prodigal on eukaryotic genomes greatly overpredicts tiny ORFs,
        # drowning out real genes and inflating synteny denominators.
        if extracted_genes:
            original_count = len(extracted_genes)
            extracted_genes.sort(key=lambda g: g['end'] - g['start'], reverse=True)
            keep_fraction = args.pred_keep_pct
            if keep_fraction <= 0 or keep_fraction > 1:
                keep_fraction = 0.10
            keep_count = max(args.n_flank * 2 + 1, int(len(extracted_genes) * keep_fraction))
            kept_genes = set()
            for g in extracted_genes[:keep_count]:
                kept_genes.add((g['chrom'], g['start'], g['end']))
            # CRITICAL: Always re-inject genes overlapping the GOI hit region
            # The GOI itself (e.g. a tiny 67aa peptide) may be filtered out by
            # the top-10% size filter. Without it, the home track won't show GOI.
            for g in extracted_genes:
                for region in target_regions:
                    if (g['chrom'] == region['chrom'] and
                        g['start'] < region['end'] and
                        g['end'] > region['start']):
                        kept_genes.add((g['chrom'], g['start'], g['end']))
                        break
            extracted_genes = [g for g in extracted_genes
                              if (g['chrom'], g['start'], g['end']) in kept_genes]
            # FALLBACK: If no predicted gene overlaps a target region,
            # inject the target region itself as a GOI pseudo-gene.
            # This handles cases where Prodigal can't predict the GOI at all.
            for region in target_regions:
                has_overlap = any(
                    g['chrom'] == region['chrom'] and
                    g['start'] < region['end'] and
                    g['end'] > region['start']
                    for g in extracted_genes
                )
                if not has_overlap and region['chrom'] in genome_seqs:
                    goi_start = region['start']
                    goi_end = region['end']
                    goi_dna = genome_seqs[region['chrom']][goi_start:goi_end]
                    goi_prot = translate(goi_dna)
                    if '*' in goi_prot:
                        goi_prot = goi_prot.split('*')[0]
                    goi_id = f"GOI_{region['chrom']}_{goi_start}"
                    extracted_genes.append({
                        'chrom': region['chrom'],
                        'start': goi_start,
                        'end': goi_end,
                        'strand': '+',
                        'attrs': {'ID': goi_id},
                        'seq': goi_prot
                    })
                    print(f"Injected GOI pseudo-gene at {region['chrom']}:{goi_start}-{goi_end}")
            # Re-sort by genomic position for proper flanking order
            extracted_genes.sort(key=lambda g: (g['chrom'], g['start']))
            print(f"Prodigal filter: kept {len(extracted_genes)} genes from {original_count} predictions (top {keep_count} longest + GOI overlaps)")

    else:
        # EXISTING LOGIC FOR GFF
        all_genes = parse_gff_for_genes(args.gff)
        
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
            
            # Add genes from window
            window_genes = chrom_genes[start_idx:end_idx]
            extracted_genes.extend(window_genes)
            
            # Ensure genes overlapping the target region are included
            overlap_genes = [
                g for g in chrom_genes
                if g['start'] < region['end'] and g['end'] > region['start']
            ]
            for g in overlap_genes:
                if g not in extracted_genes:
                    extracted_genes.append(g)
            
            # If no overlap genes found, inject a GOI pseudo-gene from the region
            if not overlap_genes and region['chrom'] in genome_seqs:
                goi_start = region['start']
                goi_end = region['end']
                goi_dna = genome_seqs[region['chrom']][goi_start:goi_end]
                goi_prot = translate(goi_dna)
                if '*' in goi_prot:
                    goi_prot = goi_prot.split('*')[0]
                goi_id = f"GOI_{region['chrom']}_{goi_start}"
                extracted_genes.append({
                    'chrom': region['chrom'],
                    'start': goi_start,
                    'end': goi_end,
                    'strand': '+',
                    'attrs': {'ID': goi_id},
                    'seq': goi_prot
                })
                print(f"Injected GOI pseudo-gene at {region['chrom']}:{goi_start}-{goi_end}")

    # Write Outputs
    all_fasta_records = []  # Collect all FASTA records
    
    with open(args.out_bed, 'w') as bed_out:
        seen = set()
        for gene in extracted_genes:
            raw_id = gene['attrs'].get('ID', '')
            # Include chromosome and position in dedup key so genes with same ID
            # on different contigs are not silently dropped
            gid = raw_id if raw_id else f"{gene['chrom']}_{gene['start']}"
            dedup_key = f"{gene['chrom']}:{gene['start']}-{gene['end']}:{gid}"
            if dedup_key in seen: continue
            seen.add(dedup_key)
            display_label = _preferred_gene_label(gene).replace('\t', ' ').strip() or gid
            seq_header = gid if display_label == gid else f"{gid} label={display_label}"
            
            # Write BED (always the full gene)
            bed_out.write(
                f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t{gid}\t.\t{gene['strand']}\t{display_label}\n"
            )
            
            # Write FASTA - depends on exon_mode
            if 'seq' in gene:
                # From prediction - no exon info available, use whole sequence
                prot_seq = str(gene['seq'])
                if len(prot_seq) * 3 >= args.min_size:
                    all_fasta_records.append((seq_header, prot_seq))
            else:
                if gene['chrom'] not in genome_seqs:
                    continue
                    
                seq_record = genome_seqs[gene['chrom']]  # Now a string, not SeqRecord
                
                if gene.get('cds_parts'):
                    # We have CDS parts (exons)
                    cds_parts = sorted(gene['cds_parts'], key=lambda x: x['start'])
                    
                    if exon_mode:
                        # EXON MODE: Output each exon as separate sequence
                        # ID format: gene_id|exon_N (1-indexed)
                        for exon_idx, part in enumerate(cds_parts, start=1):
                            exon_id = f"{gid}|exon_{exon_idx}"
                            
                            # Extract exon DNA sequence
                            exon_dna = seq_record[part['start']:part['end']]
                            
                            # Handle strand for this exon
                            if gene['strand'] == '-':
                                exon_dna = reverse_complement(exon_dna)
                            
                            # Translate with correct frame (use phase from GFF if available)
                            phase = int(part.get('phase', 0)) if part.get('phase', '.') != '.' else 0
                            
                            # Skip phase nucleotides at start
                            if phase > 0:
                                exon_dna = exon_dna[phase:]
                            
                            # Pad if needed for translation
                            remainder = len(exon_dna) % 3
                            if remainder:
                                exon_dna = exon_dna[:-remainder]
                            
                            if len(exon_dna) < 9:  # Skip very short exons (< 3 amino acids)
                                continue
                                
                            exon_prot = translate(exon_dna)
                            # Remove stop codons
                            exon_prot = exon_prot.replace('*', '')
                            
                            # Write exon sequence with metadata in header
                            exon_header = f"{exon_id} parent={gid} exon={exon_idx}/{len(cds_parts)} coords={part['start']}-{part['end']} strand={gene['strand']}"
                            all_fasta_records.append((exon_header, exon_prot))

                        # ALSO emit the full-length protein (CDS DNA
                        # concatenated in genomic order, then translated).
                        # Per-exon records drive sensitive MMseqs search;
                        # the full-length protein is used by miniprot for
                        # multi-exon gene modeling in target genomes.
                        full_dna = ""
                        for part in cds_parts:
                            full_dna += seq_record[part['start']:part['end']]
                        if gene['strand'] == '-':
                            full_dna = reverse_complement(full_dna)
                        remainder = len(full_dna) % 3
                        if remainder:
                            full_dna = full_dna[:-remainder]
                        full_prot = translate(full_dna)
                        if '*' in full_prot:
                            full_prot = full_prot.split('*')[0]
                        if len(full_prot) >= 10:
                            full_header = f"{gid} full_length_protein exons={len(cds_parts)} strand={gene['strand']}"
                            all_fasta_records.append((full_header, full_prot))
                    else:
                        # WHOLE PROTEIN MODE: Concatenate all exons
                        dna_seq = ""
                        for part in cds_parts:
                            part_seq = seq_record[part['start']:part['end']]
                            dna_seq += part_seq
                        
                        if gene['strand'] == '-':
                            dna_seq = reverse_complement(dna_seq)
                            
                        remainder = len(dna_seq) % 3
                        if remainder:
                            dna_seq = dna_seq[:-remainder]
                            
                        prot_seq = translate(dna_seq)
                        # Stop at first stop codon
                        if '*' in prot_seq:
                            prot_seq = prot_seq.split('*')[0]
                        
                        if len(prot_seq) * 3 >= args.min_size:
                            all_fasta_records.append((seq_header, prot_seq))
                        
                else:
                    # NAIVE FALLBACK (No CDS parts, use full genomic)
                    feature_seq = seq_record[gene['start']:gene['end']]
                    if gene['strand'] == '-':
                        feature_seq = reverse_complement(feature_seq)
                    prot_seq = translate(feature_seq)
                    # Stop at first stop codon
                    if '*' in prot_seq:
                        prot_seq = prot_seq.split('*')[0]
                    
                    if len(prot_seq) * 3 >= args.min_size:
                        all_fasta_records.append((seq_header, prot_seq))
    
    # Write all FASTA records at once
    write_fasta(all_fasta_records, args.out_faa)

if __name__ == "__main__":
    main()
