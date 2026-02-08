#!/usr/bin/env python3
"""
borrow_annotations.py — Borrow gene annotations from an annotated target genome

When the home genome lacks GFF annotations but a related target genome has them,
this script uses the annotated target's CDS to improve home genome annotation:

Flow:
  1. Load Prodigal-predicted proteins from home genome (de novo)
  2. Find target genomes with GFF annotations (CDS features)
  3. Search home proteins against annotated target → locate the syntenic region
  4. Extract CDS genes from that region in the target
  5. Search those CDS exon-wise against the home genome
  6. Annotate flanking genes around the GOI in the home genome

This gives higher-quality gene models than raw Prodigal, especially for
eukaryotic genomes where Prodigal misses intron-containing genes.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    from sequence_utils import (
        parse_fasta, write_fasta, load_genome, reverse_complement,
        translate, parse_gff, get_feature_id
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, load_genome, reverse_complement,
        translate, parse_gff, get_feature_id
    )


def find_annotated_target(genomes_dir):
    """
    Find a target genome in the genomes directory that has a matching GFF file
    with actual CDS features.
    
    Returns: (genome_fna, genome_gff) or (None, None)
    """
    genomes_path = Path(genomes_dir)
    
    fna_files = sorted(genomes_path.glob("*.fna"))
    
    best = None
    best_cds_count = 0
    
    for fna in fna_files:
        gff = fna.with_suffix('.gff')
        if not gff.exists():
            continue
        
        # Count CDS features to verify quality
        cds_count = 0
        gene_count = 0
        try:
            with open(gff) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        if parts[2] == 'CDS':
                            cds_count += 1
                        elif parts[2] == 'gene':
                            gene_count += 1
        except Exception:
            continue
        
        if cds_count >= 10:  # Minimum quality: at least 10 CDS features
            if cds_count > best_cds_count:
                best = (str(fna), str(gff))
                best_cds_count = cds_count
                print(f"  Candidate: {fna.name} — {cds_count} CDS, {gene_count} genes")
    
    if best:
        print(f"  Selected: {Path(best[0]).name} ({best_cds_count} CDS features)")
    else:
        print("  No annotated target genomes found")
    
    return best if best else (None, None)


def extract_region_cds(gff_file, genome_file, chrom, region_start, region_end,
                       flanking_bp=50000):
    """
    Extract CDS protein sequences from a GFF region.
    
    Returns: list of (gene_id, protein_seq) tuples
    """
    genome_seqs = load_genome(genome_file)
    
    if chrom not in genome_seqs:
        print(f"  Warning: {chrom} not in genome")
        return []
    
    chrom_seq = genome_seqs[chrom]
    
    # Expand region with flanking
    start = max(0, region_start - flanking_bp)
    end = min(len(chrom_seq), region_end + flanking_bp)
    
    # Parse GFF for genes in region
    genes_in_region = {}  # gene_id -> {'chrom', 'start', 'end', 'strand', 'cds_parts': []}
    gene_parents = {}  # child_id -> gene_id
    
    for feature in parse_gff(gff_file):
        f_chrom = feature.get('seqid', '')
        f_start = feature.get('start', 0)
        f_end = feature.get('end', 0)
        f_type = feature.get('type', '')
        attrs = feature.get('attributes', {})
        
        if f_chrom != chrom:
            continue
        
        # Check if feature overlaps our region
        if f_end < start or f_start > end:
            continue
        
        fid = get_feature_id(feature)
        
        if f_type in ('gene', 'mRNA', 'transcript'):
            parent = attrs.get('Parent', '')
            if f_type == 'gene':
                genes_in_region[fid] = {
                    'chrom': f_chrom, 'start': f_start, 'end': f_end,
                    'strand': feature.get('strand', '+'),
                    'cds_parts': [],
                    'name': attrs.get('Name', fid)
                }
            else:
                # mRNA/transcript — link to parent gene
                gene_parents[fid] = parent if parent else fid
        
        elif f_type == 'CDS':
            parent = attrs.get('Parent', '')
            # Find the gene this CDS belongs to
            gene_id = gene_parents.get(parent, parent)
            # If gene not in our dict, try direct parent
            if gene_id not in genes_in_region:
                # Create ad-hoc gene entry
                genes_in_region[gene_id] = {
                    'chrom': f_chrom, 'start': f_start, 'end': f_end,
                    'strand': feature.get('strand', '+'),
                    'cds_parts': [],
                    'name': attrs.get('Name', gene_id)
                }
            
            phase = feature.get('phase', '.')
            genes_in_region[gene_id]['cds_parts'].append({
                'start': f_start, 'end': f_end,
                'phase': phase
            })
            # Expand gene bounds
            genes_in_region[gene_id]['start'] = min(genes_in_region[gene_id]['start'], f_start)
            genes_in_region[gene_id]['end'] = max(genes_in_region[gene_id]['end'], f_end)
    
    # Extract protein sequences from CDS
    proteins = []
    for gid, gene in genes_in_region.items():
        if not gene['cds_parts']:
            continue
        
        cds_parts = sorted(gene['cds_parts'], key=lambda x: x['start'])
        
        # Concatenate CDS exons
        dna = ""
        for part in cds_parts:
            dna += chrom_seq[part['start']:part['end']]
        
        if gene['strand'] == '-':
            dna = reverse_complement(dna)
        
        # Trim to codon boundary
        dna = dna[:len(dna) - len(dna) % 3]
        if len(dna) < 30:
            continue
        
        prot = translate(dna)
        if '*' in prot:
            prot = prot.split('*')[0]
        
        if len(prot) < 10:
            continue
        
        proteins.append((gene['name'], prot))
    
    return proteins


def search_proteins_against_genome(query_faa, genome_fna, output_bed, output_gff,
                                    goi_bed, n_flanking=10, min_size=100):
    """
    Use MMseqs2 to search borrowed proteins against home genome.
    Filter results to keep genes near the GOI region.
    Output a GFF-like annotation and BED file for the flanking genes.
    """
    # Parse GOI location
    goi_chroms = {}
    goi_positions = []
    with open(goi_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                goi_chroms[chrom] = True
                goi_positions.append((chrom, start, end))
    
    if not goi_positions:
        print("  Warning: No GOI positions found")
        return
    
    goi_chrom = goi_positions[0][0]
    goi_center = (goi_positions[0][1] + goi_positions[0][2]) // 2
    
    # Run MMseqs2 search
    with tempfile.TemporaryDirectory() as tmpdir:
        hits_file = os.path.join(tmpdir, "hits.m8")
        
        cmd = [
            'mmseqs', 'easy-search',
            query_faa, genome_fna, hits_file, os.path.join(tmpdir, 'tmp'),
            '--search-type', '2',  # protein vs translated DNA
            '-s', '7.5',
            '-e', '0.001',
            '--format-output',
            'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"  Warning: MMseqs2 search failed: {e}")
            return
        
        if not os.path.exists(hits_file) or os.path.getsize(hits_file) == 0:
            print("  No hits found searching borrowed proteins against home genome")
            return
        
        # Parse hits
        hits = []
        with open(hits_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                try:
                    chrom = parts[1]
                    tstart = int(parts[8])
                    tend = int(parts[9])
                    strand = '+' if tstart < tend else '-'
                    gstart = min(tstart, tend) - 1
                    gend = max(tstart, tend)
                    
                    hits.append({
                        'gene_name': parts[0],
                        'chrom': chrom,
                        'start': gstart,
                        'end': gend,
                        'strand': strand,
                        'pident': float(parts[2]),
                        'evalue': float(parts[10]),
                        'bitscore': float(parts[11])
                    })
                except (ValueError, IndexError):
                    continue
        
        if not hits:
            print("  No parseable hits")
            return
        
        # Filter: only keep hits on the GOI chromosome and near GOI
        nearby_hits = [h for h in hits if h['chrom'] == goi_chrom]
        if not nearby_hits:
            # Fall back to best chromosome
            nearby_hits = hits
        
        # Deduplicate by gene name (keep best hit per gene)
        best_by_gene = {}
        for h in nearby_hits:
            gname = h['gene_name']
            if gname not in best_by_gene or h['evalue'] < best_by_gene[gname]['evalue']:
                best_by_gene[gname] = h
        
        # Sort by distance to GOI
        gene_list = sorted(best_by_gene.values(),
                          key=lambda h: abs((h['start'] + h['end']) // 2 - goi_center)
                          if h['chrom'] == goi_chrom else float('inf'))
        
        # Write outputs
        with open(output_bed, 'w') as bed_out, open(output_gff, 'w') as gff_out:
            gff_out.write("##gff-version 3\n")
            
            for gene in gene_list[:n_flanking * 4]:  # Take more than needed, let downstream filter
                bed_out.write(f"{gene['chrom']}\t{gene['start']}\t{gene['end']}\t"
                            f"{gene['gene_name']}\t{gene['pident']:.1f}\t{gene['strand']}\n")
                
                gff_out.write(f"{gene['chrom']}\tborrowed_annotation\tgene\t"
                            f"{gene['start'] + 1}\t{gene['end']}\t"
                            f"{gene['pident']:.1f}\t{gene['strand']}\t.\t"
                            f"ID={gene['gene_name']};Name={gene['gene_name']};"
                            f"Note=Borrowed from annotated target genome\n")
        
        print(f"  Wrote {min(len(gene_list), n_flanking * 4)} borrowed annotations")


def main():
    parser = argparse.ArgumentParser(
        description="Borrow gene annotations from annotated target genomes"
    )
    parser.add_argument("--home_genome", required=True,
                       help="Home genome FASTA")
    parser.add_argument("--home_proteins", required=True,
                       help="Prodigal-predicted proteins from home genome")
    parser.add_argument("--genomes_dir", required=True,
                       help="Directory containing target genomes (and optional GFFs)")
    parser.add_argument("--goi_bed", required=True,
                       help="BED file with GOI location in home genome")
    parser.add_argument("--output_gff", required=True,
                       help="Output GFF with borrowed annotations")
    parser.add_argument("--output_proteins", required=True,
                       help="Output FASTA with borrowed protein sequences")
    parser.add_argument("--n_flanking", type=int, default=10,
                       help="Number of flanking genes to consider")
    
    args = parser.parse_args()
    
    print("[borrow_annotations] Searching for annotated target genomes...")
    
    # Step 1: Find an annotated target genome
    target_fna, target_gff = find_annotated_target(args.genomes_dir)
    
    if not target_fna or not target_gff:
        print("[borrow_annotations] No annotated targets available. Using Prodigal only.")
        # Create empty outputs
        Path(args.output_gff).write_text("##gff-version 3\n")
        write_fasta([], args.output_proteins)
        return
    
    print(f"[borrow_annotations] Using annotated target: {Path(target_fna).name}")
    
    # Step 2: Search home proteins against target genome to find the region
    print("[borrow_annotations] Searching home proteins against annotated target...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        hits_file = os.path.join(tmpdir, "home_vs_target.m8")
        
        cmd = [
            'mmseqs', 'easy-search',
            args.home_proteins, target_fna, hits_file,
            os.path.join(tmpdir, 'tmp'),
            '--search-type', '2',
            '-s', '7.5',
            '-e', '0.001',
            '--format-output',
            'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"[borrow_annotations] MMseqs2 search failed: {e}")
            Path(args.output_gff).write_text("##gff-version 3\n")
            write_fasta([], args.output_proteins)
            return
        
        if not os.path.exists(hits_file) or os.path.getsize(hits_file) == 0:
            print("[borrow_annotations] No hits — can't locate region in target")
            Path(args.output_gff).write_text("##gff-version 3\n")
            write_fasta([], args.output_proteins)
            return
        
        # Parse: find best region in target genome
        hits_by_chrom = {}
        with open(hits_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                chrom = parts[1]
                tstart = int(parts[8])
                tend = int(parts[9])
                gstart = min(tstart, tend) - 1
                gend = max(tstart, tend)
                evalue = float(parts[10])
                
                if chrom not in hits_by_chrom:
                    hits_by_chrom[chrom] = []
                hits_by_chrom[chrom].append({
                    'start': gstart, 'end': gend, 'evalue': evalue
                })
        
        if not hits_by_chrom:
            print("[borrow_annotations] No parseable hits")
            Path(args.output_gff).write_text("##gff-version 3\n")
            write_fasta([], args.output_proteins)
            return
        
        # Pick chromosome with most hits
        best_chrom = max(hits_by_chrom.keys(),
                        key=lambda c: len(hits_by_chrom[c]))
        chrom_hits = hits_by_chrom[best_chrom]
        
        region_start = min(h['start'] for h in chrom_hits)
        region_end = max(h['end'] for h in chrom_hits)
        
        print(f"[borrow_annotations] Target region: {best_chrom}:{region_start}-{region_end} "
              f"({len(chrom_hits)} hits)")
    
    # Step 3: Extract CDS from target region
    print("[borrow_annotations] Extracting CDS from annotated target region...")
    proteins = extract_region_cds(target_gff, target_fna, best_chrom,
                                  region_start, region_end, flanking_bp=100000)
    
    if not proteins:
        print("[borrow_annotations] No CDS found in target region")
        Path(args.output_gff).write_text("##gff-version 3\n")
        write_fasta([], args.output_proteins)
        return
    
    print(f"[borrow_annotations] Extracted {len(proteins)} CDS from target region")
    
    # Step 4: Write borrowed proteins
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as tmp:
        borrowed_faa = tmp.name
        write_fasta(proteins, borrowed_faa)
    
    # Step 5: Search borrowed proteins against home genome
    print("[borrow_annotations] Searching borrowed CDS against home genome...")
    search_proteins_against_genome(
        borrowed_faa, args.home_genome, 
        args.output_proteins.replace('.faa', '.bed'),  # temp bed
        args.output_gff,
        args.goi_bed,
        n_flanking=args.n_flanking
    )
    
    # Write the borrowed proteins as output
    write_fasta(proteins, args.output_proteins)
    
    # Cleanup
    try:
        os.unlink(borrowed_faa)
    except Exception:
        pass
    
    print(f"[borrow_annotations] Done. Wrote {len(proteins)} borrowed annotations.")


if __name__ == "__main__":
    main()
