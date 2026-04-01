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
            # Prefer a transcript that has CDS features (skips misc_RNA / lncRNA
            # isoforms listed first, which would leave cds_parts empty for
            # protein-coding genes like THEM6 that have both misc_RNA and mRNA
            # transcripts in NCBI RefSeq GFFs).
            best_t = next(
                (t for t in transcripts if t in cds_by_parent),
                transcripts[0]
            )
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


def _lcs_similarity(seq_a: str, seq_b: str) -> float:
    """
    Longest Common Subsequence (LCS) similarity (0-100), normalized by the
    longer sequence length.  More sensitive than k-mer Jaccard for detecting
    distant protein family relationships (e.g. LY6/3FTx members vs LY6E).
    O(m*n) — fast enough for typical flanking gene protein lengths (<1 kb aa).
    """
    m, n = len(seq_a), len(seq_b)
    if not m or not n:
        return 0.0
    prev = [0] * (n + 1)
    for c in seq_a:
        curr = [0] * (n + 1)
        for j, d in enumerate(seq_b, 1):
            curr[j] = prev[j - 1] + 1 if c == d else max(prev[j], curr[j - 1])
        prev = curr
    return 100.0 * prev[n] / max(m, n)


def _gene_protein(gene: dict, genome_seqs: dict) -> str:
    """
    Translate the CDS of a gene dict to a protein string.
    Returns empty string if CDS parts are unavailable or too short.
    """
    cds_parts = gene.get('cds_parts', [])
    if not cds_parts:
        return ''
    chrom = gene.get('chrom', '')
    if chrom not in genome_seqs:
        return ''
    strand = gene.get('strand', '+')
    parts_sorted = sorted(cds_parts, key=lambda p: p['start'])
    cds_dna = ''
    for part in parts_sorted:
        seg = genome_seqs[chrom][part['start']:part['end']]
        cds_dna += seg
    if not cds_dna:
        return ''
    if strand == '-':
        cds_dna = reverse_complement(cds_dna)
    cds_dna = cds_dna[:len(cds_dna) - (len(cds_dna) % 3)]
    if len(cds_dna) < 30:
        return ''
    protein = translate(cds_dna)
    # Take longest ORF (drop stop codons)
    return protein.split('*')[0]


def _load_goi_sequences(goi_faa: str) -> list:
    """Parse a FASTA file and return list of protein sequences."""
    seqs = []
    current = []
    with open(goi_faa) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current))
                    current = []
            else:
                current.append(line)
    if current:
        seqs.append(''.join(current))
    # Keep only full-length sequences; exon fragments (<50 aa) cause false-positive
    # similarity matches against unrelated proteins.
    full_seqs = [s for s in seqs if len(s) >= 50]
    return full_seqs if full_seqs else [s for s in seqs if s]


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
    parser.add_argument("--goi_faa", default=None,
                        help="GOI protein FASTA — flanking candidates too similar to the GOI "
                             "will be filtered out (avoids picking LY6/3FTx family members as anchors)")
    parser.add_argument("--max_goi_similarity", type=float, default=35.0,
                        help="Max allowed LCS similarity (0-100) to GOI; genes above this threshold "
                             "are excluded from the flanking set (default: 35.0)")
    parser.add_argument("--max_flanking_distance", type=int, default=0,
                        help="Max distance (bp) from GOI center to walk for flanking genes. "
                             "0 = unlimited (legacy). Prevents reaching into gene deserts when "
                             "the GOI neighbours a tandem array that fills one side (default: 0)")
    parser.add_argument("--expand_goi_similar", type=str, default="false",
                        help="If true, GOI-similar genes within --expand_goi_similar_distance "
                             "are emitted as additional GOI queries (GOI_NEIGHBOR_ prefix) so "
                             "all paralogs/copies are searched and end up in the IQ-TREE "
                             "(default: false)")
    parser.add_argument("--expand_goi_similar_distance", type=int, default=300000,
                        help="Max bp from GOI center for GOI-similar genes to be treated as "
                             "additional GOI queries when --expand_goi_similar is set (default: 300000)")
    parser.add_argument("--out_bed", required=True, help="Output BED")
    parser.add_argument("--out_faa", required=True, help="Output FASTA")
    
    args = parser.parse_args()
    prefer_large = args.prefer_large.lower() == 'true'
    exon_mode = args.exon_mode.lower() == 'true'
    expand_goi_similar = args.expand_goi_similar.lower() == 'true'

    if exon_mode:
        print("Exon mode enabled: extracting individual CDS exon sequences")

    # Load GOI sequences for similarity filtering
    goi_seqs = []
    if args.goi_faa and args.goi_faa != 'NO_GOI' and os.path.exists(args.goi_faa):
        goi_seqs = _load_goi_sequences(args.goi_faa)
        print(f"Loaded {len(goi_seqs)} GOI sequence(s) for flanking similarity filter "
              f"(max_similarity={args.max_goi_similarity}%)")

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
    goi_similar_neighbors = []  # populated in the GFF path when expand_goi_similar is on
    
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
            # If GOI-similarity filtering is active, expand the window until we
            # collect n_flank non-similar genes on each side.  This prevents
            # the flanking set being filled entirely with GOI-family members
            # (e.g. LY6/3FTx proteins) that are useless as synteny anchors.
            if goi_seqs:
                def _is_goi_similar(g: dict) -> bool:
                    seq = _gene_protein(g, genome_seqs)
                    if not seq:
                        return False
                    return any(
                        _lcs_similarity(seq, gs) >= args.max_goi_similarity
                        for gs in goi_seqs
                    )

                # Walk outward from the GOI center collecting genes into the window.
                # - GOI-similar genes are skipped as flanking anchors (useless for synteny).
                #   When --expand_goi_similar is set AND the gene is within
                #   --expand_goi_similar_distance, it is instead queued as an additional
                #   GOI query (GOI_NEIGHBOR_ prefix in FAA) so all paralogs/copies are
                #   searched in target genomes and end up in IQ-TREE.
                # - Non-coding genes (no CDS) are added to the window for BED
                #   visualization but do NOT count toward the n_flank protein-coding
                #   quota — this prevents lncRNAs and pseudogenes from pushing
                #   useful synteny anchors (e.g. TOP1MT) out of the flanking set.
                # - When --max_flanking_distance > 0, the walk stops once the gene
                #   is beyond that many bp from the GOI center (prevents reaching
                #   1Mb+ away when the GOI neighbours a tandem array that eats one side).
                goi_similar_neighbors = []  # genes to emit as GOI_NEIGHBOR_ queries

                upstream_all = []   # all non-GOI-similar genes visited (closest→farthest)
                upstream_slots = 0  # count of protein-coding genes (toward n_flank limit)
                max_expand = min(len(chrom_genes), args.n_flank * 5)
                for offset in range(1, max_expand + 1):
                    if upstream_slots >= args.n_flank:
                        break
                    idx = center_idx - offset
                    if idx < 0:
                        break
                    g = chrom_genes[idx]
                    # Distance cap: chrom_genes is sorted by position; upstream genes
                    # move monotonically farther from GOI as offset increases.
                    if args.max_flanking_distance > 0:
                        gene_center = (g['start'] + g['end']) / 2
                        if abs(gene_center - reg_center) > args.max_flanking_distance:
                            print(f"  [flank-dist] Upstream walk stopped at offset {offset}: "
                                  f"gene {_preferred_gene_label(g)} is "
                                  f"{abs(gene_center - reg_center)/1000:.0f}kb away "
                                  f"(limit {args.max_flanking_distance/1000:.0f}kb)")
                            break
                    if _is_goi_similar(g):
                        name = _preferred_gene_label(g)
                        gene_center = (g['start'] + g['end']) / 2
                        if expand_goi_similar and abs(gene_center - reg_center) <= args.expand_goi_similar_distance:
                            goi_similar_neighbors.append(g)
                            print(f"  [goi-expand] Upstream GOI-similar gene queued as additional GOI: {name}")
                        else:
                            print(f"  [flank-filter] Skipping GOI-similar upstream gene: {name}")
                        continue
                    upstream_all.append(g)
                    if g.get('cds_parts'):
                        upstream_slots += 1

                downstream_all = []
                downstream_slots = 0
                for offset in range(1, max_expand + 1):
                    if downstream_slots >= args.n_flank:
                        break
                    idx = center_idx + offset
                    if idx >= len(chrom_genes):
                        break
                    g = chrom_genes[idx]
                    if args.max_flanking_distance > 0:
                        gene_center = (g['start'] + g['end']) / 2
                        if abs(gene_center - reg_center) > args.max_flanking_distance:
                            print(f"  [flank-dist] Downstream walk stopped at offset {offset}: "
                                  f"gene {_preferred_gene_label(g)} is "
                                  f"{abs(gene_center - reg_center)/1000:.0f}kb away "
                                  f"(limit {args.max_flanking_distance/1000:.0f}kb)")
                            break
                    if _is_goi_similar(g):
                        name = _preferred_gene_label(g)
                        gene_center = (g['start'] + g['end']) / 2
                        if expand_goi_similar and abs(gene_center - reg_center) <= args.expand_goi_similar_distance:
                            goi_similar_neighbors.append(g)
                            print(f"  [goi-expand] Downstream GOI-similar gene queued as additional GOI: {name}")
                        else:
                            print(f"  [flank-filter] Skipping GOI-similar downstream gene: {name}")
                        continue
                    downstream_all.append(g)
                    if g.get('cds_parts'):
                        downstream_slots += 1

                window_genes = list(reversed(upstream_all)) + downstream_all
                print(
                    f"  [flank-filter] Selected {upstream_slots} upstream + "
                    f"{downstream_slots} downstream protein-coding flanking genes "
                    f"(+{len(upstream_all)-upstream_slots+len(downstream_all)-downstream_slots} "
                    f"non-coding neighbours; expanded up to {max_expand} candidates per side)."
                )
                if goi_similar_neighbors:
                    print(f"  [goi-expand] {len(goi_similar_neighbors)} GOI-similar neighbor(s) "
                          f"will be emitted as additional GOI queries in FAA.")
            else:
                start_idx = max(0, center_idx - args.n_flank)
                end_idx = min(len(chrom_genes), center_idx + args.n_flank + 1)
                window_genes = chrom_genes[start_idx:end_idx]

            # Add genes from window
            extracted_genes.extend(window_genes)

            # Ensure genes overlapping the target region are included,
            # but respect the GOI-similarity filter: do NOT re-inject genes
            # that would have been excluded by _is_goi_similar(). Without
            # this gate the overlap injection undoes the walking-loop filter.
            overlap_genes = [
                g for g in chrom_genes
                if g['start'] < region['end'] and g['end'] > region['start']
            ]
            for g in overlap_genes:
                if g not in extracted_genes:
                    if goi_seqs and _is_goi_similar(g):
                        name = _preferred_gene_label(g)
                        print(f"  [overlap-filter] NOT re-injecting GOI-similar overlap gene: {name}")
                        continue
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
                        # GFF phase = bases to skip at the coding 5' end.
                        # + strand: first exon (index 0) carries the phase.
                        # - strand: last exon genomically carries the phase;
                        #           apply after reverse-complementing.
                        coding_first_phase = int(cds_parts[0].get('phase', 0) or 0) \
                            if gene['strand'] == '+' \
                            else int(cds_parts[-1].get('phase', 0) or 0)
                        dna_seq = ""
                        for part in cds_parts:
                            part_seq = seq_record[part['start']:part['end']]
                            dna_seq += part_seq
                        
                        if gene['strand'] == '-':
                            dna_seq = reverse_complement(dna_seq)

                        # Trim phase at coding 5' end
                        if coding_first_phase > 0:
                            dna_seq = dna_seq[coding_first_phase:]
                            
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
                    # No CDS parts — non-coding gene (lncRNA, pseudogene, etc.).
                    # Include in BED (already written above) but do NOT emit a
                    # protein sequence: translating genomic DNA in frame 0 would
                    # produce junk that pollutes the MMseqs2 search DB.
                    pass
    
    # Emit GOI-similar neighbors as additional GOI queries (GOI_NEIGHBOR_ prefix).
    # These are NOT in the BED (don't affect synteny denominator), but ARE in the
    # FAA so they are searched in target genomes and end up in IQ-TREE alongside the
    # main GOI, allowing paralogs vs orthologs to be resolved by the tree.
    if expand_goi_similar and goi_similar_neighbors:
        genome_seqs_loaded = genome_seqs  # already loaded above
        emitted = 0
        for g in goi_similar_neighbors:
            raw_id = g['attrs'].get('ID', '') or f"{g['chrom']}_{g['start']}"
            neighbor_id = f"GOI_NEIGHBOR_{raw_id}"

            if g.get('cds_parts'):
                cds_parts = sorted(g['cds_parts'], key=lambda x: x['start'])
                if exon_mode:
                    for exon_idx, part in enumerate(cds_parts, start=1):
                        exon_id = f"{neighbor_id}|exon_{exon_idx}"
                        exon_dna = genome_seqs_loaded[g['chrom']][part['start']:part['end']]
                        if g['strand'] == '-':
                            exon_dna = reverse_complement(exon_dna)
                        phase = int(part.get('phase', 0)) if part.get('phase', '.') != '.' else 0
                        if phase > 0:
                            exon_dna = exon_dna[phase:]
                        remainder = len(exon_dna) % 3
                        if remainder:
                            exon_dna = exon_dna[:-remainder]
                        if len(exon_dna) < 9:
                            continue
                        exon_prot = translate(exon_dna).replace('*', '')
                        exon_header = (f"{exon_id} parent={neighbor_id} "
                                       f"exon={exon_idx}/{len(cds_parts)} "
                                       f"coords={part['start']}-{part['end']} "
                                       f"strand={g['strand']}")
                        all_fasta_records.append((exon_header, exon_prot))
                    # Full-length protein
                    full_dna = ""
                    for part in cds_parts:
                        full_dna += genome_seqs_loaded[g['chrom']][part['start']:part['end']]
                    if g['strand'] == '-':
                        full_dna = reverse_complement(full_dna)
                    remainder = len(full_dna) % 3
                    if remainder:
                        full_dna = full_dna[:-remainder]
                    full_prot = translate(full_dna)
                    if '*' in full_prot:
                        full_prot = full_prot.split('*')[0]
                    if len(full_prot) >= 10:
                        all_fasta_records.append((
                            f"{neighbor_id} full_length_protein exons={len(cds_parts)} strand={g['strand']}",
                            full_prot
                        ))
                else:
                    coding_first_phase = int(cds_parts[0].get('phase', 0) or 0) \
                        if g['strand'] == '+' \
                        else int(cds_parts[-1].get('phase', 0) or 0)
                    dna_seq = ""
                    for part in cds_parts:
                        dna_seq += genome_seqs_loaded[g['chrom']][part['start']:part['end']]
                    if g['strand'] == '-':
                        dna_seq = reverse_complement(dna_seq)
                    if coding_first_phase > 0:
                        dna_seq = dna_seq[coding_first_phase:]
                    remainder = len(dna_seq) % 3
                    if remainder:
                        dna_seq = dna_seq[:-remainder]
                    prot_seq = translate(dna_seq)
                    if '*' in prot_seq:
                        prot_seq = prot_seq.split('*')[0]
                    if prot_seq:
                        all_fasta_records.append((neighbor_id, prot_seq))
                emitted += 1
        if emitted:
            print(f"  [goi-expand] Emitted {emitted} GOI-similar neighbor(s) as GOI_NEIGHBOR_ sequences in FAA.")

    # Write all FASTA records at once
    write_fasta(all_fasta_records, args.out_faa)

if __name__ == "__main__":
    main()
