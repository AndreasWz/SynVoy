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
        load_genome, reverse_complement, translate, parse_gff as parse_gff_base,
        sw_pair_evalue, set_genetic_code,
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, reverse_complement, translate, parse_gff as parse_gff_base,
        sw_pair_evalue, set_genetic_code,
    )

try:
    from gene_predictor import predict_genes
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from gene_predictor import predict_genes

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
    feature_ids = set()   # IDs of every feature that is not a CDS (to find unattached CDS)
    
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
            if feature_type != 'CDS' and attributes.get('ID'):
                feature_ids.add(attributes['ID'])
            
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
    # can still select neighboring models around the GOI region. In a file that also
    # has gene rows (the borrowed annotations of a home genome without its own GFF,
    # merged with the predictions), only CDS that belong to no feature of the file count
    # here: no Parent, or a Parent no row defines (batch HG-5a: one borrowed gene row
    # used to switch the fallback off, and every prediction was lost).
    mixed = bool(processed_genes)
    if mixed:
        cds_by_parent = {pid: rows for pid, rows in cds_by_parent.items() if pid not in feature_ids}
    if cds_by_parent or cds_orphans:
        if not mixed:
            print("GFF contains no gene features; falling back to CDS-derived pseudo genes")
        else:
            print(f"GFF mixes gene features with {len(cds_by_parent) + len(cds_orphans)} CDS "
                  "model(s) that belong to no gene; adding them as CDS-derived pseudo genes")

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

        processed_genes = processed_genes + pseudo_genes if mixed else pseudo_genes

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


LCS_FALLBACK_THRESHOLD = 35.0   # only used when no aligner is available, see is_goi_relative


# E-value at or below which the annotated gene at the GOI's own position counts as the
# GOI itself, whatever --goi_relative_evalue says about the neighbours.
GOI_SELF_EVALUE = 1e-5
BEST_INTERVAL = 'best_interval'   # split_loci.py's name for the interval a locus walk starts from


def is_goi_relative(protein: str, goi_seqs, max_evalue: float, statistic: str = 'sw') -> bool:
    """True when ``protein`` aligns to one of the GOI sequences with E <= ``max_evalue``.

    ``statistic``: 'sw' (Smith-Waterman, BLOSUM62 11/1, Karlin-Altschul) or 'blastp_cbs'
    (the label's family test, synvoy_specificity.pairwise_evalue: blastp with
    composition-based statistics), with which a low-complexity neighbour of a
    low-complexity GOI is no longer a "relative" by its composition alone (preset plan
    PR-8, the low_complexity module; melittin's PEPEP pro-region passes the plain test at
    E 8.4e-7).

    A gene related to the GOI is not used as an anchor: its matches in a target genome
    fall on the GOI's own relatives, so it would count the same sequence evidence a
    second time as "neighbourhood".

    The test is a local alignment judged by its E-value. Until 2026-10-08 it was the
    longest common subsequence divided by the longer length, with a threshold of 35.
    Two unrelated proteins of similar length score 35 on that measure by chance, so about
    half of all genes of the GOI's size were removed, whatever their sequence (measured on
    five benchmark loci: 22 of 33 removed genes had no significant alignment to the GOI),
    and the nearest neighbours of the gene were lost as anchors.

    ``max_evalue <= 0`` switches the test off. Without parasail the old measure is used,
    with a warning, so that the step still runs in a degraded environment.
    """
    if not protein or not goi_seqs or max_evalue <= 0:
        return False
    if statistic == 'blastp_cbs':
        # The GOI is the query, as in the label's family test and home_context.py (review fix
        # RF-2: the two directions differ; this call had the neighbour as the query).
        from synvoy_specificity import pairwise_evalue
        return any(pairwise_evalue(gs, protein) <= max_evalue for gs in goi_seqs)
    try:
        return any(sw_pair_evalue(protein, gs)[0] <= max_evalue for gs in goi_seqs)
    except RuntimeError as exc:
        if not getattr(is_goi_relative, '_warned', False):
            print(f"WARNING: no Smith-Waterman aligner ({exc}). The GOI-relative test "
                  f"falls back to the common-subsequence measure, which sits at chance level "
                  f"for proteins of similar length.", file=sys.stderr)
            is_goi_relative._warned = True
        return any(_lcs_similarity(protein, gs) >= LCS_FALLBACK_THRESHOLD for gs in goi_seqs)


def _lcs_similarity(seq_a: str, seq_b: str) -> float:
    """
    Longest Common Subsequence (LCS) similarity (0-100), normalized by the
    longer sequence length. Kept only as the fallback of is_goi_relative: for two
    unrelated proteins of similar length it is about 35 by chance.
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
    parser.add_argument("--exon_mode", type=str, default="false", 
                        help="If true, output individual exon CDS sequences instead of full protein")
    parser.add_argument("--pred_flank_window", type=int, default=50000,
                        help="Flanking window around GOI hits for gene prediction")
    parser.add_argument("--pred_keep_pct", type=float, default=0.10,
                        help="Fraction of longest gene predictions to keep")
    parser.add_argument("--gene_predictor", type=str, default="auto",
                        help="Gene predictor: auto, augustus, or prodigal (default: auto)")
    parser.add_argument("--augustus_species", type=str, default="fly",
                        help="Augustus species model (default: fly)")
    parser.add_argument("--goi_faa", default=None,
                        help="GOI protein FASTA — flanking candidates too similar to the GOI "
                             "will be filtered out (avoids picking LY6/3FTx family members as anchors)")
    parser.add_argument("--goi_relative_evalue", type=float, default=1e-5,
                        help="A gene whose protein aligns to the GOI with an E-value at or below "
                             "this (Smith-Waterman, BLOSUM62, pairwise) is a relative of the GOI "
                             "and is not used as an anchor. 0 switches the test off "
                             "(default: 1e-5)")
    parser.add_argument("--goi_relative_statistic", choices=["sw", "blastp_cbs"], default="sw",
                        help="Statistic of that test: Smith-Waterman E-value (default) or blastp "
                             "with composition-based statistics (the label's family test). The "
                             "gene at the GOI's own position is always tested with sw")
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
    
    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the home genome (default 1, the standard code)")
    args = parser.parse_args()
    set_genetic_code(args.genetic_code)
    exon_mode = args.exon_mode.lower() == 'true'
    expand_goi_similar = args.expand_goi_similar.lower() == 'true'

    if exon_mode:
        print("Exon mode enabled: extracting individual CDS exon sequences")

    # Load GOI sequences for similarity filtering
    goi_seqs = []
    if args.goi_faa and args.goi_faa != 'NO_GOI' and os.path.exists(args.goi_faa):
        goi_seqs = _load_goi_sequences(args.goi_faa)
        statistic = ("blastp, composition-based statistics" if args.goi_relative_statistic == "blastp_cbs"
                     else "Smith-Waterman")
        print(f"Loaded {len(goi_seqs)} GOI sequence(s); genes related to them "
              f"(alignment E <= {args.goi_relative_evalue:g}, {statistic}) are not used as anchors")

    # Load INPUT Genes: the hit intervals of one home locus (split_loci.py). A row that
    # split_loci snapped to an annotated gene carries its hit span in columns 8-9.
    target_regions = []
    with open(args.bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                region = {
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': parts[3] if len(parts) > 3 else '',
                }
                if len(parts) >= 9 and parts[7].isdigit() and parts[8].isdigit():
                    region['hit_start'], region['hit_end'] = int(parts[7]), int(parts[8])
                target_regions.append(region)

    if not target_regions:
        print("No target regions found.")
        open(args.out_bed, 'w').close()
        open(args.out_faa, 'w').close()
        return
    # One walk per home locus (batch HG-1B), from its best interval: the row split_loci
    # names best_interval, else the first row (a BED with one row, or one written before).
    # The locus's other intervals are copies or relatives of the gene: genes there are not
    # anchors, and each is written as a GOI_ row without a protein. Until 2026-10-16 every
    # row had a walk of its own and a GOI pseudo-gene.
    walk_region = next((r for r in target_regions if r.get('name') == BEST_INTERVAL),
                       target_regions[0])
    other_regions = [r for r in target_regions if r is not walk_region]

    def _at_other_interval(g):
        return any(g['chrom'] == r['chrom'] and g['start'] < r['end'] and g['end'] > r['start']
                   for r in other_regions)

    genome_seqs = load_genome(args.genome)

    extracted_genes = []
    goi_similar_neighbors = []  # populated in the GFF path when expand_goi_similar is on
    
    # MODE SWITCH
    if args.gff == "NO_GFF" or not os.path.exists(args.gff):
        predictor_name = getattr(args, 'gene_predictor', 'auto')
        augustus_sp = getattr(args, 'augustus_species', 'fly')
        print(f"No GFF provided. Running gene prediction on flanking regions "
              f"(predictor={predictor_name}, augustus_species={augustus_sp})...")

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

            # Run gene prediction (Augustus for eukaryotes, Prodigal for prokaryotes)
            temp_out_faa = f"temp_{sub_id}.faa"
            temp_out_gff = f"temp_{sub_id}.gff"

            try:
                predict_genes(
                    temp_fa, temp_out_faa, temp_out_gff,
                    predictor=predictor_name,
                    augustus_species=augustus_sp,
                )

                if not os.path.exists(temp_out_faa):
                    continue

                # Parse output FAA using our utility
                # Both Augustus and Prodigal output Prodigal-style headers:
                #   >id # start # end # strand # attrs
                for header, clean_id, seq in parse_fasta(temp_out_faa):
                    parts = header.split(" # ")
                    if len(parts) >= 4:
                        local_start = int(parts[1]) # 1-based
                        local_end = int(parts[2])   # 1-based inclusive
                        strand_code = parts[3].strip() # 1 or -1
                        strand = "+" if strand_code in ("1", "+") else "-"

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
                if os.path.exists(temp_out_gff): os.remove(temp_out_gff)

        # FILTER: Keep only top N% longest gene predictions.
        # Prodigal on eukaryotic genomes greatly overpredicts tiny ORFs;
        # Augustus is better but may still predict partial genes.
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
                        'seq': goi_prot,
                        'placeholder': True,
                    })
                    print(f"Injected GOI pseudo-gene at {region['chrom']}:{goi_start}-{goi_end}")
            # Re-sort by genomic position for proper flanking order
            extracted_genes.sort(key=lambda g: (g['chrom'], g['start']))
            print(f"Prodigal filter: kept {len(extracted_genes)} genes from {original_count} predictions (top {keep_count} longest + GOI overlaps)")

    else:
        # EXISTING LOGIC FOR GFF
        all_genes = parse_gff_for_genes(args.gff)

        for region in [walk_region]:
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
                _relative_cache = {}

                def _is_goi_similar(g: dict) -> bool:
                    key = (g.get('chrom'), g.get('start'), g.get('end'), g.get('id'))
                    if key not in _relative_cache:
                        _relative_cache[key] = is_goi_relative(
                            _gene_protein(g, genome_seqs), goi_seqs, args.goi_relative_evalue,
                            args.goi_relative_statistic)
                    return _relative_cache[key]

                def _is_goi_itself(g: dict) -> bool:
                    # At the GOI's own position the test cannot be switched off or made
                    # stricter than its default: the gene there that aligns to the query IS
                    # the gene of interest, and a gene is never its own neighbour. (With
                    # --goi_relative_evalue 0 it came back as an anchor and raised the
                    # flanking count of every true call by one.)
                    return is_goi_relative(_gene_protein(g, genome_seqs), goi_seqs,
                                           max(args.goi_relative_evalue, GOI_SELF_EVALUE))

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

                # The walk starts one gene out because the gene nearest the locus is
                # normally the gene of interest itself. When the locus overlaps no gene
                # (an unannotated seed such as an ORF, or a locus placed by alignment
                # between genes) that nearest gene is a real neighbour and was lost; it now
                # opens the walk on its own side (bug hunt 2026-10-13, A7).
                centre = chrom_genes[center_idx]
                centre_free = not (centre['start'] < region['end'] and centre['end'] > region['start'])
                centre_up = centre_free and (centre['start'] + centre['end']) / 2 <= reg_center
                centre_down = centre_free and not centre_up

                upstream_all = []   # all non-GOI-similar genes visited (closest→farthest)
                upstream_slots = 0  # count of protein-coding genes (toward n_flank limit)
                max_expand = min(len(chrom_genes), args.n_flank * 5)
                for offset in range(0 if centre_up else 1, max_expand + 1):
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
                    if _at_other_interval(g):
                        print(f"  [flank-filter] Skipping upstream gene at another hit interval "
                              f"of this locus: {_preferred_gene_label(g)}")
                        continue
                    upstream_all.append(g)
                    if g.get('cds_parts'):
                        upstream_slots += 1

                downstream_all = []
                downstream_slots = 0
                for offset in range(0 if centre_down else 1, max_expand + 1):
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
                    if _at_other_interval(g):
                        print(f"  [flank-filter] Skipping downstream gene at another hit interval "
                              f"of this locus: {_preferred_gene_label(g)}")
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
                window_genes = [g for g in chrom_genes[start_idx:end_idx] if not _at_other_interval(g)]

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
            overlap_added_count = 0
            for g in overlap_genes:
                if g not in extracted_genes:
                    if goi_seqs and _is_goi_itself(g):
                        name = _preferred_gene_label(g)
                        print(f"  [overlap-filter] NOT re-injecting GOI-similar overlap gene: {name}")
                        continue
                    extracted_genes.append(g)
                    overlap_added_count += 1
                else:
                    # Already in extracted_genes (e.g. it was a flanking neighbor).
                    overlap_added_count += 1

            # Inject a GOI pseudo-gene either when no overlap genes exist OR
            # when every overlap gene was filtered out as GOI-similar (this is
            # the home-locus case: the GOI itself looks like the GOI). Without
            # this, the home BED ends up with only flanking anchors and the
            # synteny plot's home track has no GOI marker.
            if overlap_added_count == 0 and region['chrom'] in genome_seqs:
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
                    'seq': goi_prot,
                    'placeholder': True,
                })
                print(f"Injected GOI pseudo-gene at {region['chrom']}:{goi_start}-{goi_end}")

            # The locus's other intervals: the annotated gene at each (the largest overlap) is
            # written as a GOI_ row without a protein, so the region scoring counts it as the
            # GOI and not as an anchor (cluster_grs.py) and the figures draw it as a copy. An
            # interval without a gene gets no row.
            for other in other_regions:
                at = [(min(g['end'], other['end']) - max(g['start'], other['start']), i, g)
                      for i, g in enumerate(all_genes)
                      if g['chrom'] == other['chrom'] and g['start'] < other['end']
                      and g['end'] > other['start']]
                if not at:
                    continue
                gene = max(at, key=lambda t: (t[0], -t[1]))[2]
                goi_row = dict(gene)
                goi_row['attrs'] = dict(gene.get('attrs') or {}, ID=f"GOI_{gene['chrom']}_{gene['start']}")
                goi_row['placeholder'] = True
                extracted_genes.append(goi_row)
                print(f"GOI row for the gene at another hit interval of this locus: "
                      f"{_preferred_gene_label(gene)} {gene['chrom']}:{gene['start']}-{gene['end']}")

    # Write Outputs
    all_fasta_records = []  # Collect all FASTA records
    # Home coordinates of every anchor and the length of the home genome travel with the
    # proteins (FASTA header fields `home=` and `genome=`): the search derives its
    # distances from them (bin/synvoy_scale.py).
    genome_bp = sum(len(seq) for seq in genome_seqs.values())
    # The home locus itself (`goi=`), so the search knows on which side of the gene each
    # anchor lies (its slot between the nearest neighbours, bin/synvoy_evalue.py).
    # Over every interval of the locus, from the hits (a row snapped to a gene keeps its hit
    # span in columns 8-9), so snapping leaves the search's slots as they were.
    goi_span = {}
    for region in target_regions:
        h0, h1 = region.get('hit_start', region['start']), region.get('hit_end', region['end'])
        lo, hi = goi_span.get(region['chrom'], (h0, h1))
        goi_span[region['chrom']] = (min(lo, h0), max(hi, h1))

    def _home_fields(gene):
        fields = f"home={gene['start']}-{gene['end']} genome={genome_bp}"
        if gene.get('chrom') in goi_span:
            lo, hi = goi_span[gene['chrom']]
            fields += f" goi={lo}-{hi}"
        return fields

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
            if gene.get('placeholder'):
                # The marker of the GOI in the home BED, not a query: its "protein" is the
                # translation of the hit interval on the plus strand up to the first stop.
                # For a single-exon gene on the plus strand that is the gene once more; it
                # then entered the search as a second GOI and collected most of the
                # low-confidence rows (yeast STE2: 16 of 19). For any other gene it is noise.
                pass
            elif 'seq' in gene:
                # From prediction - no exon info available, use whole sequence
                prot_seq = str(gene['seq'])
                if len(prot_seq) * 3 >= args.min_size:
                    all_fasta_records.append((f"{seq_header} {_home_fields(gene)}", prot_seq))
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
                            full_header = (f"{gid} full_length_protein exons={len(cds_parts)} "
                                           f"strand={gene['strand']} {_home_fields(gene)}")
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
                            all_fasta_records.append((f"{seq_header} {_home_fields(gene)}", prot_seq))
                        
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
