#!/usr/bin/env python3
"""
annotate_goi_exons.py - Annotate GOI exons in home genome

Core idea:
  We want to find and annotate individual exons of the Gene of Interest (GOI)
  in the home genome. This gives us exon-level protein sequences that can be
  searched individually (protein → DNA) in target genomes.

Two modes:
  A) GFF available: Match GOI to annotated gene, extract CDS/exons directly
  B) No GFF: Use tblastn/MMseqs2 hits to identify exon boundaries via
     splice site detection (GT-AG rule), start/stop codon analysis

Output:
  - goi_exons.faa: Individual exon protein sequences + full GOI protein
  - goi_annotation.bed: Exon locations in home genome (BED format)
  - goi_info.json: Metadata about the annotation

IMPORTANT: Always uses protein → DNA search, never DNA → DNA
"""

import argparse
import json
import os
import re
import subprocess
import sys
import urllib.request
import urllib.error

try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, reverse_complement, translate, parse_gff, get_feature_id,
        parse_gff_attributes
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, reverse_complement, translate, parse_gff, get_feature_id,
        parse_gff_attributes
    )


# =============================================================================
# UNIPROT NAME LOOKUP
# =============================================================================

def fetch_uniprot_names(query_id):
    """
    Query UniProt REST API for gene names and synonyms.
    Returns list of gene name strings (primary + synonyms).
    Returns empty list on failure (no internet, invalid ID, etc.)
    """
    names = []

    # Check if query_id looks like a UniProt accession
    if not re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', query_id):
        return names

    url = f"https://rest.uniprot.org/uniprotkb/{query_id}.json"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())

        # Extract gene names
        for gene in data.get('genes', []):
            if 'geneName' in gene and 'value' in gene['geneName']:
                names.append(gene['geneName']['value'])
            for syn in gene.get('synonyms', []):
                if 'value' in syn:
                    names.append(syn['value'])
            for oln in gene.get('orderedLocusNames', []):
                if 'value' in oln:
                    names.append(oln['value'])
            for orf in gene.get('orfNames', []):
                if 'value' in orf:
                    names.append(orf['value'])

        # Also extract protein name
        protein_desc = data.get('proteinDescription', {})
        rec_name = protein_desc.get('recommendedName', {})
        if 'fullName' in rec_name and 'value' in rec_name['fullName']:
            names.append(rec_name['fullName']['value'])

        print(f"[UniProt] Found names for {query_id}: {names}")
    except Exception as e:
        print(f"[UniProt] Could not fetch names for {query_id}: {e}")

    return names


# =============================================================================
# GFF-BASED GOI MATCHING (Scenario B)
# =============================================================================

def match_goi_in_gff(query_id, query_seq, gff_file, genome_file, blast_hits, mmseqs_hits,
                     gff_search_window):
    """
    Try to match the GOI to an annotated gene in the GFF.

    Strategy (in order):
    1. If query_id is a UniProt ID → fetch gene names from UniProt API
       Search GFF for those names
    2. If query_id is an NCBI gene name → search GFF for that name
    3. If only sequence or name not found:
       a. Use LOCATE_GENE hits to know approximate genomic location
       b. Extract proteins from GFF genes in that region only
       c. Compare those regional proteins to GOI (fast, small search)

    Returns:
        dict with 'gene_id', 'chrom', 'strand', 'cds_parts' (list of
        {'start', 'end', 'phase'}) or None if not found
    """
    print("[GFF Match] Attempting to match GOI to annotated gene...")

    # Collect candidate gene names
    candidate_names = []

    # Try UniProt lookup
    uniprot_names = fetch_uniprot_names(query_id)
    candidate_names.extend(uniprot_names)

    # Also add the raw query_id itself as a candidate name
    candidate_names.append(query_id)

    # Add common transformations
    for name in list(candidate_names):
        candidate_names.append(name.lower())
        candidate_names.append(name.upper())
        # NCBI-style: gene-Name or gene-LOCxxxxxx
        candidate_names.append(f"gene-{name}")

    # Deduplicate while preserving order
    seen = set()
    unique_names = []
    for n in candidate_names:
        if n.lower() not in seen:
            seen.add(n.lower())
            unique_names.append(n)
    candidate_names = unique_names

    print(f"[GFF Match] Searching GFF for names: {candidate_names[:10]}...")

    # --- Strategy 1: Name-based search in GFF ---
    match = _search_gff_by_name(gff_file, candidate_names)
    if match:
        print(f"[GFF Match] Found by name: {match['gene_id']}")
        return match

    # --- Strategy 2: Region-based protein search ---
    # Use LOCATE_GENE hits to narrow the search region
    hit_regions = _parse_hit_regions(blast_hits, mmseqs_hits)
    if hit_regions:
        print(f"[GFF Match] Name not found. Searching proteins in {len(hit_regions)} hit region(s)...")
        match = _search_gff_by_region_and_sequence(
            gff_file, genome_file, query_seq, hit_regions, gff_search_window
        )
        if match:
            print(f"[GFF Match] Found by regional sequence search: {match['gene_id']}")
            return match

    print("[GFF Match] GOI not found in GFF annotations.")
    return None


def _search_gff_by_name(gff_file, candidate_names):
    """
    Search GFF for genes matching any of the candidate names.
    
    Scoring system:
    - Exact Name match: 100 points
    - Exact ID match: 90 points
    - Name/description contains full candidate name (>4 chars): 50 points
    - Partial word match (>5 chars): 20 points
    
    Collects ALL matches, returns the highest-scoring one.
    """
    # Build a set of lowercased names for fast lookup
    name_set = set(n.lower() for n in candidate_names)

    # Build word patterns for partial matching (only longer words, skip generic ones)
    SKIP_WORDS = {'gene', 'protein', 'like', 'family', 'domain', 'type', 'group'}
    name_words = set()
    for n in candidate_names:
        for word in re.split(r'[^a-zA-Z0-9]+', n):
            if len(word) > 4 and word.lower() not in SKIP_WORDS:
                name_words.add(word.lower())

    matched_genes = {}  # gene_id -> {'score': N, 'info': {...}}

    gene_features = {}  # gene_id -> {'chrom', 'start', 'end', 'strand'}
    transcript_to_gene = {}  # transcript_id -> gene_id
    cds_by_parent = {}  # parent_id -> [cds_parts]

    for feature in parse_gff(gff_file):
        attrs = feature.get('attributes', {})

        if feature['type'] == 'gene':
            gene_id = attrs.get('ID', '')
            gene_name = attrs.get('Name', '')
            gene_desc = attrs.get('description', attrs.get('product', ''))
            dbxref = attrs.get('Dbxref', '')
            gene_sym = attrs.get('gene', '')  # NCBI uses 'gene' attribute

            # Score this gene against our candidates
            match_score = 0

            # Check exact Name matches (highest priority)
            name_lower = gene_name.lower() if gene_name else ''
            gene_sym_lower = gene_sym.lower() if gene_sym else ''
            desc_lower = gene_desc.lower() if gene_desc else ''

            for candidate in name_set:
                # Exact gene name/symbol match
                if candidate and (candidate == name_lower or candidate == gene_sym_lower):
                    match_score = max(match_score, 100)
                # Exact ID match (e.g., gene-MELT)
                elif candidate and candidate == gene_id.lower():
                    match_score = max(match_score, 90)
                # Description contains the full candidate name (if >4 chars)
                elif candidate and len(candidate) > 4 and candidate in desc_lower:
                    match_score = max(match_score, 50)
                # Name contains candidate (e.g. gene name "Melt" matches "melt")
                elif candidate and len(candidate) > 3 and candidate in name_lower:
                    match_score = max(match_score, 60)

            # Partial word matching (lower priority, only for longer words)
            if match_score == 0:
                searchable_lower = f"{name_lower} {desc_lower} {gene_sym_lower}"
                for word in name_words:
                    if word in searchable_lower:
                        match_score = max(match_score, 20)
                        break

            if match_score > 0:
                gene_features[gene_id] = {
                    'chrom': feature['seqid'],
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand'],
                    'gene_name': gene_name or gene_sym or gene_id
                }
                matched_genes[gene_id] = match_score
                print(f"[GFF Match]   Candidate: {gene_id} (Name={gene_name}, "
                      f"desc={gene_desc[:40]}...) score={match_score}")

        elif feature['type'] in ['mRNA', 'transcript']:
            tid = attrs.get('ID', '')
            parent = attrs.get('Parent', '')
            if tid and parent:
                transcript_to_gene[tid] = parent

        elif feature['type'] == 'CDS':
            parent = attrs.get('Parent', '')
            if parent:
                if parent not in cds_by_parent:
                    cds_by_parent[parent] = []
                cds_by_parent[parent].append({
                    'chrom': feature['seqid'],
                    'start': feature['start'] - 1,  # GFF is 1-based, convert to 0-based
                    'end': feature['end'],
                    'strand': feature['strand'],
                    'phase': feature.get('phase', 0) or 0
                })

    if not matched_genes:
        return None

    # Pick the highest-scoring gene
    best_gene = max(matched_genes, key=matched_genes.get)
    best_score = matched_genes[best_gene]
    print(f"[GFF Match] Best match: {best_gene} (score={best_score})")

    if best_gene not in gene_features:
        return None

    # Find CDS for the best matching gene
    gene_info = gene_features[best_gene]
    cds_parts = []

    # Check direct CDS children of the gene
    if best_gene in cds_by_parent:
        cds_parts = cds_by_parent[best_gene]

    # Check via transcript
    if not cds_parts:
        for tid, gid in transcript_to_gene.items():
            if gid == best_gene and tid in cds_by_parent:
                cds_parts = cds_by_parent[tid]
                break  # Use first transcript (TODO: pick longest)

    if not cds_parts:
        print(f"[GFF Match] Gene {best_gene} found but no CDS features. Using gene span.")
        # Fall back to gene span
        cds_parts = [{
            'chrom': gene_info['chrom'],
            'start': gene_info['start'] - 1,
            'end': gene_info['end'],
            'strand': gene_info['strand'],
            'phase': 0
        }]

    return {
        'gene_id': best_gene,
        'gene_name': gene_info.get('gene_name', best_gene),
        'chrom': gene_info['chrom'],
        'start': gene_info['start'],
        'end': gene_info['end'],
        'strand': gene_info['strand'],
        'cds_parts': sorted(cds_parts, key=lambda x: x['start'])
    }


def _parse_hit_regions(blast_hits_file, mmseqs_hits_file):
    """
    Parse BLAST/MMseqs raw hit files (format 6) to get genomic regions.
    Returns list of {'chrom', 'start', 'end', 'qstart', 'qend', 'strand'}.
    """
    regions = []

    for hits_file in [blast_hits_file, mmseqs_hits_file]:
        if not hits_file or not os.path.exists(hits_file) or os.path.getsize(hits_file) == 0:
            continue
        try:
            with open(hits_file) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue
                    try:
                        qstart = int(parts[6])
                        qend = int(parts[7])
                        tstart = int(parts[8])
                        tend = int(parts[9])
                        evalue = float(parts[10])

                        strand = '+' if tstart < tend else '-'
                        # Convert 1-based BLAST/mmseqs coords to 0-based half-open
                        gstart = min(tstart, tend) - 1  # 0-based start
                        gend = max(tstart, tend)         # exclusive end

                        regions.append({
                            'chrom': parts[1],
                            'start': gstart,
                            'end': gend,
                            'qstart': min(qstart, qend),
                            'qend': max(qstart, qend),
                            'strand': strand,
                            'evalue': evalue
                        })
                    except (ValueError, IndexError):
                        continue
        except Exception:
            continue

    return regions


def _search_gff_by_region_and_sequence(gff_file, genome_file, query_seq, hit_regions,
                                       gff_search_window):
    """
    Search GFF genes near hit regions and compare protein sequences.

    Instead of searching the entire proteome (slow), we:
    1. Identify GFF genes overlapping with or near our hit regions
    2. Extract their protein sequences
    3. Compare to GOI protein
    """
    SEARCH_WINDOW = max(0, int(gff_search_window))

    # Build search windows
    search_windows = []
    for region in hit_regions:
        search_windows.append({
            'chrom': region['chrom'],
            'start': max(0, region['start'] - SEARCH_WINDOW),
            'end': region['end'] + SEARCH_WINDOW
        })

    # Merge overlapping windows
    search_windows.sort(key=lambda x: (x['chrom'], x['start']))
    merged_windows = []
    for w in search_windows:
        if merged_windows and merged_windows[-1]['chrom'] == w['chrom'] and w['start'] <= merged_windows[-1]['end']:
            merged_windows[-1]['end'] = max(merged_windows[-1]['end'], w['end'])
        else:
            merged_windows.append(dict(w))

    # Find GFF genes in these windows
    candidate_genes = {}  # gene_id -> gene_info
    transcript_to_gene = {}
    cds_by_parent = {}

    for feature in parse_gff(gff_file):
        attrs = feature.get('attributes', {})

        if feature['type'] == 'gene':
            gene_id = attrs.get('ID', '')
            # Check if gene overlaps any search window
            for window in merged_windows:
                if (feature['seqid'] == window['chrom'] and
                        feature['start'] <= window['end'] and
                        feature['end'] >= window['start']):
                    candidate_genes[gene_id] = {
                        'chrom': feature['seqid'],
                        'start': feature['start'],
                        'end': feature['end'],
                        'strand': feature['strand'],
                        'gene_name': attrs.get('Name', gene_id)
                    }
                    break

        elif feature['type'] in ['mRNA', 'transcript']:
            tid = attrs.get('ID', '')
            parent = attrs.get('Parent', '')
            if tid and parent:
                transcript_to_gene[tid] = parent

        elif feature['type'] == 'CDS':
            parent = attrs.get('Parent', '')
            if parent:
                if parent not in cds_by_parent:
                    cds_by_parent[parent] = []
                cds_by_parent[parent].append({
                    'chrom': feature['seqid'],
                    'start': feature['start'] - 1,
                    'end': feature['end'],
                    'strand': feature['strand'],
                    'phase': feature.get('phase', 0) or 0
                })

    if not candidate_genes:
        return None

    print(f"[GFF Match] Found {len(candidate_genes)} candidate genes in hit regions")

    # Load genome for protein extraction
    genome_seqs = load_genome(genome_file)

    # Extract proteins for candidate genes and compare to GOI
    best_match = None
    best_identity = 0

    for gene_id, gene_info in candidate_genes.items():
        cds_parts = []

        # Direct CDS
        if gene_id in cds_by_parent:
            cds_parts = cds_by_parent[gene_id]

        # Via transcript
        if not cds_parts:
            for tid, gid in transcript_to_gene.items():
                if gid == gene_id and tid in cds_by_parent:
                    cds_parts = cds_by_parent[tid]
                    break

        if not cds_parts:
            continue

        # Extract and translate protein
        chrom = gene_info['chrom']
        if chrom not in genome_seqs:
            continue

        cds_parts_sorted = sorted(cds_parts, key=lambda x: x['start'])
        dna_seq = ""
        for part in cds_parts_sorted:
            dna_seq += genome_seqs[chrom][part['start']:part['end']]

        if gene_info['strand'] == '-':
            dna_seq = reverse_complement(dna_seq)

        # Trim to codon boundary
        dna_seq = dna_seq[:len(dna_seq) - len(dna_seq) % 3]
        if len(dna_seq) < 9:
            continue

        protein = translate(dna_seq).split('*')[0]  # Stop at first stop

        # Compare to GOI (simple identity calculation)
        identity = _quick_protein_identity(query_seq, protein)

        if identity > best_identity and identity > 30:
            best_identity = identity
            best_match = {
                'gene_id': gene_id,
                'gene_name': gene_info.get('gene_name', gene_id),
                'chrom': chrom,
                'start': gene_info['start'],
                'end': gene_info['end'],
                'strand': gene_info['strand'],
                'cds_parts': cds_parts_sorted,
                'identity': identity,
                'protein': protein
            }

    if best_match:
        print(f"[GFF Match] Best protein match: {best_match['gene_id']} "
              f"({best_match['identity']:.1f}% identity)")
    return best_match


def _quick_protein_identity(seq1, seq2):
    """Quick and dirty protein identity using global alignment heuristic."""
    if not seq1 or not seq2:
        return 0.0

    # Use k-mer based similarity (fast, no alignment needed)
    k = 3
    if len(seq1) < k or len(seq2) < k:
        return 0.0

    kmers1 = set(seq1[i:i + k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i + k] for i in range(len(seq2) - k + 1))

    if not kmers1 or not kmers2:
        return 0.0

    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)

    # Jaccard-like similarity scaled to percentage
    similarity = (intersection / union) * 100

    # Also factor in length ratio
    len_ratio = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))
    similarity *= len_ratio

    return similarity


def extract_exons_from_gff_match(match, genome_file):
    """
    Extract individual CDS/exon protein sequences from a GFF match.

    Returns list of dicts:
    [{'id': 'GOI_genename|exon_1', 'seq': 'MKKV...', 'coords': (start, end)}, ...]
    """
    genome_seqs = load_genome(genome_file)
    chrom = match['chrom']

    if chrom not in genome_seqs:
        print(f"[GFF Exons] Chromosome {chrom} not found in genome!")
        return []

    exons = []
    cds_parts = sorted(match['cds_parts'], key=lambda x: x['start'])

    # If on minus strand, exons are numbered from the 3' end of DNA (= 5' of mRNA)
    # For minus strand: last CDS part (highest coord) is exon 1
    if match['strand'] == '-':
        cds_parts = list(reversed(cds_parts))

    gene_name = match.get('gene_name', match['gene_id'])

    # Extract each exon individually
    for i, part in enumerate(cds_parts, start=1):
        exon_dna = genome_seqs[chrom][part['start']:part['end']]

        if match['strand'] == '-':
            exon_dna = reverse_complement(exon_dna)

        # Handle phase/frame
        phase = int(part.get('phase', 0)) if part.get('phase') not in [None, '.'] else 0
        if phase > 0:
            exon_dna = exon_dna[phase:]

        # Trim to codon boundary
        remainder = len(exon_dna) % 3
        if remainder:
            exon_dna = exon_dna[:len(exon_dna) - remainder]

        if len(exon_dna) < 9:  # < 3 amino acids
            continue

        exon_prot = translate(exon_dna).replace('*', '')

        if len(exon_prot) < 3:
            continue

        exons.append({
            'id': f"GOI_{gene_name}|exon_{i}",
            'seq': exon_prot,
            'coords': (part['start'], part['end']),
            'exon_num': i,
            'total_exons': len(cds_parts),
            'strand': match['strand'],
            'chrom': chrom
        })

    # Also build the full concatenated protein
    all_dna = ""
    for part in sorted(match['cds_parts'], key=lambda x: x['start']):
        all_dna += genome_seqs[chrom][part['start']:part['end']]

    if match['strand'] == '-':
        all_dna = reverse_complement(all_dna)

    all_dna = all_dna[:len(all_dna) - len(all_dna) % 3]
    full_protein = translate(all_dna).split('*')[0]

    return exons, full_protein


# =============================================================================
# TANDEM DUPLICATION DETECTION
# =============================================================================

def detect_tandem_duplications(hits, query_seq, chrom_seq, chrom_name,
                                max_intergenic_distance=50000):
    """
    Detect tandem duplications among search hits.
    
    Distinguishes between:
    - Exons: Hits covering DIFFERENT parts of the query (non-overlapping in query space)
    - Tandem copies: Hits each covering a LARGE/SIMILAR part of the query 
      (overlapping in query space, non-overlapping in genome space)
    
    Criteria for tandem duplication:
    1. Multiple hits on same chromosome within max_intergenic_distance
    2. Hits overlap significantly in QUERY space (each covers >40% of query)
       OR hits that individually align well to the query
    3. Translated products from each hit region have high identity to query
    
    Returns:
        (is_tandem, copies) where copies is a list of dicts with:
        {'id', 'seq', 'chrom', 'gstart', 'gend', 'strand', 'pident', 'qstart', 'qend'}
        Returns (False, []) if not tandem duplication.
    """
    if not hits or len(hits) < 2:
        return False, []
    
    query_len = len(query_seq)
    if query_len == 0:
        return False, []
    
    # Only consider deduplicated hits
    hits = _deduplicate_hits(hits, overlap_threshold=0.8)
    
    # Check how many hits cover a large fraction of the query
    # An "exon" would cover a SMALL fraction; a tandem copy covers a LARGE fraction
    large_coverage_hits = []
    for h in hits:
        q_cov = (h['qend'] - h['qstart'] + 1) / query_len
        if q_cov >= 0.35:  # hit covers >=35% of query → could be full copy
            large_coverage_hits.append(h)
    
    # Also check: do most hits overlap in query space?
    # Sort by qstart
    sorted_hits = sorted(hits, key=lambda h: h['qstart'])
    query_overlap_count = 0
    for i in range(len(sorted_hits) - 1):
        for j in range(i + 1, len(sorted_hits)):
            q_overlap = _calc_overlap(
                sorted_hits[i]['qstart'], sorted_hits[i]['qend'],
                sorted_hits[j]['qstart'], sorted_hits[j]['qend']
            )
            shorter = min(
                sorted_hits[i]['qend'] - sorted_hits[i]['qstart'] + 1,
                sorted_hits[j]['qend'] - sorted_hits[j]['qstart'] + 1
            )
            if shorter > 0 and q_overlap / shorter > 0.5:
                query_overlap_count += 1
    
    total_pairs = len(sorted_hits) * (len(sorted_hits) - 1) / 2
    overlap_fraction = query_overlap_count / total_pairs if total_pairs > 0 else 0
    
    # Decision: if most hits cover large parts of query AND overlap in query space
    # → tandem duplication
    is_tandem = False
    
    if len(large_coverage_hits) >= 2:
        # Multiple hits each covering >35% of query → likely tandem copies
        is_tandem = True
        print(f"[Tandem] {len(large_coverage_hits)} hits each cover >35% of query → tandem copies")
    elif overlap_fraction > 0.5 and len(hits) >= 3:
        # Most hit pairs overlap in query space → likely tandem copies, not exons
        is_tandem = True
        print(f"[Tandem] {overlap_fraction:.0%} of hit pairs overlap in query space → tandem copies")
    
    if not is_tandem:
        return False, []
    
    # Extract each copy as a separate protein
    # Determine strand consensus
    strand_votes = {'+': 0, '-': 0}
    for h in hits:
        strand_votes[h.get('strand', '+')] += 1
    consensus_strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'
    
    copies = []
    # Use the large-coverage hits as the basis
    candidate_hits = large_coverage_hits if large_coverage_hits else hits
    
    # Sort by genomic position
    candidate_hits.sort(key=lambda h: h['gstart'])
    
    # Cluster nearby hits into individual gene copies
    gene_clusters = []
    current_cluster = [candidate_hits[0]]
    
    for h in candidate_hits[1:]:
        prev = current_cluster[-1]
        # If hits are very close in genome space (<5kb), they're likely same gene
        if h['gstart'] - prev['gend'] < 5000:
            current_cluster.append(h)
        else:
            gene_clusters.append(current_cluster)
            current_cluster = [h]
    gene_clusters.append(current_cluster)
    
    for copy_num, cluster in enumerate(gene_clusters, start=1):
        # Merge cluster hits into one region
        gstart = min(h['gstart'] for h in cluster)
        gend = max(h['gend'] for h in cluster)
        best_hit = min(cluster, key=lambda h: h['evalue'])
        
        # Extract and translate the DNA
        exon_dna = chrom_seq[gstart:gend]
        strand = best_hit.get('strand', consensus_strand)
        if strand == '-':
            exon_dna = reverse_complement(exon_dna)
        
        # Trim to codon boundary
        exon_dna = exon_dna[:len(exon_dna) - len(exon_dna) % 3]
        if len(exon_dna) < 9:
            continue
        
        prot = translate(exon_dna)
        # Stop at first stop codon
        if '*' in prot:
            prot = prot.split('*')[0]
        
        if len(prot) < 10:
            continue
        
        copies.append({
            'id': f"GOI_copy_{copy_num}",
            'seq': prot,
            'chrom': chrom_name,
            'gstart': gstart,
            'gend': gend,
            'strand': strand,
            'pident': best_hit['pident'],
            'evalue': best_hit['evalue'],
            'qstart': best_hit['qstart'],
            'qend': best_hit['qend'],
            'exon_num': copy_num,
            'has_start_codon': False,
            'has_stop_codon': False,
            'splice_donor': None,
            'splice_acceptor': None,
            'coords': (gstart, gend),
        })
    
    if len(copies) >= 2:
        print(f"[Tandem] Identified {len(copies)} tandem copies "
              f"spanning {copies[0]['gstart']}-{copies[-1]['gend']} on {chrom_name}")
        return True, copies
    
    return False, []


# =============================================================================
# HIT-BASED EXON ANNOTATION (Scenario A - No GFF)
# =============================================================================

def _filter_exon_hits(hits, query_len, min_query_cov, min_alnlen):
    """Filter exon hits by query coverage and alignment length."""
    if not hits:
        return hits
    if query_len <= 0:
        return hits
    kept = []
    for h in hits:
        qstart = h.get('qstart')
        qend = h.get('qend')
        if not qstart or not qend:
            continue
        qspan = abs(qend - qstart) + 1
        alnlen = h.get('alnlen', qspan)
        qcov = qspan / query_len if query_len > 0 else 0
        if min_alnlen and alnlen < min_alnlen:
            continue
        if min_query_cov and qcov < min_query_cov:
            continue
        kept.append(h)
    return kept



# =============================================================================
# MINIPROT-BASED GENE MODELING
# =============================================================================

def _check_miniprot():
    """Check if miniprot is available."""
    try:
        result = subprocess.run(["miniprot", "--version"],
                                capture_output=True, text=True, timeout=5)
        return result.returncode == 0 or result.stderr.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

MINIPROT_AVAILABLE = _check_miniprot()


def _parse_miniprot_gff(gff_text):
    """
    Parse miniprot GFF3 output to extract gene models.
    
    Returns list of gene models, each with:
      - rank, identity, score
      - list of CDS features with (start, end, strand, phase, target_info)
    """
    models = {}  # ID -> model dict
    
    for line in gff_text.strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 9:
            continue
        
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts
        
        # Parse attributes
        attr_dict = {}
        for kv in attrs.split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                attr_dict[k] = v
        
        if ftype in ('mRNA', 'transcript'):
            model_id = attr_dict.get('ID', '')
            rank = int(attr_dict.get('Rank', 1))
            identity = float(attr_dict.get('Identity', 0))
            score_val = float(score) if score != '.' else 0
            models[model_id] = {
                'id': model_id,
                'rank': rank,
                'identity': identity,
                'score': score_val,
                'strand': strand,
                'cds_list': []
            }
        elif ftype == 'CDS':
            parent = attr_dict.get('Parent', '')
            if parent in models:
                # Keep CDS coordinates local to the provided target sequence.
                # The caller applies region offset exactly once when converting to global coords.
                cds_start = int(start) - 1  # GFF 1-based -> 0-based
                cds_end = int(end)          # inclusive -> exclusive
                cds_phase = int(phase) if phase.isdigit() else 0
                
                # Parse Target attribute (query protein coords)
                target = attr_dict.get('Target', '')
                qstart, qend = 0, 0
                if target:
                    target_parts = target.split()
                    if len(target_parts) >= 3:
                        qstart = int(target_parts[1])
                        qend = int(target_parts[2])
                
                models[parent]['cds_list'].append({
                    'gstart': cds_start,
                    'gend': cds_end,
                    'strand': strand,
                    'phase': cds_phase,
                    'qstart': qstart,
                    'qend': qend
                })
    
    # Sort models by rank (best first), then by score descending
    result = sorted(models.values(), key=lambda m: (m['rank'], -m['score']))
    # Sort CDS within each model by genomic position
    for model in result:
        model['cds_list'].sort(key=lambda c: c['gstart'])
    
    return result


def annotate_using_miniprot(query_seq, chrom_seq, chrom_name,
                             strand=None, max_intron=200000,
                             sensitive=False, region_offset=0):
    """
    Use miniprot to align query protein to genomic region.
    
    Args:
        query_seq: Query protein sequence (string)
        chrom_seq: Target genomic DNA sequence (string)
        chrom_name: Chromosome/scaffold name
        strand: Unused (miniprot searches both strands)
        max_intron: Max intron size (-G flag)
        sensitive: If True, use ultra-sensitive params (for GOI)
        region_offset: Offset to add to coordinates (if chrom_seq is a sub-region)
    
    Returns:
        (exon_list, full_protein) matching the old function signature.
        exon_list: list of dicts with keys:
            id, seq, exon_num, qstart, qend, gstart, gend, strand,
            chrom, pident, has_start_codon, has_stop_codon,
            splice_donor, splice_acceptor
        full_protein: concatenated translated CDS
    """
    if not query_seq or not chrom_seq:
        return [], query_seq or ''
    
    if not MINIPROT_AVAILABLE:
        print("[annotate_goi_exons] WARNING: miniprot not available, returning raw query",
              file=sys.stderr)
        return [], query_seq
    
    pid = os.getpid()
    query_file = f"/tmp/synterra_mp_query_{pid}.faa"
    target_file = f"/tmp/synterra_mp_target_{pid}.fna"
    
    try:
        # Write temp files
        write_fasta([("query", query_seq)], query_file)
        write_fasta([(chrom_name, chrom_seq)], target_file)
        
        # Build miniprot command
        cmd = [
            "miniprot",
            "--gff",                    # GFF3 output
            "--trans",                  # Also output translated protein
            "-G", str(max_intron),      # Max intron size
            "-t", "1",                  # Single thread (called per gene)
            "--outc=0.01",               # Min 1% query coverage for output
        ]
        
        query_len = len(query_seq)
        
        if sensitive:
            # Ultra-sensitive for GOI and divergent sequences
            cmd.extend([
                "-n", "2",              # Lower min syncmers for seeding
                "-p", "0.3",            # Lower secondary threshold
                "-N", "50",             # More secondary alignments
                "--outs=0.3",           # Output if >= 30% of best score
            ])
        
        if query_len < 60:
            # Short peptide mode (e.g., toxins)
            cmd.extend([
                "-L", "10",             # Min ORF length 10aa (default 30)
            ])
            if query_len < 30:
                cmd.extend([
                    "-S",               # No splicing for very short peptides
                    "-G", "2000",       # Override max intron to small value
                ])
        
        cmd.extend([target_file, query_file])
        
        # Run miniprot
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        if result.returncode != 0:
            print(f"[miniprot] Warning: non-zero exit ({result.returncode}): "
                  f"{result.stderr[:200]}", file=sys.stderr)
            return [], query_seq
        
        # Parse GFF output
        gff_output = result.stdout
        if not gff_output.strip():
            print("[miniprot] No alignments found", file=sys.stderr)
            return [], query_seq
        
        models = _parse_miniprot_gff(gff_output)
        
        if not models:
            print("[miniprot] No gene models parsed from GFF", file=sys.stderr)
            return [], query_seq
        
        # Take the best model (rank 1)
        best = models[0]
        
        if not best['cds_list']:
            return [], query_seq
        
        # Extract CDS DNA and build exon list
        exons = []
        dna_fragments = []
        model_strand = best['strand']
        
        cds_ordered = best['cds_list']
        if model_strand == '-':
            cds_ordered = list(reversed(cds_ordered))
        
        for i, cds in enumerate(cds_ordered, 1):
            gs = cds['gstart']
            ge = cds['gend']
            
            # Extract DNA
            exon_dna = chrom_seq[gs:ge]
            if model_strand == '-':
                exon_dna = reverse_complement(exon_dna)
            
            # Handle phase (trim N bases from start if phase > 0)
            phase = cds['phase']
            if phase > 0 and i == 1:
                exon_dna = exon_dna[phase:]
            
            dna_fragments.append(exon_dna)
            
            # Check splice sites
            splice_donor = None
            splice_acceptor = None
            if gs >= 2:
                if model_strand == '+':
                    splice_acceptor = chrom_seq[gs-2:gs].upper()
                else:
                    splice_donor = reverse_complement(chrom_seq[gs-2:gs]).upper()
            if ge + 2 <= len(chrom_seq):
                if model_strand == '+':
                    splice_donor = chrom_seq[ge:ge+2].upper()
                else:
                    splice_acceptor = reverse_complement(chrom_seq[ge:ge+2]).upper()
            
            exon_prot = translate(exon_dna)
            
            exons.append({
                'id': f'exon_{i}',
                # Remove all stop codons; internal stops can appear in noisy models
                # and degrade downstream protein->DNA searches.
                'seq': exon_prot.replace('*', ''),
                'exon_num': i,
                'qstart': cds['qstart'],
                'qend': cds['qend'],
                'gstart': gs + region_offset,
                'gend': ge + region_offset,
                'strand': model_strand,
                'chrom': chrom_name,
                'pident': best['identity'] * 100,  # miniprot reports 0-1
                'has_start_codon': (i == 1 and exon_prot.startswith('M')),
                'has_stop_codon': (i == len(cds_ordered) and
                                   exon_dna[-3:].upper() in ('TAA', 'TAG', 'TGA')),
                'splice_donor': splice_donor,
                'splice_acceptor': splice_acceptor,
                'method': 'miniprot'
            })
        
        # Assemble full protein from DNA
        full_cds = ''.join(dna_fragments)
        # Trim to codon boundary
        full_cds = full_cds[:len(full_cds) - len(full_cds) % 3]
        # Keep a stop-free protein for downstream querying.
        full_protein = translate(full_cds).replace('*', '')
        
        if not full_protein:
            full_protein = query_seq
        
        # Parse --trans output for verification (appears as PAF comments)
        # The GFF-based extraction above is authoritative
        
        print(f"[miniprot] Best model: {len(exons)} exon(s), "
              f"{len(full_protein)} aa, identity={best['identity']:.2f}, "
              f"score={best['score']:.0f}", file=sys.stderr)
        
        return exons, full_protein
        
    except subprocess.TimeoutExpired:
        print("[miniprot] Timed out after 120s", file=sys.stderr)
        return [], query_seq
    except Exception as e:
        print(f"[miniprot] Error: {e}", file=sys.stderr)
        return [], query_seq
    finally:
        for f in [query_file, target_file]:
            if os.path.exists(f):
                os.remove(f)


def annotate_exons_from_hit_list(hits, query_seq, chrom_seq, chrom_name,
                                  search_missing=True,
                                  gap_min_size=10,
                                  gap_search_window=50000,
                                  gap_evalue=10.0,
                                  gap_min_identity=25.0,
                                  gap_min_alnlen=10,
                                  gap_max_hits=5,
                                  exon_query_mode=False,
                                  min_exon_query_cov=0.25,
                                  min_exon_alnlen=30,
                                  sensitive=False):
    """
    Entry point for exon annotation - uses miniprot for gene modeling.
    
    Kept for backward compatibility with iterative_search_runner.py.
    The 'hits' parameter is used to determine the region and strand,
    then miniprot does the actual gene modeling.
    """
    if not hits or not query_seq or not chrom_seq:
        return [], query_seq or ''

    # Determine strand consensus from hits (informational only)
    strand_votes = {'+': 0, '-': 0}
    for h in hits:
        strand_votes[h.get('strand', '+')] += 1
    consensus_strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'
    
    # Determine search region from hits
    hit_starts = [h.get('gstart', h.get('start', 0)) for h in hits]
    hit_ends = [h.get('gend', h.get('end', 0)) for h in hits]
    
    if not hit_starts or not hit_ends:
        return [], query_seq
    
    region_start = max(0, min(hit_starts) - gap_search_window)
    region_end = min(len(chrom_seq), max(hit_ends) + gap_search_window)
    
    # Extract sub-region for miniprot
    sub_seq = chrom_seq[region_start:region_end]
    
    # Detect if this is likely a GOI query (use sensitive mode)
    is_goi = sensitive or any('GOI' in str(h.get('query', '')) for h in hits)
    
    # Calculate max intron from hit spread
    max_intron = max(20000, region_end - region_start)
    
    # Run miniprot
    exons, full_protein = annotate_using_miniprot(
        query_seq, sub_seq, chrom_name,
        strand=consensus_strand,
        max_intron=max_intron,
        sensitive=is_goi,
        region_offset=region_start
    )
    
    return exons, full_protein


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Annotate GOI exons in home genome"
    )
    parser.add_argument("--query", required=True,
                        help="GOI protein sequence (FASTA)")
    parser.add_argument("--genome", required=True,
                        help="Home genome (FASTA)")
    parser.add_argument("--gff", default="NO_GFF",
                        help="Home GFF annotation file (or 'NO_GFF')")
    parser.add_argument("--blast_hits", default=None,
                        help="Raw tblastn hits from LOCATE_GENE (format 6)")
    parser.add_argument("--mmseqs_hits", default=None,
                        help="Raw MMseqs2 hits from LOCATE_GENE (format 6)")
    parser.add_argument("--query_id", default="",
                        help="Query ID (UniProt accession or gene name)")
    parser.add_argument("--gff_search_window", type=int, default=100000,
                        help="Window around hits for GFF-based search")
    parser.add_argument("--gap_search_window", type=int, default=50000,
                        help="Window around existing exons for gap search")
    parser.add_argument("--gap_min_size", type=int, default=10,
                        help="Minimum gap size (aa) to search")
    parser.add_argument("--gap_evalue", type=float, default=10.0,
                        help="E-value threshold for gap search")
    parser.add_argument("--gap_min_identity", type=float, default=25.0,
                        help="Minimum identity for gap hits")
    parser.add_argument("--gap_min_alnlen", type=int, default=10,
                        help="Minimum alignment length for gap hits")
    parser.add_argument("--gap_max_hits", type=int, default=5,
                        help="Maximum hits to consider per gap search")
    parser.add_argument("--min_exon_query_cov", type=float, default=0.25,
                        help="Minimum query coverage for exon hits")
    parser.add_argument("--min_exon_alnlen", type=int, default=30,
                        help="Minimum alignment length for exon hits")
    parser.add_argument("--output_exons", required=True,
                        help="Output: exon protein FASTA")
    parser.add_argument("--output_bed", required=True,
                        help="Output: exon locations BED")
    parser.add_argument("--output_info", required=True,
                        help="Output: annotation info JSON")

    args = parser.parse_args()

    # Load query sequence
    query_records = list(parse_fasta(args.query))
    if not query_records:
        print("ERROR: No sequences in query file!", file=sys.stderr)
        sys.exit(1)

    query_header, query_clean_id, query_seq = query_records[0]
    query_id = args.query_id if args.query_id else query_clean_id

    print(f"[annotate_goi_exons] GOI: {query_id} ({len(query_seq)} aa)")

    gff_available = (args.gff != "NO_GFF" and
                     os.path.exists(args.gff) and
                     os.path.getsize(args.gff) > 0)

    exons = []
    full_protein = query_seq
    method_used = "none"

    # ========== Try GFF-based approach first ==========
    if gff_available:
        print(f"[annotate_goi_exons] GFF available, trying annotation-based approach...")

        match = match_goi_in_gff(
            query_id,
            query_seq,
            args.gff,
            args.genome,
            args.blast_hits,
            args.mmseqs_hits,
            args.gff_search_window
        )

        if match:
            result = extract_exons_from_gff_match(match, args.genome)
            if result:
                exons, full_protein = result
                method_used = "gff_annotation"
                print(f"[annotate_goi_exons] Extracted {len(exons)} exons from GFF")

    # ========== Fall back to hit-based annotation ==========
    if not exons:
        print(f"[annotate_goi_exons] Using hit-based exon annotation...")

        if args.blast_hits or args.mmseqs_hits:
            # Parse hit regions to determine where to search
            hit_regions = _parse_hit_regions(args.blast_hits, args.mmseqs_hits)
            if hit_regions:
                # Load target genome
                genome_data = load_genome(args.genome)
                
                # Group hits by chromosome
                chroms = {}
                for hr in hit_regions:
                    c = hr['chrom']
                    if c not in chroms:
                        chroms[c] = []
                    chroms[c].append(hr)
                
                # For each chromosome with hits, run miniprot
                for chrom, chrom_hits in chroms.items():
                    if chrom not in genome_data:
                        continue
                    chrom_seq = genome_data[chrom]
                    
                    # Determine region
                    region_start = max(0, min(h['start'] for h in chrom_hits) - args.gap_search_window)
                    region_end = min(len(chrom_seq), max(h['end'] for h in chrom_hits) + args.gap_search_window)
                    sub_seq = chrom_seq[region_start:region_end]
                    
                    exons, full_protein = annotate_using_miniprot(
                        query_seq, sub_seq, chrom,
                        max_intron=max(20000, region_end - region_start),
                        sensitive=True,  # GOI always gets sensitive
                        region_offset=region_start
                    )
                    if exons:
                        break  # Take first chromosome with results
            # Detect if tandem duplication was found (copies have id "GOI_copy_N")
            if exons and any(e['id'].startswith('GOI_copy_') for e in exons):
                method_used = "tandem_duplication"
            else:
                method_used = "hit_annotation"
        else:
            print("[annotate_goi_exons] No hits provided, using full protein only")
            method_used = "full_protein_only"

    # ========== Build output ==========
    fasta_records = []
    bed_records = []

    # Always include the full GOI protein
    goi_full_id = f"GOI_{query_id}"
    fasta_records.append((goi_full_id, query_seq))

    # Add individual exon sequences
    # Normalize IDs so all GOI-derived entries are explicitly prefixed with GOI_.
    for idx, exon in enumerate(exons, start=1):
        raw_exon_id = exon.get('id', '') or f"exon_{idx}"
        exon_num = exon.get('exon_num', idx)
        exon_id = raw_exon_id
        if not (raw_exon_id.startswith("GOI_") or raw_exon_id.startswith("GOI_copy_")):
            if raw_exon_id.startswith("exon_"):
                exon_id = f"{goi_full_id}|{raw_exon_id}"
            else:
                exon_id = f"{goi_full_id}|exon_{exon_num}"
        exon['id'] = exon_id
        fasta_records.append((exon_id, exon['seq']))

        # BED record
        if 'gstart' in exon:
            bed_records.append(
                f"{exon['chrom']}\t{exon['gstart']}\t{exon['gend']}\t"
                f"{exon_id}\t{exon.get('pident', 0):.1f}\t{exon['strand']}"
            )
        elif 'coords' in exon:
            bed_records.append(
                f"{exon['chrom']}\t{exon['coords'][0]}\t{exon['coords'][1]}\t"
                f"{exon_id}\t0\t{exon.get('strand', '.')}"
            )

    # Write FASTA output
    write_fasta(fasta_records, args.output_exons)

    # Write BED output
    with open(args.output_bed, 'w') as f:
        for rec in bed_records:
            f.write(rec + '\n')

    # Write info JSON
    info = {
        'query_id': query_id,
        'query_length': len(query_seq),
        'method': method_used,
        'num_exons': len(exons),
        'exons': [{
            'id': e['id'],
            'length_aa': len(e['seq']),
            'exon_num': e.get('exon_num', 0),
            'pident': e.get('pident', 0),
            'has_start_codon': e.get('has_start_codon', False),
            'has_stop_codon': e.get('has_stop_codon', False),
            'splice_donor': e.get('splice_donor', None),
            'splice_acceptor': e.get('splice_acceptor', None)
        } for e in exons],
        'full_protein_length': len(full_protein)
    }

    with open(args.output_info, 'w') as f:
        json.dump(info, f, indent=2)

    print(f"\n[annotate_goi_exons] SUMMARY:")
    print(f"  Method: {method_used}")
    print(f"  Exons found: {len(exons)}")
    print(f"  Full protein: {len(full_protein)} aa")
    print(f"  Output records: {len(fasta_records)} (1 full + {len(exons)} exons)")

    for exon in exons:
        splice = ""
        if exon.get('splice_acceptor'):
            splice += f" acceptor={exon['splice_acceptor']}"
        if exon.get('splice_donor'):
            splice += f" donor={exon['splice_donor']}"
        start_stop = ""
        if exon.get('has_start_codon'):
            start_stop += " [START]"
        if exon.get('has_stop_codon'):
            start_stop += " [STOP]"

        print(f"  Exon {exon.get('exon_num', '?')}: {len(exon['seq'])} aa, "
              f"query {exon.get('qstart', '?')}-{exon.get('qend', '?')}"
              f"{start_stop}{splice}")


if __name__ == "__main__":
    main()
