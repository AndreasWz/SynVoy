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

def annotate_exons_from_hit_list(hits, query_seq, chrom_seq, chrom_name,
                                  search_missing=True,
                                  gap_min_size=10,
                                  gap_search_window=50000,
                                  gap_evalue=10.0,
                                  gap_min_identity=25.0,
                                  gap_min_alnlen=10,
                                  gap_max_hits=5):
    """
    Core exon annotation from pre-parsed hit dicts.

    Reusable entry point for both home-genome (annotate_goi_exons.py) and
    target-genome (iterative_search_runner.py) annotation.

    Each hit dict should have at minimum:
        qstart, qend:  query protein positions (1-based)
        gstart, gend:  genomic positions within chrom_seq (0-based)
        evalue:        e-value
        pident:        percent identity
    Optional hit fields:
        alnlen, strand, chrom, bitscore

    Args:
        hits:           list of hit dicts
        query_seq:      full query protein sequence
        chrom_seq:      DNA sequence of the chromosome/region
        chrom_name:     chromosome/region name (for output IDs)
        search_missing: if True, run tblastn to search for missing exons

    Returns:
        (annotated_exons, full_protein) where annotated_exons is a list of
        exon dicts and full_protein is the original query sequence
    """
    if not hits or not query_seq or not chrom_seq:
        return [], query_seq or ''

    # Sort by e-value (best first)
    hits = sorted(hits, key=lambda h: h.get('evalue', 999))

    # Deduplicate overlapping hits
    hits = _deduplicate_hits(hits)

    if not hits:
        return [], query_seq

    print(f"[Exon Annotate] Working with {len(hits)} unique hits")

    query_len = len(query_seq)

    # Determine strand consensus
    strand_votes = {'+': 0, '-': 0}
    for h in hits:
        strand_votes[h.get('strand', '+')] += 1
    consensus_strand = '+' if strand_votes['+'] >= strand_votes['-'] else '-'

    # Sort by query position (tells us exon order)
    hits.sort(key=lambda h: h.get('qstart', 0))

    # Annotate each hit as a candidate exon
    annotated_exons = []
    for i, hit in enumerate(hits, start=1):
        exon = _annotate_single_exon(
            hit, chrom_seq, query_seq, consensus_strand,
            i, len(hits), chrom_name
        )
        if exon:
            annotated_exons.append(exon)

    if not annotated_exons:
        return [], query_seq

    # Check query coverage
    covered = [False] * query_len
    for exon in annotated_exons:
        for j in range(exon['qstart'] - 1, min(exon['qend'], query_len)):
            covered[j] = True

    coverage = sum(covered) / query_len * 100 if query_len > 0 else 0
    print(f"[Exon Annotate] Query coverage: {coverage:.1f}%")

    # Search for missing exons if requested
    if search_missing:
        missing_regions = _find_gaps(covered, min_gap_size=gap_min_size)
        if missing_regions:
            print(f"[Exon Annotate] {len(missing_regions)} gap(s), searching nearby...")
            new_exons = _search_missing_exons(
                query_seq, missing_regions, annotated_exons,
                None, chrom_name, chrom_seq, consensus_strand,
                gap_search_window=gap_search_window,
                gap_evalue=gap_evalue,
                gap_min_identity=gap_min_identity,
                gap_min_alnlen=gap_min_alnlen,
                gap_max_hits=gap_max_hits
            )
            annotated_exons.extend(new_exons)
            annotated_exons.sort(key=lambda e: e.get('qstart', 0))

    return annotated_exons, query_seq


def annotate_exons_from_hits(query_seq, genome_file, blast_hits_file, mmseqs_hits_file,
                             gap_min_size=10,
                             gap_search_window=50000,
                             gap_evalue=10.0,
                             gap_min_identity=25.0,
                             gap_min_alnlen=10,
                             gap_max_hits=5):
    """
    Use tblastn/MMseqs2 hit files to annotate individual exons of the GOI.

    Parses format-6 hit files, picks the best chromosome, then delegates
    to annotate_exons_from_hit_list() for the core annotation logic.

    Returns:
        list of exon dicts, full_protein string
    """
    print("[Hit Annotate] Annotating exons from protein→DNA hits...")

    all_hits = _parse_format6_hits(blast_hits_file, mmseqs_hits_file)
    if not all_hits:
        print("[Hit Annotate] No hits found!")
        return [], query_seq

    # Load genome
    genome_seqs = load_genome(genome_file)

    # Group hits by chromosome
    hits_by_chrom = {}
    for h in all_hits:
        chrom = h['chrom']
        if chrom not in hits_by_chrom:
            hits_by_chrom[chrom] = []
        hits_by_chrom[chrom].append(h)

    # Pick the chromosome with the most/best hits
    best_chrom = max(hits_by_chrom.keys(),
                     key=lambda c: (len(hits_by_chrom[c]),
                                    -min(h['evalue'] for h in hits_by_chrom[c])))

    chrom_seq = genome_seqs.get(best_chrom)
    if not chrom_seq:
        print(f"[Hit Annotate] Chromosome {best_chrom} not found!")
        return [], query_seq

    print(f"[Hit Annotate] Best chromosome: {best_chrom} with {len(hits_by_chrom[best_chrom])} hits")

    # Check for tandem duplications BEFORE exon annotation
    is_tandem, copies = detect_tandem_duplications(
        hits_by_chrom[best_chrom], query_seq, chrom_seq, best_chrom
    )
    
    if is_tandem:
        print(f"[Hit Annotate] TANDEM DUPLICATION detected: {len(copies)} copies")
        return copies, query_seq

    return annotate_exons_from_hit_list(
        hits_by_chrom[best_chrom],
        query_seq,
        chrom_seq,
        best_chrom,
        search_missing=True,
        gap_min_size=gap_min_size,
        gap_search_window=gap_search_window,
        gap_evalue=gap_evalue,
        gap_min_identity=gap_min_identity,
        gap_min_alnlen=gap_min_alnlen,
        gap_max_hits=gap_max_hits
    )


def _parse_format6_hits(blast_file, mmseqs_file):
    """Parse BLAST/MMseqs format 6 output files."""
    hits = []

    for hits_file in [blast_file, mmseqs_file]:
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
                        pident = float(parts[2])
                        alnlen = int(parts[3])

                        strand = '+' if tstart < tend else '-'
                        # Convert 1-based BLAST/mmseqs coords to 0-based half-open
                        gstart = min(tstart, tend) - 1  # 0-based start
                        gend = max(tstart, tend)         # exclusive end

                        hits.append({
                            'query': parts[0],
                            'chrom': parts[1],
                            'pident': pident,
                            'alnlen': alnlen,
                            'qstart': min(qstart, qend),
                            'qend': max(qstart, qend),
                            'gstart': gstart,
                            'gend': gend,
                            'strand': strand,
                            'evalue': evalue,
                            'bitscore': float(parts[11]),
                            'source': os.path.basename(hits_file)
                        })
                    except (ValueError, IndexError):
                        continue
        except Exception:
            continue

    return hits


def _deduplicate_hits(hits, overlap_threshold=0.5):
    """Remove overlapping hits, keeping the one with best e-value.
    
    Two-pass deduplication:
    1. Remove hits redundant in BOTH query AND genomic space
    2. Merge hits that overlap significantly in GENOMIC space only
       (from fragment-based search producing many overlapping hits)
    """
    if len(hits) <= 1:
        return hits

    # Sort by e-value
    hits.sort(key=lambda h: h['evalue'])

    # Pass 1: Remove hits redundant in both query AND genomic space
    kept = []
    for hit in hits:
        is_redundant = False
        for existing in kept:
            # Check query overlap
            q_overlap = _calc_overlap(
                hit['qstart'], hit['qend'],
                existing['qstart'], existing['qend']
            )
            q_span = hit['qend'] - hit['qstart'] + 1

            if q_span > 0 and q_overlap / q_span > overlap_threshold:
                # Also check genomic overlap
                g_overlap = _calc_overlap(
                    hit['gstart'], hit['gend'],
                    existing['gstart'], existing['gend']
                )
                g_span = hit['gend'] - hit['gstart'] + 1

                if g_span > 0 and g_overlap / g_span > overlap_threshold:
                    is_redundant = True
                    break

        if not is_redundant:
            kept.append(hit)

    # Pass 2: Merge hits that overlap significantly in GENOMIC space
    # This handles fragment-based search producing many hits covering
    # different query ranges but the same genomic region
    if len(kept) > 1:
        kept.sort(key=lambda h: h['evalue'])
        merged = []
        for hit in kept:
            is_genomic_dup = False
            for existing in merged:
                g_overlap = _calc_overlap(
                    hit['gstart'], hit['gend'],
                    existing['gstart'], existing['gend']
                )
                smaller_span = min(
                    hit['gend'] - hit['gstart'] + 1,
                    existing['gend'] - existing['gstart'] + 1
                )
                if smaller_span > 0 and g_overlap / smaller_span > 0.6:
                    # Merge: expand the existing hit to cover both ranges
                    existing['gstart'] = min(existing['gstart'], hit['gstart'])
                    existing['gend'] = max(existing['gend'], hit['gend'])
                    existing['qstart'] = min(existing['qstart'], hit['qstart'])
                    existing['qend'] = max(existing['qend'], hit['qend'])
                    is_genomic_dup = True
                    break
            if not is_genomic_dup:
                merged.append(hit)
        kept = merged

    return kept


def _calc_overlap(s1, e1, s2, e2):
    """Calculate overlap between two intervals."""
    return max(0, min(e1, e2) - max(s1, s2))


def _annotate_single_exon(hit, chrom_seq, query_seq, strand, exon_num, total_hits, chrom):
    """
    Annotate a single exon from a tblastn/MMseqs hit.

    Checks:
    - Splice sites (GT-AG at intron boundaries)
    - Start codon (first exon)
    - Stop codon (last exon)
    - Reading frame consistency
    """
    gstart = hit['gstart']
    gend = hit['gend']

    # Expand slightly to capture full codons at boundaries
    # tblastn coordinates are nucleotide positions
    # Add a small buffer for splice site checking
    buffer = 10
    check_start = max(0, gstart - buffer)
    check_end = min(len(chrom_seq), gend + buffer)

    region = chrom_seq[check_start:check_end]

    # Check splice sites
    splice_info = _check_splice_sites(chrom_seq, gstart, gend, strand)

    # Extract the hit region DNA
    exon_dna = chrom_seq[gstart:gend]
    if strand == '-':
        exon_dna = reverse_complement(exon_dna)

    # Translate
    exon_dna_trimmed = exon_dna[:len(exon_dna) - len(exon_dna) % 3]
    if len(exon_dna_trimmed) < 9:
        return None

    exon_prot = translate(exon_dna_trimmed).replace('*', '')
    if len(exon_prot) < 3:
        return None

    # Check start codon (first exon should have ATG)
    has_start = False
    if exon_num == 1:
        if strand == '+':
            first_codon = chrom_seq[gstart:gstart + 3].upper()
        else:
            first_codon = reverse_complement(chrom_seq[gend - 3:gend]).upper()
        has_start = first_codon == 'ATG'

    # Check stop codon (last exon should end with stop)
    has_stop = False
    if exon_num == total_hits:
        if strand == '+':
            last_codon = chrom_seq[gend:gend + 3].upper()
        else:
            last_codon = reverse_complement(chrom_seq[gstart - 3:gstart]).upper()
        has_stop = last_codon in ['TAA', 'TAG', 'TGA']

    return {
        'id': f"GOI_hit|exon_{exon_num}",
        'seq': exon_prot,
        'exon_num': exon_num,
        'total_exons': total_hits,
        'chrom': chrom,
        'gstart': gstart,
        'gend': gend,
        'strand': strand,
        'qstart': hit['qstart'],
        'qend': hit['qend'],
        'pident': hit['pident'],
        'evalue': hit['evalue'],
        'has_start_codon': has_start,
        'has_stop_codon': has_stop,
        'splice_donor': splice_info.get('donor', None),
        'splice_acceptor': splice_info.get('acceptor', None),
        'coords': (gstart, gend)
    }


def _check_splice_sites(chrom_seq, gstart, gend, strand):
    """
    Check for canonical splice sites flanking the exon.

    Eukaryotic introns:
    - Donor (5' end of intron): GT (>98% of cases)
    - Acceptor (3' end of intron): AG (>98% of cases)

    For + strand exon at [gstart, gend]:
    - Splice acceptor: 2bp before gstart (should be AG)
    - Splice donor: 2bp after gend (should be GT)

    For - strand, it's reversed.
    """
    result = {}

    if strand == '+':
        # Check acceptor site (AG before exon start)
        if gstart >= 2:
            acceptor = chrom_seq[gstart - 2:gstart].upper()
            result['acceptor'] = acceptor
            result['acceptor_canonical'] = acceptor == 'AG'

        # Check donor site (GT after exon end)
        if gend + 2 <= len(chrom_seq):
            donor = chrom_seq[gend:gend + 2].upper()
            result['donor'] = donor
            result['donor_canonical'] = donor == 'GT'
    else:
        # Minus strand: complement everything
        # Acceptor (AG) appears as CT on + strand, after the exon end
        if gend + 2 <= len(chrom_seq):
            acceptor_raw = chrom_seq[gend:gend + 2].upper()
            result['acceptor'] = reverse_complement(acceptor_raw).upper()
            result['acceptor_canonical'] = result['acceptor'] == 'AG'

        # Donor (GT) appears as AC on + strand, before the exon start
        if gstart >= 2:
            donor_raw = chrom_seq[gstart - 2:gstart].upper()
            result['donor'] = reverse_complement(donor_raw).upper()
            result['donor_canonical'] = result['donor'] == 'GT'

    return result


def _find_gaps(covered, min_gap_size=10):
    """Find gaps in coverage array."""
    gaps = []
    in_gap = False
    gap_start = 0

    for i, c in enumerate(covered):
        if not c and not in_gap:
            in_gap = True
            gap_start = i
        elif c and in_gap:
            in_gap = False
            if i - gap_start >= min_gap_size:
                gaps.append((gap_start, i))

    # Final gap
    if in_gap and len(covered) - gap_start >= min_gap_size:
        gaps.append((gap_start, len(covered)))

    return gaps


def _search_missing_exons(query_seq, missing_regions, existing_exons,
                          genome_file, chrom, chrom_seq, strand,
                          gap_search_window=50000,
                          gap_evalue=10.0,
                          gap_min_identity=25.0,
                          gap_min_alnlen=10,
                          gap_max_hits=5):
    """
    Search for missing parts of the GOI near existing exon hits.

    For each gap in query coverage:
    1. Extract the missing protein subsequence
    2. Search nearby genomic regions (within 50kb of existing exons)
    3. Use tblastn (protein → DNA)
    4. Annotate any new hits

    Always uses protein → DNA search!
    """
    new_exons = []

    if not existing_exons:
        return new_exons

    # Define search region: around existing exons
    all_gstarts = [e['gstart'] for e in existing_exons]
    all_gends = [e['gend'] for e in existing_exons]
    region_start = max(0, min(all_gstarts) - gap_search_window)
    region_end = min(len(chrom_seq), max(all_gends) + gap_search_window)

    next_exon_num = max(e['exon_num'] for e in existing_exons) + 1

    for gap_start, gap_end in missing_regions:
        gap_protein = query_seq[gap_start:gap_end]
        if len(gap_protein) < gap_min_alnlen:
            continue

        print(f"[Missing Exon] Searching for query positions {gap_start + 1}-{gap_end} "
              f"({len(gap_protein)} aa) near existing exons...")

        # Write gap protein to temp file
        gap_query_file = f"/tmp/synterra_gap_query_{os.getpid()}.faa"
        gap_target_file = f"/tmp/synterra_gap_target_{os.getpid()}.fna"
        gap_hits_file = f"/tmp/synterra_gap_hits_{os.getpid()}.txt"

        try:
            write_fasta([("gap_query", gap_protein)], gap_query_file)

            # Extract nearby genomic region
            search_region = chrom_seq[region_start:region_end]
            write_fasta([(f"{chrom}_{region_start}_{region_end}", search_region)], gap_target_file)

            # Run tblastn (protein → DNA search!)
            cmd = [
                "tblastn",
                "-query", gap_query_file,
                "-subject", gap_target_file,
                "-out", gap_hits_file,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                "-evalue", str(gap_evalue),  # Relaxed for short fragments
                "-seg", "no",  # Don't mask low complexity (short queries)
                "-max_target_seqs", str(gap_max_hits)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

            if result.returncode == 0 and os.path.exists(gap_hits_file):
                gap_hits = _parse_format6_hits(gap_hits_file, None)

                for hit in gap_hits:
                    if hit['pident'] < gap_min_identity or hit['alnlen'] < gap_min_alnlen:
                        continue

                    # Convert coordinates back to global
                    actual_gstart = region_start + hit['gstart']
                    actual_gend = region_start + hit['gend']

                    # Check if this overlaps with existing exons
                    overlaps = False
                    for existing in existing_exons:
                        if (_calc_overlap(actual_gstart, actual_gend,
                                          existing['gstart'], existing['gend']) > 10):
                            overlaps = True
                            break

                    if overlaps:
                        continue

                    # Extract and annotate
                    exon_dna = chrom_seq[actual_gstart:actual_gend]
                    if strand == '-':
                        exon_dna = reverse_complement(exon_dna)

                    exon_dna = exon_dna[:len(exon_dna) - len(exon_dna) % 3]
                    if len(exon_dna) < 9:
                        continue

                    exon_prot = translate(exon_dna).replace('*', '')
                    if len(exon_prot) < 3:
                        continue

                    # Check splice sites
                    splice_info = _check_splice_sites(
                        chrom_seq, actual_gstart, actual_gend, strand
                    )

                    new_exons.append({
                        'id': f"GOI_hit|exon_{next_exon_num}",
                        'seq': exon_prot,
                        'exon_num': next_exon_num,
                        'total_exons': -1,  # Unknown
                        'chrom': chrom,
                        'gstart': actual_gstart,
                        'gend': actual_gend,
                        'strand': strand,
                        'qstart': gap_start + hit['qstart'],
                        'qend': gap_start + hit['qend'],
                        'pident': hit['pident'],
                        'evalue': hit['evalue'],
                        'has_start_codon': False,
                        'has_stop_codon': False,
                        'splice_donor': splice_info.get('donor', None),
                        'splice_acceptor': splice_info.get('acceptor', None),
                        'coords': (actual_gstart, actual_gend),
                        'method': 'gap_search'
                    })

                    print(f"[Missing Exon] Found new exon at {chrom}:{actual_gstart}-{actual_gend} "
                          f"({hit['pident']:.1f}% identity)")
                    next_exon_num += 1
                    break  # Take first good hit per gap

        except subprocess.TimeoutExpired:
            print(f"[Missing Exon] tblastn timed out for gap query")
        except FileNotFoundError:
            print(f"[Missing Exon] tblastn not found, skipping gap search")
        except Exception as e:
            print(f"[Missing Exon] Gap search failed: {e}")
        finally:
            for f in [gap_query_file, gap_target_file, gap_hits_file]:
                if os.path.exists(f):
                    os.remove(f)

    return new_exons


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
            exons, full_protein = annotate_exons_from_hits(
                query_seq,
                args.genome,
                args.blast_hits,
                args.mmseqs_hits,
                gap_min_size=args.gap_min_size,
                gap_search_window=args.gap_search_window,
                gap_evalue=args.gap_evalue,
                gap_min_identity=args.gap_min_identity,
                gap_min_alnlen=args.gap_min_alnlen,
                gap_max_hits=args.gap_max_hits
            )
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
    for exon in exons:
        exon_id = exon['id']
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
