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
import tempfile
import urllib.request
import urllib.error

try:
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, load_genome_subset, reverse_complement, translate, parse_gff, get_feature_id,
        parse_gff_attributes, get_genetic_code, set_genetic_code, stop_codons
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        parse_fasta, write_fasta, extract_id, extract_base_id,
        load_genome, load_genome_subset, reverse_complement, translate, parse_gff, get_feature_id,
        parse_gff_attributes, get_genetic_code, set_genetic_code, stop_codons
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
                     gff_search_window, gene_name=''):
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
    if gene_name:
        # Given by the user (--home_goi_gene): tried first.
        candidate_names.append(gene_name)

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
    match = _search_gff_by_name(gff_file, candidate_names,
                                query_seq=query_seq, genome_file=genome_file)
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


# A name match only POINTS at a gene; it does not show that the gene is the query.
# Every gene taken from the annotation, by name or by position, is therefore compared
# with the query protein, and the transcript that matches the query best is used.
# Without that check a query whose only UniProt name is a description ("Putative
# uncharacterized protein FLJ37218") matched 12,484 human genes on the word
# "uncharacterized", and the first one in the file -- a non-coding gene on another
# chromosome -- became the home gene: the run then searched every target genome for an
# 11-residue translation of that gene's span. Without the transcript choice, TP53 was
# read from the first transcript in the file, the 354-residue isoform.
HOME_GENE_MIN_IDENTITY = 60.0      # percent identity of the local alignment
HOME_GENE_MIN_SCORE_FRAC = 0.5     # alignment score / self-score of the QUERY: the gene explains half of the query
HOME_GENE_MAX_NAME_CANDIDATES = 50  # name matches compared with the query, best name score first


def _cds_protein(cds_parts, chrom_seq, strand):
    """Protein of one transcript: its CDS rows joined in coding order, read in the
    frame given by the phase of the first coding row, up to the first stop codon."""
    parts = sorted(cds_parts, key=lambda x: x['start'])
    if not parts:
        return ''
    first = parts[0] if strand == '+' else parts[-1]
    try:
        phase = int(first.get('phase', 0) or 0)
    except (TypeError, ValueError):
        phase = 0
    dna = ''.join(chrom_seq[part['start']:part['end']] for part in parts)
    if strand == '-':
        dna = reverse_complement(dna)
    if phase > 0:
        dna = dna[phase:]
    dna = dna[:len(dna) - len(dna) % 3]
    if len(dna) < 9:
        return ''
    return translate(dna).split('*')[0]


def _query_match(query_seq, protein):
    """How well an annotated protein agrees with the query.

    Returns ``(score_fraction, identity_percent)``: the Smith-Waterman score divided by
    the self-score of the QUERY, and the identity of the aligned columns. The fraction is
    taken of the query, not of the shorter sequence: a six-residue translation of a
    non-coding gene aligns perfectly to some six residues of any query. Without parasail
    both values come from the shared-3-mer measure, so the step still runs."""
    if not query_seq or not protein:
        return 0.0, 0.0
    try:
        from sequence_utils import sw_align
        score, ident = sw_align(query_seq, protein)
        self_q, _ = sw_align(query_seq, query_seq)
        return (score / self_q if self_q > 0 else 0.0), ident
    except (ImportError, RuntimeError):
        approx = _quick_protein_identity(query_seq, protein)
        return approx / 100.0, approx


def _gene_transcripts(gene_id, transcript_to_gene, cds_by_parent):
    """``[(transcript_id, cds_parts)]`` of a gene, in a fixed order: CDS rows attached
    to the gene itself first, then its transcripts sorted by id."""
    found = []
    if cds_by_parent.get(gene_id):
        found.append((gene_id, cds_by_parent[gene_id]))
    for tid in sorted(t for t, g in transcript_to_gene.items() if g == gene_id):
        if cds_by_parent.get(tid):
            found.append((tid, cds_by_parent[tid]))
    return found


def _best_transcript(gene_id, strand, chrom_seq, transcript_to_gene, cds_by_parent,
                     query_seq, measure=None):
    """The transcript of ``gene_id`` whose protein agrees best with the query.

    ``measure(query, protein)`` returns a tuple that is compared as a whole; by default
    it is ``_query_match``. Ties go to the protein whose length is closest to the
    query's (a query that is itself a short isoform is contained in the long transcript
    just as well), then to the first transcript in the fixed order of
    ``_gene_transcripts``. Returns
    ``(cds_parts, protein, measure_value)`` or ``None`` when the gene has no CDS."""
    measure = measure or _query_match
    best = None
    for tid, parts in _gene_transcripts(gene_id, transcript_to_gene, cds_by_parent):
        protein = _cds_protein(parts, chrom_seq, strand)
        if not protein:
            continue
        value = measure(query_seq, protein)
        key = (value, -abs(len(protein) - len(query_seq or '')))
        if best is None or key > best[0]:
            best = (key, parts, protein, value)
    if best is None:
        return None
    return best[1], best[2], best[3]


def _is_home_gene(score_frac, identity):
    return identity >= HOME_GENE_MIN_IDENTITY and score_frac >= HOME_GENE_MIN_SCORE_FRAC


def _search_gff_by_name(gff_file, candidate_names, query_seq=None, genome_file=None):
    """
    Search GFF for genes matching any of the candidate names.
    
    Scoring system:
    - Exact Name match: 100 points
    - Exact ID match: 90 points
    - Name/description contains full candidate name (>4 chars): 50 points
    - Partial word match (>5 chars): 20 points
    
    Collects ALL matches. With ``query_seq`` and ``genome_file`` (the pipeline always
    passes them) the matches are tried from the best name score down, and the first gene
    whose protein agrees with the query is returned, with the transcript that agrees
    best (see ``_is_home_gene``); a gene without CDS can only pass when the translation
    of its span does. Without them the highest-scoring gene is returned unchecked.
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

    def _result(gene_id, cds_parts):
        info = gene_features[gene_id]
        return {
            'gene_id': gene_id,
            'gene_name': info.get('gene_name', gene_id),
            'chrom': info['chrom'],
            'start': info['start'],
            'end': info['end'],
            'strand': info['strand'],
            'cds_parts': sorted(cds_parts, key=lambda x: x['start'])
        }

    def _span_as_cds(gene_id):
        info = gene_features[gene_id]
        return [{
            'chrom': info['chrom'],
            'start': info['start'] - 1,
            'end': info['end'],
            'strand': info['strand'],
            'phase': 0
        }]

    # Best name score first; equal scores in the order of the file.
    file_order = {gid: n for n, gid in enumerate(matched_genes)}
    ranked = sorted(matched_genes, key=lambda g: (-matched_genes[g], file_order[g]))

    if not query_seq or not genome_file:
        # No query to compare with: the highest-scoring gene, as given by its name.
        best_gene = ranked[0]
        print(f"[GFF Match] Best match: {best_gene} (score={matched_genes[best_gene]}), "
              f"not compared with the query")
        transcripts = _gene_transcripts(best_gene, transcript_to_gene, cds_by_parent)
        if transcripts:
            return _result(best_gene, transcripts[0][1])
        print(f"[GFF Match] Gene {best_gene} found but no CDS features. Using gene span.")
        return _result(best_gene, _span_as_cds(best_gene))

    genome_seqs = load_genome_subset(
        genome_file,
        {gene_features[g]['chrom'] for g in ranked[:HOME_GENE_MAX_NAME_CANDIDATES]})
    if len(ranked) > HOME_GENE_MAX_NAME_CANDIDATES:
        print(f"[GFF Match] {len(ranked)} genes match a name; comparing the "
              f"{HOME_GENE_MAX_NAME_CANDIDATES} with the best name score with the query.")
    for gene_id in ranked[:HOME_GENE_MAX_NAME_CANDIDATES]:
        info = gene_features[gene_id]
        chrom_seq = genome_seqs.get(info['chrom'])
        if not chrom_seq:
            continue
        chosen = _best_transcript(gene_id, info['strand'], chrom_seq,
                                  transcript_to_gene, cds_by_parent, query_seq)
        if chosen is None:
            # A gene row without CDS rows: its span is the coding sequence only if the
            # translation of the span is the query.
            span = _span_as_cds(gene_id)
            protein = _cds_protein(span, chrom_seq, info['strand'])
            chosen = (span, protein, _query_match(query_seq, protein))
        cds_parts, protein, (score_frac, identity) = chosen
        if _is_home_gene(score_frac, identity):
            print(f"[GFF Match] Best match: {gene_id} (score={matched_genes[gene_id]}); "
                  f"its protein ({len(protein)} aa) agrees with the query: "
                  f"{identity:.1f}% identity, score fraction {score_frac:.2f}")
            return _result(gene_id, cds_parts)
        print(f"[GFF Match]   Rejected {gene_id} (name score={matched_genes[gene_id]}): "
              f"its protein does not agree with the query "
              f"({identity:.1f}% identity, score fraction {score_frac:.2f})")

    print("[GFF Match] No gene found by name agrees with the query sequence.")
    return None


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
        chrom = gene_info['chrom']
        if chrom not in genome_seqs:
            continue

        # Every transcript of the gene is compared with the query and the best one is
        # kept: the first transcript in the file can be a short isoform.
        chosen = _best_transcript(
            gene_id, gene_info['strand'], genome_seqs[chrom],
            transcript_to_gene, cds_by_parent, query_seq,
            measure=lambda q, p: (_quick_protein_identity(q, p),))
        if chosen is None:
            continue
        cds_parts, protein, (identity,) = chosen
        cds_parts_sorted = sorted(cds_parts, key=lambda x: x['start'])
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


def compute_genomic_footprint(exons):
    """Home GOI genomic footprint from annotated exons.

    Returns ``(genomic_span_bp, max_intron_bp)`` — the gene's start-to-end span
    and its largest intron (gap between consecutive exons). Used downstream to
    derive a data-driven max-intron for the GOI fallback hit-chain. Reads each
    exon's genomic interval from ``gstart``/``gend`` or ``coords`` (mirrors the
    BED-writing logic). Returns ``(None, None)`` when no exon carries coords.
    """
    intervals = []
    for e in exons:
        if 'gstart' in e and 'gend' in e:
            s, t = int(e['gstart']), int(e['gend'])
        elif e.get('coords'):
            s, t = int(e['coords'][0]), int(e['coords'][1])
        else:
            continue
        intervals.append((min(s, t), max(s, t)))
    if not intervals:
        return None, None
    intervals.sort()
    genomic_span_bp = intervals[-1][1] - intervals[0][0]
    max_intron_bp = 0
    for (_s0, e0), (s1, _e1) in zip(intervals, intervals[1:]):
        max_intron_bp = max(max_intron_bp, s1 - e0)
    return genomic_span_bp, max_intron_bp


def extract_exons_from_gff_match(match, genome_file, query_seq=None):
    """
    Extract individual CDS/exon protein sequences from a GFF match.

    When query_seq is provided, computes per-exon percent identity (pident)
    by aligning the extracted exon protein against the corresponding region
    of the query. Also detects start/stop codons and splice sites.

    Returns tuple: (exons_list, full_protein)
    """
    genome_seqs = load_genome(genome_file)
    chrom = match['chrom']

    if chrom not in genome_seqs:
        print(f"[GFF Exons] Chromosome {chrom} not found in genome!")
        return []

    chrom_seq = genome_seqs[chrom]
    exons = []
    cds_parts = sorted(match['cds_parts'], key=lambda x: x['start'])

    # If on minus strand, exons are numbered from the 3' end of DNA (= 5' of mRNA)
    # For minus strand: last CDS part (highest coord) is exon 1
    if match['strand'] == '-':
        cds_parts = list(reversed(cds_parts))

    gene_name = match.get('gene_name', match['gene_id'])
    n_exons = len(cds_parts)

    # First pass: translate all exons and build full protein to compute identity
    exon_proteins = []
    exon_dna_seqs = []
    for i, part in enumerate(cds_parts, start=1):
        exon_dna = chrom_seq[part['start']:part['end']]

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

        exon_dna_seqs.append(exon_dna)
        exon_proteins.append(translate(exon_dna) if len(exon_dna) >= 3 else '')

    # Build full concatenated protein
    # GFF phase = bases to skip at the 5' coding end (first CDS in coding direction).
    cds_parts_for_full = sorted(match['cds_parts'], key=lambda x: x['start'])
    coding_first_phase = int(cds_parts_for_full[0].get('phase', 0) or 0) \
        if match['strand'] == '+' \
        else int(cds_parts_for_full[-1].get('phase', 0) or 0)
    all_dna = ""
    for part in cds_parts_for_full:
        all_dna += chrom_seq[part['start']:part['end']]
    if match['strand'] == '-':
        all_dna = reverse_complement(all_dna)
    if coding_first_phase > 0:
        all_dna = all_dna[coding_first_phase:]
    all_dna = all_dna[:len(all_dna) - len(all_dna) % 3]
    full_protein = translate(all_dna).split('*')[0]

    # Compute per-exon identity against query protein
    # Map each exon to its approximate query region by cumulative AA offset
    cumulative_aa = 0
    for i, (part, exon_prot, exon_dna) in enumerate(
            zip(cds_parts, exon_proteins, exon_dna_seqs), start=1):
        clean_prot = exon_prot.replace('*', '')
        if len(clean_prot) < 3:
            cumulative_aa += len(clean_prot)
            continue

        # Compute pident: find best matching region in query via sliding window
        pident = 0.0
        best_qstart = cumulative_aa
        best_qend = cumulative_aa + len(clean_prot)
        if query_seq and len(clean_prot) >= 3:
            pident, best_qstart, best_qend = _best_match_pident(
                clean_prot, query_seq, hint_start=cumulative_aa)

        # Update cumulative offset using best match position for next exon
        # This self-corrects any drift from split codons at exon boundaries
        if best_qend > 0:
            cumulative_aa = best_qend
        else:
            cumulative_aa += len(clean_prot)

        # Detect splice sites from genomic DNA
        gs = part['start']
        ge = part['end']
        splice_donor = None
        splice_acceptor = None

        if i > 1:  # Not the first exon
            if match['strand'] == '+':
                if gs >= 2:
                    splice_acceptor = chrom_seq[gs - 2:gs].upper()
            else:
                if ge + 2 <= len(chrom_seq):
                    splice_acceptor = reverse_complement(
                        chrom_seq[ge:ge + 2]).upper()

        if i < n_exons:  # Not the last exon
            if match['strand'] == '+':
                if ge + 2 <= len(chrom_seq):
                    splice_donor = chrom_seq[ge:ge + 2].upper()
            else:
                if gs >= 2:
                    splice_donor = reverse_complement(
                        chrom_seq[gs - 2:gs]).upper()

        # Detect start/stop codons
        has_start_codon = (i == 1 and exon_prot.startswith('M'))
        has_stop_codon = False
        if i == n_exons and len(exon_dna) >= 3:
            last_codon = exon_dna[-3:].upper()
            # Check if original (before stop removal) protein ended with stop
            has_stop_codon = last_codon in stop_codons()

        exons.append({
            'id': f"GOI_{gene_name}|exon_{i}",
            'seq': clean_prot,
            'coords': (gs, ge),
            'exon_num': i,
            'total_exons': n_exons,
            'strand': match['strand'],
            'chrom': chrom,
            'gstart': gs,
            'gend': ge,
            'qstart': best_qstart,
            'qend': best_qend,
            'pident': pident,
            'has_start_codon': has_start_codon,
            'has_stop_codon': has_stop_codon,
            'splice_donor': splice_donor,
            'splice_acceptor': splice_acceptor,
            'method': 'gff_annotation'
        })

    return exons, full_protein


def _compute_pident(seq1, seq2):
    """
    Compute percent identity between two protein sequences.
    Uses a simple positional comparison (global alignment style).
    For equal-length or near-equal-length sequences from exon extraction,
    this gives an accurate pident.
    """
    if not seq1 or not seq2:
        return 0.0

    # If lengths differ significantly, use the shorter for comparison
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))

    if min_len == 0:
        return 0.0

    # Simple positional identity (no gaps)
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)

    # Penalize length difference
    pident = (matches / max_len) * 100.0
    return round(pident, 1)


def _best_match_pident(exon_prot, query_seq, hint_start=0, search_margin=20):
    """
    Find the best-matching window in query_seq for exon_prot using a sliding
    window approach.  Returns (pident, best_qstart, best_qend).

    hint_start is the expected start position (from cumulative AA counting).
    search_margin is how many positions to search around the hint to handle
    boundary offsets from split codons.
    """
    exon_len = len(exon_prot)
    query_len = len(query_seq)

    if exon_len == 0 or query_len == 0:
        return 0.0, hint_start, hint_start

    best_pident = 0.0
    best_start = hint_start
    best_end = min(hint_start + exon_len, query_len)

    # Define search window around hint_start
    search_start = max(0, hint_start - search_margin)
    search_end = min(query_len - exon_len + 1, hint_start + search_margin + 1)

    # Ensure we have at least some range to search
    if search_end <= search_start:
        search_end = search_start + 1

    for start in range(search_start, min(search_end, query_len)):
        end = min(start + exon_len, query_len)
        window = query_seq[start:end]
        compare_len = min(len(window), exon_len)
        if compare_len == 0:
            continue
        matches = sum(1 for a, b in zip(exon_prot[:compare_len],
                                         window[:compare_len]) if a == b)
        denom = max(exon_len, len(window))
        pid = (matches / denom) * 100.0
        if pid > best_pident:
            best_pident = pid
            best_start = start
            best_end = end

    return round(best_pident, 1), best_start, best_end


# =============================================================================
# TANDEM DUPLICATION DETECTION
# =============================================================================


def _calc_overlap(start1, end1, start2, end2):
    """Calculate the overlap length between two 1-based closed intervals."""
    return max(0, min(end1, end2) - max(start1, start2) + 1)


def _deduplicate_hits(hits, overlap_threshold=0.8):
    """
    Remove redundant hits that overlap heavily in genomic space.

    Keeps the hit with the higher identity (pident) when two hits overlap
    by >= *overlap_threshold* of the shorter hit's length.

    Parameters
    ----------
    hits : list of dict
        Each dict must have 'gstart', 'gend', 'chrom', and 'pident'.
    overlap_threshold : float
        Minimum fraction of the shorter hit covered by overlap to consider
        two hits redundant.

    Returns
    -------
    list of dict – deduplicated hits, sorted by genomic start.
    """
    if len(hits) <= 1:
        return list(hits)
    # Sort by descending identity so the best hit is kept first
    sorted_hits = sorted(hits, key=lambda h: -h.get('pident', 0))
    kept = []
    for h in sorted_hits:
        is_dup = False
        for k in kept:
            if h.get('chrom') != k.get('chrom'):
                continue
            ov = _calc_overlap(h['gstart'], h['gend'], k['gstart'], k['gend'])
            shorter = min(h['gend'] - h['gstart'] + 1,
                          k['gend'] - k['gstart'] + 1)
            if shorter > 0 and ov / shorter >= overlap_threshold:
                is_dup = True
                break
        if not is_dup:
            kept.append(h)
    kept.sort(key=lambda h: (h.get('chrom', ''), h['gstart']))
    return kept


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

def _check_miniprot():
    """Check if miniprot is available."""
    try:
        result = subprocess.run(["miniprot", "--version"],
                                capture_output=True, text=True, timeout=5)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False

MINIPROT_AVAILABLE = _check_miniprot()
# Queries shorter than this are aligned with miniprot -L 10 (ORFs from 10 aa, default 30).
# The search sets it per run (iterative_search_runner --miniprot_short_orf_max_len; the
# short_query preset module raises it to 100); the home gene model always uses 60.
MINIPROT_SHORT_ORF_MAX_LEN = 60
_MINIPROT_T_OK = None


def _miniprot_genetic_code_args():
    """['-T', n] for miniprot when the process's genetic code is not the standard one (an
    index-building option: miniprot builds its index from the target on every call here).
    -T arrived in miniprot 0.13-r248; an older miniprot would refuse the option, so that is
    an error, not a silent standard-code run."""
    global _MINIPROT_T_OK
    code = get_genetic_code()
    if code == 1:
        return []
    if _MINIPROT_T_OK is None:
        try:
            out = subprocess.run(["miniprot", "--version"], capture_output=True, text=True,
                                 timeout=5).stdout.strip()
            m = re.match(r"(\d+)\.(\d+)", out)
            _MINIPROT_T_OK = bool(m) and (int(m.group(1)), int(m.group(2))) >= (0, 13)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            _MINIPROT_T_OK = False
    if not _MINIPROT_T_OK:
        raise RuntimeError(f"genetic code {code} needs miniprot >= 0.13 (option -T)")
    return ["-T", str(code)]


def _parse_miniprot_gff(gff_text):
    """
    Parse miniprot GFF3 output to extract gene models.

    Returns list of gene models, each with:
      - rank, identity, score
      - list of CDS features with (start, end, strand, phase, target_info)
      - 'protein': the translated model from miniprot --trans (the ##STA line that
        precedes the model's mRNA record), or '' when --trans was not used
    """
    models = {}  # ID -> model dict
    pending_sta = ''

    for line in gff_text.strip().split('\n'):
        if line.startswith('##STA'):
            # miniprot --trans writes ##PAF, then ##STA<TAB>protein, then that
            # alignment's mRNA/CDS records.
            parts = line.split('\t')
            pending_sta = parts[1].strip() if len(parts) > 1 else ''
            continue
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
                'cds_list': [],
                'protein': pending_sta,
            }
            pending_sta = ''
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


# =============================================================================
# TERMINAL REFINEMENT OF MINIPROT MODELS
# =============================================================================
# miniprot reports only the part of the query it aligns, so a model normally starts a few
# codons inside the gene (without the initiator Met) and ends before the stop codon. On the
# 2026-09 family benchmark this was the most common disagreement with curated annotations:
# every apyrase model was short by exactly 3 codons at the N terminus and 2 codons plus the
# stop at the C terminus, while all ten internal splice junctions were byte-identical to the
# curated gene. Walking in frame outwards from the terminal CDS recovers those ends.

TERMINUS_MAX_START_CODONS = 30
TERMINUS_MAX_STOP_CODONS = 20
_TERMINUS_UNKNOWN_QUERY_BOUND = 10


def _terminus_bound(unaligned_query_residues, hard_cap):
    """How many codons the walk may cross.

    The unaligned query residues at that terminus are the natural bound: they are what the
    model is missing. With none unaligned (the query itself may be truncated) a small fixed
    budget is allowed instead, so the walk can complete a start/stop codon but cannot invent
    a long terminus.
    """
    if unaligned_query_residues > 0:
        return max(1, min(hard_cap, unaligned_query_residues + 3))
    return min(hard_cap, _TERMINUS_UNKNOWN_QUERY_BOUND)


def refine_model_termini(cds_ordered, chrom_seq, strand,
                         max_start_codons=TERMINUS_MAX_START_CODONS,
                         max_stop_codons=TERMINUS_MAX_STOP_CODONS):
    """Extend a miniprot model's terminal CDS to the gene's own start and stop codon.

    ``cds_ordered`` is the model's CDS list in CODING order (the caller reverses it on the
    minus strand), with ``gstart`` 0-based inclusive and ``gend`` exclusive in ``chrom_seq``
    coordinates. It is modified in place.

    The walk is deliberately timid, because extending a terminus is inventing sequence the
    aligner did not support:
      * it refuses to start when the model's terminal codon is already ATG / a stop codon;
      * it refuses to cross a splice site (``AG`` immediately 5' of the first coding base, or
        ``GT`` immediately 3' of the last one) — there the model is missing a neighbouring
        EXON, and walking into the intron would fabricate residues;
      * the upstream walk stops at the first in-frame stop codon rather than crossing it;
      * it only extends when it finds an ATG (upstream) or a stop codon (downstream) inside
        the bound, so a terminus it cannot resolve is left exactly as miniprot reported it;
      * a first CDS with a non-zero phase, or a last CDS whose coding length is not a whole
        number of codons, is left alone (the model does not begin or end on a codon border).

    The splice-site guard costs some legitimate extensions, since ``GT``/``AG`` also occur
    inside coding sequence: of four real Apis genes checked on 2026-09-15 (query trimmed at
    both ends), three were recovered to the curated boundary exactly and in the fourth a
    chance ``GT`` in a valine codon blocked a 5-codon walk to the real stop. Refusing to
    invent residues is the safer error.

    Returns ``(prefix_aa, suffix_aa)``: the residues the extension adds to miniprot's own
    translation. The stop codon itself is not translated, so ``suffix_aa`` holds only the
    coding residues in front of it.
    """
    stops = stop_codons()       # of the process's genetic code (set_genetic_code)
    if not cds_ordered or not chrom_seq:
        return '', ''

    first, last = cds_ordered[0], cds_ordered[-1]
    prefix = suffix = ''
    n = len(chrom_seq)

    # --- start codon: walk upstream, in frame, from the first coding base
    if int(first.get('phase', 0) or 0) == 0:
        if strand == '+':
            edge = int(first['gstart'])
            at_splice = edge >= 2 and chrom_seq[edge - 2:edge].upper() == 'AG'
            already = chrom_seq[edge:edge + 3].upper() == 'ATG'
            steps = min(max_start_codons, edge // 3)
            if not at_splice and not already:
                for k in range(1, steps + 1):
                    s = edge - 3 * k
                    codon = chrom_seq[s:s + 3].upper()
                    if codon in stops:
                        break
                    if codon == 'ATG':
                        first['gstart'] = s
                        prefix = translate(chrom_seq[s:edge]).replace('*', '')
                        break
        else:
            edge = int(first['gend'])
            at_splice = (edge + 2 <= n and
                         reverse_complement(chrom_seq[edge:edge + 2]).upper() == 'AG')
            already = (edge >= 3 and
                       reverse_complement(chrom_seq[edge - 3:edge]).upper() == 'ATG')
            steps = min(max_start_codons, (n - edge) // 3)
            if not at_splice and not already:
                for k in range(1, steps + 1):
                    e = edge + 3 * k
                    codon = reverse_complement(chrom_seq[e - 3:e]).upper()
                    if codon in stops:
                        break
                    if codon == 'ATG':
                        first['gend'] = e
                        prefix = translate(
                            reverse_complement(chrom_seq[edge:e])).replace('*', '')
                        break

    # --- stop codon: walk downstream, in frame, from the last coding base
    coding_len = (int(last['gend']) - int(last['gstart'])) - int(last.get('phase', 0) or 0)
    if coding_len > 0 and coding_len % 3 == 0:
        if strand == '+':
            edge = int(last['gend'])
            at_splice = edge + 2 <= n and chrom_seq[edge:edge + 2].upper() == 'GT'
            already = (edge >= 3 and chrom_seq[edge - 3:edge].upper() in stops)
            steps = min(max_stop_codons, (n - edge) // 3)
            if not at_splice and not already:
                for k in range(1, steps + 1):
                    e = edge + 3 * k
                    if chrom_seq[e - 3:e].upper() in stops:
                        last['gend'] = e
                        suffix = translate(chrom_seq[edge:e]).replace('*', '')
                        break
        else:
            edge = int(last['gstart'])
            at_splice = (edge >= 2 and
                         reverse_complement(chrom_seq[edge - 2:edge]).upper() == 'GT')
            already = (edge + 3 <= n and
                       reverse_complement(chrom_seq[edge:edge + 3]).upper() in stops)
            steps = min(max_stop_codons, edge // 3)
            if not at_splice and not already:
                for k in range(1, steps + 1):
                    s = edge - 3 * k
                    if reverse_complement(chrom_seq[s:s + 3]).upper() in stops:
                        last['gstart'] = s
                        suffix = translate(
                            reverse_complement(chrom_seq[s:edge])).replace('*', '')
                        break

    return prefix, suffix


def _run_miniprot(query_seq, chrom_seq, chrom_name, max_intron=200000, sensitive=False):
    """Align a protein to a stretch of DNA with miniprot; the parsed models, best first.

    An empty list when miniprot is missing, fails, times out or finds nothing.
    """
    if not query_seq or not chrom_seq or not MINIPROT_AVAILABLE:
        if query_seq and chrom_seq:
            print("[annotate_goi_exons] WARNING: miniprot not available, returning raw query",
                  file=sys.stderr)
        return []

    # mkstemp, not a pid-derived name: pids repeat across containers/PID namespaces
    # that share /tmp, and a collision would silently swap one task's model for another's.
    _qfd, query_file = tempfile.mkstemp(prefix="synvoy_mp_query_", suffix=".faa")
    _tfd, target_file = tempfile.mkstemp(prefix="synvoy_mp_target_", suffix=".fna")
    os.close(_qfd)
    os.close(_tfd)

    try:
        write_fasta([("query", query_seq)], query_file)
        write_fasta([(chrom_name, chrom_seq)], target_file)

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

        if query_len < MINIPROT_SHORT_ORF_MAX_LEN:
            # Short peptide mode (e.g., toxins)
            cmd.extend([
                "-L", "10",             # Min ORF length 10aa (default 30)
            ])
            if query_len < 30:
                cmd.extend([
                    "-S",               # No splicing for very short peptides
                    "-G", "2000",       # Override max intron to small value
                ])

        cmd.extend(_miniprot_genetic_code_args())
        cmd.extend([target_file, query_file])

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        if result.returncode != 0:
            print(f"[miniprot] Warning: non-zero exit ({result.returncode}): "
                  f"{result.stderr[:200]}", file=sys.stderr)
            return []
        if not result.stdout.strip():
            print("[miniprot] No alignments found", file=sys.stderr)
            return []

        models = _parse_miniprot_gff(result.stdout)
        if not models:
            print("[miniprot] No gene models parsed from GFF", file=sys.stderr)
        return models

    except subprocess.TimeoutExpired:
        print("[miniprot] Timed out after 120s", file=sys.stderr)
        return []
    except Exception as e:
        print(f"[miniprot] Error: {e}", file=sys.stderr)
        return []
    finally:
        for f in [query_file, target_file]:
            if os.path.exists(f):
                os.remove(f)


def _build_model(best, query_seq, chrom_seq, chrom_name, region_offset=0):
    """Exon list and protein of ONE parsed miniprot model: ``(exons, full_protein)``.

    ``exons`` are dicts with id, seq, exon_num, qstart, qend, gstart, gend, strand, chrom,
    pident, has_start_codon, has_stop_codon, splice_donor, splice_acceptor; coordinates
    are those of ``chrom_seq`` plus ``region_offset``. A model without CDS gives ([], '').
    """
    if not best['cds_list']:
        return [], ''
    
    # Extract CDS DNA and build exon list
    exons = []
    model_strand = best['strand']

    cds_ordered = best['cds_list']
    if model_strand == '-':
        cds_ordered = list(reversed(cds_ordered))

    # miniprot reports only the aligned span: recover the model's own start and stop
    # codon where the sequence supports them (see refine_model_termini). The codon count
    # BEFORE the walk is what miniprot's own translation covers.
    cds_codons_before = sum(int(c['gend']) - int(c['gstart'])
                            for c in cds_ordered) // 3
    term_prefix, term_suffix = refine_model_termini(
        cds_ordered, chrom_seq, model_strand,
        max_start_codons=_terminus_bound(
            max(0, int(cds_ordered[0].get('qstart', 0) or 0) - 1),
            TERMINUS_MAX_START_CODONS),
        max_stop_codons=_terminus_bound(
            max(0, len(query_seq) - int(cds_ordered[-1].get('qend', 0) or 0)),
            TERMINUS_MAX_STOP_CODONS),
    )
    if term_prefix or term_suffix:
        print(f"[miniprot] Terminal refinement: +{len(term_prefix)} aa at the N "
              f"terminus, +{len(term_suffix)} aa + stop at the C terminus",
              file=sys.stderr)

    for i, cds in enumerate(cds_ordered, 1):
        gs = cds['gstart']
        ge = cds['gend']

        # Extract DNA
        exon_dna = chrom_seq[gs:ge]
        if model_strand == '-':
            exon_dna = reverse_complement(exon_dna)

        # Every CDS is read in ITS OWN frame: GFF3 phase = bases to skip at the
        # coding-direction 5' end to reach the first complete codon (after the
        # reverse-complement above, that end is on the left). Trimming only the
        # first CDS read every exon behind a phase-1/2 intron in the wrong frame:
        # a 64 %-identity Harpegnathos hyaluronidase model came out at 47.6 %, and
        # every multi-exon GOI model of the 2026-09-11 family benchmark was affected.
        # The residue of a codon split by an intron is dropped here; the model
        # protein below comes from miniprot's own translation, which keeps it.
        phase = cds['phase']
        coding_dna = exon_dna[phase:] if phase > 0 else exon_dna

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
        
        exon_prot = translate(coding_dna)

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
            'phase': phase,
            'chrom': chrom_name,
            'pident': best['identity'] * 100,  # miniprot reports 0-1
            'has_start_codon': (i == 1 and exon_prot.startswith('M')),
            'has_stop_codon': (i == len(cds_ordered) and
                               exon_dna[-3:].upper() in stop_codons()),
            'splice_donor': splice_donor,
            'splice_acceptor': splice_acceptor,
            'method': 'miniprot'
        })

    # Model protein: miniprot's own translation (--trans, the ##STA record). It is
    # the only construction that stays in frame across a frameshift miniprot has
    # modelled (Frameshift=1) and keeps intron-split codons; concatenating the CDS
    # DNA runs out of frame from the first frameshift to the end of the model.
    # Fallback, without --trans output: the per-exon in-frame translations.
    full_protein = (best.get('protein') or '').replace('*', '')
    if full_protein:
        # miniprot's translation covers the aligned span only; the terminal refinement
        # above widened the terminal CDS, so those residues are added here. The per-exon
        # fallback below translates the widened CDS and already contains them.
        # An alignment that ends mid-codon leaves miniprot's translation one residue
        # longer than its own CDS, and that residue is the first one the C-terminal walk
        # re-translates — drop the overlap instead of duplicating it.
        overlap = max(0, len(full_protein) - cds_codons_before)
        full_protein = term_prefix + full_protein + term_suffix[overlap:]
    else:
        full_protein = ''.join(e['seq'] for e in exons)
    print(f"[miniprot] Model rank {best.get('rank', 1)}: {len(exons)} exon(s), "
          f"{len(full_protein)} aa, identity={best['identity']:.2f}, "
          f"score={best['score']:.0f}", file=sys.stderr)
    
    return exons, full_protein


def annotate_using_miniprot(query_seq, chrom_seq, chrom_name,
                             strand=None, max_intron=200000,
                             sensitive=False, region_offset=0):
    """
    Use miniprot to align query protein to genomic region: the BEST gene model.

    Args:
        query_seq: Query protein sequence (string)
        chrom_seq: Target genomic DNA sequence (string)
        chrom_name: Chromosome/scaffold name
        strand: Unused (miniprot searches both strands)
        max_intron: Max intron size (-G flag)
        sensitive: If True, use ultra-sensitive params (for GOI)
        region_offset: Offset to add to coordinates (if chrom_seq is a sub-region)

    Returns:
        (exon_list, full_protein); ([], query_seq) when there is no model.
    """
    if not query_seq or not chrom_seq:
        return [], query_seq or ''
    models = _run_miniprot(query_seq, chrom_seq, chrom_name, max_intron, sensitive)
    if not models:
        return [], query_seq
    try:
        exons, full_protein = _build_model(models[0], query_seq, chrom_seq, chrom_name, region_offset)
    except Exception as e:
        print(f"[miniprot] Error: {e}", file=sys.stderr)
        return [], query_seq
    if not exons:
        return [], query_seq
    return exons, full_protein or query_seq


# A further alignment is a model of its own when it scores at least this share of the best
# one. It is the limit miniprot applies itself (--outs: 0.3 in the sensitive call used
# for the gene of interest, 0.99 otherwise) and is applied again to the masked rounds
# of miniprot_gene_models.
MINIPROT_OUTS_SENSITIVE = 0.3
MINIPROT_OUTS_DEFAULT = 0.99


def miniprot_gene_models(query_seq, chrom_seq, chrom_name, max_intron=200000,
                         sensitive=False, region_offset=0, max_models=10, evidence=None):
    """Every gene model miniprot finds in the region that does not lie on a better one.

    miniprot chains the query to each place it aligns and reports the alignments in order
    of score. In a tandem array these are the copies of the gene. A model is kept when
    less than half of its coding sequence overlaps the coding sequence of a model kept
    before it, on either strand, so a second alignment to the same gene is dropped and a
    neighbouring copy is not.

    One run does not report every copy: with a long maximum intron (here the length of
    the region) miniprot leaves out a neighbouring copy that is more diverged than the
    best one (measured: of two copies 2 kb apart, the one at 80 % identity is reported
    with a maximum intron up to 6 kb and not from 10 kb on). ``evidence`` closes that
    gap: intervals ``(start, end)`` in the coordinates of ``chrom_seq`` where the query
    is known to align (the hits of the window search). As long as one of them lies on no
    coding exon of a kept model, the coding exons of the kept models are masked and
    miniprot runs again on what is left. A locus whose hits are all explained by the
    first model costs one run, as before.

    Returns ``[(exons, full_protein), ...]``, best first, at most ``max_models``.
    """
    if not query_seq or not chrom_seq:
        return []
    limit = max(1, int(max_models))
    min_ratio = MINIPROT_OUTS_SENSITIVE if sensitive else MINIPROT_OUTS_DEFAULT
    kept, kept_cds = [], []
    best_score = None
    target = chrom_seq

    def unexplained():
        for start, end in evidence or ():
            if not any(min(end, ke) - max(start, kb) > 0 for cds in kept_cds for kb, ke in cds):
                return True
        return False

    for _round in range(limit):
        added = 0
        for model in _run_miniprot(query_seq, target, chrom_name, max_intron, sensitive):
            if len(kept) >= limit:
                break
            score = float(model.get('score', 0) or 0)
            if best_score is None:
                best_score = score
            elif score < min_ratio * best_score:
                continue
            cds = [(c['gstart'], c['gend']) for c in model.get('cds_list', [])]
            coding = sum(e - b for b, e in cds)
            if coding <= 0:
                continue
            shared = sum(max(0, min(e, ke) - max(b, kb))
                         for b, e in cds for other in kept_cds for kb, ke in other)
            if shared > 0.5 * coding:
                continue
            try:
                # Always built on the unmasked sequence: the ends of a model are refined
                # on the DNA next to it.
                exons, protein = _build_model(model, query_seq, chrom_seq, chrom_name, region_offset)
            except Exception as e:
                print(f"[miniprot] Error: {e}", file=sys.stderr)
                continue
            if not exons:
                continue
            kept.append((exons, protein or ''))
            kept_cds.append(cds)
            added += 1
        if not added or len(kept) >= limit or not unexplained():
            break
        masked = list(chrom_seq)
        for cds in kept_cds:
            for b, e in cds:
                masked[b:e] = 'N' * (e - b)
        target = ''.join(masked)
    return kept


def _hit_region(hits, chrom_seq, gap_search_window):
    """The stretch of DNA a model may use: the hits plus ``gap_search_window`` on each
    side. Returns ``(start, end)`` or None without usable hits."""
    hit_starts = [h.get('gstart', h.get('start', 0)) for h in hits]
    hit_ends = [h.get('gend', h.get('end', 0)) for h in hits]
    if not hit_starts or not hit_ends:
        return None
    return (max(0, min(hit_starts) - gap_search_window),
            min(len(chrom_seq), max(hit_ends) + gap_search_window))


def annotate_exons_from_hit_list(hits, query_seq, chrom_seq, chrom_name,
                                  gap_search_window=50000, sensitive=False):
    """
    The best gene model for a set of hits: one miniprot run over the hits plus
    ``gap_search_window`` on each side. The hits only say where to look.

    Returns ``(exons, full_protein)``; ``([], query_seq)`` when there is no model.
    """
    if not hits or not query_seq or not chrom_seq:
        return [], query_seq or ''
    region = _hit_region(hits, chrom_seq, gap_search_window)
    if region is None:
        return [], query_seq
    region_start, region_end = region
    is_goi = sensitive or any('GOI' in str(h.get('query', '')) for h in hits)
    return annotate_using_miniprot(
        query_seq, chrom_seq[region_start:region_end], chrom_name,
        max_intron=max(20000, region_end - region_start),
        sensitive=is_goi,
        region_offset=region_start
    )


def annotate_models_from_hit_list(hits, query_seq, chrom_seq, chrom_name,
                                   gap_search_window=50000, sensitive=False, max_models=10):
    """
    Every gene model for a set of hits (see miniprot_gene_models), best first:
    ``[(exons, full_protein), ...]``, empty when there is none. The hits are also the
    evidence for further copies: while one of them lies on no model, miniprot is run
    again on the region with the models found so far masked.
    """
    if not hits or not query_seq or not chrom_seq:
        return []
    region = _hit_region(hits, chrom_seq, gap_search_window)
    if region is None:
        return []
    region_start, region_end = region
    is_goi = sensitive or any('GOI' in str(h.get('query', '')) for h in hits)
    evidence = []
    for h in hits:
        h_start = h.get('gstart', h.get('start'))
        h_end = h.get('gend', h.get('end'))
        if h_start is None or h_end is None:
            continue
        lo, hi = min(h_start, h_end) - region_start, max(h_start, h_end) - region_start
        if hi > 0 and lo < region_end - region_start:
            evidence.append((max(0, lo), hi))
    return miniprot_gene_models(
        query_seq, chrom_seq[region_start:region_end], chrom_name,
        max_intron=max(20000, region_end - region_start),
        sensitive=is_goi,
        region_offset=region_start,
        max_models=max_models,
        evidence=evidence if int(max_models) > 1 else None,
    )


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
    parser.add_argument("--gene_name", default="",
                        help="Name of the home gene as it stands in the annotation, for a "
                             "query whose own name is not found there. It is tried first; "
                             "like every name match it must agree with the query sequence")
    parser.add_argument("--gap_search_window", type=int, default=50000,
                        help="DNA added on each side of the hits before miniprot builds the "
                             "home gene model (only when no annotated gene matches)")
    parser.add_argument("--output_exons", required=True,
                        help="Output: exon protein FASTA")
    parser.add_argument("--output_bed", required=True,
                        help="Output: exon locations BED")
    parser.add_argument("--output_info", required=True,
                        help="Output: annotation info JSON")

    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the home genome (default 1, the standard code)")
    args = parser.parse_args()
    set_genetic_code(args.genetic_code)

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
            args.gff_search_window,
            gene_name=args.gene_name.strip(),
        )

        if match:
            result = extract_exons_from_gff_match(match, args.genome,
                                                  query_seq=query_seq)
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

    # Enforce a consistent prefix across the full sequence and all exon sequences.
    base_goi_name = query_id
    if exons and exons[0].get('id', '').startswith('GOI_'):
        import re
        m = re.match(r'^GOI_([A-Za-z0-9_.-]+)(?:\||$)', exons[0]['id'])
        if m:
            base_goi_name = m.group(1)
            
    # Fallback if base_goi_name is empty (e.g. Pro Mode with no query_id)
    if not base_goi_name:
        base_goi_name = "target"
        
    goi_full_id = f"GOI_{base_goi_name}"
    fasta_records.append((goi_full_id, full_protein))

    # Add individual exon sequences
    for idx, exon in enumerate(exons, start=1):
        raw_exon_id = exon.get('id', '') or f"exon_{idx}"
        exon_num = exon.get('exon_num', idx)
        
        if raw_exon_id.startswith('GOI_copy_'):
            # Tandem copies maintain their own ID
            exon_id = raw_exon_id
        else:
            # Overwrite all exon IDs to strictly match the full sequence's parent ID
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

    # Home GOI genomic footprint. Downstream (iterative_search_runner) scales the
    # home GOI's OWN largest intron by a margin to bound the GOI fallback hit-chain,
    # so two far-apart spurious hits can't be stitched into one gene — instead of a
    # hand-tuned --max_intron.
    genomic_span_bp, max_intron_bp = compute_genomic_footprint(exons)

    # Write info JSON
    info = {
        'query_id': query_id,
        'query_length': len(query_seq),
        'method': method_used,
        'genomic_span_bp': genomic_span_bp,
        'max_intron_bp': max_intron_bp,
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
