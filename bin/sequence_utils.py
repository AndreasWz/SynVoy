#!/usr/bin/env python3
"""
sequence_utils.py - Robust sequence handling without external dependencies

Philosophy: Assume nothing. Work with what's there. If nothing is there, number sequentially.

Handles:
- FASTA parsing (no BioPython)
- GFF parsing (no dependencies)
- ID extraction from any format (NCBI, Ensembl, custom, unknown)
- Fallback to sequential numbering
"""

import argparse
import gzip
import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union
from urllib.parse import unquote


# =============================================================================
# FASTA PARSING (No BioPython)
# =============================================================================

def parse_fasta(filepath: str, user_header: bool = False) -> Iterator[Tuple[str, str, str]]:
    """
    Parse FASTA file yielding (raw_header, id, sequence) tuples.
    
    Works with any FASTA format. No assumptions about header structure. The id is the header's
    first word, verbatim (first_word): for a file SynVoy wrote that is the id it wrote and the id
    the aligners report, so two records never read as one. Until 2026-10-31 every reader took
    extract_id, which read FBgn0000001|exon_2 and FBgn0000001|exon_1 as FBgn0000001 (fix-phase
    item 25).

    user_header: the file comes from a user, so the id goes through extract_id, which knows the
    forms of public databases (P12345|LY6E_HUMAN and sp|P12345|LY6E_HUMAN read P12345).
    """
    read_id = extract_id if user_header else first_word
    current_header = None
    current_seq = []
    
    _open = gzip.open if str(filepath).endswith('.gz') else open
    
    with _open(filepath, 'rt') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if current_header is not None:
                    seq = ''.join(current_seq)
                    yield current_header, read_id(current_header), seq
                current_header = line[1:]  # Remove '>'
                current_seq = []
            elif line and current_header is not None:
                # Skip empty lines, handle sequence
                current_seq.append(line.strip())
        
        # Don't forget the last sequence
        if current_header is not None:
            seq = ''.join(current_seq)
            yield current_header, read_id(current_header), seq


def parse_fasta_to_dict(filepath: str) -> Dict[str, str]:
    """Parse FASTA to dict {id: sequence}."""
    return {clean_id: seq for _, clean_id, seq in parse_fasta(filepath)}


def write_fasta(records: List[Tuple[str, str]], filepath: str, wrap: int = 80):
    """
    Write FASTA file.
    
    Args:
        records: List of (header, sequence) tuples
        filepath: Output path
        wrap: Line wrap width (0 = no wrap)
    """
    with open(filepath, 'w') as f:
        for header, seq in records:
            # Ensure header doesn't start with '>'
            header = header.lstrip('>')
            f.write(f'>{header}\n')
            
            if wrap > 0:
                for i in range(0, len(seq), wrap):
                    f.write(seq[i:i+wrap] + '\n')
            else:
                f.write(seq + '\n')


def count_sequences(filepath: str) -> int:
    """Count sequences in FASTA file efficiently."""
    count = 0
    _open = gzip.open if str(filepath).endswith('.gz') else open
    with _open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


# =============================================================================
# ID EXTRACTION - The heart of robust parsing
# =============================================================================

# Known ID patterns (order matters - more specific first)
ID_PATTERNS = [
    # The pipeline's own GOI records: GOI_<name>, GOI_<name>|exon_1, GOI_copy_1|... The id is
    # the whole first word, whatever <name> starts with. Since item 25 no pattern below can cut
    # one and first_word gives the same id: a guard, in front of any pattern that may come later.
    (r'^(GOI_\S+)', 'goi'),

    # NCBI RefSeq proteins: XP_123456.1, NP_123456.1, YP_123456.1, WP_123456.1
    (r'^([XNYWZ]P_\d+(?:\.\d+)?)(?=$|[\s|])', 'ncbi_refseq'),
    
    # NCBI GenBank proteins: AAA12345.1
    # Require token boundary after accession so long IDs like LOC143834063 are
    # not incorrectly truncated to LOC14383.
    (r'^([A-Z]{3}\d{5}(?:\.\d+)?)(?=$|[\s|])', 'ncbi_genbank'),
    
    # NCBI GI (legacy): gi|12345|
    (r'^gi\|(\d+)\|', 'ncbi_gi'),
    
    # NCBI accession with version: NC_123456.1, NZ_CP123456.1
    (r'^((?:NC|NZ|NW|NT|AC|NG)_[A-Z]*\d+(?:\.\d+)?)(?=$|[\s|])', 'ncbi_nucl'),
    
    # NCBI's "FASTA Nucleotide CDS" and "FASTA Protein" downloads of a nucleotide record:
    # lcl|NC_123456.1_cds_XP_123456.1_1 [gene=...] and lcl|NC_123456.1_prot_XP_123456.1_1 -> the
    # protein, a RefSeq or a GenBank one (KAF1234567.1). The pattern stood here since c6bd07a with
    # a lookahead that can never hold, so such a header kept its whole first word (fix item 28).
    (r'^lcl\|\S*?_(?:cds|prot)_((?:[XNYWZ]P_\d+|[A-Z]{3}(?:\d{5}|\d{7}))(?:\.\d+)?)_\d+(?=$|\s)', 'ncbi_cds'),

    # A local id, lcl|X -> X; NCBI's database tags, ref|X|, gb|X|LOCUS, emb|X|, dbj|X| -> X; ENA's
    # ENA|X|X.1 -> X.1. The query id becomes the GOI's name, and SynVoy reads a GOI's name up to its
    # first '|', so such an id read as its tag ('ref', 'ENA') (freeze fix e2); an id that still
    # holds a '|' stops normalize_query.py.
    (r'^lcl\|([^|\s]+)(?=$|\s)', 'local_id'),
    (r'^(?:ref|gb|emb|dbj)\|([^|\s]+)(?:\|\S*)?(?=$|\s)', 'ncbi_tagged'),
    (r'^ENA\|[^|\s]+\|([^|\s]+)(?=$|\s)', 'ena'),
    
    # UniProt: P12345, Q9NXB0, A0A1B2C3D4
    (r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?=$|[\s|])', 'uniprot'),
    
    # UniProt with prefix: sp|P12345|NAME, tr|A0A123|NAME
    (r'^(?:sp|tr)\|([^\|]+)\|', 'uniprot_full'),
    
    # Ensembl proteins: ENSP00000123456, ENSMUSP00000123456
    (r'^(ENS[A-Z]*P\d+(?:\.\d+)?)(?=$|[\s|])', 'ensembl_protein'),
    
    # Ensembl genes: ENSG00000123456
    (r'^(ENS[A-Z]*G\d+(?:\.\d+)?)(?=$|[\s|])', 'ensembl_gene'),
    
    # Ensembl transcripts: ENST00000123456
    (r'^(ENS[A-Z]*T\d+(?:\.\d+)?)(?=$|[\s|])', 'ensembl_transcript'),
    
    # FlyBase: FBgn0000001, FBpp0000001
    (r'^(FB[a-z]{2}\d+)(?=$|[\s|])', 'flybase'),
    
    # WormBase: WBGene00000001
    (r'^(WBGene\d+)(?=$|[\s|])', 'wormbase'),
    
    # TAIR (Arabidopsis): AT1G01010.1
    (r'^(AT[1-5MC]G\d+(?:\.\d+)?)(?=$|[\s|])', 'tair'),
    
    # Prodigal format: >1_1 # start # end # strand # ID=1_1;...
    (r'^(\d+_\d+)\s*#', 'prodigal'),
    
    # No generic patterns here (fix-phase item 25). ^([A-Z]{2,6}_\d+), ^(gene_\d+) and
    # ^(prot(?:ein)?_\d+) stood here without a boundary after the digits, so they cut the first
    # word: NODE_1_length_5000_cov_1.0 read NODE_1, SUPER_10_unloc_1 read SUPER_10, and
    # FUN_000001|exon_2 read FUN_000001. With a boundary of whitespace or the end they would
    # return what first_word returns, so leaving them out is the same rule in less code.
    
    # Just take first word if nothing else matches
    (r'^(\S+)', 'first_word'),
]


_FIRST_WORD = re.compile(r'\S+', re.ASCII)


def first_word(header: str) -> str:
    """A FASTA record's name: its header up to the first whitespace, verbatim. samtools, miniprot
    and tblastn call a sequence by this word, and SynVoy reads MMseqs2's tables by it too (the
    header columns, mmseqs_first_words). ASCII whitespace, the C library's: samtools, miniprot and
    MMseqs2 keep a no-break space in the word. "unknown" for an empty header, as extract_id."""
    m = _FIRST_WORD.search(header or '')
    return m.group(0) if m else "unknown"


# Every MMseqs2 table SynVoy parses: the 12 BLAST-tabular columns, then the query's and the target's
# FASTA header. MMseqs2's own query and target columns rewrite some names (measured on MMseqs2
# 18-8cc5c, fix-phase item 25): a word that starts with a database tag by the accession inside
# (lcl|X, ref|X|, sp|X|Y, gnl|db|X, gi|1|ref|X| as X; gnl|a|1 and gnl|b|1 both as 1), a leading uc
# or consensus_ dropped (ucon_5 as on_5, consensus_24 as 24, which can be another record's name).
# The headers hold each record's name as its FASTA does, so the tables are read by their first word
# (fix item 29; until then such names stopped the run, since 2026-11-02). The target header comes
# last: a tab inside it stays part of it. A query header must hold none: SynVoy writes the queries.
MMSEQS_FORMAT_OUTPUT = ("query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,"
                        "qheader,theader")


def mmseqs_first_words(line: str) -> List[str]:
    """One row of an MMSEQS_FORMAT_OUTPUT table as the 12 BLAST-tabular columns, the query and the
    target named by the first word of their FASTA headers; StepFailed on a row without the headers,
    and on a row whose own names are no part of those words. MMseqs2 names a record by a part of
    its header's first word (a field of a tagged word, or the word without a leading uc or
    consensus_; an empty name included), so such a row is misaligned (a tab inside the query's
    header would make the target the first word of the query's description), not renamed."""
    p = line.rstrip('\n').split('\t')
    if len(p) < 14:
        raise StepFailed(f"an MMseqs2 table row has {len(p)} columns, not the 12 BLAST-tabular ones and "
                         f"the two headers (--format-output {MMSEQS_FORMAT_OUTPUT}): {line[:200]!r}")
    query, target = first_word(p[12]), first_word('\t'.join(p[13:]))
    if p[0] not in query or p[1] not in target:
        raise StepFailed(f"an MMseqs2 table row names its query {p[0]!r} and its target {p[1]!r}, but its "
                         f"header columns read {query!r} and {target!r}: the row is misaligned (a tab inside "
                         f"the query's header?): {line[:200]!r}")
    return [query, target] + p[2:12]


def mmseqs_table_first_words(raw_path: str, out_path: str) -> int:
    """Rewrite an MMSEQS_FORMAT_OUTPUT table into the 12-column table every reader of the pipeline
    parses, rows in their order, columns 3-12 verbatim: the number of rows."""
    n = 0
    with open(raw_path, errors='surrogateescape') as src, \
            open(out_path, 'w', errors='surrogateescape') as out:
        for line in src:
            if line.strip():
                out.write('\t'.join(mmseqs_first_words(line)) + '\n')
                n += 1
    return n


def id_problem(name: str, whitespace: bool = True, pipe: bool = False, gff: bool = False) -> Optional[str]:
    """Why SynVoy cannot use ``name`` as a sequence's or gene's id, or None (fix-phase d4: the
    names that lost results in silence). ``whitespace``: a FASTA record of it is read by its first
    word, so ``ORF 1`` and ``ORF 2`` are one id (item 32). ``pipe``: SynVoy reads a gene id up to
    its first '|' (a predicted home gene's id holds its scaffold's name, item 30; a flanking gene
    id, fix item 29), so two such ids can read as one. ``gff``: the id is written into GFF column
    9, where ',' makes a Parent a list, ';' ends the attribute, '=' and '&' are reserved and '%'
    is decoded (fix item 28)."""
    if gff:
        bad = sorted(set(re.findall(r"[,;=&%]", name)))
        if bad:
            return ("holds " + " and ".join(f"'{c}'" for c in bad) + ", reserved in GFF3: SynVoy writes "
                    "the GOI's name into the ID, Parent and SynVoy_Parent of its models, where ',' makes a "
                    "list, ';' ends the attribute and '%' is decoded")
    if whitespace and re.search(r"\s", name):
        return ("holds whitespace: a FASTA record of it is read by its first word, so two such "
                "ids can read as one")
    if pipe and "|" in name:
        return ("holds '|': SynVoy reads a gene id up to its first '|' (the ids of the home genes "
                "predicted on a sequence hold the sequence's name), so two such ids can read as one")
    return None


def check_gene_ids(ids, source: str, written: str, pipe: bool = False) -> None:
    """Stop (StepFailed) on the first gene id of ``ids`` that id_problem rejects: ``source`` is the
    file the id comes from (the GFF), ``written`` the file SynVoy would write it to. ``pipe``: also
    a '|' (a flanking gene's id is read up to it)."""
    for rid in ids:
        problem = id_problem(rid, pipe=pipe)
        if problem:
            raise StepFailed(f"{source}: the gene id {rid!r} {problem} (SynVoy writes it to {written}). "
                             "Rename the gene's ID in that file: a word without whitespace"
                             + (" or '|'." if pipe else "."))


def extract_id(header: str, fallback_counter: Optional[int] = None) -> str:
    """
    Extract a clean, usable ID from any FASTA header a user wrote (parse_fasta's user_header).
    A file SynVoy wrote is read by its first word (first_word); this function can read two of
    its records as one.
    
    Philosophy:
    1. Try known patterns (NCBI, Ensembl, UniProt, etc.)
    2. Fall back to first non-whitespace token
    3. If still nothing, use counter-based ID
    
    Args:
        header: FASTA header (without '>')
        fallback_counter: Optional counter for generating fallback IDs
        
    Returns:
        Clean ID string
    """
    if not header or not header.strip():
        if fallback_counter is not None:
            return f"seq_{fallback_counter}"
        return "unknown"
    
    header = header.strip()
    
    # Try each pattern
    for pattern, id_type in ID_PATTERNS:
        match = re.search(pattern, header)
        if match:
            extracted = match.group(1)
            # Clean up common issues
            extracted = extracted.strip('|').strip()
            if extracted:
                return extracted
    
    # Last resort: hash the header
    if fallback_counter is not None:
        return f"seq_{fallback_counter}"
    
    # Generate deterministic ID from header
    hash_id = hashlib.md5(header.encode()).hexdigest()[:8]
    return f"seq_{hash_id}"


# The ids fragment_query.generate_fragments writes: <query id>|frag_<start>_<end>_<total> and
# <query id>|slide_<start>_<end>_<total>; <start> and <end> are 1-based in the whole protein. An id
# that holds |frag_ or |slide_ in another form is a fragment too, without positions.
_FRAGMENT_ID_RE = re.compile(r"\|(frag|slide)_(?:(\d+)_(\d+)_(\d+)$)?")


def fragment_of(query_id: str):
    """(kind, start, end, total) of a fragment's id, kind 'frag' or 'slide' (the numbers None when
    the id holds none); None for any other id.

    The one reader of these ids: a hit's query positions start at the fragment's own residue 1, so
    a reader that knows only one kind reads the other as a full-length protein with positions
    near its start (fix-phase item 12)."""
    m = _FRAGMENT_ID_RE.search(query_id or "")
    if not m:
        return None
    numbers = (int(m.group(2)), int(m.group(3)), int(m.group(4))) if m.group(2) else (None, None, None)
    return (m.group(1), *numbers)


# =============================================================================
# GFF PARSING (No dependencies)
# =============================================================================

def parse_gff(filepath: str, feature_types: Optional[List[str]] = None) -> Iterator[Dict[str, Any]]:
    """
    Parse GFF/GFF3/GTF file yielding feature dictionaries.
    
    Works with any GFF-like format. Handles:
    - GFF3 (ID=...; Parent=...)
    - GTF (gene_id "..."; transcript_id "...")
    - Edge cases and malformed files
    
    Args:
        filepath: Path to GFF file
        feature_types: Optional list of feature types to include (e.g., ['CDS', 'gene'])
        
    Yields:
        Dict with keys: seqid, source, type, start, end, score, strand, phase, attributes
    """
    _open = gzip.open if str(filepath).endswith('.gz') else open
    with _open(filepath, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n\r')
            
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 8:
                continue  # Malformed line, skip silently
            
            feature_type = parts[2]
            
            # Filter by feature type if specified
            if feature_types and feature_type not in feature_types:
                continue
            
            # Parse coordinates
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue  # Invalid coordinates
            
            # Parse score
            score = None
            if parts[5] != '.':
                try:
                    score = float(parts[5])
                except ValueError:
                    pass
            
            # Parse phase
            phase = None
            if len(parts) > 7 and parts[7] != '.':
                try:
                    phase = int(parts[7])
                except ValueError:
                    pass
            
            # Parse attributes
            attributes = {}
            if len(parts) > 8:
                attributes = parse_gff_attributes(parts[8])
            
            yield {
                'seqid': parts[0],
                'source': parts[1],
                'type': feature_type,
                'start': start,
                'end': end,
                'score': score,
                'strand': parts[6] if parts[6] in ['+', '-'] else '.',
                'phase': phase,
                'attributes': attributes,
                'line_num': line_num
            }


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF3/GTF attribute string into dictionary.
    
    Handles both:
    - GFF3: ID=gene1;Name=ABC;Parent=transcript1
    - GTF: gene_id "gene1"; transcript_id "tx1";
    """
    attributes = {}
    
    if not attr_string or attr_string == '.':
        return attributes
    
    # Try GFF3 format first (key=value;)
    if '=' in attr_string:
        for pair in attr_string.split(';'):
            pair = pair.strip()
            if '=' in pair:
                key, value = pair.split('=', 1)
                # Full URL decode (covers %3B %3D %26 %2C and any other percent-escapes
                # that GFF3 producers may emit, e.g. %20 spaces in product names).
                attributes[key.strip()] = unquote(value.strip())
    
    # Try GTF format (key "value";)
    elif '"' in attr_string:
        for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
            key, value = match.groups()
            attributes[key] = value
    
    return attributes


def get_feature_id(feature: Dict[str, Any], fallback_counter: Optional[int] = None) -> str:
    """
    Get a usable ID from a GFF feature.
    
    Tries in order:
    1. ID attribute
    2. gene_id, transcript_id, protein_id
    3. Name attribute
    4. locus_tag
    5. Generate from position
    """
    attrs = feature.get('attributes', {})
    
    # Try standard ID attributes
    for key in ['ID', 'gene_id', 'transcript_id', 'protein_id', 'Name', 'locus_tag', 'gene']:
        if key in attrs and attrs[key]:
            return attrs[key]
    
    # Generate from position
    if fallback_counter is not None:
        return f"{feature['type']}_{fallback_counter}"
    
    seqid = feature.get('seqid', 'unknown')
    start = feature.get('start', 0)
    return f"{feature['type']}_{seqid}_{start}"


# =============================================================================
# SEQUENCE OPERATIONS
# =============================================================================

def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


# NCBI genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi): the
# amino acid of each of the 64 codons in TCAG order (TTT, TTC, TTA, TTG, TCT, ...), as NCBI
# prints them ('*' = stop). Codons that are a stop only in context (tables 27, 28, 31)
# carry their amino acid, as in NCBI's string.
NCBI_GENETIC_CODES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # standard
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",   # vertebrate mito
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # yeast mito
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # mold mito, Mycoplasma
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",   # invertebrate mito
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # ciliate nuclear
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",   # echinoderm mito
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # euplotid nuclear
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # bacterial
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # alternative yeast nuclear (CUG-Ser1)
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",  # ascidian mito
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # alt. flatworm mito
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Blepharisma
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # chlorophycean mito
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # trematode mito
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Scenedesmus mito
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Thraustochytrium mito
    24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Pterobranchia mito
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # candidate division SR1
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Pachysolen (CUG-Ala)
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # karyorelict nuclear
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Condylostoma nuclear
    29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Mesodinium nuclear
    30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # peritrich nuclear
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Blastocrithidia nuclear
    33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Cephalodiscidae mito
}
# The tables every engine of the search translates with: BLAST+ 2.17 (-db_gencode / -query
# _gencode: 1-6, 9-16, 21-31, 33), MMseqs2 18 (--translation-table: 1-6, 9-16, 21-31) and
# miniprot 0.18 (-T: 1-33, of which 7, 8, 17-20 do not exist). 33 is missing in MMseqs2.
SUPPORTED_GENETIC_CODES = frozenset(list(range(1, 7)) + list(range(9, 17)) + list(range(21, 32)))
_CODON_ORDER = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
_CODON_TABLES: Dict[int, Dict[str, str]] = {}
_GENETIC_CODE = 1     # the process default (set_genetic_code)


def check_genetic_code(table: Any) -> int:
    """The table as an int; ValueError when not every engine of the search supports it."""
    try:
        n = int(str(table).strip())
    except (TypeError, ValueError):
        raise ValueError(f"genetic code '{table}' is not a number") from None
    if n not in SUPPORTED_GENETIC_CODES:
        raise ValueError(
            f"genetic code {n} is not supported by every engine SynVoy uses (BLAST+, MMseqs2, "
            f"miniprot); use one of {', '.join(map(str, sorted(SUPPORTED_GENETIC_CODES)))}")
    return n


def set_genetic_code(table: Any) -> int:
    """Set the process default genetic code that translate() and stop_codons() use when
    no table is named. One per process: the search sets it per target genome at the start
    of each genome (one genome per worker task)."""
    global _GENETIC_CODE
    _GENETIC_CODE = check_genetic_code(table)
    return _GENETIC_CODE


def get_genetic_code() -> int:
    """The process default genetic code (1 unless set_genetic_code changed it)."""
    return _GENETIC_CODE


def codon_table(table: Optional[int] = None) -> Dict[str, str]:
    """{codon: amino acid} of an NCBI genetic code (the process default for None)."""
    n = _GENETIC_CODE if table is None else check_genetic_code(table)
    if n not in _CODON_TABLES:
        _CODON_TABLES[n] = dict(zip(_CODON_ORDER, NCBI_GENETIC_CODES[n]))
    return _CODON_TABLES[n]


def stop_codons(table: Optional[int] = None) -> frozenset:
    """The stop codons of a genetic code (the process default for None)."""
    return frozenset(c for c, aa in codon_table(table).items() if aa == '*')


def translate(seq: str, table: Optional[int] = None) -> str:
    """
    Translate DNA to protein with an NCBI genetic code.

    Args:
        seq: DNA sequence (U is read as T; a codon with N or another letter is X)
        table: NCBI table number; None = the process default (set_genetic_code, 1 unless
               set), so every caller that names no table follows the target's code
    """
    codons = codon_table(table)
    seq = seq.upper().replace('U', 'T')
    protein = []
    
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if 'N' in codon:
            protein.append('X')
        else:
            protein.append(codons.get(codon, 'X'))
    
    return ''.join(protein)


def extract_sequence(genome_dict: Dict[str, str], chrom: str, start: int, end: int, strand: str = '+') -> str:
    """
    Extract sequence from genome dictionary.
    
    Args:
        genome_dict: Dict of {chrom: sequence}
        chrom: Chromosome/contig name
        start: 1-based start position
        end: 1-based end position (inclusive)
        strand: '+' or '-'
    """
    if chrom not in genome_dict:
        return ""
    
    seq = genome_dict[chrom]
    # Convert to 0-based for Python slicing
    extracted = seq[start-1:end]
    
    if strand == '-':
        extracted = reverse_complement(extracted)
    
    return extracted


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def _genome_name(seen: Dict[str, str], header: str, filepath: str) -> str:
    """The sequence name of a genome's header, its first word: the name the aligners, samtools and
    the GFF use (fix-phase item 25; until 2026-10-31 extract_id's cut, so a block on
    NODE_1_length_5000_cov_1.0 found no sequence and SUPER_10_unloc_1 overwrote SUPER_10). Two
    records with one name stop the task: no tool can tell them apart, and samtools faidx keeps
    only one of them."""
    name = first_word(header)
    if name in seen:
        raise StepFailed(f"{filepath}: two records have the name {name}: >{seen[name]} and >{header}. "
                         "A genome's sequence names must be unique up to the first whitespace: "
                         "samtools faidx and a name lookup keep only one of them, so rename one.")
    seen[name] = header
    return name


def iter_genome(filepath: str) -> Iterator[Tuple[str, str]]:
    """(name, sequence) of each record of a genome FASTA, in file order; ``_genome_name``."""
    seen: Dict[str, str] = {}
    for header, _id, seq in parse_fasta(filepath):
        yield _genome_name(seen, header, filepath), seq


def genome_window(genome_path: str, chrom: str, start: int, end: int,
                  samtools_bin: str = "samtools") -> str:
    """genome[chrom:start..end] (1-based inclusive): ``samtools faidx`` when it can (a plain
    FASTA or a bgzip file), else read from the file (a plain-gzip genome, which samtools refuses,
    or no samtools). "" when the genome holds no such sequence. The rescue passes' one reader
    (fix-phase d5, item 35: the strong-synteny rescue skipped every block of a plain-gzip
    target)."""
    try:
        res = subprocess.run([samtools_bin, "faidx", genome_path, f"{chrom}:{start}-{end}"],
                             check=True, capture_output=True, text=True)
        return "".join(l for l in res.stdout.splitlines() if not l.startswith(">"))
    except (FileNotFoundError, subprocess.CalledProcessError):
        return load_genome_subset(genome_path, [chrom]).get(chrom, "")[max(0, start - 1):end]


def load_genome(filepath: str) -> Dict[str, str]:
    """Load genome FASTA as dict {chrom: sequence}."""
    return dict(iter_genome(filepath))


def load_genome_subset(filepath: str, names) -> Dict[str, str]:
    """Like ``load_genome``, but holds only the sequences named in ``names``: one chromosome of
    a 3 Gb genome, not the genome. The whole file is read, one record at a time, so that a name
    used twice further down still stops the task."""
    want = set(names)
    found: Dict[str, str] = {}
    if not want:
        return found
    for name, seq in iter_genome(filepath):
        if name in want:
            found[name] = seq
    return found


def get_sequence_lengths(filepath: str) -> Dict[str, int]:
    """Get lengths of all sequences in FASTA."""
    return {name: len(seq) for name, seq in iter_genome(filepath)}


_COPY_PARENT_RE = re.compile(r'_b\d+_[lr]\d+(?:m\d+)?_by_')


def strip_copy_parent(suffix: str) -> str:
    """A model id's part behind the first '|' without a seed-form tandem copy's parent:
    ``<genome>_b<k>_<l|r><j>_by_<parent>`` -> ``<genome>_b<k>_<l|r><j>``. The parent comes last
    and can hold a block-like token of its own (a GOI named Ly6e_b3_l1; fix item 26)."""
    m = _COPY_PARENT_RE.search(suffix or '')
    return suffix[:m.end() - len('_by_')] if m else (suffix or '')


def block_index_of(model_id: str, genome_tag: str) -> Optional[int]:
    """The search block a model of this genome was found in: N of
    ``<parent>|<genome_tag>_b<N>_...`` (``genome_tag`` = the genome name with '.', '-' and
    ' ' replaced by '_', as the model ids carry it). None for a model without a block
    number (cross-chromosome flanking models) or of another genome. Here, not in the
    search, so that the figures read it without importing the search."""
    rest = (model_id or '').split('|', 1)
    prefix = f"{genome_tag}_b"
    if len(rest) < 2 or not rest[1].startswith(prefix):
        return None
    digits = rest[1][len(prefix):].split('_', 1)[0]
    return int(digits) if digits.isdigit() else None


# =============================================================================
# BED PARSING
# =============================================================================

def parse_bed(path: Union[str, Path]) -> List[Dict[str, Any]]:
    """
    Parse a BED file into a list of dicts.

    Each row carries: chrom, start, end, name, strand, display_name.
    Missing optional columns get sensible defaults:
      - name: '{chrom}:{start}-{end}' if BED has fewer than 4 columns
      - strand: '+'
      - display_name: '' (BED col 7 used by extract_flanking_genes.py)

    Returns an empty list for missing, empty, or unreadable paths so callers
    can treat absence the same as "no rows".
    """
    rows: List[Dict[str, Any]] = []
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            chrom = parts[0]
            name = parts[3] if len(parts) > 3 and parts[3] else f"{chrom}:{start}-{end}"
            strand = parts[5] if len(parts) > 5 and parts[5] in {'+', '-'} else '+'
            display_name = parts[6] if len(parts) > 6 else ''
            rows.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'strand': strand,
                'display_name': display_name,
            })
    return rows


# =============================================================================
# FASTA HEADERS-ONLY (cheap iterator)
# =============================================================================

def iter_fasta_headers(path: Union[str, Path]) -> Iterator[str]:
    """
    Yield raw FASTA headers (without the leading '>') from a file.

    Cheap alternative to `parse_fasta` for scripts that only need IDs and
    never touch the sequence (e.g. counting, indexing, dedup).
    """
    _open = gzip.open if str(path).endswith('.gz') else open
    with _open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                yield line[1:].rstrip('\n\r')


# =============================================================================
# CLI / PROCESS HELPERS
# =============================================================================

_SPECIES_QUALIFIER_RE = re.compile(
    r"\s*\((?:strain|isolate|serotype|serovar|var\.|subsp\.|biotype|genotype)\b[^)]*\)",
    re.IGNORECASE)


def strip_species_qualifier(species_name):
    """Species name without a strain/isolate qualifier in parentheses.

    UniProt reports organisms with the qualifier attached: D6VTK4 (yeast STE2) comes
    back as ``Saccharomyces cerevisiae (strain ATCC 204508 / S288c)``. NCBI knows the
    species as ``Saccharomyces cerevisiae``, so the long form is rejected by
    ``datasets``, has no taxonomy entry, never equals the species name of an assembly
    record, and is noise in a figure label. Returns the input unchanged when there is
    nothing to strip, so a name whose real form contains other parentheses is kept.
    """
    if not species_name:
        return species_name
    out = re.sub(r"\s+", " ", _SPECIES_QUALIFIER_RE.sub("", species_name)).strip()
    return out or species_name


DEFAULT_TAXA_DB = os.path.expanduser("~/.etetoolkit/taxa.sqlite")


def offline_ncbi_taxa(dbfile: Optional[str] = None):
    """ete3's NCBITaxa on the NCBI taxonomy database already on disk (``~/.etetoolkit/
    taxa.sqlite``, the one synvoy_tree reads), never a download: None without ete3 or the
    database. ete3's web plugin imports the stdlib ``cgi`` module, which Python 3.13 removed,
    so a stub stands in for it during the import, as in synvoy_tree."""
    path = dbfile or DEFAULT_TAXA_DB
    if not os.path.exists(path):
        return None
    try:
        import sys as _sys
        import types as _types
        import warnings as _warnings
        if "cgi" not in _sys.modules:
            try:
                with _warnings.catch_warnings():
                    _warnings.simplefilter("ignore", DeprecationWarning)   # 3.11-3.12
                    import cgi  # noqa: F401
            except ModuleNotFoundError:
                _shim = _types.ModuleType("cgi")
                _shim.escape = lambda s, quote=True: s
                _sys.modules["cgi"] = _shim
        from ete3 import NCBITaxa  # type: ignore
        return NCBITaxa(dbfile=path)
    except Exception:
        return None


def str2bool(value: Any) -> bool:
    """
    Argparse-friendly boolean parser.

    Accepts: True/False, 'true'/'false', '1'/'0', 'yes'/'no', 'y'/'n', 't'/'f'
    (case-insensitive). Raises argparse.ArgumentTypeError on anything else,
    which gives a clean argparse error message.
    """
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    val = str(value).strip().lower()
    if val in {'true', '1', 'yes', 'y', 't'}:
        return True
    if val in {'false', '0', 'no', 'n', 'f'}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def run_command(
    cmd: Sequence[str],
    *,
    check: bool = True,
    capture: bool = False,
    text: bool = True,
    timeout: Optional[float] = None,
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Dict[str, str]] = None,
    input: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """
    Thin, consistent subprocess.run wrapper used across SynVoy scripts.

    Defaults match the common case: check=True, no capture (output streams to
    the task log), text mode. Pass capture=True when you need stdout/stderr
    on the returned CompletedProcess.

    Raises subprocess.CalledProcessError if check=True and the command fails;
    subprocess.TimeoutExpired if timeout elapses.
    """
    return subprocess.run(
        list(cmd),
        check=check,
        capture_output=capture,
        text=text,
        timeout=timeout,
        cwd=cwd,
        env=env,
        input=input,
    )


class StepFailed(BaseException):
    """A tool call that did not finish within its guard, or failed where no fallback may
    stand in for it (fix-phase item 18: the wall clock may decide whether a result exists,
    never what it contains). A BaseException on purpose: the search's ``except Exception``
    fallbacks (a window kept without its Smith-Waterman hits, a gene without its model) must
    not swallow it. The script that owns the task reports it and exits non-zero."""


class MemoryKilled(StepFailed):
    """A tool, or a search worker, killed for memory (fix-phase d6). The script exits
    MEMORY_KILL_EXIT, which Nextflow retries with more memory, as item 18 did for the
    Smith-Waterman child: a killed miniprot read as "no model", and a killed Augustus or search
    worker exited 1, which is not retried. The message says "ran out of memory"."""


MEMORY_KILL_EXIT = 137


def memory_killed(returncode: int, stderr: str = "") -> bool:
    """Whether a tool's exit means it was killed for memory: SIGKILL (-9 from subprocess, 137
    from a shell) or an allocation failure it reported itself."""
    text = stderr or ""
    return returncode in (-9, 137) or "MemoryError" in text or "std::bad_alloc" in text


def guard_failed(what: str, guard_seconds: float, param: Optional[str] = None) -> StepFailed:
    """The StepFailed of a tool that did not finish within its guard: ``what`` names the tool
    and the place, ``param`` the pipeline parameter that sets the guard."""
    return StepFailed(f"{what} did not finish within {guard_seconds:g} s"
                      + (f" (--{param})" if param else "")
                      + ". The task stops instead of going on without this result; on a slow "
                        "or busy machine raise the guard")


def sequence_not_in_genome(chrom: str, genome: str, names) -> StepFailed:
    """The StepFailed of a hit or search block on a sequence the loaded genome does not hold
    (fix-phase item 25). The hits come from this genome's own FASTA and every table is read by the
    first word of its headers (fix item 29), so it can only mean that a tool and the FASTA name the
    sequence differently. ``names``: the genome's sequence names."""
    names = sorted(names)
    shown = ", ".join(names[:3]) + (", ..." if len(names) > 3 else "")
    return StepFailed(
        f"[{genome}] a hit or search block lies on the sequence {chrom!r}, which this genome's FASTA "
        f"does not hold (its {len(names)} sequence names: {shown}). A tool and the FASTA name the "
        "sequence differently; every tool of the search names a sequence by the first word of its "
        "header. Please report this with the genome's first header lines.")


def run_guarded(cmd: Sequence[str], guard_seconds: float, what: str, param: Optional[str] = None,
                **popen_kwargs: Any) -> subprocess.CompletedProcess:
    """subprocess.run(capture_output=True, text=True) in a process group of its own, with a
    guard: when ``cmd`` has not finished after ``guard_seconds``, the whole group is killed
    (the tool's own children too, e.g. ssearch36 under smith_waterman_search.py) and
    StepFailed is raised, naming ``what`` (the tool and the place) and the pipeline parameter
    ``param`` that sets the guard. The return code is not checked."""
    import signal
    proc = subprocess.Popen(list(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            start_new_session=True, **popen_kwargs)
    try:
        out, err = proc.communicate(timeout=guard_seconds)
    except BaseException as exc:
        # The guard, and also an exception in the caller while it waits (Ctrl-C's
        # KeyboardInterrupt): the child has a session of its own, so a signal to the caller's
        # process group does not reach it. A SIGTERM or SIGKILL of the caller runs no Python
        # code; Nextflow then kills the task's whole process tree, and Slurm or a container
        # its cgroup, which hold the child too.
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        proc.wait()
        if isinstance(exc, subprocess.TimeoutExpired):
            raise guard_failed(what, guard_seconds, param) from None
        raise
    return subprocess.CompletedProcess(list(cmd), proc.returncode, out, err)


def require_binary(name: str, install_hint: Optional[str] = None) -> str:
    """
    Resolve an external tool on PATH or exit with a clear, consistent error.

    Returns the absolute path to the binary on success. On miss, prints a
    one-line error (plus install_hint if given) to stderr and raises
    FileNotFoundError so the calling script's main() can decide whether to
    sys.exit(1) or recover.
    """
    found = shutil.which(name)
    if found:
        return found
    msg = f"required binary not found on PATH: {name!r}"
    if install_hint:
        msg += f" — {install_hint}"
    print(f"ERROR: {msg}", file=sys.stderr)
    raise FileNotFoundError(msg)


# The exit code of a script whose StepFailed must end the run although its module ignores
# other failures (the rescue passes: errorStrategy maps this code to 'terminate').
STEP_FAILED_EXIT = 3


def exit_step_failed(main_fn) -> None:
    """Run a script's main() and exit with its return value; a StepFailed (a stop: a name, a guard,
    a failed tool) is printed and ends the script with STEP_FAILED_EXIT, without a traceback, a
    MemoryKilled with MEMORY_KILL_EXIT (retried). Exit 3 ends the run under the global
    errorStrategy, and run_step.py passes it through (the rescue passes)."""
    try:
        rc = main_fn()
    except MemoryKilled as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(MEMORY_KILL_EXIT)
    except StepFailed as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(STEP_FAILED_EXIT)
    sys.exit(rc)



def setup_logging(
    level: Union[int, str] = logging.INFO,
    *,
    name: Optional[str] = None,
    with_date: bool = True,
) -> logging.Logger:
    """
    One-call logging setup used across SynVoy bin/ scripts.

    Idempotent: if the root logger already has handlers (e.g. set up by a
    library or a parent script) we leave them alone and just adjust the
    requested logger's level. This avoids the duplicate-line spam that
    happens when basicConfig is called twice.
    """
    root = logging.getLogger()
    if not root.handlers:
        fmt = '%(asctime)s [%(levelname)s] %(message)s'
        kwargs: Dict[str, Any] = {'level': level, 'format': fmt}
        if with_date:
            kwargs['datefmt'] = '%Y-%m-%d %H:%M:%S'
        logging.basicConfig(**kwargs)
    logger = logging.getLogger(name) if name else root
    logger.setLevel(level)
    return logger


# =============================================================================
# JSON IO
# =============================================================================

def read_json(path: Union[str, Path]) -> Any:
    """Read a JSON file. Raises on missing or malformed input."""
    with open(path) as fh:
        return json.load(fh)


def write_json(
    path: Union[str, Path],
    obj: Any,
    *,
    indent: int = 2,
    atomic: bool = True,
) -> None:
    """
    Write `obj` as JSON to `path`.

    With atomic=True (default), the JSON is written to a sibling temp file
    and then os.replace'd into place — readers never see a half-written
    file if the writer is interrupted. Set atomic=False only when the
    destination is on a different filesystem than its parent dir.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if atomic:
        tmp = path.with_suffix(path.suffix + '.tmp')
        with open(tmp, 'w') as fh:
            json.dump(obj, fh, indent=indent)
        os.replace(tmp, path)
    else:
        with open(path, 'w') as fh:
            json.dump(obj, fh, indent=indent)


# =============================================================================
# ALIGNMENT (parasail Smith-Waterman)
# =============================================================================

def sw_align(query_seq: str, subject_seq: str,
             gap_open: int = 11, gap_extend: int = 1) -> Tuple[float, float]:
    """Smith-Waterman align two protein sequences, returning (score, percent_identity).

    The raw parasail SW score doubles as a bitscore proxy — we only ever use it to
    *rank* subjects against each other for a fixed query (paralog RBH / locus
    ownership), so no e-value calibration is needed. Identity is read off the
    traceback CIGAR when parasail exposes match/mismatch ops (=/X); otherwise it
    falls back to a score/length approximation.

    parasail is imported lazily so the many bin/ scripts that import sequence_utils
    but never align don't take a hard dependency on the bioconda parasail-python
    package. Raises RuntimeError with an actionable message when it is missing.
    """
    try:
        import parasail  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without parasail
        raise RuntimeError(
            "parasail-python is required for Smith-Waterman alignment "
            "(install via the synvoy_env conda env: `conda install -c bioconda parasail-python`)"
        ) from exc

    if not query_seq or not subject_seq:
        return 0.0, 0.0
    res = parasail.sw_trace_striped_16(query_seq, subject_seq, gap_open, gap_extend,
                                       parasail.blosum62)
    score = float(res.score)
    cigar = getattr(res, "cigar", None)
    cigar_ops = getattr(cigar, "seq", None) if cigar is not None else None
    if cigar_ops is None or len(cigar_ops) == 0:
        return score, 0.0
    matches = 0
    aligned = 0
    for op in cigar_ops:
        length = op >> 4
        code = op & 0xF
        # parasail cigar ops: 0=M (match-or-mismatch), 7='=', 8='X', 1=I, 2=D
        if code in (0, 7, 8):
            aligned += length
            if code == 7:
                matches += length
    if aligned == 0 or matches == 0:
        # blosum alphabet didn't tag =/X — approximate identity from the SW score
        # relative to the query length, capped at 100%.
        ident = min(1.0, score / max(1, len(query_seq))) * 100.0
    else:
        ident = 100.0 * matches / aligned
    return score, ident


# ---------------------------------------------------------------- home hits (S0)
# A hit continues the query after another when it starts at most this many residues before
# the other ends (the overlap one_copy_score tolerates between the hits of one copy).
QUERY_CONTINUATION_TOL = 10


def load_home_hits(paths, max_evalue, with_duplicates: bool = False):
    """Hits of the home search (12-column BLAST / MMseqs2 tables) with E <= ``max_evalue``,
    one per stretch of DNA: of hits that overlap by more than half of the shorter one the
    best-scoring stays, whatever their strands (both search programs report every exon).

    Dicts with ``chrom``, ``start`` (0-based), ``end``, ``strand`` ('+', or '-' when the
    subject coordinates run backwards), ``evalue``, ``bits``, ``qstart``, ``qend``. Moved from
    split_loci.load_raw_hits (batch of 2026-10-15, S0) with the same rule and order; strand
    and E-value are new keys, for the callers that chain hits. ``with_duplicates``: each kept
    hit also carries ``duplicates``, the hits it replaced (the first kept hit that a dropped
    one overlaps), so a caller can report the whole stretch they cover and their best E-value
    (the two programs' E-values are not on one scale: the best-scoring hit need not have the
    best E-value)."""
    hits = []
    for path in paths or []:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) < 12:
                    continue
                try:
                    a, b = int(p[8]), int(p[9])
                    evalue, bits = float(p[10]), float(p[11])
                except ValueError:
                    continue
                if evalue > max_evalue:
                    continue
                try:
                    qa, qb = int(p[6]), int(p[7])
                except ValueError:
                    continue
                hits.append({'chrom': p[1], 'start': min(a, b) - 1, 'end': max(a, b),
                             'strand': '+' if a <= b else '-', 'evalue': evalue, 'bits': bits,
                             'qstart': min(qa, qb), 'qend': max(qa, qb)})
    kept = []
    by_chrom = {}
    for h in sorted(hits, key=lambda x: (-x['bits'], x['chrom'], x['start'])):
        length = h['end'] - h['start']
        clash = None
        for k in by_chrom.get(h['chrom'], []):
            overlap = min(h['end'], k['end']) - max(h['start'], k['start'])
            if overlap > 0.5 * min(length, k['end'] - k['start']):
                clash = k
                break
        if clash is None:
            if with_duplicates:
                h['duplicates'] = []
            by_chrom.setdefault(h['chrom'], []).append(h)
            kept.append(h)
        elif with_duplicates:
            clash['duplicates'].append(h)
    return kept


def continues_query(a, b, tol: int = QUERY_CONTINUATION_TOL) -> bool:
    """Whether hit ``b`` continues the query after hit ``a``, the two taken in coding order
    (in_coding_order: on the minus strand the later one lies upstream): ``b`` starts at most
    ``tol`` residues before ``a`` ends, and it ends after ``a``. The end test refuses a hit
    that lies inside the end of ``a``, which would otherwise extend a chain with up to
    ``tol + 1`` residues already covered. ``a`` and ``b`` are dicts with ``qstart`` and
    ``qend`` (query positions, ascending on both strands)."""
    return b['qstart'] >= a['qend'] - tol and b['qend'] > a['qend']


def in_coding_order(hits, strand: str):
    """``hits`` (dicts with ``start``) in coding order: ascending genomic start on the plus
    strand, descending on the minus strand, so continues_query can read them pairwise."""
    return sorted(hits, key=lambda h: (h['start'], h.get('end', 0)), reverse=(strand == '-'))


def one_copy_score(hits, max_query_overlap: int = 10) -> float:
    """Best total bit score of hits that together cover the query at most once.

    The score of ONE copy of the gene at a locus: the heaviest set of hits whose query
    intervals do not overlap (by more than a few residues). ``hits`` are dicts with
    ``qstart``, ``qend`` and ``bits``. Summing all hits instead
    would add up the copies of a tandem array, or several weak relatives in a row, into
    one strong locus.
    """
    ordered = sorted(hits, key=lambda h: (h['qend'], h['qstart']))
    best = []
    for i, h in enumerate(ordered):
        prev = 0.0
        for j in range(i - 1, -1, -1):
            if ordered[j]['qend'] - max_query_overlap <= h['qstart']:
                prev = best[j]
                break
        best.append(max(best[i - 1] if i else 0.0, prev + h['bits']))
    return best[-1] if best else 0.0


def sw_align_span(query_seq: str, subject_seq: str) -> Dict[str, float]:
    """Local alignment of two proteins with the part of the QUERY it covers.

    Returns ``{'score', 'identity', 'q_start', 'q_end', 'query_cov'}``: identity in
    percent of the alignment columns between the first and the last aligned residue pair,
    gap columns included (as BLAST, MMseqs2 and miniprot count it; counting only the
    gap-free columns read several points higher on a gapped alignment and was compared
    with the same thresholds), ``q_start``/``q_end`` 1-based inclusive positions in
    the query, ``query_cov`` the share of the query between them (0-1). BLOSUM62 at BLAST
    gap costs 11/1. Raises RuntimeError without parasail, like ``sw_align``.
    """
    empty = {'score': 0.0, 'identity': 0.0, 'q_start': 0, 'q_end': 0, 'query_cov': 0.0}
    if not query_seq or not subject_seq:
        return empty
    try:
        import parasail  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without parasail
        raise RuntimeError("parasail-python is required for Smith-Waterman alignment") from exc
    res = parasail.sw_trace_striped_32(query_seq, subject_seq, 12, 1, parasail.blosum62)
    cigar = res.cigar
    # parasail writes the unaligned head of either sequence as leading I/D operations, so
    # the query position is tracked through the whole string and the aligned part is the
    # stretch from the first to the last column that holds a residue of both sequences.
    q_pos = int(cigar.beg_query)
    matches = aligned = 0
    columns = pending_gaps = 0         # gap columns count once a pair follows them
    q_first = q_last = None
    for op in cigar.seq:
        length, code = int(op >> 4), int(op & 0xF)
        if code in (0, 7, 8):          # M, =, X: columns with a residue on both sides
            if q_first is None:
                q_first = q_pos
            q_last = q_pos + length - 1
            aligned += length
            columns += pending_gaps + length
            pending_gaps = 0
            q_pos += length
            if code == 7:
                matches += length
        else:
            if q_first is not None:    # an inner gap; leading and trailing ones are not columns
                pending_gaps += length
            if code == 1:              # I: residues of the query only
                q_pos += length
    if not aligned:
        return empty
    return {
        'score': float(res.score),
        'identity': 100.0 * matches / columns,
        'q_start': q_first + 1,
        'q_end': q_last + 1,
        'query_cov': min(1.0, (q_last - q_first + 1) / len(query_seq)),
    }


# Karlin-Altschul parameters of gapped BLOSUM62 alignments at BLAST gap costs 11/1, the
# pair smith_waterman_search.py uses for the window search.
SW_LAMBDA_11_1 = 0.267
SW_K_11_1 = 0.041


def sw_bits(query_seq: str, subject_seq: str) -> float:
    """Bit score of the best local alignment of two proteins: BLOSUM62, BLAST gap costs
    11/1 (parasail 12/1), ``(lambda * S - ln K) / ln 2`` with lambda 0.267, K 0.041.

    A score in bits does not depend on the lengths or the database, so models of one gene
    built from different evidence (a spliced model, a chain of hits, a single exon) can be
    ranked by it. 0.0 for an empty sequence. Raises RuntimeError without parasail.
    """
    if not query_seq or not subject_seq:
        return 0.0
    try:
        import parasail  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without parasail
        raise RuntimeError("parasail-python is required for Smith-Waterman alignment") from exc
    score = parasail.sw_striped_32(query_seq, subject_seq, 12, 1, parasail.blosum62).score
    return (SW_LAMBDA_11_1 * score - math.log(SW_K_11_1)) / math.log(2)


def sw_pair_evalue(seq_a: str, seq_b: str) -> Tuple[float, float, float]:
    """E-value of the best local alignment of two proteins: ``(evalue, score, identity)``.

    The number of alignments with at least this score expected by chance between two
    unrelated sequences of these lengths (search space = len(a) x len(b), no edge
    correction). It answers "are these two proteins related" in a unit that does not
    depend on their lengths, which a percent identity or a common-subsequence share does.
    """
    if not seq_a or not seq_b:
        return float('inf'), 0.0, 0.0
    # parasail charges `open` for the first gap residue: 12/1 there is BLAST 11/1.
    score, ident = sw_align(seq_a, seq_b, gap_open=12, gap_extend=1)
    evalue = SW_K_11_1 * len(seq_a) * len(seq_b) * math.exp(-SW_LAMBDA_11_1 * score)
    return evalue, score, ident


# =============================================================================
# TESTING
# =============================================================================

def longest_collinear_run(ranks):
    """Longest monotonic (non-decreasing OR non-increasing) subsequence of home ranks.

    Single source of truth for "how collinear is this set of anchors". `None` entries
    (GOI proxies, unmapped genes) are ignored. Returns (length, direction) with '+' for
    increasing and '-' for decreasing; ties favour '+'.

    Both halves of the pipeline must agree on this. Before it lived here, the search
    core scored collinearity with a proper LIS while the region scorer used a fraction
    of adjacent monotonic PAIRS — two different definitions of one concept, so a block
    could be judged collinear when seeding and non-collinear when ranking. Adjacent-pair
    counting is also fragile: one transposed gene in the middle of an otherwise perfect
    run breaks two pairs, while an LIS absorbs it as a single skip.
    """
    vals = [r for r in ranks if r is not None]
    n = len(vals)
    if n <= 1:
        return n, '+'

    def _lis_non_decreasing(seq):
        m = len(seq)
        dp = [1] * m
        best = 1
        for i in range(m):
            for j in range(i):
                if seq[j] <= seq[i] and dp[j] + 1 > dp[i]:
                    dp[i] = dp[j] + 1
            if dp[i] > best:
                best = dp[i]
        return best

    inc = _lis_non_decreasing(vals)
    dec = _lis_non_decreasing(list(reversed(vals)))
    return (inc, '+') if inc >= dec else (dec, '-')


# =============================================================================
# PROTEIN SEQUENCE FEATURES (query_features.py)
# =============================================================================

# Kyte & Doolittle (1982) hydropathy scale. Residues outside the 20 (X, B, Z, U, *) are 0.
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4,
    'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
    'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}


def mean_hydropathy(seq: str) -> float:
    """Mean Kyte-Doolittle hydropathy of ``seq`` (0.0 for an empty one)."""
    if not seq:
        return 0.0
    return sum(KYTE_DOOLITTLE.get(c, 0.0) for c in seq.upper()) / len(seq)


def clopper_pearson(k: int, n: int, conf: float = 0.95) -> Tuple[float, float]:
    """Exact two-sided (Clopper-Pearson) interval of a proportion of ``k`` in ``n``, by
    bisection on the binomial distribution; (nan, nan) for n = 0."""
    if n <= 0:
        return float("nan"), float("nan")
    tail = (1.0 - conf) / 2.0

    def cdf(x: int, p: float) -> float:
        if p <= 0.0:
            return 1.0
        if p >= 1.0:
            return 0.0 if x < n else 1.0
        return sum(math.comb(n, i) * p ** i * (1.0 - p) ** (n - i) for i in range(x + 1))

    def solve(pred) -> float:
        lo, hi = 0.0, 1.0
        for _ in range(100):
            mid = (lo + hi) / 2.0
            if pred(mid):
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2.0

    lower = 0.0 if k == 0 else solve(lambda p: 1.0 - cdf(k - 1, p) >= tail)
    upper = 1.0 if k >= n else solve(lambda p: cdf(k, p) <= tail)
    return lower, upper


def shannon_entropy(seq: str) -> float:
    """Shannon entropy in bits of the residue composition of ``seq``."""
    if not seq:
        return 0.0
    counts: Dict[str, int] = {}
    for c in seq:
        counts[c] = counts.get(c, 0) + 1
    n = len(seq)
    return -sum((k / n) * math.log2(k / n) for k in counts.values())


def low_complexity_mask(seq: str, window: int = 12, max_entropy: float = 2.2) -> List[bool]:
    """True for every residue inside some window of ``window`` residues whose composition
    has a Shannon entropy of at most ``max_entropy`` bits. SEG's trigger window and
    complexity (K1) of BLAST's `-seg 12 2.2 2.5`, without SEG's extension and trimming, so
    an approximation of what SEG masks. A sequence shorter than the window is one window."""
    s = seq.upper()
    n = len(s)
    mask = [False] * n
    if n == 0:
        return mask
    w = min(window, n)
    for i in range(0, n - w + 1):
        if shannon_entropy(s[i:i + w]) <= max_entropy:
            for j in range(i, i + w):
                mask[j] = True
    return mask


def repeat_kmer_mask(seq: str, k: int = 4, min_offset: int = 4) -> List[bool]:
    """True for every residue covered by an exact ``k``-mer that occurs again in ``seq`` at
    an offset of at least ``min_offset`` (so a k-mer does not count as repeated by
    overlapping itself, e.g. AAAAA)."""
    s = seq.upper()
    n = len(s)
    mask = [False] * n
    where: Dict[str, List[int]] = {}
    for i in range(0, n - k + 1):
        where.setdefault(s[i:i + k], []).append(i)
    for starts in where.values():
        if len(starts) < 2:
            continue
        for i in starts:
            if any(abs(i - j) >= min_offset for j in starts):
                for p in range(i, i + k):
                    mask[p] = True
    return mask


# =============================================================================
# CODING MODELS FROM A GFF (gff_to_faa.py, home_context.py)
# =============================================================================

def cds_models_by_gene(gff_file: str, chroms: Optional[set] = None) -> Dict[str, Dict[str, Any]]:
    """The longest coding model of every gene of a GFF3: ``{gene_id: {'chrom', 'strand',
    'exons': [(start, end, phase), ...]}}`` (1-based CDS rows in file order).

    CDS rows are grouped by their Parent (the mRNA); an mRNA is assigned to its gene by
    walking up the Parent links until an ID starting with ``gene-``/``gene:`` (NCBI,
    Ensembl) or the top of the chain; per gene the mRNA with the most CDS bases wins (the
    first one on a tie). ``chroms`` restricts the models to those sequences. Missing file
    -> {}."""
    relationships: Dict[str, str] = {}
    mrnas: Dict[str, Dict[str, Any]] = {}
    if not os.path.exists(gff_file):
        return {}

    def attrs_of(text: str) -> Dict[str, str]:
        out = {}
        for item in text.strip().split(';'):
            if '=' in item:
                k, v = item.strip().split('=', 1)
                out[k] = v
        return out

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            attrs = attrs_of(parts[8])
            fid, parent = attrs.get('ID'), attrs.get('Parent')
            if fid and parent:
                relationships[fid] = parent
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            if chroms is not None and parts[0] not in chroms:
                continue
            mrna_id = attrs_of(parts[8]).get('Parent')
            if not mrna_id:
                continue
            if mrna_id not in mrnas:
                mrnas[mrna_id] = {'chrom': parts[0], 'strand': parts[6], 'exons': []}
            mrnas[mrna_id]['exons'].append((int(parts[3]), int(parts[4]),
                                            int(parts[7]) if parts[7] != '.' else 0))
    longest: Dict[str, Dict[str, Any]] = {}
    for mrna_id, info in mrnas.items():
        ancestor = mrna_id
        visited = set()
        while ancestor in relationships and ancestor not in visited:
            visited.add(ancestor)
            if ancestor.startswith('gene-') or ancestor.startswith('gene:'):
                break
            ancestor = relationships[ancestor]
        length = sum(e[1] - e[0] + 1 for e in info['exons'])
        if ancestor not in longest or length > longest[ancestor]['length']:
            longest[ancestor] = {'id': mrna_id, 'length': length, 'info': info}
    return {gid: d['info'] for gid, d in longest.items()}


def cds_model_protein(info: Dict[str, Any], chrom_seq: str) -> str:
    """Protein of a coding model (``cds_models_by_gene``) on its chromosome's sequence:
    CDS joined in coding order, the first coding CDS's phase trimmed, translated, cut at
    the first stop. '' when nothing is left."""
    exons = sorted(info['exons'], key=lambda x: x[0])
    coding_first_phase = exons[0][2] if info['strand'] == '+' else exons[-1][2]
    dna_seq = ''.join(chrom_seq[s - 1:e] for s, e, _p in exons)
    if info['strand'] == '-':
        dna_seq = reverse_complement(dna_seq)
    if coding_first_phase > 0:
        dna_seq = dna_seq[coding_first_phase:]
    prot = translate(dna_seq)
    return prot.split('*')[0] if '*' in prot else prot


def stream_fasta_lengths(filepath: str) -> Dict[str, int]:
    """{sequence name: length} of every record of a genome FASTA (plain or .gz), without holding
    a sequence in memory: a 3 Gb genome is read line by line. Names as ``_genome_name``."""
    lengths: Dict[str, int] = {}
    seen: Dict[str, str] = {}
    name = None
    _open = gzip.open if str(filepath).endswith('.gz') else open
    with _open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                name = _genome_name(seen, line[1:].rstrip('\n\r'), filepath)
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def sw_aligned_pairs(query_seq: str, subject_seq: str) -> List[Tuple[int, int]]:
    """0-based (query position, subject position) of every column of the best local
    alignment of two proteins that holds a residue of both sequences (BLOSUM62, BLAST gap
    costs 11/1). [] for an empty sequence; raises RuntimeError without parasail."""
    if not query_seq or not subject_seq:
        return []
    try:
        import parasail  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised only without parasail
        raise RuntimeError("parasail-python is required for Smith-Waterman alignment") from exc
    res = parasail.sw_trace_striped_32(query_seq, subject_seq, 12, 1, parasail.blosum62)
    cigar = res.cigar
    q, s = int(cigar.beg_query), int(cigar.beg_ref)
    pairs: List[Tuple[int, int]] = []
    for op in cigar.seq:
        length, code = int(op >> 4), int(op & 0xF)
        if code in (0, 7, 8):          # M, =, X: a residue of both
            pairs.extend((q + k, s + k) for k in range(length))
            q += length
            s += length
        elif code == 1:                # I: query only
            q += length
        elif code == 2:                # D: subject only
            s += length
    return pairs


# =============================================================================
# GENETIC CODE OF A TARGET GENOME (lib/GeneticCodes.groovy hands the option over)
# =============================================================================

def parse_genetic_code_spec(spec: Optional[str]) -> Dict[str, int]:
    """{genome: table} of a 'genome=table,genome=table' option ('' or None: {})."""
    out: Dict[str, int] = {}
    for item in (spec or "").split(","):
        if "=" not in item:
            continue
        name, table = item.rsplit("=", 1)
        if name.strip():
            out[name.strip()] = check_genetic_code(table)
    return out


def genetic_code_for(genome: str, spec: Dict[str, int], default: int = 1) -> int:
    """The table of a genome: named in ``spec`` by its file name (as given, e.g. with
    '.fna'), or by that name without '.gz' and its last extension; else ``default``."""
    name = os.path.basename(str(genome))
    if name in spec:
        return spec[name]
    stem = name[:-3] if name.endswith(".gz") else name
    stem = stem.rsplit(".", 1)[0] if "." in stem else stem
    return spec.get(stem, check_genetic_code(default))


def genetic_code_flag(tool: str, table: Optional[int] = None) -> List[str]:
    """The command-line option that gives an engine a genetic code, [] for the standard
    code (so a default run's commands do not change). ``tool``: 'mmseqs', 'tblastn',
    'miniprot', 'synvoy' (a SynVoy script's --genetic_code)."""
    n = get_genetic_code() if table is None else check_genetic_code(table)
    if n == 1:
        return []
    return {"mmseqs": ["--translation-table", str(n)],
            "tblastn": ["-db_gencode", str(n)],
            "miniprot": ["-T", str(n)],
            "synvoy": ["--genetic_code", str(n)]}[tool]
