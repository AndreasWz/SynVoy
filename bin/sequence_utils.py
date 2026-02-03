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

import re
import hashlib
from pathlib import Path
from typing import Iterator, Dict, List, Tuple, Optional, Any


# =============================================================================
# FASTA PARSING (No BioPython)
# =============================================================================

def parse_fasta(filepath: str) -> Iterator[Tuple[str, str, str]]:
    """
    Parse FASTA file yielding (raw_header, clean_id, sequence) tuples.
    
    Works with any FASTA format. No assumptions about header structure.
    """
    current_header = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if current_header is not None:
                    seq = ''.join(current_seq)
                    yield current_header, extract_id(current_header), seq
                current_header = line[1:]  # Remove '>'
                current_seq = []
            elif line and current_header is not None:
                # Skip empty lines, handle sequence
                current_seq.append(line.strip())
        
        # Don't forget the last sequence
        if current_header is not None:
            seq = ''.join(current_seq)
            yield current_header, extract_id(current_header), seq


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
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


# =============================================================================
# ID EXTRACTION - The heart of robust parsing
# =============================================================================

# Known ID patterns (order matters - more specific first)
ID_PATTERNS = [
    # NCBI RefSeq proteins: XP_123456.1, NP_123456.1, YP_123456.1, WP_123456.1
    (r'^([XNYWZ]P_\d+(?:\.\d+)?)', 'ncbi_refseq'),
    
    # NCBI GenBank proteins: AAA12345.1
    (r'^([A-Z]{3}\d{5}(?:\.\d+)?)', 'ncbi_genbank'),
    
    # NCBI GI (legacy): gi|12345|
    (r'^gi\|(\d+)\|', 'ncbi_gi'),
    
    # NCBI accession with version: NC_123456.1, NZ_CP123456.1
    (r'^((?:NC|NZ|NW|NT|AC|NG)_[A-Z]*\d+(?:\.\d+)?)', 'ncbi_nucl'),
    
    # NCBI assembly accession in sequence name: lcl|NC_123456.1_cds_XP_123456.1_1
    (r'cds_([XNYWZ]P_\d+(?:\.\d+)?)_', 'ncbi_cds'),
    
    # UniProt: P12345, Q9NXB0, A0A1B2C3D4
    (r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', 'uniprot'),
    
    # UniProt with prefix: sp|P12345|NAME, tr|A0A123|NAME
    (r'^(?:sp|tr)\|([^\|]+)\|', 'uniprot_full'),
    
    # Ensembl proteins: ENSP00000123456, ENSMUSP00000123456
    (r'^(ENS[A-Z]*P\d+(?:\.\d+)?)', 'ensembl_protein'),
    
    # Ensembl genes: ENSG00000123456
    (r'^(ENS[A-Z]*G\d+(?:\.\d+)?)', 'ensembl_gene'),
    
    # Ensembl transcripts: ENST00000123456
    (r'^(ENS[A-Z]*T\d+(?:\.\d+)?)', 'ensembl_transcript'),
    
    # FlyBase: FBgn0000001, FBpp0000001
    (r'^(FB[a-z]{2}\d+)', 'flybase'),
    
    # WormBase: WBGene00000001
    (r'^(WBGene\d+)', 'wormbase'),
    
    # TAIR (Arabidopsis): AT1G01010.1
    (r'^(AT[1-5MC]G\d+(?:\.\d+)?)', 'tair'),
    
    # Prodigal format: >1_1 # start # end # strand # ID=1_1;...
    (r'^(\d+_\d+)\s*#', 'prodigal'),
    
    # Generic locus_tag pattern: XXXX_12345
    (r'^([A-Z]{2,6}_\d+)', 'locus_tag'),
    
    # Generic gene ID with underscore: gene_12345
    (r'^(gene_\d+)', 'generic_gene'),
    
    # Generic protein ID: protein_12345, prot_12345
    (r'^(prot(?:ein)?_\d+)', 'generic_protein'),
    
    # Just take first word if nothing else matches
    (r'^(\S+)', 'first_word'),
]


def extract_id(header: str, fallback_counter: Optional[int] = None) -> str:
    """
    Extract a clean, usable ID from any FASTA header.
    
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


def extract_base_id(id_string: str) -> str:
    """
    Extract the base gene/protein ID, stripping variants, fragments, versions.
    
    Examples:
        "XP_123456.1" -> "XP_123456"
        "gene|var_1" -> "gene"
        "gene|frag_half_1" -> "gene"
        "ENSP00000123456.2" -> "ENSP00000123456"
        "lcl|NC_123_cds_XP_456.1_7" -> "XP_456"
    """
    if not id_string:
        return id_string
    
    original = id_string
    
    # Remove our own variant/fragment suffixes
    if '|' in id_string:
        parts = id_string.split('|')
        # Check if last part is a variant/fragment marker
        if len(parts) >= 2:
            last = parts[-1]
            if any(last.startswith(p) for p in ['var_', 'frag_', 'half_', 'third_', 'quarter_']):
                id_string = '|'.join(parts[:-1])
            elif re.match(r'^(var|frag|half|third|quarter)', last):
                id_string = '|'.join(parts[:-1])
    
    # Handle lcl|..._cds_XP_123.1_N format
    cds_match = re.search(r'cds_([XNYWZ]P_\d+)', id_string)
    if cds_match:
        id_string = cds_match.group(1)
    
    # Remove version numbers (.1, .2, etc.) from the end
    id_string = re.sub(r'\.\d+$', '', id_string)
    
    # Remove trailing underscore + number (like _1, _2 for CDS copies)
    id_string = re.sub(r'_\d+$', '', id_string)
    
    return id_string if id_string else original


def normalize_id(id_string: str) -> str:
    """
    Normalize ID for comparison (lowercase, no special chars).
    """
    if not id_string:
        return ""
    return re.sub(r'[^a-z0-9]', '', id_string.lower())


def ids_match(id1: str, id2: str) -> bool:
    """
    Check if two IDs refer to the same entity (fuzzy matching).
    """
    if not id1 or not id2:
        return False
    
    # Exact match
    if id1 == id2:
        return True
    
    # Base ID match
    base1 = extract_base_id(id1)
    base2 = extract_base_id(id2)
    if base1 == base2:
        return True
    
    # Normalized match
    if normalize_id(base1) == normalize_id(base2):
        return True
    
    # One contains the other (for complex IDs)
    if base1 in id2 or base2 in id1:
        return True
    
    return False


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
    with open(filepath, 'r') as f:
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
                # URL decode common escapes
                value = value.replace('%3B', ';').replace('%3D', '=')
                value = value.replace('%26', '&').replace('%2C', ',')
                attributes[key.strip()] = value.strip()
    
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


def translate(seq: str, table: int = 1) -> str:
    """
    Translate DNA to protein.
    
    Args:
        seq: DNA sequence
        table: Genetic code table (1 = standard, 11 = bacterial)
    """
    # Standard genetic code
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    seq = seq.upper().replace('U', 'T')
    protein = []
    
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if 'N' in codon:
            protein.append('X')
        else:
            protein.append(codon_table.get(codon, 'X'))
    
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
# BATCH ID ASSIGNMENT
# =============================================================================

class IdAssigner:
    """
    Assigns clean IDs to sequences, ensuring uniqueness.
    
    Use when you have many sequences and need consistent, unique IDs.
    """
    
    def __init__(self, prefix: str = "seq"):
        self.prefix = prefix
        self.seen_ids = {}
        self.counter = 0
    
    def get_id(self, original_header: str) -> str:
        """
        Get a unique ID for a sequence.
        
        If the extracted ID was seen before, appends a number.
        """
        # Try to extract ID
        base_id = extract_id(original_header)
        
        if base_id in self.seen_ids:
            # Already seen, add suffix
            self.seen_ids[base_id] += 1
            return f"{base_id}_{self.seen_ids[base_id]}"
        else:
            self.seen_ids[base_id] = 0
            return base_id
    
    def get_sequential_id(self) -> str:
        """Get next sequential ID."""
        self.counter += 1
        return f"{self.prefix}_{self.counter}"


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def load_genome(filepath: str) -> Dict[str, str]:
    """Load genome FASTA as dict {chrom: sequence}."""
    return {clean_id: seq for _, clean_id, seq in parse_fasta(filepath)}


def get_sequence_lengths(filepath: str) -> Dict[str, int]:
    """Get lengths of all sequences in FASTA."""
    return {clean_id: len(seq) for _, clean_id, seq in parse_fasta(filepath)}


def filter_fasta(input_path: str, output_path: str, keep_ids: set):
    """Filter FASTA to keep only specified IDs."""
    records = []
    for header, clean_id, seq in parse_fasta(input_path):
        if clean_id in keep_ids or extract_base_id(clean_id) in keep_ids:
            records.append((header, seq))
    write_fasta(records, output_path)


def merge_fastas(input_paths: List[str], output_path: str, deduplicate: bool = True):
    """Merge multiple FASTA files."""
    seen = set()
    records = []
    
    for path in input_paths:
        for header, clean_id, seq in parse_fasta(path):
            if deduplicate:
                if clean_id in seen:
                    continue
                seen.add(clean_id)
            records.append((header, seq))
    
    write_fasta(records, output_path)


# =============================================================================
# TESTING
# =============================================================================

if __name__ == '__main__':
    import sys
    
    # Test ID extraction
    test_headers = [
        "XP_123456.1 hypothetical protein",
        "sp|P12345|ABC_HUMAN Some protein",
        "ENSP00000123456.2 description",
        "lcl|NC_001234.5_cds_XP_789012.1_15 gene description",
        "gene|var_1",
        "1_234 # 1000 # 2000 # 1 # ID=1_234;partial=00",
        "random_header_with no known pattern",
        "ABC_12345 locus tag style",
        "",
        "   ",
    ]
    
    print("ID Extraction Tests:")
    print("-" * 60)
    for header in test_headers:
        extracted = extract_id(header)
        base = extract_base_id(extracted)
        print(f"  '{header[:50]}...' -> '{extracted}' (base: '{base}')")
    
    # Test file if provided
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
        print(f"\nParsing: {filepath}")
        print("-" * 60)
        for i, (header, clean_id, seq) in enumerate(parse_fasta(filepath)):
            if i >= 5:
                print(f"  ... and more")
                break
            print(f"  {clean_id}: {len(seq)} bp")
