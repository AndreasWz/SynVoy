#!/usr/bin/env python3
"""
Resolve gene input: accepts UniProt ID, NCBI protein accession, or FASTA file.
Outputs JSON with resolved FASTA path + species name.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path


# ─── Input Type Detection ─────────────────────────────────────────────────────

# UniProt accession patterns (https://www.uniprot.org/help/accession_numbers)
UNIPROT_RE = re.compile(
    r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|'
    r'^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'
)

# NCBI RefSeq protein accessions
NCBI_REFSEQ_RE = re.compile(r'^[XNWY]P_\d+\.\d+$')

# NCBI GenBank protein accessions  (e.g., KAF1234567.1, AAB12345.1)
NCBI_GENBANK_RE = re.compile(r'^[A-Z]{3}\d{5,}\.\d+$')


def detect_input_type(gene_input):
    """Detect whether input is a file, UniProt ID, or NCBI accession."""
    # Check if it's a file on disk
    if os.path.isfile(gene_input):
        return 'file', gene_input
    
    # Strip version suffix for matching (e.g., P01501.2 → P01501)
    clean = gene_input.strip().split('.')[0] if '.' in gene_input else gene_input.strip()
    
    if UNIPROT_RE.match(clean):
        return 'uniprot', gene_input.strip()
    
    if NCBI_REFSEQ_RE.match(gene_input.strip()):
        return 'ncbi', gene_input.strip()
    
    if NCBI_GENBANK_RE.match(gene_input.strip()):
        return 'ncbi', gene_input.strip()
    
    # If nothing matched, and it's NOT a file, treat as symbol
    return 'symbol', gene_input


# ─── UniProt Fetching ─────────────────────────────────────────────────────────

def fetch_uniprot(uniprot_id, output_dir):
    """Fetch sequence and species from UniProt REST API."""
    base_url = "https://rest.uniprot.org/uniprotkb"
    
    # Fetch JSON metadata (species, protein name, etc.)
    json_url = f"{base_url}/{uniprot_id}.json"
    print(f"  Fetching metadata from UniProt: {json_url}")
    
    try:
        req = urllib.request.Request(json_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  ERROR: UniProt ID '{uniprot_id}' not found", file=sys.stderr)
            return None
        raise
    
    # Extract species
    organism = data.get('organism', {})
    species = organism.get('scientificName', '')
    taxid = organism.get('taxonId', '')
    
    # Extract protein name
    protein_desc = data.get('proteinDescription', {})
    rec_name = protein_desc.get('recommendedName', {})
    protein_name = rec_name.get('fullName', {}).get('value', '')
    if not protein_name:
        # Try submittedName
        sub_names = protein_desc.get('submissionNames', [])
        if sub_names:
            protein_name = sub_names[0].get('fullName', {}).get('value', '')
    
    # Fetch FASTA
    fasta_url = f"{base_url}/{uniprot_id}.fasta"
    print(f"  Fetching FASTA from UniProt: {fasta_url}")
    
    try:
        with urllib.request.urlopen(fasta_url, timeout=30) as resp:
            fasta_content = resp.read().decode()
    except urllib.error.HTTPError:
        print(f"  ERROR: Could not fetch FASTA for '{uniprot_id}'", file=sys.stderr)
        return None
    
    if not fasta_content.strip():
        print(f"  ERROR: Empty FASTA for '{uniprot_id}'", file=sys.stderr)
        return None
    
    # Write FASTA
    fasta_path = Path(output_dir) / f"{uniprot_id}.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_path.write_text(fasta_content)
    
    print(f"  Species: {species}")
    print(f"  Protein: {protein_name}")
    print(f"  TaxID:   {taxid}")
    print(f"  FASTA:   {fasta_path}")
    
    return {
        'fasta_path': str(fasta_path),
        'species': species,
        'taxid': str(taxid),
        'protein_name': protein_name,
        'source': 'uniprot',
        'input_id': uniprot_id
    }


def search_uniprot(gene_name, species, output_dir):
    """Search UniProt for a gene name + species."""
    print(f"  Searching UniProt for gene '{gene_name}' in '{species}'...")
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # URL encode query
    query = f"gene_exact:{gene_name} AND organism_name:\"{species}\" AND reviewed:true"
    params = urllib.parse.urlencode({'query': query, 'format': 'json', 'size': 1})
    search_url = f"{base_url}?{params}"
    
    try:
        req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  ERROR: UniProt search failed: {e}", file=sys.stderr)
        return None
        
    results = data.get('results', [])
    if not results:
        # Try unreviewed
        print("  No Swiss-Prot hit, trying TrEMBL...")
        query = f"gene_exact:{gene_name} AND organism_name:\"{species}\""
        params = urllib.parse.urlencode({'query': query, 'format': 'json', 'size': 1})
        search_url = f"{base_url}?{params}"
        try:
            req = urllib.request.Request(search_url, headers={'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            results = data.get('results', [])
        except Exception:
            pass
            
    if not results:
        print(f"  ERROR: No results found for gene '{gene_name}' in '{species}'", file=sys.stderr)
        return None
        
    # Take first hit
    hit = results[0]
    accession = hit.get('primaryAccession')
    print(f"  Found accession: {accession}")
    
    return fetch_uniprot(accession, output_dir)


# ─── NCBI Fetching ───────────────────────────────────────────────────────────

def fetch_ncbi(accession, output_dir):
    """Fetch sequence and species from NCBI."""
    print(f"  Fetching from NCBI: {accession}")
    
    # Fetch FASTA
    try:
        result = subprocess.run(
            ['efetch', '-db', 'protein', '-id', accession, '-format', 'fasta'],
            capture_output=True, text=True, check=True, timeout=30
        )
        fasta_content = result.stdout.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"  ERROR: Could not fetch FASTA from NCBI: {e}", file=sys.stderr)
        return None
    
    if not fasta_content:
        print(f"  ERROR: Empty FASTA for '{accession}'", file=sys.stderr)
        return None
    
    # Fetch metadata (organism name)
    species = ''
    protein_name = ''
    taxid = ''
    try:
        result = subprocess.run(
            ['efetch', '-db', 'protein', '-id', accession, '-format', 'gpc'],
            capture_output=True, text=True, check=True, timeout=30
        )
        gpc_output = result.stdout
        
        # Parse organism from GPC XML
        import xml.etree.ElementTree as ET
        try:
            root = ET.fromstring(gpc_output)
            # Look for INSDSeq_organism
            for elem in root.iter():
                if elem.tag == 'INSDSeq_organism' and elem.text:
                    species = elem.text.strip()
                elif elem.tag == 'INSDSeq_definition' and elem.text:
                    protein_name = elem.text.strip()
                elif elem.tag == 'INSDSeq_taxonomy' and elem.text:
                    pass  # Full taxonomy string, not needed
        except ET.ParseError:
            pass
        
        # Get TaxID via esearch if we have species
        if species:
            try:
                tax_result = subprocess.run(
                    ['esearch', '-db', 'taxonomy', '-query', species],
                    capture_output=True, text=True, timeout=15,
                    stdin=subprocess.DEVNULL
                )
                efetch_result = subprocess.run(
                    ['efetch', '-format', 'uid'],
                    input=tax_result.stdout, capture_output=True, text=True, timeout=15
                )
                taxid = efetch_result.stdout.strip()
            except Exception:
                pass
                
    except Exception as e:
        print(f"  Warning: Could not fetch metadata: {e}", file=sys.stderr)
    
    # Write FASTA
    safe_name = accession.replace('.', '_')
    fasta_path = Path(output_dir) / f"{safe_name}.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    fasta_path.write_text(fasta_content + '\n')
    
    print(f"  Species: {species}")
    print(f"  Protein: {protein_name}")
    print(f"  TaxID:   {taxid}")
    print(f"  FASTA:   {fasta_path}")
    
    return {
        'fasta_path': str(fasta_path),
        'species': species,
        'taxid': taxid,
        'protein_name': protein_name,
        'source': 'ncbi',
        'input_id': accession
    }


# ─── FASTA File Handling ──────────────────────────────────────────────────────

def handle_fasta_file(fasta_path, output_dir, species=None):
    """Handle a local FASTA file."""
    if not os.path.isfile(fasta_path):
        print(f"ERROR: File not found: {fasta_path}", file=sys.stderr)
        return None

    # Keep all resolved inputs in one place with a predictable extension so
    # Nextflow output collection (`resolved_query/*.fasta`) always works.
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    copied_fasta = outdir / "input_query.fasta"
    shutil.copyfile(fasta_path, copied_fasta)
    
    return {
        'fasta_path': str(copied_fasta.resolve()),
        'species': species or '',
        'taxid': '',
        'protein_name': '',
        'source': 'file',
        'input_id': Path(fasta_path).name
    }


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Resolve gene input: UniProt ID, NCBI accession, or FASTA file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input P01501                        # UniProt ID
  %(prog)s --input XP_006565763.1                # NCBI RefSeq protein
  %(prog)s --input my_protein.fasta --species "Apis mellifera"  # FASTA file
"""
    )
    
    parser.add_argument("--input", required=True,
                       help="Gene input: UniProt ID, NCBI protein accession, or path to FASTA file")
    parser.add_argument("--species", default=None,
                       help="Override species name (required for FASTA files without --query_id)")
    parser.add_argument("--outdir", default="resolved_query",
                       help="Output directory for resolved FASTA (default: resolved_query)")
    
    args = parser.parse_args()
    
    # Detect input type
    input_type, clean_input = detect_input_type(args.input)
    print(f"Input: '{args.input}' → detected as: {input_type}")
    
    # Resolve
    result = None
    
    if input_type == 'uniprot':
        result = fetch_uniprot(clean_input, args.outdir)
        
    elif input_type == 'ncbi':
        result = fetch_ncbi(clean_input, args.outdir)
        
    elif input_type == 'file':
        result = handle_fasta_file(clean_input, output_dir=args.outdir, species=args.species)

    elif input_type == 'symbol':
        # If it's a symbol, we need a species to search
        if not args.species:
            # Fallback: Maybe it WAS a file but missing?
            print(f"ERROR: Input '{clean_input}' is not a file and no species provided for search.", file=sys.stderr)
            sys.exit(1)
        
        result = search_uniprot(clean_input, args.species, args.outdir)
    
    if result is None:
        print("ERROR: Failed to resolve gene input", file=sys.stderr)
        sys.exit(1)
    
    # Override species if provided on command line
    if args.species:
        result['species'] = args.species
    
    # Validate: if no species, warn (but don't fail — main.nf handles this)
    if not result['species']:
        print("WARNING: No species detected. --home_species will be required.", file=sys.stderr)
    
    # Write JSON result
    output_path = Path(args.outdir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    json_path = output_path / "resolved_input.json"
    with open(json_path, 'w') as f:
        json.dump(result, f, indent=2)
    
    # Also write plain-text outputs for easy Nextflow consumption
    (output_path / "resolved_species.txt").write_text(result['species'])
    (output_path / "resolved_fasta.txt").write_text(result['fasta_path'])
    
    print(f"\n{'='*60}")
    print(f"Resolved: {result['source']} → {result['species'] or '(no species)'}")
    print(f"  FASTA: {result['fasta_path']}")
    print(f"  JSON:  {json_path}")
    print(f"{'='*60}")
    
    # Print JSON to stdout for piping
    print(json.dumps(result))


if __name__ == "__main__":
    main()
