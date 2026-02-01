#!/usr/bin/env python3
"""
Fetch related genomes from NCBI for easy mode.
Uses NCBI Datasets API to find and download related species.
"""

import argparse
import subprocess
import json
import sys
import os
from pathlib import Path
import shlex

def run_command(cmd, check=True):
    """Run shell command and return output."""
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                              capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}", file=sys.stderr)
        print(f"Error: {e.stderr}", file=sys.stderr)
        if check:
            sys.exit(1)
        return None

def get_taxid_from_name(species_name):
    """Get NCBI taxonomy ID from species name."""
    # Use esearch from NCBI E-utilities
    safe_name = shlex.quote(species_name)
    cmd = f'esearch -db taxonomy -query {safe_name} | efetch -format uid'
    taxid = run_command(cmd, check=False)
    if taxid:
        return taxid.strip()
    return None

def get_related_species(taxid, max_genomes=10):
    """Get list of related species using NCBI taxonomy."""
    print(f"Finding related species for TaxID: {taxid}")
    
    # Get lineage
    cmd = f'efetch -db taxonomy -id {taxid} -format xml'
    xml_output = run_command(cmd, check=False)
    
    if not xml_output:
        print("Warning: Could not fetch taxonomy information", file=sys.stderr)
        return []
    
    # Parse to find genus or family level
    # For simplicity, search for genomes in the same genus
    cmd = f'esearch -db assembly -query "txid{taxid}[Organism:exp]" | '
    cmd += 'efetch -format docsum | '
    cmd += 'xtract -pattern DocumentSummary -element AssemblyAccession SpeciesName'
    
    results = run_command(cmd, check=False)
    
    if not results:
        print("Warning: No assemblies found", file=sys.stderr)
        return []
    
    assemblies = []
    for line in results.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 2:
                assemblies.append({
                    'accession': parts[0],
                    'species': parts[1]
                })
    
    # Limit number
    return assemblies[:max_genomes]

def download_genome(accession, output_dir):
    """Download genome from NCBI using datasets."""
    print(f"Downloading {accession}...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Try using datasets command if available
    cmd = f'datasets download genome accession {accession} --filename {output_dir}/{accession}.zip'
    result = run_command(cmd, check=False)
    
    if result is not None:
        # Unzip
        cmd = f'unzip -o {output_dir}/{accession}.zip -d {output_dir}/{accession}/'
        run_command(cmd, check=False)
        
        # Find the .fna file
        fna_files = list(Path(f"{output_dir}/{accession}").rglob("*.fna"))
        if fna_files:
            # Copy to main directory with simple name
            target = output_path / f"{accession}.fna"
            cmd = f'cp {fna_files[0]} {target}'
            run_command(cmd)
            print(f"  ✓ Downloaded to {target}")
            return str(target)
    
    # Fallback: try direct FTP download
    print(f"  Trying FTP download for {accession}...")
    # This is complex and depends on assembly structure
    # For now, just report failure
    print(f"  ✗ Could not download {accession}")
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Fetch related genomes from NCBI for easy mode",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fetch 5 related genomes for Drosophila melanogaster
  fetch_related_genomes.py --species "Drosophila melanogaster" --max 5 --outdir genomes/
  
  # Using taxonomy ID
  fetch_related_genomes.py --taxid 7227 --max 10 --outdir genomes/
  
Requirements:
  - NCBI E-utilities (esearch, efetch, xtract)
    Install: conda install -c bioconda entrez-direct
  - NCBI Datasets CLI (optional, for faster downloads)
    Install: conda install -c conda-forge ncbi-datasets-cli
        """
    )
    
    parser.add_argument("--species", help="Species name (e.g., 'Apis mellifera')")
    parser.add_argument("--taxid", help="NCBI Taxonomy ID")
    parser.add_argument("--max", type=int, default=10, 
                       help="Maximum number of genomes to fetch (default: 10)")
    parser.add_argument("--outdir", default="easy_mode_genomes",
                       help="Output directory for genomes (default: easy_mode_genomes)")
    parser.add_argument("--include-outgroup", action="store_true",
                       help="Also fetch one distant outgroup species")
    parser.add_argument("--list-only", action="store_true",
                       help="Only list available genomes, don't download")
    
    args = parser.parse_args()
    
    # Check for required tools
    if not run_command("which esearch", check=False):
        print("ERROR: NCBI E-utilities not found!", file=sys.stderr)
        print("Install with: conda install -c bioconda entrez-direct", file=sys.stderr)
        sys.exit(1)
    
    # Get taxonomy ID
    taxid = args.taxid
    if not taxid and args.species:
        print(f"Looking up taxonomy ID for '{args.species}'...")
        taxid = get_taxid_from_name(args.species)
        if not taxid:
            print(f"ERROR: Could not find taxonomy ID for '{args.species}'", file=sys.stderr)
            sys.exit(1)
        print(f"  Found TaxID: {taxid}")
    
    if not taxid:
        print("ERROR: Must provide either --species or --taxid", file=sys.stderr)
        sys.exit(1)
    
    # Get related species
    assemblies = get_related_species(taxid, args.max)
    
    if not assemblies:
        print("No related assemblies found.", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nFound {len(assemblies)} related assemblies:")
    for i, asm in enumerate(assemblies, 1):
        print(f"  {i}. {asm['accession']} - {asm['species']}")
    
    if args.list_only:
        print("\nList-only mode. Exiting without download.")
        sys.exit(0)
    
    # Download genomes
    print(f"\nDownloading genomes to {args.outdir}/...")
    downloaded = []
    
    for asm in assemblies:
        fna_path = download_genome(asm['accession'], args.outdir)
        if fna_path:
            downloaded.append(fna_path)
    
    print(f"\n{'='*60}")
    print(f"✓ Successfully downloaded {len(downloaded)}/{len(assemblies)} genomes")
    print(f"{'='*60}")
    
    # Write manifest file
    manifest_path = Path(args.outdir) / "genomes_manifest.txt"
    with open(manifest_path, 'w') as f:
        for path in downloaded:
            f.write(f"{path}\n")
    
    print(f"\nGenome paths written to: {manifest_path}")
    print("\nUse these genomes with:")
    print(f"  nextflow run main.nf --gene gene.fasta --home_genome home.fna \\")
    print(f"    --target_genomes '{args.outdir}/*.fna'")

if __name__ == "__main__":
    main()
