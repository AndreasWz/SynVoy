#!/usr/bin/env python3
"""
Fetch the reference genome for a species from NCBI.
Downloads both the genome FASTA and GFF annotation.
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path
import shutil


def run_command(cmd, check=True):
    """Run command and return output."""
    try:
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}", file=sys.stderr)
        if check:
            raise
        return None


def get_reference_genome(species_name):
    """Get the reference/representative genome accession for a species."""
    print(f"Finding reference genome for '{species_name}'...")
    
    # Use datasets to find reference genome
    # datasets summary genome taxon "Apis mellifera" --reference --as-json-lines
    cmd = [
        'datasets', 'summary', 'genome', 'taxon', species_name,
        '--reference', '--as-json-lines'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        import json
        
        for line in result.stdout.strip().split('\n'):
            if line:
                data = json.loads(line)
                
                # JSON structure: accession is at root level with --as-json-lines
                if 'accession' in data:
                    accession = data.get('accession')
                    org_name = data.get('organism', {}).get('organism_name', species_name)
                    annotation = data.get('annotation_info', {})
                    has_annotation = annotation.get('status') == 'Full annotation' if annotation else False
                    
                    print(f"  Found: {accession} ({org_name})")
                    print(f"  Has annotation: {has_annotation}")
                    return accession, has_annotation
                
                # Fallback: check for 'reports' array (older API format)
                elif 'reports' in data:
                    for report in data['reports']:
                        accession = report.get('accession')
                        org_name = report.get('organism', {}).get('organism_name', species_name)
                        annotation = report.get('annotation_info', {})
                        has_annotation = annotation.get('status') == 'Full annotation' if annotation else False
                        
                        print(f"  Found: {accession} ({org_name})")
                        print(f"  Has annotation: {has_annotation}")
                        return accession, has_annotation
                        
    except subprocess.CalledProcessError as e:
        print(f"Error querying NCBI: {e.stderr}", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
    
    return None, False


def download_genome_with_annotation(accession, output_dir):
    """Download genome and GFF annotation using NCBI datasets."""
    print(f"\nDownloading {accession} with annotation...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    zip_file = output_path / f"{accession}.zip"
    
    # Download with annotation
    cmd = [
        'datasets', 'download', 'genome', 'accession', accession,
        '--include', 'genome,gff3',
        '--filename', str(zip_file)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Unzip
        extract_dir = output_path / 'extracted'
        subprocess.run(['unzip', '-o', str(zip_file), '-d', str(extract_dir)], 
                      check=True, capture_output=True)
        
        # Find files
        fna_files = list(extract_dir.rglob("*.fna"))
        gff_files = list(extract_dir.rglob("*.gff"))
        
        genome_path = None
        gff_path = None
        
        if fna_files:
            genome_path = output_path / "home_genome.fna"
            shutil.copy(fna_files[0], genome_path)
            print(f"  Genome: {genome_path}")
        
        if gff_files:
            gff_path = output_path / "home_genome.gff"
            shutil.copy(gff_files[0], gff_path)
            print(f"  Annotation: {gff_path}")
        
        # Cleanup
        if zip_file.exists():
            zip_file.unlink()
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        
        return genome_path, gff_path
        
    except subprocess.CalledProcessError as e:
        print(f"Error downloading: {e.stderr}", file=sys.stderr)
        return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Fetch reference genome for a species from NCBI"
    )
    
    parser.add_argument("--species", required=True,
                       help="Species name (e.g., 'Apis mellifera')")
    parser.add_argument("--outdir", default="home_genome",
                       help="Output directory (default: home_genome)")
    
    args = parser.parse_args()
    
    # Check for datasets tool
    try:
        subprocess.run(['which', 'datasets'], check=True, capture_output=True)
    except:
        print("ERROR: NCBI datasets tool not found!", file=sys.stderr)
        print("Install with: conda install -c conda-forge ncbi-datasets-cli", file=sys.stderr)
        sys.exit(1)
    
    # Find reference genome
    accession, has_annotation = get_reference_genome(args.species)
    
    if not accession:
        print(f"ERROR: Could not find reference genome for '{args.species}'", file=sys.stderr)
        sys.exit(1)
    
    # Download
    genome_path, gff_path = download_genome_with_annotation(accession, args.outdir)
    
    if not genome_path:
        print("ERROR: Failed to download genome", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n{'='*60}")
    print(f"SUCCESS: Downloaded reference genome for {args.species}")
    print(f"  Accession: {accession}")
    print(f"  Genome:    {genome_path}")
    if gff_path:
        print(f"  GFF:       {gff_path}")
    print(f"{'='*60}")
    
    # Write paths to files for Nextflow to pick up
    (Path(args.outdir) / "genome_path.txt").write_text(str(genome_path))
    if gff_path:
        (Path(args.outdir) / "gff_path.txt").write_text(str(gff_path))


if __name__ == "__main__":
    main()
