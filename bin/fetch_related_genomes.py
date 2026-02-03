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

def run_safe_command(cmd, check=True):
    """Run command with list arguments safely."""
    try:
        # print(f"DEBUG: running command: {cmd}")
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Error: {e.stderr}", file=sys.stderr)
        if check:
            # Propagate error if check=True
            raise e
        return None

def run_piped_command(cmds):
    """
    Run a chain of commands connected by pipes.
    cmds: List of command lists. E.g. [['esearch', ...], ['efetch', ...]]
    """
    procs = []
    # Start first process
    try:
        p1 = subprocess.Popen(cmds[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procs.append(p1)
        
        # Chain subsequent processes
        prev_stdout = p1.stdout
        for i in range(1, len(cmds)):
            p_next = subprocess.Popen(cmds[i], stdin=prev_stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            procs.append(p_next)
            # Allow p1 to receive SIGPIPE if p2 exits
            prev_stdout.close()
            prev_stdout = p_next.stdout
            
        # Get output from last process
        last_proc = procs[-1]
        output, error = last_proc.communicate()
        
        if last_proc.returncode != 0:
            print(f"Error in piped command chain: {cmds[-1]}", file=sys.stderr)
            print(f"Stderr: {error.decode() if error else ''}", file=sys.stderr)
            return None
            
        return output.decode().strip()
            
    except Exception as e:
        print(f"Exception running pipe chain: {e}", file=sys.stderr)
        return None
    finally:
        # Ensure cleanup
        for p in procs:
            if p.poll() is None:
                p.kill()

def get_taxid_from_name(species_name):
    """Get NCBI taxonomy ID from species name."""
    # esearch -db taxonomy -query name | efetch -format uid
    cmds = [
        ['esearch', '-db', 'taxonomy', '-query', species_name],
        ['efetch', '-format', 'uid']
    ]
    return run_piped_command(cmds)

def get_related_species(taxid, max_genomes=10):
    """Get list of related species using NCBI taxonomy."""
    print(f"Finding related species for TaxID: {taxid}")
    
    # Get lineage check (optional, can skip directly to esearch)
    
    # esearch -db assembly -query "txid{taxid}[Organism:exp]" | efetch -format docsum | xtract ...
    query = f"txid{taxid}[Organism:exp]"
    
    cmds = [
        ['esearch', '-db', 'assembly', '-query', query],
        ['efetch', '-format', 'docsum'],
        ['xtract', '-pattern', 'DocumentSummary', '-element', 'AssemblyAccession', 'SpeciesName']
    ]
    
    results = run_piped_command(cmds)
    
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
    
    return assemblies[:max_genomes]

def download_genome(accession, output_dir):
    """Download genome from NCBI using datasets."""
    print(f"Downloading {accession}...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    zip_file = output_path / f"{accession}.zip"
    
    # Try using datasets command
    cmd = ['datasets', 'download', 'genome', 'accession', accession, '--filename', str(zip_file)]
    
    try:
        run_safe_command(cmd)
        
        # Unzip
        # unzip -o file.zip -d output_dir/result
        extract_dir = output_path / accession
        cmd_unzip = ['unzip', '-o', str(zip_file), '-d', str(extract_dir)]
        run_safe_command(cmd_unzip)
        
        # Find .fna file
        fna_files = list(extract_dir.rglob("*.fna"))
        if fna_files:
            target = output_path / f"{accession}.fna"
            # Move/Copy
            os.rename(fna_files[0], target)
            print(f"  ✓ Downloaded to {target}")
            
            # Cleanup zip and extracted dir
            if zip_file.exists(): zip_file.unlink()
            import shutil
            if extract_dir.exists(): shutil.rmtree(extract_dir)
                
            return str(target)
            
    except Exception as e:
        print(f"  ✗ Could not download {accession}: {e}")
    
    return None

def main():
    parser = argparse.ArgumentParser(
        description="Fetch related genomes from NCBI for easy mode",
        formatter_class=argparse.RawDescriptionHelpFormatter
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
    try:
        run_safe_command(['which', 'esearch'])
    except:
        print("ERROR: NCBI E-utilities (esearch) not found!", file=sys.stderr)
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

if __name__ == "__main__":
    main()
