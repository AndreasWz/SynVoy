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

def get_related_species(taxid, max_genomes=10, refseq_only=True, exclude_species=None):
    """Get list of related species using NCBI taxonomy.
    
    Args:
        taxid: NCBI taxonomy ID
        max_genomes: Maximum number of genomes to return
        refseq_only: Only return RefSeq assemblies (GCF_ prefix) - higher quality
        exclude_species: Species name to exclude (e.g., home species)
    """
    print(f"Finding related species for TaxID: {taxid}")
    
    # Query NCBI and get assembly info including RefSeq category
    query = f"txid{taxid}[Organism:exp]"
    print(f"  Query: {query}")
    
    cmds = [
        ['esearch', '-db', 'assembly', '-query', query],
        ['efetch', '-format', 'docsum'],
        ['xtract', '-pattern', 'DocumentSummary', '-element', 
         'AssemblyAccession', 'SpeciesName', 'RefSeq_category']
    ]
    
    results = run_piped_command(cmds)
    
    if not results:
        print("Warning: No assemblies found", file=sys.stderr)
        return []
    
    exclude_lower = exclude_species.lower().strip() if exclude_species else None
    
    # Parse results and categorize
    reference_genomes = []
    other_refseq = []
    genbank_only = []
    
    for line in results.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            accession = parts[0]
            species = parts[1]
            refseq_cat = parts[2] if len(parts) > 2 else 'na'
            
            # Skip excluded species
            if exclude_lower and species.lower().strip() == exclude_lower:
                continue
            
            entry = {'accession': accession, 'species': species, 'category': refseq_cat}
            
            is_refseq = accession.startswith('GCF_')
            is_reference = refseq_cat in ('reference genome', 'representative genome')
            
            if is_refseq and is_reference:
                reference_genomes.append(entry)
            elif is_refseq:
                other_refseq.append(entry)
            elif not refseq_only:
                genbank_only.append(entry)
    
    # Prioritize: reference genomes > other RefSeq > GenBank
    assemblies = reference_genomes + other_refseq + genbank_only
    
    # Deduplicate by species (keep first/best per species)
    seen_species = set()
    unique_assemblies = []
    for asm in assemblies:
        sp = asm['species'].lower().strip()
        if sp not in seen_species:
            seen_species.add(sp)
            unique_assemblies.append(asm)
            if len(unique_assemblies) >= max_genomes:
                break
    
    print(f"  Found {len(unique_assemblies)} unique species (excluding {exclude_species or 'none'})")
    if unique_assemblies:
        print(f"  Reference/representative genomes: {len([a for a in unique_assemblies if a['category'] in ('reference genome', 'representative genome')])}")
    
    return unique_assemblies

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
    
    parser.add_argument("--home-species", dest="home_species", required=True,
                       help="Home species name (e.g., 'Apis mellifera') - searches genus, excludes this species")
    parser.add_argument("--max", type=int, default=10, 
                       help="Maximum number of genomes to fetch (default: 10)")
    parser.add_argument("--outdir", default="easy_mode_genomes",
                       help="Output directory for genomes (default: easy_mode_genomes)")
    parser.add_argument("--list-only", action="store_true",
                       help="Only list available genomes, don't download")
    
    args = parser.parse_args()
    
    # Check for required tools
    try:
        run_safe_command(['which', 'esearch'])
    except:
        print("ERROR: NCBI E-utilities (esearch) not found!", file=sys.stderr)
        sys.exit(1)
    
    # Extract genus from home species (first word)
    genus = args.home_species.split()[0]
    print(f"Home species: {args.home_species}")
    print(f"Searching genus: {genus}")
    
    # Get taxonomy ID for genus
    print(f"Looking up taxonomy ID for '{genus}'...")
    taxid = get_taxid_from_name(genus)
    if not taxid:
        print(f"ERROR: Could not find taxonomy ID for '{genus}'", file=sys.stderr)
        sys.exit(1)
    print(f"  Found TaxID: {taxid}")
    
    # Get related species - pass exclude_species to filter early
    assemblies = get_related_species(taxid, args.max, refseq_only=True, exclude_species=args.home_species)
    
    # If no assemblies found at genus level, try family level
    # This is important for cases like Homo (only has extinct relatives)
    # or genera with few species
    if not assemblies:
        print(f"No non-home assemblies found at genus level. Trying broader search...")
        
        # Known family mappings for common problematic cases
        family_mappings = {
            'homo': 'Hominidae',  # Great apes
            'pan': 'Hominidae',   # Chimps
            'gorilla': 'Hominidae',
            'pongo': 'Hominidae',  # Orangutans
        }
        
        family = family_mappings.get(genus.lower())
        if family:
            print(f"Searching family: {family}")
            family_taxid = get_taxid_from_name(family)
            if family_taxid:
                print(f"  Found family TaxID: {family_taxid}")
                # Fetch more for family-level search
                family_fetch_count = max(50, args.max * 5)
                assemblies = get_related_species(family_taxid, family_fetch_count, 
                                                  refseq_only=True, exclude_species=args.home_species)
    
    if not assemblies:
        print("WARNING: No related assemblies found.", file=sys.stderr)
        # Create empty manifest to avoid downstream errors
        output_path = Path(args.outdir)
        output_path.mkdir(parents=True, exist_ok=True)
        manifest_path = output_path / "genomes_manifest.txt"
        manifest_path.touch()
        print(f"Created empty manifest: {manifest_path}")
        sys.exit(0)  # Exit gracefully, not with error
    
    # assemblies already filtered above; limit to requested max
    assemblies = assemblies[:args.max]
    
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
    
    # Ensure output directory exists
    output_path = Path(args.outdir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Write manifest file
    manifest_path = output_path / "genomes_manifest.txt"
    with open(manifest_path, 'w') as f:
        for path in downloaded:
            f.write(f"{path}\n")
    
    print(f"\nGenome paths written to: {manifest_path}")

if __name__ == "__main__":
    main()
