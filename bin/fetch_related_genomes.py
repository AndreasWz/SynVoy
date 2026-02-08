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

def get_parent_taxa(taxid):
    """
    Walk up the NCBI taxonomy tree from a given taxid.
    Returns list of (rank_label, taxid) for interesting ranks:
      family, order, class  (skipping genus since we already tried it).
    """
    target_ranks = ['family', 'order', 'class']
    results = []
    try:
        # efetch the full lineage for this taxid
        cmds = [
            ['efetch', '-db', 'taxonomy', '-id', str(taxid), '-format', 'xml'],
            ['xtract', '-pattern', 'Taxon', '-block', 'LineageEx/Taxon',
             '-element', 'Rank', '-element', 'TaxId', '-element', 'ScientificName'],
        ]
        raw = run_piped_command(cmds)
        if not raw:
            return results
        # Output: lines of "rank\ttaxid\tname" repeated for every ancestor
        parts = raw.split('\t')
        # Group into triples
        triples = []
        for i in range(0, len(parts) - 2, 3):
            triples.append((parts[i].strip(), parts[i+1].strip(), parts[i+2].strip()))
        for rank, tid, name in triples:
            if rank.lower() in target_ranks:
                results.append((f"{rank} ({name})", tid))
        # Sort by target_ranks order (family first, then order, then class)
        rank_order = {r: i for i, r in enumerate(target_ranks)}
        results.sort(key=lambda x: rank_order.get(x[0].split()[0].lower(), 99))
    except Exception as e:
        print(f"Warning: Could not retrieve lineage for TaxID {taxid}: {e}")
    return results

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
    """Download genome from NCBI using datasets. Also attempts to download GFF if available."""
    print(f"Downloading {accession}...")
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    zip_file = output_path / f"{accession}.zip"
    
    # Try using datasets command — include GFF annotations if available
    cmd = ['datasets', 'download', 'genome', 'accession', accession,
           '--include', 'genome,gff3', '--filename', str(zip_file)]
    
    try:
        run_safe_command(cmd)
        
        # Unzip
        extract_dir = output_path / accession
        cmd_unzip = ['unzip', '-o', str(zip_file), '-d', str(extract_dir)]
        run_safe_command(cmd_unzip)
        
        # Find .fna file
        fna_files = list(extract_dir.rglob("*.fna"))
        fna_path = None
        if fna_files:
            target = output_path / f"{accession}.fna"
            os.rename(fna_files[0], target)
            fna_path = str(target)
            print(f"  ✓ Genome: {target}")
        
        # Find .gff file (annotations)
        gff_files = list(extract_dir.rglob("*.gff"))
        if gff_files:
            # Check if GFF has actual CDS features (not just scaffold entries)
            best_gff = None
            for gff in gff_files:
                try:
                    with open(gff) as f:
                        has_cds = False
                        for line_num, line in enumerate(f):
                            if line_num > 500:
                                break
                            if '\tCDS\t' in line or '\tgene\t' in line:
                                has_cds = True
                                break
                        if has_cds:
                            best_gff = gff
                            break
                except Exception:
                    continue
            
            if best_gff:
                gff_target = output_path / f"{accession}.gff"
                os.rename(best_gff, gff_target)
                print(f"  ✓ GFF annotations: {gff_target}")
            else:
                print(f"  ○ GFF found but no CDS features (scaffold-only)")
        else:
            print(f"  ○ No GFF annotations available")
        
        # Cleanup zip and extracted dir
        if zip_file.exists(): zip_file.unlink()
        import shutil
        if extract_dir.exists(): shutil.rmtree(extract_dir)
            
        return fna_path
            
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
    parser.add_argument("--target-species", dest="target_species", default=None,
                       help="Comma-separated list of species names to fetch (bypasses automatic search)")
    
    args = parser.parse_args()
    
    # Check for required tools
    try:
        run_safe_command(['which', 'esearch'])
    except:
        print("ERROR: NCBI E-utilities (esearch) not found!", file=sys.stderr)
        sys.exit(1)
    
    # --- TARGET SPECIES MODE ---
    # If user provided specific species names, fetch those directly
    if args.target_species:
        species_list = [s.strip() for s in args.target_species.split(',') if s.strip()]
        print(f"Target species mode: fetching {len(species_list)} specified species")
        
        assemblies = []
        for sp_name in species_list:
            print(f"\nLooking up '{sp_name}'...")
            sp_taxid = get_taxid_from_name(sp_name)
            if not sp_taxid:
                print(f"  WARNING: Could not find taxonomy ID for '{sp_name}', skipping")
                continue
            
            # Search for assemblies of this specific species
            sp_assemblies = get_related_species(sp_taxid, max_genomes=3, refseq_only=False,
                                                 exclude_species=args.home_species)
            # Also try exact species name search
            if not sp_assemblies:
                sp_assemblies = get_related_species(sp_taxid, max_genomes=3, refseq_only=False,
                                                     exclude_species=None)
            
            if sp_assemblies:
                # Keep best assembly for this species
                best = sp_assemblies[0]
                assemblies.append(best)
                print(f"  Found: {best['accession']} ({best['species']}, {best['category']})")
            else:
                print(f"  WARNING: No assemblies found for '{sp_name}'")
        
        if not assemblies:
            print("ERROR: No assemblies found for any of the specified species", file=sys.stderr)
            output_path = Path(args.outdir)
            output_path.mkdir(parents=True, exist_ok=True)
            (output_path / "genomes_manifest.txt").touch()
            (output_path / "species_mapping.tsv").touch()
            sys.exit(0)
        
        print(f"\nFound {len(assemblies)} assemblies for specified species:")
        for i, asm in enumerate(assemblies, 1):
            print(f"  {i}. {asm['accession']} - {asm['species']}")
        
        if not args.list_only:
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
            
            # Write species mapping file
            species_map_path = output_path / "species_mapping.tsv"
            with open(species_map_path, 'w') as f:
                for asm in assemblies:
                    f.write(f"{asm['accession']}\t{asm['species']}\n")
            print(f"Species mapping written to: {species_map_path}")
        
        return
    
    # --- AUTOMATIC TAXONOMIC SEARCH MODE ---
    
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
    
    # If no assemblies found at genus level, try broader taxonomic levels
    # Walk up the taxonomy tree: genus -> family -> order
    if not assemblies:
        print(f"No non-home assemblies found at genus level. Trying broader search...")
        
        # Get lineage from NCBI taxonomy automatically
        parent_ranks = get_parent_taxa(taxid)
        
        for rank_name, rank_taxid in parent_ranks:
            print(f"Searching {rank_name} (TaxID: {rank_taxid})...")
            broader_count = max(50, args.max * 5)
            assemblies = get_related_species(rank_taxid, broader_count,
                                              refseq_only=True, exclude_species=args.home_species)
            if assemblies:
                print(f"  Found {len(assemblies)} assemblies at {rank_name} level")
                break
            print(f"  No assemblies at {rank_name} level either")
    
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

    # Write species mapping file (accession -> species name)
    species_map_path = output_path / "species_mapping.tsv"
    with open(species_map_path, 'w') as f:
        for asm in assemblies:
            f.write(f"{asm['accession']}\t{asm['species']}\n")
    print(f"Species mapping written to: {species_map_path}")

if __name__ == "__main__":
    main()
