#!/usr/bin/env python3

import sys
import os
import argparse
try:
    from ete3 import NCBITaxa
except ImportError:
    NCBITaxa = None

def main():
    parser = argparse.ArgumentParser(description="Sort genomes by phylogenetic distance")
    parser.add_argument("--home", required=True, help="Home genome name or TaxID")
    parser.add_argument("--targets", nargs='+', help="List of target genome files or names")
    parser.add_argument("--targets_dir", help="Directory containing target genomes")
    parser.add_argument("--img_ext", default=".fna", help="Extension for targets in directory (default: .fna)")
    parser.add_argument("--taxdb", required=True, help="Path to directory containing taxdump (nodes.dmp, names.dmp) or sqlite db")
    parser.add_argument("--output", required=True, help="Output sorted list of genomes")
    
    args = parser.parse_args()
    
    target_list = []
    if args.targets:
        target_list = args.targets
    elif args.targets_dir:
        # List files in directory
        if os.path.isdir(args.targets_dir):
             for f in os.listdir(args.targets_dir):
                 if f.endswith(args.img_ext) or f.endswith(".fasta") or f.endswith(".fa"):
                     target_list.append(os.path.join(args.targets_dir, f))
    
    if not target_list:
        print("Error: No targets provided via --targets or --targets_dir")
        sys.exit(1)
        
    # Standardize to basenames for output sorting consistency if needed, 
    # BUT we need full paths or at least consistent identifiers.
    # The downstream tools assume paths or names.
    # If we output just basenames, iterative_search needs to know the dir.
    # Let's output Basenames in the sorted list?
    # Original script wrote what it got.
    # We should stick to writing Basenames if we use targets_dir, to allow re-construction.
    
    # If NCBITaxa is not available, fallback immediately
    if NCBITaxa is None:
        print("ETE3 not installed. Falling back to alphabetical sorting.")
        fallback_sort(args)
        return

    # 1. Initialize NCBITaxa
    # ... (rest of code)
    # ETE3 uses a sqlite database. If args.taxdb serves a directory with dump files, 
    # we might need to specify the db file location or let ete3 build it.
    # However, standard NCBITaxa() uses ~/.etetoolkit/taxa.sqlite.
    # If the user provides a custom path, we should try to use it.
    # If args.taxdb is a file, assume it's the sqlite.
    # If it's a directory, we might need to update/create the db.
    # For this script to be robust in a pipeline, we'll try to use the provided DB or standard one.
    
    print(f"Initializing Taxonomy from {args.taxdb}...")
    
    try:
        if os.path.isfile(args.taxdb):
            ncbi = NCBITaxa(dbfile=args.taxdb)
        elif os.path.isdir(args.taxdb):
            # Check for taxdump
            db_path = os.path.join(args.taxdb, 'taxa.sqlite')
            if os.path.exists(db_path):
                 ncbi = NCBITaxa(dbfile=db_path)
            else:
                 # If we only have dmp files, ete3 constructs the DB.
                 # This might take time.
                 # Let's assume standard usage for now or fallback if fails.
                 # If the user passed the taxdump folder, we can point NCBITaxa to it?
                 # NCBITaxa(taxdump_file=...) takes a tar.gz usually.
                 # If uncompressed, we might need to rely on default behavior or custom handling.
                 # Let's try to load standard if provided arg is just a holder, 
                 # but arguably we should respect the input.
                 
                 # Optimization: If the user provides the taxdump FOLDER, we can't easily tell ETE3 to use it 
                 # without rebuilding the SQL. 
                 # Let's assume the user has set up ETE3 or provided a valid SQlite DB path.
                 # If not, and they provided the dump folder, we warn them.
                 pass
            # Fallback for now:
            ncbi = NCBITaxa(dbfile=db_path if os.path.exists(db_path) else None)
        else:
             # Try default
             ncbi = NCBITaxa()
    except Exception as e:
        print(f"Error initializing NCBITaxa: {e}")
        fallback_sort(args)
        return

    # 2. Resolve Home TaxID
    home_taxid = get_taxid(args.home, ncbi)
    if not home_taxid:
        print(f"Could not find TaxID for home: {args.home}")
        # Fallback?? Or Exit? Exit is safer for a "sort" tool but maybe fallback is better for pipeline continuity.
        fallback_sort(args)
        return

    print(f"Home TaxID: {home_taxid}")

    # 3. Process Targets
    scored_targets = []
    
    # Pre-fetch linege for home
    try:
        home_lineage = ncbi.get_lineage(home_taxid)
    except ValueError:
        print(f"TaxID {home_taxid} not found in DB")
        fallback_sort(args)
        return

    for target in target_list:
        tid = get_taxid(target, ncbi)
        if not tid:
            print(f"Warning: Could not resolve TaxID for {target}")
            # If target is a path, use basename for output if we want clean list
            scored_targets.append((float('inf'), os.path.basename(target)))
            continue
            
        try:
            # Simple lineage overlap measure:
            target_lineage = ncbi.get_lineage(tid)
            shared = 0
            for a, b in zip(home_lineage, target_lineage):
                if a == b:
                    shared += 1
                else:
                    break
            
            dist = (len(home_lineage) - shared) + (len(target_lineage) - shared)
            scored_targets.append((dist, os.path.basename(target)))
            
        except Exception as e:
            print(f"Error calculating distance for {target}: {e}")
            scored_targets.append((float('inf'), os.path.basename(target)))
            
    # 4. Sort and Write
    scored_targets.sort(key=lambda x: x[0])
    
    with open(args.output, 'w') as f:
        for dist, target in scored_targets:
            f.write(f"{target}\t{dist}\n")
            
def get_taxid(name, ncbi):
    # If name is integer-like, assume it's a taxid
    # Clean up name first
    clean = os.path.basename(name).split('.')[0].replace('_', ' ')
    # ... (rest of function)
    
    if clean.isdigit():
        return int(clean)
        
    # Lookup by name
    try:
        name2taxid = ncbi.get_name_translator([clean])
        if clean in name2taxid:
            return name2taxid[clean][0]
    except:
        pass
        
    # Try fuzzy or synonyms? ETE3 name translator is exact.
    # Try parts of name?
    # e.g. "Drosophila melanogaster release 6" -> "Drosophila melanogaster"
    parts = clean.split()
    for i in range(len(parts), 0, -1):
        subname = " ".join(parts[:i])
        try:
            name2taxid = ncbi.get_name_translator([subname])
            if subname in name2taxid:
                return name2taxid[subname][0]
        except:
            continue
            
    return None

def fallback_sort(args):
    print("Falling back to alphabetical sorting")
    # We need to re-generate target list because it's not passed to this function?
    # It is better to move target listing to main scope or pass it.
    # But args doesn't have it.
    # Let's verify if we can access target_list.
    # Changing function signature is needed or quick fix:
    
    target_list = []
    if args.targets:
        target_list = args.targets
    elif args.targets_dir:
        if os.path.isdir(args.targets_dir):
             for f in os.listdir(args.targets_dir):
                 if f.endswith(args.img_ext) or f.endswith(".fasta") or f.endswith(".fa"):
                     target_list.append(os.path.basename(f)) # Use basename for safety
    
    with open(args.output, 'w') as out:
        for target in sorted(target_list):
            out.write(f"{os.path.basename(target)}\t0\n")

if __name__ == "__main__":
    main()
