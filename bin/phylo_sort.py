#!/usr/bin/env python3
"""
phylo_sort.py - Sort target genomes by phylogenetic distance from home genome.

Strategies (tried in order):
  1. ETE3 + local taxonomy database (offline, fast)
  2. NCBI Datasets CLI - extract taxonomy from GCF_/GCA_ accessions (online)
  3. Alphabetical fallback with distance=0

The key insight: in pro mode, genome filenames typically contain NCBI assembly
accessions (e.g. GCF_029169275.1.fna). We can extract taxonomy lineages
from these accessions via the NCBI Datasets CLI, which is already a pipeline
dependency.
"""

import hashlib
import os
import re
import json
import subprocess
import argparse
try:
    from ete3 import NCBITaxa
except ImportError:
    NCBITaxa = None


# =============================================================================
# ACCESSION-BASED TAXONOMY (via NCBI Datasets CLI)
# =============================================================================

def _extract_accession(filename):
    """
    Extract an NCBI assembly accession (GCF_/GCA_ prefix) from a filename.
    Returns accession string or None.
    Examples:
        GCF_029169275.1.fna -> GCF_029169275.1
        GCA_928718305.1.fasta.gz -> GCA_928718305.1
    """
    basename = os.path.basename(filename)
    m = re.match(r'(GC[AF]_\d{9}\.\d+)', basename)
    return m.group(1) if m else None


def _datasets_cli_available():
    """Check if NCBI datasets CLI is available."""
    try:
        result = subprocess.run(
            ['datasets', '--version'],
            capture_output=True, text=True, timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _fetch_lineage_from_accessions(accessions):
    """
    Fetch taxonomy lineage for a list of assembly accessions using NCBI
    datasets CLI.

    Two-step approach:
      1. `datasets summary genome accession` → get tax_id per accession
      2. `datasets summary taxonomy taxon` → get full lineage (parents list)

    Returns dict: accession -> list of parent tax IDs (from root to species).
    """
    if not accessions:
        return {}

    lineages = {}

    # Step 1: Resolve accessions to tax IDs
    acc_to_taxid = {}
    acc_to_name = {}
    try:
        cmd = ['datasets', 'summary', 'genome', 'accession'] + list(accessions)
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            print(f"[phylo_sort] datasets genome summary error: {result.stderr[:200]}")
            return {}

        data = json.loads(result.stdout)
        for report in data.get('reports', []):
            acc = report.get('accession', '')
            org = report.get('organism', {})
            tax_id = org.get('tax_id', 0)
            species_name = org.get('organism_name', '')
            if acc and tax_id:
                acc_to_taxid[acc] = tax_id
                acc_to_name[acc] = species_name

    except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError) as e:
        print(f"[phylo_sort] Error in genome summary: {e}")
        return {}

    if not acc_to_taxid:
        return {}

    # Step 2: Get full lineage for each tax_id via taxonomy summary
    taxids = list(set(acc_to_taxid.values()))
    taxid_to_parents = {}  # tax_id -> [parent_tax_ids from root]

    try:
        cmd = ['datasets', 'summary', 'taxonomy', 'taxon'] + \
              [str(t) for t in taxids]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            print(f"[phylo_sort] datasets taxonomy error: {result.stderr[:200]}")
            # Fall back to just species names
            for acc, name in acc_to_name.items():
                if name:
                    lineages[acc] = [name]
            return lineages

        data = json.loads(result.stdout)
        for report in data.get('reports', []):
            tax = report.get('taxonomy', {})
            tax_id = tax.get('tax_id', 0)
            parents = tax.get('parents', [])

            if tax_id and parents:
                # parents list goes from root (1) to immediate parent
                # Append the species tax_id itself
                full_lineage = parents + [tax_id]
                taxid_to_parents[tax_id] = full_lineage

    except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError) as e:
        print(f"[phylo_sort] Error in taxonomy summary: {e}")

    # Map back to accessions
    for acc, tax_id in acc_to_taxid.items():
        name = acc_to_name.get(acc, '?')
        if tax_id in taxid_to_parents:
            lineages[acc] = taxid_to_parents[tax_id]
            print(f"[phylo_sort] {acc} -> {name} "
                  f"(lineage depth: {len(lineages[acc])})")
        elif name:
            lineages[acc] = [name]
            print(f"[phylo_sort] {acc} -> {name} (name only, no lineage)")

    return lineages


def _lineage_distance(lin1, lin2):
    """
    Compute a simple topological distance between two lineage lists
    (lists of tax IDs or names, ordered root → species).
    Distance = number of unshared nodes (symmetric difference).
    Lower values indicate more closely related organisms.
    """
    if not lin1 or not lin2:
        return float('inf')

    # Convert to sets for overlap, then also check prefix
    set1 = set(lin1)
    set2 = set(lin2)
    shared = len(set1 & set2)

    # Distance is total unique nodes minus shared (= symmetric difference)
    return (len(set1) + len(set2)) - 2 * shared


def _fetch_species_lineage(species_name):
    """
    Fetch taxonomy lineage for a species name using NCBI Datasets CLI.
    Returns list of parent tax IDs (root → species) or None.
    """
    if not species_name:
        return None
    # Clean: strip path components, extensions
    clean = os.path.basename(species_name)
    clean = re.sub(r'\.(fna|fasta|fa)(\.gz)?$', '', clean)
    clean = clean.replace('_', ' ').strip()
    if not clean or clean.lower() in ('home genome', 'home'):
        return None

    try:
        cmd = ['datasets', 'summary', 'taxonomy', 'taxon', clean]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return None

        data = json.loads(result.stdout)
        for report in data.get('reports', []):
            tax = report.get('taxonomy', {})
            tax_id = tax.get('tax_id', 0)
            parents = tax.get('parents', [])
            if tax_id and parents:
                return parents + [tax_id]

    except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError):
        pass

    return None


def sort_by_accession_taxonomy(home_identifier, target_list):
    """
    Sort targets by phylogenetic distance using NCBI Datasets CLI to
    resolve taxonomy from assembly accessions.

    Returns list of (distance, basename) or None if method fails.
    """
    if not _datasets_cli_available():
        print("[phylo_sort] NCBI datasets CLI not available, skipping accession-based sort")
        return None

    # Collect target accessions
    acc_map = {}  # accession -> filename(s)
    for f in target_list:
        acc = _extract_accession(f)
        if acc:
            acc_map.setdefault(acc, []).append(f)

    # Also try to extract home accession
    home_acc = _extract_accession(home_identifier)
    if home_acc:
        acc_map.setdefault(home_acc, []).append(home_identifier)

    if not acc_map:
        print("[phylo_sort] No GCF_/GCA_ accessions found in filenames")
        return None

    # Fetch lineages for all accessions
    accessions = list(acc_map.keys())
    print(f"[phylo_sort] Fetching taxonomy for {len(accessions)} accessions...")
    lineages = _fetch_lineage_from_accessions(accessions)

    if not lineages:
        print("[phylo_sort] No lineage data retrieved")
        return None

    # Determine home lineage
    home_lineage = None
    if home_acc and home_acc in lineages:
        home_lineage = lineages[home_acc]
        print(f"[phylo_sort] Home lineage from accession {home_acc} "
              f"(depth: {len(home_lineage)})")
    else:
        # Home identifier might be a species name — look it up via taxonomy API
        home_lineage = _fetch_species_lineage(home_identifier)
        if home_lineage:
            print(f"[phylo_sort] Home lineage from species name '{home_identifier}' "
                  f"(depth: {len(home_lineage)})")

    if not home_lineage:
        print(f"[phylo_sort] Could not determine home lineage for '{home_identifier}'")
        # Use the most detailed lineage as reference (imperfect but better than nothing)
        if lineages:
            home_lineage = max(lineages.values(), key=len)
            print(f"[phylo_sort] Using longest lineage as fallback reference")
        else:
            return None

    # Score each target
    scored = []
    for target in target_list:
        acc = _extract_accession(target)
        basename = os.path.basename(target)
        if acc and acc in lineages:
            dist = _lineage_distance(home_lineage, lineages[acc])
            scored.append((dist, basename))
        else:
            # Unknown taxonomy → push to end
            scored.append((float('inf'), basename))

    # Ensure we got real distances for at least some targets
    real_distances = [d for d, _ in scored if d < float('inf')]
    if not real_distances:
        print("[phylo_sort] No taxonomic distances computed")
        return None

    scored.sort(key=lambda x: (x[0], x[1]))
    return scored


# =============================================================================
# ETE3-BASED TAXONOMY SORTING (original method)
# =============================================================================

def get_taxid(name, ncbi):
    """Resolve a name/accession to an NCBI TaxID via ETE3."""
    clean = os.path.basename(name).split('.')[0].replace('_', ' ')

    if clean.isdigit():
        return int(clean)

    try:
        name2taxid = ncbi.get_name_translator([clean])
        if clean in name2taxid:
            return name2taxid[clean][0]
    except Exception:
        pass

    # Try progressively shorter prefixes
    parts = clean.split()
    for i in range(len(parts), 0, -1):
        subname = " ".join(parts[:i])
        try:
            name2taxid = ncbi.get_name_translator([subname])
            if subname in name2taxid:
                return name2taxid[subname][0]
        except Exception:
            continue

    return None


def sort_by_ete3(args, target_list):
    """Sort using ETE3 NCBITaxa. Returns list of (distance, basename) or None."""
    if NCBITaxa is None:
        return None

    print(f"[phylo_sort] Initializing ETE3 Taxonomy from {args.taxdb}...")

    try:
        if os.path.isfile(args.taxdb):
            ncbi = NCBITaxa(dbfile=args.taxdb)
        elif os.path.isdir(args.taxdb):
            db_path = os.path.join(args.taxdb, 'taxa.sqlite')
            if os.path.exists(db_path):
                ncbi = NCBITaxa(dbfile=db_path)
            else:
                ncbi = NCBITaxa()
        else:
            ncbi = NCBITaxa()
    except Exception as e:
        print(f"[phylo_sort] ETE3 init failed: {e}")
        return None

    home_taxid = get_taxid(args.home, ncbi)
    if not home_taxid:
        print(f"[phylo_sort] Could not find TaxID for home: {args.home}")
        return None

    print(f"[phylo_sort] Home TaxID: {home_taxid}")

    try:
        home_lineage = ncbi.get_lineage(home_taxid)
    except ValueError:
        print(f"[phylo_sort] TaxID {home_taxid} not found in DB")
        return None

    scored = []
    for target in target_list:
        tid = get_taxid(target, ncbi)
        if not tid:
            scored.append((float('inf'), os.path.basename(target)))
            continue

        try:
            target_lineage = ncbi.get_lineage(tid)
            shared = 0
            for a, b in zip(home_lineage, target_lineage):
                if a == b:
                    shared += 1
                else:
                    break
            dist = (len(home_lineage) - shared) + (len(target_lineage) - shared)
            scored.append((dist, os.path.basename(target)))
        except Exception as e:
            print(f"[phylo_sort] Error for {target}: {e}")
            scored.append((float('inf'), os.path.basename(target)))

    real_distances = [d for d, _ in scored if d < float('inf')]
    if not real_distances:
        return None

    scored.sort(key=lambda x: (x[0], x[1]))
    return scored


# =============================================================================
# FALLBACK SORT
# =============================================================================

def parse_fasta_simple(filepath):
    import gzip
    _open = gzip.open if str(filepath).endswith('.gz') else open
    current_header = None
    current_seq = []
    with _open(filepath, 'rt') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if current_header is not None:
                    yield ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            elif line and current_header is not None:
                current_seq.append(line.strip())
        if current_header is not None:
            yield ''.join(current_seq)

def minhash_sketch(filepath, k=21, sketch_size=10000):
    import heapq
    sketch = []
    try:
        for seq in parse_fasta_simple(filepath):
            seq = seq.upper()
            stride = max(1, len(seq) // 500000)
            for i in range(0, len(seq) - k + 1, stride):
                kmer = seq[i:i+k]
                if 'N' in kmer: continue
                h = int(hashlib.md5(kmer.encode()).hexdigest(), 16)
                if len(sketch) < sketch_size:
                    heapq.heappush(sketch, -h)
                elif -h > sketch[0]:
                    heapq.heappushpop(sketch, -h)
    except Exception as e:
        pass
    return set([-x for x in sketch])

def jaccard_distance(s1, s2):
    if not s1 or not s2: return 1.0
    union = len(s1 | s2)
    return 1.0 - (len(s1 & s2) / union) if union > 0 else 1.0

def _collect_targets(args):
    target_list = []
    if args.targets:
        target_list = args.targets
    elif args.targets_dir:
        if os.path.isdir(args.targets_dir):
            for f in os.listdir(args.targets_dir):
                if (f.endswith(args.img_ext) or f.endswith(".fasta") or
                    f.endswith(".fa") or f.endswith(args.img_ext + ".gz") or
                    f.endswith(".fasta.gz") or f.endswith(".fa.gz")):
                    target_list.append(os.path.join(args.targets_dir, f))
    return target_list

def fallback_sort(args, target_list):
    """Fallback sort using MinHash Jaccard distance or alphabetical if fasta missing."""
    if not hasattr(args, 'home_fasta') or not args.home_fasta or not os.path.exists(args.home_fasta):
        print("[phylo_sort] Falling back to alphabetical sorting (no home_fasta)")
        with open(args.output, 'w') as out:
            for target in sorted([os.path.basename(t) for t in target_list]):
                out.write(f"{target}\t0\n")
        return

    print(f"[phylo_sort] Falling back to MinHash sorting against {args.home_fasta}")
    home_sketch = minhash_sketch(args.home_fasta)
    results = []
    for tgt in target_list:
        tgt_path = tgt if os.path.exists(tgt) else os.path.join(args.targets_dir or "", tgt)
        if not os.path.exists(tgt_path):
            results.append((1.0, os.path.basename(tgt)))
            continue
        tgt_sketch = minhash_sketch(tgt_path)
        dist = jaccard_distance(home_sketch, tgt_sketch)
        # Store scaled distance so sorting and integers work, max 1000
        dist_scaled = dist * 1000.0
        results.append((dist_scaled, os.path.basename(tgt)))
    
    results.sort(key=lambda x: x[0])
    with open(args.output, 'w') as f:
        for dist, target in results:
            f.write(f"{target}\t{int(dist)}\n")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Sort genomes by phylogenetic distance"
    )
    parser.add_argument("--home", required=True,
                        help="Home genome name, species name, or TaxID")
    parser.add_argument("--home_fasta", default="",
                        help="Path to home genome FASTA for MinHash fallback")
    parser.add_argument("--targets", nargs='+',
                        help="List of target genome files or names")
    parser.add_argument("--targets_dir",
                        help="Directory containing target genomes")
    parser.add_argument("--img_ext", default=".fna",
                        help="Extension for targets in directory (default: .fna)")
    parser.add_argument("--taxdb", required=True,
                        help="Path to taxonomy DB (directory or sqlite file)")
    parser.add_argument("--output", required=True,
                        help="Output sorted list of genomes")

    args = parser.parse_args()

    # Collect target files
    target_list = []
    if args.targets:
        target_list = args.targets
    elif args.targets_dir:
        if os.path.isdir(args.targets_dir):
            for f in os.listdir(args.targets_dir):
                if (f.endswith(args.img_ext) or f.endswith(".fasta") or
                    f.endswith(".fa") or f.endswith(args.img_ext + ".gz") or
                    f.endswith(".fasta.gz") or f.endswith(".fa.gz")):
                    target_list.append(os.path.join(args.targets_dir, f))

    if not target_list:
        print("[phylo_sort] Warning: No targets found")
        with open(args.output, 'w') as f:
            pass
        return

    # Strategy 1: ETE3 + local taxonomy DB
    result = sort_by_ete3(args, target_list)
    if result:
        print(f"[phylo_sort] Sorted by ETE3 taxonomy ({len(result)} genomes)")
        _write_result(args.output, result)
        return

    # Strategy 2: NCBI Datasets CLI (accession-based)
    home_identifier = args.home if args.home else ''
    # If home is a file path, also pass it for accession extraction
    if not home_identifier and hasattr(args, 'home'):
        home_identifier = args.home
    result = sort_by_accession_taxonomy(home_identifier, target_list)
    if result:
        print(f"[phylo_sort] Sorted by accession taxonomy ({len(result)} genomes)")
        _write_result(args.output, result)
        return

    # Strategy 3: MinHash / Alphabetical fallback
    fallback_sort(args, target_list)


def _write_result(output_file, scored_targets):
    """Write sorted results to output file."""
    with open(output_file, 'w') as f:
        for dist, target in scored_targets:
            dist_val = int(dist) if dist != float('inf') else 0
            f.write(f"{target}\t{dist_val}\n")


if __name__ == "__main__":
    main()
