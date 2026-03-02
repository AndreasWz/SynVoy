#!/usr/bin/env python3
"""
Fetch the reference genome for a species from NCBI.
Downloads both the genome FASTA and GFF annotation.
"""

import argparse
import json
import select
import subprocess
import sys
import time
import zipfile
from pathlib import Path
import shutil

LARGE_RANK = 10**18


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


def extract_zip_archive(zip_path: Path, extract_dir: Path):
    """
    Extract a ZIP archive using Python stdlib so Docker/Conda runs do not
    depend on external `unzip` being installed.
    """
    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)


def parse_int(value):
    if value is None:
        return None
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)):
        return int(value)
    if isinstance(value, str):
        txt = value.strip().replace(",", "")
        if not txt or txt.lower() in {"na", "n/a", "none", "null"}:
            return None
        try:
            return int(float(txt))
        except ValueError:
            return None
    return None


def parse_float(value):
    if value is None:
        return None
    if isinstance(value, bool):
        return float(int(value))
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        txt = value.strip().replace(",", "")
        if not txt or txt.lower() in {"na", "n/a", "none", "null"}:
            return None
        try:
            return float(txt)
        except ValueError:
            return None
    return None


def refseq_priority(entry):
    cat = (entry.get("category") or "").strip().lower()
    acc = entry.get("accession", "")
    if "reference genome" in cat:
        return 0
    if "representative genome" in cat:
        return 1
    if acc.startswith("GCF_"):
        return 2
    if acc.startswith("GCA_"):
        return 3
    return 4


def assembly_level_priority(level):
    txt = (level or "").strip().lower()
    if "chromosome" in txt:
        return 0
    if "complete genome" in txt or txt == "complete":
        return 1
    if "scaffold" in txt:
        return 2
    if "contig" in txt:
        return 3
    return 4


def asc_rank(value):
    if value is None:
        return LARGE_RANK
    return value


def desc_rank(value):
    if value is None:
        return LARGE_RANK
    return -value


def assembly_rank_tuple(entry, ranking_mode):
    ref_rank = refseq_priority(entry)
    level_rank = assembly_level_priority(entry.get("assembly_status"))
    chrom = entry.get("chromosome_count")
    scaff = entry.get("scaffold_count")
    contigs = entry.get("contig_count")
    scaf_n50 = entry.get("scaffold_n50")
    cont_n50 = entry.get("contig_n50")
    scaf_n80 = entry.get("scaffold_n80")
    cont_n80 = entry.get("contig_n80")
    acc = entry.get("accession", "")

    if ranking_mode == "counts":
        return (
            ref_rank,
            level_rank,
            asc_rank(chrom),
            asc_rank(scaff),
            asc_rank(contigs),
            desc_rank(scaf_n50),
            desc_rank(cont_n50),
            desc_rank(scaf_n80),
            desc_rank(cont_n80),
            acc,
        )
    if ranking_mode == "nstats":
        return (
            ref_rank,
            level_rank,
            desc_rank(scaf_n80),
            desc_rank(cont_n80),
            desc_rank(scaf_n50),
            desc_rank(cont_n50),
            asc_rank(contigs),
            asc_rank(scaff),
            asc_rank(chrom),
            acc,
        )
    return (
        ref_rank,
        level_rank,
        asc_rank(contigs),
        asc_rank(scaff),
        desc_rank(cont_n50),
        desc_rank(scaf_n50),
        desc_rank(cont_n80),
        desc_rank(scaf_n80),
        asc_rank(chrom),
        acc,
    )


def get_assembly_quality(accession):
    """Fetch quality fields for one assembly accession from NCBI assembly docsum."""
    query = accession
    cmds_search = ["esearch", "-db", "assembly", "-query", query]
    cmds_fetch = ["efetch", "-format", "docsum"]
    cmds_extract = [
        "xtract",
        "-pattern",
        "DocumentSummary",
        "-element",
        "AssemblyAccession",
        "SpeciesName",
        "RefSeq_category",
        "AssemblyStatus",
        "ScaffoldCount",
        "ContigCount",
        "ScaffoldN50",
        "ContigN50",
        "ScaffoldN80",
        "ContigN80",
    ]
    try:
        p1 = subprocess.Popen(cmds_search, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(cmds_fetch, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(cmds_extract, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2.stdout.close()
        output, _ = p3.communicate(timeout=30)
        line = output.decode().strip().split("\n")[0] if output else ""
        if not line:
            return {}
        parts = line.split("\t")
        parts += [""] * (10 - len(parts))
        return {
            "accession": parts[0].strip(),
            "species": parts[1].strip(),
            "category": parts[2].strip() if len(parts) > 2 else None,
            "assembly_status": parts[3].strip() if len(parts) > 3 else None,
            "scaffold_count": parse_int(parts[4]) if len(parts) > 4 else None,
            "contig_count": parse_int(parts[5]) if len(parts) > 5 else None,
            "scaffold_n50": parse_float(parts[6]) if len(parts) > 6 else None,
            "contig_n50": parse_float(parts[7]) if len(parts) > 7 else None,
            "scaffold_n80": parse_float(parts[8]) if len(parts) > 8 else None,
            "contig_n80": parse_float(parts[9]) if len(parts) > 9 else None,
        }
    except Exception:
        return {}


def format_quality(entry):
    return (
        f"level={entry.get('assembly_status') or 'NA'}, "
        f"chr={entry.get('chromosome_count') or 'NA'}, "
        f"scaf={entry.get('scaffold_count') or 'NA'}, "
        f"contigs={entry.get('contig_count') or 'NA'}, "
        f"N50(contig/scaf)={entry.get('contig_n50') or 'NA'}/{entry.get('scaffold_n50') or 'NA'}"
    )


def is_bad_quality(entry, args):
    # Chromosome-level and complete-genome assemblies are always acceptable.
    # NCBI DocSum counts include alternate haplotypes and unplaced sequences,
    # which inflate contig/scaffold numbers even for excellent assemblies.
    level_rank = assembly_level_priority(entry.get("assembly_status"))
    if level_rank <= 1:  # Chromosome or Complete Genome
        return False, []

    reasons = []
    contigs = entry.get("contig_count")
    scaff = entry.get("scaffold_count")
    contig_n50 = entry.get("contig_n50")
    scaffold_n50 = entry.get("scaffold_n50")
    best_n50 = max(v for v in (contig_n50, scaffold_n50) if v is not None) if (contig_n50 is not None or scaffold_n50 is not None) else None

    if contigs is not None and contigs > args.bad_max_contigs:
        reasons.append(f"contig_count={contigs} > {args.bad_max_contigs}")
    if scaff is not None and scaff > args.bad_max_scaffolds:
        reasons.append(f"scaffold_count={scaff} > {args.bad_max_scaffolds}")
    if best_n50 is not None and best_n50 < args.bad_min_n50:
        reasons.append(f"best_N50={int(best_n50)} < {args.bad_min_n50}")
    if (
        best_n50 is None
        and contigs is None
        and scaff is None
        and level_rank >= 3
    ):
        reasons.append("assembly level is contig-like and no quality metrics are available")

    return len(reasons) > 0, reasons


def ask_keep_bad_quality(entry, reasons, timeout_seconds):
    msg = [
        "",
        f"WARNING: Home genome candidate {entry.get('accession', 'unknown')} looks low quality.",
        f"  Details: {format_quality(entry)}",
        "  Reasons: " + "; ".join(reasons),
    ]
    question = f"  Keep this genome? [y/N] (auto-NO after {timeout_seconds}s): "

    try:
        with open("/dev/tty", "r", encoding="utf-8", errors="ignore") as tty_in, open(
            "/dev/tty", "w", encoding="utf-8", errors="ignore"
        ) as tty_out:
            tty_out.write("\n".join(msg) + "\n")
            tty_out.write(question)
            tty_out.flush()
            ready, _, _ = select.select([tty_in], [], [], timeout_seconds)
            if ready:
                answer = tty_in.readline().strip().lower()
                return answer in {"y", "yes"}
            tty_out.write("\n  No response. Proceeding with NO.\n")
            tty_out.flush()
            return False
    except Exception:
        if sys.stdin is not None and sys.stdin.isatty():
            print("\n".join(msg))
            print(question, end="", flush=True)
            ready, _, _ = select.select([sys.stdin], [], [], timeout_seconds)
            if ready:
                answer = sys.stdin.readline().strip().lower()
                return answer in {"y", "yes"}
            print("\n  No response. Proceeding with NO.")
            return False

    print("\n".join(msg))
    print(f"  No interactive terminal available. Waiting {timeout_seconds}s, then proceeding with NO.")
    if timeout_seconds > 0:
        time.sleep(timeout_seconds)
    return False


def get_reference_genome(species_name):
    """Get the reference/representative genome accession for a species.
    Returns (accession, has_annotation, actual_species) tuple.
    """
    print(f"Finding reference genome for '{species_name}'...")
    
    cmd = [
        'datasets', 'summary', 'genome', 'taxon', species_name,
        '--reference', '--as-json-lines'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        for line in result.stdout.strip().split('\n'):
            if line:
                data = json.loads(line)
                
                if 'accession' in data:
                    accession = data.get('accession')
                    org_name = data.get('organism', {}).get('organism_name', species_name)
                    annotation = data.get('annotation_info', {})
                    has_annotation = annotation.get('status') == 'Full annotation' if annotation else False
                    
                    print(f"  Found reference: {accession} ({org_name})")
                    return accession, has_annotation, org_name
                
                elif 'reports' in data:
                    for report in data['reports']:
                        accession = report.get('accession')
                        org_name = report.get('organism', {}).get('organism_name', species_name)
                        annotation = report.get('annotation_info', {})
                        has_annotation = annotation.get('status') == 'Full annotation' if annotation else False
                        
                        print(f"  Found reference: {accession} ({org_name})")
                        return accession, has_annotation, org_name
                        
    except subprocess.CalledProcessError as e:
        print(f"  No reference genome found via datasets: {e.stderr.strip()}", file=sys.stderr)
    except Exception as e:
        print(f"  Error: {e}", file=sys.stderr)
    
    return None, False, species_name


def find_any_genome(species_name, ranking_mode="hybrid"):
    """Find any genome assembly (RefSeq or GenBank) for a species.
    Falls back from reference/representative and ranks by assembly quality.
    Returns (entry_dict, has_annotation, actual_species) tuple.
    """
    print(f"  Searching for any genome assembly for '{species_name}'...")
    
    # Use esearch + efetch to find all assemblies
    try:
        query = f'"{species_name}"[Organism]'
        cmds_search = ['esearch', '-db', 'assembly', '-query', query]
        cmds_fetch = ['efetch', '-format', 'docsum']
        cmds_extract = [
            'xtract', '-pattern', 'DocumentSummary', '-element',
            'AssemblyAccession', 'SpeciesName', 'RefSeq_category',
            'AssemblyStatus', 'ScaffoldCount', 'ContigCount',
            'ScaffoldN50', 'ContigN50', 'ScaffoldN80', 'ContigN80'
        ]
        
        p1 = subprocess.Popen(cmds_search, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(cmds_fetch, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(cmds_extract, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2.stdout.close()
        output, _ = p3.communicate(timeout=30)
        
        results = output.decode().strip()
        if not results:
            return None, False, species_name
        
        candidates = []
        for line in results.split('\n'):
            parts = line.split('\t')
            if len(parts) >= 2:
                parts += [''] * (10 - len(parts))
                entry = {
                    'accession': parts[0].strip(),
                    'species': parts[1].strip(),
                    'category': parts[2].strip() if len(parts) > 2 else None,
                    'assembly_status': parts[3].strip() if len(parts) > 3 else None,
                    'scaffold_count': parse_int(parts[4]) if len(parts) > 4 else None,
                    'contig_count': parse_int(parts[5]) if len(parts) > 5 else None,
                    'scaffold_n50': parse_float(parts[6]) if len(parts) > 6 else None,
                    'contig_n50': parse_float(parts[7]) if len(parts) > 7 else None,
                    'scaffold_n80': parse_float(parts[8]) if len(parts) > 8 else None,
                    'contig_n80': parse_float(parts[9]) if len(parts) > 9 else None,
                }
                candidates.append(entry)

        if candidates:
            candidates.sort(key=lambda x: assembly_rank_tuple(x, ranking_mode))
            best = candidates[0]
            src = 'RefSeq' if best['accession'].startswith('GCF_') else 'GenBank'
            print(f"  Found {src} assembly: {best['accession']} ({best['species']})")
            print(f"  Quality: {format_quality(best)}")
            return best, True, best['species']  # annotation may exist
        
    except Exception as e:
        print(f"  Error searching assemblies: {e}", file=sys.stderr)
    
    return None, False, species_name


def find_closest_relative_genome(species_name):
    """Walk up taxonomy to find closest relative with a genome.
    Returns (accession, has_annotation, actual_species) tuple.
    """
    print(f"  No genome for '{species_name}'. Searching closest relatives...")
    
    # Get taxonomy ID for species
    try:
        p1 = subprocess.Popen(
            ['esearch', '-db', 'taxonomy', '-query', species_name],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p2 = subprocess.Popen(
            ['efetch', '-format', 'uid'],
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p1.stdout.close()
        output, _ = p2.communicate(timeout=15)
        taxid = output.decode().strip()
    except Exception:
        return None, False, species_name
    
    if not taxid:
        return None, False, species_name
    
    # Get lineage
    try:
        p1 = subprocess.Popen(
            ['efetch', '-db', 'taxonomy', '-id', taxid, '-format', 'xml'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p2 = subprocess.Popen(
            ['xtract', '-pattern', 'Taxon', '-block', 'LineageEx/Taxon',
             '-element', 'Rank', '-element', 'TaxId', '-element', 'ScientificName'],
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p1.stdout.close()
        output, _ = p2.communicate(timeout=15)
        raw = output.decode().strip()
    except Exception:
        return None, False, species_name
    
    if not raw:
        return None, False, species_name
    
    # Parse lineage into (rank, taxid, name) triples
    parts = raw.split('\t')
    triples = []
    for i in range(0, len(parts) - 2, 3):
        triples.append((parts[i].strip(), parts[i+1].strip(), parts[i+2].strip()))
    
    # Try genus first, then family, then order
    target_ranks = ['genus', 'family', 'order']
    for target_rank in target_ranks:
        for rank, tid, name in triples:
            if rank.lower() == target_rank:
                print(f"  Trying {rank} '{name}' (TaxID: {tid})...")
                
                # Use datasets to find reference genome in this taxon
                cmd = [
                    'datasets', 'summary', 'genome', 'taxon', tid,
                    '--reference', '--as-json-lines'
                ]
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=30)
                    for line in result.stdout.strip().split('\n'):
                        if line:
                            data = json.loads(line)
                            if 'accession' in data:
                                acc = data['accession']
                                org = data.get('organism', {}).get('organism_name', name)
                                print(f"  ✓ Using closest relative: {acc} ({org})")
                                return acc, True, org
                            elif 'reports' in data:
                                for report in data['reports']:
                                    acc = report['accession']
                                    org = report.get('organism', {}).get('organism_name', name)
                                    print(f"  ✓ Using closest relative: {acc} ({org})")
                                    return acc, True, org
                except Exception:
                    continue
    
    return None, False, species_name


def filter_chromosomes_only(fna_path: Path) -> None:
    """
    If a FASTA contains both chromosome sequences (NC_/CM_ accession prefixes)
    and unlocalized/placed scaffolds (NW_, NZ_, or other prefixes), rewrite the
    file keeping only the chromosome sequences.

    If NO chromosome sequences are found (pure scaffold assembly), the file is
    left untouched so we do not discard all sequence data.
    """
    CHROM_PREFIXES = ("NC_", "CM")

    chrom_ids = []
    non_chrom_ids = []
    try:
        with open(fna_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                seq_id = line[1:].split()[0] if line.strip() else ""
                if any(seq_id.startswith(p) for p in CHROM_PREFIXES):
                    chrom_ids.append(seq_id)
                else:
                    non_chrom_ids.append(seq_id)
    except Exception as e:
        print(f"  [chr-filter] Warning: could not scan {fna_path}: {e}")
        return

    if not chrom_ids:
        print(
            f"  [chr-filter] No NC_/CM_ chromosomes found — retaining all "
            f"{len(non_chrom_ids)} sequences as-is (scaffold assembly)."
        )
        return

    if not non_chrom_ids:
        print(
            f"  [chr-filter] {len(chrom_ids)} chromosome sequence(s) — no scaffolds to remove."
        )
        return

    keep_set = set(chrom_ids)
    tmp_path = Path(str(fna_path) + ".chr_only.tmp")
    kept = 0
    skipping = False
    try:
        with open(fna_path) as src, open(tmp_path, "w") as dst:
            for line in src:
                if line.startswith(">"):
                    seq_id = line[1:].split()[0] if line.strip() else ""
                    skipping = seq_id not in keep_set
                    if not skipping:
                        kept += 1
                if not skipping:
                    dst.write(line)
        tmp_path.replace(fna_path)
        print(
            f"  [chr-filter] Kept {kept} chromosome sequence(s), "
            f"removed {len(non_chrom_ids)} scaffold/unlocalized sequence(s)."
        )
    except Exception as e:
        print(f"  [chr-filter] Warning: filter failed for {fna_path}: {e}")
        if tmp_path.exists():
            tmp_path.unlink()


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
        
        # Extract archive (stdlib; avoids external unzip dependency)
        extract_dir = output_path / 'extracted'
        extract_zip_archive(zip_file, extract_dir)
        
        # Find files
        fna_files = list(extract_dir.rglob("*.fna"))
        gff_files = list(extract_dir.rglob("*.gff"))
        
        genome_path = None
        gff_path = None
        
        if fna_files:
            genome_path = output_path / "home_genome.fna"
            shutil.copy(fna_files[0], genome_path)
            # Filter to chromosome sequences only when the assembly contains
            # both chromosomes (NC_/CM_) and unlocalized scaffolds (NW_/NZ_).
            filter_chromosomes_only(genome_path)
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
    parser.add_argument(
        "--assembly-ranking",
        choices=["hybrid", "counts", "nstats"],
        default="hybrid",
        help="Ranking strategy for non-reference fallback selection (default: hybrid)",
    )
    parser.add_argument(
        "--bad-quality-policy",
        choices=["ask", "drop", "keep"],
        default="ask",
        help="What to do if only low-quality home assembly is available (default: ask)",
    )
    parser.add_argument(
        "--bad-quality-timeout",
        type=int,
        default=300,
        help="Prompt timeout in seconds when bad-quality-policy=ask (default: 300)",
    )
    parser.add_argument(
        "--bad-max-contigs",
        type=int,
        default=100000,
        help="Assemblies above this contig count are flagged low quality (default: 100000)",
    )
    parser.add_argument(
        "--bad-max-scaffolds",
        type=int,
        default=50000,
        help="Assemblies above this scaffold count are flagged low quality (default: 50000)",
    )
    parser.add_argument(
        "--bad-min-n50",
        type=int,
        default=20000,
        help="Assemblies with best N50 below this are flagged low quality (default: 20000)",
    )
    
    args = parser.parse_args()
    
    # Check for datasets tool
    try:
        subprocess.run(['which', 'datasets'], check=True, capture_output=True)
    except:
        print("ERROR: NCBI datasets tool not found!", file=sys.stderr)
        print("Install with: conda install -c conda-forge ncbi-datasets-cli", file=sys.stderr)
        sys.exit(1)
    
    # 3-tier fallback:
    # 1. Reference/representative genome (highest quality)
    # 2. Any genome for this species (RefSeq > GenBank)
    # 3. Closest relative's genome (walk up taxonomy)
    
    accession, has_annotation, actual_species = get_reference_genome(args.species)
    selected_entry = None
    
    if not accession:
        print(f"  No reference genome. Trying any assembly...")
        selected_entry, has_annotation, actual_species = find_any_genome(
            args.species, ranking_mode=args.assembly_ranking
        )
        accession = selected_entry.get("accession") if selected_entry else None
    else:
        selected_entry = get_assembly_quality(accession)
        if selected_entry and not selected_entry.get("species"):
            selected_entry["species"] = actual_species
    
    if not accession:
        print(f"  No genome at all. Searching closest relatives...")
        accession, has_annotation, actual_species = find_closest_relative_genome(args.species)
        selected_entry = get_assembly_quality(accession) if accession else None
        if selected_entry and not selected_entry.get("species"):
            selected_entry["species"] = actual_species
    
    if not accession:
        print(f"ERROR: Could not find any genome for '{args.species}' or its relatives", file=sys.stderr)
        sys.exit(1)
    
    surrogate = actual_species.lower().strip() != args.species.lower().strip()
    if surrogate:
        print(f"  WARNING: Using surrogate genome from '{actual_species}' instead of '{args.species}'")

    if selected_entry:
        bad, reasons = is_bad_quality(selected_entry, args)
        if bad:
            action = args.bad_quality_policy
            if action == "keep":
                print(
                    "  [keep] Using low-quality home assembly by policy: "
                    + "; ".join(reasons)
                )
            elif action == "drop":
                print(
                    "ERROR: Home assembly rejected by low-quality policy: "
                    + "; ".join(reasons),
                    file=sys.stderr,
                )
                sys.exit(1)
            else:
                keep = ask_keep_bad_quality(selected_entry, reasons, args.bad_quality_timeout)
                if not keep:
                    print(
                        "ERROR: Home assembly rejected by user/timeout due to low quality.",
                        file=sys.stderr,
                    )
                    sys.exit(1)
    
    # Download
    genome_path, gff_path = download_genome_with_annotation(accession, args.outdir)
    
    if not genome_path:
        print("ERROR: Failed to download genome", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n{'='*60}")
    if surrogate:
        print(f"SUCCESS: Downloaded surrogate genome ({actual_species})")
        print(f"  Requested species: {args.species}")
        print(f"  Surrogate species: {actual_species}")
    else:
        print(f"SUCCESS: Downloaded genome for {args.species}")
    print(f"  Accession: {accession}")
    print(f"  Genome:    {genome_path}")
    if gff_path:
        print(f"  GFF:       {gff_path}")
    print(f"{'='*60}")
    
    # Write paths to files for Nextflow to pick up
    (Path(args.outdir) / "genome_path.txt").write_text(str(genome_path))
    if gff_path:
        (Path(args.outdir) / "gff_path.txt").write_text(str(gff_path))
    # Write actual species name (may differ if surrogate)
    (Path(args.outdir) / "actual_species.txt").write_text(actual_species)


if __name__ == "__main__":
    main()
