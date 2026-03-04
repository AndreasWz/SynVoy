#!/usr/bin/env python3
"""
Fetch related genomes from NCBI for easy mode.
Uses NCBI E-utilities + Datasets and ranks assemblies by quality.
"""

import argparse
import gzip
import json
import os
import select
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path


LARGE_RANK = 10**18


def run_safe_command(cmd, check=True):
    """Run command with list arguments safely."""
    try:
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Error: {e.stderr}", file=sys.stderr)
        if check:
            raise e
        return None


def extract_zip_archive(zip_path: Path, extract_dir: Path):
    """
    Extract datasets ZIP archives without requiring external `unzip`.
    """
    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)


def run_piped_command(cmds):
    """
    Run a chain of commands connected by pipes.
    cmds: List of command lists. E.g. [['esearch', ...], ['efetch', ...]]
    """
    procs = []
    try:
        p1 = subprocess.Popen(cmds[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        procs.append(p1)

        prev_stdout = p1.stdout
        for i in range(1, len(cmds)):
            p_next = subprocess.Popen(
                cmds[i], stdin=prev_stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            procs.append(p_next)
            prev_stdout.close()
            prev_stdout = p_next.stdout

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
        for p in procs:
            if p.poll() is None:
                p.kill()


def normalize_species(name):
    return (name or "").strip().lower()


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
    if isinstance(value, list):
        return len(value)
    if isinstance(value, dict):
        for k in ("count", "value", "number", "total"):
            if k in value:
                parsed = parse_int(value[k])
                if parsed is not None:
                    return parsed
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
    if isinstance(value, dict):
        for k in ("count", "value", "number", "total"):
            if k in value:
                parsed = parse_float(value[k])
                if parsed is not None:
                    return parsed
    return None


def normalize_key(key):
    return "".join(ch for ch in str(key).lower() if ch.isalnum())


def extract_metric_from_json(data, aliases, parser):
    alias_set = {normalize_key(a) for a in aliases}
    stack = [data]
    while stack:
        node = stack.pop()
        if isinstance(node, dict):
            for key, value in node.items():
                if normalize_key(key) in alias_set:
                    parsed = parser(value)
                    if parsed is not None:
                        return parsed
                if isinstance(value, (dict, list)):
                    stack.append(value)
        elif isinstance(node, list):
            stack.extend(node)
    return None


def load_datasets_quality(accession, metadata_cache):
    if accession in metadata_cache:
        return metadata_cache[accession]

    cmd = [
        "datasets",
        "summary",
        "genome",
        "accession",
        accession,
        "--as-json-lines",
    ]
    out = None
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=45)
        out = result.stdout
    except Exception:
        metadata_cache[accession] = {}
        return {}

    record = None
    for raw in out.strip().splitlines():
        raw = raw.strip()
        if not raw:
            continue
        try:
            payload = json.loads(raw)
        except json.JSONDecodeError:
            continue

        if isinstance(payload, dict) and payload.get("accession") == accession:
            record = payload
            break

        reports = payload.get("reports") if isinstance(payload, dict) else None
        if isinstance(reports, list):
            exact = None
            for rep in reports:
                if isinstance(rep, dict) and rep.get("accession") == accession:
                    exact = rep
                    break
            if exact is not None:
                record = exact
                break
            if record is None and reports:
                first = reports[0]
                if isinstance(first, dict):
                    record = first

        if record is None and isinstance(payload, dict):
            record = payload

    if record is None:
        metadata_cache[accession] = {}
        return {}

    metadata = {
        "assembly_status": extract_metric_from_json(
            record, {"assembly_status", "assembly_level"}, str
        ),
        "category": extract_metric_from_json(record, {"refseq_category"}, str),
        "chromosome_count": extract_metric_from_json(
            record,
            {
                "number_of_chromosomes",
                "chromosome_count",
                "num_chromosomes",
                "total_number_of_chromosomes",
            },
            parse_int,
        ),
        "scaffold_count": extract_metric_from_json(
            record,
            {"number_of_scaffolds", "scaffold_count", "num_scaffolds", "total_scaffold_count"},
            parse_int,
        ),
        "contig_count": extract_metric_from_json(
            record,
            {"number_of_contigs", "contig_count", "num_contigs", "total_contig_count"},
            parse_int,
        ),
        "scaffold_n50": extract_metric_from_json(record, {"scaffold_n50"}, parse_float),
        "contig_n50": extract_metric_from_json(record, {"contig_n50"}, parse_float),
        "scaffold_n80": extract_metric_from_json(record, {"scaffold_n80"}, parse_float),
        "contig_n80": extract_metric_from_json(record, {"contig_n80"}, parse_float),
    }
    metadata_cache[accession] = metadata
    return metadata


def enrich_quality_metadata(entry, metadata_cache):
    needs_enrichment = (
        entry.get("assembly_status") is None
        or (
            entry.get("chromosome_count") is None
            and entry.get("scaffold_count") is None
            and entry.get("contig_count") is None
            and entry.get("scaffold_n50") is None
            and entry.get("contig_n50") is None
        )
    )
    if not needs_enrichment:
        return entry

    metadata = load_datasets_quality(entry["accession"], metadata_cache)
    if not metadata:
        return entry
    for key, value in metadata.items():
        if entry.get(key) is None and value is not None:
            entry[key] = value
    return entry


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
    # Hybrid: prioritize reference/level, then contiguity counts, then N50/N80.
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


def format_quality(entry):
    return (
        f"level={entry.get('assembly_status') or 'NA'}, "
        f"chr={entry.get('chromosome_count') or 'NA'}, "
        f"scaf={entry.get('scaffold_count') or 'NA'}, "
        f"contigs={entry.get('contig_count') or 'NA'}, "
        f"N50(contig/scaf)={entry.get('contig_n50') or 'NA'}/{entry.get('scaffold_n50') or 'NA'}, "
        f"N80(contig/scaf)={entry.get('contig_n80') or 'NA'}/{entry.get('scaffold_n80') or 'NA'}"
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
    if best_n50 is None and contigs is None and scaff is None and level_rank >= 3:
        reasons.append("assembly level is contig-like and no quality metrics are available")

    return len(reasons) > 0, reasons


def ask_keep_bad_quality(entry, reasons, timeout_seconds):
    msg = [
        "",
        f"WARNING: {entry['species']} ({entry['accession']}) looks low quality.",
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


def apply_bad_quality_policy(assemblies, args):
    kept = []
    dropped = []
    for asm in assemblies:
        bad, reasons = is_bad_quality(asm, args)
        asm["bad_quality"] = bad
        asm["bad_quality_reasons"] = "; ".join(reasons)
        if not bad:
            kept.append(asm)
            continue

        if args.bad_quality_policy == "keep":
            print(
                f"  [keep] Retaining low-quality assembly {asm['accession']} for {asm['species']}: "
                + "; ".join(reasons)
            )
            kept.append(asm)
            continue
        if args.bad_quality_policy == "drop":
            print(
                f"  [drop] Excluding low-quality assembly {asm['accession']} for {asm['species']}: "
                + "; ".join(reasons)
            )
            dropped.append(asm)
            continue

        if ask_keep_bad_quality(asm, reasons, args.bad_quality_timeout):
            kept.append(asm)
        else:
            dropped.append(asm)

    if dropped:
        print(f"  Excluded {len(dropped)} low-quality assembly(ies) by policy '{args.bad_quality_policy}'.")
    return kept


def get_taxid_from_name(species_name):
    """Get NCBI taxonomy ID from species name."""
    cmds = [
        ["esearch", "-db", "taxonomy", "-query", species_name],
        ["efetch", "-format", "uid"],
    ]
    return run_piped_command(cmds)


def get_parent_taxa(taxid, include_genus=False):
    """
    Walk up the NCBI taxonomy tree from a given taxid.
    Returns list of (rank_label, taxid) for interesting ranks.
    """
    target_ranks = ["genus", "family", "order", "class"] if include_genus else ["family", "order", "class"]
    results = []
    try:
        cmds = [
            ["efetch", "-db", "taxonomy", "-id", str(taxid), "-format", "xml"],
            [
                "xtract",
                "-pattern",
                "Taxon",
                "-block",
                "LineageEx/Taxon",
                "-element",
                "Rank",
                "-element",
                "TaxId",
                "-element",
                "ScientificName",
            ],
        ]
        raw = run_piped_command(cmds)
        if not raw:
            return results
        parts = raw.split("\t")
        triples = []
        for i in range(0, len(parts) - 2, 3):
            triples.append((parts[i].strip(), parts[i + 1].strip(), parts[i + 2].strip()))
        for rank, tid, name in triples:
            if rank.lower() in target_ranks:
                results.append((f"{rank} ({name})", tid))
        rank_order = {r: i for i, r in enumerate(target_ranks)}
        results.sort(key=lambda x: rank_order.get(x[0].split()[0].lower(), 99))
    except Exception as e:
        print(f"Warning: Could not retrieve lineage for TaxID {taxid}: {e}")
    return results


def get_related_species(
    taxid,
    max_genomes=10,
    exclude_species=None,
    tax_level="",
    ranking_mode="hybrid",
    metadata_cache=None,
):
    """
    Get list of related species using NCBI taxonomy and select best assembly per species.
    """
    if metadata_cache is None:
        metadata_cache = {}

    print(f"  Searching TaxID {taxid} ({tax_level})...")

    query = f"txid{taxid}[Organism:exp]"
    cmds = [
        ["esearch", "-db", "assembly", "-query", query],
        ["efetch", "-format", "docsum"],
        [
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
        ],
    ]
    results = run_piped_command(cmds)
    if not results:
        return []

    exclude_lower = normalize_species(exclude_species) if exclude_species else None
    candidates_by_species = {}
    max_scan_species = max(100, max_genomes * 8)
    max_scan_lines = max(2000, max_genomes * 200)
    line_count = 0

    for line in results.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        parts += [""] * (10 - len(parts))
        accession = parts[0].strip()
        species = parts[1].strip()
        if not accession or not species:
            continue
        if exclude_lower and normalize_species(species) == exclude_lower:
            continue

        entry = {
            "accession": accession,
            "species": species,
            "category": parts[2].strip() if len(parts) > 2 else None,
            "assembly_status": parts[3].strip() if len(parts) > 3 else None,
            "scaffold_count": parse_int(parts[4]) if len(parts) > 4 else None,
            "contig_count": parse_int(parts[5]) if len(parts) > 5 else None,
            "scaffold_n50": parse_float(parts[6]) if len(parts) > 6 else None,
            "contig_n50": parse_float(parts[7]) if len(parts) > 7 else None,
            "scaffold_n80": parse_float(parts[8]) if len(parts) > 8 else None,
            "contig_n80": parse_float(parts[9]) if len(parts) > 9 else None,
            "tax_level": tax_level,
        }
        sp_key = normalize_species(species)
        candidates_by_species.setdefault(sp_key, []).append(entry)

        line_count += 1
        if len(candidates_by_species) >= max_scan_species and line_count >= max_scan_lines:
            break

    best_per_species = []
    for species_entries in candidates_by_species.values():
        species_entries.sort(key=lambda x: assembly_rank_tuple(x, ranking_mode))
        enrich_n = min(len(species_entries), 4)
        for i in range(enrich_n):
            enrich_quality_metadata(species_entries[i], metadata_cache)
        for asm in species_entries:
            asm["_rank"] = assembly_rank_tuple(asm, ranking_mode)
        best = min(species_entries, key=lambda x: x["_rank"])
        best_per_species.append(best)

    best_per_species.sort(key=lambda x: x["_rank"])
    unique_assemblies = best_per_species[:max_genomes]

    n_refseq = len([a for a in unique_assemblies if a["accession"].startswith("GCF_")])
    n_genbank = len([a for a in unique_assemblies if a["accession"].startswith("GCA_")])
    print(f"    Found {len(unique_assemblies)} species ({n_refseq} RefSeq, {n_genbank} GenBank)")

    return unique_assemblies


def _extract_best_gff(extract_dir: Path, gff_target: Path) -> str:
    """
    Extract best available annotation GFF from an extracted datasets directory.

    Returns:
      "ok"          -> valid GFF copied to gff_target
      "no_features" -> GFF files exist but no gene/CDS found in quick scan
      "absent"      -> no GFF-like files found
    """
    gff_files = []
    for pat in ("*.gff", "*.gff3", "*.gff.gz", "*.gff3.gz"):
        gff_files.extend(extract_dir.rglob(pat))
    if not gff_files:
        return "absent"

    best_gff = None
    for gff in gff_files:
        try:
            _open = gzip.open if str(gff).endswith(".gz") else open
            with _open(gff, "rt") as f:
                has_features = False
                for line_num, line in enumerate(f):
                    if line_num > 500:
                        break
                    if "\tCDS\t" in line or "\tgene\t" in line:
                        has_features = True
                        break
                if has_features:
                    best_gff = gff
                    break
        except Exception:
            continue

    if not best_gff:
        return "no_features"

    if str(best_gff).endswith(".gz"):
        with gzip.open(best_gff, "rt") as src, open(gff_target, "w") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.move(str(best_gff), str(gff_target))
    return "ok"


def filter_chromosomes_only(fna_path: Path) -> None:
    """
    If a FASTA contains both chromosome sequences (NC_/CM_ accession prefixes)
    and unlocalized/placed scaffolds (NW_, NZ_, or other prefixes), rewrite the
    file keeping only the chromosome sequences.

    Rationale:
      Chromosome-level assemblies are distributed with the assembled chromosomes
      *plus* unlocalized/unplaced scaffolds in the same .fna.  Including NW_
      scaffolds alongside NC_ chromosomes adds noise to the synteny search:
      flanking genes may falsely appear multiple times (once on the chromosome,
      again on a scaffold that overlaps the same region).

    If NO chromosome sequences are found (e.g. pure scaffold assembly), the file
    is left untouched so we do not discard all sequence data.
    """
    # Prefixes that indicate assembled (placed) chromosomes in NCBI nomenclature.
    # NC_ = RefSeq chromosome; CM_ = GenBank chromosome/complete genomic molecule.
    CHROM_PREFIXES = ("NC_", "CM")

    # First pass: scan only FASTA headers (fast — avoids reading sequence content)
    chrom_ids = []
    non_chrom_ids = []
    try:
        with open(fna_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                # The sequence ID is the first whitespace-delimited token after ">"
                seq_id = line[1:].split()[0] if line.strip() else ""
                if any(seq_id.startswith(p) for p in CHROM_PREFIXES):
                    chrom_ids.append(seq_id)
                else:
                    non_chrom_ids.append(seq_id)
    except Exception as e:
        print(f"  [chr-filter] Warning: could not scan {fna_path}: {e}")
        return

    if not chrom_ids:
        # Pure scaffold/contig assembly — keep as-is
        print(
            f"  [chr-filter] No NC_/CM_ chromosomes found — retaining all "
            f"{len(non_chrom_ids)} sequences as-is (scaffold assembly)."
        )
        return

    if not non_chrom_ids:
        # Already chromosomes only
        print(
            f"  [chr-filter] {len(chrom_ids)} chromosome sequence(s) — no scaffold sequences to remove."
        )
        return

    # Mixed assembly: chromosomes + scaffolds. Filter to chromosomes only.
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
        # Atomic rename
        tmp_path.replace(fna_path)
        print(
            f"  [chr-filter] Kept {kept} chromosome sequence(s), "
            f"removed {len(non_chrom_ids)} scaffold/unlocalized sequence(s)."
        )
    except Exception as e:
        print(f"  [chr-filter] Warning: filter failed for {fna_path}: {e}")
        if tmp_path.exists():
            tmp_path.unlink()


def download_genome(accession, output_dir, max_retries=3):
    """Download genome from NCBI using datasets with retry. Also attempts to download GFF if available."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    zip_file = output_path / f"{accession}.zip"
    extract_dir = output_path / accession
    gff_target = output_path / f"{accession}.gff"

    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        "genome,gff3",
        "--filename",
        str(zip_file),
    ]

    for attempt in range(1, max_retries + 1):
        print(f"Downloading {accession} (attempt {attempt}/{max_retries})...")

        # Clean up artefacts from a previous failed attempt
        if zip_file.exists():
            zip_file.unlink()
        if extract_dir.exists():
            shutil.rmtree(extract_dir)

        try:
            run_safe_command(cmd)
            extract_zip_archive(zip_file, extract_dir)

            fna_files = list(extract_dir.rglob("*.fna"))
            fna_path = None

            if fna_files:
                target = output_path / f"{accession}.fna"
                shutil.move(str(fna_files[0]), str(target))
                fna_path = str(target)
                # Filter to chromosome sequences only when both chromosomes
                # and unlocalized scaffolds are present in the same assembly.
                filter_chromosomes_only(target)
                print(f"  ✓ Genome: {target}")

            gff_status = _extract_best_gff(extract_dir, gff_target)
            if gff_status == "ok":
                print(f"  ✓ GFF annotations: {gff_target}")
            elif gff_status == "no_features":
                print("  ○ GFF found but no CDS/gene features (scaffold-only)")
            else:
                print("  ○ No GFF annotations available")

            # Fallback: if GenBank accession has no useful GFF, try RefSeq counterpart.
            if gff_status != "ok" and accession.startswith("GCA_"):
                refseq_acc = "GCF_" + accession[4:]
                ref_zip = output_path / f"{refseq_acc}.zip"
                ref_extract = output_path / refseq_acc
                try:
                    print(f"  ○ Trying RefSeq annotation fallback: {refseq_acc}")
                    run_safe_command(
                        [
                            "datasets",
                            "download",
                            "genome",
                            "accession",
                            refseq_acc,
                            "--include",
                            "gff3",
                            "--filename",
                            str(ref_zip),
                        ]
                    )
                    extract_zip_archive(ref_zip, ref_extract)
                    ref_status = _extract_best_gff(ref_extract, gff_target)
                    if ref_status == "ok":
                        print(f"  ✓ GFF annotations (RefSeq fallback): {gff_target}")
                    else:
                        print("  ○ RefSeq fallback did not yield usable GFF")
                except Exception:
                    print("  ○ RefSeq annotation fallback unavailable")
                finally:
                    if ref_zip.exists():
                        ref_zip.unlink()
                    if ref_extract.exists():
                        shutil.rmtree(ref_extract)

            return fna_path

        except (subprocess.CalledProcessError, zipfile.BadZipFile) as e:
            msg = e.stderr if hasattr(e, 'stderr') else str(e)
            print(f"  ✗ Attempt {attempt} failed for {accession}: {msg}", file=sys.stderr)
            if attempt < max_retries:
                wait = 10 * attempt
                print(f"  Retrying in {wait}s...", file=sys.stderr)
                import time; time.sleep(wait)
            else:
                print(f"  ✗ All {max_retries} attempts exhausted for {accession}.", file=sys.stderr)
                return None
        finally:
            if zip_file.exists():
                zip_file.unlink()
            if extract_dir.exists():
                shutil.rmtree(extract_dir)


def write_quality_report(assemblies, output_dir):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    report_path = output_path / "assembly_quality.tsv"
    with open(report_path, "w") as out:
        out.write(
            "accession\tspecies\ttax_level\trefseq_category\tassembly_level\tchromosomes\t"
            "scaffolds\tcontigs\tcontig_n50\tscaffold_n50\tcontig_n80\tscaffold_n80\t"
            "bad_quality\tbad_reasons\n"
        )
        for asm in assemblies:
            out.write(
                f"{asm.get('accession', '')}\t{asm.get('species', '')}\t{asm.get('tax_level', '')}\t"
                f"{asm.get('category', '')}\t{asm.get('assembly_status', '')}\t"
                f"{asm.get('chromosome_count', '')}\t{asm.get('scaffold_count', '')}\t{asm.get('contig_count', '')}\t"
                f"{asm.get('contig_n50', '')}\t{asm.get('scaffold_n50', '')}\t"
                f"{asm.get('contig_n80', '')}\t{asm.get('scaffold_n80', '')}\t"
                f"{str(asm.get('bad_quality', False)).lower()}\t{asm.get('bad_quality_reasons', '')}\n"
            )
    print(f"Assembly quality report written to: {report_path}")


def write_outputs(assemblies, downloaded_paths, output_dir):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    manifest_path = output_path / "genomes_manifest.txt"
    with open(manifest_path, "w") as f:
        for path in downloaded_paths:
            f.write(f"{path}\n")
    print(f"Genome paths written to: {manifest_path}")

    downloaded_acc = {Path(p).stem for p in downloaded_paths}
    selected = [a for a in assemblies if a.get("accession") in downloaded_acc]
    species_map_path = output_path / "species_mapping.tsv"
    with open(species_map_path, "w") as f:
        for asm in selected:
            tax_level = asm.get("tax_level", "unknown")
            f.write(f"{asm['accession']}\t{asm['species']}\t{tax_level}\n")
    print(f"Species mapping written to: {species_map_path}")

    write_quality_report(assemblies, output_dir)


def print_selected_assemblies(assemblies, title):
    print(f"\n{title}:")
    for i, asm in enumerate(assemblies, 1):
        cat = asm.get("category", "na")
        print(f"  {i}. {asm['accession']} - {asm['species']} [{cat}]")
        print(f"     {format_quality(asm)}")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch related genomes from NCBI for easy mode",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--home-species",
        dest="home_species",
        required=True,
        help="Home species name (e.g., 'Apis mellifera') - searches genus, excludes this species",
    )
    parser.add_argument("--max", type=int, default=10, help="Maximum number of genomes to fetch (default: 10)")
    parser.add_argument(
        "--outdir",
        default="easy_mode_genomes",
        help="Output directory for genomes (default: easy_mode_genomes)",
    )
    parser.add_argument("--list-only", action="store_true", help="Only list available genomes, don't download")
    parser.add_argument(
        "--target-species",
        dest="target_species",
        default=None,
        help="Comma-separated list of species names to fetch (bypasses automatic search)",
    )
    parser.add_argument(
        "--assembly-ranking",
        choices=["hybrid", "counts", "nstats"],
        default="hybrid",
        help="Ranking strategy for assembly quality selection (default: hybrid)",
    )
    parser.add_argument(
        "--bad-quality-policy",
        choices=["ask", "drop", "keep"],
        default="ask",
        help="What to do when only low-quality assembly is available for a species (default: ask)",
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
        default=500000,
        help="Assemblies above this contig count are flagged low quality (default: 500000)",
    )
    parser.add_argument(
        "--bad-max-scaffolds",
        type=int,
        default=500000,
        help="Assemblies above this scaffold count are flagged low quality (default: 500000)",
    )
    parser.add_argument(
        "--bad-min-n50",
        type=int,
        default=5000,
        help="Assemblies with best N50 below this are flagged low quality (default: 5000)",
    )
    args = parser.parse_args()

    try:
        run_safe_command(["which", "esearch"])
    except Exception:
        print("ERROR: NCBI E-utilities (esearch) not found!", file=sys.stderr)
        sys.exit(1)

    metadata_cache = {}

    # --- TARGET SPECIES MODE ---
    if args.target_species:
        species_list = [s.strip() for s in args.target_species.split(",") if s.strip()]
        print(f"Target species mode: fetching {len(species_list)} specified species")

        assemblies = []
        for sp_name in species_list:
            print(f"\nLooking up '{sp_name}'...")
            sp_taxid = get_taxid_from_name(sp_name)
            if not sp_taxid:
                print(f"  WARNING: Could not find taxonomy ID for '{sp_name}', skipping")
                continue

            sp_assemblies = get_related_species(
                sp_taxid,
                max_genomes=10,
                exclude_species=args.home_species,
                tax_level=f"target ({sp_name})",
                ranking_mode=args.assembly_ranking,
                metadata_cache=metadata_cache,
            )
            if not sp_assemblies:
                sp_assemblies = get_related_species(
                    sp_taxid,
                    max_genomes=10,
                    exclude_species=None,
                    tax_level=f"target ({sp_name})",
                    ranking_mode=args.assembly_ranking,
                    metadata_cache=metadata_cache,
                )

            if sp_assemblies:
                exact = [a for a in sp_assemblies if normalize_species(a["species"]) == normalize_species(sp_name)]
                best = exact[0] if exact else sp_assemblies[0]
                assemblies.append(best)
                print(
                    f"  Best assembly: {best['accession']} ({best['species']}, {best.get('category', 'na')})"
                )
                print(f"    {format_quality(best)}")
            else:
                print(f"  WARNING: No assemblies found for '{sp_name}'")

        assemblies = apply_bad_quality_policy(assemblies, args)

        if not assemblies:
            print("ERROR: No acceptable assemblies found for specified species", file=sys.stderr)
            output_path = Path(args.outdir)
            output_path.mkdir(parents=True, exist_ok=True)
            (output_path / "genomes_manifest.txt").touch()
            (output_path / "species_mapping.tsv").touch()
            write_quality_report([], args.outdir)
            sys.exit(0)

        print_selected_assemblies(assemblies, f"Selected {len(assemblies)} assembly(ies)")

        if args.list_only:
            print("\nList-only mode. Exiting without download.")
            write_quality_report(assemblies, args.outdir)
            return

        print(f"\nDownloading genomes to {args.outdir}/...")
        downloaded = []
        for asm in assemblies:
            fna_path = download_genome(asm["accession"], args.outdir)
            if fna_path:
                downloaded.append(fna_path)

        print(f"\n{'=' * 60}")
        print(f"✓ Successfully downloaded {len(downloaded)}/{len(assemblies)} genomes")
        print(f"{'=' * 60}")
        write_outputs(assemblies, downloaded, args.outdir)
        return

    # --- AUTOMATIC TAXONOMIC SEARCH MODE ---
    genus = args.home_species.split()[0]
    print(f"Home species: {args.home_species}")

    print(f"Looking up taxonomy for '{args.home_species}'...")
    species_taxid = get_taxid_from_name(args.home_species)
    if not species_taxid:
        species_taxid = get_taxid_from_name(genus)
    if not species_taxid:
        print(f"ERROR: Could not find taxonomy ID for '{args.home_species}'", file=sys.stderr)
        sys.exit(1)
    print(f"  Species TaxID: {species_taxid}")

    genus_taxid = get_taxid_from_name(genus)
    search_levels = []
    if genus_taxid:
        search_levels.append((f"genus ({genus})", genus_taxid))
    parent_ranks = get_parent_taxa(species_taxid)
    search_levels.extend(parent_ranks)

    if not search_levels:
        print("ERROR: Could not determine taxonomy levels", file=sys.stderr)
        sys.exit(1)
    print(f"  Search levels: {', '.join(name for name, _ in search_levels)}")

    max_genomes = args.max
    if max_genomes <= 0:
        n_levels = len(search_levels)
        max_genomes = min(n_levels * 3, 20)
        print(f"  Auto genome count: {n_levels} levels x 3 = {max_genomes} genomes")
    print(f"  Target: {max_genomes} genomes\n")

    n_levels = len(search_levels)
    base_budget = max(2, max_genomes // n_levels)
    remainder = max_genomes - (base_budget * n_levels)
    level_budgets = [base_budget] * n_levels
    for i in range(remainder):
        mid_idx = min(1 + i, n_levels - 1)
        level_budgets[mid_idx] += 1
    budget_str = ", ".join(f"{name}: {b}" for (name, _), b in zip(search_levels, level_budgets))
    print(f"  Budget allocation: {budget_str}\n")

    collected = []
    seen_species = {normalize_species(args.home_species)}
    carry_over = 0

    for i, (level_name, level_taxid) in enumerate(search_levels):
        budget = level_budgets[i] + carry_over
        carry_over = 0
        if len(collected) >= max_genomes:
            break
        budget = min(budget, max_genomes - len(collected))
        print(f"Level: {level_name} (budget: {budget})")

        broader_count = max(50, budget * 10)
        candidates = get_related_species(
            level_taxid,
            broader_count,
            exclude_species=args.home_species,
            tax_level=level_name,
            ranking_mode=args.assembly_ranking,
            metadata_cache=metadata_cache,
        )

        added = 0
        for c in candidates:
            sp_key = normalize_species(c["species"])
            if sp_key not in seen_species:
                seen_species.add(sp_key)
                collected.append(c)
                added += 1
                if added >= budget:
                    break

        if added < budget:
            carry_over = budget - added
            print(f"    Added {added}/{budget} (carrying {carry_over} to next level)")
        else:
            print(f"    Added {added} new species (total: {len(collected)}/{max_genomes})")

    assemblies = apply_bad_quality_policy(collected, args)

    if not assemblies:
        print("WARNING: No acceptable related assemblies found.", file=sys.stderr)
        output_path = Path(args.outdir)
        output_path.mkdir(parents=True, exist_ok=True)
        (output_path / "genomes_manifest.txt").touch()
        (output_path / "species_mapping.tsv").touch()
        write_quality_report([], args.outdir)
        sys.exit(0)

    print_selected_assemblies(assemblies, f"Found {len(assemblies)} related assemblies")

    if args.list_only:
        print("\nList-only mode. Exiting without download.")
        write_quality_report(assemblies, args.outdir)
        sys.exit(0)

    print(f"\nDownloading genomes to {args.outdir}/...")
    downloaded = []
    for asm in assemblies:
        fna_path = download_genome(asm["accession"], args.outdir)
        if fna_path:
            downloaded.append(fna_path)

    print(f"\n{'=' * 60}")
    print(f"✓ Successfully downloaded {len(downloaded)}/{len(assemblies)} genomes")
    print(f"{'=' * 60}")
    write_outputs(assemblies, downloaded, args.outdir)


if __name__ == "__main__":
    main()
