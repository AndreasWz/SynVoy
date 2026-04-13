#!/usr/bin/env python3
"""
gene_predictor.py — Unified gene/ORF prediction for SynVoy

Dispatches to the appropriate gene predictor based on organism domain:
  - Prokaryotes: Prodigal (metagenomic mode)
  - Eukaryotes:  Augustus (with species model)

SynVoy primarily works with eukaryotic genomes (insects, vertebrates, etc.),
so Augustus is the default.  Prodigal is a prokaryotic gene finder — it cannot
predict intron-containing genes and will miss most eukaryotic genes.

Usage:
    from gene_predictor import predict_genes, predict_orfs

    # Full gene prediction (for flanking gene extraction, home proteome)
    genes = predict_genes(fasta_in, faa_out, gff_out,
                          predictor="augustus", augustus_species="fly")

    # ORF prediction on small regions (for PLM/structural search)
    orfs = predict_orfs(region_fasta, output_dir,
                        predictor="augustus", augustus_species="fly")
"""

import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------

def check_prodigal_available() -> bool:
    return shutil.which("prodigal") is not None


def check_augustus_available() -> bool:
    return shutil.which("augustus") is not None


def check_predictor_available(predictor: str) -> bool:
    """Check if the requested predictor binary is on PATH."""
    if predictor == "prodigal":
        return check_prodigal_available()
    if predictor == "augustus":
        return check_augustus_available()
    return False


# ---------------------------------------------------------------------------
# Augustus species model helpers
# ---------------------------------------------------------------------------

# Common Augustus species models grouped by clade.
# Used by auto-detection and as documentation for users.
AUGUSTUS_SPECIES_MAP = {
    # Insects
    "fly": "fly",
    "drosophila": "fly",
    "honeybee": "honeybee1",
    "apis": "honeybee1",
    "beetle": "tribolium2012",
    "tribolium": "tribolium2012",
    "mosquito": "aedes",
    "aedes": "aedes",
    "bombus": "bombus_terrestris1",
    "butterfly": "heliconius_melpomene1",
    "silkworm": "bombyx_mori",
    "ant": "nasonia",  # nasonia is a close hymenoptera model
    "wasp": "nasonia",
    "nasonia": "nasonia",
    # Vertebrates
    "human": "human",
    "mouse": "mouse",
    "chicken": "chicken",
    "zebrafish": "zebrafish",
    # Other eukaryotes
    "nematode": "caenorhabditis",
    "c_elegans": "caenorhabditis",
    "yeast": "saccharomyces",
    "arabidopsis": "arabidopsis",
    "rice": "rice",
    "maize": "maize5",
    # Generic
    "generic": "fly",  # fly is a reasonable general eukaryote model
}


def resolve_augustus_species(species_hint: str) -> str:
    """
    Resolve a user-friendly species hint to an Augustus --species parameter.

    Accepts:
      - Direct Augustus species names (e.g., "honeybee1", "nasonia")
      - Common names from AUGUSTUS_SPECIES_MAP (e.g., "fly", "honeybee")
      - Falls back to the hint as-is (user may have a custom model)
    """
    if not species_hint:
        return "fly"  # safe default for eukaryotes
    hint_lower = species_hint.strip().lower()
    return AUGUSTUS_SPECIES_MAP.get(hint_lower, species_hint)


def list_augustus_species() -> List[str]:
    """List available Augustus species models on this system."""
    try:
        result = subprocess.run(
            ["augustus", "--species=help"],
            capture_output=True, text=True, timeout=10,
        )
        # Augustus prints species list to stderr
        output = result.stderr + result.stdout
        species = []
        for line in output.splitlines():
            line = line.strip()
            if line and not line.startswith("#") and not line.startswith("usage"):
                # Augustus lists species one per line, sometimes with descriptions
                name = line.split()[0] if line.split() else ""
                if name and not name.startswith("-"):
                    species.append(name)
        return species
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return []


# ---------------------------------------------------------------------------
# Prodigal gene prediction
# ---------------------------------------------------------------------------

def run_prodigal(
    fasta_in: str,
    faa_out: str,
    gff_out: str,
    meta_mode: bool = True,
    quiet: bool = True,
) -> bool:
    """
    Run Prodigal on a FASTA file.  Best for prokaryotic genomes or
    metagenomic fragments.

    Returns True on success.
    """
    cmd = ["prodigal", "-i", fasta_in, "-a", faa_out, "-f", "gff", "-o", gff_out]
    if meta_mode:
        cmd.extend(["-p", "meta"])
    if quiet:
        cmd.append("-q")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as exc:
        logger.debug(f"Prodigal failed: {exc}")
        return False


# ---------------------------------------------------------------------------
# Augustus gene prediction
# ---------------------------------------------------------------------------

def run_augustus(
    fasta_in: str,
    faa_out: str,
    gff_out: str,
    species: str = "fly",
    softmasking: bool = True,
) -> bool:
    """
    Run Augustus on a FASTA file.  Designed for eukaryotic genomes —
    handles introns, splice sites, and multi-exon gene structures.

    Returns True on success.
    """
    species = resolve_augustus_species(species)

    cmd = [
        "augustus",
        f"--species={species}",
        "--gff3=on",
        "--protein=on",
        "--codingseq=on",
    ]
    if softmasking:
        cmd.append("--softmasking=1")
    cmd.append(fasta_in)

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            logger.debug(f"Augustus exited {result.returncode}: {result.stderr[:300]}")
            return False

        # Augustus writes GFF3 to stdout
        gff_content = result.stdout

        # Parse GFF and extract protein sequences
        _parse_augustus_output(gff_content, faa_out, gff_out)
        return True

    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        logger.debug(f"Augustus failed: {exc}")
        return False


def _parse_augustus_output(
    gff_content: str,
    faa_out: str,
    gff_out: str,
) -> None:
    """
    Parse Augustus GFF3+protein stdout into separate GFF and FAA files.

    Augustus emits GFF3 lines followed by protein sequences in FASTA-like
    comment blocks:  # protein sequence = [MXXX...XXX]
    """
    gff_lines = []
    proteins: Dict[str, str] = {}
    current_gene_id = ""
    current_protein_lines: List[str] = []
    in_protein_block = False

    for line in gff_content.splitlines():
        if line.startswith("#"):
            # Check for protein sequence blocks
            # Augustus format: # protein sequence = [MSEQ...]
            if "protein sequence" in line:
                # Start of protein block
                match = re.search(r"protein sequence = \[(.+)", line)
                if match:
                    seq_part = match.group(1)
                    if seq_part.endswith("]"):
                        # Single-line protein
                        seq_part = seq_part.rstrip("]")
                        if current_gene_id and seq_part:
                            proteins[current_gene_id] = seq_part
                    else:
                        in_protein_block = True
                        current_protein_lines = [seq_part]
            elif in_protein_block:
                # Continuation of multi-line protein
                seq_part = line.lstrip("# ").strip()
                if seq_part.endswith("]"):
                    seq_part = seq_part.rstrip("]")
                    current_protein_lines.append(seq_part)
                    if current_gene_id:
                        proteins[current_gene_id] = "".join(current_protein_lines)
                    current_protein_lines = []
                    in_protein_block = False
                else:
                    current_protein_lines.append(seq_part)
            continue

        if not line.strip():
            continue

        parts = line.split("\t")
        if len(parts) >= 9:
            gff_lines.append(line)

            # Track gene IDs for protein association
            feature_type = parts[2]
            if feature_type == "gene":
                attrs = parts[8]
                gene_match = re.search(r"ID=([^;]+)", attrs)
                if gene_match:
                    current_gene_id = gene_match.group(1)

    # Write GFF
    with open(gff_out, "w") as f:
        f.write("##gff-version 3\n")
        for line in gff_lines:
            f.write(line + "\n")

    # Write proteins in Prodigal-compatible FAA format
    # We match Prodigal's header format: >gene_id # start # end # strand # attrs
    with open(faa_out, "w") as f:
        for gene_id, seq in proteins.items():
            if not seq or len(seq) < 5:
                continue
            # Find the gene's coordinates from GFF lines
            start, end, strand = _find_gene_coords(gff_lines, gene_id)
            strand_code = "1" if strand == "+" else "-1"
            f.write(f">{gene_id} # {start} # {end} # {strand_code} # ID={gene_id}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def _find_gene_coords(
    gff_lines: List[str], gene_id: str,
) -> Tuple[int, int, str]:
    """Find start, end, strand for a gene ID from GFF lines."""
    for line in gff_lines:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        if parts[2] not in ("gene", "mRNA", "transcript"):
            continue
        if f"ID={gene_id}" in parts[8] or f"Parent={gene_id}" in parts[8]:
            try:
                return int(parts[3]), int(parts[4]), parts[6]
            except (ValueError, IndexError):
                pass
    return 0, 0, "+"


# ---------------------------------------------------------------------------
# Unified prediction interface
# ---------------------------------------------------------------------------

def predict_genes(
    fasta_in: str,
    faa_out: str,
    gff_out: str,
    predictor: str = "augustus",
    augustus_species: str = "fly",
) -> bool:
    """
    Predict genes in a FASTA file using the specified predictor.

    Args:
        fasta_in:        input genome/region FASTA
        faa_out:         output protein FASTA
        gff_out:         output GFF3
        predictor:       'prodigal', 'augustus', or 'auto'
        augustus_species: Augustus species model (ignored for prodigal)

    Returns:
        True on success
    """
    predictor = _resolve_predictor(predictor)

    if predictor == "prodigal":
        return run_prodigal(fasta_in, faa_out, gff_out)
    elif predictor == "augustus":
        success = run_augustus(fasta_in, faa_out, gff_out, species=augustus_species)
        if not success and check_prodigal_available():
            logger.warning("Augustus failed; falling back to Prodigal")
            return run_prodigal(fasta_in, faa_out, gff_out)
        return success
    else:
        logger.error(f"Unknown predictor: {predictor}")
        return False


def predict_orfs(
    region_fasta: str,
    output_dir: str,
    predictor: str = "augustus",
    augustus_species: str = "fly",
    min_aa: int = 20,
) -> List[Dict[str, Any]]:
    """
    Predict ORFs/genes in a genomic region.  Replacement for the old
    predict_orfs_prodigal() that only supported prokaryotic prediction.

    Args:
        region_fasta:    input region FASTA
        output_dir:      working directory for temp files
        predictor:       'prodigal', 'augustus', or 'auto'
        augustus_species: Augustus species model
        min_aa:          minimum protein length

    Returns:
        list of dicts: {id, seq, start (0-based), end (exclusive), strand}
    """
    predictor = _resolve_predictor(predictor)

    proteins_file = os.path.join(output_dir, "predicted_orfs.faa")
    gff_file = os.path.join(output_dir, "predicted_orfs.gff")

    if predictor == "prodigal":
        success = run_prodigal(region_fasta, proteins_file, gff_file)
    elif predictor == "augustus":
        success = run_augustus(
            region_fasta, proteins_file, gff_file,
            species=augustus_species,
        )
        if not success and check_prodigal_available():
            logger.debug("Augustus failed on region; falling back to Prodigal")
            success = run_prodigal(region_fasta, proteins_file, gff_file)
    else:
        success = False

    if not success:
        return []

    if not os.path.exists(proteins_file) or os.path.getsize(proteins_file) == 0:
        return []

    return _parse_predicted_faa(proteins_file, min_aa=min_aa)


def _parse_predicted_faa(
    faa_path: str,
    min_aa: int = 20,
) -> List[Dict[str, Any]]:
    """
    Parse a predicted protein FASTA (from either Prodigal or Augustus)
    with Prodigal-style headers: >id # start # end # strand # attrs

    Returns list of dicts: {id, seq, start, end, strand}
    """
    orfs: List[Dict[str, Any]] = []
    current_id: Optional[str] = None
    current_seq: List[str] = []
    current_meta: Dict[str, Any] = {}

    with open(faa_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    seq = "".join(current_seq).rstrip("*")
                    if len(seq) >= min_aa:
                        orfs.append({"id": current_id, "seq": seq, **current_meta})

                parts = line[1:].split(" # ")
                current_id = parts[0].strip()
                current_seq = []
                current_meta = {}
                if len(parts) >= 4:
                    try:
                        current_meta["start"] = int(parts[1]) - 1  # 0-based
                        current_meta["end"] = int(parts[2])        # exclusive
                        strand_raw = parts[3].strip()
                        if strand_raw in ("1", "+"):
                            current_meta["strand"] = "+"
                        elif strand_raw in ("-1", "-"):
                            current_meta["strand"] = "-"
                        else:
                            current_meta["strand"] = "+"
                    except (ValueError, IndexError):
                        pass
            else:
                current_seq.append(line)

    # Last record
    if current_id and current_seq:
        seq = "".join(current_seq).rstrip("*")
        if len(seq) >= min_aa:
            orfs.append({"id": current_id, "seq": seq, **current_meta})

    return orfs


def _resolve_predictor(predictor: str) -> str:
    """
    Resolve 'auto' to the best available predictor.
    Prefers Augustus (eukaryote-capable) over Prodigal.
    """
    if predictor == "auto":
        if check_augustus_available():
            return "augustus"
        if check_prodigal_available():
            logger.warning(
                "Augustus not available; falling back to Prodigal. "
                "Prodigal is designed for prokaryotes and will miss "
                "intron-containing eukaryotic genes. "
                "Install Augustus: conda install -c bioconda augustus"
            )
            return "prodigal"
        logger.error("No gene predictor available (neither Augustus nor Prodigal)")
        return "prodigal"  # will fail gracefully at runtime
    return predictor
