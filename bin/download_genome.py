#!/usr/bin/env python3
"""
Robust genome downloader for SynTerra.

Strategy (best-of-both-worlds hybrid):
  1. Use `datasets summary` to resolve accession → direct NCBI FTP/HTTPS URL.
  2. Download the .fna.gz directly with wget -c (resume on failure, progress bar).
  3. Fall back to `datasets download` (streaming) if no FTP path is available.

Why not pure `datasets download`?
  - No --resume flag: every failure restarts from 0 bytes.
  - Uses HTTP/2 multiplexed streams that NCBI servers occasionally reset mid-download.
  - ZIP wrapper adds one more failure mode (bad zip on truncation).

Why not pure wget?
  - Doesn't handle accession → URL resolution.
  - Doesn't handle the genome+gff3 packaging logic.
  - Doesn't know about best-assembly selection.

Hybrid gives us:
  ✓ NCBI official metadata/accession resolution (datasets)
  ✓ Resumable downloads (-c flag)
  ✓ Automatic retry with backoff (--tries, --retry-connrefused)
  ✓ Live progress bar
  ✓ Direct .fna.gz (no zip extraction step = less temp disk)
  ✓ GFF fallback still via datasets (small download)
"""

import gzip
import json
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def check_disk_space(path: Path, required_gb: float = 5.0) -> None:
    """Raise RuntimeError if free disk space at *path* is below required_gb."""
    st = shutil.disk_usage(path)
    free_gb = st.free / (1024 ** 3)
    if free_gb < required_gb:
        raise RuntimeError(
            f"Insufficient disk space at {path}: {free_gb:.1f} GB free, "
            f"{required_gb:.1f} GB required. Free up space and retry."
        )


def run_streaming(cmd: list, label: str = "") -> None:
    """
    Run *cmd* with live stdout/stderr forwarded to the console.
    Raises subprocess.CalledProcessError on non-zero exit (with stderr captured).
    """
    stderr_lines: list = []
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    ) as proc:
        for line in proc.stdout:  # type: ignore[union-attr]
            print(line, end="", flush=True)
        stderr_lines = proc.stderr.read().splitlines()  # type: ignore[union-attr]
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(
                proc.returncode, cmd, stderr="\n".join(stderr_lines)
            )


def run_safe(cmd: list, timeout: int = 60) -> str:
    """Run command, capture and return stdout. Raises on failure."""
    result = subprocess.run(
        cmd, capture_output=True, text=True, check=True, timeout=timeout
    )
    return result.stdout.strip()


# ---------------------------------------------------------------------------
# FTP path resolution via E-utilities
# ---------------------------------------------------------------------------

def get_ftp_path(accession: str) -> str | None:
    """
    Resolve an assembly accession to its NCBI FTP base directory.
    Returns a string like:
      https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14
    or None if not resolvable.
    """
    try:
        procs = []
        p1 = subprocess.Popen(
            ["esearch", "-db", "assembly", "-query", accession],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        procs.append(p1)
        p2 = subprocess.Popen(
            ["efetch", "-format", "docsum"],
            stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        p1.stdout.close()
        procs.append(p2)
        p3 = subprocess.Popen(
            ["xtract", "-pattern", "DocumentSummary", "-element",
             "FtpPath_RefSeq", "FtpPath_GenBank"],
            stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        p2.stdout.close()
        procs.append(p3)
        out, _ = p3.communicate(timeout=45)
        line = out.decode().strip().split("\n")[0] if out else ""
        if not line:
            return None
        parts = line.split("\t")
        for p in parts:
            p = p.strip()
            # Prefer HTTPS; efetch returns ftp:// but NCBI supports https://
            if p and p not in ("na", ""):
                return p.replace("ftp://", "https://")
        return None
    except Exception as e:
        print(f"  [ftp-resolve] Warning: {e}")
        return None


def ftp_url_to_fna_url(ftp_base: str) -> tuple[str, str]:
    """
    Given the FTP base directory URL, return (fna_url, gff_url).
    NCBI naming convention: the last path component is the assembly name.
    e.g. .../GCF_000001405.40_GRCh38.p14  →
         .../GCF_000001405.40_GRCh38.p14_genomic.fna.gz
         .../GCF_000001405.40_GRCh38.p14_genomic.gff.gz
    """
    base = ftp_base.rstrip("/")
    assembly_name = base.split("/")[-1]
    fna_url = f"{base}/{assembly_name}_genomic.fna.gz"
    gff_url = f"{base}/{assembly_name}_genomic.gff.gz"
    return fna_url, gff_url


# ---------------------------------------------------------------------------
# wget-based download with resume
# ---------------------------------------------------------------------------

def wget_download(url: str, dest: Path, max_retries: int = 5) -> bool:
    """
    Download *url* to *dest* using wget with resume (-c) and retry.
    Returns True on success.
    """
    cmd = [
        "wget",
        "--continue",           # -c: resume partial downloads
        "--tries", str(max_retries),
        "--retry-connrefused",
        "--timeout", "60",      # per-read timeout
        "--waitretry", "15",    # wait between retries
        "--show-progress",
        "--progress", "bar:force",
        "-O", str(dest),
        url,
    ]
    print(f"  wget: {url}")
    try:
        # wget prints progress to stderr; use Popen to forward it live
        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # merge so progress bar shows
            text=True,
            bufsize=1,
        ) as proc:
            for line in proc.stdout:  # type: ignore[union-attr]
                print(line, end="", flush=True)
            proc.wait()
            if proc.returncode != 0:
                print(f"  wget exited {proc.returncode}", file=sys.stderr)
                return False
        return True
    except FileNotFoundError:
        print("  wget not found — falling back to datasets", file=sys.stderr)
        return False


# ---------------------------------------------------------------------------
# datasets fallback (streaming, no resume)
# ---------------------------------------------------------------------------

def datasets_download(accession: str, output_path: Path,
                      fna_target: Path, gff_target: Path,
                      max_retries: int = 3) -> tuple[bool, bool]:
    """
    Fall back to `datasets download ... --include genome,gff3`.
    Returns (fna_ok, gff_ok).
    """
    zip_file = output_path / f"{accession}.zip"
    extract_dir = output_path / f"{accession}_extracted"
    cmd = [
        "datasets", "download", "genome", "accession", accession,
        "--include", "genome,gff3",
        "--filename", str(zip_file),
    ]
    for attempt in range(1, max_retries + 1):
        print(f"  datasets download (attempt {attempt}/{max_retries})...")
        if zip_file.exists():
            zip_file.unlink()
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        try:
            run_streaming(cmd, label=accession)
            extract_dir.mkdir(exist_ok=True)
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(extract_dir)
            fna_files = list(extract_dir.rglob("*.fna"))
            gff_files = [f for pat in ("*.gff", "*.gff3") for f in extract_dir.rglob(pat)]
            fna_ok = False
            gff_ok = False
            if fna_files:
                shutil.copy(fna_files[0], fna_target)
                fna_ok = True
            if gff_files:
                shutil.copy(gff_files[0], gff_target)
                gff_ok = True
            return fna_ok, gff_ok
        except (subprocess.CalledProcessError, zipfile.BadZipFile, OSError) as e:
            msg = e.stderr if hasattr(e, "stderr") else str(e)
            print(f"  datasets attempt {attempt} failed: {msg}", file=sys.stderr)
            if zip_file.exists():
                zip_file.unlink()
            if extract_dir.exists():
                shutil.rmtree(extract_dir)
            if attempt < max_retries:
                wait = 15 * attempt
                print(f"  Retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
        finally:
            if zip_file.exists():
                zip_file.unlink()
            if extract_dir.exists():
                shutil.rmtree(extract_dir)
    return False, False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def download_genome_robust(
    accession: str,
    output_path: Path,
    fna_name: str = "genome.fna",
    gff_name: str = "genome.gff",
    max_retries: int = 5,
) -> tuple[Path | None, Path | None]:
    """
    Download genome FASTA (and optionally GFF) for *accession*.

    Strategy:
      1. Resolve FTP URL via E-utilities.
      2. If resolved: wget -c .fna.gz → decompress → filter chromosomes.
                      wget -c .gff.gz → decompress.
      3. If FTP not available: fall back to datasets download (streaming).

    Returns (fna_path, gff_path) — either may be None if unavailable.
    """
    output_path.mkdir(parents=True, exist_ok=True)

    # Pre-flight disk check
    check_disk_space(output_path, required_gb=5.0)

    fna_target = output_path / fna_name
    gff_target = output_path / gff_name

    print(f"\n[download] Resolving FTP path for {accession}...")
    ftp_base = get_ftp_path(accession)

    if ftp_base:
        fna_gz_url, gff_gz_url = ftp_url_to_fna_url(ftp_base)
        fna_gz = output_path / f"{accession}_genomic.fna.gz"
        gff_gz = output_path / f"{accession}_genomic.gff.gz"

        print(f"  FTP base: {ftp_base}")

        # --- Download FASTA ---
        fna_ok = False
        for attempt in range(1, max_retries + 1):
            check_disk_space(output_path, required_gb=5.0)
            if wget_download(fna_gz_url, fna_gz):
                try:
                    print(f"  Decompressing {fna_gz.name}...")
                    with gzip.open(fna_gz, "rb") as f_in, open(fna_target, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    fna_gz.unlink()
                    fna_ok = True
                    break
                except (OSError, EOFError) as e:
                    print(f"  Decompression failed: {e} — retrying", file=sys.stderr)
                    if fna_gz.exists():
                        fna_gz.unlink()  # Remove corrupt gz so wget restarts
            if attempt < max_retries:
                time.sleep(15 * attempt)

        # --- Download GFF ---
        gff_ok = False
        if wget_download(gff_gz_url, gff_gz):
            try:
                print(f"  Decompressing {gff_gz.name}...")
                with gzip.open(gff_gz, "rb") as f_in, open(gff_target, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                gff_gz.unlink()
                gff_ok = True
            except (OSError, EOFError) as e:
                print(f"  GFF decompression failed: {e} (non-fatal)", file=sys.stderr)
                if gff_gz.exists():
                    gff_gz.unlink()

        if fna_ok:
            return fna_target, (gff_target if gff_ok else None)

        print("  FTP download failed — falling back to datasets...", file=sys.stderr)

    # Fallback: datasets download (no resume, but reliable for small genomes)
    print(f"  Using datasets download fallback for {accession}")
    fna_ok, gff_ok = datasets_download(
        accession, output_path, fna_target, gff_target, max_retries=3
    )
    return (fna_target if fna_ok else None), (gff_target if gff_ok else None)
