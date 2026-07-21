#!/usr/bin/env python3
"""Stage user-supplied target-genome annotations next to their genome FASTAs.

Pro-Mode counterpart to Easy Mode's automatic annotation download. In Easy Mode
``fetch_related_genomes.py`` writes ``{accession}.gff`` alongside
``{accession}.fna``, which is exactly where the search layer looks for it:

  * ``iterative_search_runner.find_native_annotation_path()`` probes
    ``<genome>.gff`` / ``.gff3`` (+ ``.gz``) next to each genome FASTA and, when
    found, stamps ``TargetGene`` / ``TargetProduct`` / ``TargetID`` onto every
    emitted model and feeds the family-consistency gate.
  * ``borrow_annotations.find_annotated_target()`` looks for ``<genome>.gff``
    when the HOME genome has no annotation of its own.

In Pro Mode the user passes plain FASTAs (``--target_genomes`` is FASTA-only by
design), so neither consumer ever saw an annotation. This script closes that
gap: given a separate folder/glob of GFFs (``--target_gffs``), it matches each
annotation to its genome and links it into the staged genomes directory under
the ONE filename both consumers probe first.

Naming is not cosmetic -- it is a correctness constraint:

    ``modules/cluster_regions.nf`` resolves the genome with
    ``find -L <dir> -name "<genome_name>*" -type f | head -n 1`` and
    ``<genome_name>`` carries the extension (e.g. ``cow.fna``). Naming the
    annotation ``cow.fna.gff`` would let that glob return the GFF INSTEAD of the
    genome. ``Path(genome).with_suffix('.gff')`` (-> ``cow.gff``) is both the
    first probe of ``find_native_annotation_path`` and immune to that glob.

Matching is tiered and refuses to guess: exact stem, then assembly accession,
then unambiguous one-sided prefix. Anything ambiguous or unmatched is reported
loudly (and is fatal unless ``--allow_unmatched``) rather than silently
producing a run with no annotations -- the failure mode this pipeline has been
bitten by repeatedly.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)
from sequence_utils import parse_gff  # noqa: E402

FASTA_EXTS = (".fna", ".fa", ".fasta", ".fas", ".seq")
GFF_EXTS = (".gff", ".gff3", ".gtf")
# Sentinel placeholders the workflow passes when the user supplied nothing.
SENTINELS = {"NO_GFFS", "NO_GFF"}
ACCESSION_RE = re.compile(r"(GC[AF]_\d+(?:\.\d+)?)")
# Feature types load_native_annotation_index() will actually index.
USABLE_TYPES = {"gene", "mRNA", "mrna", "transcript"}


def strip_gz(name: str) -> str:
    """Drop a trailing .gz so extension logic sees the real file type."""
    return name[:-3] if name.lower().endswith(".gz") else name


def file_stem(name: str, exts: Tuple[str, ...]) -> str:
    """Filename minus .gz and minus one known extension (case-insensitive)."""
    base = strip_gz(name)
    lowered = base.lower()
    for ext in exts:
        if lowered.endswith(ext):
            return base[: -len(ext)]
    return base


def accession_of(name: str) -> Optional[str]:
    """Extract a GCA_/GCF_ assembly accession from a filename, if present."""
    m = ACCESSION_RE.search(name)
    return m.group(1) if m else None


def native_annotation_name(genome_filename: str) -> str:
    """The filename an annotation must have to be found next to this genome.

    Mirrors the FIRST candidate of
    ``iterative_search_runner.find_native_annotation_path()`` so the search
    layer picks it up, while staying clear of cluster_regions.nf's
    ``<genome_name>*`` glob (see module docstring).
    """
    return Path(genome_filename).with_suffix(".gff").name


def list_genomes(genomes_dir: Path) -> List[Path]:
    """Genome FASTAs in the staged directory (same extensions the rest of the
    pipeline recognises). Symlinks are expected -- STAGE_GENOMES links them."""
    out = []
    for p in sorted(genomes_dir.iterdir()):
        if not p.is_file() and not p.is_symlink():
            continue
        base = strip_gz(p.name).lower()
        if base.endswith(FASTA_EXTS):
            out.append(p)
    return out


def collect_gffs(paths: List[str]) -> List[Path]:
    """Filter the raw --gffs list down to real annotation files."""
    out = []
    for raw in paths:
        p = Path(raw)
        if p.name in SENTINELS:
            continue
        if not p.exists() or not p.is_file():
            continue
        if not strip_gz(p.name).lower().endswith(GFF_EXTS):
            continue
        out.append(p)
    return out


def match_gffs_to_genomes(genomes: List[Path],
                          gffs: List[Path]) -> Tuple[Dict[str, Path], List[Tuple[Path, str]]]:
    """Assign at most one GFF to each genome.

    Returns ``(by_genome_name, problems)`` where ``problems`` is a list of
    ``(gff_path, reason)`` for annotations that could not be assigned.
    """
    genome_by_stem: Dict[str, List[Path]] = defaultdict(list)
    genome_by_acc: Dict[str, List[Path]] = defaultdict(list)
    for g in genomes:
        genome_by_stem[file_stem(g.name, FASTA_EXTS).lower()].append(g)
        acc = accession_of(g.name)
        if acc:
            genome_by_acc[acc.lower()].append(g)

    assigned: Dict[str, Path] = {}
    claimed_by: Dict[str, Path] = {}
    problems: List[Tuple[Path, str]] = []

    def claim(genome: Path, gff: Path) -> bool:
        prev = claimed_by.get(genome.name)
        if prev is not None:
            problems.append(
                (gff, f"genome '{genome.name}' already annotated by "
                      f"'{prev.name}' -- refusing to overwrite")
            )
            return False
        assigned[genome.name] = gff
        claimed_by[genome.name] = gff
        return True

    for gff in gffs:
        stem = file_stem(gff.name, GFF_EXTS)
        stem_l = stem.lower()

        # Tier 1: exact stem match (cow.gff <-> cow.fna).
        cands = genome_by_stem.get(stem_l, [])
        if len(cands) == 1:
            claim(cands[0], gff)
            continue
        if len(cands) > 1:
            problems.append((gff, f"stem '{stem}' matches {len(cands)} genomes"))
            continue

        # Tier 2: assembly accession (GCF_000001405.40_XYZ.gff <-> GCF_000001405.40.fna).
        acc = accession_of(gff.name)
        if acc:
            cands = genome_by_acc.get(acc.lower(), [])
            if len(cands) == 1:
                claim(cands[0], gff)
                continue
            if len(cands) > 1:
                problems.append((gff, f"accession '{acc}' matches {len(cands)} genomes"))
                continue

        # Tier 3: unambiguous one-sided prefix (longer name extends the shorter).
        cands = [
            g for g in genomes
            if (lambda gs: gs.startswith(stem_l) or stem_l.startswith(gs))(
                file_stem(g.name, FASTA_EXTS).lower())
        ]
        if len(cands) == 1:
            claim(cands[0], gff)
            continue
        if len(cands) > 1:
            problems.append(
                (gff, f"prefix '{stem}' is ambiguous across {len(cands)} genomes: "
                      f"{', '.join(sorted(c.name for c in cands)[:4])}")
            )
            continue

        problems.append((gff, f"no genome matches stem '{stem}'"))

    return assigned, problems


def count_usable_features(path: Path) -> Tuple[int, int]:
    """(indexable features, CDS features) -- cheap sanity check on an annotation.

    ``load_native_annotation_index`` only indexes gene/mRNA/transcript rows, and
    ``borrow_annotations`` wants >= 10 CDS. A GFF with neither is staged but
    warned about, because it will silently contribute nothing.
    """
    usable = 0
    cds = 0
    try:
        for feat in parse_gff(str(path)):
            ftype = feat.get("type")
            if ftype in USABLE_TYPES:
                usable += 1
            elif ftype == "CDS":
                cds += 1
    except Exception as exc:  # noqa: BLE001 - malformed GFF must not kill staging
        print(f"  [WARN] could not parse {path.name}: {exc}")
        return (0, 0)
    return (usable, cds)


def place(gff: Path, dest: Path) -> None:
    """Link the annotation into place; decompress .gz so every consumer reads it.

    ``borrow_annotations.find_annotated_target()`` only probes plain ``.gff``,
    so gzipped input is expanded rather than linked -- otherwise it would work
    for the search layer but silently not for annotation borrowing.
    """
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if gff.name.lower().endswith(".gz"):
        with gzip.open(gff, "rt") as src, open(dest, "w") as out:
            shutil.copyfileobj(src, out)
    else:
        os.symlink(os.path.realpath(gff), dest)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Stage target-genome GFFs next to their genome FASTAs (Pro mode)."
    )
    ap.add_argument("--genomes_dir", required=True,
                    help="Staged genomes directory (annotations are linked into it)")
    ap.add_argument("--gffs", nargs="*", default=[],
                    help="Target-genome GFF/GFF3/GTF files (plain or .gz)")
    ap.add_argument("--allow_unmatched", action="store_true",
                    help="Warn instead of failing when a GFF matches no genome")
    args = ap.parse_args()

    genomes_dir = Path(args.genomes_dir)
    if not genomes_dir.is_dir():
        print(f"[stage_target_gffs] ERROR: {genomes_dir} is not a directory")
        return 1

    gffs = collect_gffs(args.gffs)
    if not gffs:
        print("[stage_target_gffs] No target GFFs supplied — "
              "targets will be searched without native annotations.")
        return 0

    genomes = list_genomes(genomes_dir)
    if not genomes:
        print(f"[stage_target_gffs] ERROR: no genome FASTAs in {genomes_dir}")
        return 1

    print(f"[stage_target_gffs] Matching {len(gffs)} annotation(s) "
          f"to {len(genomes)} genome(s)...")

    assigned, problems = match_gffs_to_genomes(genomes, gffs)

    staged = 0
    for genome in genomes:
        gff = assigned.get(genome.name)
        if gff is None:
            continue
        dest_name = native_annotation_name(genome.name)
        if dest_name == genome.name:
            print(f"  [WARN] {genome.name}: cannot derive a distinct annotation "
                  f"name — skipping (rename the genome to carry an extension)")
            continue
        dest = genomes_dir / dest_name
        place(gff, dest)
        usable, cds = count_usable_features(dest)
        note = f"{usable} gene/mRNA, {cds} CDS"
        if usable == 0 and cds == 0:
            note += "  [WARN] no usable features — will contribute nothing"
        print(f"  {genome.name}  <-  {gff.name}   ({dest_name}; {note})")
        staged += 1

    missing = [g.name for g in genomes if g.name not in assigned]
    if missing:
        print(f"[stage_target_gffs] {len(missing)} genome(s) without annotation: "
              f"{', '.join(sorted(missing)[:6])}"
              f"{' ...' if len(missing) > 6 else ''}")

    if problems:
        print(f"[stage_target_gffs] {len(problems)} annotation(s) could not be assigned:")
        for gff, reason in problems:
            print(f"  - {gff.name}: {reason}")
        if not args.allow_unmatched:
            print("[stage_target_gffs] ERROR: refusing to continue with unassigned "
                  "annotations. Name each GFF after its genome (cow.fna -> cow.gff), "
                  "or pass --allow_unmatched to ignore them.")
            return 1

    print(f"[stage_target_gffs] Staged {staged} annotation(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
