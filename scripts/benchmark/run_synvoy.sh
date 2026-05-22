#!/usr/bin/env bash
# Run SynVoy in Pro Mode on the 18 benchmark genomes and emit calls.tsv.
# Reuses the existing ground-truth Apis mellifera home + the 5 Koludarov
# target genomes plus the 13 newly-fetched targets.
#
# Usage: bash scripts/benchmark/run_synvoy.sh

set -euo pipefail

QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
# Use the freshly-fetched RefSeq home (the original symlinks under
# local_data/ground_truth/melettin/home/ point to a deleted local_runs path).
HOME_FA="${HOME_FA:-local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.fna}"
HOME_GFF="${HOME_GFF:-local_data/benchmark/genomes/home/Apis_mellifera/Apis_mellifera.gff}"
GENOMES_DIR="${GENOMES_DIR:-local_data/benchmark/genomes}"
EXISTING_TARGETS="${EXISTING_TARGETS:-local_data/ground_truth/melettin/targets}"
OUTDIR="${OUTDIR:-benchmark_results/synvoy}"
RESULTS_DIR="${RESULTS_DIR:-results/melittin_benchmark}"

mkdir -p "${OUTDIR}"

# Stage all 17 target genomes (excluding Apis mellifera home) into a single dir
TARGETS_STAGED="${OUTDIR}/staged_targets"
mkdir -p "${TARGETS_STAGED}"

# Existing 5
for sp in Colletes_gigas Euglossa_dilemma Melipona_beecheii Tetragonula_carbonaria Xylocopa_violacea; do
    [[ -e "${EXISTING_TARGETS}/${sp}.fa" ]] && cp -n "${EXISTING_TARGETS}/${sp}.fa" "${TARGETS_STAGED}/"
done

# Newly fetched 12 (skip home)
for tier in tier1 tier2 tier3; do
    for species_dir in "${GENOMES_DIR}/${tier}"/*/; do
        sp=$(basename "${species_dir}")
        fa=$(ls "${species_dir}"*.fna 2>/dev/null | head -1)
        if [[ -n "${fa}" ]]; then
            cp -n "${fa}" "${TARGETS_STAGED}/${sp}.fna"
        fi
    done
done

n_targets=$(ls "${TARGETS_STAGED}/" | wc -l)
echo "Staged ${n_targets} target genomes"

if [[ "${n_targets}" -lt 17 ]]; then
    echo "WARNING: fewer than 17 targets staged. Continuing anyway." >&2
fi

# Run SynVoy. Nextflow needs Java; both live in synvoy_env. Calling the
# nextflow binary directly inherits the calling shell's PATH (which usually
# lacks java). `mamba run -n synvoy_env` properly activates the env so java
# is on PATH for nextflow's startup.
NEXTFLOW_ENV="${NEXTFLOW_ENV:-synvoy_env}"
if ! mamba run -n "${NEXTFLOW_ENV}" which nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow not found in conda env '${NEXTFLOW_ENV}'." >&2
    echo "  install: mamba install -n ${NEXTFLOW_ENV} -c bioconda nextflow openjdk" >&2
    echo "  or set NEXTFLOW_ENV to a different env name." >&2
    exit 1
fi
echo "Using nextflow from conda env: ${NEXTFLOW_ENV}"

echo "Running SynVoy (Pro Mode, this can take 6-12 hours)..."
mamba run -n "${NEXTFLOW_ENV}" nextflow run main.nf -profile standard --mode pro \
    --query "${QUERY}" \
    --home_genome "${HOME_FA}" \
    --home_gff "${HOME_GFF}" \
    --target_genomes "${TARGETS_STAGED}/*" \
    --n_flanking_genes 5 \
    --outdir "${RESULTS_DIR}" \
    -resume 2>&1 | tail -40 || true

# Parse calls from synvoy_report.json + per-genome regions BEDs.
python3 - "${RESULTS_DIR}" "${OUTDIR}/calls.tsv" <<'PY'
import json
import re
import sys
from pathlib import Path

results_dir = Path(sys.argv[1])
calls_path = Path(sys.argv[2])

# Locate the report
report = None
for cand in (results_dir / "synvoy_report.json",
             results_dir / "logs" / "GENERATE_REPORT__synvoy_report.json"):
    if cand.is_file():
        report = json.loads(cand.read_text())
        break

if report is None:
    print(f"WARNING: synvoy_report.json not found under {results_dir}", file=sys.stderr)
    report = {}

def canon_species(name):
    """Strip .fa / .fna / .fa.gz / .fna.gz to canonicalize species name.
    The Koludarov genomes carry the .fa suffix into the report's 'genome' field,
    while NCBI-fetched ones are extension-stripped — normalize so they merge.
    """
    return re.sub(r"\.f(n)?a(\.gz)?$", "", name)

ann_per = {}
for entry in report.get("annotations", {}).get("per_genome", []) or []:
    ann_per[canon_species(entry["genome"])] = entry

reg_per = {}
for entry in report.get("regions", {}).get("per_genome", []) or []:
    reg_per[canon_species(entry["genome"])] = entry

# Map each region BED file to its top region (best score, highest confidence)
regions_dir = results_dir / "regions"
top_region = {}  # species -> (chrom, start, end, strand, confidence_label, score)
if regions_dir.exists():
    for bed in regions_dir.glob("*.regions.bed"):
        sp = re.sub(r"\.f(n)?a(\.gz)?\.regions\.bed$", "", bed.name)
        best = None
        with bed.open() as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                chrom, start, end, name, score, strand = fields[:6]
                # name format: Species|Reg1_G10_CHIGH_S0.70
                m = re.search(r"_C([A-Z]+)_S([0-9.]+)", name)
                conf = m.group(1) if m else "?"
                sc = float(m.group(2)) if m else float(score or 0)
                rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}.get(conf, 0)
                if best is None or (rank, sc) > (best[5], best[6]):
                    best = (chrom, start, end, strand, conf, rank, sc)
        if best:
            top_region[sp] = best

# Union all species we know about (annotations or regions or BEDs)
species_seen = set(ann_per) | set(reg_per) | set(top_region)
if not species_seen:
    print("WARNING: no SynVoy data found — empty report and no region BEDs", file=sys.stderr)

calls = ["species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra"]

for sp in sorted(species_seen):
    a = ann_per.get(sp, {})
    r = reg_per.get(sp, {})
    tr = top_region.get(sp)

    resolved = a.get("resolved_goi_annotations", 0)
    ambiguous = a.get("ambiguous_goi_annotations", 0)
    total_goi = a.get("goi_annotations", 0)
    high_med_regions = (r.get("confidence_counts", {}).get("HIGH", 0)
                        + r.get("confidence_counts", {}).get("MEDIUM", 0))
    goi_anchor_regions = r.get("goi_anchor_regions", 0)

    # PRESENT = at least one resolved GOI annotation AND a HIGH/MEDIUM region
    # AMBIGUOUS = only LOW region (synteny didn't resolve cleanly), or
    #             only ambiguous-class annotations, or
    #             best_score below 0.30 (= region barely scored)
    # ABSENT = no GOI annotations at all
    #
    # The LOW-region rule rejects 3 v1 false positives (Polistes/Vespa LOW@~0.16,
    # Vollenhovia LOW@0.15) where SynVoy hit something but the region itself
    # had no flanking-gene support to anchor it. These are the cases the
    # paper also marks ABSENT/PSEUDOGENE.
    top_score = (tr[6] if tr else 0.0)
    top_conf = (tr[4] if tr else "-")
    has_high_med = (top_conf in ("HIGH", "MEDIUM")) or high_med_regions >= 1
    if resolved >= 1 and has_high_med and top_score >= 0.30:
        status = "PRESENT"
    elif resolved >= 1 or ambiguous >= 1 or total_goi >= 1:
        status = "AMBIGUOUS"
    else:
        status = "ABSENT"

    if tr:
        chrom, start, end, strand, conf_label, _, score = tr
    else:
        chrom = start = end = strand = "-"
        conf_label = "-"
        score = 0.0

    extra_bits = []
    if a:
        extra_bits.append(f"goi={total_goi}")
        extra_bits.append(f"resolved={resolved}")
        if ambiguous:
            extra_bits.append(f"ambig={ambiguous}")
    if r:
        extra_bits.append(f"regs={r.get('total_regions', 0)}")
        if score:
            extra_bits.append(f"best_score={score:.2f}")
    extra = ";".join(extra_bits) if extra_bits else "-"

    calls.append(
        f"{sp}\t-\t{status}\t{chrom}\t{start}\t{end}\t{strand}\t{conf_label}\t{extra}"
    )

calls_path.write_text("\n".join(calls) + "\n")
print(f"wrote {calls_path}")
PY

echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${OUTDIR}/calls.tsv" | sort | uniq -c
