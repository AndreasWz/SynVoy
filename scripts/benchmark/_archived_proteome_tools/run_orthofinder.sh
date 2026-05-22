#!/usr/bin/env bash
# OrthoFinder v3 on the 18 benchmark proteomes.
# Identifies whether melittin (query) clusters with each target's bee venom
# proteins, and emits calls.tsv.
#
# Usage: bash scripts/benchmark/run_orthofinder.sh

set -euo pipefail

PROTEOMES_DIR="${PROTEOMES_DIR:-local_data/benchmark/proteomes}"
QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
OUTDIR="${OUTDIR:-benchmark_results/orthofinder}"
# Lowered from 8 to 4 to fit a 16 GB laptop. OrthoFinder + DIAMOND fan out
# pairwise comparisons; 8 threads on this hardware caused SonicParanoid to
# OOM-kill in a sibling run.
THREADS="${THREADS:-4}"

mkdir -p "${OUTDIR}"

if ! command -v orthofinder >/dev/null 2>&1; then
    echo "ERROR: orthofinder not found. conda env create -f scripts/benchmark/competitors.yml" >&2
    exit 1
fi

# OrthoFinder needs all proteomes in one dir, one fasta per species.
# Inject the query as a "species" so we can ask which orthogroup it lands in.
WORK="${OUTDIR}/work"
mkdir -p "${WORK}/proteomes_with_query"
for f in "${PROTEOMES_DIR}"/*.faa; do
    cp "${f}" "${WORK}/proteomes_with_query/"
done
cp "${QUERY}" "${WORK}/proteomes_with_query/Apis_mellifera_query.faa"

# Resume detection: if a previous run left a WorkingDirectory with cached
# Blast*.txt.gz files, skip the (slow) blast phase and go straight to ortholog
# inference via OrthoFinder's `-b <WorkingDirectory>` mode.
BEST_WD=""
BEST_N=0
for wd in "${WORK}"/run_*/Results_*/WorkingDirectory; do
    [[ -d "${wd}" ]] || continue
    n=$(find "${wd}" -maxdepth 1 -name 'Blast*.txt.gz' 2>/dev/null | wc -l)
    if [[ "${n}" -gt "${BEST_N}" ]]; then
        BEST_N=${n}
        BEST_WD="${wd}"
    fi
done

# Detect already-completed orthogroup data: Orthogroups.txt is written before
# the (OOM-prone) species-tree-rooting step, so it usually survives crashes.
# If we have it, skip the OF re-run entirely and parse what's on disk.
EXISTING_OG_TXT=$(find "${WORK}" -name 'Orthogroups.txt' -size +1c 2>/dev/null | head -1)

if [[ -n "${EXISTING_OG_TXT}" ]]; then
    echo "Found existing Orthogroups.txt: ${EXISTING_OG_TXT}"
    echo "Skipping OrthoFinder re-run; parsing cached results directly."
    # Walk back up to the Results_* dir
    RESULTS_DIR=$(echo "${EXISTING_OG_TXT}" | sed 's|/Orthogroups/Orthogroups\.txt||' | xargs dirname | xargs dirname)
    # Actually we want the Results_* dir whose WorkingDirectory contains
    # SpeciesIDs.txt + SequenceIDs.txt.
    INNER_RESULTS_DIR=$(echo "${EXISTING_OG_TXT}" | sed 's|/Orthogroups/Orthogroups\.txt||')
    RESULTS_DIR="${INNER_RESULTS_DIR}"
elif [[ -n "${BEST_WD}" && "${BEST_N}" -ge 100 ]]; then
    echo "Resuming OrthoFinder from cached blast results:"
    echo "  WorkingDirectory: ${BEST_WD}"
    echo "  Cached blasts:    ${BEST_N}/324"
    orthofinder -b "${BEST_WD}" -t "${THREADS}" 2>&1 | tail -50 || true
    # Resumed run writes Results_* alongside the WorkingDirectory's parent
    RESULTS_DIR=$(find "$(dirname "${BEST_WD}")/.." -maxdepth 2 -name 'Results_*' -type d 2>/dev/null | sort | tail -1)
else
    # Fresh run with timestamped output dir.
    RUN_TAG="$(date +%Y%m%d_%H%M%S)"
    RUN_OUT="${WORK}/run_${RUN_TAG}"
    echo "Running OrthoFinder fresh (this can take 30-90 min at ${THREADS} threads)..."
    orthofinder -f "${WORK}/proteomes_with_query" -t "${THREADS}" -o "${RUN_OUT}" \
        2>&1 | tail -50 || true
    RESULTS_DIR=$(find "${RUN_OUT}" -maxdepth 2 -name 'Results_*' -type d | head -1)
fi

if [[ -z "${RESULTS_DIR}" ]]; then
    echo "ERROR: OrthoFinder produced no Results_* directory" >&2
    exit 1
fi
echo "OrthoFinder results in: ${RESULTS_DIR}"

# OrthoFinder may have crashed during the species-tree-rooting step (see e.g.
# FastTree OOM); the orthogroup clustering still completes earlier and writes
# Orthogroups.txt + the per-OG Sequences_ids/OG*.fa files. We parse those
# directly using SpeciesIDs.txt + SequenceIDs.txt for species attribution,
# which is robust against missing Orthogroups.tsv.

# SpeciesIDs.txt + SequenceIDs.txt live in the OUTER WorkingDirectory of the
# OrthoFinder run, while Orthogroups.txt + Sequences_ids/ live deeper in an
# inner Results_*/. Search the whole workspace to locate both.
INNER_WD=$(find "${WORK}" -name 'SpeciesIDs.txt' -size +1c -printf '%h\n' 2>/dev/null | head -1)
SEQUENCES_IDS_DIR=$(find "${WORK}" -type d -name 'Sequences_ids' 2>/dev/null | head -1)
SEQ_IDS="${INNER_WD}/SequenceIDs.txt"
SPECIES_IDS="${INNER_WD}/SpeciesIDs.txt"

if [[ ! -s "${SEQ_IDS}" || ! -s "${SPECIES_IDS}" || ! -d "${SEQUENCES_IDS_DIR}" ]]; then
    echo "ERROR: cannot locate SpeciesIDs.txt / SequenceIDs.txt / Sequences_ids/ under ${WORK}" >&2
    echo "  SpeciesIDs.txt: ${SPECIES_IDS}" >&2
    echo "  SequenceIDs.txt: ${SEQ_IDS}" >&2
    echo "  Sequences_ids dir: ${SEQUENCES_IDS_DIR}" >&2
    exit 1
fi
echo "Using ID maps from: ${INNER_WD}"
echo "Using OG sequences from: ${SEQUENCES_IDS_DIR}"

QUERY_HEADER=$(grep '^>' "${QUERY}" | head -1 | sed 's/^>//; s/ .*//')

python3 - "${SPECIES_IDS}" "${SEQ_IDS}" "${SEQUENCES_IDS_DIR}" "${QUERY_HEADER}" "${OUTDIR}/calls.tsv" <<'PY'
import sys
from pathlib import Path

species_ids_path, seq_ids_path, seq_dir, query_header, calls_path = sys.argv[1:]
seq_dir = Path(seq_dir)

# Map species_id -> species_name (canonical, _query/_home stripped)
species_id_to_name = {}
with open(species_ids_path) as fh:
    for line in fh:
        line = line.strip()
        if not line or ":" not in line:
            continue
        sid, fname = line.split(":", 1)
        sid = sid.strip()
        sp = fname.strip().replace(".faa", "").replace("_query", "").replace("_home", "")
        species_id_to_name[sid] = sp

# Build reverse: species_id -> set of (internal_id, gene_name)
# Plus: gene_name + species_id -> internal_id (for finding the home gene-Melt)
seq_id_to_gene = {}
gene_to_seq_ids = {}  # gene_name -> list[(species_id, internal_id)]
with open(seq_ids_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line or ":" not in line:
            continue
        iid, gene = line.split(":", 1)
        iid = iid.strip()
        gene = gene.strip()
        seq_id_to_gene[iid] = gene
        sp_id = iid.split("_")[0]
        gene_to_seq_ids.setdefault(gene, []).append((sp_id, iid))

# Find the home Apis_mellifera gene-Melt: species id whose name == "Apis_mellifera"
# (NOT the _query variant). The query was named gene-Melt too, so disambiguate by species.
home_sid = None
query_sid = None
for sid, sp in species_id_to_name.items():
    if sp == "Apis_mellifera":
        # Could be either the host proteome or the injected query single-seq file.
        # Distinguish: query species file contains exactly 1 sequence.
        # Look up via species's seq IDs.
        n_seqs = sum(1 for iid in seq_id_to_gene if iid.startswith(f"{sid}_"))
        if n_seqs <= 5:
            query_sid = sid
        else:
            home_sid = sid
print(f"home Apis_mellifera species_id={home_sid}, query species_id={query_sid}",
      file=sys.stderr)

# Find the internal ID of the home's gene-Melt
home_melt_iid = None
for sp_id, iid in gene_to_seq_ids.get(query_header, []):
    if sp_id == home_sid:
        home_melt_iid = iid
        break

if home_melt_iid is None:
    # Fall back to the query species
    for sp_id, iid in gene_to_seq_ids.get(query_header, []):
        if sp_id == query_sid:
            home_melt_iid = iid
            break

print(f"melittin internal ID: {home_melt_iid}", file=sys.stderr)

# Find the orthogroup file containing that internal ID
melittin_og = None
melittin_members = []  # list of (species_id, internal_id)
for og_file in sorted(seq_dir.glob("OG*.fa")):
    members = []
    found = False
    with og_file.open() as fh:
        for line in fh:
            if line.startswith(">"):
                iid = line[1:].strip().split()[0]
                sp_id = iid.split("_")[0]
                members.append((sp_id, iid))
                if iid == home_melt_iid:
                    found = True
    if found:
        melittin_og = og_file.stem
        melittin_members = members
        break

print(f"melittin orthogroup: {melittin_og}, n_members={len(melittin_members)}",
      file=sys.stderr)

# species_id -> list of internal_ids in the melittin OG
sid_to_members = {}
for sp_id, iid in melittin_members:
    sid_to_members.setdefault(sp_id, []).append(iid)

calls = ["species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra"]

# Skip the home and query "species" — they are not target genomes
skip_species = {"Apis_mellifera"}  # canonical name; matches both host and _query

# Build a per-canonical-species presence map (in case there are duplicate
# canonical names like Apis_mellifera vs Apis_mellifera_query).
canonical_present = {}  # canonical_name -> list[internal_id]
for sid, sp_canonical in species_id_to_name.items():
    if sp_canonical in skip_species:
        continue
    members_for_sid = sid_to_members.get(sid, [])
    if members_for_sid:
        canonical_present.setdefault(sp_canonical, []).extend(members_for_sid)
    else:
        canonical_present.setdefault(sp_canonical, [])

if melittin_og is None:
    # No OG with home_melt found — every target counted as ABSENT
    for sp in canonical_present:
        calls.append(f"{sp}\t-\tABSENT\t-\t-\t-\t-\torthogroup=not_found\thome_melt_iid_missing")
else:
    for sp, members in canonical_present.items():
        if members:
            first_gene = seq_id_to_gene.get(members[0], members[0])
            calls.append(
                f"{sp}\t-\tPRESENT\t-\t-\t-\t-\t"
                f"orthogroup={melittin_og}\tn_members={len(members)};first={first_gene}"
            )
        else:
            calls.append(f"{sp}\t-\tABSENT\t-\t-\t-\t-\torthogroup={melittin_og}\tno_member_in_OG")

with open(calls_path, "w") as fh:
    fh.write("\n".join(calls) + "\n")
print(f"wrote {calls_path}")
PY

echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${OUTDIR}/calls.tsv" | sort | uniq -c
