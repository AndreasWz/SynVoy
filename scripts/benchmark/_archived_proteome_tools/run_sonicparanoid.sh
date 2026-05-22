#!/usr/bin/env bash
# SonicParanoid2 on the 18 benchmark proteomes.
# Same orthogroup-membership question as OrthoFinder.
#
# Usage: bash scripts/benchmark/run_sonicparanoid.sh

set -euo pipefail

PROTEOMES_DIR="${PROTEOMES_DIR:-local_data/benchmark/proteomes}"
QUERY="${QUERY:-local_data/ground_truth/melettin/query_melittin.faa}"
OUTDIR="${OUTDIR:-benchmark_results/sonicparanoid}"
# Lowered from 8 to 4 to stay within 16 GB RAM + 4 GB swap. SP+DIAMOND
# concurrently fan out one worker per pair (~ N^2/2 = 153 pairs at N=18),
# and at 8 threads on a 16 GB laptop it OOM-killed the prior run.
THREADS="${THREADS:-4}"

mkdir -p "${OUTDIR}"

if ! command -v sonicparanoid >/dev/null 2>&1; then
    echo "ERROR: sonicparanoid not found. conda env create -f scripts/benchmark/competitors.yml" >&2
    exit 1
fi

WORK="${OUTDIR}/work"
mkdir -p "${WORK}/input"
for f in "${PROTEOMES_DIR}"/*.faa; do
    cp "${f}" "${WORK}/input/"
done
cp "${QUERY}" "${WORK}/input/Apis_mellifera_query.faa"

# SonicParanoid refuses to use a project name whose runs/<NAME> dir already
# exists. Auto-rotate: pick the next available name. Per-pair alignments live
# in sp_results/orthologs_db/ and ARE reused across runs even with new names.
SP_RUN_NAME="melittin_benchmark"
i=2
while [[ -d "${WORK}/sp_results/runs/${SP_RUN_NAME}" ]]; do
    SP_RUN_NAME="melittin_benchmark_v${i}"
    i=$((i + 1))
done
echo "Using SonicParanoid project name: ${SP_RUN_NAME}"

echo "Running SonicParanoid with ${THREADS} threads (this can take 30-90 min on 16 GB)..."
echo "Cached pairwise alignments in orthologs_db/ are reused across project names."
# -ow forces SP to recompute even when it detects identical inputs to a previous
# run (it hashes input and refuses to re-run otherwise). Without -ow we get:
# "WARNING: the input is the same as the one used in the last run, ..."
sonicparanoid -i "${WORK}/input" -o "${WORK}/sp_results" \
    -p "${SP_RUN_NAME}" -t "${THREADS}" -ow 2>&1 | tail -50 || true

# SonicParanoid emits ortholog_groups.tsv
OG_TSV=$(find "${WORK}/sp_results/runs/${SP_RUN_NAME}/ortholog_groups" -name 'ortholog_groups.tsv' 2>/dev/null | head -1)
if [[ -z "${OG_TSV}" || ! -s "${OG_TSV}" ]]; then
    echo "ERROR: SonicParanoid ortholog_groups.tsv not found under ${WORK}/sp_results/" >&2
    echo "  This usually means SP was killed (OOM) before the inference step finished." >&2
    echo "  Check: $(du -sh ${WORK}/sp_results/runs/melittin_benchmark/aux 2>/dev/null) of alignment matrices exist;" >&2
    echo "  re-running with THREADS=2 should resume and complete in ~10-20 min." >&2
    exit 1
fi

QUERY_HEADER=$(grep '^>' "${QUERY}" | head -1 | sed 's/^>//; s/ .*//')

python3 - "${OG_TSV}" "${QUERY_HEADER}" "${OUTDIR}/calls.tsv" <<'PY'
import sys
og_tsv, query_header, calls_path = sys.argv[1:]

with open(og_tsv) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    # SonicParanoid columns: group_id, group_size, sp1, sp1_genes, sp2, sp2_genes, ...
    sp_cols = []
    for i, h in enumerate(header):
        if h.endswith(".faa") or h in ("group_id", "group_size"):
            continue
        # Could be 'species_name' before the corresponding '_genes' col
    # Simpler approach: look at every other column starting from idx 4 (per docs)
    # group_id\tgroup_size\tsp1_id\tsp1_genes\tsp2_id\tsp2_genes\t...
    # But header layout varies; use index-based parse: from col 4 onward, pairs
    # are (species, members).
    # Find the line containing the query
    melittin_row = None
    for line in fh:
        if query_header in line:
            melittin_row = line.rstrip("\n").split("\t")
            break

calls = ["species\taccession\tcalled_status\tlocus_chrom\tlocus_start\tlocus_end\tstrand\tconfidence\textra"]

if melittin_row is None:
    # Query is singleton or absent
    seen = set()
    for h in header:
        if h.endswith(".faa"):
            sp = h.replace(".faa", "").replace("_query", "").replace("_home", "")
            if sp == "Apis_mellifera":
                continue
            if sp in seen:
                continue
            seen.add(sp)
            calls.append(f"{sp}\t-\tABSENT\t-\t-\t-\t-\torthogroup=singleton\tquery_unclustered")
else:
    group_id = melittin_row[0]
    # Pair up species cols with members; SonicParanoid puts members in
    # columns whose header ends with .faa
    seen = set()
    for i, h in enumerate(header):
        if not h.endswith(".faa"):
            continue
        sp = h.replace(".faa", "").replace("_query", "").replace("_home", "")
        if sp in ("Apis_mellifera",):
            continue
        if sp in seen:
            continue
        seen.add(sp)
        members = melittin_row[i] if i < len(melittin_row) else "*"
        if members and members != "*":
            calls.append(f"{sp}\t-\tPRESENT\t-\t-\t-\t-\torthogroup={group_id}\tmembers={members}")
        else:
            calls.append(f"{sp}\t-\tABSENT\t-\t-\t-\t-\torthogroup={group_id}\tno_member_in_OG")

with open(calls_path, "w") as fh:
    fh.write("\n".join(calls) + "\n")
print(f"wrote {calls_path}")
PY

echo "Summary:"
awk -F'\t' 'NR>1 {print $3}' "${OUTDIR}/calls.tsv" | sort | uniq -c
