#!/bin/bash
# =============================================================================
#  run_uniprot.sh — LIVE DEMO: audience picks a UniProt ID, SynVoy runs on it
# =============================================================================
#
#     demo/run_uniprot.sh P36192
#
#  Fetches that protein from UniProt, then runs SynVoy in PRO mode against three
#  pre-downloaded Drosophila genomes (no genome download at demo time).
#
#  The ID must be a *Drosophila melanogaster* protein (OX=7227) — the home genome
#  is D. melanogaster, so a protein from any other organism has no home locus to
#  anchor on. The script refuses non-Dmel IDs rather than producing garbage.
#
#  Known-good IDs to fall back on if the audience's pick misbehaves. All four
#  verified to resolve, to be D. melanogaster, and to be single-copy:
#     P36192  Def     92 aa   Defensin, antimicrobial peptide  (VALIDATED, ~8 min)
#     P61851  Sod1   153 aa   Cu/Zn superoxide dismutase
#     P00334  Adh    256 aa   alcohol dehydrogenase
#     P20477  Gs1    399 aa   glutamine synthetase 1
#  Avoid P10987 (Act5C) — Drosophila has six actins, so it pulls in the family.
# =============================================================================
set -o pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

ACC="${1:?usage: demo/run_uniprot.sh <UNIPROT_ID>   e.g. demo/run_uniprot.sh P36192}"
shift 2>/dev/null

DEMO=local_data/demo_insect
[ -d "$DEMO/genomes" ] || { echo "ERROR: $DEMO missing. Run demo/fetch_data_insect.sh first."; exit 1; }

# A VS Code-activated .venv shadows the conda env inside Nextflow tasks and makes
# parasail (Smith-Waterman) disappear. Strip it so the demo behaves identically
# from a terminal and from the IDE.
if [ -n "${VIRTUAL_ENV:-}" ]; then
    echo "[demo] dropping active virtualenv $VIRTUAL_ENV"
    PATH="$(echo "$PATH" | tr ':' '\n' | grep -v "^${VIRTUAL_ENV}/bin$" | paste -sd:)"
    unset VIRTUAL_ENV; export PATH
fi

mkdir -p "$DEMO/query"
FAA="$DEMO/query/${ACC}.faa"

echo "==> fetching $ACC from UniProt"
curl -fsSL --max-time 30 "https://rest.uniprot.org/uniprotkb/${ACC}.fasta" -o "$FAA" \
    || { echo "ERROR: could not fetch $ACC from UniProt (bad ID, or no network)."; exit 1; }
[ -s "$FAA" ] || { echo "ERROR: UniProt returned an empty record for $ACC — check the accession."; rm -f "$FAA"; exit 1; }

HDR=$(head -1 "$FAA")
LEN=$(grep -v '^>' "$FAA" | tr -d '\n \t' | wc -c)
SYM=$(sed -nE 's/.*GN=([^ ]+).*/\1/p' <<<"$HDR"); SYM=${SYM:-$ACC}
echo "    $HDR"
echo "    gene symbol: $SYM | length: ${LEN} aa"

# Guard 1: wrong organism. The home genome is D. melanogaster (taxid 7227).
if ! grep -q "OX=7227" <<<"$HDR"; then
    echo
    echo "REFUSING TO RUN: $ACC is not Drosophila melanogaster (no OX=7227 in the header)."
    echo "The home genome for this demo is D. melanogaster, so there would be no home"
    echo "locus to anchor the synteny search on. Pick a Dmel protein — e.g. search"
    echo "UniProt for the gene name with 'AND organism_id:7227'."
    exit 1
fi

# Guard 2: too short. NORMALIZE_QUERY enforces min_query_length BEFORE the preset
# logic can adjust it, so a <30 aa query must be handled at launch (docs/TODO.md
# §1f caveat). Lower the floor rather than aborting.
EXTRA=()
if [ "$LEN" -lt 30 ]; then
    echo "    [note] ${LEN} aa is below the 30 aa default floor — passing --min_query_length ${LEN}"
    EXTRA+=(--min_query_length "$LEN")
fi

# Guard 3: runtime. Search cost scales with query length; the validated 92 aa run is
# ~6 min. This does NOT block — it warns, so you can pick again before committing to
# a long wait in front of an audience.
if [ "$LEN" -gt 1200 ]; then
    echo
    echo "    [WARNING] ${LEN} aa is a long query — expect well over 10 minutes."
    echo "              Ctrl-C now and pick a shorter one if you are on the clock."
    echo "              (P36192 Def 92 aa is the validated ~6 min run.)"
    sleep 5
fi

OUTDIR="results/demo_${SYM}"

# Explicit two-genome list, NOT a glob over $DEMO/genomes/. D. simulans and
# D. yakuba are still on disk as spares; naming the targets here keeps the demo
# at three genomes without moving files around.
TARGETS="$DEMO/genomes/Drosophila_pseudoobscura.fna,$DEMO/genomes/Anopheles_gambiae.fna"
TGFFS="$DEMO/target_gffs/Drosophila_pseudoobscura.gff,$DEMO/target_gffs/Anopheles_gambiae.gff"

echo
echo "=============================================================="
echo " SynVoy live demo — $SYM ($ACC, ${LEN} aa)"
echo " home:    Drosophila melanogaster    (144 Mb)"
echo " target:  Drosophila pseudoobscura   (163 Mb, ~25 Mya)"
echo " target:  Anopheles gambiae          (264 Mb, ~250 Mya, outside Drosophila)"
echo " output:  $OUTDIR"
echo "=============================================================="

./run_synvoy.sh \
    -profile low_mem \
    -c demo/demo_insect.config \
    --mode pro \
    --query          "$FAA" \
    --home_genome    "$DEMO/home/Drosophila_melanogaster.fna" \
    --home_gff       "$DEMO/home/Drosophila_melanogaster.gff" \
    --home_species   "Drosophila melanogaster" \
    --target_genomes "$TARGETS" \
    --target_gffs    "$TGFFS" \
    --outdir         "$OUTDIR" \
    "${EXTRA[@]}" "$@"

echo
echo "Figure:  $OUTDIR/synteny_block_locus_1_anchor_grid.html"
echo "Report:  $OUTDIR/synvoy_report.json"
