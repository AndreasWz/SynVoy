#!/bin/bash
# Re-fetch the INSECT demo inputs into local_data/demo_insect/ (gitignored, ~600 MB).
# Needs the NCBI `datasets` CLI and network access.
#
# Why Drosophila and not something smaller: the two genuinely smaller annotated
# insect genomes are unusable for a synteny demo. Belgica antarctica (90 Mb, the
# smallest insect genome known) has NO annotation at all; Clunio marinus (85 Mb)
# is annotated but sits in 24,952 scaffolds, which fragments the flanking signal.
# Drosophila is the smallest insect clade that is chromosome-level AND annotated
# across several species, which is what a synteny demo actually needs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

DEMO=local_data/demo_insect
mkdir -p "$DEMO"/{home,genomes,target_gffs,query,_dl}

declare -A NAME=(
  [GCF_000001215.4]=Drosophila_melanogaster   # home,    144 Mb, gold-standard annotation
  [GCF_016746395.2]=Drosophila_simulans       # target,  132 Mb,  97 seqs, ~5 Mya
  [GCF_016746365.2]=Drosophila_yakuba         # target,  148 Mb,  64 seqs, ~10 Mya
  [GCF_009870125.1]=Drosophila_pseudoobscura  # target,  163 Mb,  71 seqs, ~25 Mya
)

for acc in "${!NAME[@]}"; do
    n=${NAME[$acc]}
    echo "==> $n ($acc)"
    datasets download genome accession "$acc" --include genome,gff3 \
        --filename "$DEMO/_dl/$acc.zip" --no-progressbar
    rm -rf "$DEMO/_dl/x" && mkdir -p "$DEMO/_dl/x"
    unzip -q "$DEMO/_dl/$acc.zip" -d "$DEMO/_dl/x"
    fna=$(find "$DEMO/_dl/x" -name "*_genomic.fna" | head -1)
    gff=$(find "$DEMO/_dl/x" -name "genomic.gff" | head -1)
    if [ "$acc" = "GCF_000001215.4" ]; then
        cp "$fna" "$DEMO/home/$n.fna"; cp "$gff" "$DEMO/home/$n.gff"
    else
        cp "$fna" "$DEMO/genomes/$n.fna"; cp "$gff" "$DEMO/target_gffs/$n.gff"
    fi
done
rm -rf "$DEMO/_dl"

# Query proteins (UniProt, reviewed D. melanogaster entries)
curl -sL "https://rest.uniprot.org/uniprotkb/P36192.fasta" -o "$DEMO/query/Defensin.faa"  # AMP, 92 aa, fast-evolving
curl -sL "https://rest.uniprot.org/uniprotkb/P00334.fasta" -o "$DEMO/query/Adh.faa"       # alcohol dehydrogenase, 256 aa, the control

echo
echo "Done. Inputs in $DEMO ($(du -sh "$DEMO" | cut -f1))"
