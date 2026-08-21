#!/bin/bash
# Re-fetch the demo inputs into local_data/demo/ (gitignored, ~76 MB).
# Needs the NCBI `datasets` CLI and network access.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

DEMO=local_data/demo
mkdir -p "$DEMO"/{home,genomes,target_gffs,query,_dl}

declare -A NAME=(
  [GCF_000146045.2]=Saccharomyces_cerevisiae     # home, S288C R64
  [GCF_000237345.1]=Naumovozyma_castellii        # target, post-WGD
  [GCF_000002515.2]=Kluyveromyces_lactis         # target, pre-WGD
  [GCF_000142805.1]=Lachancea_thermotolerans     # target, pre-WGD
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
    if [ "$acc" = "GCF_000146045.2" ]; then
        cp "$fna" "$DEMO/home/$n.fna"; cp "$gff" "$DEMO/home/$n.gff"
    else
        cp "$fna" "$DEMO/genomes/$n.fna"; cp "$gff" "$DEMO/target_gffs/$n.gff"
    fi
done
rm -rf "$DEMO/_dl"

# Query proteins (UniProt, reviewed S. cerevisiae entries)
curl -sL "https://rest.uniprot.org/uniprotkb/D6VTK4.fasta" -o "$DEMO/query/STE2.faa"  # pheromone receptor, 431 aa
curl -sL "https://rest.uniprot.org/uniprotkb/P00560.fasta" -o "$DEMO/query/PGK1.faa"  # phosphoglycerate kinase, 416 aa

echo
echo "Done. Inputs in $DEMO ($(du -sh "$DEMO" | cut -f1))"
