process STAGE_GENOMES {
    tag "staging"
    
    input:
    path genomes
    
    output:
    path "staged_genomes", emit: dir
    path "species_mapping.tsv", emit: species_map
    
    script:
    """
    mkdir staged_genomes
    # Copy files into directory
    # Use cp -r to handle potential directories or multiple files
    cp -f -L -r $genomes staged_genomes/

    # Extract organism binomial from FASTA headers to build a species map.
    # Typical NCBI header: >NC_039506.1 Ooceraea biroi isolate clonal line ...
    # We use the file basename (minus extension) as the accession key and
    # grab the first two words after the sequence ID as genus + species.
    > species_mapping.tsv
    for f in staged_genomes/*; do
        [ -f "\$f" ] || continue
        acc=\$(basename "\$f")
        acc="\${acc%%.*}"          # strip all extensions  e.g. GCF_003672135
        # Re-add the version suffix if the original filename had one (GCF_xxx.1)
        bn=\$(basename "\$f")
        ver=\$(echo "\$bn" | grep -oP '\\.[0-9]+(?=\\.[a-zA-Z])' || true)
        acc="\${acc}\${ver}"

        # Read first header line (handle gzipped files too)
        case "\$f" in
            *.gz) hdr=\$(zcat "\$f" | head -1) ;;
            *)    hdr=\$(head -1 "\$f") ;;
        esac

        # Extract organism: skip the sequence accession (first word), take next two words
        species=\$(echo "\$hdr" | sed 's/^>[^ ]* //' | awk '{print \$1, \$2}')
        if [ -n "\$species" ] && [ "\$species" != " " ]; then
            printf '%s\\t%s\\tgenome\\n' "\$acc" "\$species" >> species_mapping.tsv
        fi
    done
    """
}
