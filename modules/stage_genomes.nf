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

    # Build species_mapping.tsv (key = filename-derived accession, value = species).
    # Strategy: prefer organism binomial from FASTA header (NCBI-style:
    #   ">NC_039506.1 Ooceraea biroi isolate ..."). Fall back to the filename
    # itself (with underscores → spaces) when the header has no binomial —
    # this is the common Pro-Mode case where users pass files named after the
    # species (e.g. Colletes_gigas.fa with header ">WUUM01000001.1").
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

        # Extract organism: skip the sequence accession (first word), take next two words.
        species=\$(echo "\$hdr" | sed 's/^>[^ ]* //' | awk '{print \$1, \$2}')
        # Trim whitespace
        species=\$(echo "\$species" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*\$//')

        # Reject bogus extractions: starts with '>' (no organism in header), or
        # is a single token (also no organism — awk took just the seq ID).
        if [ -z "\$species" ] || [ "\${species:0:1}" = ">" ] || ! echo "\$species" | grep -q ' '; then
            # Fall back to the filename (acc), converting underscores to spaces so
            # 'Colletes_gigas' becomes the binomial 'Colletes gigas'.
            species=\$(echo "\$acc" | tr '_' ' ')
        fi

        if [ -n "\$species" ] && [ "\$species" != " " ]; then
            printf '%s\\t%s\\tgenome\\n' "\$acc" "\$species" >> species_mapping.tsv
        fi
    done
    """
}
