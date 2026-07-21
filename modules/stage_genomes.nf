process STAGE_GENOMES {
    tag "staging"
    
    input:
    path genomes
    // Pro-Mode target annotations (may be the NO_GFFS sentinel when unset).
    // Easy Mode gets these for free: fetch_related_genomes.py already writes
    // {accession}.gff next to {accession}.fna, which is exactly where
    // iterative_search_runner.find_native_annotation_path() probes.
    path target_gffs

    output:
    path "staged_genomes", emit: dir
    path "species_mapping.tsv", emit: species_map
    
    script:
    """
    mkdir staged_genomes
    # Link (don't copy) the input genomes into one directory. Nextflow already
    # staged them into this task dir, so the old `cp -f -L -r` dereferenced those
    # stage-symlinks into full byte copies — duplicating tens of GB of FASTA and,
    # on a slow shared filesystem, blowing past the per-process wall-clock limit
    # (the reported "staging exceeded the 1h cutoff"). `readlink -f` resolves to
    # the real user file so the link is stable; downstream consumers use
    # `find -L`, so symlinks are transparent to them.
    for src in $genomes; do
        if [ -d "\$src" ]; then
            # Defensive: a directory slipped through (resolveTargetGenomeFiles
            # normally flattens folders to files). Link each FASTA inside it.
            for inner in "\$src"/*; do
                [ -f "\$inner" ] && ln -sf "\$(readlink -f "\$inner")" "staged_genomes/\$(basename "\$inner")"
            done
        elif [ -f "\$src" ]; then
            ln -sf "\$(readlink -f "\$src")" "staged_genomes/\$(basename "\$src")"
        fi
    done

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

    # Link user-supplied target annotations in AFTER the species map is built —
    # the loop above walks staged_genomes/* untyped, so a .gff present at that
    # point would be recorded as a bogus "genome" row.
    #
    # stage_target_gffs.py names each annotation Path(genome).with_suffix('.gff')
    # (cow.fna -> cow.gff). That is deliberate: cluster_regions.nf resolves the
    # genome with `find -name "<genome_name>*"` WITHOUT an extension filter, and
    # <genome_name> carries the extension, so cow.fna.gff would be returned in
    # place of the genome. cow.gff is invisible to that glob and is also the
    # first path find_native_annotation_path() probes.
    stage_target_gffs.py \\
        --genomes_dir staged_genomes \\
        ${params.allow_unmatched_target_gffs ? '--allow_unmatched' : ''} \\
        --gffs $target_gffs
    """
}
