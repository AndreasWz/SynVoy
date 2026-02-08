process PREPARE_INITIAL_DB {
    tag "prepare_db_${locus_id}"
    label 'process_low'

    input:
    tuple val(locus_id), path(flanking_faa)
    path goi_exons  // From ANNOTATE_GOI: contains full GOI + individual exon sequences

    output:
    tuple val(locus_id), path("initial_db_${locus_id}.faa"), emit: db

    script:
    """
    #!/usr/bin/env python3
    
    # Build the initial search database from:
    # 1. Flanking genes (from GFF or Prodigal prediction)
    # 2. GOI exon sequences (from ANNOTATE_GOI - real exons, not arbitrary fragments)
    # 3. Fallback fragments (halves/thirds) only if no real exons were found
    
    import sys
    import os
    
    sys.path.insert(0, "${projectDir}/bin")
    from sequence_utils import parse_fasta, write_fasta
    from fragment_query import generate_fragments
    
    output_faa = "initial_db_${locus_id}.faa"
    all_records = []
    
    # 1. Add all flanking genes
    print("Loading flanking genes from ${flanking_faa}...")
    for header, clean_id, seq in parse_fasta("${flanking_faa}"):
        all_records.append((clean_id, seq))
    flanking_count = len(all_records)
    print(f"  Loaded {flanking_count} flanking genes")
    
    # 2. Add GOI sequences (full protein + individual exons from ANNOTATE_GOI)
    print("Loading GOI exon sequences from ${goi_exons}...")
    goi_records = list(parse_fasta("${goi_exons}"))
    if not goi_records:
        print("ERROR: No sequences found in GOI exons file!", file=sys.stderr)
        sys.exit(1)
    
    goi_full_count = 0
    exon_count = 0
    tandem_count = 0
    full_goi_seq = None
    full_goi_id = None
    
    for header, clean_id, seq in goi_records:
        all_records.append((clean_id, seq))
        
        if '|exon_' in clean_id:
            exon_count += 1
            print(f"  Added exon: {clean_id} ({len(seq)} aa)")
        elif clean_id.startswith('GOI_copy_'):
            tandem_count += 1
            print(f"  Added tandem copy: {clean_id} ({len(seq)} aa)")
        else:
            goi_full_count += 1
            full_goi_seq = seq
            full_goi_id = clean_id
            print(f"  Added full GOI: {clean_id} ({len(seq)} aa)")
    
    # 3. Fallback: Generate arbitrary fragments ONLY if no real exons or tandem copies found
    fragment_count = 0
    if exon_count == 0 and tandem_count == 0 and full_goi_seq:
        print("  No real exons found, generating fallback fragments...")
        min_fragment_size = 20  # amino acids
        fragments = generate_fragments(full_goi_seq, full_goi_id, min_size=min_fragment_size)
        
        for frag_id, frag_seq, frag_desc in fragments:
            if "fragment_type=full" not in frag_desc:
                all_records.append((frag_id, frag_seq))
                fragment_count += 1
        
        print(f"  Generated {fragment_count} fallback fragments")
    elif exon_count > 0:
        print(f"  Using {exon_count} real exon sequences (no arbitrary fragments needed)")
    elif tandem_count > 0:
        print(f"  Using {tandem_count} tandem copies as queries (no arbitrary fragments needed)")
    
    # 4. Write combined database
    write_fasta(all_records, output_faa)
    
    print(f"\\nInitial database created: {output_faa}")
    print(f"  Total sequences: {len(all_records)}")
    print(f"  - Flanking genes: {flanking_count}")
    print(f"  - GOI (full): {goi_full_count}")
    print(f"  - GOI exons: {exon_count}")
    print(f"  - GOI tandem copies: {tandem_count}")
    print(f"  - Fallback fragments: {fragment_count}")
    """
}
