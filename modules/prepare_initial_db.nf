process PREPARE_INITIAL_DB {
    tag "prepare_db_${locus_id}"
    label 'process_low'

    input:
    tuple val(locus_id), path(flanking_faa)
    path query_gene

    output:
    tuple val(locus_id), path("initial_db_${locus_id}.faa"), emit: db

    script:
    """
    #!/usr/bin/env python3
    
    # This critical fix ensures the Query Gene of Interest (GOI) is included
    # in the iterative search database, along with its fragments and exons
    
    import sys
    import os
    
    # Import our sequence utilities
    sys.path.insert(0, "${projectDir}/bin")
    from sequence_utils import parse_fasta, write_fasta
    from fragment_query import generate_fragments
    
    output_faa = "initial_db_${locus_id}.faa"
    all_records = []
    
    # 1. Add all flanking genes
    print("Loading flanking genes from ${flanking_faa}...")
    for header, clean_id, seq in parse_fasta("${flanking_faa}"):
        all_records.append((clean_id, seq))
    print(f"  Loaded {len(all_records)} flanking genes")
    
    # 2. Add the Query Gene of Interest (GOI) - CRITICAL FIX!
    print("Loading query gene from ${query_gene}...")
    goi_records = list(parse_fasta("${query_gene}"))
    if not goi_records:
        print("ERROR: No sequences found in query gene file!", file=sys.stderr)
        sys.exit(1)
    
    goi_count = 0
    for header, clean_id, seq in goi_records:
        # Mark as GOI for special handling in iterative search
        goi_id = f"GOI_{clean_id}"
        all_records.append((goi_id, seq))
        goi_count += 1
        
        # 3. Generate systematic fragments (halves, thirds, quarters)
        # This allows detection of partial genes, truncations, pseudogenes
        min_fragment_size = 20  # amino acids
        fragments = generate_fragments(seq, goi_id, min_size=min_fragment_size)
        
        for frag_id, frag_seq, frag_desc in fragments:
            # Skip the full-length duplicate (already added)
            if "fragment_type=full" not in frag_desc:
                all_records.append((frag_id, frag_seq))
        
        print(f"  Added GOI: {goi_id} ({len(seq)} aa)")
        print(f"  Generated {len(fragments)-1} fragments for progressive search")
    
    # 4. Write combined database
    write_fasta(all_records, output_faa)
    
    print(f"\\nInitial database created: {output_faa}")
    print(f"  Total sequences: {len(all_records)}")
    print(f"  - Flanking genes: {len(all_records) - goi_count - len(fragments) + 1}")
    print(f"  - Query genes (GOI): {goi_count}")
    print(f"  - Query fragments: {len(fragments) - 1}")
    print(f"\\nCRITICAL: GOI is now included in iterative search!")
    """
}
