process COMPUTE_TREE {
    tag "tree"
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(locus_id), path(fasta_files)

    output:
    tuple val(locus_id), path("*.nwk"), emit: tree

    script:
    """
    # Concatenate all fasta files
    cat ${fasta_files} > all_sequences.faa
    
    # Filter to GOI sequences: keep ALL GOI hits across all genomes (multiple
    # per genome are expected when expand_goi_similar finds paralogs).  Only
    # exon fragments and exact-header duplicates are removed.
    python3 -c "
import sys

# ---- Parse FASTA ----
entries = []
current_name = None
current_seq = []
for line in open('all_sequences.faa'):
    if line.startswith('>'):
        if current_name and current_seq:
            entries.append((current_name, ''.join(current_seq)))
        current_name = line.strip()[1:]
        current_seq = []
    else:
        current_seq.append(line.strip())
if current_name and current_seq:
    entries.append((current_name, ''.join(current_seq)))

# ---- Filter to GOI only, skip exon fragments ----
goi_entries = [(n, s) for n, s in entries
               if ('GOI_' in n or 'GOI|' in n) and '|exon_' not in n]

# ---- Deduplicate by exact header (keep first occurrence) ----
seen = set()
unique = []
for name, seq in goi_entries:
    if name not in seen:
        seen.add(name)
        unique.append((name, seq))

# ---- Write output ----
for name, seq in unique:
    sys.stdout.write('>' + name + chr(10))
    for i in range(0, len(seq), 80):
        sys.stdout.write(seq[i:i+80] + chr(10))

print(f'Tree filter: {len(entries)} total -> {len(goi_entries)} GOI -> {len(unique)} unique seqs (all per genome)', file=sys.stderr)
" > goi_only.faa 2> goi_dedup.log
    head -5 goi_dedup.log || true
    
    # Check if we have sequences
    count=\$(grep -c '^>' goi_only.faa 2>/dev/null || echo 0)
    echo "Filtered to \$count GOI sequences for tree"
    
    if [ "\$count" -lt 3 ]; then
        echo "Not enough GOI sequences (<3) for tree, creating placeholder"
        echo "(GOI_placeholder:0.0);" > ${locus_id}_tree.nwk
    else
        compute_tree.py \\
            --input goi_only.faa \\
            --output ${locus_id}_tree.nwk \\
            --threads ${task.cpus}
    fi
    """
}
