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
    
    # CRITICAL: Filter to only GOI sequences, then deduplicate by sequence content
    # Identical sequences from different query proteins should be collapsed
    python3 -c "
import sys
from collections import OrderedDict

# Pass 1: collect all GOI sequences
entries = []  # list of (name, seq)
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

# Filter to GOI only
goi_entries = [(n, s) for n, s in entries if 'GOI_' in n or 'GOI|' in n]

# Pass 2: deduplicate by sequence content
# For identical sequences, keep one representative name
seq_to_name = OrderedDict()  # seq -> first name seen
name_counts = {}  # handle name collisions

for name, seq in goi_entries:
    if seq in seq_to_name:
        continue  # skip identical sequence
    # Handle name collisions (same name, different sequence)
    if name in name_counts:
        name_counts[name] += 1
        unique_name = name + '_var' + str(name_counts[name])
    else:
        name_counts[name] = 1
        unique_name = name
    seq_to_name[seq] = unique_name

# Write deduplicated output
for seq, name in seq_to_name.items():
    sys.stdout.write('>' + name + chr(10))
    # Write sequence in 80-char lines
    for i in range(0, len(seq), 80):
        sys.stdout.write(seq[i:i+80] + chr(10))

print(f'Deduplicated: {len(goi_entries)} GOI entries -> {len(seq_to_name)} unique sequences', file=sys.stderr)
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
