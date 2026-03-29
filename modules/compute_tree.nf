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
    
    # CRITICAL: Filter to GOI sequences, then keep ONE representative per genome.
    # The tree is only used for genome-level ordering & colouring, so we do not
    # need hundreds of GOI copies – one longest rep per genome is sufficient.
    python3 -c "
import sys, re
from collections import OrderedDict

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

# ---- Extract genome_id from header ----
_gcf_re = re.compile(r'(GC[FA]_[0-9]+_[0-9]+)')
def genome_id(name):
    m = _gcf_re.search(name)
    if m:
        parts = m.group(1).split('_')          # ['GCF','012345','6']
        return f'{parts[0]}_{parts[1]}.{parts[2]}'
    # Non-GCF genomes: extract the genome name encoded between '|' and '_b<N>_'
    # e.g. 'GOI_Melt|Colletes_gigas_fa_b0_l1_fallback' -> 'Colletes_gigas_fa'
    m2 = re.search(r'[|](.+?)_b[0-9]+_', name)
    if m2:
        return m2.group(1)
    return 'home'

# ---- Keep longest representative per genome ----
best = {}  # genome_id -> (name, seq)
for name, seq in goi_entries:
    gid = genome_id(name)
    if gid not in best or len(seq) > len(best[gid][1]):
        best[gid] = (name, seq)

# ---- Write output ----
for gid, (name, seq) in best.items():
    sys.stdout.write('>' + name + chr(10))
    for i in range(0, len(seq), 80):
        sys.stdout.write(seq[i:i+80] + chr(10))

print(f'Tree filter: {len(entries)} total -> {len(goi_entries)} GOI -> {len(best)} representative seqs (1 per genome)', file=sys.stderr)
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
