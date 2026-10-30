#!/usr/bin/env python3
"""
smith_waterman_search.py - Smith-Waterman local alignment for GOI search

Uses parasail library for vectorized Smith-Waterman alignment.
This provides more sensitive local alignment than MMseqs2 for divergent sequences.

Hits are written as BLAST m8 with real Karlin-Altschul bit scores and E-values for the
window's six-frame search space (BLOSUM62, BLAST gap costs 11/1 by default).

Every piece of a translated frame is scored first without a traceback (parasail's score-only
kernel, several times faster); the traceback runs only where that score can give a hit the
runner keeps (min_reportable_score). The hits with E <= 1 are the same as with a traceback
on every piece; the time per window is what the inputs make it, so the runner's guard
(--sw_timeout_seconds) can sit far above it. A failed alignment ends the program with an
error (exit 1) instead of skipping the piece.

Usage:
    python smith_waterman_search.py --query goi.faa --target region.fna \\
        --output hits.tsv --matrix BLOSUM62 --gap_open 11 --gap_extend 1
"""

import argparse
import math
import sys
import os
try:
    from sequence_utils import parse_fasta, translate, reverse_complement, set_genetic_code
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import parse_fasta, translate, reverse_complement, set_genetic_code

try:
    import parasail
    HAS_PARASAIL = True
except ImportError:
    HAS_PARASAIL = False
    print("WARNING: parasail is not installed: the Smith-Waterman search can use ssearch36 only.",
          file=sys.stderr)


# parasail picks its SIMD code at run time, the first of these the CPU supports.
SIMD_ORDER = ('avx2', 'sse41', 'sse2', 'altivec', 'neon')


def backend_info(method='auto'):
    """What ``--method`` runs on this machine, for the run's record (fix-phase item 18):
    {'backend', 'parasail_version', 'simd'}. 'parasail' is the main path (a score-only pass,
    then a traceback where a kept hit can come out); 'ssearch36' with parasail installed runs
    run_parasail_sw ('parasail_best_hit_rounds'); without parasail it is ssearch36 itself."""
    if method == 'auto':
        method = 'parasail' if HAS_PARASAIL else 'ssearch36'
    if method == 'ssearch36':
        backend = 'parasail_best_hit_rounds' if HAS_PARASAIL else 'ssearch36'
    else:
        backend = 'parasail'
    version = simd = ''
    if HAS_PARASAIL:
        # the Python binding's version and the C library's, which holds the alignment code
        version = str(getattr(parasail, '__version__', '') or '')
        try:
            version += f" (libparasail {'.'.join(str(v) for v in parasail.version())})"
        except Exception:      # an older binding without version()
            pass
        simd = next((isa for isa in SIMD_ORDER
                     if getattr(parasail, f'can_use_{isa}', lambda: False)()), 'none')
    return {'backend': backend, 'parasail_version': version, 'simd': simd}


# Gapped Karlin-Altschul parameters (lambda, K) for BLOSUM62, keyed by BLAST gap costs
# (open, extend), where a gap of length k costs open + k*extend. Values as printed by
# NCBI blastp 2.17.0 for each setting.
#
# Until 2026-09-11 every SW hit was written with a placeholder E-value (0.001) and its
# RAW score in the bit-score column. SW reports the best local alignment in each of the
# six reading frames of every window, homolog or not, so each window contributed ~6 hits
# that cleared every E-value filter downstream and, carrying raw scores (about twice a
# bit score), won the spatial dedup against real tblastn/MMseqs2 hits of the same gene.
# Scored properly those hits sit at E 14-93 (median) on shuffled queries.
KARLIN_ALTSCHUL_BLOSUM62 = {
    (12, 1): (0.283, 0.059),
    (11, 1): (0.267, 0.041),   # BLAST default
    (10, 1): (0.243, 0.024),
    (9, 1): (0.206, 0.010),
}


def to_parasail_gaps(gap_open, gap_extend):
    """parasail charges `open` for a gap's first residue and `extend` for each further one
    (gap of length k = open + (k-1)*extend); BLAST charges open + k*extend. BLAST 11/1 is
    therefore parasail (12, 1)."""
    return gap_open + gap_extend, gap_extend


def karlin_altschul(score, query_len, search_space_len, gap_open, gap_extend,
                    matrix_name='BLOSUM62'):
    """(bit score, E-value) of a raw SW score. The search space is the query length times
    the total translated length searched, with no edge-effect correction (as in SSEARCH);
    for six-frame genomic windows this is conservative, because stop codons cut the target
    into short open frames (measured on shuffled queries: 2 % of per-frame best hits reach
    E <= 1 at 11/1, where a calibrated E would give ~10 %)."""
    params = KARLIN_ALTSCHUL_BLOSUM62.get((gap_open, gap_extend)) if matrix_name == 'BLOSUM62' else None
    if params is None:
        raise ValueError(
            f"No Karlin-Altschul parameters for {matrix_name} gap {gap_open}/{gap_extend}; "
            f"supported (BLAST convention): {sorted(KARLIN_ALTSCHUL_BLOSUM62)}")
    lam, k = params
    bits = (lam * score - math.log(k)) / math.log(2)
    evalue = k * query_len * search_space_len * math.exp(-lam * score)
    return bits, evalue


def min_reportable_score(query_len, search_space_len, gap_open, gap_extend,
                         matrix_name='BLOSUM62'):
    """The lowest raw score that can give a hit the runner keeps: one below the lowest score
    with E <= 1. The runner reads the E-value back from the m8 line, printed with two
    decimals ('%.2e'), and keeps E <= 1, so a hit with E up to 1.005 still counts; one raw
    score lower multiplies E by exp(lambda) (1.31 at BLAST 11/1), so every score that can be
    printed as 1.00e+00 or less is at or above this value."""
    def evalue(score):
        return karlin_altschul(score, query_len, search_space_len, gap_open, gap_extend,
                               matrix_name)[1]
    lam = KARLIN_ALTSCHUL_BLOSUM62[(gap_open, gap_extend)][0]
    s = max(0, math.ceil(math.log(max(evalue(0), 1e-300)) / lam))
    while s > 0 and evalue(s - 1) <= 1.0:
        s -= 1
    while evalue(s) > 1.0:
        s += 1
    return s - 1


def best_local_score(query_seq, target_seq, matrix_name='BLOSUM62', gap_open=11, gap_extend=1):
    """parasail's score-only Smith-Waterman (same scoring as smith_waterman_parasail): the
    result's ``score`` and ``end_ref``. On equal scores both kernels report the alignment that
    ends first in the target."""
    matrix = parasail.blosum62 if matrix_name == 'BLOSUM62' else parasail.pam100
    p_open, p_extend = to_parasail_gaps(gap_open, gap_extend)
    return parasail.sw_striped_32(query_seq, target_seq, p_open, p_extend, matrix)


def check_same_alignment(best, result, where):
    """The traceback must find the alignment the score-only pass scored; a difference would
    mean the score-only pass can skip a piece the traceback would report from."""
    if (best.score, best.end_ref) != (result['score'], result['end_ref']):
        raise RuntimeError(
            f"parasail's score-only and traceback kernels disagree on {where}: score "
            f"{best.score} / {result['score']}, end {best.end_ref} / {result['end_ref']}")


def _traceback_counts(tb):
    """(identities, mismatches, gap openings) from a parasail traceback. parasail marks
    identities '|', positive substitutions ':', others '.' and gap columns ' '."""
    comp = tb.comp or ''
    identities = comp.count('|')
    mismatches = comp.count(':') + comp.count('.')
    gapopen = 0
    for s in (tb.ref or '', tb.query or ''):
        in_gap = False
        for c in s:
            if c == '-':
                if not in_gap:
                    gapopen += 1
                    in_gap = True
            else:
                in_gap = False
    return identities, mismatches, gapopen


def translate_in_six_frames(dna_seq):
    """
    Translate DNA sequence in all 6 reading frames.
    Returns list of (frame_id, protein_seq, start_offset) tuples.
    """
    frames = []
    
    # Forward strand (frames +1, +2, +3)
    for frame in range(3):
        prot = translate(dna_seq[frame:])
        frames.append((f"+{frame+1}", prot, frame))
    
    # Reverse strand (frames -1, -2, -3)
    rev_seq = reverse_complement(dna_seq)
    for frame in range(3):
        prot = translate(rev_seq[frame:])
        frames.append((f"-{frame+1}", prot, frame))
    
    return frames


# An alignment with traceback needs memory in proportion to query length x target length:
# measured 4.2 bytes per cell (550 aa against the frames of a 2 Mb window: 1.5 GB). One
# window of 13.7 Mb once took 10 GB and the kernel ended the run. A translated frame
# longer than this many cells per query is therefore aligned in overlapping pieces, each
# about 0.8 GB at most. Up to that size (1.1 Mb of DNA for a 550-aa query, 8.6 Mb for
# 70 aa) nothing changes.
SW_MAX_CELLS = 200_000_000
# A piece cannot be shorter than four times the query (two overlaps), so for a very long
# query the cells per piece grow with the square of its length: 10,000 aa already need
# 1.7 GB, a 34,000-aa titin 20 GB. Such a query is not searched with Smith-Waterman; the
# other two search programs find a protein of that length without it.
SW_MAX_QUERY_LEN = 10_000


def six_frame_length(dna_len):
    """Residues of the six translated frames of a DNA sequence of this length: translate()
    reads (length - o) // 3 codons from offset o = 0, 1, 2, on each strand. This is the
    search space of a window (search_space_len in main())."""
    return 2 * sum(max(0, int(dna_len) - o) // 3 for o in range(3))


def smith_waterman_cells(dna_len, query_lens):
    """The cells the search of one window scores: each query of at most SW_MAX_QUERY_LEN
    residues against the window's six frames. The score-only pass scores each of them once; a
    traceback (only where a kept hit can come out) scores its piece once more. Fixed by the
    inputs, so the iterative search logs it per window and names it when a window fails."""
    frames = six_frame_length(dna_len)
    return sum(int(q) * frames for q in query_lens if int(q) <= SW_MAX_QUERY_LEN)


def protein_chunks(protein, query_len, max_cells=SW_MAX_CELLS):
    """[(start, piece)] covering ``protein``. One piece (the whole sequence) while
    query length x protein length stays within ``max_cells``; otherwise pieces of that
    many cells that overlap by twice the query length, so that every local alignment of
    the query lies entirely inside at least one piece."""
    query_len = max(1, int(query_len))
    if len(protein) * query_len <= max_cells:
        return [(0, protein)]
    overlap = 2 * query_len
    size = max(2 * overlap, max_cells // query_len)
    step = size - overlap
    return [(start, protein[start:start + size])
            for start in range(0, max(1, len(protein) - overlap), step)]


def smith_waterman_parasail(query_seq, target_seq, matrix_name='BLOSUM62',
                            gap_open=11, gap_extend=1):
    """
    Perform Smith-Waterman alignment using parasail (fast vectorized implementation).

    gap_open/gap_extend follow the BLAST convention (gap of length k costs
    open + k*extend) and are converted for parasail here.

    Returns dict with alignment details.
    """
    if not HAS_PARASAIL:
        raise ImportError("parasail library required for Smith-Waterman search")

    # Get substitution matrix
    matrix = parasail.blosum62 if matrix_name == 'BLOSUM62' else parasail.pam100
    p_open, p_extend = to_parasail_gaps(gap_open, gap_extend)

    # Run Smith-Waterman with traceback
    result = parasail.sw_trace_striped_32(
        query_seq, target_seq, p_open, p_extend, matrix
    )

    traceback = result.traceback
    aln_len = len(traceback.comp) if traceback.comp else 0
    identities, mismatches, gapopen = _traceback_counts(traceback)
    identity = (identities / aln_len) * 100 if aln_len else 0
    # Count how many reference/query positions are consumed (non-gap)
    ref_consumed = sum(1 for c in traceback.ref if c != '-') if traceback.ref else aln_len
    query_consumed = sum(1 for c in traceback.query if c != '-') if traceback.query else aln_len

    return {
        'score': result.score,
        'end_query': result.end_query,
        'end_ref': result.end_ref,
        'identity': identity,
        'length': aln_len,
        'mismatch': mismatches,
        'gapopen': gapopen,
        'ref_consumed': ref_consumed,
        'query_consumed': query_consumed,
    }


# A reported alignment is masked with this residue before the piece is aligned again.
# '*' scores -4 against every residue in BLOSUM62, so no alignment starts or ends on a
# masked residue; 'X' scores 0 against A, S and T, and an alignment could end on a masked
# residue and be refused as an overlap.
MASK_RESIDUE = '*'


def several_alignments(query_seq, frame_seq, max_alignments, evalue_of, matrix_name='BLOSUM62',
                       gap_open=11, gap_extend=1, max_cells=None, skip_below=None):
    """Up to ``max_alignments`` local alignments of the query in each piece of one translated
    frame, as [(result, piece_start)] with ``result`` as smith_waterman_parasail returns it.

    The first alignment of a piece is the one the single-alignment loop reports, unless an
    earlier piece masked part of it. After each alignment its residues are masked in FRAME
    coordinates, so a later piece that overlaps this one does not find it again, and the
    piece is aligned again: a score-only pass first, and a traceback only when the E-value of
    that score (``evalue_of``) is at most 1. A piece stops after ``max_alignments``
    alignments, at a further alignment with E > 1, after a first alignment with E > 1 (no
    later one can score higher), or when an alignment runs through a masked span (not
    reported; counted in the second return value). An alignment that ends on the last residue
    of a piece with a successor is masked in this piece only: the next piece holds it whole.

    ``skip_below`` (min_reportable_score): a piece whose best score is lower is skipped
    without a traceback when that best alignment reaches no later piece (it ends before the
    next piece starts, or on the piece's last residue, where it is masked in this piece only).
    Such a piece reports nothing the runner keeps (E > 1), and the loop above would stop after
    that alignment, whose mask could reach no other piece; so the alignments with E <= 1 are
    the same. Otherwise the piece is aligned as before.
    """
    if max_cells is None:
        max_cells = SW_MAX_CELLS
    matrix = parasail.blosum62 if matrix_name == 'BLOSUM62' else parasail.pam100
    p_open, p_extend = to_parasail_gaps(gap_open, gap_extend)
    masked = bytearray(len(frame_seq))
    found, refused = [], 0
    pieces = protein_chunks(frame_seq, len(query_seq), max_cells)
    for idx, (piece_start, piece) in enumerate(pieces):
        has_next = idx < len(pieces) - 1
        local_mask = bytearray(masked[piece_start:piece_start + len(piece)])
        residues = [MASK_RESIDUE if m else aa for aa, m in zip(piece, local_mask)]
        first = None
        if skip_below is not None:
            first = parasail.sw_striped_32(query_seq, ''.join(residues), p_open, p_extend, matrix)
            if first.score < skip_below and (
                    not has_next or first.end_ref < pieces[idx + 1][0] - piece_start
                    or first.end_ref == len(piece) - 1):
                continue
        kept, deferred = 0, False
        while kept < max_alignments:
            text = ''.join(residues)
            if kept or deferred:
                score = parasail.sw_striped_32(query_seq, text, p_open, p_extend, matrix).score
                if score <= 0 or evalue_of(score) > 1:
                    break
            result = smith_waterman_parasail(query_seq, text, matrix_name, gap_open, gap_extend)
            if first is not None and not (kept or deferred):
                check_same_alignment(first, result, f"a piece at {piece_start}")
            if result['score'] <= 0:
                break
            end = result['end_ref']
            start = end - result['ref_consumed'] + 1
            if any(local_mask[start:end + 1]):
                refused += 1
                break
            for i in range(start, end + 1):
                local_mask[i] = 1
                residues[i] = MASK_RESIDUE
            if has_next and end == len(piece) - 1 and not deferred:
                deferred = True
                continue
            masked[piece_start + start:piece_start + end + 1] = b'\x01' * (end - start + 1)
            found.append((result, piece_start))
            kept += 1
            if evalue_of(result['score']) > 1:
                break
    return found, refused


def sw_hit(q_id, t_id, result, piece_start, frame_id, offset, target_len,
           query_len, search_space_len, args):
    """One m8 record of a parasail result in a piece of a translated frame."""
    end_ref = result['end_ref'] + piece_start
    start_ref = end_ref - result['ref_consumed'] + 1
    start_query = result['end_query'] - result['query_consumed'] + 1
    if frame_id.startswith('+'):
        # Forward frame: protein pos p -> DNA pos offset + p*3
        dna_start = offset + start_ref * 3 + 1  # 1-based
        dna_end = offset + (end_ref + 1) * 3  # 1-based inclusive
    else:
        # Reverse frame: positions are in reverse-complement space
        # Convert back to forward strand coordinates.
        # BLAST m8 convention: tstart > tend signals minus strand.
        rc_start = offset + start_ref * 3
        rc_end = offset + (end_ref + 1) * 3 - 1
        dna_start = target_len - rc_start     # larger coord (1-based)
        dna_end = target_len - rc_end         # smaller coord (1-based)
    bits, evalue = karlin_altschul(result['score'], query_len, search_space_len,
                                   args.gap_open, args.gap_extend, args.matrix)
    return {
        'query': q_id,
        'target': t_id,
        'pident': result['identity'],
        'alnlen': result['length'],
        'mismatch': result['mismatch'],
        'gapopen': result['gapopen'],
        'qstart': start_query + 1,
        'qend': result['end_query'] + 1,
        'tstart': dna_start,
        'tend': dna_end,
        'evalue': evalue,
        'bits': bits,
    }


def smith_waterman_ssearch(query_faa, target_fna, output_tsv, threads=1):
    """
    Perform iterative Smith-Waterman search using Parasail (if available) or ssearch36.
    
    If Parasail is available:
        - Translates target DNA in 6 frames.
        - Runs striped SW (protein-protein).
        - Iteratively masks hits to find secondary alignments (tandem dupes).
        - Maps coordinates back to genomic DNA.
        
    If Parasail is missing:
        - Falls back to ssearch36 (external binary) with iterative masking loop.
    """
    import subprocess
    import shutil
    
    if HAS_PARASAIL:
        run_parasail_sw(query_faa, target_fna, output_tsv)
        return

    # Check if ssearch36 is available
    try:
        subprocess.run(['ssearch36', '-h'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print("ERROR: ssearch36 not found and parasail-python not installed.", file=sys.stderr)
        print("  Please install one of them:", file=sys.stderr)
        print("  pip install parasail", file=sys.stderr)
        print("  conda install -c bioconda fasta3", file=sys.stderr)
        sys.exit(1)
    
    # Read original target sequence
    targets = list(parse_fasta(target_fna))
    if not targets:
        return
        
    # We assume one target sequence for simplicity in this context
    # (Checking one region against one query set)
    t_header, t_id, t_seq_orig = targets[0]
    t_seq_mutable = list(t_seq_orig) # Mutable list of chars
    
    all_hits_lines = []
    
    # Iteration loop
    max_iter = 20
    tmp_target = f"{output_tsv}.tmp.fna"
    
    for i in range(max_iter):
        # Write current masked target
        with open(tmp_target, 'w') as f:
            f.write(f">{t_id}\n{''.join(t_seq_mutable)}\n")
            
        # Run ssearch36
        cmd = [
            'ssearch36',
            '-m', '8',  # Tabular output
            '-T', str(threads),
            '-E', '20000',  # Relaxed E-value
            '-3',  # Query is protein, target is DNA
            query_faa,
            tmp_target
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout
        except subprocess.CalledProcessError as e:
            # Ends the program (exit 1): the runner stops the genome's search instead of
            # keeping the hits of the iterations before (fix-phase item 18).
            raise RuntimeError(f"ssearch36 failed at iteration {i}: {e}; {(e.stderr or '')[-300:]}") from e
            
        if not output.strip():
            break # No more hits
            
        # Parse hits from this iteration
        new_hits_found = False
        lines = output.strip().split('\n')
        for line in lines:
            if not line.strip() or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts) < 12: continue
            
            # Hit found!
            new_hits_found = True
            all_hits_lines.append(line)
            
            # Mask region
            try:
                ts = int(parts[8])
                te = int(parts[9])
                start = min(ts, te)
                end = max(ts, te)
                
                # 1-based to 0-based
                start_0 = max(0, start - 1)
                end_0 = min(len(t_seq_mutable), end)
                
                # Mask with N
                for k in range(start_0, end_0):
                    t_seq_mutable[k] = 'N'
                    
            except ValueError:
                continue
        
        if not new_hits_found:
            break
            
    # Cleanup
    if os.path.exists(tmp_target):
        os.remove(tmp_target)
        
    # Write aggregated hits
    with open(output_tsv, 'w') as out:
        for line in all_hits_lines:
            out.write(line + '\n')
    
    print(f"Smith-Waterman search complete (ssearch36): {output_tsv} ({len(all_hits_lines)} hits)", file=sys.stderr)


def run_parasail_sw(query_faa, target_fna, output_tsv):
    """
    In-memory iterative Smith-Waterman using Parasail.
    Uses 'Best Hit per Iteration' strategy:
    1. Iterate searching all frames.
    2. Pick BEST hit globally.
    3. Mask that region in ALL frames.
    4. Repeat.
    """
    # 1. Load Sequences
    queries = list(parse_fasta(query_faa))
    targets = list(parse_fasta(target_fna))
    
    if not queries or not targets:
        return

    # Assume single target region
    t_header, t_id, t_seq_dna = targets[0]
    
    # 2. Translate Target (6 frames)
    # Store: {'seq': mutable_list_aa, 'frame': int, 'strand': str, 'offset': int}
    frames = []
    
    # Forward (+1, +2, +3)
    for i in range(3):
        seq_trans = list(translate(t_seq_dna[i:]))
        frames.append({'seq': seq_trans, 'frame': i+1, 'strand': '+', 'offset': i})
        
    # Reverse (-1, -2, -3)
    # RevComp DNA first
    rc_dna = reverse_complement(t_seq_dna)
    for i in range(3):
        seq_trans = list(translate(rc_dna[i:]))
        frames.append({'seq': seq_trans, 'frame': -(i+1), 'strand': '-', 'offset': i})
        
    all_hits = []
    # BLAST 11/1, as the main() path; real Karlin-Altschul statistics instead of the old
    # 1e-10 / raw-score placeholders. Search space = all six frames of the region.
    blast_open, blast_extend = 11, 1
    p_open, p_extend = to_parasail_gaps(blast_open, blast_extend)
    search_space_len = sum(len(f['seq']) for f in frames)

    # 3. For each query
    for q_head, q_id, q_seq in queries:
        # Use 32-bit profile; a failure ends the program (fix-phase item 18: no query is
        # left out in silence)
        profile = parasail.profile_create_32(q_seq, parasail.blosum62)
            
        # Iteration Loop (Find multiple non-overlapping hits)
        for iteration in range(20): 
            best_score = 0
            best_result = None
            best_frame_idx = -1
            
            # Search all frames
            for f_idx, frame in enumerate(frames):
                seq_str = "".join(frame['seq'])
                # sw_trace_striped_profile_32
                result = parasail.sw_trace_striped_profile_32(profile, seq_str, p_open, p_extend)
                
                if result.score > best_score:
                    best_score = result.score
                    best_result = result
                    best_frame_idx = f_idx
            
            # Check threshold
            if best_score < 40: # Raw score threshold
                break
                
            # Process Best Hit
            result = best_result
            frame = frames[best_frame_idx]
            
            # CIGAR Parsing for coordinates (a failure ends the program, fix-phase item 18)
            cigar_decoded = result.cigar.decode.decode()
                
            # Parse CIGAR to find start and alignment length
            # We must handle leading 'D's which shift the logical start
            cigar_ops = []
            curr_num = ""
            for char in cigar_decoded:
                if char.isdigit():
                    curr_num += char
                else:
                    n = int(curr_num) if curr_num else 1
                    curr_num = ""
                    cigar_ops.append((char, n))
            
            # Calculate total trace length
            t_len_aln = 0
            q_len_aln = 0
            aln_len = 0
            
            for op, n in cigar_ops:
                aln_len += n
                if op in ['M', '=', 'X', 'D']:
                    t_len_aln += n
                if op in ['M', '=', 'X', 'I']:
                    q_len_aln += n

            t_end_aa = result.end_ref
            q_end_aa = result.end_query
            
            t_start_aa = t_end_aa - t_len_aln + 1
            q_start_aa = q_end_aa - q_len_aln + 1
            
            # Adjust start for leading gaps (D or I)
            # Leading D means we skipped Ref bases -> t_start shifts forward
            # Leading I means we skipped Query bases -> q_start shifts forward
            
            # Trim leading D (Deletion from Ref)
            # Usually Parasail S-W trace starts with D if it skips beginning of ref?
            # We iterate ops to adjust starts
            for op, n in cigar_ops:
                if op == 'D':
                    t_start_aa += n # Advance match start on Ref
                elif op == 'I':
                    q_start_aa += n # Advance match start on Query
                else:
                    break # First match/sub stops the trimming

            
            # DEBUG
            # print(f"DEBUG_HIT: Iter={iteration} Frame={frame['frame']} Score={result.score}")
            # print(f"DEBUG_HIT: Cigar={cigar_decoded} T_Range=[{t_start_aa}, {t_end_aa}]")
            
            # DNA Coordinates of the Hit
            if frame['strand'] == '+':
                t_start_dna = t_start_aa * 3 + frame['offset']
                t_end_dna = t_end_aa * 3 + frame['offset'] + 2
                ts = t_start_dna + 1
                te = t_end_dna + 1
            else:
                rc_start = t_start_aa * 3 + frame['offset']
                rc_end = t_end_aa * 3 + frame['offset'] + 2
                L = len(t_seq_dna)
                te_fwd = L - rc_start
                ts_fwd = L - rc_end
                # BLAST m8 convention: tstart > tend signals minus strand
                ts = te_fwd   # larger coord first → minus strand
                te = ts_fwd   # smaller coord second
                
                # For masking, we need the internal "frame coordinates" for every frame.
                # Simplest: Define Genomic Interval [min_dna, max_dna]
                t_start_dna = min(ts, te) - 1
                t_end_dna = max(ts, te) - 1 + 1 # exclusive upper bound?
            
            # print(f"DEBUG_DNA: [{t_start_dna}, {t_end_dna}]")

            # Compute real percent identity, mismatch, and gap-open from parasail traceback
            try:
                tb = result.traceback
                if tb and tb.comp:
                    real_pident = tb.comp.count('|') / len(tb.comp) * 100
                    mismatch = tb.comp.count('.')
                    # Count gap-open events: each run of '-' in ref or query is one gap-open
                    gapopen = 0
                    for s in (tb.ref or '', tb.query or ''):
                        in_gap = False
                        for c in s:
                            if c == '-':
                                if not in_gap:
                                    gapopen += 1
                                    in_gap = True
                            else:
                                in_gap = False
                else:
                    real_pident = 0.0
                    mismatch = 0
                    gapopen = 0
            except Exception:
                real_pident = 0.0
                mismatch = 0
                gapopen = 0
            bits, evalue = karlin_altschul(result.score, len(q_seq), search_space_len,
                                           blast_open, blast_extend)
            hit_line = f"{q_id}\t{t_id}\t{real_pident:.1f}\t{aln_len}\t{mismatch}\t{gapopen}\t{q_start_aa+1}\t{q_end_aa+1}\t{ts}\t{te}\t{evalue:.2e}\t{bits:.1f}"
            all_hits.append(hit_line)
            
            # MASKING: Mask this genomic region in ALL frames
            # Genomic Interval: (ts-1) to (te) (0-based)
            g_start = min(ts, te) - 1
            g_end = max(ts, te)
            
            # print(f"DEBUG_MASKING: Interval [{g_start}, {g_end}]")
            
            for f_i, f in enumerate(frames):
                # Map genomic [g_start, g_end] to frame AA coords
                # Frame Fwd: AA_idx = (DNA_idx - offset) / 3
                # We mask conservatively: any AA that overlaps the genomic region
                
                if f['strand'] == '+':
                    start_i = int((g_start - f['offset'] - 2)/3)
                    end_i = int((g_end - f['offset'])/3) + 1 # Exclusive
                    
                else:
                    rc_g_start = len(t_seq_dna) - g_end
                    rc_g_end = len(t_seq_dna) - g_start
                    
                    start_i = int((rc_g_start - f['offset'] - 2)/3)
                    end_i = int((rc_g_end - f['offset'])/3) + 1
                
                # Clamp
                start_i = max(0, start_i)
                end_i = min(len(f['seq']), end_i)
                
                # DEBUG BEFORE
                # if f_i == 1: # Check Frame 2 specifically
                #     sample = "".join(f['seq'][30:40])
                #     print(f"DEBUG_F2_BEFORE: ...{sample}...")

                # Apply Mask
                changed = False
                for k in range(start_i, end_i):
                    f['seq'][k] = 'X'
                    changed = True
                
                # if f_i == 1:
                #      sample = "".join(f['seq'][30:40])
                #      print(f"DEBUG_F2_AFTER: ...{sample}... Changed={changed} Range=[{start_i}, {end_i}]")
                     
    # Write output
                    
    # Write output
    with open(output_tsv, 'w') as f:
        for line in all_hits:
            f.write(line + "\n")
    
    print(f"Smith-Waterman search complete (Parasail): {output_tsv} ({len(all_hits)} hits)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Smith-Waterman search for GOI in genomic regions")
    parser.add_argument("--query", required=True, help="Query protein FASTA (GOI)")
    parser.add_argument("--target", required=True, help="Target DNA FASTA (genomic region)")
    parser.add_argument("--output", required=True, help="Output TSV (BLAST m8 format)")
    parser.add_argument("--matrix", default="BLOSUM62", choices=["BLOSUM62", "PAM100"],
                       help="Substitution matrix (E-values need BLOSUM62)")
    parser.add_argument("--gap_open", type=int, default=11,
                        help="Gap opening cost, BLAST convention: a gap of length k costs "
                             "gap_open + k*gap_extend (default 11/1, BLAST's default)")
    parser.add_argument("--gap_extend", type=int, default=1, help="Gap extension cost (BLAST convention)")
    parser.add_argument("--min_score", type=int, default=50, help="Minimum alignment score")
    parser.add_argument("--min_identity", type=float, default=20.0, help="Minimum percent identity")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--method", default="auto", choices=["auto", "parasail", "ssearch36"],
                       help="Smith-Waterman implementation")
    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the target DNA, for the six-frame "
                             "translation (default 1, the standard code)")
    parser.add_argument("--max_alignments_per_piece", type=int, default=1,
                        help="Local alignments reported per piece of a translated frame "
                             "(parasail only). 1 (default): the best one. Above 1: after each "
                             "alignment its residues are masked and the piece is aligned again, "
                             "while the E-value stays at most 1 (several exons or copies in one "
                             "frame)")

    args = parser.parse_args()
    set_genetic_code(args.genetic_code)
    if args.max_alignments_per_piece < 1:
        parser.error("--max_alignments_per_piece must be at least 1")
    
    # Decide method
    method = args.method
    if method == "auto":
        method = "parasail" if HAS_PARASAIL else "ssearch36"
    
    print(f"Using Smith-Waterman implementation: {method}", file=sys.stderr)
    
    if method == "ssearch36":
        if args.max_alignments_per_piece > 1:
            print("Warning: --max_alignments_per_piece applies to the parasail path only; "
                  "ssearch36 reports its own alignments.", file=sys.stderr)
        # Use external ssearch36 program
        smith_waterman_ssearch(args.query, args.target, args.output, args.threads)
    
    elif method == "parasail":
        # Use Python parasail library (faster, in-memory)
        if not HAS_PARASAIL:
            print("ERROR: parasail not installed but requested. Install with:", file=sys.stderr)
            print("  pip install parasail", file=sys.stderr)
            sys.exit(1)
        
        # Load query proteins
        queries = list(parse_fasta(args.query))
        
        # Load target DNA and translate in 6 frames
        targets = list(parse_fasta(args.target))
        
        hits = []
        refused_total = 0
        cells = tracebacks = 0
        # Validate the statistics up front (fails loud on an unsupported matrix/gap pair).
        karlin_altschul(0, 1, 1, args.gap_open, args.gap_extend, args.matrix)
        target_frames = [(t_id, t_seq, translate_in_six_frames(t_seq))
                         for _t_header, t_id, t_seq in targets]
        # Search space = every translated residue the query was aligned against.
        search_space_len = sum(len(prot) for _t_id, _t_seq, frames in target_frames
                               for _f, prot, _o in frames)

        for q_header, q_id, q_seq in queries:
            if len(q_seq) > SW_MAX_QUERY_LEN:
                print(f"Warning: Smith-Waterman skipped for {q_id}: {len(q_seq)} aa is longer than "
                      f"{SW_MAX_QUERY_LEN} aa (the traceback would need too much memory).", file=sys.stderr)
                continue
            print(f"Searching with query: {q_id}", file=sys.stderr)
            seen = set()        # an alignment in the overlap of two pieces is found twice
            reportable = min_reportable_score(len(q_seq), search_space_len,
                                              args.gap_open, args.gap_extend, args.matrix)
            cells += len(q_seq) * search_space_len

            for t_id, t_seq, frames in target_frames:
                for frame_id, prot_seq, offset in frames:
                    if args.max_alignments_per_piece > 1:
                        # A failed alignment ends the program (exit 1), so the runner stops the
                        # genome's search instead of keeping the window without this frame.
                        found, refused = several_alignments(
                            q_seq, prot_seq, args.max_alignments_per_piece,
                            lambda score: karlin_altschul(
                                score, len(q_seq), search_space_len,
                                args.gap_open, args.gap_extend, args.matrix)[1],
                            args.matrix, args.gap_open, args.gap_extend,
                            skip_below=reportable)
                        refused_total += refused
                        for result, piece_start in found:
                            if result['score'] >= args.min_score and result['identity'] >= args.min_identity:
                                hits.append(sw_hit(q_id, t_id, result, piece_start, frame_id, offset,
                                                   len(t_seq), len(q_seq), search_space_len, args))
                        continue
                    for piece_start, piece in protein_chunks(prot_seq, len(q_seq), SW_MAX_CELLS):
                        # Score-only first; the traceback only where a kept hit can come out
                        # (score >= min_score and E low enough). A failed alignment ends the
                        # program (exit 1) instead of skipping the piece.
                        best = best_local_score(q_seq, piece, args.matrix, args.gap_open, args.gap_extend)
                        if best.score < max(args.min_score, reportable):
                            continue
                        tracebacks += 1
                        result = smith_waterman_parasail(
                            q_seq, piece,
                            args.matrix, args.gap_open, args.gap_extend
                        )
                        check_same_alignment(best, result, f"{q_id} vs {t_id} frame {frame_id}, piece at {piece_start}")

                        # Filter by thresholds
                        if result['score'] >= args.min_score and result['identity'] >= args.min_identity:
                            hit = sw_hit(q_id, t_id, result, piece_start, frame_id, offset,
                                         len(t_seq), len(q_seq), search_space_len, args)
                            if (t_id, hit['tstart'], hit['tend']) in seen:
                                continue
                            seen.add((t_id, hit['tstart'], hit['tend']))
                            hits.append(hit)
        if refused_total:
            print(f"Smith-Waterman: {refused_total} further alignment(s) ran through an already "
                  f"reported one; those pieces were searched no further.", file=sys.stderr)

        # Write output
        with open(args.output, 'w') as f:
            for hit in hits:
                f.write(f"{hit['query']}\t{hit['target']}\t{hit['pident']:.1f}\t"
                       f"{hit['alnlen']}\t{hit['mismatch']}\t{hit['gapopen']}\t"
                       f"{hit['qstart']}\t{hit['qend']}\t{hit['tstart']}\t{hit['tend']}\t"
                       f"{hit['evalue']:.2e}\t{hit['bits']:.1f}\n")
        
        work = f"{cells} cells scored" + (f", {tracebacks} tracebacks"
                                          if args.max_alignments_per_piece == 1 else "")
        print(f"Found {len(hits)} Smith-Waterman alignments above threshold ({work})",
              file=sys.stderr)


if __name__ == "__main__":
    main()
