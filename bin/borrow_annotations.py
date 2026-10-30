#!/usr/bin/env python3
"""
borrow_annotations.py — gene models for a home genome without an annotation, transferred from
an annotated target genome (batch of 2026-10-15, HG-5b; spec_5).

Without a home GFF the home genome has only gene predictions around the query's hits
(PREPARE_HOME_PROTEOME: Prodigal or AUGUSTUS). This script adds the genes of the closest
annotated target as real gene models, so they can be anchors:

  1. The home locus: the primary locus BED of split_loci (its best_interval row, else the row
     with the best bit score; never just the first row), and a window around it of
     --pred_flank_window bp, widened until it holds 2 x n_flanking predicted genes or 500 kb.
  2. The target: the closest annotated target (>= 10 CDS rows) by the distances of
     --sorted_genomes; between several at that distance, or without distances, the one whose
     hits of the window's predicted proteins have the highest median identity.
  3. The target region: the window's predicted proteins searched against the target (MMseqs2,
     E <= 1e-3); their hits clustered per scaffold with the 50 kb gap of split_loci; the
     cluster with the most distinct home proteins, +- --target_flank bp.
  4. Each annotated protein of that region aligned to the home window with miniprot (the home
     genome's genetic code); a model is kept at identity >= 30 % and query coverage >= 50 %
     and only when the target region's best protein for it is the one it came from
     (reciprocal check, Smith-Waterman).
  5. One model per place: a borrowed model replaces a prediction it overlaps by more than half
     of the shorter one (it has a name and introns); of two borrowed models at one place the
     one with the higher identity stays.
  6. One merged GFF (predictions + borrowed models, gene/mRNA/CDS rows with phases,
     ID=borrowed_<target gene>, Note=transferred from <target>) and the home proteome built
     from it (the longest coding model per gene, as gff_to_faa.py does), without the genes
     over the query's own home model (goi_annotation.bed), with home_genes.tsv and
     home_annotation.json (source, counts, partial: true).

Before 2026-10-16 the borrowed rows were gene rows without CDS (the walk read them as
non-coding: drawn, never anchors), placed by the first row of the hit BED (the alphabetically
first scaffold with a hit), from the best-annotated target rather than the closest one, with a
target region of min..max of all hits.
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    from sequence_utils import (
        load_genome, parse_gff, get_feature_id, parse_bed, write_fasta, translate,
        reverse_complement, cds_models_by_gene, cds_model_protein, sw_bits,
        genetic_code_flag, genetic_code_for, parse_genetic_code_spec, set_genetic_code,
    )
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sequence_utils import (
        load_genome, parse_gff, get_feature_id, parse_bed, write_fasta, translate,
        reverse_complement, cds_models_by_gene, cds_model_protein, sw_bits,
        genetic_code_flag, genetic_code_for, parse_genetic_code_spec, set_genetic_code,
    )

MIN_TARGET_CDS = 10            # an annotated target needs at least this many CDS rows
WINDOW_MAX_BP = 500_000        # the home window grows up to this
CLUSTER_GAP_BP = 50_000        # target hits closer than this are one cluster (split_loci's gap)
MIN_IDENTITY = 0.30            # a borrowed model's miniprot identity ...
MIN_QUERY_COVERAGE = 0.50      # ... and its coverage of the target protein
BEST_INTERVAL = 'best_interval'


# ------------------------------------------------------------------ the home locus
def goi_position(bed_path):
    """``(chrom, start, end)`` of the query's home locus: the best_interval row of the primary
    locus BED (split_loci, batch HG-1B), else its row with the best bit score (column 7), else
    the first row. None for an empty BED."""
    rows = parse_bed(bed_path)
    if not rows:
        return None
    marked = [r for r in rows if r.get('name') == BEST_INTERVAL]
    if marked:
        r = marked[0]
        return r['chrom'], r['start'], r['end']

    def bits(r):
        try:
            return float(r.get('display_name') or 'nan')
        except ValueError:
            return float('nan')
    scored = [(bits(r), i, r) for i, r in enumerate(rows)]
    scored = [t for t in scored if t[0] == t[0]]          # drop NaN
    r = max(scored, key=lambda t: (t[0], -t[1]))[2] if scored else rows[0]
    return r['chrom'], r['start'], r['end']


def predicted_models(gff):
    """{gene_id: model} of a GFF (cds_models_by_gene, the longest coding model per gene)."""
    if not gff or not os.path.exists(gff) or os.path.getsize(gff) == 0:
        return {}
    return cds_models_by_gene(gff)


def model_span(info):
    return (min(s for s, _e, _p in info['exons']) - 1, max(e for _s, e, _p in info['exons']))


def home_window(goi, models, chrom_len, flank, n_flanking):
    """The window around the home locus: +- ``flank`` bp, widened (doubled) until it holds
    2 x n_flanking predicted genes or reaches WINDOW_MAX_BP; returns (start, end)."""
    chrom, start, end = goi
    on_chrom = [model_span(m) for m in models.values() if m['chrom'] == chrom]
    pad = max(1, int(flank))
    while True:
        lo, hi = max(0, start - pad), min(chrom_len, end + pad)
        n = sum(1 for s, e in on_chrom if s >= lo and e <= hi)
        if n >= 2 * n_flanking or (hi - lo) >= WINDOW_MAX_BP or (lo == 0 and hi == chrom_len):
            return lo, hi
        pad *= 2


# ------------------------------------------------------------------ the target
def annotated_targets(genomes_dir):
    """[(genome fna, its gff, CDS rows)] of the targets with an annotation of >= 10 CDS rows."""
    out = []
    for fna in sorted(Path(genomes_dir).glob("*.fna")):
        gff = fna.with_suffix('.gff')
        if not gff.exists():
            continue
        n = 0
        with open(gff) as fh:
            for line in fh:
                if not line.startswith('#') and line.count('\t') >= 2 and line.split('\t', 3)[2] == 'CDS':
                    n += 1
        if n >= MIN_TARGET_CDS:
            out.append((str(fna), str(gff), n))
    return out


def read_sorted_genomes(path):
    """``{genome file name: distance}`` from PHYLO_SORT's sorted_genomes.txt (name, tab,
    distance; taxonomy steps or a MinHash distance). {} without the file, and when every
    distance is the same (the alphabetical fallback writes 0 for all): no order then."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if not f[0].strip():
                continue
            try:
                out[os.path.basename(f[0].strip())] = float(f[1])
            except (IndexError, ValueError):
                continue
    return out if len(set(out.values())) > 1 else {}


def choose_target(targets, distance, median_identity=None):
    """The annotated target to borrow from: the one at the smallest phylogenetic
    ``distance``; several at that distance, or no distance known, are told apart by
    ``median_identity(fna)`` (the median identity of the home proteins' hits; ties: the
    first by name). ``(fna, gff, why)`` or None."""
    if not targets:
        return None
    known = [t for t in targets if os.path.basename(t[0]) in distance]
    pool, why = targets, "no phylogenetic distance known"
    if known:
        d = min(distance[os.path.basename(t[0])] for t in known)
        pool = [t for t in known if distance[os.path.basename(t[0])] == d]
        why = f"the closest annotated target (distance {d:g})"
        if len(pool) == 1:
            return pool[0][0], pool[0][1], why
        why += f", {len(pool)} at that distance"
    if median_identity is None:
        return pool[0][0], pool[0][1], why + "; the first by name"
    scored = [(median_identity(t[0]), -i, t) for i, t in enumerate(pool)]
    ident, _i, best = max(scored, key=lambda x: (x[0], x[1]))
    return best[0], best[1], (f"{why}; the highest median identity of the home proteins' hits "
                              f"({ident:.1f} %)")


def mmseqs_hits(query_faa, genome_fna, evalue=1e-3):
    """Hits of proteins against a genome (MMseqs2 translated search): dicts with query,
    chrom, start (0-based), end, pident, evalue, bits; [] when the search fails."""
    with tempfile.TemporaryDirectory(prefix="synvoy_borrow_") as tmp:
        out = os.path.join(tmp, "hits.m8")
        cmd = ['mmseqs', 'easy-search', query_faa, genome_fna, out, os.path.join(tmp, 'tmp'),
               '--search-type', '2', '-s', '7.5', '-e', str(evalue), '-v', '1',
               '--format-output',
               'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'
               ] + genetic_code_flag('mmseqs')
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print(f"  WARNING: MMseqs2 search against {Path(genome_fna).name} failed: {exc}")
            return []
        hits = []
        if os.path.exists(out):
            with open(out) as fh:
                for line in fh:
                    p = line.rstrip('\n').split('\t')
                    if len(p) < 12:
                        continue
                    a, b = int(p[8]), int(p[9])
                    hits.append({'query': p[0], 'chrom': p[1], 'start': min(a, b) - 1,
                                 'end': max(a, b), 'pident': float(p[2]),
                                 'evalue': float(p[10]), 'bits': float(p[11])})
        return hits


def target_region(hits, flank):
    """The target region: hits clustered per scaffold with CLUSTER_GAP_BP, the cluster with
    the most distinct home proteins (ties: more hits, then position), +- ``flank`` bp:
    ``(chrom, start, end, n proteins)`` or None."""
    best = None
    by_chrom = {}
    for h in hits:
        by_chrom.setdefault(h['chrom'], []).append(h)
    for chrom in sorted(by_chrom):
        cl = []
        for h in sorted(by_chrom[chrom], key=lambda x: (x['start'], x['end'])) + [None]:
            if h is not None and cl and h['start'] - max(x['end'] for x in cl) <= CLUSTER_GAP_BP:
                cl.append(h)
                continue
            if cl:
                key = (len({x['query'] for x in cl}), len(cl), -min(x['start'] for x in cl))
                if best is None or key > best[0]:
                    best = (key, chrom, min(x['start'] for x in cl), max(x['end'] for x in cl))
            cl = [h] if h is not None else []
    if best is None:
        return None
    return best[1], max(0, best[2] - flank), best[3] + flank, best[0][0]


def region_proteins(gff, genome_fna, chrom, start, end):
    """[(gene id, name, protein)] of the annotated genes of a target region (the longest coding
    model per gene)."""
    genome = load_genome(genome_fna)
    if chrom not in genome:
        return []
    names = {}
    for f in parse_gff(gff):
        if f.get('type') == 'gene' and f.get('seqid') == chrom:
            a = f.get('attributes', {})
            names[get_feature_id(f)] = a.get('Name') or a.get('gene') or get_feature_id(f)
    out = []
    for gid, info in sorted(cds_models_by_gene(gff, {chrom}).items()):
        s, e = model_span(info)
        if e <= start or s >= end:
            continue
        prot = cds_model_protein(info, genome[chrom])
        if len(prot) >= 10:
            out.append((gid, names.get(gid, gid), prot))
    return out


# ------------------------------------------------------------------ back to home
def borrowed_models(proteins, window_seq, chrom, offset, max_intron):
    """Each target protein aligned to the home window (miniprot): the kept models as dicts
    with source id, name, exons (global, 0-based half-open, phase), strand, identity,
    coverage and protein."""
    import annotate_goi_exons as ag
    out = []
    for gid, name, prot in proteins:
        models = ag._run_miniprot(prot, window_seq, chrom, max_intron=max_intron)
        if not models:
            continue
        best = models[0]
        exons, protein = ag._build_model(best, prot, window_seq, chrom, region_offset=offset)
        if not exons:
            continue
        covered = sum(max(0, int(e['qend']) - int(e['qstart']) + 1) for e in exons)
        coverage = min(1.0, covered / max(1, len(prot)))
        if best['identity'] < MIN_IDENTITY or coverage < MIN_QUERY_COVERAGE:
            continue
        out.append({'source': gid, 'name': name, 'strand': best['strand'],
                    'exons': [(int(e['gstart']), int(e['gend']), int(e.get('phase', 0) or 0)) for e in exons],
                    'identity': best['identity'], 'coverage': coverage,
                    'protein': protein.rstrip('*')})
    return out


def reciprocal(models, proteins):
    """The models whose best protein of the target region (Smith-Waterman bits) is the one
    they came from."""
    kept = []
    for m in models:
        if not m['protein']:
            continue
        best = max(proteins, key=lambda p: (sw_bits(m['protein'], p[2]), p[0] == m['source']))
        if best[0] == m['source']:
            kept.append(m)
    return kept


def _span(exons):
    return min(s for s, _e, _p in exons), max(e for _s, e, _p in exons)


def _overlaps_half(a, b):
    ov = min(a[1], b[1]) - max(a[0], b[0])
    return ov > 0 and ov > 0.5 * min(a[1] - a[0], b[1] - b[0])


def one_model_per_place(borrowed, predicted):
    """``(borrowed kept, predicted gene ids dropped)``: of two borrowed models at one place
    (overlap > half of the shorter) the higher identity stays; a prediction overlapping a
    kept borrowed model that much is dropped."""
    kept = []
    for m in sorted(borrowed, key=lambda m: (-m['identity'], m['source'])):
        span = _span(m['exons'])
        if not any(_overlaps_half(span, _span(k['exons'])) for k in kept):
            kept.append(m)
    dropped = set()
    for gid, info in predicted.items():
        ps = model_span(info)
        for m in kept:
            if info['chrom'] == m['chrom'] and _overlaps_half(ps, _span(m['exons'])):
                dropped.add(gid)
                break
    return kept, dropped


# ------------------------------------------------------------------ outputs
def borrowed_gff_rows(m, target_name):
    s, e = _span(m['exons'])
    gid = f"borrowed_{m['source']}"
    note = f"transferred from {target_name}"
    attrs = f"ID={gid};Name={m['name']};gene_biotype=protein_coding;Note={note}"
    rows = [f"{m['chrom']}\tborrowed_annotation\tgene\t{s + 1}\t{e}\t.\t{m['strand']}\t.\t{attrs}",
            f"{m['chrom']}\tborrowed_annotation\tmRNA\t{s + 1}\t{e}\t{m['identity'] * 100:.1f}\t"
            f"{m['strand']}\t.\tID={gid}.t1;Parent={gid};Identity={m['identity']:.3f};"
            f"QueryCoverage={m['coverage']:.3f}"]
    for i, (cs, ce, ph) in enumerate(sorted(m['exons']), 1):
        rows.append(f"{m['chrom']}\tborrowed_annotation\tCDS\t{cs + 1}\t{ce}\t.\t{m['strand']}\t{ph}\t"
                    f"ID={gid}.cds{i};Parent={gid}.t1")
    return rows


def write_merged_gff(path, predicted_gff, dropped, borrowed, target_name):
    """The predictions (without the dropped genes and their children) and the borrowed models."""
    drop_ids = set(dropped)
    lines = ["##gff-version 3"]
    if predicted_gff and os.path.exists(predicted_gff):
        rows = [ln.rstrip('\n') for ln in open(predicted_gff) if ln.strip() and not ln.startswith('#')]
        parent = {}
        for ln in rows:
            a = dict(kv.split('=', 1) for kv in ln.split('\t')[8].split(';') if '=' in kv) if ln.count('\t') >= 8 else {}
            if a.get('ID') and a.get('Parent'):
                parent[a['ID']] = a['Parent']

        def top(fid):
            seen = set()
            while fid in parent and fid not in seen:
                seen.add(fid)
                fid = parent[fid]
            return fid
        for ln in rows:
            f = ln.split('\t')
            if len(f) < 9:
                continue
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            fid = a.get('ID') or a.get('Parent', '')
            if top(fid) in drop_ids or top(a.get('Parent', '')) in drop_ids:
                continue
            lines.append(ln)
    for m in borrowed:
        lines.extend(borrowed_gff_rows(m, target_name))
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def write_home_proteome(merged_gff, genome, goi_rows, out_dir, record):
    """The home proteome of a run without a home annotation: the longest coding model of every
    gene of the merged GFF, translated (as gff_to_faa.py), without the genes over the query's
    own home model; home_genes.tsv and home_annotation.json beside it."""
    os.makedirs(out_dir, exist_ok=True)
    proteins, coords, left_out = [], [], []
    for gid, info in sorted(predicted_models(merged_gff).items()):
        if info['chrom'] not in genome:
            continue
        s, e = model_span(info)
        if any(c == info['chrom'] and s < ge and e > gs for c, gs, ge in goi_rows):
            left_out.append(gid)
            continue
        prot = cds_model_protein(info, genome[info['chrom']])
        if prot:
            proteins.append((gid, prot))
            coords.append((gid, info['chrom'], s + 1, e, info['strand']))
    write_fasta(proteins, os.path.join(out_dir, 'home_proteome.faa'))
    with open(os.path.join(out_dir, 'home_genes.tsv'), 'w') as fh:
        fh.write("gene_id\tchrom\tstart\tend\tstrand\n")
        for row in coords:
            fh.write('\t'.join(map(str, row)) + '\n')
    record.update(proteins=len(proteins), left_out_over_goi=left_out)
    with open(os.path.join(out_dir, 'home_annotation.json'), 'w') as fh:
        json.dump(record, fh, indent=1)
    return len(proteins)


def main():
    parser = argparse.ArgumentParser(
        description="Gene models from an annotated target genome for a home genome without annotation")
    parser.add_argument("--home_genome", required=True, help="Home genome FASTA")
    parser.add_argument("--predicted_gff", required=True,
                        help="The home gene predictions (PREPARE_HOME_PROTEOME's home_predicted.gff)")
    parser.add_argument("--genomes_dir", required=True,
                        help="Directory with the target genomes and, for some, <genome>.gff")
    parser.add_argument("--goi_bed", required=True, help="The primary home locus BED (split_loci)")
    parser.add_argument("--sorted_genomes", default=None,
                        help="PHYLO_SORT's sorted_genomes.txt: targets closest first")
    parser.add_argument("--goi_annotation", default=None,
                        help="goi_annotation.bed of ANNOTATE_GOI: genes over it are left out of "
                             "the home proteome")
    parser.add_argument("--n_flanking", type=int, default=10, help="Flanking genes per side")
    parser.add_argument("--window", type=int, default=50000,
                        help="Home window around the locus, widened as needed (pred_flank_window)")
    parser.add_argument("--target_flank", type=int, default=100000,
                        help="bp added on each side of the target region (borrow_target_flank)")
    parser.add_argument("--output_gff", required=True, help="The merged home annotation")
    parser.add_argument("--output_proteins", required=True,
                        help="The target proteins of the region (FASTA)")
    parser.add_argument("--output_db_dir", required=True,
                        help="Folder for home_proteome.faa, home_genes.tsv, home_annotation.json")
    parser.add_argument("--genetic_code", type=int, default=1,
                        help="NCBI genetic code of the home genome (default 1)")
    parser.add_argument("--target_genetic_codes", default="",
                        help="Genetic codes of target genomes, 'genome=table,...'")
    parser.add_argument("--default_target_genetic_code", type=int, default=1,
                        help="NCBI genetic code of a target genome not named there (default 1)")
    parser.add_argument("--miniprot_timeout", type=int, default=3600,
                        help="Guard of each miniprot run that transfers a target protein to the "
                             "home window, in seconds (pipeline: flanking_miniprot_timeout): a run "
                             "that does not finish stops the task (until 2026-10-28 a timeout of "
                             "120 s left that protein out)")
    args = parser.parse_args()
    set_genetic_code(args.genetic_code)
    if args.miniprot_timeout < 1:
        parser.error("--miniprot_timeout must be at least 1")
    import annotate_goi_exons as ag
    ag.FLANKING_MINIPROT_TIMEOUT = args.miniprot_timeout

    genome = load_genome(args.home_genome)
    predicted = predicted_models(args.predicted_gff)
    goi_rows = [(r['chrom'], r['start'], r['end']) for r in parse_bed(args.goi_annotation)] \
        if args.goi_annotation else []
    record = {'source': 'predicted', 'partial': True, 'target': None, 'borrowed': 0,
              'predicted': len(predicted), 'predictions_dropped': 0}

    def finish(borrowed, dropped, target_name, proteins):
        write_merged_gff(args.output_gff, args.predicted_gff, dropped, borrowed, target_name)
        write_fasta([(f"{gid} {name}", p) for gid, name, p in proteins], args.output_proteins)
        record.update(borrowed=len(borrowed), predictions_dropped=len(dropped),
                      predicted=len(predicted) - len(dropped))
        if borrowed:
            record['source'] = 'predicted+borrowed'
        n = write_home_proteome(args.output_gff, genome, goi_rows, args.output_db_dir, record)
        print(f"[borrow_annotations] {len(borrowed)} borrowed gene model(s), "
              f"{len(predicted) - len(dropped)} prediction(s) kept, {len(dropped)} replaced; "
              f"home proteome {n} protein(s) (partial: genes near the query's hits only).")

    goi = goi_position(args.goi_bed)
    if goi is None or goi[0] not in genome:
        print("[borrow_annotations] No usable home locus: predictions only.")
        return finish([], set(), '', [])
    lo, hi = home_window(goi, predicted, len(genome[goi[0]]), args.window, args.n_flanking)
    in_window = {gid: m for gid, m in predicted.items()
                 if m['chrom'] == goi[0] and model_span(m)[0] >= lo and model_span(m)[1] <= hi}
    print(f"[borrow_annotations] Home locus {goi[0]}:{goi[1]}-{goi[2]}; window {lo}-{hi} "
          f"({len(in_window)} predicted gene(s)).")
    window_proteins = [(gid, cds_model_protein(m, genome[goi[0]])) for gid, m in sorted(in_window.items())]
    window_proteins = [(g, p) for g, p in window_proteins if len(p) >= 30]
    targets = annotated_targets(args.genomes_dir)
    if not targets or not window_proteins:
        print("[borrow_annotations] No annotated target genome, or no predicted protein near "
              "the locus: predictions only.")
        return finish([], set(), '', [])

    with tempfile.TemporaryDirectory(prefix="synvoy_borrow_") as tmp:
        qfa = os.path.join(tmp, "window.faa")
        write_fasta(window_proteins, qfa)
        codes = parse_genetic_code_spec(args.target_genetic_codes)
        hits_by = {}

        def hits_of(fna):
            if fna not in hits_by:
                set_genetic_code(genetic_code_for(fna, codes, args.default_target_genetic_code))
                hits_by[fna] = mmseqs_hits(qfa, fna)
            return hits_by[fna]

        def median_identity(fna):
            hits = hits_of(fna)
            return statistics.median(h['pident'] for h in hits) if hits else 0.0

        choice = choose_target(targets, read_sorted_genomes(args.sorted_genomes), median_identity)
        target_fna, target_gff, why = choice
        target_name = Path(target_fna).stem
        print(f"[borrow_annotations] Target {Path(target_fna).name}: {why}.")
        region = target_region(hits_of(target_fna), args.target_flank)
    set_genetic_code(args.genetic_code)
    if region is None:
        print("[borrow_annotations] The home proteins do not align to the target: predictions only.")
        return finish([], set(), target_name, [])
    chrom, rs, re_, n_prot = region
    proteins = region_proteins(target_gff, target_fna, chrom, rs, re_)
    print(f"[borrow_annotations] Target region {chrom}:{rs}-{re_} ({n_prot} home protein(s) "
          f"align there; {len(proteins)} annotated protein(s)).")
    window_seq = genome[goi[0]][lo:hi]
    # miniprot's intron bound: the window itself, at least 20 kb (as for the home gene model
    # built from hits, annotate_goi_exons.py).
    models = borrowed_models(proteins, window_seq, goi[0], lo, max(20000, hi - lo))
    for m in models:
        m['chrom'] = goi[0]
    recip = reciprocal(models, proteins)
    kept, dropped = one_model_per_place(recip, predicted)
    print(f"[borrow_annotations] {len(models)} model(s) at >= {MIN_IDENTITY:.0%} identity and "
          f">= {MIN_QUERY_COVERAGE:.0%} coverage, {len(recip)} reciprocal, {len(kept)} kept "
          f"(one per place).")
    record['target'] = target_name
    record['target_choice'] = why
    record['target_region'] = f"{chrom}:{rs}-{re_}"
    return finish(kept, dropped, target_name, proteins)


if __name__ == "__main__":
    main()
