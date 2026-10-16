// The gene gap of the home loci (batch of 2026-10-15, HG-1C): two places where the query aligns
// at home are one locus when at most this many annotated genes lie between them, since their
// flanking genes would largely be the same. 'auto' (the default) is half the flanking genes per
// side, floor(n_flanking_genes / 2): 5 at the default of 10, as the fixed default was. Resolved
// here, once for split_loci.py and profile_hits.py, which take a number (a string fails in
// both, and profile_hits.py's failure is swallowed by LOCATE_GENE, so hit_profile.json and the
// preset choice would vanish without a word). It counts only with locus_gap = 'auto' and a home
// annotation: a number for locus_gap (and -profile legacy_distances) measures the gap in bp.
class LocusGap {

    static boolean isAuto(params) {
        return params.locus_gap_genes.toString().trim() == 'auto'
    }

    static int genes(params) {
        if (isAuto(params)) {
            return new BigDecimal(params.n_flanking_genes.toString().trim()).intValue().intdiv(2)
        }
        return new BigDecimal(params.locus_gap_genes.toString().trim()).intValue()
    }
}
