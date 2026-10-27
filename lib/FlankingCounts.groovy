// Counts tied to n_flanking_genes (profile handoff, stages 6 and 4: constraint C9). Each is 'auto'
// by default and then follows n = n_flanking_genes; at the default of 10 it is the fixed value it
// had before, so such a run passes the same numbers to the same scripts:
//   synteny_bridge_max_gap_genes        2n               20  one home neighbourhood of n genes per side
//   strong_synteny_min_flanking         n / 2 rounded up  5  a count of at least n / 2 (0 would switch
//                                                            the region class off)
//   goi_hull_min_flanking               round(0.4 n)      4  (0 acts as 1: a hull holds at least one)
//   strong_synteny_no_goi_min_flanking  2n               20  the report's strong_synteny_no_goi flag
// A number is used as given. Resolved here, as LocusGap resolves locus_gap_genes (floor(n / 2)),
// because the runner, cluster_grs.py, rescue_goi_hull.py and generate_report.py take numbers.
class FlankingCounts {

    static final List<String> NAMES = ['synteny_bridge_max_gap_genes', 'strong_synteny_min_flanking',
                                       'goi_hull_min_flanking', 'strong_synteny_no_goi_min_flanking']

    static boolean isAuto(params, String name) {
        return params[name].toString().trim() == 'auto'
    }

    static int auto(String name, int n) {
        switch (name) {
            case 'synteny_bridge_max_gap_genes':
                return 2 * n
            case 'strong_synteny_min_flanking':
                return (n + 1).intdiv(2)
            case 'goi_hull_min_flanking':
                return (0.4 * n).setScale(0, java.math.RoundingMode.HALF_UP).intValue()
            case 'strong_synteny_no_goi_min_flanking':
                return 2 * n
        }
        throw new IllegalArgumentException("not a count tied to n_flanking_genes: ${name}")
    }

    static Object value(params, String name) {
        if (!isAuto(params, name)) {
            return params[name]
        }
        return auto(name, new BigDecimal(params.n_flanking_genes.toString().trim()).intValue())
    }
}
