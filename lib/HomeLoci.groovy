// The home loci BEDs of SPLIT_LOCI (batch of 2026-10-15, HG-5b). split_loci.py writes the
// primary locus first; its BED is taken by its index (the number at the end of the file name),
// not by the name locus_1.bed, which a renaming would break (a .filter on the name would then
// never emit, and the run would wait for nothing).
class HomeLoci {

    static Object primary(beds) {
        def list = (beds instanceof Collection) ? (beds as List) : [beds]
        return list.min { f -> index(f) }
    }

    static int index(f) {
        def name = f.getFileName().toString().replaceFirst(/\.[^.]*$/, '')
        def m = (name =~ /(\d+)$/)
        return m.find() ? (m.group(1) as int) : Integer.MAX_VALUE
    }
}
