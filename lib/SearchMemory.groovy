import groovy.transform.CompileStatic
import java.lang.management.ManagementFactory
import nextflow.util.MemoryUnit

// The memory request of ITERATIVE_SEARCH (profile stage 5b; it changes no result). One task
// searches every target genome of a home locus. The genomes of a wave are searched at the same
// time, each by a worker that runs its own MMseqs2 (at most mmseqs_split_memory_limit) and loads the
// whole genome (bin/iterative_search_runner.py: process_single_genome, load_genome), so the peak
// follows the largest genomes that can be searched at once. Measured (tmp/profile_2026-10-27/checks/
// memory_measurements.md): LRZ, split 16G, one 1.04 Gb genome at a time 13.9 GiB for the whole job,
// three genomes of 1.4-2.3 Gb at a time 37-38 GiB; the laptop, split 2G, at most 2.7 GiB.
// iterative_search_memory 'auto' (the default): base = max(FLOOR, need), need = 1 GiB + for each of
// the W largest targets min(split limit, MMSEQS_PER_GB x genome) + genome + 1 GiB, W = the most
// genomes one wave searches at once (min(cpus, the largest wave)). A size the user (or a profile)
// sets is the base. Attempt n asks for n x base, after an out-of-memory kill (137) as after any
// retried exit. Under the local executor the request is at most the memory the executor has
// (executor.memory, else the JVM's total memory, a cgroup limit included): Nextflow refuses a task
// that asks for more.
class SearchMemory {

    static final MemoryUnit FLOOR = MemoryUnit.of('10 GB')
    static final long GIB = 1L << 30
    // MMseqs2 per Gb of target genome below the split limit: 11.4 GiB measured for 1.04 Gb (LRZ)
    static final double MMSEQS_PER_GB = 12.0
    static final List<String> GENOME_SUFFIXES = ['.fna', '.fa', '.fasta', '.fas']

    // Bytes of sequence in a target genome file (0 for anything else); a gzip file holds about 4
    // times its size (bin/iterative_search_runner.py, genome_scale_args).
    static long sequenceBytes(File f) {
        def name = f.name.toLowerCase()
        def gz = name.endsWith('.gz')
        def stem = gz ? name.substring(0, name.length() - 3) : name
        if (!f.isFile() || !GENOME_SUFFIXES.any { s -> stem.endsWith(s) }) {
            return 0L
        }
        return gz ? f.length() * 4 : f.length()
    }

    // A path of the task's inputs as a File. In a directive, an input path is a TaskPath whose
    // toString() is the name staged in the work folder; its target is the real path.
    static File realFile(g) {
        def p = g
        if (g != null && g.metaClass.hasProperty(g, 'target') && g.target != null) {
            p = g.target
        }
        return p instanceof File ? p : (p instanceof java.nio.file.Path ? p.toFile() : new File(p.toString()))
    }

    // The target genomes of the task's `genomes` input (a folder, a file or a list), largest first.
    static List<Long> targetBytes(genomes) {
        def files = []
        def items = genomes instanceof Collection ? genomes : [genomes]
        items.each { g ->
            def f = realFile(g)
            if (f.isDirectory()) {
                f.listFiles()?.each { c -> files << c }
            } else {
                files << f
            }
        }
        return files.collect { f -> sequenceBytes(f) }.findAll { b -> b > 0 }.sort { a, b -> b <=> a }
    }

    // The largest wave for n genomes: bin/iterative_search_runner.py assign_waves_by_rank (the
    // default), one wave with disable_wavefront, at most 5 with the old distance binning.
    static int largestWave(int n, boolean rankBinning, boolean wavefrontOff) {
        if (n <= 0) {
            return 0
        }
        if (wavefrontOff) {
            return n
        }
        if (!rankBinning) {
            return Math.min(5, n)
        }
        int i = 0
        int largest = 0
        while (i < n) {
            double q = i / (double) Math.max(1, n - 1)
            int size = q < 0.10 ? 1 : (q < 0.35 ? 2 : (q < 0.70 ? 3 : 5))
            largest = Math.max(largest, Math.min(size, n - i))
            i += size
        }
        return largest
    }

    // MMseqs2's --split-memory-limit in bytes ('8G', '512M'; binary units); 0 = no limit.
    static long splitBytes(value) {
        def m = (value ?: '0').toString().trim().toUpperCase() =~ /^([0-9.]+)\s*([BKMGT]?)B?$/
        if (!m.matches()) {
            return 0L
        }
        def power = ['', 'B', 'K', 'M', 'G', 'T'].indexOf(m.group(2))
        return (long) (new BigDecimal(m.group(1)) * (1L << (10 * Math.max(0, power - 1))))
    }

    // What one genome of `bytes` needs while it is searched: its MMseqs2 run, the genome, 1 GiB.
    static long workerBytes(long bytes, long split) {
        long mmseqs = (long) (MMSEQS_PER_GB * GIB * bytes / 1e9)
        return (split > 0 ? Math.min(split, mmseqs) : mmseqs) + bytes + GIB
    }

    static Map need(params, genomes, int cpus) {
        def sizes = targetBytes(genomes)
        int w = Math.min(Math.max(1, cpus), largestWave(sizes.size(),
                params.rank_wave_binning.toString().toBoolean(),
                params.disable_wavefront.toString().toBoolean()))
        long split = splitBytes(params.mmseqs_split_memory_limit)
        long bytes = GIB + ((sizes.take(w).collect { b -> workerBytes(b, split) }.sum() ?: 0L) as long)
        return [memory: new MemoryUnit(bytes), at_once: w,
                largest: sizes ? new MemoryUnit(sizes[0]) : new MemoryUnit(0)]
    }

    @CompileStatic
    static long machineMemory() {
        def os = (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        return os.getTotalMemorySize()
    }

    // The local executor's memory: executor.memory (or executor.$local.memory), else the machine's.
    static MemoryUnit localMemory() {
        try {
            def ex = nextflow.Global.session?.config?.executor
            if (ex instanceof Map) {
                def scoped = ex['$local']
                def v = (scoped instanceof Map && scoped.memory) ? scoped.memory : ex.memory
                if (v) {
                    return v instanceof MemoryUnit ? v : MemoryUnit.of(v.toString())
                }
            }
        } catch (Exception ignored) {
            // an unreadable setting: the machine's memory below
        }
        return new MemoryUnit(machineMemory())
    }

    // [memory: the request, line: what the console prints about it]
    static Map request(params, genomes, task) {
        int attempt = (task.attempt ?: 1) as int
        def set = params.iterative_search_memory.toString().trim()
        MemoryUnit base
        String how
        if (set == 'auto') {
            def n = need(params, genomes, ((task.cpus ?: params.iterative_search_cpus) ?: 1) as int)
            base = n.memory > FLOOR ? n.memory : FLOOR
            how = "max(${FLOOR}, ${n.memory} for ${n.at_once} genome(s) at once, the largest ${n.largest})"
        } else {
            base = MemoryUnit.of(set)
            how = "iterative_search_memory ${set}"
        }
        def mem = new MemoryUnit(base.toBytes() * attempt)
        def capped = ''
        if ((task.executor ?: 'local').toString() == 'local') {
            def local = localMemory()
            if (mem > local) {
                mem = local
                capped = ", at most ${local}, the local executor's memory"
            }
        }
        return [memory: mem,
                line  : "ITERATIVE_SEARCH memory: ${mem} (attempt ${attempt} x ${base}; ${how}${capped})"]
    }
}
