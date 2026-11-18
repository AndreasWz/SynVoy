// What the end summary of a run lists (main.nf, workflow.onComplete): the files of the output
// folder. Which of them belong to this run, the caller decides (main.nf fromThisRun).
class RunSummary {

    // The region tables of this run: CLUSTER_REGIONS publishes regions/<locus>/<genome>.regions.bed,
    // one folder per home locus. Sorted by folder, then name.
    static List<File> regionBeds(File outdir, Closure<Boolean> keep) {
        def folders = (new File(outdir, 'regions').listFiles()?.findAll { d -> d.isDirectory() } ?: [])
            .sort { d -> d.name }
        return folders.collectMany { d ->
            (d.listFiles()?.findAll { f -> f.isFile() && f.name.endsWith('.regions.bed') && keep.call(f) } ?: [])
                .sort { f -> f.name }
        }
    }

    // 'tree' or 'placeholder'. compute_tree.py writes <locus>_tree.info.json beside a tree; a
    // placeholder tree has none (skip_tree, fewer than 3 GOI sequences), or one whose bootstrap
    // status is placeholder_tree (an alignment without information). An unreadable info counts as
    // a placeholder.
    // The info must be this run's: publishDir overwrites what a run writes and leaves the rest, so
    // a placeholder written into a reused output folder can sit next to an earlier run's info. The
    // caller's keep (main.nf fromThisRun) refuses an info older than this run's start on a fresh
    // run; on -resume it keeps every file, so an info that is older than its tree is refused here
    // too. COMPUTE_TREE publishes the two files together, within SAME_TASK_MS of each other.
    static final long SAME_TASK_MS = 2000

    static String treeKind(File nwk, Closure<Boolean> keep = { File f -> true }) {
        def info = new File(nwk.parentFile, nwk.name.replaceFirst(/\.nwk$/, '.info.json'))
        if (!info.isFile() || !keep.call(info) || info.lastModified() < nwk.lastModified() - SAME_TASK_MS) {
            return 'placeholder'
        }
        def status = null
        try {
            status = new groovy.json.JsonSlurper().parse(info)?.bootstrap?.status
        } catch (Exception ignored) {
            return 'placeholder'
        }
        return status == 'placeholder_tree' ? 'placeholder' : 'tree'
    }
}
