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
    static String treeKind(File nwk) {
        def info = new File(nwk.parentFile, nwk.name.replaceFirst(/\.nwk$/, '.info.json'))
        if (!info.isFile()) {
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
