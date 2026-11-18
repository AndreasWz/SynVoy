// The task logs of a run (main.nf, workflow.onComplete): stdout, stderr and the script of each task
// this session ran, copied from Nextflow's work folder into <outdir>/logs/.
class TaskLogs {

    // The task folders of this session: <workDir>/<2 hex>/<30 hex>/ whose .command.run is not older
    // than sinceMs (this session wrote it; a task taken from the cache on -resume keeps its old one).
    // Two levels are listed, so no other folder of the work folder is entered (the conda and
    // Singularity caches, a task's own subfolders). Sorted by path.
    static List<File> ofThisRun(File workDir, long sinceMs) {
        def out = []
        def tops = (workDir.listFiles()?.findAll { d -> d.isDirectory() && d.name ==~ /[0-9a-f]{2}/ } ?: [])
            .sort { d -> d.name }
        tops.each { top ->
            def tasks = (top.listFiles()?.findAll { d -> d.isDirectory() && d.name ==~ /[0-9a-f]{30}/ } ?: [])
                .sort { d -> d.name }
            tasks.each { task ->
                def run = new File(task, '.command.run')
                if (run.isFile() && run.lastModified() >= sinceMs) {
                    out << task
                }
            }
        }
        return out
    }

    // <PROCESS>__<tag>__<hash>: the process (with its tag) from the header of .command.run, its
    // first output, and the first 8 characters of the task hash (the two folder names), so two
    // tasks never share a folder and a rerun of the same tasks gives the same names.
    static String folderName(File task) {
        def processName = 'unknown'
        def tag = ''
        try {
            task.toPath().resolve('.command.run').withReader { reader ->
                for (int i = 0; i < 12; i++) {
                    def line = reader.readLine()
                    if (line == null) {
                        break
                    }
                    def m = (line =~ /name:\s*'([^']+)'/)
                    if (m.find()) {
                        processName = m.group(1).replaceAll(/[^a-zA-Z0-9_.-]/, '_')
                    }
                    def out = (line =~ /- '([^']+)'/)
                    if (out.find() && !tag) {
                        tag = out.group(1).replaceAll(/.*\//, '').replaceAll(/[^a-zA-Z0-9_.-]/, '_')
                    }
                }
            }
        } catch (Exception ignored) {
            // an unreadable header: the folder is named 'unknown' with the task hash
        }
        def name = (tag ? "${processName}__${tag}" : processName).take(120)
        return "${name}__${task.parentFile.name}${task.name.take(6)}"
    }

    // Copies the non-empty .command.log, .command.err and .command.sh of each task of this session
    // into logsDir/<folderName>/ as stdout.log, stderr.log and script.sh; a task whose log and error
    // output are both empty gets no folder. Returns the number of folders written.
    static int collect(File workDir, File logsDir, long sinceMs) {
        if (!workDir.isDirectory()) {
            return 0
        }
        int n = 0
        ofThisRun(workDir, sinceMs).each { task ->
            def log = new File(task, '.command.log')
            def err = new File(task, '.command.err')
            if (!(log.isFile() && log.length() > 0) && !(err.isFile() && err.length() > 0)) {
                return
            }
            def dest = new File(logsDir, folderName(task))
            dest.mkdirs()
            [[log, 'stdout.log'], [err, 'stderr.log'], [new File(task, '.command.sh'), 'script.sh']].each { pair ->
                File src = pair[0] as File
                if (src.isFile() && src.length() > 0) {
                    new File(dest, pair[1] as String).bytes = src.bytes
                }
            }
            n++
        }
        return n
    }
}
