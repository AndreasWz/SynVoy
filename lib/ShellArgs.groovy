// A value quoted for the bash of a process's script: single quotes keep every byte literal; a ' is
// closed, given in double quotes and opened again (the idiom of HomeGenes.resolverArgs). Groovy fills
// in ${...} before bash reads the script, so the quoting has to happen here (fix item 28: ANNOTATE_GOI
// passed --query_id unquoted and --gene_name in plain single quotes, so 14-3-3(epsilon) or
// 5'-nucleotidase ended the task with a shell syntax error, and P04637;x ran a second command;
// RESOLVE_GENE_INPUT gave the same easy-mode input in double quotes, where $(...) and ` run).
// A value without ' renders as 'value': as before where it was in single quotes (--gene_name, the
// home genes, --run_resource), in single quotes now where it was not (--query_id, RESOLVE_GENE_INPUT's
// input and --species), so those tasks get a new script.
class ShellArgs {

    static String q(Object v) {
        return "'" + v.toString().replace("'", "'\"'\"'") + "'"
    }
}
