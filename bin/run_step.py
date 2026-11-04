#!/usr/bin/env python3
"""Run one adjudication step and record its failure instead of losing its result in silence.

    run_step.py --step NAME [--locus L] [--genome G] --record FILE [--record FILE ...]
                [--header LINE] [--empty FILE ...] -- command [args ...]

The steps that a run may finish without (the paralog check, locus ownership, the two rescue
passes) ran under Nextflow's errorStrategy 'ignore': a failed task left no output and only the
console's count of ignored tasks said so, while a missing ownership check leaves calls
unrelabelled (fix-phase d5, item 21). The command's own output goes through unchanged. Exit 0:
exit 0. Exit 3 (sequence_utils.STEP_FAILED_EXIT: a guard) and a kill (137, 140, 143; SIGKILL and
SIGTERM read as 137 and 143): the same code, so the module's errorStrategy ends the run or retries the task as
before. Any other exit: every --record file is written as the --header lines (if any) and one
line

    # SYNVOY_STEP_FAILED<TAB>step=NAME<TAB>locus=L<TAB>genome=G<TAB>exit=N<TAB>reason=<last stderr line>

every --empty file is emptied, and this script exits 0. A crash by another signal (139, 134) is such a
failure: only 3 and the kill codes keep their code. generate_report.py collects these lines
(collect_step_failures): a warning names the genome and the step, and an ownership check that
failed marks every record it would have checked (ownership_not_checked).
"""
import argparse
import subprocess
import sys

KILL_CODES = {3, 137, 140, 143}
MARKER = "# SYNVOY_STEP_FAILED"


def marker_line(step, locus, genome, code, reason):
    reason = " ".join((reason or "").split())[:300]
    return (f"{MARKER}\tstep={step}\tlocus={locus}\tgenome={genome}\texit={code}\t"
            f"reason={reason}\n")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--" not in argv:
        sys.exit("run_step.py: the command follows '--'")
    cut = argv.index("--")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--step", required=True)
    ap.add_argument("--locus", default="")
    ap.add_argument("--genome", default="")
    ap.add_argument("--record", action="append", required=True,
                    help="an output file of the step that receives the failure line")
    ap.add_argument("--header", action="append", default=[],
                    help="a line written before the failure line (e.g. ##gff-version 3)")
    ap.add_argument("--empty", action="append", default=[],
                    help="an output file of the step that is emptied on a failure")
    args = ap.parse_args(argv[:cut])
    cmd = argv[cut + 1:]
    if not cmd:
        sys.exit("run_step.py: no command after '--'")
    proc = subprocess.run(cmd, stderr=subprocess.PIPE, text=True, errors="replace")
    sys.stderr.write(proc.stderr or "")
    code = proc.returncode
    if code == 0:
        return 0
    if code < 0:                         # killed by a signal: the shell's code for it
        code = 128 - code
    if code in KILL_CODES:
        return code
    # Any other exit, a crash by another signal (139 SIGSEGV, 134 SIGABRT) included, is recorded:
    # the module's errorStrategy retries only the kill codes (first review of the session).
    lines = [l for l in (proc.stderr or "").splitlines() if l.strip()]
    line = marker_line(args.step, args.locus, args.genome, code, lines[-1] if lines else "")
    for path in args.record:
        with open(path, "w") as fh:
            fh.writelines(h + "\n" for h in args.header)
            fh.write(line)
    for path in args.empty:
        open(path, "w").close()
    sys.stderr.write(f"[run_step] {args.step} failed (exit {code})"
                     + (f" for {args.genome}" if args.genome else "")
                     + (f", {args.locus}" if args.locus else "")
                     + "; recorded for the report.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
