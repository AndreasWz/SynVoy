#!/usr/bin/env python3

import argparse
import json
import os
import re
import shutil
import subprocess
import sys


def read_fasta_entries(fasta_path):
    """Return FASTA entries as a list of (name, sequence)."""
    entries = []
    current_name = None
    current_seq = []

    with open(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    entries.append((current_name, "".join(current_seq)))
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        entries.append((current_name, "".join(current_seq)))

    return entries


def sanitize_newick_label(name):
    """Convert FASTA headers into safe unquoted Newick labels."""
    safe = []
    for ch in name:
        if ch.isalnum() or ch in "._-|":
            safe.append(ch)
        else:
            safe.append("_")
    label = "".join(safe).strip("_")
    return label or "placeholder"


def write_placeholder_tree(input_fasta, output_nwk, reason):
    """
    Write a deterministic placeholder tree preserving sequence labels.
    """
    names = [sanitize_newick_label(name) for name, _ in read_fasta_entries(input_fasta)]
    if not names:
        tree = "(placeholder:0.0);"
    elif len(names) == 1:
        tree = f"({names[0]}:0.0);"
    else:
        leaves = ",".join(f"{name}:0.0" for name in names)
        tree = f"({leaves});"

    with open(output_nwk, "w") as handle:
        handle.write(tree + "\n")

    print(f"Wrote placeholder tree ({reason}) to {output_nwk}")


def is_data_limited_iqtree_error(text):
    """
    Detect known IQ-TREE failures caused by tiny/uninformative alignments.
    """
    lowered = (text or "").lower()
    markers = [
        "less than 4 sequences",
        "no parsimony-informative",
        "all sites are constant",
        "not enough variation",
        "too few taxa",
    ]
    return any(marker in lowered for marker in markers)


def run_mafft(input_fasta, output_aln, threads=1):
    """Run MAFFT alignment."""
    cmd = ["mafft", "--amino", "--auto", "--thread", str(threads), input_fasta]
    with open(output_aln, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)


# IQ-TREE's random seed: fixed, so the same alignment gives the same tree run after run (fix-phase
# d6; without --seed IQ-TREE seeds from the clock). The seed and the threads are recorded beside
# the tree (write_info), as the thread count can change IQ-TREE's search.
IQTREE_SEED = 1


def write_info(path, **info):
    """<locus>_tree.info.json: how the tree was made (seed, threads, sequences, bootstrap) and its
    leaves."""
    with open(path, "w") as fh:
        json.dump(info, fh, indent=1, sort_keys=True)


def run_iqtree(input_aln, output_nwk, threads=1, bootstrap_reps=1000, seed=IQTREE_SEED):
    """
    Run IQ-TREE with automatic model selection and optional ultrafast bootstrap.
    """
    prefix = output_nwk.replace(".nwk", "")
    cmd = [
        "-s",
        input_aln,
        "--prefix",
        prefix,
        "-m",
        "MFP",
        "-T",
        str(threads),
        "--quiet",
        "--redo",
        "--seed",
        str(seed),
    ]
    if bootstrap_reps and bootstrap_reps > 0:
        cmd.extend(["-B", str(bootstrap_reps)])

    try:
        subprocess.run(["iqtree2"] + cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        # Without iqtree2: `iqtree` (IQ-TREE 1, or IQ-TREE 3's own name for itself; environment.yml pins 3.1.2).
        subprocess.run(["iqtree"] + cmd, check=True, capture_output=True, text=True)

    treefile = prefix + ".treefile"
    if os.path.exists(treefile):
        shutil.move(treefile, output_nwk)
    else:
        raise RuntimeError(f"IQ-TREE completed but did not write {treefile}")

    # Cleanup IQ-TREE auxiliary files.
    for ext in [
        ".iqtree",
        ".log",
        ".mldist",
        ".model.gz",
        ".bionj",
        ".ckp.gz",
        ".contree",
        ".splits.nex",
        ".uniqueseq.phy",
    ]:
        aux = prefix + ext
        if os.path.exists(aux):
            os.remove(aux)


def main():
    parser = argparse.ArgumentParser(description="Compute phylogenetic tree from FASTA")
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument("--output", required=True, help="Output Newick tree file")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=IQTREE_SEED,
                        help="IQ-TREE's random seed (fixed by default, so a rerun gives the same tree)")
    parser.add_argument("--info", default="",
                        help="Where to write how the tree was made (default <output stem>.info.json): "
                             "seed, threads, sequences, and the bootstrap, also when it was dropped "
                             "after a failure")
    args = parser.parse_args()
    info_path = args.info or re.sub(r"\.nwk$", "", args.output) + ".info.json"
    locus = re.sub(r"_tree$", "", os.path.basename(re.sub(r"\.nwk$", "", args.output)))

    if not os.path.exists(args.input):
        print(f"Input file {args.input} not found.")
        sys.exit(1)

    entries = read_fasta_entries(args.input)
    count = len(entries)
    if count < 3:
        print("Not enough sequences to build a tree (<3).")
        write_placeholder_tree(args.input, args.output, reason="<3 sequences")
        sys.exit(0)

    aln_file = args.input + ".aln"
    print(f"Aligning {count} sequences with MAFFT...")
    try:
        run_mafft(args.input, aln_file, args.threads)
    except Exception as e:
        print(f"MAFFT failed: {e}")
        sys.exit(1)

    bootstrap_reps = 1000 if count >= 4 else None
    if bootstrap_reps:
        print("Inferring tree with IQ-TREE (auto model selection + ultrafast bootstrap)...")
    else:
        print("Inferring tree with IQ-TREE (auto model selection, no bootstrap for <4 taxa)...")

    # The leaves, as the tree names them (a header's first word): the report marks those it does
    # not count, and the tree page draws them so (fix item 46).
    leaves = [name.split()[0] for name, _seq in entries if name.split()]

    def info(status, reps, reason=""):
        write_info(info_path, locus=locus, seed=args.seed, threads=args.threads, sequences=count,
                   bootstrap={"status": status, "replicates": reps, "reason": reason}, leaves=leaves)

    try:
        run_iqtree(aln_file, args.output, args.threads, bootstrap_reps=bootstrap_reps, seed=args.seed)
        info("ultrafast" if bootstrap_reps else "none_fewer_than_4_sequences", bootstrap_reps or 0)
    except subprocess.CalledProcessError as e:
        details = ((e.stdout or "") + "\n" + (e.stderr or "")).strip()
        if bootstrap_reps:
            # Retry without bootstrap for tiny alignments or short loci.
            print("IQ-TREE failed with bootstrap; retrying without bootstrap...")
            first = (details.splitlines() or [""])[-1][:200]
            try:
                run_iqtree(aln_file, args.output, args.threads, bootstrap_reps=None, seed=args.seed)
                # Recorded (fix-phase d6): the tree page and the report say the tree has no support.
                info("dropped_after_failure", 0, f"IQ-TREE with {bootstrap_reps} ultrafast replicates "
                     f"failed (exit {e.returncode}: {first}); the tree was inferred without them")
            except subprocess.CalledProcessError as e2:
                details2 = ((e2.stdout or "") + "\n" + (e2.stderr or "")).strip()
                if is_data_limited_iqtree_error(details2):
                    write_placeholder_tree(
                        args.input,
                        args.output,
                        reason="low-information alignment",
                    )
                    info("placeholder_tree", 0, "low-information alignment")
                    sys.exit(0)
                print(f"IQ-TREE failed: {e2}")
                if details2:
                    print(details2)
                sys.exit(1)
        elif is_data_limited_iqtree_error(details):
            write_placeholder_tree(args.input, args.output, reason="low-information alignment")
            info("placeholder_tree", 0, "low-information alignment")
            sys.exit(0)
        else:
            print(f"IQ-TREE failed: {e}")
            if details:
                print(details)
            sys.exit(1)
    except Exception as e:
        print(f"IQ-TREE failed: {e}")
        sys.exit(1)

    print(f"Tree written to {args.output}")


if __name__ == "__main__":
    main()
