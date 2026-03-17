#!/usr/bin/env python3
"""
Filter a phylogenetically sorted genome list using the aggregated QC JSON.

Input sorted list format:
    genome_name<TAB>distance

QC JSON format:
    [
      {"genome": "GCF_xxx.fna", "status": "PASS", ...},
      ...
    ]
"""

import argparse
import json
import os
import sys


def _basename_key(value: str) -> str:
    return os.path.basename((value or "").strip())


def _load_qc_statuses(qc_json_path: str) -> dict[str, dict]:
    with open(qc_json_path) as fh:
        data = json.load(fh)

    statuses = {}
    for rec in data:
        genome = _basename_key(rec.get("genome", ""))
        if genome:
            statuses[genome] = rec
    return statuses


def main():
    parser = argparse.ArgumentParser(
        description="Filter sorted genome list using SynTerra QC summary"
    )
    parser.add_argument("--sorted", required=True, help="Sorted genomes TSV")
    parser.add_argument("--qc_json", required=True, help="Aggregated QC JSON")
    parser.add_argument("--output", required=True, help="Filtered output TSV")
    parser.add_argument(
        "--policy",
        choices=["drop", "keep"],
        default="drop",
        help="How to handle QC failures (default: drop)",
    )
    args = parser.parse_args()

    qc_statuses = _load_qc_statuses(args.qc_json)

    kept_lines = []
    dropped = []
    unknown = []
    original_lines = []

    with open(args.sorted) as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            original_lines.append(line)
            parts = line.split("\t")
            genome_name = parts[0]
            key = _basename_key(genome_name)
            qc_rec = qc_statuses.get(key)
            if qc_rec is None:
                kept_lines.append(line)
                unknown.append(key)
                continue

            status = str(qc_rec.get("status", "PASS")).upper()
            if status == "FAIL" and args.policy == "drop":
                dropped.append(
                    {
                        "genome": key,
                        "msg": qc_rec.get("msg", "QC failure"),
                    }
                )
                continue

            kept_lines.append(line)

    # Avoid converting QC into a hard pipeline failure when every target fails.
    # Keep the original ordering in that case, but log loudly.
    if original_lines and not kept_lines and args.policy == "drop":
        kept_lines = list(original_lines)
        print(
            "[qc-filter] All genomes failed QC; keeping the original sorted list "
            "to avoid a full no-target abort.",
            file=sys.stderr,
        )

    with open(args.output, "w") as out:
        for line in kept_lines:
            out.write(line + "\n")

    print(
        f"[qc-filter] policy={args.policy} kept={len(kept_lines)} "
        f"dropped={len(dropped)} unknown={len(unknown)}",
        file=sys.stderr,
    )
    for rec in dropped:
        print(
            f"[qc-filter] dropped {rec['genome']}: {rec['msg']}",
            file=sys.stderr,
        )
    if unknown:
        print(
            f"[qc-filter] kept {len(unknown)} genome(s) without QC record: "
            + ", ".join(sorted(set(unknown))),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
