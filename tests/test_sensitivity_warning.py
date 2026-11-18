"""The search log's MMseqs2 sensitivity warning uses the start's bounds (bug list item 28, 2026-11-16).

The search warned about every sensitivity above 9, the pipeline's own default of 9.5 included, while
the start (main.nf) warns only outside 1-12. Both now warn outside 1-12, with the same text.
"""
import os
import re
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))

import iterative_search_runner as isr  # noqa: E402


@pytest.mark.parametrize("sens", [1.0, 7.0, 7.5, 9.5, 12.0])
def test_no_warning_inside_the_bounds(sens):
    assert isr.mmseqs_sensitivity_warning(sens) is None


@pytest.mark.parametrize("sens", [0.5, 12.5, 13.0])
def test_a_warning_outside_them(sens):
    assert isr.mmseqs_sensitivity_warning(sens) == f"Unusual MMseqs sensitivity {sens}. Typical range is 1-9.5."


def test_the_bounds_are_the_starts():
    with open(os.path.join(ROOT, "main.nf")) as fh:
        main = fh.read()
    m = re.search(r"mmseqsSensitivity < ([0-9.]+) \|\| mmseqsSensitivity > ([0-9.]+)", main)
    assert m, "main.nf's sensitivity check not found"
    assert isr.MMSEQS_SENSITIVITY_TYPICAL == (float(m.group(1)), float(m.group(2)))
    assert "Typical range is 1-9.5." in main


def test_the_default_does_not_warn():
    with open(os.path.join(ROOT, "nextflow.config")) as fh:
        config = fh.read()
    default = float(re.search(r"^\s*mmseqs_sensitivity\s*=\s*([0-9.]+)", config, re.M).group(1))
    assert isr.mmseqs_sensitivity_warning(default) is None


def test_main_uses_the_helper():
    src = open(os.path.join(ROOT, "bin", "iterative_search_runner.py")).read()
    body = src.split("\ndef main(", 1)[1]
    assert "mmseqs_sensitivity_warning(args.mmseqs_sens)" in body
    assert "outside typical range" not in src
