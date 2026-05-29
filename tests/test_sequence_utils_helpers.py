#!/usr/bin/env python3
"""
Tests for the shared helpers in sequence_utils.py: parse_bed, str2bool,
run_command, require_binary, setup_logging, read_json/write_json,
iter_fasta_headers, and the upgraded parse_gff_attributes (full URL unquote).
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))

from sequence_utils import (  # noqa: E402
    iter_fasta_headers,
    parse_bed,
    parse_gff_attributes,
    read_json,
    require_binary,
    run_command,
    setup_logging,
    str2bool,
    write_json,
)


class TestParseBed(unittest.TestCase):
    def test_missing_file_returns_empty(self):
        self.assertEqual(parse_bed('/nonexistent/path.bed'), [])
        self.assertEqual(parse_bed(''), [])
        self.assertEqual(parse_bed(None), [])

    def test_empty_file_returns_empty(self):
        with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
            path = f.name
        try:
            self.assertEqual(parse_bed(path), [])
        finally:
            os.unlink(path)

    def test_skips_comments_and_blanks(self):
        with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
            f.write("# header comment\n")
            f.write("\n")
            f.write("chr1\t100\t200\tgene1\t.\t+\n")
            path = f.name
        try:
            rows = parse_bed(path)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['chrom'], 'chr1')
            self.assertEqual(rows[0]['start'], 100)
            self.assertEqual(rows[0]['end'], 200)
            self.assertEqual(rows[0]['name'], 'gene1')
            self.assertEqual(rows[0]['strand'], '+')
            self.assertEqual(rows[0]['display_name'], '')
        finally:
            os.unlink(path)

    def test_three_column_bed_synthesises_name(self):
        with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
            f.write("chr2\t500\t1500\n")
            path = f.name
        try:
            rows = parse_bed(path)
            self.assertEqual(rows[0]['name'], 'chr2:500-1500')
            self.assertEqual(rows[0]['strand'], '+')
        finally:
            os.unlink(path)

    def test_seven_column_keeps_display_name(self):
        with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
            f.write("chrX\t10\t20\tg1\t.\t-\tFancyGeneName\n")
            path = f.name
        try:
            rows = parse_bed(path)
            self.assertEqual(rows[0]['strand'], '-')
            self.assertEqual(rows[0]['display_name'], 'FancyGeneName')
        finally:
            os.unlink(path)

    def test_invalid_coords_are_skipped(self):
        with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as f:
            f.write("chr1\tnope\t100\tg\n")
            f.write("chr1\t10\t20\tgood\n")
            path = f.name
        try:
            rows = parse_bed(path)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['name'], 'good')
        finally:
            os.unlink(path)


class TestStr2Bool(unittest.TestCase):
    def test_truthy_values(self):
        for v in ['true', 'True', 'TRUE', '1', 'yes', 'y', 't', True]:
            self.assertTrue(str2bool(v), f"expected True for {v!r}")

    def test_falsy_values(self):
        for v in ['false', 'False', '0', 'no', 'n', 'f', False, None]:
            self.assertFalse(str2bool(v), f"expected False for {v!r}")

    def test_invalid_raises_argparse_error(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            str2bool('maybe')

    def test_whitespace_tolerated(self):
        self.assertTrue(str2bool('  yes  '))
        self.assertFalse(str2bool('  no\n'))


class TestRunCommand(unittest.TestCase):
    def test_success_returns_completed_process(self):
        result = run_command(['true'])
        self.assertEqual(result.returncode, 0)

    def test_failure_raises_when_check_true(self):
        with self.assertRaises(subprocess.CalledProcessError):
            run_command(['false'])

    def test_failure_returns_when_check_false(self):
        result = run_command(['false'], check=False)
        self.assertNotEqual(result.returncode, 0)

    def test_capture_returns_stdout(self):
        result = run_command(['echo', 'hello'], capture=True)
        self.assertEqual(result.stdout.strip(), 'hello')


class TestRequireBinary(unittest.TestCase):
    def test_present_binary_returns_path(self):
        # 'sh' is on every POSIX system
        path = require_binary('sh')
        self.assertTrue(os.path.isabs(path))
        self.assertTrue(os.path.exists(path))

    def test_missing_binary_raises(self):
        with self.assertRaises(FileNotFoundError):
            require_binary('definitely_not_a_real_binary_xyz_123')


class TestSetupLogging(unittest.TestCase):
    def test_returns_logger_at_level(self):
        log = setup_logging(logging.DEBUG, name='synvoy.test.setup_logging')
        self.assertEqual(log.level, logging.DEBUG)

    def test_idempotent_no_duplicate_handlers(self):
        before = len(logging.getLogger().handlers)
        setup_logging()
        setup_logging()
        after = len(logging.getLogger().handlers)
        # Handler count must not grow beyond what the first call established
        self.assertLessEqual(after, max(1, before))


class TestJsonIO(unittest.TestCase):
    def test_round_trip(self):
        obj = {'a': 1, 'b': [1, 2, 3], 'c': {'nested': True}}
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'out.json'
            write_json(path, obj)
            self.assertEqual(read_json(path), obj)

    def test_atomic_write_no_partial_file_on_replace(self):
        # After write_json completes, no .tmp sibling should remain.
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'sub' / 'out.json'
            write_json(path, {'x': 1})
            self.assertTrue(path.exists())
            self.assertFalse(path.with_suffix(path.suffix + '.tmp').exists())

    def test_creates_parent_dirs(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'deep' / 'tree' / 'out.json'
            write_json(path, {'ok': True})
            self.assertEqual(json.loads(path.read_text()), {'ok': True})


class TestIterFastaHeaders(unittest.TestCase):
    def test_yields_clean_headers(self):
        with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as f:
            f.write(">seq1 description one\n")
            f.write("AAAA\n")
            f.write(">seq2\n")
            f.write("CCCC\n")
            path = f.name
        try:
            headers = list(iter_fasta_headers(path))
            self.assertEqual(headers, ['seq1 description one', 'seq2'])
        finally:
            os.unlink(path)


class TestParseGffAttributesUnquote(unittest.TestCase):
    """Full urllib unquote should handle %20 (spaces) and other escapes that
    the old replace-loop missed."""

    def test_handles_space_escape(self):
        attrs = parse_gff_attributes("ID=foo;product=DNA%20polymerase%20III")
        self.assertEqual(attrs['product'], 'DNA polymerase III')

    def test_handles_semicolon_escape(self):
        attrs = parse_gff_attributes("Name=geneA%3BgeneB")
        self.assertEqual(attrs['Name'], 'geneA;geneB')

    def test_empty_input(self):
        self.assertEqual(parse_gff_attributes(''), {})
        self.assertEqual(parse_gff_attributes('.'), {})


if __name__ == '__main__':
    unittest.main()
