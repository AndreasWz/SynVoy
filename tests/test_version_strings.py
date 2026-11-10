"""Fix item 7 (9 Nov 2026): the run banner printed "SynVoy v2.0" while VERSION and the manifest said 1.0.0 (since
6dd1d18, 2026-10-13); it reads the manifest's version now, and VERSION and the manifest must agree."""
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_the_banner_reads_the_manifest_version():
    main = open(os.path.join(ROOT, "main.nf")).read()
    body = main[main.index("def printHeader()"):]
    body = body[:body.index("\n}\n")]
    assert "workflow.manifest.version" in body
    assert not re.search(r"\bv\d+\.\d+", body)


def test_the_version_file_equals_the_manifest():
    version = open(os.path.join(ROOT, "VERSION")).read().strip()
    cfg = open(os.path.join(ROOT, "nextflow.config")).read()
    m = re.search(r"manifest\s*\{[^}]*?\bversion\s*=\s*'([^']+)'", cfg, re.S)
    assert m and m.group(1) == version
