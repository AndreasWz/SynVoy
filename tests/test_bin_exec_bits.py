r"""Regression guard: every bin/*.py script invoked PATH-relative from a
modules/*.nf process MUST have the +x bit set.

Why: Nextflow process scripts call helpers as e.g. `profile_hits.py …`,
relying on `beforeScript = "export PATH=${projectDir}/bin:\$PATH"` + the
file's exec bit. A new helper created via the editor without `chmod +x`
fails the invocation with ``Permission denied`` — and the surrounding
``errorStrategy 'ignore'`` / ``|| echo "…skipped"`` patterns SWALLOW the
error, so the step silently no-ops in production.

This bit us four times: the §1f advisory (profile_hits.py +
auto_select_preset.py from the 2026-05-27 landing) never actually ran on
any pipeline invocation until the missing exec bit was caught on
2026-05-29; rescue_strong_synteny.py (§1e) and
reciprocal_best_paralog_check.py (§1j Phase B) shipped the same way in
the 2026-05-29 batch.

Scripts loaded via ``python3 bin/foo.py`` or imported as Python modules
(``import gene_predictor``) don't need +x — only those invoked
PATH-relative inside a process ``script:`` block. We discover that list
by grepping modules/*.nf for ``name.py`` calls; anything in bin/ with the
same basename must be executable.
"""
from __future__ import annotations

import os
import re
import stat
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
MODULES_DIR = ROOT / "modules"


def _path_relative_script_invocations() -> set[str]:
    """Scan modules/*.nf for bare-name script calls like ``foo.py …``.

    We pick up names that appear at the start of a script-line (after
    leading whitespace) or after a continuation `\\` or `&&`. We
    explicitly exclude lines where the script is called via ``python3
    bin/foo.py`` or ``${projectDir}/bin/foo.py`` — those don't need +x.
    """
    pat = re.compile(r"(?<![/\w.-])(\w[\w_]*\.py)\b")
    excluded = set()
    found: set[str] = set()
    for nf in MODULES_DIR.glob("*.nf"):
        for raw in nf.read_text().splitlines():
            line = raw.strip()
            # Skip comments, output declarations, etc.
            if not line or line.startswith("//"):
                continue
            for m in pat.finditer(raw):
                # If preceded immediately by "python3 ", "python ",
                # "bin/", or "${projectDir}/bin/", treat as explicit path
                # (no +x needed).
                start = m.start()
                preceding = raw[max(0, start - 30):start]
                if ("python3 " in preceding[-15:] or
                        "python " in preceding[-15:] or
                        "/bin/" in preceding[-10:]):
                    excluded.add(m.group(1))
                    continue
                found.add(m.group(1))
    return found - excluded


PATH_INVOKED = sorted(_path_relative_script_invocations())


def test_path_invoked_scripts_are_discovered():
    """Sanity check: we should find at least the known broken-on-2026-05-29
    helpers (profile_hits, rescue_strong_synteny, reciprocal_best_paralog_check)
    in the discovery sweep. If this assertion fails it means the discovery
    regex was tightened in a way that drops real cases — fix it instead of
    relaxing the test."""
    must_include = {
        "profile_hits.py", "auto_select_preset.py",
        "rescue_strong_synteny.py", "reciprocal_best_paralog_check.py",
    }
    missing = must_include - set(PATH_INVOKED)
    assert not missing, (
        f"Discovery regex missed known PATH-relative invocations: {sorted(missing)}. "
        f"Found: {PATH_INVOKED}"
    )


@pytest.mark.parametrize("script_name", PATH_INVOKED)
def test_path_invoked_script_is_executable(script_name: str):
    """For each helper that a Nextflow process invokes PATH-relative, the
    file in ``bin/`` must have the +x bit set for owner. Without it the
    process fails with ``Permission denied`` — and most §1f-style
    advisory/§1e rescue/§1j paralog processes silently swallow that
    failure via ``errorStrategy 'ignore'``."""
    path = BIN_DIR / script_name
    if not path.exists():
        # Some helpers live in conf/ or are dynamically generated;
        # only enforce the +x rule for scripts actually shipped in bin/.
        pytest.skip(f"{script_name} not present in bin/ — likely external/dynamic")
    mode = path.stat().st_mode
    assert mode & stat.S_IXUSR, (
        f"bin/{script_name} is invoked PATH-relative from a Nextflow process "
        f"but lacks the +x bit (mode={oct(mode)}). Fix: chmod +x bin/{script_name}. "
        f"Symptom in production: the calling process logs 'Permission denied' but "
        f"errorStrategy 'ignore' or '|| echo \"…skipped\"' swallows the failure, so "
        f"the step silently no-ops."
    )


def test_no_orphaned_executable_bit_on_unused_scripts():
    """Soft inverse: scripts NOT invoked PATH-relative and NOT imported as
    Python modules don't need +x. This isn't enforced (a script may grow
    callers later) — it's just a heads-up if you see surprising +x bits.
    Currently a no-op assertion; kept here so future audits have a hook."""
    # No-op for now; presence-of-test documents the intent.
    assert True
