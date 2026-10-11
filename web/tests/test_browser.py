"""The browser run-through (``tests/browser/run.mjs``) as a test. Off by default: it
needs Chrome and Node 22 or newer and takes about a minute.

    SYNVOY_WEB_BROWSER=1 python -m pytest tests/test_browser.py -s
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

RUN = Path(__file__).with_name("browser") / "run.mjs"


def _node():
    node = os.environ.get("SYNVOY_WEB_NODE") or shutil.which("node")
    if not node:
        return None
    try:
        version = subprocess.run([node, "--version"], capture_output=True, text=True,
                                 timeout=20).stdout.strip().lstrip("v")
        return node if int(version.split(".")[0]) >= 22 else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


@pytest.mark.skipif(os.environ.get("SYNVOY_WEB_BROWSER") != "1",
                    reason="set SYNVOY_WEB_BROWSER=1 to drive a real browser")
def test_the_site_in_a_real_browser():
    node = _node()
    chrome = os.environ.get("SYNVOY_WEB_CHROME") or shutil.which("google-chrome")
    if not node or not chrome:
        pytest.skip("needs Node 22 or newer (SYNVOY_WEB_NODE) and Chrome (SYNVOY_WEB_CHROME)")
    env = dict(os.environ, SYNVOY_WEB_PYTHON=sys.executable, SYNVOY_WEB_CHROME=chrome)
    env.pop("VIRTUAL_ENV", None)
    env.pop("SYNVOY_WEB_CONFIG", None)
    done = subprocess.run([node, str(RUN)], capture_output=True, text=True, timeout=900, env=env)
    print(done.stdout[-6000:])
    assert done.returncode == 0, done.stdout[-3000:] + done.stderr[-2000:]
