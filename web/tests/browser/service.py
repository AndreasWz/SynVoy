"""A throw-away installation for the browser run-through: the real pages, the
stand-in pipeline, a made-up genome library, canned accession look-ups, and a stand-in
NCBI that has two more genomes and hands them out slowly.

    python service.py <folder> <port>

Prints "ready" when it can be used and runs until it is killed.
"""

import os
import sys
import threading
import time
from pathlib import Path

TESTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TESTS.parent))
sys.path.insert(0, str(TESTS))
os.environ.setdefault("FAKE_SLOW_SECONDS", "900")

from conftest import Env  # noqa: E402
from synvoy_web.api import App  # noqa: E402
from synvoy_web.server import make_server  # noqa: E402


def main():
    root, port = Path(sys.argv[1]), int(sys.argv[2])
    root.mkdir(parents=True, exist_ok=True)
    env = Env(root, site_dir=str(TESTS.parent / "site"), max_running=2, port=port,
              genome_store_gb=0.001)
    # Two genomes that are at "NCBI" and not on the server. The download takes some
    # seconds, so the page of a waiting search can be seen saying what it waits for.
    env.ncbi.add("GCA_000000101.1", "Formica nova", size=400000, lineage=[1, 2, 3, 15])
    env.ncbi.add("GCF_000000102.1", "Gigas maximus", reported=5000000000)
    env.ncbi.chunk, env.ncbi.delay = 16, 0.15
    worker = env.worker(log=lambda message: None)
    server = make_server(App(env.ctx, log=lambda message: None), "127.0.0.1", port)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    worker.tick()
    print("ready", flush=True)
    while True:
        worker.tick()
        time.sleep(0.2)


if __name__ == "__main__":
    main()
