"""SynVoy web service: HTTP API, job queue and worker.

Standard library only, so it runs wherever the pipeline's Python runs and its tests
need nothing beyond pytest. The pipeline itself is never imported: it is started as a
separate process and read through the files it publishes.
"""

__version__ = "0.1.0"
API_VERSION = "v1"
