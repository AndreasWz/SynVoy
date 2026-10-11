"""Putting the parts together, and the built-in HTTP server.

The built-in server (wsgiref with one thread per request) is enough for trying the
site and for a small installation behind a reverse proxy. It speaks plain HTTP only:
TLS, request-rate limits and slow-client protection belong in the proxy
(``deploy/Caddyfile``).
"""

import socketserver
import threading
from wsgiref.simple_server import WSGIRequestHandler, WSGIServer

from .api import App
from .config import load_settings
from .db import Store
from .jobs import Context
from .library import LibraryWatcher


def build_context(settings=None):
    settings = settings or load_settings()
    store = Store(settings.db_path)
    return Context(settings, store, LibraryWatcher(settings.library_path))


def build_app(settings=None, ctx=None):
    return App(ctx or build_context(settings))


class _Handler(WSGIRequestHandler):
    timeout = 30                      # a client that stops sending is dropped

    def log_message(self, *args):     # the app writes its own access log, without addresses
        pass


class ThreadedServer(socketserver.ThreadingMixIn, WSGIServer):
    daemon_threads = True
    request_queue_size = 64


def make_server(app, host, port):
    server = ThreadedServer((host, port), _Handler)
    server.set_app(app)
    return server


def serve_in_thread(app, host="127.0.0.1", port=0):
    """Start the server on a background thread; returns ``(server, thread)``. With
    ``port=0`` the system picks a free port (``server.server_address[1]``)."""
    server = make_server(app, host, port)
    thread = threading.Thread(target=server.serve_forever, name="synvoy-web-http", daemon=True)
    thread.start()
    return server, thread
