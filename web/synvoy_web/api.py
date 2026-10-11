"""The HTTP side: a small JSON API under /api/v1/ and the static pages.

Written against WSGI with the standard library only. ``python -m synvoy_web serve``
runs it with the built-in threaded server; behind a reverse proxy any WSGI server
(gunicorn, waitress) can run ``synvoy_web.wsgi:application`` instead.

    GET  /api/v1/status                     is the service taking searches, limits
    GET  /api/v1/library[?home=<id>]        genomes on the server (+ targets by relatedness;
                                            ?taxon=<NCBI id> orders them for a species
                                            whose genome is not here)
    POST /api/v1/lookup                     {"accession"} or {"sequence"} -> what it is
    POST /api/v1/genomes/find               {"term", "purpose"?} -> assemblies at NCBI
    POST /api/v1/searches                   start a search
    GET  /api/v1/searches/<id>              state and progress
    POST /api/v1/searches/<id>/cancel       {"token"}
    GET  /api/v1/searches/<id>/result       the result (when done)
    GET  /api/v1/searches/<id>/files/<name> one figure or download of the result
    GET  /healthz                           "ok" when the database answers

Errors are ``{"error": {"code", "message", "field"}}`` with a matching status code;
``message`` is written for the visitor and shown as it is.

The id of a search is the only key to it, so it must not leak: responses are never
cached by shared caches, pages send no referrer, and the access log shortens ids.
"""

import json
import mimetypes
import os
import re
import sys
import time
import traceback
from pathlib import Path
from urllib.parse import parse_qs

from . import db, jobs, results
from .query import QueryError

MAX_BODY_BYTES = 128 * 1024
_ID = r"(?P<id>[A-Za-z0-9][A-Za-z0-9_-]{21})"
_FILE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,120}$")

_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
    ".woff2": "font/woff2",
    ".txt": "text/plain; charset=utf-8",
    ".tsv": "text/tab-separated-values; charset=utf-8",
    ".faa": "text/plain; charset=utf-8",
    ".nwk": "text/plain; charset=utf-8",
    ".gff": "text/plain; charset=utf-8",
    ".zip": "application/zip",
    ".webmanifest": "application/manifest+json",
}
_INLINE = (".html", ".svg", ".png")

# Pages of the site itself: only our own scripts and styles, nothing from elsewhere.
SITE_CSP = ("default-src 'none'; script-src 'self'; style-src 'self'; img-src 'self' data:; "
            "font-src 'self'; connect-src 'self'{connect}; frame-src 'self'; base-uri 'none'; "
            "form-action 'self'; frame-ancestors 'none'")
# Figures written by the pipeline. They are self-contained pages with inline code, so
# they are shown in a sandbox: no access to this site's origin, no network.
FIGURE_PAGE_CSP = ("sandbox allow-scripts allow-downloads; default-src 'none'; "
                   "script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; "
                   "frame-ancestors 'self'")
FIGURE_IMAGE_CSP = "sandbox; default-src 'none'; style-src 'unsafe-inline'; img-src data:"

_STATUS_TEXT = {
    200: "OK", 204: "No Content", 304: "Not Modified", 400: "Bad Request", 403: "Forbidden",
    404: "Not Found", 405: "Method Not Allowed", 409: "Conflict", 410: "Gone",
    411: "Length Required", 413: "Payload Too Large", 415: "Unsupported Media Type",
    429: "Too Many Requests", 500: "Internal Server Error", 503: "Service Unavailable",
}


class HttpError(Exception):
    def __init__(self, status, code, message, headers=()):
        Exception.__init__(self, message)
        self.status = status
        self.code = code
        self.message = message
        self.headers = list(headers)


class Response:
    def __init__(self, status=200, body=b"", headers=None, file_path=None):
        self.status = status
        self.body = body
        self.headers = dict(headers or {})
        self.file_path = file_path


def json_response(data, status=200):
    body = json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return Response(status, body, {"Content-Type": "application/json; charset=utf-8",
                                   "Cache-Control": "no-store",
                                   "X-Robots-Tag": "noindex, nofollow"})


def error_response(status, code, message, field=None):
    return json_response({"error": {"code": code, "message": message, "field": field}}, status)


def bundle_headers(name, private):
    """Headers for a figure or download of a result."""
    suffix = Path(name).suffix.lower()
    headers = {"Content-Type": _TYPES.get(suffix, "application/octet-stream"),
               "X-Robots-Tag": "noindex, nofollow"}
    if suffix == ".html":
        headers["Content-Security-Policy"] = FIGURE_PAGE_CSP
    else:
        headers["Content-Security-Policy"] = FIGURE_IMAGE_CSP
    if suffix not in _INLINE:
        headers["Content-Disposition"] = 'attachment; filename="%s"' % Path(name).name
    headers["Cache-Control"] = "private, max-age=3600" if private else "no-cache"
    return headers


class App:
    def __init__(self, ctx, log=None):
        self.ctx = ctx
        self.settings = ctx.settings
        self.site_dir = Path(self.settings.site_dir).resolve()
        self.log = log if log is not None else self._default_log
        connect = "".join(" " + origin for origin in self.settings.cors_origins)
        self.site_csp = SITE_CSP.format(connect=connect)
        self.routes = [
            ("GET", re.compile(r"^/api/v1/status$"), self.get_status),
            ("GET", re.compile(r"^/api/v1/library$"), self.get_library),
            ("POST", re.compile(r"^/api/v1/lookup$"), self.post_lookup),
            ("POST", re.compile(r"^/api/v1/genomes/find$"), self.post_find_genomes),
            ("POST", re.compile(r"^/api/v1/searches$"), self.post_search),
            ("GET", re.compile(r"^/api/v1/searches/" + _ID + r"$"), self.get_search),
            ("POST", re.compile(r"^/api/v1/searches/" + _ID + r"/cancel$"), self.post_cancel),
            ("GET", re.compile(r"^/api/v1/searches/" + _ID + r"/result$"), self.get_result),
            ("GET", re.compile(r"^/api/v1/searches/" + _ID + r"/files/(?P<name>[^/]+)$"),
             self.get_file),
        ]

    @staticmethod
    def _default_log(message):
        sys.stderr.write("%s api: %s\n" % (time.strftime("%Y-%m-%dT%H:%M:%S"), message))
        sys.stderr.flush()

    # ── WSGI entry ──────────────────────────────────────────────────────────
    def __call__(self, environ, start_response):
        started = time.monotonic()
        method = environ.get("REQUEST_METHOD", "GET").upper()
        path = environ.get("PATH_INFO", "") or "/"
        try:
            response = self.handle(method, path, environ)
        except QueryError as exc:
            response = error_response(exc.status, exc.code, exc.message, exc.field)
        except HttpError as exc:
            response = error_response(exc.status, exc.code, exc.message)
            response.headers.update(dict(exc.headers))
        except Exception:
            incident = "%08x" % (int(time.time() * 1000) & 0xFFFFFFFF)
            self.log("incident %s on %s %s\n%s" % (incident, method, _short(path),
                                                   traceback.format_exc()))
            response = error_response(
                500, "internal_error",
                "Something went wrong on the server (reference %s). Please try again." % incident)
        self._finish_headers(response, environ, path)
        body = response.body
        send_file = response.file_path is not None and response.status == 200 and method != "HEAD"
        if response.file_path is not None and response.status == 200:
            response.headers["Content-Length"] = str(os.path.getsize(str(response.file_path)))
        else:
            response.headers["Content-Length"] = str(len(body))
        status_line = "%d %s" % (response.status, _STATUS_TEXT.get(response.status, "Status"))
        start_response(status_line, list(response.headers.items()))
        self.log("%s %s -> %d (%d ms)" % (method, _short(path), response.status,
                                          int((time.monotonic() - started) * 1000)))
        if method == "HEAD" or response.status in (204, 304):
            return [b""]
        if send_file:
            handle = open(str(response.file_path), "rb")
            wrapper = environ.get("wsgi.file_wrapper")
            return wrapper(handle, 65536) if wrapper else _iter_file(handle)
        return [body]

    def _finish_headers(self, response, environ, path):
        headers = response.headers
        headers.setdefault("X-Content-Type-Options", "nosniff")
        headers.setdefault("Referrer-Policy", "no-referrer")
        headers.setdefault("Content-Security-Policy", self.site_csp)
        headers.setdefault("Cross-Origin-Opener-Policy", "same-origin")
        headers.setdefault("Permissions-Policy", "geolocation=(), camera=(), microphone=()")
        origin = environ.get("HTTP_ORIGIN")
        if origin and origin in self.settings.cors_origins and path.startswith("/api/"):
            headers["Access-Control-Allow-Origin"] = origin
            headers["Vary"] = "Origin"

    # ── dispatch ────────────────────────────────────────────────────────────
    def handle(self, method, path, environ):
        if any(ord(ch) < 32 or ord(ch) > 126 for ch in path) or "\\" in path:
            raise HttpError(400, "bad_path", "This address is not valid.")
        if path == "/healthz":
            if self.ctx.store.quick_check() != "ok":
                raise HttpError(503, "database", "The database does not answer.")
            return Response(200, b"ok\n", {"Content-Type": "text/plain; charset=utf-8",
                                           "Cache-Control": "no-store"})
        if path.startswith("/api/"):
            if method == "OPTIONS":
                return self.preflight(environ)
            allowed = []
            for route_method, pattern, handler in self.routes:
                match = pattern.match(path)
                if not match:
                    continue
                if route_method == method or (method == "HEAD" and route_method == "GET"):
                    return handler(environ, **match.groupdict())
                allowed.append(route_method)
            if allowed:
                raise HttpError(405, "method_not_allowed", "This address does not take %s." % method,
                                [("Allow", ", ".join(sorted(set(allowed))))])
            raise HttpError(404, "not_found", "There is nothing at this address.")
        if method not in ("GET", "HEAD"):
            raise HttpError(405, "method_not_allowed", "This address does not take %s." % method,
                            [("Allow", "GET, HEAD")])
        return self.static(path, environ)

    def preflight(self, environ):
        origin = environ.get("HTTP_ORIGIN")
        headers = {"Cache-Control": "no-store"}
        if origin and origin in self.settings.cors_origins:
            headers.update({"Access-Control-Allow-Origin": origin,
                            "Access-Control-Allow-Methods": "GET, POST",
                            "Access-Control-Allow-Headers": "Content-Type",
                            "Access-Control-Max-Age": "600", "Vary": "Origin"})
        return Response(204, b"", headers)

    # ── helpers ─────────────────────────────────────────────────────────────
    def client_address(self, environ):
        if self.settings.trust_proxy:
            forwarded = environ.get("HTTP_X_FORWARDED_FOR", "")
            last = forwarded.split(",")[-1].strip()
            if last:
                return last
        return environ.get("REMOTE_ADDR", "") or "unknown"

    @staticmethod
    def read_json(environ):
        content_type = (environ.get("CONTENT_TYPE") or "").split(";")[0].strip().lower()
        if content_type != "application/json":
            raise HttpError(415, "unsupported_type", "Send the request as application/json.")
        try:
            length = int(environ.get("CONTENT_LENGTH") or "")
        except ValueError:
            raise HttpError(411, "length_required", "The request has no length.")
        if length < 0:
            raise HttpError(400, "bad_request", "The request is not in the expected form.")
        if length > MAX_BODY_BYTES:
            raise HttpError(413, "too_large", "The request is too large.")
        raw = environ["wsgi.input"].read(length) if length else b""
        if len(raw) != length:
            raise HttpError(400, "bad_request", "The request arrived incomplete.")
        try:
            return json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, ValueError):
            raise HttpError(400, "bad_json", "The request is not valid JSON.")

    def job_or_404(self, job_id):
        job = self.ctx.store.get_job(job_id)
        if job is None:
            raise HttpError(404, "not_found", "There is no search with this link.")
        return job

    def result_dir(self, job):
        """Folder of the finished result, or an HttpError that says why there is none."""
        if job["state"] == db.EXPIRED:
            raise HttpError(410, "expired", "This result was deleted after %d days. Start the "
                            "search again to get it back." % self.settings.retention_days)
        if job["state"] != db.DONE:
            raise HttpError(409, "not_ready", "This search has no result (state: %s)." % job["state"])
        folder = self.settings.jobs_dir / job["id"] / "result"
        if not (folder / "result.json").is_file():
            raise HttpError(410, "expired", "The files of this result are no longer on the server.")
        return folder

    # ── API handlers ────────────────────────────────────────────────────────
    def get_status(self, environ):
        return json_response(jobs.status_view(self.ctx))

    def get_library(self, environ):
        library = self.ctx.library()
        data = library.public()
        query = parse_qs(environ.get("QUERY_STRING", ""), max_num_fields=4)
        home_id = (query.get("home") or [""])[0]
        if home_id:
            home = library.get(home_id)
            if home is None:
                raise HttpError(404, "home_unknown", "This genome is not on the server.")
            data["home"] = home.public()
            targets = library.targets_for(home)
            data["targets"] = [g.id for g in targets]
            # ... and which of them the form ticks by default.
            data["related"] = [g.id for g in library.related(home, targets)]
        taxon = (query.get("taxon") or [""])[0]
        if taxon and not home_id:
            if not re.match(r"^[1-9][0-9]{0,9}$", taxon):
                raise HttpError(400, "bad_request", "taxon must be an NCBI taxon id.")
            data["targets"], data["related"] = jobs.targets_by_taxon(
                self.ctx, library, int(taxon), self.client_address(environ))
        response = json_response(data)
        response.headers["Cache-Control"] = "no-cache"
        return response

    def post_lookup(self, environ):
        payload = self.read_json(environ)
        return json_response(jobs.lookup_view(self.ctx, payload, self.client_address(environ)))

    def post_find_genomes(self, environ):
        payload = self.read_json(environ)
        return json_response(jobs.find_view(self.ctx, payload, self.client_address(environ)))

    def post_search(self, environ):
        payload = self.read_json(environ)
        outcome = jobs.submit(self.ctx, payload, self.client_address(environ))
        job = self.ctx.store.get_job(outcome["id"])
        view = jobs.job_view(self.ctx, job)
        view.update({"token": outcome["token"], "reused": outcome["reused"]})
        return json_response(view, 200 if outcome["reused"] else 201)

    def get_search(self, environ, id):
        return json_response(jobs.job_view(self.ctx, self.job_or_404(id)))

    def post_cancel(self, environ, id):
        payload = self.read_json(environ)
        token = payload.get("token") if isinstance(payload, dict) else None
        outcome = jobs.cancel(self.ctx, id, token)
        view = jobs.job_view(self.ctx, self.job_or_404(id))
        view.update(outcome)
        return json_response(view)

    def get_result(self, environ, id):
        folder = self.result_dir(self.job_or_404(id))
        response = Response(200, b"", {"Content-Type": "application/json; charset=utf-8",
                                       "Cache-Control": "private, no-cache",
                                       "X-Robots-Tag": "noindex, nofollow"},
                            file_path=folder / "result.json")
        return response

    def get_file(self, environ, id, name):
        if not _FILE_NAME_RE.match(name):
            raise HttpError(404, "not_found", "There is no such file.")
        folder = self.result_dir(self.job_or_404(id))
        try:
            with open(str(folder / "result.json")) as fh:
                listed = results.bundle_files(json.load(fh))
        except (OSError, ValueError):
            raise HttpError(410, "expired", "The files of this result are no longer on the server.")
        relative = "files/" + name
        path = folder / relative
        if relative not in listed or not path.is_file():
            raise HttpError(404, "not_found", "There is no such file.")
        return Response(200, b"", bundle_headers(name, private=True), file_path=path)

    # ── static pages ────────────────────────────────────────────────────────
    def static(self, path, environ):
        parts = [p for p in path.split("/") if p]
        if any(p in (".", "..") or p.startswith(".") for p in parts):
            return self.not_found()
        candidate = self.site_dir.joinpath(*parts) if parts else self.site_dir
        try:
            resolved = candidate.resolve()
            resolved.relative_to(self.site_dir)
        except (OSError, ValueError):
            return self.not_found()
        if resolved.is_dir():
            resolved = resolved / "index.html"
        if not resolved.is_file():
            if parts == ["data", "gallery.json"]:
                # The examples of this site have not been built (yet): that is an
                # empty list, not an error in the visitor's browser.
                return Response(200, b'{"examples": []}',
                                {"Content-Type": "application/json; charset=utf-8",
                                 "Cache-Control": "no-cache"})
            return self.not_found()
        suffix = resolved.suffix.lower()
        stat = resolved.stat()
        etag = '"%x-%x"' % (stat.st_mtime_ns, stat.st_size)
        in_bundle = len(parts) > 1 and parts[0] == "data"
        if in_bundle and suffix in (".html", ".svg"):
            headers = bundle_headers(resolved.name, private=False)
        else:
            headers = {"Content-Type": _TYPES.get(suffix) or mimetypes.guess_type(resolved.name)[0]
                       or "application/octet-stream", "Cache-Control": "no-cache"}
        headers["ETag"] = etag
        if environ.get("HTTP_IF_NONE_MATCH") == etag:
            return Response(304, b"", headers)
        return Response(200, b"", headers, file_path=resolved)

    def not_found(self):
        page = self.site_dir / "404.html"
        try:
            body = page.read_bytes()
        except OSError:
            body = b"<!doctype html><title>Not found</title><p>There is nothing at this address.</p>"
        return Response(404, body, {"Content-Type": "text/html; charset=utf-8",
                                    "Cache-Control": "no-store"})


def _short(path):
    """Path for the log with search ids cut to six characters: an id is a key."""
    return re.sub(r"(/searches/[A-Za-z0-9_-]{6})[A-Za-z0-9_-]{16}", r"\1..", path)


def _iter_file(handle, block=65536):
    try:
        while True:
            data = handle.read(block)
            if not data:
                break
            yield data
    finally:
        handle.close()
