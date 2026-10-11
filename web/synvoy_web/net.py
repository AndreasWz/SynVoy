"""Outbound HTTP, used for three things only: looking up a protein accession, reading
the taxonomy of a genome, and getting a genome assembly from NCBI.

Every request goes to one of the hosts below over https, has a timeout and (except for
the genome files, whose size is checked by the caller) a size cap, and a redirect to
any other host is refused. The accession that ends up in a URL has been matched
against a strict pattern before, so a visitor cannot make the server fetch an address
of their choosing.
"""

import socket
import threading
import time
import urllib.error
import urllib.request
from urllib.parse import urlsplit

ALLOWED_HOSTS = ("rest.uniprot.org", "eutils.ncbi.nlm.nih.gov", "api.ncbi.nlm.nih.gov",
                 "ftp.ncbi.nlm.nih.gov")
# Hosts whose answers are small look-ups: NCBI asks to keep these to a few per second.
_THROTTLED_HOSTS = ("eutils.ncbi.nlm.nih.gov", "api.ncbi.nlm.nih.gov")
USER_AGENT = "SynVoy-web (+https://github.com/AndreasWz/SynVoy)"

# NCBI asks for at most 3 requests per second without an API key.
_NCBI_MIN_INTERVAL = 0.4
_ncbi_lock = threading.Lock()
_ncbi_last = [0.0]


class FetchError(Exception):
    """A look-up did not give a usable answer.

    ``kind`` is one of: not_found, timeout, unavailable, too_large, refused.
    """

    def __init__(self, kind, detail=""):
        Exception.__init__(self, "%s: %s" % (kind, detail) if detail else kind)
        self.kind = kind
        self.detail = detail


def _check_url(url):
    parts = urlsplit(url)
    if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:
        raise FetchError("refused", "not an allowed address: %s" % url)
    return parts


class _SameHostsOnly(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _check_url(newurl)
        return urllib.request.HTTPRedirectHandler.redirect_request(
            self, req, fp, code, msg, headers, newurl)


_opener = urllib.request.build_opener(_SameHostsOnly())


def http_get(url, timeout=15, max_bytes=2000000, accept="*/*"):
    """Body of ``url`` as bytes, or FetchError."""
    parts = _check_url(url)
    if parts.hostname in _THROTTLED_HOSTS:
        with _ncbi_lock:
            wait = _NCBI_MIN_INTERVAL - (time.monotonic() - _ncbi_last[0])
            if wait > 0:
                time.sleep(wait)
            _ncbi_last[0] = time.monotonic()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    try:
        with _opener.open(request, timeout=timeout) as response:
            data = response.read(max_bytes + 1)
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 404, 410):
            raise FetchError("not_found", "HTTP %d" % exc.code)
        raise FetchError("unavailable", "HTTP %d" % exc.code)
    except (socket.timeout, TimeoutError):
        raise FetchError("timeout", "no answer within %d s" % timeout)
    except (urllib.error.URLError, ConnectionError, OSError) as exc:
        reason = getattr(exc, "reason", exc)
        if isinstance(reason, (socket.timeout, TimeoutError)):
            raise FetchError("timeout", "no answer within %d s" % timeout)
        raise FetchError("unavailable", str(reason))
    if len(data) > max_bytes:
        raise FetchError("too_large", "more than %d bytes" % max_bytes)
    return data


def http_stream(url, offset=0, timeout=60):
    """Open ``url`` for reading from byte ``offset``. Returns ``(response, total)``:
    ``response.read(n)`` gives the bytes, ``response.partial`` says whether the server
    honoured the offset, ``total`` is the size of the whole file (None when the server
    does not say). The caller closes the response. Raises FetchError."""
    _check_url(url)
    headers = {"User-Agent": USER_AGENT, "Accept": "*/*", "Accept-Encoding": "identity"}
    if offset:
        headers["Range"] = "bytes=%d-" % offset
    request = urllib.request.Request(url, headers=headers)
    try:
        response = _opener.open(request, timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 416:                      # nothing after the offset: all is here
            raise FetchError("complete", "HTTP 416")
        if exc.code in (400, 404, 410):
            raise FetchError("not_found", "HTTP %d" % exc.code)
        raise FetchError("unavailable", "HTTP %d" % exc.code)
    except (socket.timeout, TimeoutError):
        raise FetchError("timeout", "no answer within %d s" % timeout)
    except (urllib.error.URLError, ConnectionError, OSError) as exc:
        raise FetchError("unavailable", str(getattr(exc, "reason", exc)))
    total = None
    response.partial = response.status == 206
    content_range = response.headers.get("Content-Range") or ""
    if response.partial and "/" in content_range:
        try:
            total = int(content_range.rsplit("/", 1)[1])
        except ValueError:
            total = None
    elif response.headers.get("Content-Length"):
        try:
            total = int(response.headers["Content-Length"])
        except ValueError:
            total = None
    return response, total
