"""Entry point for a WSGI server, e.g. ``gunicorn synvoy_web.wsgi:application``.

Settings come from ``$SYNVOY_WEB_CONFIG`` and ``SYNVOY_WEB_*`` variables.
"""

from .server import build_app

application = build_app()
