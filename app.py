"""
The Genome Browser service for Trellis Tools.

A Flask application, backed by its own SQLite database, that lets lab members
load reference genomes, build a library of aligned data tracks, and view them
in an IGV-style browser. It follows the suite's conventions: it keeps no
accounts or labs of its own, reading the shared directory (MoClo's ``users``
table) and ``suite.db`` instead; it registers in ``landing/tools.json`` and
exposes a ``/health`` endpoint; and it runs behind the proxy with
``SUITE_AUTH=proxy``, trusting the ``X-Suite-User`` header as the signed-in
member.

The proxy mounts this service under ``/genome`` and strips that prefix before
passing a request on, so the routes here live at ``/health``, ``/`` and
``/api/...``; ``URL_PREFIX`` is kept for redirects and links that have to put
the prefix back.

This module is the application factory and its configuration. The data model,
the readers, and the endpoints live in their own modules; only the wiring is
here.
"""

import os
from pathlib import Path

from flask import Flask, jsonify

import db
from figures import register_figures
from genomes import register_genomes
from tracks import register_tracks
from regions import register_regions
from notebook import register_notebook
from tags import register_tags
from views import register_views

APP_VERSION = '2.34.4'


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#
# Everything the module needs to find the shared stores, reach the other tools,
# and bound conversion comes from the environment, so compose is the single
# place it is set. The defaults match the compose file in the design, so the
# module also runs with nothing set for development and tests.

def _flag(name: str, default: str = '') -> bool:
    return os.environ.get(name, default).strip().lower() in ('1', 'true', 'yes')


class Config:
    """The module's settings, read once from the environment."""

    # Where the proxy mounts this service. It strips the prefix before passing a
    # request on, so redirects and generated links have to put it back.
    URL_PREFIX = '/' + os.environ.get('URL_PREFIX', '/genome').strip('/')

    # Inside the suite the proxy says who is signed in (X-Suite-User); the
    # module never implements sign-in of its own.
    SUITE_AUTH = os.environ.get('SUITE_AUTH', '').strip().lower() == 'proxy'

    # The shared stores, read-only. suite.db holds labs and chair groups
    # (labs.py); the MoClo users table is the account directory (directory.py).
    SUITE_DATABASE_PATH = Path(os.environ.get('SUITE_DATABASE_PATH', '/suite/suite.db'))
    MOCLO_DB = Path(os.environ.get('MOCLO_DB', '/moclo/moclo.db'))

    # The other tools, for tag validation, links, and direct-to-notebook export.
    NOTEBOOK_URL = os.environ.get('NOTEBOOK_URL', '/notebook')
    MOCLO_URL = os.environ.get('MOCLO_URL', '/')
    EDITOR_URL = os.environ.get('EDITOR_URL', '/editor')
    SEEDS_URL = os.environ.get('SEEDS_URL', '/seeds')

    # This module's own data root (the genome-data volume): the SQLite
    # database, the per-genome and per-track files, and the conversion scratch.
    DATA_ROOT = Path(os.environ.get('GENOME_DATA_ROOT', '/genome'))

    # Conversion guardrails. The upload ceiling is enforced in the app before a
    # job is created; the free-disk floor gates a job before heavy work; the
    # thread cap keeps samtools from starving the other tools.
    MAX_BAM_GB = int(os.environ.get('MAX_BAM_GB', '50'))
    CONVERT_MIN_FREE_GB = int(os.environ.get('CONVERT_MIN_FREE_GB', '50'))
    SAMTOOLS_THREADS = int(os.environ.get('SAMTOOLS_THREADS', '2'))


# ---------------------------------------------------------------------------
# Application factory
# ---------------------------------------------------------------------------

def create_app(config: Config = None) -> Flask:
    """
    Build the Flask application.

    The data model is created on startup (tables if absent), so a fresh volume
    comes up ready. Endpoints are registered by their own modules as they are
    added; for now the factory wires configuration, the data root, and health.
    """
    config = config or Config()
    app = Flask(__name__, static_folder=None)
    app.config.from_object(config)

    # Make sure the data root and its sub-directories exist before anything
    # tries to open the database or stream an upload into scratch.
    for sub in ('', 'genomes', 'tracks', 'scratch'):
        (config.DATA_ROOT / sub).mkdir(parents=True, exist_ok=True)

    # Create the module's own tables if they are not there yet.
    db.init_db(config.DATA_ROOT / 'genome.db')

    register_health(app, config)
    register_genomes(app)
    register_tracks(app)
    register_figures(app)
    register_regions(app)
    register_notebook(app)
    register_tags(app)
    register_views(app)
    return app


def register_health(app: Flask, config: Config) -> None:
    """
    The health endpoint, reachable without sign-in (nginx leaves it public).

    200 when the module can reach its SQLite database, 503 when it cannot, so
    the proxy and any monitoring can tell the service apart from its store.
    (Requirements 1.1, 1.2, 1.3.)
    """

    @app.route('/health')
    def health():
        reachable = db.database_reachable(config.DATA_ROOT / 'genome.db')
        body = {'status': 'healthy' if reachable else 'unhealthy',
                'version': APP_VERSION, 'database': reachable}
        return jsonify(body), (200 if reachable else 503)


# The WSGI entry point gunicorn runs (``app:app``).
app = create_app()


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5008)
