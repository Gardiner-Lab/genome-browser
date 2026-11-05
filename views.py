"""
The browser page, the vendored static assets, and Browser_View persistence.

This module owns the three things the interactive browser needs that are not a
data feed: the HTML page that boots igv.js (``GET /``, served at ``/genome/``
through the proxy), the static route that serves the vendored client assets
(``/static/<path>``, the browser-visible ``/genome/static/<path>``), and the
Browser_View persistence endpoint (``/api/view``) that remembers, per member
and genome, which locus is shown, which tracks are enabled and in what order,
and each track's rendering options.

Three conventions from the rest of the module shape it:

- **The proxy strips ``/genome``.** nginx mounts the service under ``/genome``
  and removes that prefix before the request reaches Flask, so the routes here
  are ``/``, ``/static/<path>`` and ``/api/view`` — *without* the prefix. The
  browser, on the other hand, still sees the prefix, so every asset URL the
  page emits (the igv.js script, anything it loads) is written with the
  ``/genome/static/...`` prefix put back. ``Config.URL_PREFIX`` is the single
  source of that prefix.

- **The app was built with ``static_folder=None``.** ``app.py`` deliberately
  disables Flask's built-in static handler, so this module registers its own
  ``/static/<path>`` route over the vendored ``genome/static/`` directory, with
  a path-traversal guard so a crafted ``..`` cannot escape that directory.

- **Access and identity are never reinvented.** ``/api/view`` requires a
  signed-in member (:func:`identity.require_member`, proxy-trusted) and only
  touches a genome the member may view (:func:`access.can_view_genome`, else a
  404 that reveals nothing). The persisted enabled-track set is filtered to the
  tracks that actually belong to the genome and are viewable by the member, so
  a view can never name a track out of reach. The locus is validated and
  clamped against the genome's ``chrom_sizes`` through the pure :mod:`locus`
  math, so an unknown chromosome is rejected naming it and an out-of-range
  window snaps to the chromosome bounds — exactly the region endpoint's rules.

The Browser_View row is a single upsert keyed by ``(username, genome_id)``
(Requirements 6.1, 6.2, 6.3, 8.8): enabling a track adds it to ``enabled_json``
(order preserved), disabling it removes it, reordering rewrites the order, and
the per-track rendering options are stored in ``options_json`` for the client
and the SVG renderer to read back. The endpoint stores and returns the
normalised view; it does not itself draw anything.

Wired by :func:`register_views`, mirroring the sibling ``register_*``
registrars so the app factory stays a thin wiring layer.
"""

import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional

from flask import (Flask, Response, abort, current_app, jsonify, request,
                   send_from_directory)

import access
import db
import identity
import locus as locus_math

# The vendored client assets live beside this module in ``genome/static/``.
# app.py built the app with static_folder=None, so this module serves them.
STATIC_DIR = Path(__file__).with_name('static')


def register_views(app: Flask) -> None:
    """
    Wire the browser page, the static route, and ``/api/view`` onto the app.

    Called once from ``create_app`` with a single line, mirroring the other
    ``register_*`` registrars. Keeping the routes here lets this module own the
    page, its assets, and the Browser_View persistence without touching the
    factory body or the sibling endpoint modules. The shared
    ``identity.Unauthenticated`` 401 handler registered in ``genomes.py`` covers
    the unauthenticated case for ``/api/view``.
    """

    @app.route('/')
    def browser_page():
        """
        The igv.js browser page (served at ``/genome/`` through the proxy).

        Returns the self-contained HTML shell that boots igv.js from the
        vendored ``/genome/static/igv.min.js`` and mounts an IGV container. The
        page is open (no per-member gate) the way a tool's landing page is; the
        data it loads goes through the authenticated API endpoints.
        """
        return Response(_browser_html(_static_prefix()), mimetype='text/html')

    @app.route('/static/<path:asset>')
    def static_asset(asset):
        """
        Serve a vendored client asset from ``genome/static/`` (``/genome/static/``).

        app.py disabled Flask's built-in static handler, so this module serves
        the directory itself. ``send_from_directory`` resolves the path against
        ``STATIC_DIR`` and refuses anything that escapes it, so a crafted
        ``..`` cannot read files outside the vendored assets. A missing asset is
        a plain 404.
        """
        if not _safe_asset(asset):
            abort(404)
        try:
            return send_from_directory(STATIC_DIR, asset)
        except Exception:
            abort(404)

    @app.route('/api/view', methods=['GET'])
    def get_view():
        """
        Read back the persisted Browser_View for a genome (helper read).

        Requires a signed-in member and a viewable genome (else 404). Returns
        the stored locus, enabled track ids (in order), and rendering options;
        when the member has no saved view for the genome yet, returns an empty
        view (no locus, no enabled tracks, no options) rather than a 404, so the
        client can treat "never saved" and "saved empty" the same way.
        """
        member = identity.require_member()
        genome_id = (request.args.get('genome_id') or '').strip()
        if not genome_id:
            return jsonify({'error': 'A genome_id query parameter is required.'}), 400

        conn = db.connect(_data_root() / 'genome.db')
        try:
            genome = _load_genome(conn, genome_id)
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404

            row = _load_view(conn, member, genome_id)
            return jsonify(_view_body(genome_id, row))
        finally:
            conn.close()

    @app.route('/api/view', methods=['POST'])
    def save_view():
        """
        Persist (upsert) the Browser_View for the signed-in member and genome.

        The JSON body carries ``genome_id`` (required), and optionally
        ``locus`` (a display string or ``{chrom,start,end}``), ``enabled`` (an
        ordered list of track ids), and ``options`` (a per-track map of
        rendering options). The genome must be viewable (else 404). The locus is
        validated and clamped against the genome's ``chrom_sizes`` — an unknown
        chromosome is rejected naming it (Requirement 8.2 behaviour); the
        enabled set is filtered to tracks that belong to the genome and are
        viewable by the member, with order preserved and duplicates dropped
        (Requirements 6.1, 6.2, 6.3); the options are stored verbatim for the
        client and renderer to apply (Requirement 8.8). The row is upserted
        keyed by ``(username, genome_id)`` and the stored view is returned.
        """
        member = identity.require_member()
        body = request.get_json(silent=True)
        if not isinstance(body, dict):
            return jsonify({'error': 'A JSON object body is required.'}), 400

        genome_id = str(body.get('genome_id') or '').strip()
        if not genome_id:
            return jsonify({'error': 'A genome_id field is required.'}), 400

        conn = db.connect(_data_root() / 'genome.db')
        try:
            genome = _load_genome(conn, genome_id)
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404

            chrom_sizes = _chrom_sizes(conn, genome_id)
            try:
                locus_str = _normalise_locus(body.get('locus'), chrom_sizes)
            except locus_math.UnknownChrom as exc:
                return jsonify({'error': str(exc), 'chrom': exc.chrom}), 400
            except locus_math.LocusParseError as exc:
                return jsonify({'error': str(exc)}), 400

            enabled = _validate_enabled(conn, member, genome_id,
                                        body.get('enabled'))
            options = _normalise_options(body.get('options'), enabled)

            _upsert_view(conn, member, genome_id, locus_str, enabled, options)
            conn.commit()

            row = _load_view(conn, member, genome_id)
            return jsonify(_view_body(genome_id, row))
        finally:
            conn.close()


# ---------------------------------------------------------------------------
# Locus / enabled / options normalisation
# ---------------------------------------------------------------------------


def _normalise_locus(raw, chrom_sizes: Dict[str, int]) -> Optional[str]:
    """
    Validate a requested locus and return its clamped display string, or None.

    Accepts either a display string (``'Chr1:1,000-5,000'`` or a bare
    ``'Chr1'``) parsed by :func:`locus.parse`, or an object
    ``{chrom, start, end}`` with 0-based half-open coordinates clamped by
    :func:`locus.clamp`. Both run through the pure locus math, so an unknown
    chromosome raises :class:`locus.UnknownChrom` (naming it) and an
    out-of-range window is clamped to the chromosome bounds (Requirement 8.3).
    An absent/empty locus returns ``None`` — a saved view need not pin a locus.
    """
    if raw is None:
        return None
    if isinstance(raw, str):
        text = raw.strip()
        if not text:
            return None
        return locus_math.parse(text, chrom_sizes).format()
    if isinstance(raw, dict):
        chrom = str(raw.get('chrom') or '').strip()
        if not chrom:
            return None
        if chrom not in chrom_sizes:
            raise locus_math.UnknownChrom(chrom)
        start = _int_or_zero(raw.get('start'))
        end = _int_or_zero(raw.get('end'), default=int(chrom_sizes[chrom]))
        clamped = locus_math.clamp(
            locus_math.Locus(chrom=chrom, start=start, end=end), chrom_sizes)
        return clamped.format()
    # Any other type is treated as "no locus" rather than an error.
    return None


def _validate_enabled(conn, member: str, genome_id: str,
                      raw) -> List[str]:
    """
    Filter the requested enabled-track ids to viewable tracks of this genome.

    The order the member sent is preserved (Requirement 6.3) and duplicates are
    dropped (first occurrence wins). A track id is kept only when it names a
    track that belongs to ``genome_id`` and is viewable by the member
    (:func:`access.can_view_track`, which also requires ``track_status='ready'``),
    so a persisted view can never enable a track that is out of reach, not part
    of the genome, or not yet ready. A non-list body yields an empty enabled set.
    """
    if not isinstance(raw, list):
        return []

    kept: List[str] = []
    seen = set()
    for item in raw:
        track_id = str(item or '').strip()
        if not track_id or track_id in seen:
            continue
        track = _load_track(conn, track_id)
        if track is None:
            continue
        if track['genome_id'] != genome_id:
            continue
        if not access.can_view_track(member, track, conn):
            continue
        kept.append(track_id)
        seen.add(track_id)
    return kept


def _normalise_options(raw, enabled: List[str]) -> Dict[str, Any]:
    """
    Keep the rendering options for enabled tracks, as a JSON-serialisable map.

    ``options`` is a per-track map (``track_id -> {colour, height, ...}``) the
    client and the SVG renderer read back to apply each track's rendering
    choices (Requirement 8.8). Only entries for currently-enabled tracks are
    kept, so options never accumulate for tracks that have been disabled or
    removed. Each entry must be an object; anything else is dropped. The values
    are stored verbatim — this module validates colours/heights only at render
    time, not here.
    """
    if not isinstance(raw, dict):
        return {}
    enabled_set = set(enabled)
    options: Dict[str, Any] = {}
    for track_id, value in raw.items():
        key = str(track_id or '').strip()
        if key and key in enabled_set and isinstance(value, dict):
            options[key] = value
    return options


def _int_or_zero(value, default: int = 0) -> int:
    """Parse an integer coordinate, falling back to ``default`` on bad input."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Database reads / the upsert
# ---------------------------------------------------------------------------


def _upsert_view(conn, username: str, genome_id: str, locus_str: Optional[str],
                 enabled: List[str], options: Dict[str, Any]) -> None:
    """
    Insert or replace the Browser_View row keyed by ``(username, genome_id)``.

    A single ``INSERT ... ON CONFLICT ... DO UPDATE`` upsert so the member's
    view for a genome is remembered and overwritten in place (Requirements 6.1,
    6.2, 6.3, 8.8). ``enabled``/``options`` are stored as JSON text in
    ``enabled_json``/``options_json`` and ``updated_at`` is refreshed.
    """
    conn.execute(
        'INSERT INTO browser_views '
        '  (username, genome_id, locus, enabled_json, options_json, updated_at) '
        'VALUES (?, ?, ?, ?, ?, datetime(\'now\')) '
        'ON CONFLICT(username, genome_id) DO UPDATE SET '
        '  locus = excluded.locus, '
        '  enabled_json = excluded.enabled_json, '
        '  options_json = excluded.options_json, '
        '  updated_at = excluded.updated_at',
        (username, genome_id, locus_str,
         json.dumps(enabled), json.dumps(options)))


def _load_view(conn, username: str, genome_id: str):
    """The browser_views row for this member + genome, or None when absent."""
    return conn.execute(
        'SELECT username, genome_id, locus, enabled_json, options_json, '
        'updated_at FROM browser_views WHERE username = ? AND genome_id = ?',
        (username, genome_id)).fetchone()


def _view_body(genome_id: str, row) -> Dict[str, Any]:
    """
    Shape a browser_views row (or its absence) into the JSON view response.

    A missing row reads as an empty view — no locus, no enabled tracks, no
    options — so the client treats "never saved" and "saved empty" the same.
    The ``enabled`` list and ``options`` map are decoded from their JSON text;
    malformed stored JSON degrades to empty rather than erroring.
    """
    if row is None:
        return {'genome_id': genome_id, 'locus': None,
                'enabled': [], 'options': {}, 'updated_at': None}
    return {
        'genome_id': genome_id,
        'locus': row['locus'],
        'enabled': _json_or(row['enabled_json'], []),
        'options': _json_or(row['options_json'], {}),
        'updated_at': row['updated_at'],
    }


def _json_or(text: Optional[str], default):
    """Decode JSON text, falling back to ``default`` on absent/invalid input."""
    if not text:
        return default
    try:
        return json.loads(text)
    except (TypeError, ValueError):
        return default


def _load_genome(conn: sqlite3.Connection, genome_id: str):
    """The reference_genomes row for this id, or None when there is none."""
    return conn.execute(
        'SELECT id, name, assembly_id, owner_lab, genome_visibility '
        'FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()


def _load_track(conn: sqlite3.Connection, track_id: str):
    """The tracks row needed to check genome membership and viewability."""
    return conn.execute(
        'SELECT id, genome_id, name, type, owner_lab, visibility, '
        'track_status, file_path FROM tracks WHERE id = ?',
        (track_id,)).fetchone()


def _chrom_sizes(conn: sqlite3.Connection, genome_id: str) -> Dict[str, int]:
    """The genome's ``chrom -> length_bp`` map, for validating the locus."""
    rows = conn.execute(
        'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?',
        (genome_id,)).fetchall()
    return {row['chrom']: int(row['length_bp']) for row in rows}


def _data_root() -> Path:
    """The configured data root (the genome-data volume), as the siblings do."""
    return Path(current_app.config['DATA_ROOT'])


# ---------------------------------------------------------------------------
# Static-asset safety + the browser page
# ---------------------------------------------------------------------------


def _safe_asset(asset: str) -> bool:
    """
    Whether a requested static path stays inside the vendored assets dir.

    A belt-and-braces guard on top of ``send_from_directory``'s own check:
    reject absolute paths, parent references, and drive-letter/UNC escapes
    before touching the filesystem, so no crafted path can read outside
    ``genome/static/``.
    """
    if not asset:
        return False
    normalised = asset.replace('\\', '/')
    if normalised.startswith('/') or '..' in normalised.split('/'):
        return False
    # Reject Windows drive letters / absolute-ish paths defensively.
    if ':' in normalised:
        return False
    return True


def _static_prefix() -> str:
    """
    The browser-visible URL prefix for vendored assets (``/genome/static``).

    The proxy strips ``/genome`` before the request reaches Flask, but the
    browser still sees it, so asset URLs the page emits must put the prefix
    back. ``Config.URL_PREFIX`` is the single source of the prefix; a bare
    ``/`` prefix (no mount) degrades to ``/static``.
    """
    prefix = str(current_app.config.get('URL_PREFIX', '/genome')).rstrip('/')
    return f'{prefix}/static' if prefix else '/static'


def _browser_html(static_prefix: str) -> str:
    """
    The self-contained igv.js browser page HTML.

    References the vendored ``{static_prefix}/igv.min.js`` (browser-visible
    ``/genome/static/igv.min.js``) rather than a CDN, so the suite has no
    outbound dependency at runtime. It mounts an IGV container and, once the
    vendored script is present, boots a browser; a small inline note explains
    that ``igv.min.js`` must be vendored into ``static/`` when the page is built
    (the static README records the same). The data the browser loads comes
    through the authenticated API endpoints under the same ``/genome`` mount.
    """
    igv_src = f'{static_prefix}/igv.min.js'
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Genome Browser</title>
  <style>
    :root {{ color-scheme: light dark; }}
    html, body {{ margin: 0; height: 100%; font-family: system-ui, sans-serif; }}
    header {{ padding: 0.75rem 1rem; border-bottom: 1px solid #ccc; }}
    header h1 {{ font-size: 1.1rem; margin: 0; }}
    #igv-container {{ min-height: 70vh; }}
    .vendor-note {{
      margin: 1rem; padding: 0.75rem 1rem; border: 1px solid #c9a227;
      background: #fff8e1; color: #5f4b00; border-radius: 6px; font-size: 0.9rem;
    }}
    .vendor-note code {{ font-family: ui-monospace, monospace; }}
    .vendor-note[hidden] {{ display: none; }}
  </style>
</head>
<body>
  <header>
    <h1>Genome Browser</h1>
  </header>

  <noscript>
    <p class="vendor-note">This browser requires JavaScript to render tracks.</p>
  </noscript>

  <div id="vendor-note" class="vendor-note" hidden>
    The IGV viewer (<code>igv.min.js</code>) has not been vendored into
    <code>static/</code> yet. Add the built <code>igv.min.js</code> there to
    enable the interactive browser; the API endpoints under
    <code>/genome/api/</code> work regardless.
  </div>

  <div id="igv-container"
       role="application"
       aria-label="Interactive genome browser"></div>

  <script src="{igv_src}"
          onerror="document.getElementById('vendor-note').hidden = false;"></script>
  <script>
    (function () {{
      var container = document.getElementById('igv-container');
      // igv.js is vendored at {igv_src}; when it is present it exposes a global
      // `igv`. Until then, show the vendor note and leave the container empty.
      if (typeof igv === 'undefined' || !igv || typeof igv.createBrowser !== 'function') {{
        document.getElementById('vendor-note').hidden = false;
        return;
      }}
      igv.createBrowser(container, {{ genome: null }}).then(function (browser) {{
        window.__genomeBrowser = browser;
      }}).catch(function (err) {{
        console.error('Failed to start the genome browser', err);
      }});
    }})();
  </script>
</body>
</html>
"""
