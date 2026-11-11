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
from html import escape
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional

from flask import (Flask, Response, abort, current_app, jsonify, request,
                   send_from_directory)

import access
import db
import identity
import locus as locus_math
import pages

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

    @app.route('/help')
    def help_page():
        """How to use the Genome Browser (``/genome/help``), like every tool's."""
        prefix = _static_prefix()
        return Response(pages.help_html(prefix, _mount(prefix)), mimetype='text/html')

    @app.route('/about')
    def about_page():
        """What the Genome Browser is (``/genome/about``), with its version."""
        from app import APP_VERSION  # the factory's module; imported here, not at load
        prefix = _static_prefix()
        return Response(pages.about_html(prefix, _mount(prefix), APP_VERSION),
                        mimetype='text/html')

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


def _mount(static_prefix: str) -> str:
    """Where the tool is mounted (``/genome``), from its asset prefix."""
    return static_prefix[:-len('/static')] if static_prefix.endswith('/static') else ''


def _browser_html(static_prefix: str) -> str:
    """
    The Genome Browser page: the landing page and the browser in one shell.

    With no query it shows the landing page, where a member chooses a genome,
    the tracks to show and where to start; ``?genome=<id>`` opens igv.js on
    that genome with the member's saved view. The behaviour is in
    ``static/browser.js`` and the layout in ``static/browser.css``, over the
    suite's look (``static/trellis/``). igv.js is the vendored
    ``{static_prefix}/igv.min.js``, never a CDN, so the suite has no outbound
    dependency at runtime; if it does not load, a note says so.
    """
    prefix = escape(_mount(static_prefix))
    static = escape(static_prefix)
    editor = escape(str(current_app.config.get('EDITOR_URL') or '').rstrip('/'))
    notebook = escape(str(current_app.config.get('NOTEBOOK_URL') or '').rstrip('/'))
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Genome Browser · Trellis Tools</title>
  <link rel="icon" href="{static}/trellis/genome.svg" type="image/svg+xml">
  <link rel="stylesheet" href="{static}/trellis/trellis.css">
  <link rel="stylesheet" href="{static}/browser.css">
</head>
<body class="tt" data-prefix="{prefix}" data-editor="{editor}" data-notebook="{notebook}">
  <header class="tt-bar">
    <a class="tt-home" href="/"><img src="{static}/trellis/trellis.svg" alt=""><span>Trellis Tools</span></a>
    <a class="tt-module" href="{prefix}/"><img src="{static}/trellis/genome.svg" alt="">Genome Browser</a>
    <span class="tt-sub" id="crumb" hidden></span>
    <nav class="tt-links">
      <a href="{prefix}/help">Help</a>
      <a href="{prefix}/about">About</a>
    </nav>
    <span class="tt-user" id="user"></span>
    <a class="tt-logout" href="/account/logout">Log out</a>
  </header>

  <noscript>
    <p class="vendor-note">The Genome Browser needs JavaScript to show genomes and tracks.</p>
  </noscript>

  <main id="choose" class="gb-page" hidden>
    <div class="gb-intro">
      <h1>Choose what to show</h1>
      <p>Pick a reference genome, the tracks to lay over it, and where to start. The browser
         remembers your view of each genome, so it opens where you left off next time.</p>
    </div>
    <p class="gb-error" id="choose-note" role="alert" hidden></p>

    <aside class="gb-guide" id="intro-guide" hidden aria-labelledby="intro-guide-title">
      <h2 id="intro-guide-title">New to the Genome Browser?</h2>
      <ol>
        <li><b>Choose a genome</b> below, tick the <b>tracks</b> to show on it, and <b>open</b> it.</li>
        <li>In the browser, <b>drag</b> to move along and use <b>+ / &minus;</b> to zoom. <b>Click</b> a gene
            for its details.</li>
        <li><b>Double-click</b> a gene, or <b>Shift-drag</b> across any stretch, to select its sequence; then
            open it in the Sequence Editor, BLAST it, or download it as FASTA.</li>
        <li>Your own data goes in with <b>+ Add tracks</b>, one file or many at once (or <b>+ Add a genome</b> for a new assembly).</li>
      </ol>
      <div class="gb-guide-actions">
        <a class="tt-btn" href="{prefix}/help">Read the help</a>
        <button type="button" class="tt-btn primary" id="intro-dismiss">Got it</button>
      </div>
    </aside>

    <form id="choose-form" novalidate>
      <section class="gb-step" aria-labelledby="step-genome-title">
        <div class="gb-step-head">
          <span class="gb-step-num" aria-hidden="true">1</span>
          <h2 id="step-genome-title">Genome</h2>
          <a class="gb-help-link" href="{prefix}/help#genome" title="About genomes and adding one">?</a>
          <span class="gb-actions">
            <button type="button" class="tt-btn" id="share-genome" hidden
                    title="See and change which labs see this genome">Share genome</button>
            <button type="button" class="tt-btn" id="add-genome">+ Add a genome</button>
          </span>
        </div>
        <div id="genome-list" class="gb-genomes" role="group" aria-labelledby="step-genome-title"></div>
      </section>

      <section class="gb-step" id="step-tracks" aria-labelledby="step-tracks-title" hidden>
        <div class="gb-step-head">
          <span class="gb-step-num" aria-hidden="true">2</span>
          <h2 id="step-tracks-title">Tracks</h2>
          <a class="gb-help-link" href="{prefix}/help#track" title="About tracks and adding one">?</a>
          <span class="gb-hint">shown in this order, top to bottom</span>
          <span class="gb-actions">
            <span id="track-actions">
              <button type="button" class="gb-link" id="tracks-all">Select all</button>
              <button type="button" class="gb-link" id="tracks-none">Select none</button>
            </span>
            <button type="button" class="tt-btn" id="add-track">+ Add tracks</button>
          </span>
        </div>
        <div class="gb-search" id="track-search-row">
          <input id="track-search" type="search" autocomplete="off" spellcheck="false"
                 placeholder="Search tracks by name or kind: rnaseq rep2, bigwig, peaks…"
                 aria-label="Search tracks" aria-controls="track-list">
          <label class="gb-sort">Sort by
            <select id="track-sort" aria-label="Sort tracks by">
              <option value="experiment">Experiment</option>
              <option value="name">Name</option>
              <option value="newest">Newest first</option>
            </select>
          </label>
          <span class="gb-hint" id="track-search-count" aria-live="polite"></span>
          <button type="button" class="gb-link" id="tracks-shown" hidden>Tick shown</button>
        </div>
        <ul id="track-list" class="gb-tracks"></ul>
        <div class="gb-elsewhere" id="track-elsewhere" hidden>
          <h3>On your other genomes</h3>
          <ul id="track-elsewhere-list"></ul>
        </div>
        <ul id="pending-list" class="gb-pending" aria-live="polite"></ul>
      </section>

      <section class="gb-step" id="step-locus" aria-labelledby="step-locus-title" hidden>
        <div class="gb-step-head">
          <span class="gb-step-num" aria-hidden="true">3</span>
          <h2 id="step-locus-title">Start at</h2>
          <span class="gb-hint">a chromosome, or a region such as Chr1:3,600-9,000</span>
        </div>
        <div class="gb-locus">
          <input id="locus" name="locus" type="text" autocomplete="off" spellcheck="false"
                 placeholder="Chr1:1-100,000" aria-labelledby="step-locus-title"
                 aria-describedby="locus-error">
        </div>
        <p id="locus-error" class="gb-error" role="alert" hidden></p>
        <div id="chrom-list" class="gb-chroms" aria-label="Chromosomes"></div>
      </section>

      <div class="gb-open" id="open-row" hidden>
        <button type="submit" class="tt-btn primary">Open in browser</button>
        <span class="gb-hint" id="open-hint"></span>
      </div>
    </form>
  </main>

  <section id="view" class="gb-view" hidden>
    <div class="gb-toolbar">
      <h1 id="view-title"></h1>
      <span class="gb-hint" id="view-hint"></span>
      <span class="gb-spacer"></span>
      <span class="gb-saved" id="saved" aria-live="polite"></span>
      <button type="button" class="tt-btn" id="controls-toggle" aria-expanded="false" aria-controls="controls">
        <span aria-hidden="true">?</span> Controls</button>
      <button type="button" class="tt-btn" id="change" title="Back to choosing the genome, tracks and where to start">Change genome or tracks</button>
    </div>
    <aside class="gb-controls" id="controls" hidden aria-labelledby="controls-title">
      <div class="gb-controls-head">
        <h2 id="controls-title">Controls</h2>
        <a href="{prefix}/help#move">Full help</a>
        <button type="button" class="gb-link" id="controls-close">Got it, hide this</button>
      </div>
      <div class="gb-controls-grid">
        <section>
          <h3>Moving around</h3>
          <dl>
            <dt>Drag the tracks</dt><dd>move along the genome</dd>
            <dt><b>+</b> / <b>&minus;</b>, slider</dt><dd>zoom in and out</dd>
            <dt>Box in the toolbar</dt><dd>go to <code>Chr1:3,600-9,000</code> or a chromosome</dd>
            <dt>Click a feature</dt><dd>its details</dd>
            <dt>Bulk track settings</dt><dd>height, colour, display or scale of several tracks at once</dd>
            <dt>⚙ on a track</dt><dd>the same, for one track</dd>
          </dl>
        </section>
        <section>
          <h3>Selecting sequence</h3>
          <dl>
            <dt>Double-click a feature</dt><dd>select exactly that gene or transcript</dd>
            <dt><kbd>Shift</kbd> + drag</dt><dd>select any stretch</dd>
            <dt>Select region</dt><dd>while on, a plain drag selects</dd>
            <dt>Type a region</dt><dd>in the bar, then <b>Select</b></dd>
            <dt><kbd>Esc</kbd>, or click outside it</dt><dd>clear the selection (Esc first leaves Select region)</dd>
          </dl>
        </section>
        <section>
          <h3>Taking it away</h3>
          <dl>
            <dt>Open in Sequence Editor</dt><dd>up to 30 kb, with its annotations</dd>
            <dt>BLAST at NCBI</dt><dd>up to 2 Mb; sends it to NCBI</dd>
            <dt>Copy sequence</dt><dd>up to 2 Mb</dd>
            <dt>Download FASTA</dt><dd>any length</dd>
          </dl>
          <p class="gb-hint">Minus-strand features are sent reverse-complemented, as they read.</p>
          <h3>Figures</h3>
          <dl>
            <dt>Save to notebook</dt><dd>the view as a figure on a page under today, or as SVG</dd>
          </dl>
        </section>
      </div>
    </aside>
    <div id="vendor-note" class="vendor-note" hidden>
      The IGV viewer (<code>igv.min.js</code>) did not load, so the interactive
      browser cannot start. It is served from <code>{static}/</code> and
      never fetched from the internet, so this means the file is missing from the
      image. The API endpoints under <code>{prefix}/api/</code> work regardless.
    </div>
    <div class="gb-selection" id="selection">
      <div class="gb-select-tools">
        <button type="button" class="tt-btn gb-bulk" id="tracks-toggle" aria-expanded="false" aria-controls="tracks-panel"
                title="Change the height, colour, display or scale of several tracks at once">
          <span class="gb-bulk-icon" aria-hidden="true"></span>Bulk track settings</button>
        <button type="button" class="tt-btn" id="figure-open"
                title="Save this view as a figure: to your Lab Notebook as a page under today, or as an SVG file">Save to notebook</button>
        <button type="button" class="tt-btn" id="analyse-open"
                title="Profiles and heatmaps around features, differences between tracks, replicate checks, tables, bar charts and dot plots">Analyse</button>
        <button type="button" class="tt-btn" id="circos-open"
                title="Draw the genome round a circle with the tracks shown as rings: a Circos plot, as SVG, PNG or to your Lab Notebook">Circos plot</button>
        <span class="gb-tool-sep" aria-hidden="true"></span>
        <button type="button" class="tt-btn gb-mode" id="select-mode" aria-pressed="false"
                title="While on, dragging across the tracks selects sequence instead of moving along">
          <span class="gb-mode-icon" aria-hidden="true"></span>Select region</button>
        <form id="select-form" class="gb-select-form" novalidate>
          <input id="select-locus" type="text" autocomplete="off" spellcheck="false"
                 placeholder="or type Chr1:3,631-5,899" aria-label="Select a region by its coordinates">
          <button type="submit" class="tt-btn">Select</button>
        </form>
      </div>
      <div class="gb-selection-text" id="selection-text" aria-live="polite"></div>
      <div class="gb-selection-actions" id="selection-actions" hidden>
        <button type="button" class="tt-btn primary" id="send-editor"
                title="Open it as a new sequence in the Sequence Editor, with its annotations (up to 30 kb)">Open in Sequence Editor</button>
        <button type="button" class="tt-btn" id="send-blast"
                title="Open NCBI nucleotide BLAST with it filled in (up to 2 Mb). This sends the sequence to NCBI.">BLAST at NCBI</button>
        <button type="button" class="tt-btn" id="copy-seq"
                title="Copy the bases to the clipboard (up to 2 Mb)">Copy sequence</button>
        <button type="button" class="tt-btn" id="save-fasta"
                title="Download it as a FASTA file, at any length">Download FASTA</button>
        <button type="button" class="gb-link" id="clear-selection" title="Clear the selection (Esc)">Clear</button>
      </div>
      <p class="gb-note" id="selection-note" hidden></p>
      <p class="gb-error" id="selection-error" role="alert" hidden></p>
    </div>
    <div id="drag-band" class="gb-drag-band" hidden></div>
    <aside class="gb-tracks-panel" id="tracks-panel" hidden aria-labelledby="tracks-panel-title">
      <div class="gb-controls-head">
        <h2 id="tracks-panel-title">Bulk track settings</h2>
        <span class="gb-hint">Tick the tracks to change, then set them all at once. Kept with your view.</span>
        <button type="button" class="gb-link" id="tracks-panel-close">Done</button>
      </div>
      <div class="gb-tp-body">
        <div class="gb-tp-list">
          <input id="tp-search" class="gb-tp-search" type="search" autocomplete="off" spellcheck="false"
                 placeholder="Find a track…" aria-label="Find a track in this view" aria-controls="tp-tracks">
          <div class="gb-tp-list-head">
            <button type="button" class="gb-link" id="tp-all">Select all</button>
            <button type="button" class="gb-link" id="tp-none">Select none</button>
            <span class="gb-hint" id="tp-count"></span>
          </div>
          <ul id="tp-tracks"></ul>
        </div>
        <div class="gb-tp-controls" id="tp-controls">
          <div class="gb-tp-row">
            <span class="gb-tp-label">Height</span>
            <span class="gb-tp-chips">
              <button type="button" class="gb-chip" data-height="40">Short</button>
              <button type="button" class="gb-chip" data-height="80">Medium</button>
              <button type="button" class="gb-chip" data-height="160">Tall</button>
            </span>
            <form id="tp-height-form" class="gb-tp-inline">
              <input id="tp-height" type="number" min="15" max="1000" step="5" placeholder="px" aria-label="Height in pixels">
              <button type="submit" class="tt-btn">Set</button>
            </form>
          </div>
          <div class="gb-tp-row">
            <span class="gb-tp-label">Colour</span>
            <span class="gb-tp-swatches" id="tp-swatches"></span>
            <label class="gb-tp-inline" title="Pick any colour">
              <input id="tp-color" type="color" value="#2b8040" aria-label="Pick any colour">
            </label>
            <form id="tp-hex-form" class="gb-tp-inline" title="A colour as a hex code, such as #1f5fa8">
              <input id="tp-hex" type="text" maxlength="7" spellcheck="false" autocomplete="off"
                     placeholder="#1f5fa8" aria-label="Colour as a hex code" class="gb-tp-hex">
              <button type="submit" class="tt-btn">Set</button>
            </form>
            <button type="button" class="gb-link" id="tp-color-reset">Original colours</button>
          </div>
          <div class="gb-tp-row" id="tp-display-row">
            <span class="gb-tp-label">Genes</span>
            <span class="gb-tp-chips">
              <button type="button" class="gb-chip" data-display="EXPANDED" title="Every transcript on its own row">Expanded</button>
              <button type="button" class="gb-chip" data-display="SQUISHED" title="Every transcript, in thin rows">Squished</button>
              <button type="button" class="gb-chip" data-display="COLLAPSED" title="All on one row">Collapsed</button>
            </span>
          </div>
          <div class="gb-tp-row" id="tp-scale-row">
            <span class="gb-tp-label">Scale</span>
            <span class="gb-tp-chips">
              <button type="button" class="gb-chip" data-scale="auto" title="Each track fits its own data on screen">Each its own</button>
              <button type="button" class="gb-chip" data-scale="group" title="The ticked tracks share one scale, so their heights compare">Shared</button>
            </span>
            <form id="tp-range-form" class="gb-tp-inline" title="A fixed range, the same wherever you look">
              <input id="tp-min" type="number" step="any" placeholder="min" aria-label="Minimum">
              <input id="tp-max" type="number" step="any" placeholder="max" aria-label="Maximum">
              <button type="submit" class="tt-btn">Fix range</button>
            </form>
          </div>
          <div class="gb-tp-row">
            <span class="gb-tp-label"></span>
            <button type="button" class="tt-btn danger" id="tp-remove">Remove from view</button>
          </div>
          <p class="gb-error" id="tp-error" role="alert" hidden></p>
        </div>
      </div>
    </aside>
    <div id="igv-container" role="application" aria-label="Interactive genome browser"></div>
  </section>

  <dialog id="genome-dialog" class="gb-dialog" aria-labelledby="genome-dialog-title">
    <form method="dialog" id="genome-form" novalidate>
      <h2 id="genome-dialog-title">Add a genome</h2>
      <p class="gb-hint">A reference sequence as FASTA (.fa, .fasta, .fna, .fsa, gzipped or not), up
        to 4 GB. <span id="genome-owner"></span>
        <a href="{prefix}/help#genome" target="_blank" rel="noopener">How adding a genome works</a></p>
      <label class="gb-field">Name
        <input name="name" required maxlength="200" placeholder="Arabidopsis thaliana Col-0 (TAIR10)">
      </label>
      <label class="gb-field">Assembly
        <input name="assembly_id" required maxlength="100" placeholder="TAIR10">
        <span class="gb-hint">the assembly's identifier; tracks name its chromosomes</span>
      </label>
      <label class="gb-field">FASTA file
        <input name="file" type="file" required accept=".fa,.fasta,.fna,.fsa,.gz">
      </label>
      <label class="gb-check" id="suite-reference-row" hidden>
        <input name="suite_reference" type="checkbox" value="true">
        <span>Share with everyone in the suite (a suite reference genome, owned by no lab)</span>
      </label>
      <div class="gb-progress" hidden><progress max="100" value="0"></progress><span></span></div>
      <p class="gb-error" role="alert" hidden></p>
      <div class="gb-dialog-actions">
        <button type="button" class="tt-btn" value="cancel" data-close>Cancel</button>
        <button type="submit" class="tt-btn primary">Upload genome</button>
      </div>
    </form>
  </dialog>

  <dialog id="share-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="share-title">
    <form method="dialog" id="share-form" novalidate>
      <h2 id="share-title">Share</h2>
      <p class="gb-hint" id="share-intro"></p>
      <h3 class="gb-share-h">Who sees it</h3>
      <ul class="gb-share-list" id="share-list"></ul>
      <div id="share-add">
        <h3 class="gb-share-h">Share with</h3>
        <div class="gb-share-add-row">
          <select id="share-target" aria-label="A lab or chair group to share with"></select>
          <button type="button" class="tt-btn primary" id="share-go">Share</button>
        </div>
        <p class="gb-hint" id="share-note"></p>
      </div>
      <div class="gb-notice" id="share-done" hidden></div>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <button type="button" class="tt-btn" id="share-close">Close</button>
      </div>
    </form>
  </dialog>

  <dialog id="detail-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="detail-title">
    <form method="dialog" id="detail-form" novalidate>
      <h2 id="detail-title">Track</h2>
      <p class="gb-hint" id="detail-intro"></p>
      <h3 class="gb-share-h">Experiment</h3>
      <div class="gb-share-add-row">
        <input id="detail-experiment" list="experiment-names" maxlength="120" autocomplete="off"
               aria-label="Experiment" placeholder="No experiment">
        <button type="button" class="tt-btn primary" id="detail-experiment-save">Save</button>
      </div>
      <h3 class="gb-share-h">Plants, lines and seed lots</h3>
      <ul class="gb-share-list" id="detail-tags"></ul>
      <div id="detail-tag-add">
        <div class="gb-share-add-row">
          <input id="detail-stock-search" type="search" autocomplete="off"
                 aria-label="Find a line, seed lot or plant in the Seed Bank"
                 placeholder="Find a line, seed lot or plant: code or line name">
        </div>
        <ul class="gb-share-list" id="detail-stock-results"></ul>
      </div>
      <h3 class="gb-share-h">In the Lab Notebook</h3>
      <ul class="gb-share-list" id="detail-mentions"></ul>
      <p class="gb-hint">Type <b>@</b> in a notebook page and choose <i>A genome track</i> to link it there.</p>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <button type="button" class="tt-btn" id="detail-close">Close</button>
      </div>
    </form>
  </dialog>

  <dialog id="analyse-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="analyse-dialog-title">
    <form method="dialog" id="analyse-form" novalidate>
      <h2 id="analyse-dialog-title">Analyse</h2>
      <label class="gb-field">What to draw
        <select id="analyse-kind" name="kind"></select>
      </label>
      <p class="gb-hint" id="analyse-about"></p>
      <label class="gb-field">Title <span class="gb-hint">(optional)</span>
        <input name="title" maxlength="100">
      </label>
      <div data-for="tracks">
        <h3 class="gb-share-h">Coverage tracks <span class="gb-hint">(those shown are ticked)</span></h3>
        <div class="gb-analyse-checks" id="analyse-tracks"></div>
      </div>
      <div data-for="features">
        <h3 class="gb-share-h">Features</h3>
        <label class="gb-field">From the track
          <select id="analyse-feature-track"></select>
        </label>
        <p class="gb-hint" id="analyse-type-hint"></p>
        <div class="gb-analyse-checks" id="analyse-types"></div>
        <label class="gb-field">Only these <span class="gb-hint">(optional: names as in the track, separated by commas;
          for example the genes of a list)</span>
          <textarea name="names" rows="2" maxlength="100000" spellcheck="false" placeholder="AT1G01010, AT2G18790"></textarea>
        </label>
      </div>
      <div data-for="profile" class="gb-circos-row">
        <label class="gb-field">Around
          <select name="mode">
            <option value="scale" selected>The whole feature, scaled to one length</option>
            <option value="start">Its start (transcription start, for genes)</option>
            <option value="end">Its end</option>
          </select>
        </label>
        <label class="gb-field">With, either side
          <select name="flank">
            <option value="1000">1 kb</option><option value="2000" selected>2 kb</option>
            <option value="5000">5 kb</option><option value="10000">10 kb</option>
          </select>
        </label>
        <label class="gb-field">Scaled length
          <select name="body">
            <option value="1000">1 kb</option><option value="3000" selected>3 kb</option><option value="5000">5 kb</option>
          </select>
        </label>
        <label class="gb-field">Step
          <select name="bin">
            <option value="50">50 bp</option><option value="100" selected>100 bp</option><option value="200">200 bp (faster)</option>
          </select>
        </label>
      </div>
      <div data-for="profile" class="gb-circos-row">
        <label class="gb-check"><input type="checkbox" name="heatmap" checked> Heatmap too: every feature a row, sorted by signal</label>
        <label class="gb-field">Group the rows
          <select name="kmeans">
            <option value="0" selected>No: sorted only</option><option value="2">Into 2 clusters</option>
            <option value="3">Into 3 clusters</option><option value="4">Into 4 clusters</option>
          </select>
        </label>
        <label class="gb-field">Colours
          <select name="colormap">
            <option value="viridis" selected>Viridis</option><option value="magma">Magma</option>
            <option value="Reds">Reds</option><option value="Blues">Blues</option><option value="RdBu_r">Blue to red</option>
          </select>
        </label>
      </div>
      <div data-for="compare">
        <div class="gb-circos-row">
          <label class="gb-field">Track A <select id="analyse-a"></select></label>
          <label class="gb-field">Track B <select id="analyse-b"></select></label>
        </div>
        <div class="gb-circos-row">
          <label class="gb-field">Make
            <select name="operation">
              <option value="subtract" selected>A minus B</option>
              <option value="log2">log2(A / B)</option>
              <option value="ratio">A / B</option>
              <option value="mean">The mean of A and B</option>
            </select>
          </label>
          <label class="gb-field">Step
            <select name="compare_bin">
              <option value="50">50 bp</option><option value="100" selected>100 bp</option><option value="1000">1 kb</option>
            </select>
          </label>
          <label class="gb-field">Pseudocount <span class="gb-hint">(ratios: added to both, so 0 divides)</span>
            <input name="pseudocount" type="number" min="0.0001" step="any" value="1">
          </label>
        </div>
        <label class="gb-field">Name of the new track <span class="gb-hint">(optional)</span>
          <input name="track_name" maxlength="200" placeholder="drm1/2 minus WT, CHH">
        </label>
      </div>
      <div data-for="correlation" class="gb-circos-row">
        <label class="gb-field">Correlation
          <select name="method">
            <option value="pearson" selected>Pearson</option><option value="spearman">Spearman (ranks)</option>
          </select>
        </label>
        <label class="gb-field">In bins of
          <select name="corr_bin">
            <option value="5000">5 kb</option><option value="10000" selected>10 kb</option>
            <option value="50000">50 kb</option><option value="100000">100 kb</option>
          </select>
        </label>
      </div>
      <p data-for="selection" class="gb-hint">Selection: <b id="analyse-selection-text"></b></p>
      <label data-for="links" class="gb-field">Links track
        <select id="analyse-links-track"></select>
      </label>
      <div class="gb-notice gb-analyse-status" id="analyse-status" aria-live="polite"></div>
      <div id="analyse-result"></div>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <span class="gb-spacer"></span>
        <button type="button" class="tt-btn" id="analyse-close">Close</button>
        <button type="submit" class="tt-btn primary" id="analyse-run">Draw</button>
      </div>
    </form>
  </dialog>

  <dialog id="circos-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="circos-dialog-title">
    <form method="dialog" id="circos-form" novalidate>
      <h2 id="circos-dialog-title">Circos plot</h2>
      <p class="gb-hint">The genome drawn round a circle with the tracks shown, made with Circos (circos.ca). Show,
        hide, order and colour the tracks in the browser first; see Help for how each is drawn.</p>
      <label class="gb-field">Title <span class="gb-hint">(optional)</span>
        <input name="title" maxlength="120" placeholder="The genome and what is drawn">
      </label>
      <h3 class="gb-share-h">Tracks, outermost ring first</h3>
      <ul class="gb-circos-tracks" id="circos-tracks"></ul>
      <fieldset class="gb-field">
        <legend>Draw</legend>
        <label><input type="radio" name="where" value="genome" checked> The whole genome</label>
        <label class="gb-check gb-circos-indent"><input type="checkbox" name="skip_small" checked> Leave out small
          chromosomes (organelles, scaffolds: under 1% of the genome)</label>
        <label><input type="radio" name="where" value="chroms"> These chromosomes:</label>
        <div class="gb-circos-chroms gb-circos-indent" id="circos-chroms"></div>
        <label><input type="radio" name="where" value="region"> A region:
          <input name="region" class="gb-circos-region" placeholder="Chr1:1,000,000-5,000,000" spellcheck="false"></label>
      </fieldset>
      <div class="gb-circos-row">
        <label class="gb-field">Detail
          <select name="resolution">
            <option value="coarse">Coarse (about 300 steps)</option>
            <option value="medium" selected>Medium (about 800)</option>
            <option value="fine">Fine (about 2,000: slower)</option>
          </select>
        </label>
        <label class="gb-field">or steps of exactly <span class="gb-hint">(bp)</span>
          <input name="step" type="number" min="100" step="100" placeholder="e.g. 100000">
        </label>
      </div>
      <fieldset class="gb-field">
        <legend>Scale</legend>
        <label class="gb-check"><input type="checkbox" name="shared_scale"> Coverage rings share one scale, so their
          heights compare</label>
        <label class="gb-check"><input type="checkbox" name="log_scale"> Coverage on a log scale (for a few high peaks)</label>
      </fieldset>
      <label class="gb-check" id="circos-highlight-row"><input type="checkbox" name="highlight">
        Highlight the selection, <span id="circos-highlight-text"></span></label>
      <label class="gb-field">Label genes <span class="gb-hint">(optional: names as in a gene annotation track shown,
        separated by commas, up to 100)</span>
        <textarea name="labels" rows="2" maxlength="4000" placeholder="AT1G01010, AT3G24650" spellcheck="false"></textarea>
      </label>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <span class="gb-spacer"></span>
        <button type="button" class="tt-btn" id="circos-close">Cancel</button>
        <button type="submit" class="tt-btn primary" id="circos-draw">Draw</button>
      </div>
    </form>
  </dialog>

  <dialog id="figure-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="figure-dialog-title">
    <form method="dialog" id="figure-form" novalidate>
      <h2 id="figure-dialog-title">Save this view as a figure</h2>
      <p class="gb-hint" id="figure-what">Exactly what is on screen: the tracks, the position and any selection.</p>
      <div class="gb-figure-preview"><img id="figure-preview" alt="The figure as it will be saved"></div>
      <label class="gb-field">Page title
        <input name="title" maxlength="200">
      </label>
      <label class="gb-field">Note <span class="gb-hint">(optional: what it shows, which experiment)</span>
        <textarea name="note" rows="3" maxlength="4000"></textarea>
      </label>
      <p class="gb-hint" id="figure-where"></p>
      <div class="gb-notice" id="figure-done" hidden></div>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <button type="button" class="tt-btn" id="figure-svg">Download SVG</button>
        <span id="figure-downloads" class="gb-figure-downloads"></span>
        <button type="button" class="tt-btn" id="figure-circos-files" hidden
                title="Circos's own configuration and data files, to change the plot and run Circos again">Download Circos files</button>
        <span class="gb-spacer"></span>
        <button type="button" class="tt-btn" id="figure-close">Close</button>
        <button type="submit" class="tt-btn primary" id="figure-save">Save to Lab Notebook</button>
      </div>
    </form>
  </dialog>

  <dialog id="track-dialog" class="gb-dialog gb-dialog-wide" aria-labelledby="track-dialog-title">
    <form method="dialog" id="track-form" novalidate>
      <h2 id="track-dialog-title">Add tracks</h2>
      <p class="gb-hint">To <strong id="track-genome-name"></strong>. <span id="track-owner"></span></p>
      <div class="gb-drop" id="track-drop">
        <input type="file" id="track-files" multiple hidden
               accept=".gff,.gff3,.bed,.bedpe,.bw,.bigwig,.bam,.wig,.wiggle,.bedgraph,.bg,.gz">
        <p><b>Drop files here</b>, or <button type="button" class="gb-link" id="track-browse">choose files</button>:
          one or as many as you like.</p>
        <p class="gb-hint">GFF or BED (features), BigWig (coverage), or BAM, wig and bedGraph, which are
          converted to coverage on the server; gzipped is fine. Each file's chromosome names are checked
          against the genome's, and ones written another way (1, chr1, Chr1) are renamed to match.
          <a href="{prefix}/help#track" target="_blank" rel="noopener">More on tracks</a></p>
      </div>
      <label class="gb-field" id="track-experiment-row">Experiment <span class="gb-hint">(optional, for every file
        added now: the track list groups tracks by it)</span>
        <input id="track-experiment" list="experiment-names" maxlength="120" autocomplete="off"
               placeholder="Heat stress RNA-seq 2026">
      </label>
      <datalist id="experiment-names"></datalist>
      <ul class="gb-uploads" id="track-rows" aria-live="polite"></ul>
      <p class="gb-upload-summary" id="track-summary" role="status" hidden></p>
      <div class="gb-error" role="alert" hidden></div>
      <div class="gb-dialog-actions">
        <button type="button" class="tt-btn" id="track-cancel">Cancel</button>
        <button type="submit" class="tt-btn primary" id="track-submit" disabled>Upload</button>
      </div>
    </form>
  </dialog>

  <script src="{static}/igv.min.js"
          onerror="document.getElementById('vendor-note').hidden = false;"></script>
  <script src="{static}/trellis/modules.js" defer></script>
  <script src="{static}/browser.js"></script>
</body>
</html>
"""
