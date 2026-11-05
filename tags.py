"""
Plasmid and seed-stock tagging: the endpoints that link a Track or
Reference_Genome to a plasmid in the MoClo/Editor tools or a seed stock in the
Seed Bank, and that feed the Seed Bank's Family_Tree (Requirement 11).

This module owns the tagging side of the Genome Browser and is wired by
:func:`register_tags`, called once from the app factory exactly like
``register_genomes``/``register_tracks``/``register_figures`` — so the factory
stays a thin wiring layer and this module never has to edit the factory body or
its siblings. (The one app.py line this needs is noted at the bottom of this
docstring.)

Two tag kinds, two external tools, one local store
--------------------------------------------------
A tag is a reference from a Track or a Reference_Genome *out* to a record in
another tool. The Genome Browser keeps nothing of the other tool's data; it
only records the external id in its own database and validates that the id
exists, reading the other tool's store **read-only**:

- **Plasmid_Tag** — a reference to a plasmid in the MoClo Library Tool / Plasmid
  Editor. The id must exist in those tools (``plasmid_exists``). The row lands
  in ``plasmid_tags``, keyed by ``plasmid_id``.
- **Seed_Stock_Tag** — a reference to a seed stock (a line, a seed lot, or a
  plant) in the Seed Bank. The id must exist there, and its *kind*
  (``line`` | ``lot`` | ``plant``) is resolved at the same time
  (``seed_stock_lookup``). The row lands in ``seed_stock_tags``, keyed by
  ``stock_id`` — and for a **line or plant** that ``stock_id`` *is* the
  Family_Tree key: the Seed Bank reads these rows read-only by ``stock_id`` and
  surfaces the genome-browser link in that line's or plant's Family_Tree
  (Requirement 11.5; see ``CONTRACT_seed_stock_tags.md``, task 16.4).

Endpoints (all relative to the stripped ``/genome`` mount)
----------------------------------------------------------
- ``POST /api/tracks/<id>/tags`` and ``POST /api/genomes/<id>/tags`` — add a
  Plasmid_Tag or Seed_Stock_Tag. The body says which kind and the external id.
  The member must be signed in and able to *view* the target (not necessarily
  own it — the requirements do not restrict tagging to the owner lab; a member
  viewing a track/genome may annotate it). A missing external id is rejected
  naming it (Requirement 11.3), a duplicate on the same target is rejected while
  the existing reference is retained (Requirement 11.4), and a valid id writes
  one row keyed by the external id (Requirements 11.1, 11.2, 11.5).
- ``GET /api/tracks/<id>/tags`` and ``GET /api/genomes/<id>/tags`` — list the
  recorded tags, each with a navigable link to the corresponding MoClo/Editor
  record or Seed Bank record (``MOCLO_URL`` / ``EDITOR_URL`` / ``SEEDS_URL``)
  (Requirement 11.6).
- ``DELETE /api/tracks/<id>/tags/<tag_id>`` and
  ``DELETE /api/genomes/<id>/tags/<tag_id>`` — remove one tag row. For a seed
  line/plant tag this removes the Family_Tree link too, because the Seed Bank
  reads these rows (Requirements 11.7, 11.8).

Validation is isolated behind two small functions
--------------------------------------------------
The MoClo, Editor and Seed Bank databases ship as prebuilt images, so their
exact table/column names are not known from this repo. The existence checks are
therefore isolated behind :func:`plasmid_exists` and :func:`seed_stock_lookup`,
each of which tries a small set of documented candidate table/column names over
a read-only connection and tolerates a missing store or schema. This keeps the
endpoint logic independent of the external schema, makes the validators trivial
to stub in tests, and means a schema correction later touches one function.
See the ASSUMPTIONS block above each validator for the candidate names and the
read-only DB paths used.

App-factory wiring (added by the orchestrator, not here)::

    from tags import register_tags
    ...
    register_tags(app)

The 401 handler for an absent member is registered once in ``genomes.py`` for
the whole module; this module relies on it rather than re-registering it.
"""

import os
import sqlite3
import uuid
from pathlib import Path
from typing import Optional, Tuple

from flask import Flask, current_app, jsonify, request

import access
import db
import identity


# ---------------------------------------------------------------------------
# Read-only openers for the other tools' stores
# ---------------------------------------------------------------------------
#
# db.py ships read-only openers for suite.db and moclo.db; tagging additionally
# needs the Plasmid Editor's database and the Seed Bank's database. These reuse
# db.open_readonly (mode=ro, degrades to None when absent) so this module can
# never write to another tool's store, matching the module-wide convention.
#
# Paths come from the environment, defaulting to the read-only mounts the design
# gives the genome service (task 16.1): editor-data at /editor, seeds-data at
# /seeds, moclo-data at /moclo. An operator can override any of them.

def _moclo_db_path() -> str:
    """The MoClo database path (moclo-data mount), read-only."""
    return os.environ.get('MOCLO_DB', '/moclo/moclo.db')


def _parts_db_path() -> str:
    """The MoClo parts database path (moclo-data mount), read-only."""
    return os.environ.get('PARTS_DB', '/moclo/parts.db')


def _editor_db_path() -> str:
    """
    The Plasmid Editor database path, read-only.

    EDITOR_DB overrides it; the default is the editor-data mount the design
    gives the genome service. The editor stores its DB as ``plasmid_editor.db``
    (see the editor service's ``EDITOR_DATABASE_PATH`` in docker-compose.yml).
    """
    return os.environ.get('EDITOR_DB', '/editor/plasmid_editor.db')


def _seeds_db_path() -> str:
    """
    The Seed Bank database path, read-only.

    SEEDS_DB overrides it; the default is the seeds-data mount the design gives
    the genome service. The Seed Bank keeps its data under ``/seeds``; the
    database file name is an assumption documented on ``seed_stock_lookup``.
    """
    return os.environ.get('SEEDS_DB', '/seeds/seeds.db')


def _open_ro(path: str) -> Optional[sqlite3.Connection]:
    """A read-only connection to a shared store, or None when it is absent."""
    return db.open_readonly(path)


# ---------------------------------------------------------------------------
# Plasmid_Tag validation (Requirement 11.3 — the MoClo/Editor side)
# ---------------------------------------------------------------------------

def plasmid_exists(plasmid_id: str) -> bool:
    """
    Whether a plasmid with this id exists in the MoClo tools (read-only).

    Isolated behind one function so the endpoint logic is independent of the
    external schema and tests can stub it. A truthy result means the id is a
    real MoClo/Editor plasmid; a false result (including an absent/unreadable
    store) means the tag must be rejected naming the missing id (Requirement
    11.3).

    ASSUMPTIONS (the MoClo/Editor DBs ship as prebuilt images, so these are
    candidate table/column names tried in order until one resolves; a schema
    correction later edits only this function):

    - MoClo ``moclo.db``: a ``plasmids`` table keyed by ``id`` or ``name``, or a
      ``constructs``/``assemblies`` table keyed by ``id``/``name``.
    - MoClo ``parts.db``: a ``parts`` table keyed by ``id``/``name`` (a plasmid
      id may be a MoClo part/construct id).
    - Editor ``plasmid_editor.db``: a ``plasmids`` table keyed by ``id``/``name``.

    The id is matched either against the primary id column or a human-readable
    name column, because a member tags by whichever identifier the other tool
    shows them. Any lookup error on a given candidate is swallowed and the next
    candidate tried, so an unexpected schema never raises into the endpoint.
    """
    plasmid_id = (plasmid_id or '').strip()
    if not plasmid_id:
        return False

    candidates = [
        (_moclo_db_path(), ['plasmids', 'constructs', 'assemblies']),
        (_parts_db_path(), ['parts', 'plasmids']),
        (_editor_db_path(), ['plasmids', 'constructs']),
    ]
    for path, tables in candidates:
        conn = _open_ro(path)
        if conn is None:
            continue
        try:
            for table in tables:
                if _row_exists(conn, table, plasmid_id):
                    return True
        finally:
            conn.close()
    return False


# ---------------------------------------------------------------------------
# Seed_Stock_Tag validation (Requirement 11.3 — the Seed Bank side)
# ---------------------------------------------------------------------------

#: The seed-stock kinds the module records (Requirement 11.5: line/plant are the
#: Family_Tree keys; a lot is recorded too but is not a Family_Tree node).
STOCK_KINDS = ('line', 'lot', 'plant')


def seed_stock_lookup(stock_id: str) -> Tuple[bool, Optional[str]]:
    """
    Whether a seed stock with this id exists in the Seed Bank, and its kind.

    Returns ``(exists, kind)`` where ``kind`` is one of ``'line'``, ``'lot'``
    or ``'plant'`` when found, else ``(False, None)``. Isolated behind one
    function so the endpoint is independent of the Seed Bank schema and tests
    can stub it. A false result (including an absent/unreadable store) means the
    tag is rejected naming the missing id (Requirement 11.3). The ``kind`` is
    what decides, for a ``line``/``plant``, that the row is a Family_Tree key
    the Seed Bank will read (Requirement 11.5).

    ASSUMPTIONS (the Seed Bank DB ships as a prebuilt image, so these are
    candidate table names tried in order; a schema correction later edits only
    this function):

    - Read-only ``/seeds/seeds.db`` (``SEEDS_DB`` overrides the path).
    - A ``lines`` table for plant lines, a ``seed_lots`` (or ``lots``) table for
      seed lots, and a ``plants`` table for plants — each keyed by ``id`` or a
      human-readable ``name``/``code`` column. The first table in which the id
      is found decides the kind, checked in the order lines → lot → plant.

    Because kind comes from *which* table holds the id, the function checks each
    kind's table(s) in turn and returns the first match; any per-candidate
    lookup error is swallowed so an unexpected schema never raises.
    """
    stock_id = (stock_id or '').strip()
    if not stock_id:
        return False, None

    conn = _open_ro(_seeds_db_path())
    if conn is None:
        return False, None
    try:
        kind_tables = [
            ('line', ['lines', 'plant_lines']),
            ('lot', ['seed_lots', 'lots', 'seedlots']),
            ('plant', ['plants']),
        ]
        for kind, tables in kind_tables:
            for table in tables:
                if _row_exists(conn, table, stock_id):
                    return True, kind
    finally:
        conn.close()
    return False, None


# ---------------------------------------------------------------------------
# A tolerant existence check over an unknown-but-guessed schema
# ---------------------------------------------------------------------------

#: The id/name columns a record might be keyed by, tried in order. A member
#: tags by whichever identifier the other tool shows them, so both the primary
#: id and a human-readable name/code are accepted.
_ID_COLUMNS = ('id', 'name', 'code', 'identifier', 'uid')


def _row_exists(conn: sqlite3.Connection, table: str, value: str) -> bool:
    """
    Whether ``table`` has a row whose id/name column equals ``value``.

    Tolerant by design: a missing table or a missing candidate column is not an
    error (the store's real schema may differ from the guessed names), it just
    means "not found here" so the caller tries the next candidate. Only columns
    the table actually has are queried, matched case-sensitively against the
    external id exactly as the member supplied it.
    """
    present = _table_columns(conn, table)
    if not present:
        return False
    for column in _ID_COLUMNS:
        if column not in present:
            continue
        try:
            row = conn.execute(
                f'SELECT 1 FROM "{table}" WHERE "{column}" = ? LIMIT 1',
                (value,)).fetchone()
        except sqlite3.Error:
            continue
        if row is not None:
            return True
    return False


def _table_columns(conn: sqlite3.Connection, table: str) -> set:
    """The column names of ``table``, or an empty set when it does not exist."""
    try:
        rows = conn.execute(f'PRAGMA table_info("{table}")').fetchall()
    except sqlite3.Error:
        return set()
    return {r[1] for r in rows}


# ---------------------------------------------------------------------------
# Navigable links to the external records (Requirement 11.6)
# ---------------------------------------------------------------------------

def _plasmid_link(plasmid_id: str) -> str:
    """
    A navigable link to a plasmid's MoClo/Editor record (Requirement 11.6).

    Built from ``MOCLO_URL`` (the MoClo tool's mount). The link targets the
    MoClo library, which is where a plasmid is browsed in the suite; the Plasmid
    Editor (``EDITOR_URL``) opens from there. The id is URL-quoted so an id with
    spaces or slashes still forms a valid link.
    """
    base = str(current_app.config.get('MOCLO_URL') or '/')
    return _join(base, 'plasmids', plasmid_id)


def _seed_stock_link(stock_kind: str, stock_id: str) -> str:
    """
    A navigable link to a seed stock's Seed Bank record (Requirement 11.6).

    Built from ``SEEDS_URL``; the path segment reflects the kind (``lines`` /
    ``lots`` / ``plants``) so the link lands on the right Seed Bank record.
    """
    base = str(current_app.config.get('SEEDS_URL') or '/seeds')
    segment = {'line': 'lines', 'lot': 'lots', 'plant': 'plants'}.get(
        stock_kind, 'stocks')
    return _join(base, segment, stock_id)


def _join(base: str, *parts: str) -> str:
    """Join a base URL and path parts into a single ``/``-separated link."""
    from urllib.parse import quote
    pieces = [base.rstrip('/')]
    pieces.extend(quote(str(p), safe='') for p in parts if p is not None)
    return '/'.join(pieces)


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

def register_tags(app: Flask) -> None:
    """
    Wire the tagging endpoints onto the application (called by the app factory).

    One registrar, mirroring ``register_genomes``/``register_tracks``/
    ``register_figures``, so the factory wires this with a single line and this
    module owns its own routes.
    """

    # POST — add a tag to a track or a genome -------------------------------
    @app.route('/api/tracks/<target_id>/tags', methods=['POST'])
    def add_track_tag(target_id):
        return _add_tag('track', target_id)

    @app.route('/api/genomes/<target_id>/tags', methods=['POST'])
    def add_genome_tag(target_id):
        return _add_tag('genome', target_id)

    # GET — list a target's tags with navigable links -----------------------
    @app.route('/api/tracks/<target_id>/tags')
    def list_track_tags(target_id):
        return _list_tags('track', target_id)

    @app.route('/api/genomes/<target_id>/tags')
    def list_genome_tags(target_id):
        return _list_tags('genome', target_id)

    # DELETE — remove one tag row -------------------------------------------
    @app.route('/api/tracks/<target_id>/tags/<tag_id>', methods=['DELETE'])
    def delete_track_tag(target_id, tag_id):
        return _delete_tag('track', target_id, tag_id)

    @app.route('/api/genomes/<target_id>/tags/<tag_id>', methods=['DELETE'])
    def delete_genome_tag(target_id, tag_id):
        return _delete_tag('genome', target_id, tag_id)


# ---------------------------------------------------------------------------
# Endpoint bodies (shared by the track and genome variants)
# ---------------------------------------------------------------------------

def _add_tag(target_type: str, target_id: str):
    """
    Add a Plasmid_Tag or Seed_Stock_Tag to a track or genome.

    The member must be signed in and able to view the target (tagging is a
    viewer action, not an owner action — the requirements do not restrict it to
    the Owner_Lab). The body names the kind (``plasmid`` or ``seed_stock``) and
    the external id. A missing external id in the referenced tool is rejected
    naming it (11.3); a duplicate on the same target is rejected while the
    existing row is retained (11.4); a valid id writes exactly one row keyed by
    the external id and returns it (11.1, 11.2, 11.5).
    """
    member = identity.require_member()

    body = request.get_json(silent=True) or {}
    kind = (body.get('kind') or body.get('tag_type') or '').strip().lower()
    external_id = _external_id(body)

    if kind not in ('plasmid', 'seed_stock', 'seed-stock', 'seedstock'):
        return jsonify({
            'error': "A tag kind of 'plasmid' or 'seed_stock' is required."}), 400
    kind = 'plasmid' if kind == 'plasmid' else 'seed_stock'

    if not external_id:
        label = ('plasmid identifier' if kind == 'plasmid'
                 else 'seed stock identifier')
        return jsonify({'error': f'A {label} is required.'}), 400

    conn = db.connect(_db_path())
    try:
        if not _target_viewable(conn, member, target_type, target_id):
            # Same 404 whether the target is absent or merely unviewable, so a
            # member learns nothing about targets outside their reach.
            return jsonify({'error': _not_found_message(target_type)}), 404

        if kind == 'plasmid':
            return _add_plasmid_tag(conn, member, target_type, target_id,
                                    external_id)
        return _add_seed_stock_tag(conn, member, target_type, target_id,
                                   external_id)
    finally:
        conn.close()


def _add_plasmid_tag(conn, member, target_type, target_id, plasmid_id):
    """Validate against MoClo/Editor, reject missing/duplicate, else write a row."""
    # Existence in the referenced tool (Requirement 11.3). The missing id is
    # named in the error, and nothing is written.
    if not plasmid_exists(plasmid_id):
        return jsonify({
            'error': f'No plasmid with identifier "{plasmid_id}" exists in '
                     f'the MoClo/Editor tools.'}), 404

    # Duplicate on the same target (Requirement 11.4): reject, retain existing.
    existing = conn.execute(
        'SELECT id FROM plasmid_tags '
        'WHERE target_type = ? AND target_id = ? AND plasmid_id = ?',
        (target_type, target_id, plasmid_id)).fetchone()
    if existing is not None:
        return jsonify({
            'error': f'Plasmid "{plasmid_id}" is already tagged on this '
                     f'{target_type}.'}), 409

    tag_id = uuid.uuid4().hex
    try:
        conn.execute(
            'INSERT INTO plasmid_tags '
            '(id, target_type, target_id, plasmid_id, created_by) '
            'VALUES (?, ?, ?, ?, ?)',
            (tag_id, target_type, target_id, plasmid_id, member))
        conn.commit()
    except sqlite3.IntegrityError:
        # The UNIQUE(target_type, target_id, plasmid_id) constraint caught a
        # duplicate that slipped past the pre-check (a concurrent add). Same
        # duplicate response; the existing row is retained (Requirement 11.4).
        conn.rollback()
        return jsonify({
            'error': f'Plasmid "{plasmid_id}" is already tagged on this '
                     f'{target_type}.'}), 409

    row = conn.execute(
        'SELECT * FROM plasmid_tags WHERE id = ?', (tag_id,)).fetchone()
    return jsonify(_plasmid_tag_body(row)), 201


def _add_seed_stock_tag(conn, member, target_type, target_id, stock_id):
    """Validate against the Seed Bank, reject missing/duplicate, else write a row."""
    # Existence + kind in the Seed Bank (Requirement 11.3). For a line/plant the
    # stock id IS the Family_Tree key the Seed Bank reads (Requirement 11.5).
    exists, stock_kind = seed_stock_lookup(stock_id)
    if not exists:
        return jsonify({
            'error': f'No seed stock with identifier "{stock_id}" exists in '
                     f'the Seed Bank.'}), 404

    # Duplicate on the same target (Requirement 11.4): reject, retain existing.
    existing = conn.execute(
        'SELECT id FROM seed_stock_tags '
        'WHERE target_type = ? AND target_id = ? AND stock_id = ?',
        (target_type, target_id, stock_id)).fetchone()
    if existing is not None:
        return jsonify({
            'error': f'Seed stock "{stock_id}" is already tagged on this '
                     f'{target_type}.'}), 409

    tag_id = uuid.uuid4().hex
    try:
        conn.execute(
            'INSERT INTO seed_stock_tags '
            '(id, target_type, target_id, stock_kind, stock_id, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            (tag_id, target_type, target_id, stock_kind, stock_id, member))
        conn.commit()
    except sqlite3.IntegrityError:
        conn.rollback()
        return jsonify({
            'error': f'Seed stock "{stock_id}" is already tagged on this '
                     f'{target_type}.'}), 409

    row = conn.execute(
        'SELECT * FROM seed_stock_tags WHERE id = ?', (tag_id,)).fetchone()
    return jsonify(_seed_stock_tag_body(row)), 201


def _list_tags(target_type: str, target_id: str):
    """
    List the recorded tags for a target, each with a navigable link.

    Same view guard as adding: the member must be able to view the target, or
    the endpoint answers 404. Both plasmid and seed-stock tags are returned,
    each carrying a link to its MoClo/Editor or Seed Bank record (Requirement
    11.6).
    """
    member = identity.require_member()
    conn = db.connect(_db_path())
    try:
        if not _target_viewable(conn, member, target_type, target_id):
            return jsonify({'error': _not_found_message(target_type)}), 404

        plasmids = conn.execute(
            'SELECT * FROM plasmid_tags '
            'WHERE target_type = ? AND target_id = ? '
            'ORDER BY verified_at, id', (target_type, target_id)).fetchall()
        stocks = conn.execute(
            'SELECT * FROM seed_stock_tags '
            'WHERE target_type = ? AND target_id = ? '
            'ORDER BY verified_at, id', (target_type, target_id)).fetchall()
    finally:
        conn.close()

    return jsonify({
        'plasmid_tags': [_plasmid_tag_body(r) for r in plasmids],
        'seed_stock_tags': [_seed_stock_tag_body(r) for r in stocks],
    })


def _delete_tag(target_type: str, target_id: str, tag_id: str):
    """
    Remove one tag row from a target.

    The member must be able to view the target. The tag id is looked up in both
    tag tables scoped to this target; removing a ``seed_stock_tags`` row for a
    line/plant also removes the Family_Tree link, because the Seed Bank reads
    these rows by ``stock_id`` (Requirements 11.7, 11.8). Deleting a tag that is
    not there answers 404; a successful delete answers 204.
    """
    member = identity.require_member()
    conn = db.connect(_db_path())
    try:
        if not _target_viewable(conn, member, target_type, target_id):
            return jsonify({'error': _not_found_message(target_type)}), 404

        # Try both tables, scoped to this target so a tag id cannot be removed
        # from a target it does not belong to.
        for table in ('plasmid_tags', 'seed_stock_tags'):
            cur = conn.execute(
                f'DELETE FROM {table} '
                f'WHERE id = ? AND target_type = ? AND target_id = ?',
                (tag_id, target_type, target_id))
            if cur.rowcount:
                conn.commit()
                return '', 204
        conn.rollback()
        return jsonify({'error': 'Tag not found.'}), 404
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Target loading + view guard
# ---------------------------------------------------------------------------

def _target_viewable(conn, member: str, target_type: str,
                     target_id: str) -> bool:
    """
    Whether the member may view the tag target (a track or a genome).

    Tagging and listing are viewer actions: a member who can see the track or
    genome may annotate it. Visibility is resolved by ``access.can_view_track``
    / ``access.can_view_genome``, the single place the suite's two sharing
    schemes live, so tagging never reinvents the rules. A target that does not
    exist is not viewable.
    """
    if target_type == 'track':
        row = conn.execute(
            'SELECT * FROM tracks WHERE id = ?', (target_id,)).fetchone()
        return row is not None and access.can_view_track(member, row, conn)
    if target_type == 'genome':
        row = conn.execute(
            'SELECT * FROM reference_genomes WHERE id = ?',
            (target_id,)).fetchone()
        return row is not None and access.can_view_genome(member, row, conn)
    return False


def _not_found_message(target_type: str) -> str:
    """The 404 message for an absent or unviewable target."""
    if target_type == 'track':
        return 'Track not found'
    return 'Reference genome not found'


# ---------------------------------------------------------------------------
# Response bodies
# ---------------------------------------------------------------------------

def _plasmid_tag_body(row) -> dict:
    """A plasmid tag as JSON, with its navigable MoClo/Editor link (11.6)."""
    return {
        'id': row['id'],
        'kind': 'plasmid',
        'target_type': row['target_type'],
        'target_id': row['target_id'],
        'plasmid_id': row['plasmid_id'],
        'verified_at': row['verified_at'],
        'created_by': row['created_by'],
        'link': _plasmid_link(row['plasmid_id']),
    }


def _seed_stock_tag_body(row) -> dict:
    """A seed-stock tag as JSON, with its navigable Seed Bank link (11.6)."""
    return {
        'id': row['id'],
        'kind': 'seed_stock',
        'target_type': row['target_type'],
        'target_id': row['target_id'],
        'stock_kind': row['stock_kind'],
        'stock_id': row['stock_id'],
        'verified_at': row['verified_at'],
        'created_by': row['created_by'],
        'link': _seed_stock_link(row['stock_kind'], row['stock_id']),
    }


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _external_id(body: dict) -> str:
    """
    The external id a tag request names, from any of the accepted body keys.

    A plasmid tag may name ``plasmid_id``; a seed-stock tag may name
    ``stock_id``; both also accept a generic ``external_id`` / ``id`` so a
    caller can post ``{kind, id}``. Whitespace is trimmed.
    """
    for key in ('plasmid_id', 'stock_id', 'external_id', 'id'):
        value = body.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ''


def _db_path() -> Path:
    """The module's own database path, from the configured data root."""
    return Path(current_app.config['DATA_ROOT']) / 'genome.db'
