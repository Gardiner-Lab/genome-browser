"""
Reference-genome read endpoints: the genome list, a genome's chromosomes, and
a genome's Track_Library.

These are the first HTTP endpoints the module serves beyond ``/health``. They
live without the ``/genome`` prefix because the proxy strips it before passing
a request on (see ``app.py``). Identity is proxy-trusted: every endpoint here
requires the ``X-Suite-User`` header and answers 401 when it is absent, through
:func:`identity.require_member`. Visibility is never reinvented here — the
list and the Track_Library filter through ``access.can_view_genome`` and
``access.can_view_track``, the single place the suite's two sharing schemes are
resolved.

Later tasks (custom FASTA upload, genome sharing, the track endpoints) extend
this module and its sibling ``tracks.py`` rather than growing ``app.py``, so
the factory stays a thin wiring layer. This module exposes one entry point,
:func:`register_genomes`, which ``create_app`` calls.

Endpoints (all relative to the stripped ``/genome`` mount):

- ``GET /api/genomes`` — every Reference_Genome the member may view, each as
  ``{id, name, assembly_id}``, sorted alphabetically by name. (Requirement 2.1)
- ``GET /api/genomes/<id>`` — the genome's chromosome/contig names and their
  lengths in base pairs, from ``chrom_sizes``. 404 when the genome does not
  exist or the member may not view it. (Requirement 2.2)
- ``GET /api/genomes/<id>/tracks`` — the Track_Library: the genome's tracks the
  member may view and whose ``track_status='ready'``, each as
  ``{id, name, type, owner_lab, visibility}``; an empty list is a valid answer.
  (Requirements 2.7, 3.2, 3.6, 3.7)
- ``POST /api/genomes`` — upload a custom FASTA (≤4 GB) as a new Reference_Genome.
  Validates the FASTA, stores the assembly and its chromosome sizes, sets the
  Owner_Lab to the member's current lab, and records the default
  Genome_Visibility (the Owner_Lab's chair group, else owner-lab-only). Rejects
  an oversize file, an invalid FASTA (with the reason), and an assembly id that
  already exists in the same Owner_Lab. (Requirements 2.3, 2.4, 2.5, 2.6, 2.8,
  2.9, 12.1, 12.2)
- ``POST /api/genomes/<id>/share`` — share the genome with a whole chair group.
  Only a member of the Owner_Lab may change sharing (403 otherwise, nothing
  changed); the target must be a chair group, never an individual lab, and the
  chair group must exist. On success a ``genome_shares`` row is added (idempotent)
  so the labs of that chair group may view the genome on top of its existing
  default Genome_Visibility; only this genome is affected. (Requirements 12.3,
  12.4, 12.5, 12.6, 12.8)
- ``DELETE /api/genomes/<id>/share`` — revoke a previously granted chair-group
  share. Same owner-lab guard; removing the ``genome_shares`` row immediately
  stops the affected labs from viewing the genome while the default viewers keep
  it, and only this genome is affected. (Requirements 12.5, 12.6, 12.7)
"""

import os
import sqlite3
import uuid

from flask import Flask, current_app, jsonify, request

import access
import db
import identity
import labs_reader
import readers

# The FASTA upload ceiling (Requirement 2.4): a custom Reference_Genome may be
# no larger than 4 GB. Enforced here in the app so the limit holds regardless
# of any proxy body-size configuration.
MAX_FASTA_BYTES = 4 * 1024 * 1024 * 1024

# How a genome's default Genome_Visibility is recorded. A genome whose
# Owner_Lab is in a chair group defaults to that chair group; a genome whose
# Owner_Lab is in no chair group defaults to the Owner_Lab alone. (Requirements
# 2.9, 12.1, 12.2.)
VISIBILITY_CHAIR_GROUP = 'chair_group_default'
VISIBILITY_OWNER_LAB = 'owner_lab_only'


def register_genomes(app: Flask) -> None:
    """
    Wire the reference-genome read endpoints onto the application.

    Called from ``create_app``. Keeping the registration here (rather than in
    the factory) lets the later genome tasks extend one module and leaves
    ``app.py`` to wiring only.
    """

    @app.errorhandler(identity.Unauthenticated)
    def _unauthenticated(exc):
        # Behind the proxy the header is always set; this guards the
        # unproxied/misconfigured case with the 401 the design calls for.
        return jsonify({'error': str(exc) or 'Authentication required'}), 401

    @app.route('/api/genomes')
    def list_genomes():
        """
        Every genome the signed-in member may view, sorted by name.

        Returns a JSON list of ``{id, name, assembly_id}``. Membership and
        sharing are resolved by ``access.can_view_genome`` over one shared
        connection; the sort is alphabetical by name, done in SQL with a
        case-insensitive collation so ordering is stable regardless of case.
        (Requirement 2.1.)
        """
        member = identity.require_member()
        conn = db.connect(current_app.config['DATA_ROOT'] / 'genome.db')
        try:
            rows = conn.execute(
                'SELECT id, name, assembly_id, owner_lab, genome_visibility '
                'FROM reference_genomes ORDER BY name COLLATE NOCASE, id').fetchall()
            genomes = [
                {'id': row['id'], 'name': row['name'],
                 'assembly_id': row['assembly_id']}
                for row in rows
                if access.can_view_genome(member, row, conn)
            ]
        finally:
            conn.close()
        return jsonify(genomes)

    @app.route('/api/genomes', methods=['POST'])
    def upload_genome():
        """
        Upload a custom FASTA as a new Reference_Genome.

        The member's current lab is the Owner_Lab (Requirement 2.3); a member in
        no lab is rejected, because a genome has no owner otherwise. The upload
        carries a ``name`` and ``assembly_id`` as form fields and the FASTA as a
        file part; both the name and assembly id are required.

        Validation is fail-fast and persists nothing until every check passes:
        a file larger than 4 GB is rejected with a size error (2.4); a file that
        is not valid FASTA is rejected with the reason FASTA validation failed
        (2.5); an assembly id that already exists in the same Owner_Lab is
        rejected as a duplicate (2.6). On success the assembly is stored with a
        fresh immutable UUID id (2.8), its chromosome/contig sizes, the Owner_Lab
        (2.3), and the default Genome_Visibility — the Owner_Lab's chair group
        when it has one, else the Owner_Lab alone (2.9, 12.1, 12.2). The FASTA is
        kept under ``DATA_ROOT/genomes/<id>/`` per the volume layout.
        """
        member = identity.require_member()

        owner_lab = labs_reader.current_lab(member)
        if not owner_lab:
            return jsonify({
                'error': 'You are not a member of any lab, so there is no '
                         'owner lab to assign this genome to.'}), 400

        name = (request.form.get('name') or '').strip()
        assembly_id = (request.form.get('assembly_id') or '').strip()
        if not name:
            return jsonify({'error': 'A genome name is required.'}), 400
        if not assembly_id:
            return jsonify({'error': 'An assembly identifier is required.'}), 400

        upload = request.files.get('file')
        if upload is None or not (upload.filename or '').strip():
            return jsonify({'error': 'A FASTA file is required.'}), 400

        data_root = current_app.config['DATA_ROOT']
        db_path = data_root / 'genome.db'

        # Reject a duplicate assembly id for this Owner_Lab up front, with a
        # clear message, rather than relying only on the UNIQUE constraint
        # (which also guards the race) (Requirement 2.6).
        conn = db.connect(db_path)
        try:
            if _assembly_exists(conn, owner_lab, assembly_id):
                return jsonify({
                    'error': f'An assembly with identifier "{assembly_id}" '
                             f'already exists for your lab.'}), 409

            # Stream the body to a temp file while enforcing the 4 GB ceiling,
            # so an oversize file is rejected without being read into memory
            # (Requirement 2.4).
            genome_id = str(uuid.uuid4())
            genome_dir = data_root / 'genomes' / genome_id
            genome_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = genome_dir / 'reference.fasta'

            oversize = _stream_to_file(upload.stream, fasta_path, MAX_FASTA_BYTES)
            if oversize:
                _cleanup_dir(genome_dir)
                return jsonify({
                    'error': 'The uploaded file exceeds the maximum allowed '
                             'size of 4 GB.'}), 413

            # Validate the FASTA (Requirement 2.5). FASTA is ASCII; decode
            # leniently so a stray byte surfaces as a validation reason rather
            # than a decode crash.
            text = fasta_path.read_text(encoding='utf-8', errors='replace')
            result = readers.parse_fasta(text)
            if not result.valid:
                _cleanup_dir(genome_dir)
                return jsonify({'error': result.reason}), 400

            visibility = _default_visibility(owner_lab)

            try:
                _store_genome(conn, genome_id, name, assembly_id, owner_lab,
                              visibility, member, result.contigs)
            except sqlite3.IntegrityError:
                # The UNIQUE(owner_lab, assembly_id) constraint caught a
                # duplicate that slipped past the pre-check (a concurrent
                # upload). Same member-facing duplicate error (Requirement 2.6).
                conn.rollback()
                _cleanup_dir(genome_dir)
                return jsonify({
                    'error': f'An assembly with identifier "{assembly_id}" '
                             f'already exists for your lab.'}), 409
        finally:
            conn.close()

        body = {
            'id': genome_id,
            'name': name,
            'assembly_id': assembly_id,
            'owner_lab': owner_lab,
            'genome_visibility': visibility,
            'chromosomes': [
                {'name': chrom, 'length_bp': length}
                for chrom, length in result.contigs
            ],
        }
        return jsonify(body), 201

    @app.route('/api/genomes/<genome_id>')
    def genome_detail(genome_id):
        """
        A genome's chromosome/contig names and lengths in base pairs.

        404 when the genome does not exist, and the same 404 when it exists but
        the member may not view it — a member learns nothing about genomes
        outside their reach. Chromosomes come straight from ``chrom_sizes`` in a
        stable order (by name). (Requirement 2.2.)
        """
        member = identity.require_member()
        conn = db.connect(current_app.config['DATA_ROOT'] / 'genome.db')
        try:
            genome = _load_genome(conn, genome_id)
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            chroms = conn.execute(
                'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ? '
                'ORDER BY chrom COLLATE NOCASE, chrom', (genome_id,)).fetchall()
        finally:
            conn.close()
        body = {
            'id': genome['id'],
            'name': genome['name'],
            'assembly_id': genome['assembly_id'],
            'chromosomes': [
                {'name': c['chrom'], 'length_bp': c['length_bp']} for c in chroms
            ],
        }
        return jsonify(body)

    @app.route('/api/genomes/<genome_id>/tracks')
    def genome_tracks(genome_id):
        """
        The Track_Library for a genome: its viewable, ready tracks.

        404 when the genome is not visible to the member (same guard as the
        detail endpoint). Otherwise a JSON list — possibly empty — of the
        genome's tracks the member may view, each as
        ``{id, name, type, owner_lab, visibility}``. ``access.can_view_track``
        already requires ``track_status='ready'``, so converting and failed
        tracks never appear. (Requirements 2.7, 3.2, 3.6, 3.7.)
        """
        member = identity.require_member()
        conn = db.connect(current_app.config['DATA_ROOT'] / 'genome.db')
        try:
            genome = _load_genome(conn, genome_id)
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            rows = conn.execute(
                'SELECT id, name, type, owner_lab, visibility, track_status '
                'FROM tracks WHERE genome_id = ? '
                'ORDER BY name COLLATE NOCASE, id', (genome_id,)).fetchall()
            tracks = [
                {'id': row['id'], 'name': row['name'], 'type': row['type'],
                 'owner_lab': row['owner_lab'], 'visibility': row['visibility']}
                for row in rows
                if access.can_view_track(member, row, conn)
            ]
        finally:
            conn.close()
        return jsonify(tracks)

    @app.route('/api/genomes/<genome_id>/share', methods=['POST'])
    def share_genome(genome_id):
        """
        Share a genome with a whole chair group.

        The genome must exist and be viewable by the member (same 404 guard as
        the detail endpoint — a member learns nothing about genomes outside
        their reach). Only a member of the Owner_Lab may change sharing: the
        owner-lab guard runs first via ``access.require_owner_lab`` and a
        non-owner gets a 403 with the share rows untouched (Requirement 12.6).

        A Reference_Genome accepts chair groups as sharing targets only, never an
        individual lab, so an explicit ``'lab'`` target is rejected (Requirement
        12.4); and a chair group that does not exist in ``suite.db`` is rejected
        as invalid (Requirement 12.8). On success a ``genome_shares`` row is
        added — idempotently, so re-sharing the same chair group is harmless —
        letting the labs of that chair group view the genome on top of its
        existing default Genome_Visibility (Requirement 12.3). Only this genome
        is affected (Requirement 12.5).
        """
        return _mutate_genome_share(genome_id, grant=True)

    @app.route('/api/genomes/<genome_id>/share', methods=['DELETE'])
    def unshare_genome(genome_id):
        """
        Revoke a previously granted chair-group share on a genome.

        Same existence/visibility 404 guard and owner-lab guard as sharing: a
        non-owner gets a 403 with nothing changed (Requirement 12.6). Removing
        the ``genome_shares`` row immediately stops the labs of that chair group
        from viewing the genome, while the members defined by its default
        Genome_Visibility keep it (Requirement 12.7); revoking a share that was
        never granted is harmless. Only this genome is affected (Requirement
        12.5).
        """
        return _mutate_genome_share(genome_id, grant=False)


def _mutate_genome_share(genome_id, grant: bool):
    """
    Shared body of the share/unshare endpoints (grant or revoke a chair-group).

    Both verbs resolve identity, load the genome under the same 404 guard,
    enforce the owner-lab guard, and validate the chair-group target the same
    way; they differ only in whether they insert or delete the
    ``genome_shares`` row. Keeping the two in one place means the guards cannot
    drift apart.
    """
    member = identity.require_member()
    conn = db.connect(current_app.config['DATA_ROOT'] / 'genome.db')
    try:
        genome = _load_genome(conn, genome_id)
        if genome is None or not access.can_view_genome(member, genome, conn):
            return jsonify({'error': 'Reference genome not found'}), 404

        # Owner-lab guard first: a non-owner changes nothing (Requirement 12.6).
        try:
            access.require_owner_lab(member, genome['owner_lab'])
        except access.PermissionDenied as exc:
            return jsonify({'error': str(exc)}), 403

        chair_group_id, error = _share_target(request.get_json(silent=True))
        if error is not None:
            return jsonify({'error': error}), 400

        if grant:
            # Idempotent via the (genome_id, chair_group_id) primary key, so
            # re-sharing the same chair group preserves the existing row and
            # changes nothing else (Requirements 12.3, 12.5).
            conn.execute(
                'INSERT OR IGNORE INTO genome_shares '
                '(genome_id, chair_group_id, granted_by) VALUES (?, ?, ?)',
                (genome_id, chair_group_id, member))
        else:
            # Revoking a share that was never granted is harmless (Req 12.7).
            conn.execute(
                'DELETE FROM genome_shares '
                'WHERE genome_id = ? AND chair_group_id = ?',
                (genome_id, chair_group_id))
        conn.commit()
    finally:
        conn.close()

    return jsonify({
        'id': genome_id,
        'chair_group_id': chair_group_id,
        'shared': grant,
    })


def _share_target(body):
    """
    The chair-group id a share request names, or an error describing why not.

    A Reference_Genome's sharing targets are chair groups only, never individual
    labs (Requirement 12.4). The request body may name the target either as a
    bare ``chair_group_id`` or as a ``{target_type, target_id}`` pair; when a
    ``target_type`` is given it must be ``'chair_group'`` so an attempt to share
    with an individual lab is rejected explicitly. The named chair group must
    exist in ``suite.db`` or the request is rejected as invalid (Requirement
    12.8). Returns ``(chair_group_id, None)`` on success, else
    ``(None, error_message)``.
    """
    body = body or {}

    target_type = body.get('target_type')
    if target_type is not None and target_type != 'chair_group':
        # An explicit non-chair-group target — an individual lab — is never a
        # valid sharing target for a Reference_Genome (Requirement 12.4).
        return None, ('A reference genome can only be shared with a chair '
                      'group, not with an individual lab.')

    chair_group_id = (
        body.get('chair_group_id')
        or (body.get('target_id') if target_type == 'chair_group' else None))
    chair_group_id = (chair_group_id or '').strip() if isinstance(
        chair_group_id, str) else chair_group_id
    if not chair_group_id:
        return None, 'A chair group to share with is required.'

    if not labs_reader.chair_group_exists(chair_group_id):
        return None, (f'The chair group "{chair_group_id}" is invalid; '
                      f'no such chair group exists.')

    return chair_group_id, None


def _load_genome(conn: sqlite3.Connection, genome_id: str):
    """The reference_genomes row for this id, or None when there is none."""
    return conn.execute(
        'SELECT id, name, assembly_id, owner_lab, genome_visibility '
        'FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()


def _assembly_exists(conn: sqlite3.Connection, owner_lab: str,
                     assembly_id: str) -> bool:
    """Whether this Owner_Lab already has a genome with this assembly id (2.6)."""
    row = conn.execute(
        'SELECT 1 FROM reference_genomes WHERE owner_lab = ? AND assembly_id = ?',
        (owner_lab, assembly_id)).fetchone()
    return row is not None


def _default_visibility(owner_lab: str) -> str:
    """
    The default Genome_Visibility for a genome owned by this lab.

    ``chair_group_default`` when the Owner_Lab belongs to a chair group (every
    lab in that chair group may view it), otherwise ``owner_lab_only`` (the
    Owner_Lab alone). (Requirements 2.9, 12.1, 12.2.)
    """
    if labs_reader.chair_group_of(owner_lab):
        return VISIBILITY_CHAIR_GROUP
    return VISIBILITY_OWNER_LAB


def _stream_to_file(stream, dest, max_bytes: int) -> bool:
    """
    Copy ``stream`` into ``dest`` in chunks, stopping if it exceeds ``max_bytes``.

    Returns True when the stream exceeds the ceiling (the partial file is left
    for the caller to clean up), False when it was written within the limit.
    Chunked so a huge upload is never held in memory at once (Requirement 2.4).
    """
    written = 0
    chunk_size = 1024 * 1024
    with open(dest, 'wb') as handle:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                return False
            written += len(chunk)
            if written > max_bytes:
                return True
            handle.write(chunk)


def _store_genome(conn: sqlite3.Connection, genome_id: str, name: str,
                  assembly_id: str, owner_lab: str, visibility: str,
                  created_by: str, contigs) -> None:
    """
    Insert the assembly row and its chrom_sizes rows in one transaction.

    The ``reference_genomes`` insert can raise ``sqlite3.IntegrityError`` on the
    ``UNIQUE(owner_lab, assembly_id)`` constraint, which the caller turns into
    the duplicate-assembly error (Requirement 2.6). Chromosome names are stored
    exactly as parsed (case-sensitive), matching ``chrom_sizes`` elsewhere.
    """
    conn.execute(
        'INSERT INTO reference_genomes '
        '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
        'VALUES (?, ?, ?, ?, ?, ?)',
        (genome_id, name, assembly_id, owner_lab, visibility, created_by))
    conn.executemany(
        'INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
        [(genome_id, chrom, length) for chrom, length in contigs])
    conn.commit()


def _cleanup_dir(path) -> None:
    """Remove a partially-written genome directory after a rejected upload."""
    try:
        for entry in os.scandir(path):
            try:
                os.remove(entry.path)
            except OSError:
                pass
        os.rmdir(path)
    except OSError:
        pass
