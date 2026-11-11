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
  2.9, 12.1, 12.2) An administrator may send ``suite_reference=true`` to make
  a suite reference instead: owned by no lab, visible to every member.
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

import gzip
import os
import sqlite3
import uuid
from typing import List

from flask import Flask, current_app, jsonify, request

import access
import db
import identity
import labs_reader
import readers
import serving
import sharing

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
                'SELECT g.id, g.name, g.assembly_id, g.owner_lab, '
                'g.genome_visibility, '
                '(SELECT COUNT(*) FROM chrom_sizes c WHERE c.genome_id = g.id) '
                '  AS chromosome_count, '
                '(SELECT COALESCE(SUM(length_bp), 0) FROM chrom_sizes c '
                '  WHERE c.genome_id = g.id) AS length_bp '
                'FROM reference_genomes g '
                'ORDER BY g.name COLLATE NOCASE, g.id').fetchall()
            # The landing page shows what each genome is and whose it is:
            # ``suite`` marks a suite reference every member sees.
            genomes = [
                {'id': row['id'], 'name': row['name'],
                 'assembly_id': row['assembly_id'],
                 'suite': row['genome_visibility'] == access.SUITE_VISIBILITY,
                 'chromosome_count': int(row['chromosome_count']),
                 'length_bp': int(row['length_bp'])}
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

        # A suite reference belongs to the whole suite, not a lab, and only an
        # administrator may make one; it needs no lab of the uploader's.
        suite_reference = _truthy(request.form.get('suite_reference'))
        if suite_reference:
            if not identity.is_admin(member):
                return jsonify({
                    'error': 'Only an administrator can add a suite '
                             'reference genome.'}), 403
            owner_lab = access.SUITE_OWNER
        else:
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
                             f'already exists for {"the suite" if suite_reference else "your lab"}.'}), 409

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

            # A gzipped FASTA (as most are published) is unpacked, under the
            # same 4 GB ceiling for what it unpacks to.
            if _is_gzip(fasta_path):
                packed = fasta_path.with_name('upload.fasta.gz')
                fasta_path.replace(packed)
                try:
                    with gzip.open(packed, 'rb') as source:
                        oversize = _stream_to_file(source, fasta_path,
                                                   MAX_FASTA_BYTES)
                except (OSError, EOFError) as exc:
                    _cleanup_dir(genome_dir)
                    return jsonify({'error': f'The file could not be '
                                             f'unpacked: {exc}.'}), 400
                finally:
                    packed.unlink(missing_ok=True)
                if oversize:
                    _cleanup_dir(genome_dir)
                    return jsonify({
                        'error': 'The FASTA unpacks to more than the maximum '
                                 'allowed size of 4 GB.'}), 413

            # Validate the FASTA (Requirement 2.5). FASTA is ASCII; decode
            # leniently so a stray byte surfaces as a validation reason rather
            # than a decode crash.
            text = fasta_path.read_text(encoding='utf-8', errors='replace')
            result = readers.parse_fasta(text)
            if not result.valid:
                _cleanup_dir(genome_dir)
                return jsonify({'error': result.reason}), 400

            visibility = (access.SUITE_VISIBILITY if suite_reference
                          else _default_visibility(owner_lab))

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
                             f'already exists for {"the suite" if suite_reference else "your lab"}.'}), 409
        finally:
            conn.close()

        # Index the sequence for the viewer now (serving.py); if samtools is
        # missing, the first view does it instead.
        serving.prepare_reference(fasta_path)

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
                'SELECT id, name, type, owner_lab, visibility, track_status, file_path, '
                'experiment, citation, doi, geo, created_at FROM tracks WHERE genome_id = ? '
                'ORDER BY name COLLATE NOCASE, id', (genome_id,)).fetchall()
            from tracks import _export_filename
            tracks = [
                {'id': row['id'], 'name': row['name'], 'type': row['type'],
                 'owner_lab': row['owner_lab'], 'visibility': row['visibility'],
                 'experiment': row['experiment'], 'created_at': row['created_at'],
                 'citation': row['citation'], 'doi': row['doi'], 'geo': row['geo'],
                 'can_change': access.is_owner_lab_member(member, row['owner_lab']) or identity.is_admin(member),
                 # The name it downloads as (Export), which the code to regenerate a figure uses.
                 'export_name': _export_filename(row['name'], row['type'], row['file_path'] or '')}
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
    Shared body of the share/unshare endpoints: chair groups and single labs.

    Both verbs resolve identity, load the genome under the same 404 guard,
    enforce the owner-lab guard, and validate every target before changing
    anything. A share made here is a purposeful one (``auto = 0``), and takes
    over an automatic share of the same target; revoking removes only
    purposeful shares, so a lab that still has a track shared with it on this
    genome keeps seeing the genome (``sharing.sync_genome_access``), and the
    answer says which.
    """
    member = identity.require_member()
    conn = db.connect(current_app.config['DATA_ROOT'] / 'genome.db')
    try:
        genome = _load_genome(conn, genome_id)
        if genome is None or not access.can_view_genome(member, genome, conn):
            return jsonify({'error': 'Reference genome not found'}), 404
        if genome['genome_visibility'] == access.SUITE_VISIBILITY:
            return jsonify({'error': 'A suite reference genome is seen by everyone already.'}), 400

        # Owner-lab guard first: a non-owner changes nothing (Requirement 12.6).
        try:
            access.require_owner_lab(member, genome['owner_lab'])
        except access.PermissionDenied as exc:
            return jsonify({'error': str(exc)}), 403

        groups, labs, error = _share_targets(request.get_json(silent=True))
        if error is not None:
            return jsonify({'error': error}), 400
        if labs_reader.in_demo(member) and (
                set(labs) - labs_reader.demo_lab_ids()
                or any(set(labs_reader.labs_in_chair_group(g)) - labs_reader.demo_lab_ids() for g in groups)):
            return jsonify({'error': 'The demo shares with the demo labs only.'}), 400

        # The owner lab's own chair group sees a genome by default. It is not a share row but the
        # genome's visibility: taking it away narrows the genome to the owner lab, adding it back
        # restores the default, so what was set when it was uploaded can be changed later.
        named_groups = list(groups)
        owner_group = labs_reader.chair_group_of(genome['owner_lab'])
        if owner_group and owner_group in groups:
            groups = [g for g in groups if g != owner_group]
            conn.execute('UPDATE reference_genomes SET genome_visibility = ? WHERE id = ?',
                         (VISIBILITY_CHAIR_GROUP if grant else VISIBILITY_OWNER_LAB, genome_id))
            conn.execute('DELETE FROM genome_shares WHERE genome_id = ? AND chair_group_id = ? AND auto = 0',
                         (genome_id, owner_group))
        if grant:
            for group in groups:
                conn.execute(
                    'INSERT INTO genome_shares (genome_id, chair_group_id, granted_by, auto) '
                    'VALUES (?, ?, ?, 0) ON CONFLICT (genome_id, chair_group_id) DO UPDATE SET auto = 0',
                    (genome_id, group, member))
            for lab in labs:
                conn.execute(
                    'INSERT INTO genome_lab_shares (genome_id, lab_id, granted_by, auto) '
                    'VALUES (?, ?, ?, 0) ON CONFLICT (genome_id, lab_id) DO UPDATE SET auto = 0',
                    (genome_id, lab, member))
        else:
            # Revoking a share that was never granted is harmless (Req 12.7).
            for group in groups:
                conn.execute('DELETE FROM genome_shares WHERE genome_id = ? AND chair_group_id = ? '
                             'AND auto = 0', (genome_id, group))
            for lab in labs:
                conn.execute('DELETE FROM genome_lab_shares WHERE genome_id = ? AND lab_id = ? '
                             'AND auto = 0', (genome_id, lab))
        sharing.sync_genome_access(conn, genome_id, member)
        conn.commit()
        genome = _load_genome(conn, genome_id)
        body = sharing.genome_share_body(conn, member, genome)
    finally:
        conn.close()

    # Kept for older callers: the one chair group a request named.
    groups = named_groups
    body['chair_group_id'] = groups[0] if len(groups) == 1 and not labs else None
    body['shared'] = grant
    if not grant:
        # Revoked on purpose but still needed for a shared track: it stays, automatic.
        body['kept'] = [x for x in body['shares']
                        if x.get('auto') and ((x['type'] == 'lab' and x['id'] in labs)
                                              or (x['type'] == 'chair_group' and x['id'] in groups))]
    return jsonify(body)


def _share_targets(body):
    """
    The chair groups and labs a share request names, or why it cannot be done.

    The body may name ``chair_group_id`` (one chair group), ``chair_groups`` (a
    list), ``labs`` (a list of lab ids), or a ``{target_type, target_id}`` pair
    whose type is ``'chair_group'`` or ``'lab'``. Every target must exist in
    ``suite.db``; one that does not rejects the whole request, so nothing
    changes (Requirement 12.8). Returns ``(groups, labs, None)`` or
    ``(None, None, error)``.
    """
    body = body or {}
    groups: List[str] = []
    labs: List[str] = []

    def text(v):
        return v.strip() if isinstance(v, str) else ''

    target_type = body.get('target_type')
    if target_type is not None:
        if target_type == 'chair_group':
            groups.append(text(body.get('target_id')))
        elif target_type == 'lab':
            labs.append(text(body.get('target_id')))
        else:
            return None, None, 'target_type must be "chair_group" or "lab".'
    if body.get('chair_group_id'):
        groups.append(text(body.get('chair_group_id')))
    for key, into in (('chair_groups', groups), ('labs', labs)):
        raw = body.get(key)
        if raw is None:
            continue
        if not isinstance(raw, list):
            return None, None, f'{key} must be a list.'
        into.extend(text(v) for v in raw)
    groups = list(dict.fromkeys(g for g in groups if g))
    labs = list(dict.fromkeys(lab for lab in labs if lab))
    if not groups and not labs:
        return None, None, 'Name a chair group or a lab to share with.'
    for group in groups:
        if not labs_reader.chair_group_exists(group):
            return None, None, (f'The chair group "{group}" is invalid; '
                                f'no such chair group exists.')
    known = {lab['id'] for lab in labs_reader.all_labs()}
    for lab in labs:
        if lab not in known:
            return None, None, f'The lab "{lab}" is invalid; no such lab exists.'
    return groups, labs, None


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


def _is_gzip(path) -> bool:
    with open(path, 'rb') as handle:
        return handle.read(2) == b'\x1f\x8b'


def _truthy(value) -> bool:
    return (value or '').strip().lower() in ('1', 'true', 'yes', 'on')


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
