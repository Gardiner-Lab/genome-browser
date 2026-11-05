"""
Track upload: ``POST /api/genomes/<id>/tracks`` — routing and the direct-store
path.

This module owns the *write* side of a genome's Track_Library, the sibling of
``genomes.py`` (which owns the reads). It is wired by :func:`register_tracks`,
called once from the app factory, exactly like ``register_genomes`` and
``register_figures``, so the factory stays a thin wiring layer and this module
and ``genomes.py`` never have to edit the factory body at the same time.

The one endpoint here accepts a track file for a reference genome, decides what
kind of track it is from the uploaded filename, and routes it (Requirement 4.6):

- **BigWig** — the direct-add coverage path. The file is written to a scratch
  location, its chromosome header is read (``readers.read_bigwig_header``) and
  checked against the genome's ``chrom_sizes`` by ``align.check_alignment``. On
  rejection (a chromosome name that does not match, Requirement 3.4, or a
  coordinate past a chromosome's end, Requirement 3.5) **nothing is persisted**:
  no ``tracks`` row and no stored file remain, and the alignment error naming
  every offender is returned. On acceptance the BigWig is moved under
  ``tracks/<id>/`` and a ``tracks`` row is created with ``type='bigwig'`` and
  ``track_status='ready'`` (Requirements 3.3, 4.1).

- **BED / GFF / GFF3** — the Feature_Track direct-store path. The file is stored
  under ``tracks/<id>/`` and a ``tracks`` row is created with ``type='bed'`` or
  ``type='gff'`` and ``track_status='ready'`` (Requirement 4.2). These formats
  carry no coordinate system the module validates on upload; they are stored
  as-is and read on demand by the region/export endpoints.

- **BAM / wig** — the conversion inputs. These are *not* stored directly and are
  *not* converted here: the module streams the upload straight to
  ``scratch/<job_id>/input.{bam,wig}`` (never into memory), enforcing the
  ``MAX_BAM_GB`` ceiling before any row is written, then creates a
  ``converting`` track and commits a ``queued`` ``conversion_jobs`` row through
  :func:`_enqueue_conversion` and ``conversion_queue.enqueue`` before returning
  202. The separate Conversion_Worker (task 9) later claims that job
  (``conversion_queue.claim_next_job``) and drives the track to ``ready`` or
  ``failed``. An oversize upload is rejected with no track, no job, and no
  scratch left behind (Requirement 13.14).

- **Anything else** — rejected with the supported-formats list
  (``readers.supported_formats_message``), Requirements 4.4 and 4.6.

Every created track is owned by the member's current lab (``Owner_Lab``) and
starts at ``Visibility='lab_only'`` (Requirements 4.1, 4.2, 4.3). Identity is
proxy-trusted via :func:`identity.require_member`; the genome must exist and be
viewable or the endpoint answers 404, the same guard the read endpoints use so
a member learns nothing about genomes outside their reach.

The 401 for an absent member is handled by the ``identity.Unauthenticated``
error handler already registered in ``genomes.py`` (the app shares one handler
for the whole module); this module does not re-register it and instead relies on
that handler.
"""

import re
import shutil
import sqlite3
import uuid
from pathlib import Path

from flask import Flask, current_app, jsonify, request, send_file

import access
import align
import conversion_queue
import db
import identity
import labs_reader
import readers

# One gigabyte in bytes, for turning the MAX_BAM_GB ceiling into a byte limit.
_BYTES_PER_GB = 1024 ** 3


def register_tracks(app: Flask) -> None:
    """
    Wire the track-upload endpoint onto the application.

    Called from ``create_app`` with a single line, mirroring
    ``register_genomes``/``register_figures``. Keeping the registration here
    lets this module own the upload route without touching the factory body or
    ``genomes.py``.
    """

    @app.route('/api/genomes/<genome_id>/tracks', methods=['POST'])
    def upload_track(genome_id):
        """
        Upload a track to a reference genome, routing by detected format.

        Requires a signed-in member, a viewable existing genome, a track name,
        and a supported file format. BigWig is alignment-checked then stored;
        BED/GFF/GFF3 are stored directly as feature tracks; BAM/wig branch to
        the conversion hook; an unsupported format is rejected with the list of
        supported formats. (Requirements 3.3, 3.4, 3.5, 4.1–4.6.)
        """
        member = identity.require_member()

        # A track name is required; an upload with no name is rejected before
        # any file is touched (Requirement 4.5).
        name = (request.form.get('name') or '').strip()
        if not name:
            return jsonify({'error': 'A track name is required.'}), 400

        upload = request.files.get('file')
        if upload is None or not (upload.filename or '').strip():
            return jsonify({'error': 'A track file is required.'}), 400

        # Detect the format from the filename before reading the (possibly huge)
        # body. An unknown format is rejected with the supported-formats list
        # (Requirements 4.4, 4.6).
        fmt = readers.detect_format(upload.filename)
        if fmt is None:
            return jsonify({
                'error': 'Unsupported track file format. '
                         + readers.supported_formats_message(),
            }), 400

        data_root = Path(current_app.config['DATA_ROOT'])
        conn = db.connect(data_root / 'genome.db')
        try:
            genome = _load_genome(conn, genome_id)
            if genome is None or not access.can_view_genome(member, genome, conn):
                # Same 404 whether the genome is absent or merely unviewable, so
                # a member learns nothing about genomes outside their reach.
                return jsonify({'error': 'Reference genome not found'}), 404

            owner_lab = labs_reader.current_lab(member)
            if not owner_lab:
                return jsonify({
                    'error': 'You must belong to a lab to upload a track.',
                }), 400

            if fmt == 'bigwig':
                return _add_bigwig(conn, data_root, genome, member, owner_lab,
                                   name, upload)
            if fmt in readers.DIRECT_STORE_FORMATS:  # bed / gff
                return _add_feature_track(conn, data_root, genome, member,
                                          owner_lab, name, fmt, upload)
            if fmt in readers.CONVERSION_FORMATS:  # bam / wig
                return _enqueue_conversion(conn, data_root, genome, member,
                                           owner_lab, name, fmt, upload)
            # Defensive: detect_format only returns known keys, so this is
            # unreachable, but keep the supported-formats contract explicit.
            return jsonify({
                'error': 'Unsupported track file format. '
                         + readers.supported_formats_message(),
            }), 400
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>', methods=['GET'])
    def track_detail(track_id):
        """
        A track's metadata, including its status and any failure reason.

        A member may read a track's metadata when they could view it
        (``access.can_view_track``, which requires ``track_status='ready'``) OR
        when they belong to the Owner_Lab — so the owning lab can watch its own
        ``converting`` or ``failed`` track and read the recorded error for
        status reporting even though such a track is viewable by no one
        (Requirement 13.19). When neither holds the same 404 is returned whether
        the track is absent or merely out of reach, so a member learns nothing
        about tracks outside it. The body carries id, name, type, Owner_Lab,
        Visibility, Track_Status and error (Requirements 3.6, 13.19).
        """
        member = identity.require_member()
        conn = db.connect(Path(current_app.config['DATA_ROOT']) / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            if track is None or not _may_read_metadata(member, track, conn):
                return jsonify({'error': 'Track not found'}), 404
            return jsonify(_track_metadata(track))
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>/status', methods=['GET'])
    def track_status(track_id):
        """
        A track's conversion status, plus the recorded error when it failed.

        Returns ``track_status`` — ``'converting'``, ``'ready'`` or
        ``'failed'`` — so a member can poll a conversion they uploaded, and the
        recorded ``error`` when the status is ``'failed'`` (Requirements 13.6,
        13.19). Owner_Lab members may poll a ``converting``/``failed`` track;
        viewers may poll a ``ready`` one. The same 404 as the metadata endpoint
        hides a track a member may not reach.
        """
        member = identity.require_member()
        conn = db.connect(Path(current_app.config['DATA_ROOT']) / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            if track is None or not _may_read_metadata(member, track, conn):
                return jsonify({'error': 'Track not found'}), 404
            body = {'track_status': track['track_status']}
            if track['track_status'] == 'failed':
                body['error'] = track['error']
            return jsonify(body)
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>/export', methods=['GET'])
    def export_track(track_id):
        """
        Download a track's stored file, named by the track's name.

        Export is a *view* action, so unlike the metadata/status endpoints it
        gates strictly on ``access.can_view_track`` — the owner-lab metadata
        exception does not apply here. A member who may not view the track is
        denied with an HTTP 403 and no bytes (Requirement 5.3); because
        ``can_view_track`` requires ``track_status='ready'``, a converting or
        failed track is not viewable and its export is denied the same way.

        When viewable, the stored file is streamed back as an attachment named
        ``<track name>.<ext>`` — the BigWig itself for a BigWig track
        (Requirement 5.1) and the uploaded feature file in its original format
        for a Feature_Track (Requirement 5.2). The download name uses the
        track's name, sanitised for use as a filename (Requirement 5.4).
        """
        member = identity.require_member()
        conn = db.connect(Path(current_app.config['DATA_ROOT']) / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            # Absent track reads as a 404; a real but unviewable track is a 403
            # with no bytes (Req 5.3). can_view_track requires status 'ready',
            # so a converting/failed track falls into the 403 branch too.
            if track is None:
                return jsonify({'error': 'Track not found'}), 404
            if not access.can_view_track(member, track, conn):
                return jsonify({
                    'error': 'You are not permitted to export this track.'}), 403

            file_path = track['file_path']
            if not file_path or not Path(file_path).exists():
                # Should not happen for a ready track, but fail clearly rather
                # than streaming a missing file.
                return jsonify({
                    'error': 'The track file is not available for export.'}), 404

            download_name = _export_filename(track['name'], track['type'],
                                             file_path)
            return send_file(
                str(file_path),
                mimetype=_export_mimetype(track['type']),
                as_attachment=True,
                download_name=download_name)
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>/share', methods=['POST'])
    def share_track(track_id):
        """
        Share a track with specific labs and/or the Owner_Lab's chair group.

        Only a member of the track's Owner_Lab may change its sharing: the
        owner-lab guard (``access.require_owner_lab``) runs first, and a member
        outside the Owner_Lab is rejected 403 with every share row left
        unchanged (Requirement 7.6). The JSON body may name ``labs`` (a list of
        lab ids to share with, ``target_type='lab'``) and/or ``chair_group:
        true`` to share with the Owner_Lab's chair group
        (``target_type='chair_group'``).

        Validation is all-or-nothing and persists nothing on rejection: a named
        lab that does not exist is rejected with an invalid-lab error, leaving
        Visibility unchanged (Requirement 7.8); a chair-group share while the
        Owner_Lab is in no chair group is rejected with an "owner lab is in no
        chair group" error (Requirement 7.9). On success the share rows are
        inserted idempotently (a re-share is a no-op, not an error), the track's
        ``visibility`` summary is recomputed from its current shares, and only
        this track changes (Requirements 7.3, 7.4, 7.5).
        """
        member = identity.require_member()
        conn = db.connect(Path(current_app.config['DATA_ROOT']) / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            if track is None:
                return jsonify({'error': 'Track not found'}), 404

            # Owner-lab guard first: an outsider changes nothing (Req 7.6). The
            # 404-vs-403 split is deliberate — an absent track is 404, a real
            # track the member may not administer is a clear 403.
            try:
                access.require_owner_lab(member, track['owner_lab'])
            except access.PermissionDenied as exc:
                return jsonify({'error': str(exc)}), 403

            body = request.get_json(silent=True) or {}
            labs = _requested_labs(body)
            want_chair_group = bool(body.get('chair_group'))
            if not labs and not want_chair_group:
                return jsonify({
                    'error': 'Name one or more labs and/or the chair group to '
                             'share this track with.'}), 400

            # Validate every lab target before touching any row, so a single
            # invalid lab leaves Visibility unchanged (Req 7.8).
            for lab_id in labs:
                if not _lab_exists(lab_id):
                    return jsonify({
                        'error': f'The lab "{lab_id}" is invalid: it does not '
                                 f'exist or may not be shared with.'}), 400

            # Validate the chair-group target before touching any row: the
            # Owner_Lab must be in a chair group (Req 7.9).
            chair_group_id = None
            if want_chair_group:
                chair_group_id = labs_reader.chair_group_of(track['owner_lab'])
                if not chair_group_id:
                    return jsonify({
                        'error': 'The owner lab is in no chair group, so the '
                                 'track cannot be shared with a chair group.'}), 400

            # Everything validated: insert idempotently. Changes affect only
            # this track (Req 7.3, 7.4, 7.5).
            for lab_id in labs:
                _insert_share(conn, track_id, 'lab', lab_id, member)
            if chair_group_id:
                _insert_share(conn, track_id, 'chair_group', chair_group_id,
                              member)
            _refresh_visibility(conn, track_id)
            conn.commit()

            track = _load_track(conn, track_id)
            return jsonify(_track_metadata(track))
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>/share', methods=['DELETE'])
    def unshare_track(track_id):
        """
        Revoke a previously granted share: a named lab or the chair group.

        The same owner-lab guard runs first (Requirement 7.6). The JSON body
        names what to revoke — ``labs`` (a list of lab ids) and/or
        ``chair_group: true`` for the Owner_Lab's chair group. Each matching
        ``track_shares`` row is deleted, so affected members immediately lose
        view while the Owner_Lab keeps it, and the ``visibility`` summary is
        recomputed from the remaining shares. Revoking a share that was never
        granted is a no-op, not an error. Only this track changes (Requirement
        7.7).
        """
        member = identity.require_member()
        conn = db.connect(Path(current_app.config['DATA_ROOT']) / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            if track is None:
                return jsonify({'error': 'Track not found'}), 404

            try:
                access.require_owner_lab(member, track['owner_lab'])
            except access.PermissionDenied as exc:
                return jsonify({'error': str(exc)}), 403

            body = request.get_json(silent=True) or {}
            labs = _requested_labs(body)
            revoke_chair_group = bool(body.get('chair_group'))
            if not labs and not revoke_chair_group:
                return jsonify({
                    'error': 'Name one or more labs and/or the chair group to '
                             'revoke sharing from.'}), 400

            for lab_id in labs:
                _delete_share(conn, track_id, 'lab', lab_id)
            if revoke_chair_group:
                # Revoke every chair-group share on this track; the Owner_Lab's
                # chair group is the only one a track can be shared with, so
                # this clears the chair-group grant regardless of its id.
                conn.execute(
                    'DELETE FROM track_shares '
                    'WHERE track_id = ? AND target_type = ?',
                    (track_id, 'chair_group'))
            _refresh_visibility(conn, track_id)
            conn.commit()

            track = _load_track(conn, track_id)
            return jsonify(_track_metadata(track))
        finally:
            conn.close()


# ---------------------------------------------------------------------------
# BigWig direct-add (Requirements 3.3, 3.4, 3.5, 4.1) — persist-nothing on reject
# ---------------------------------------------------------------------------


def _add_bigwig(conn, data_root, genome, member, owner_lab, name, upload):
    """
    The BigWig path: alignment-check, then store a ready track or persist nothing.

    The upload is written to a scratch file first (never stored in its final
    place until it is accepted), its header is read, and it is checked against
    the genome's ``chrom_sizes`` by :func:`align.check_alignment`. On rejection
    the scratch file is removed and no ``tracks`` row is written — the
    persist-nothing contract (Requirements 3.4, 3.5). On acceptance a track id
    is minted, the BigWig is moved into ``tracks/<id>/`` and a ``track_status``
    ``ready`` row is committed (Requirements 3.3, 4.1).
    """
    scratch_dir = data_root / 'scratch' / f'upload-{uuid.uuid4().hex}'
    scratch_dir.mkdir(parents=True, exist_ok=True)
    scratch_path = scratch_dir / _safe_name(upload.filename, default='track.bw')
    try:
        upload.save(str(scratch_path))

        header = readers.read_bigwig_header(str(scratch_path))
        chrom_sizes = _chrom_sizes(conn, genome['id'])
        result = align.check_alignment(header, chrom_sizes)
        if not result.accepted:
            # Persist nothing: no row, no stored file. Scratch is cleaned in the
            # finally below, so a rejected upload leaves no trace.
            return jsonify({
                'error': 'Track does not align to the reference genome.',
                'detail': result.error_message,
                'unmatched_chroms': result.unmatched_chroms,
                'out_of_range': [
                    {'chrom': v.chrom, 'coordinate': v.coordinate,
                     'length': v.length}
                    for v in result.out_of_range
                ],
            }), 400

        track_id = uuid.uuid4().hex
        stored_path = _store_file(data_root, track_id, scratch_path)
        _insert_track(conn, track_id, genome['id'], name, 'bigwig', owner_lab,
                      member, track_status='ready', file_path=str(stored_path))
        conn.commit()
        return jsonify(_track_body(track_id, name, 'bigwig', owner_lab)), 201
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# BED / GFF / GFF3 direct store (Requirement 4.2)
# ---------------------------------------------------------------------------


def _add_feature_track(conn, data_root, genome, member, owner_lab, name, fmt,
                       upload):
    """
    The Feature_Track path: store the uploaded BED/GFF/GFF3 as-is, ready.

    A track id is minted up front so the file lands in ``tracks/<id>/``; the
    stored ``type`` is ``'bed'`` or ``'gff'`` (the key ``readers`` detection
    returns), and the row is created ``ready`` with its ``file_path`` set
    (Requirement 4.2). Feature files carry no coordinate system the module
    validates on upload, so there is no alignment gate here.
    """
    track_id = uuid.uuid4().hex
    scratch_dir = data_root / 'scratch' / f'upload-{uuid.uuid4().hex}'
    scratch_dir.mkdir(parents=True, exist_ok=True)
    scratch_path = scratch_dir / _safe_name(upload.filename, default='track')
    try:
        upload.save(str(scratch_path))
        stored_path = _store_file(data_root, track_id, scratch_path)
        _insert_track(conn, track_id, genome['id'], name, fmt, owner_lab,
                      member, track_status='ready', file_path=str(stored_path))
        conn.commit()
        return jsonify(_track_body(track_id, name, fmt, owner_lab)), 201
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# BAM / wig conversion hook (task 8.1 fills this in)
# ---------------------------------------------------------------------------


def _enqueue_conversion(conn, data_root, genome, member, owner_lab, name, fmt,
                        upload):
    """
    Accept a BAM/wig upload into the durable conversion queue, then return 202.

    These formats are never stored directly: the design routes them through the
    durable conversion queue (``conversion_queue.py``) and the separate
    Conversion_Worker (task 9), which convert them to a BigWig server-side. This
    is the single place that path is entered from the upload endpoint.

    The order here is what makes an oversize upload cost nothing and an accepted
    job durable:

    1. **Stream to scratch under the size ceiling, before any row.** A fresh
       job id is minted and the body is streamed in chunks to
       ``scratch/<job_id>/input.{bam,wig}`` — never read into memory — while
       enforcing ``MAX_BAM_GB``. If the stream exceeds the ceiling the scratch
       directory is removed and the upload is rejected with HTTP 413, creating
       **no** track row and **no** job (Requirement 13.14).
    2. **Only then persist.** Once the original is fully written within the
       ceiling, create the ``converting`` track row (``type`` the final BigWig
       form, ``file_path`` NULL) and commit a ``queued`` ``conversion_jobs`` row
       via :func:`conversion_queue.enqueue` pointing at the scratch input. The job is
       committed before the response is sent, so an accepted conversion survives
       a crash/restart (Requirements 13.1, 13.2, 13.3, 13.20).
    3. **Return 202 immediately**, without waiting for the conversion; the
       worker drives the ``converting`` track to ``ready`` or ``failed``
       (Requirement 13.3).
    """
    max_bam_bytes = int(current_app.config['MAX_BAM_GB']) * _BYTES_PER_GB

    # A job id names the scratch directory from the start, so the committed job
    # row's input_path and the directory the worker cleans always agree.
    job_id = uuid.uuid4().hex
    scratch_dir = data_root / 'scratch' / job_id
    scratch_dir.mkdir(parents=True, exist_ok=True)
    input_path = scratch_dir / f'input.{"bam" if fmt == "bam" else "wig"}'

    # Stream the body to scratch under the ceiling BEFORE creating any row, so
    # an oversize upload leaves no track, no job, and no scratch (Req 13.14).
    oversize = _stream_to_file(upload.stream, input_path, max_bam_bytes)
    if oversize:
        shutil.rmtree(scratch_dir, ignore_errors=True)
        return jsonify({
            'error': 'The uploaded file exceeds the maximum allowed BAM size '
                     f'of {current_app.config["MAX_BAM_GB"]} GB.'}), 413

    # The original is within the ceiling: create the converting track and commit
    # a queued job keyed to it. The converted output is always a BigWig track;
    # type reflects the final stored form (bam/wig -> 'bigwig').
    track_id = uuid.uuid4().hex
    track_type = readers.track_type_for_format(fmt) or 'bigwig'
    _insert_track(conn, track_id, genome['id'], name, track_type, owner_lab,
                  member, track_status='converting', file_path=None)
    conn.commit()
    # Commit the queued job (conversion_queue.enqueue commits) so the accepted
    # conversion is durable before this response returns (Req 13.3, 13.20).
    conversion_queue.enqueue(conn, track_id, fmt, str(input_path))

    body = _track_body(track_id, name, track_type, owner_lab)
    body['track_status'] = 'converting'
    # 202 Accepted: conversion proceeds asynchronously on the worker.
    return jsonify(body), 202


def _stream_to_file(stream, dest, max_bytes: int) -> bool:
    """
    Copy ``stream`` into ``dest`` in chunks, stopping if it exceeds ``max_bytes``.

    Mirrors ``genomes.py``'s streamer: returns True when the stream exceeds the
    ceiling (the partial file is left for the caller to clean up with its
    scratch directory), False when it was written within the limit. Chunked so a
    tens-of-gigabytes BAM is never held in memory at once (Requirement 13.14).
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


# ---------------------------------------------------------------------------
# Storage + persistence helpers
# ---------------------------------------------------------------------------


def _store_file(data_root: Path, track_id: str, source: Path) -> Path:
    """
    Move an accepted upload into its permanent ``tracks/<id>/`` home.

    Returns the stored file path, which goes into ``tracks.file_path``. The
    per-track directory keeps each track's data isolated, matching the design's
    on-disk layout (``/genome/tracks/<id>/``).
    """
    track_dir = data_root / 'tracks' / track_id
    track_dir.mkdir(parents=True, exist_ok=True)
    stored_path = track_dir / source.name
    shutil.move(str(source), str(stored_path))
    return stored_path


def _insert_track(conn, track_id, genome_id, name, track_type, owner_lab,
                  created_by, track_status, file_path):
    """
    Write one ``tracks`` row.

    ``visibility`` is always ``'lab_only'`` on create (Requirement 4.3);
    ``track_status`` and ``file_path`` vary by path (``ready`` with a path for
    the direct-store formats, ``converting`` with no path for the conversion
    inputs).
    """
    conn.execute(
        'INSERT INTO tracks '
        '(id, genome_id, name, type, owner_lab, visibility, track_status, '
        ' file_path, created_by) '
        'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
        (track_id, genome_id, name, track_type, owner_lab, 'lab_only',
         track_status, file_path, created_by))


def _track_body(track_id, name, track_type, owner_lab):
    """The JSON a successful create returns, matching the read endpoints' shape."""
    return {
        'id': track_id,
        'name': name,
        'type': track_type,
        'owner_lab': owner_lab,
        'visibility': 'lab_only',
    }


def _chrom_sizes(conn, genome_id: str):
    """The genome's ``chrom -> length_bp`` map, for the alignment check."""
    rows = conn.execute(
        'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?',
        (genome_id,)).fetchall()
    return {row['chrom']: int(row['length_bp']) for row in rows}


def _load_genome(conn: sqlite3.Connection, genome_id: str):
    """The reference_genomes row for this id, or None when there is none."""
    return conn.execute(
        'SELECT id, name, assembly_id, owner_lab, genome_visibility '
        'FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()


# ---------------------------------------------------------------------------
# Track metadata / status / sharing helpers (Requirements 3.6, 7, 13.6, 13.19)
# ---------------------------------------------------------------------------


def _load_track(conn: sqlite3.Connection, track_id: str):
    """The tracks row for this id, or None when there is none."""
    return conn.execute(
        'SELECT id, genome_id, name, type, owner_lab, visibility, '
        'track_status, error, file_path FROM tracks WHERE id = ?',
        (track_id,)).fetchone()


def _may_read_metadata(member: str, track, conn) -> bool:
    """
    Whether this member may read a track's metadata/status.

    A member may read it when they could view the track
    (``access.can_view_track`` — which requires ``track_status='ready'``) OR
    when they belong to the Owner_Lab. The Owner_Lab branch is what lets the
    owning lab watch its own ``converting``/``failed`` track and read the
    recorded error, even though such a track is viewable by no one
    (Requirement 13.19).
    """
    if access.can_view_track(member, track, conn):
        return True
    return access.is_owner_lab_member(member, track['owner_lab'])


def _track_metadata(track) -> dict:
    """
    The metadata body for a track (Requirements 3.6, 13.19).

    Carries id, name, type, Owner_Lab, Visibility, Track_Status and the
    recorded error (``None`` unless the status is ``failed``).
    """
    return {
        'id': track['id'],
        'name': track['name'],
        'type': track['type'],
        'owner_lab': track['owner_lab'],
        'visibility': track['visibility'],
        'track_status': track['track_status'],
        'error': track['error'],
    }


def _requested_labs(body: dict):
    """
    The lab ids named in a share/unshare body, de-duplicated and non-empty.

    Accepts ``labs`` as a list of ids (or a single id string for convenience);
    blank entries are dropped. Order is preserved so an invalid lab is reported
    predictably.
    """
    raw = body.get('labs', [])
    if isinstance(raw, str):
        raw = [raw]
    seen = []
    for lab_id in raw:
        lab_id = (lab_id or '').strip() if isinstance(lab_id, str) else ''
        if lab_id and lab_id not in seen:
            seen.append(lab_id)
    return seen


def _lab_exists(lab_id: str) -> bool:
    """
    Whether a lab with this id exists in ``suite.db`` (read-only).

    ``labs_reader`` has no direct existence check — ``chair_group_of`` returns
    ``None`` both for an unknown lab and for a real lab in no chair group — so
    this reads ``suite.db`` read-only through ``db.open_suite`` to tell the two
    apart for the invalid-lab rejection (Requirement 7.8). An absent or
    unreadable store reads as "no such lab", matching ``labs_reader``'s
    degrade-to-empty behaviour.
    """
    if not lab_id:
        return False
    conn = db.open_suite()
    if conn is None:
        return False
    try:
        row = conn.execute('SELECT 1 FROM labs WHERE id = ?',
                           (lab_id,)).fetchone()
    except sqlite3.Error:
        return False
    finally:
        conn.close()
    return row is not None


def _insert_share(conn, track_id, target_type, target_id, granted_by):
    """
    Add one ``track_shares`` row idempotently.

    ``INSERT OR IGNORE`` against the ``(track_id, target_type, target_id)``
    primary key makes a re-share a no-op rather than an error, so sharing with a
    lab that already has access succeeds quietly.
    """
    conn.execute(
        'INSERT OR IGNORE INTO track_shares '
        '(track_id, target_type, target_id, granted_by) VALUES (?, ?, ?, ?)',
        (track_id, target_type, target_id, granted_by))


def _delete_share(conn, track_id, target_type, target_id):
    """Remove one ``track_shares`` row; a no-op when the share was never granted."""
    conn.execute(
        'DELETE FROM track_shares '
        'WHERE track_id = ? AND target_type = ? AND target_id = ?',
        (track_id, target_type, target_id))


def _refresh_visibility(conn, track_id) -> None:
    """
    Recompute a track's ``visibility`` summary from its current shares.

    The summary mirrors what the track is shared with: ``'chair_group'`` when a
    chair-group share exists, ``'labs'`` when one or more specific labs are
    shared with, otherwise back to the ``'lab_only'`` default. A chair-group
    share takes precedence in the summary because it is the broader grant; both
    kinds of share rows are kept regardless, so the summary is a label, not the
    source of truth (``track_shares`` is). Kept narrow to this track id so no
    other track changes (Requirements 7.3, 7.4, 7.5, 7.7).
    """
    rows = conn.execute(
        'SELECT DISTINCT target_type FROM track_shares WHERE track_id = ?',
        (track_id,)).fetchall()
    types = {r['target_type'] for r in rows}
    if 'chair_group' in types:
        visibility = 'chair_group'
    elif 'lab' in types:
        visibility = 'labs'
    else:
        visibility = 'lab_only'
    conn.execute('UPDATE tracks SET visibility = ? WHERE id = ?',
                 (visibility, track_id))


def _safe_name(filename: str, default: str) -> str:
    """
    A safe basename for an uploaded file, stripping any path components.

    Only the final name component is kept (werkzeug already discards most path
    parts, but strip defensively), so a crafted filename cannot escape the
    scratch/track directory.
    """
    base = Path((filename or '').replace('\\', '/')).name.strip()
    return base or default


# ---------------------------------------------------------------------------
# Export naming + content type (Requirements 5.1, 5.2, 5.4)
# ---------------------------------------------------------------------------

# The extension each stored track type carries when exported. A BigWig track is
# a ``.bw``; a feature track keeps the format it was uploaded in, so its
# extension is taken from the stored file rather than a fixed value here.
_EXPORT_EXTENSIONS = {'bigwig': '.bw'}

# Content types by stored track type; feature formats have no single registered
# type, so the octet-stream default below is used for them.
_EXPORT_MIMETYPES = {'bigwig': 'application/octet-stream'}


def _export_filename(track_name: str, track_type: str, file_path: str) -> str:
    """
    The download filename for an exported track: the track name + extension.

    The base name is the track's own name sanitised for a filesystem/header
    (Requirement 5.4); the extension reflects the track's type — ``.bw`` for a
    BigWig track, and for a feature track the extension of the stored uploaded
    file so it comes back in its original format (``.bed``/``.gff3``/...,
    Requirement 5.2).
    """
    base = _safe_filename_base(track_name) or 'track'
    ext = _EXPORT_EXTENSIONS.get(track_type)
    if ext is None:
        # Feature tracks keep their uploaded extension (bed/gff/gff3/...).
        ext = Path(file_path).suffix or ''
    return f'{base}{ext}'


def _export_mimetype(track_type: str) -> str:
    """The content type for an exported track; octet-stream for feature files."""
    return _EXPORT_MIMETYPES.get(track_type, 'application/octet-stream')


def _safe_filename_base(name: str) -> str:
    """
    A filename-safe base derived from a track name (no extension).

    Path separators and characters awkward in a filename or a
    ``Content-Disposition`` header are replaced with underscores, surrounding
    whitespace/dots are trimmed, and the result is length-capped so a long
    track name cannot produce an unwieldy header. An empty result is left for
    the caller to default.
    """
    cleaned = re.sub(r'[\\/:*?"<>|\r\n\t]+', '_', (name or '').strip())
    cleaned = cleaned.strip(' .')
    return cleaned[:150]
