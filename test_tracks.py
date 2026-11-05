"""
Tests for the track-upload routing and direct-store path (task 6.1).

These exercise ``POST /api/genomes/<id>/tracks`` end to end against a seeded
temporary database and data root: a track name is required, an unknown format
is rejected with the supported-formats list, a BigWig that aligns is stored as a
``ready`` track while one that fails alignment persists nothing, and BED/GFF are
stored as feature tracks. BigWig header reads are monkeypatched so the tests run
without pyBigWig or a real BigWig file, exercising the alignment-branch logic.
"""

import io
import uuid

import pytest

import align
import app as app_module
import db
import labs_reader
import readers


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def seeded(tmp_path, monkeypatch):
    """A built app + client over a temp data root with one viewable genome."""
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True

    # Seed one genome owned by lab 'labA' with two chromosomes, and make the
    # test member a member of labA so the genome is viewable and labA is the
    # Owner_Lab for created tracks.
    genome_id = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO reference_genomes '
            '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            (genome_id, 'Test Genome', 'asm1', 'labA', 'owner_lab_only', 'alice'))
        conn.executemany(
            'INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
            [(genome_id, 'Chr1', 1000), (genome_id, 'Chr2', 2000)])
        conn.commit()
    finally:
        conn.close()

    # alice is in labA; that makes the genome viewable and labA her current lab.
    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A',
                                    'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])

    client = application.test_client()

    def post(data):
        return client.post(
            f'/api/genomes/{genome_id}/tracks', data=data,
            headers={'X-Suite-User': 'alice'},
            content_type='multipart/form-data')

    return {'app': application, 'client': client, 'genome_id': genome_id,
            'data_root': data_root, 'post': post}


def _tracks(data_root):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute('SELECT * FROM tracks').fetchall()
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Name required (Req 4.5)
# ---------------------------------------------------------------------------


def test_missing_name_rejected(seeded):
    resp = seeded['post']({
        'file': (io.BytesIO(b'chr content'), 'coverage.bw'),
    })
    assert resp.status_code == 400
    assert 'name' in resp.get_json()['error'].lower()
    assert _tracks(seeded['data_root']) == []


def test_no_file_rejected(seeded):
    resp = seeded['post']({'name': 'My Track'})
    assert resp.status_code == 400
    assert _tracks(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Unknown format rejected with supported list (Req 4.4, 4.6)
# ---------------------------------------------------------------------------


def test_unknown_format_rejected_with_supported_list(seeded):
    resp = seeded['post']({
        'name': 'Mystery',
        'file': (io.BytesIO(b'???'), 'data.xyz'),
    })
    assert resp.status_code == 400
    error = resp.get_json()['error']
    # The rejection quotes the supported-formats message verbatim.
    assert readers.supported_formats_message() in error
    assert _tracks(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# BigWig alignment accept -> ready track stored (Req 3.3, 4.1, 4.3)
# ---------------------------------------------------------------------------


def test_bigwig_aligns_stored_ready(seeded, monkeypatch):
    # Header matches the genome's chromosomes and stays within length.
    monkeypatch.setattr(
        readers, 'read_bigwig_header',
        lambda path: readers.HeaderResult(
            chrom_lengths={'Chr1': 1000}, chrom_extents={'Chr1': 500}))

    resp = seeded['post']({
        'name': 'Coverage',
        'file': (io.BytesIO(b'fake bigwig bytes'), 'coverage.bw'),
    })
    assert resp.status_code == 201
    body = resp.get_json()
    assert body['type'] == 'bigwig'
    assert body['owner_lab'] == 'labA'
    assert body['visibility'] == 'lab_only'

    rows = _tracks(seeded['data_root'])
    assert len(rows) == 1
    row = rows[0]
    assert row['track_status'] == 'ready'
    assert row['type'] == 'bigwig'
    assert row['visibility'] == 'lab_only'
    assert row['owner_lab'] == 'labA'
    # The stored file exists under tracks/<id>/.
    assert row['file_path'] is not None
    from pathlib import Path
    assert Path(row['file_path']).exists()


# ---------------------------------------------------------------------------
# BigWig alignment reject -> persist nothing (Req 3.4, 3.5)
# ---------------------------------------------------------------------------


def test_bigwig_name_mismatch_persists_nothing(seeded, monkeypatch):
    # 'chr1' (lowercase) does not match 'Chr1' (case-sensitive) -> rejected.
    monkeypatch.setattr(
        readers, 'read_bigwig_header',
        lambda path: readers.HeaderResult(
            chrom_lengths={'chr1': 1000}, chrom_extents={'chr1': 500}))

    resp = seeded['post']({
        'name': 'Bad Coverage',
        'file': (io.BytesIO(b'fake bigwig bytes'), 'coverage.bw'),
    })
    assert resp.status_code == 400
    body = resp.get_json()
    assert 'chr1' in body['detail']
    assert body['unmatched_chroms'] == ['chr1']

    # Persist-nothing contract: no track row, and no stored/scratch files.
    assert _tracks(seeded['data_root']) == []
    tracks_dir = seeded['data_root'] / 'tracks'
    assert list(tracks_dir.iterdir()) == []
    scratch_dir = seeded['data_root'] / 'scratch'
    assert list(scratch_dir.iterdir()) == []


def test_bigwig_out_of_range_persists_nothing(seeded, monkeypatch):
    # Coordinate 1500 exceeds Chr1's length of 1000 -> rejected.
    monkeypatch.setattr(
        readers, 'read_bigwig_header',
        lambda path: readers.HeaderResult(
            chrom_lengths={'Chr1': 1500}, chrom_extents={'Chr1': 1500}))

    resp = seeded['post']({
        'name': 'Too Long',
        'file': (io.BytesIO(b'fake bigwig bytes'), 'coverage.bw'),
    })
    assert resp.status_code == 400
    body = resp.get_json()
    assert body['out_of_range'] and body['out_of_range'][0]['chrom'] == 'Chr1'
    assert _tracks(seeded['data_root']) == []
    assert list((seeded['data_root'] / 'tracks').iterdir()) == []


# ---------------------------------------------------------------------------
# BED / GFF stored as feature tracks (Req 4.2)
# ---------------------------------------------------------------------------


def test_bed_stored_as_feature_track(seeded):
    resp = seeded['post']({
        'name': 'Features',
        'file': (io.BytesIO(b'Chr1\t10\t20\tfeat\n'), 'features.bed'),
    })
    assert resp.status_code == 201
    assert resp.get_json()['type'] == 'bed'
    rows = _tracks(seeded['data_root'])
    assert len(rows) == 1
    assert rows[0]['type'] == 'bed'
    assert rows[0]['track_status'] == 'ready'
    from pathlib import Path
    assert Path(rows[0]['file_path']).exists()


def test_gff_stored_as_feature_track(seeded):
    resp = seeded['post']({
        'name': 'Genes',
        'file': (io.BytesIO(b'Chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n'),
                 'genes.gff3'),
    })
    assert resp.status_code == 201
    assert resp.get_json()['type'] == 'gff'
    rows = _tracks(seeded['data_root'])
    assert len(rows) == 1
    assert rows[0]['type'] == 'gff'
    assert rows[0]['track_status'] == 'ready'


# ---------------------------------------------------------------------------
# BAM / wig upload path (task 8.1): stream to scratch, converting track +
# queued conversion_jobs row, oversize rejected (Req 4.6, 4.7, 13.1-13.3,
# 13.14, 13.20)
# ---------------------------------------------------------------------------


def _jobs(data_root):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute(
            'SELECT * FROM conversion_jobs ORDER BY created_at, id').fetchall()
    finally:
        conn.close()


def test_bam_upload_streams_to_scratch_and_queues_job(seeded):
    content = b'BAM binary reads payload'
    resp = seeded['post']({
        'name': 'Reads',
        'file': (io.BytesIO(content), 'reads.bam'),
    })
    assert resp.status_code == 202
    body = resp.get_json()
    assert body['track_status'] == 'converting'
    assert body['type'] == 'bigwig'

    # The converting track exists with no file yet.
    rows = _tracks(seeded['data_root'])
    assert len(rows) == 1
    track = rows[0]
    assert track['track_status'] == 'converting'
    assert track['file_path'] is None

    # A queued conversion_jobs row keyed to the track, pointing at the streamed
    # scratch input, which holds the uploaded bytes verbatim (Req 13.1, 13.3).
    jobs = _jobs(seeded['data_root'])
    assert len(jobs) == 1
    job = jobs[0]
    assert job['track_id'] == track['id']
    assert job['kind'] == 'bam'
    assert job['state'] == 'queued'
    from pathlib import Path
    input_path = Path(job['input_path'])
    assert input_path.name == 'input.bam'
    assert input_path.parent.parent == seeded['data_root'] / 'scratch'
    assert input_path.exists()
    assert input_path.read_bytes() == content


def test_wig_upload_streams_to_scratch_and_queues_job(seeded):
    content = b'fixedStep chrom=Chr1 start=1 step=1\n1\n2\n'
    resp = seeded['post']({
        'name': 'Signal',
        'file': (io.BytesIO(content), 'signal.wig'),
    })
    assert resp.status_code == 202
    assert resp.get_json()['track_status'] == 'converting'

    jobs = _jobs(seeded['data_root'])
    assert len(jobs) == 1
    assert jobs[0]['kind'] == 'wig'
    assert jobs[0]['state'] == 'queued'
    from pathlib import Path
    input_path = Path(jobs[0]['input_path'])
    assert input_path.name == 'input.wig'
    assert input_path.read_bytes() == content


def test_oversize_bam_rejected_no_track_no_job_no_scratch(tmp_path, monkeypatch):
    """
    A BAM over MAX_BAM_GB is rejected (413) before any row, with no track, no
    job, and no leftover scratch directory (Req 13.14).
    """
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root
        MAX_BAM_GB = 0  # any non-empty body exceeds a 0 GB ceiling

    # 0 GB -> ceiling of 0 bytes, so even one byte is oversize. (The streamer
    # treats >max_bytes as oversize; a 0-byte body would pass, so send bytes.)
    application = app_module.create_app(Config())
    application.config['TESTING'] = True

    genome_id = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO reference_genomes '
            '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            (genome_id, 'G', 'asm1', 'labA', 'owner_lab_only', 'alice'))
        conn.commit()
    finally:
        conn.close()

    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A',
                                    'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)

    client = application.test_client()
    resp = client.post(
        f'/api/genomes/{genome_id}/tracks',
        data={'name': 'Huge', 'file': (io.BytesIO(b'way too big'), 'big.bam')},
        headers={'X-Suite-User': 'alice'},
        content_type='multipart/form-data')

    assert resp.status_code == 413
    assert 'maximum' in resp.get_json()['error'].lower()
    # Persist-nothing: no track, no job, and scratch is empty.
    assert _tracks(data_root) == []
    conn = db.connect(data_root / 'genome.db')
    try:
        assert conn.execute(
            'SELECT COUNT(*) AS c FROM conversion_jobs').fetchone()['c'] == 0
    finally:
        conn.close()
    assert list((data_root / 'scratch').iterdir()) == []


def test_bam_within_ceiling_is_accepted(tmp_path, monkeypatch):
    """A BAM under MAX_BAM_GB is accepted (the ceiling only rejects oversize)."""
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root
        MAX_BAM_GB = 1  # 1 GB ceiling; a tiny body is well under it

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    genome_id = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO reference_genomes '
            '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            (genome_id, 'G', 'asm1', 'labA', 'owner_lab_only', 'alice'))
        conn.commit()
    finally:
        conn.close()
    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A',
                                    'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)

    resp = application.test_client().post(
        f'/api/genomes/{genome_id}/tracks',
        data={'name': 'Reads', 'file': (io.BytesIO(b'small bam'), 'reads.bam')},
        headers={'X-Suite-User': 'alice'},
        content_type='multipart/form-data')
    assert resp.status_code == 202
    assert len(_tracks(data_root)) == 1


# ---------------------------------------------------------------------------
# Unviewable / missing genome -> 404
# ---------------------------------------------------------------------------


def test_unknown_genome_404(seeded):
    resp = seeded['client'].post(
        '/api/genomes/does-not-exist/tracks',
        data={'name': 'X', 'file': (io.BytesIO(b'x'), 'x.bw')},
        headers={'X-Suite-User': 'alice'},
        content_type='multipart/form-data')
    assert resp.status_code == 404


def test_unauthenticated_401(seeded):
    resp = seeded['client'].post(
        f"/api/genomes/{seeded['genome_id']}/tracks",
        data={'name': 'X', 'file': (io.BytesIO(b'x'), 'x.bw')},
        content_type='multipart/form-data')
    assert resp.status_code == 401


# ---------------------------------------------------------------------------
# Track metadata, status, and sharing (task 6.2; Req 3.6, 7.3–7.9, 13.6, 13.19)
# ---------------------------------------------------------------------------
#
# These exercise the GET /api/tracks/<id>, GET /api/tracks/<id>/status, and
# POST/DELETE /api/tracks/<id>/share endpoints against a seeded temp DB, with a
# small monkeypatched lab/chair-group universe (suite.db is never mounted in
# tests). The lab universe: labA (the Owner_Lab) and labB are in chair group
# cgX; labC is in no group and belongs to nobody here. alice is in labA; bob is
# in labB; carol is in labC.

import tracks as tracks_module


def _insert_track(data_root, genome_id, *, owner_lab='labA', visibility='lab_only',
                  track_status='ready', error=None, track_type='bigwig',
                  name='A Track'):
    """Insert one tracks row directly and return its id."""
    track_id = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO tracks '
            '(id, genome_id, name, type, owner_lab, visibility, track_status, '
            ' error, file_path, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            (track_id, genome_id, name, track_type, owner_lab, visibility,
             track_status, error, '/tmp/x', 'alice'))
        conn.commit()
    finally:
        conn.close()
    return track_id


def _shares(data_root, track_id):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute(
            'SELECT target_type, target_id FROM track_shares '
            'WHERE track_id = ? ORDER BY target_type, target_id',
            (track_id,)).fetchall()
    finally:
        conn.close()


@pytest.fixture
def lab_universe(seeded, monkeypatch):
    """A small lab/chair-group universe wired through the readers (no suite.db)."""
    members = {
        'alice': [{'lab_id': 'labA', 'lab': 'Lab A', 'role': 'member'}],
        'bob': [{'lab_id': 'labB', 'lab': 'Lab B', 'role': 'member'}],
        'carol': [{'lab_id': 'labC', 'lab': 'Lab C', 'role': 'member'}],
    }
    lab_group = {'labA': 'cgX', 'labB': 'cgX', 'labC': None}
    group_labs = {'cgX': ['labA', 'labB']}
    existing_labs = {'labA', 'labB', 'labC'}

    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: members.get(u, []))
    monkeypatch.setattr(labs_reader, 'chair_group_of',
                        lambda lab: lab_group.get(lab))
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group',
                        lambda g: list(group_labs.get(g, [])))
    # suite.db is not mounted in tests; resolve lab existence from the map.
    monkeypatch.setattr(tracks_module, '_lab_exists',
                        lambda lab_id: lab_id in existing_labs)
    return seeded


def _get(client, path, user):
    return client.get(path, headers={'X-Suite-User': user})


def _post_share(client, track_id, user, body):
    return client.post(f'/api/tracks/{track_id}/share', json=body,
                       headers={'X-Suite-User': user})


def _delete_share(client, track_id, user, body):
    return client.delete(f'/api/tracks/{track_id}/share', json=body,
                        headers={'X-Suite-User': user})


# --- GET /api/tracks/<id> metadata (Req 3.6, 13.19) ------------------------


def test_metadata_ready_track_visible_to_owner(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _get(s['client'], f'/api/tracks/{tid}', 'alice')
    assert resp.status_code == 200
    body = resp.get_json()
    assert body['id'] == tid
    assert body['owner_lab'] == 'labA'
    assert body['visibility'] == 'lab_only'
    assert body['track_status'] == 'ready'
    assert body['error'] is None


def test_metadata_owner_sees_converting_track(lab_universe):
    # Req 13.19: the owner lab sees its own converting track's metadata even
    # though a converting track is viewable by no one.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'], track_status='converting')
    resp = _get(s['client'], f'/api/tracks/{tid}', 'alice')
    assert resp.status_code == 200
    assert resp.get_json()['track_status'] == 'converting'


def test_metadata_owner_sees_failed_track_error(lab_universe):
    # Req 13.19: a failed track's recorded error is visible to the owner lab.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'], track_status='failed',
                        error='chrom Chr9 not found')
    resp = _get(s['client'], f'/api/tracks/{tid}', 'alice')
    assert resp.status_code == 200
    body = resp.get_json()
    assert body['track_status'] == 'failed'
    assert body['error'] == 'chrom Chr9 not found'


def test_metadata_non_owner_cannot_see_converting_track(lab_universe):
    # A non-owner (bob, labB) may not read a converting track's metadata.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'], track_status='converting')
    resp = _get(s['client'], f'/api/tracks/{tid}', 'bob')
    assert resp.status_code == 404


def test_metadata_unknown_track_404(lab_universe):
    s = lab_universe
    resp = _get(s['client'], '/api/tracks/does-not-exist', 'alice')
    assert resp.status_code == 404


def test_metadata_unauthenticated_401(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = s['client'].get(f'/api/tracks/{tid}')
    assert resp.status_code == 401


# --- GET /api/tracks/<id>/status (Req 13.6, 13.19) -------------------------


def test_status_ready(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _get(s['client'], f'/api/tracks/{tid}/status', 'alice')
    assert resp.status_code == 200
    body = resp.get_json()
    assert body['track_status'] == 'ready'
    assert 'error' not in body


def test_status_converting(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'], track_status='converting')
    resp = _get(s['client'], f'/api/tracks/{tid}/status', 'alice')
    assert resp.status_code == 200
    assert resp.get_json()['track_status'] == 'converting'


def test_status_failed_reports_error(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'], track_status='failed',
                        error='conversion failed: bad BAM header')
    resp = _get(s['client'], f'/api/tracks/{tid}/status', 'alice')
    assert resp.status_code == 200
    body = resp.get_json()
    assert body['track_status'] == 'failed'
    assert body['error'] == 'conversion failed: bad BAM header'


# --- POST /api/tracks/<id>/share: specific labs (Req 7.3, 7.5) -------------


def test_share_with_lab_grants_view(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    # Before sharing, bob (labB) cannot see the ready track.
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 404

    resp = _post_share(s['client'], tid, 'alice', {'labs': ['labB']})
    assert resp.status_code == 200
    assert resp.get_json()['visibility'] == 'labs'
    assert [tuple(r) for r in _shares(s['data_root'], tid)] == [('lab', 'labB')]

    # Now bob may view it (Req 7.3) while alice (Owner_Lab) still can.
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 200
    assert _get(s['client'], f'/api/tracks/{tid}', 'alice').status_code == 200


def test_share_is_idempotent(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    _post_share(s['client'], tid, 'alice', {'labs': ['labB']})
    resp = _post_share(s['client'], tid, 'alice', {'labs': ['labB']})
    assert resp.status_code == 200
    assert len(_shares(s['data_root'], tid)) == 1


def test_share_locality_other_track_unchanged(lab_universe):
    # Req 7.5: sharing one track leaves another track's visibility unchanged.
    s = lab_universe
    tid1 = _insert_track(s['data_root'], s['genome_id'], name='T1')
    tid2 = _insert_track(s['data_root'], s['genome_id'], name='T2')
    _post_share(s['client'], tid1, 'alice', {'labs': ['labB']})
    other = _get(s['client'], f'/api/tracks/{tid2}', 'alice').get_json()
    assert other['visibility'] == 'lab_only'
    assert _shares(s['data_root'], tid2) == []


# --- POST share: chair group (Req 7.4, 7.9) --------------------------------


def test_share_with_chair_group_grants_group_labs(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _post_share(s['client'], tid, 'alice', {'chair_group': True})
    assert resp.status_code == 200
    assert resp.get_json()['visibility'] == 'chair_group'
    assert [tuple(r) for r in _shares(s['data_root'], tid)] == \
        [('chair_group', 'cgX')]
    # bob is in labB, a lab of cgX, so bob may now view it (Req 7.4).
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 200


def test_share_chair_group_rejected_when_owner_in_no_group(lab_universe, monkeypatch):
    # Req 7.9: owner lab in no chair group -> reject with a clear message.
    s = lab_universe
    # Put the owner lab (labA) in no chair group for this test.
    monkeypatch.setattr(labs_reader, 'chair_group_of',
                        lambda lab: None)
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _post_share(s['client'], tid, 'alice', {'chair_group': True})
    assert resp.status_code == 400
    assert 'chair group' in resp.get_json()['error'].lower()
    # Visibility unchanged, no share rows written.
    assert _get(s['client'], f'/api/tracks/{tid}', 'alice') \
        .get_json()['visibility'] == 'lab_only'
    assert _shares(s['data_root'], tid) == []


# --- POST share: invalid lab + owner guard (Req 7.6, 7.8) ------------------


def test_share_invalid_lab_rejected_visibility_unchanged(lab_universe):
    # Req 7.8: an unknown lab is rejected and leaves visibility unchanged.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _post_share(s['client'], tid, 'alice', {'labs': ['labZ']})
    assert resp.status_code == 400
    assert 'invalid' in resp.get_json()['error'].lower()
    assert _get(s['client'], f'/api/tracks/{tid}', 'alice') \
        .get_json()['visibility'] == 'lab_only'
    assert _shares(s['data_root'], tid) == []


def test_share_invalid_lab_in_batch_persists_nothing(lab_universe):
    # A valid lab alongside an invalid one must not partially apply.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _post_share(s['client'], tid, 'alice', {'labs': ['labB', 'labZ']})
    assert resp.status_code == 400
    assert _shares(s['data_root'], tid) == []


def test_share_non_owner_rejected_403(lab_universe):
    # Req 7.6: a non-owner (bob) may not change sharing; nothing changes.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    resp = _post_share(s['client'], tid, 'bob', {'labs': ['labB']})
    assert resp.status_code == 403
    assert 'not permitted' in resp.get_json()['error'].lower()
    assert _shares(s['data_root'], tid) == []


# --- DELETE /api/tracks/<id>/share: revoke (Req 7.7) -----------------------


def test_revoke_lab_share_removes_view(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    _post_share(s['client'], tid, 'alice', {'labs': ['labB']})
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 200

    resp = _delete_share(s['client'], tid, 'alice', {'labs': ['labB']})
    assert resp.status_code == 200
    assert resp.get_json()['visibility'] == 'lab_only'
    # bob immediately loses view; alice (Owner_Lab) keeps it (Req 7.7).
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 404
    assert _get(s['client'], f'/api/tracks/{tid}', 'alice').status_code == 200
    assert _shares(s['data_root'], tid) == []


def test_revoke_chair_group_share(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    _post_share(s['client'], tid, 'alice', {'chair_group': True})
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 200

    resp = _delete_share(s['client'], tid, 'alice', {'chair_group': True})
    assert resp.status_code == 200
    assert resp.get_json()['visibility'] == 'lab_only'
    assert _get(s['client'], f'/api/tracks/{tid}', 'bob').status_code == 404
    assert _shares(s['data_root'], tid) == []


def test_revoke_non_owner_rejected_403(lab_universe):
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    _post_share(s['client'], tid, 'alice', {'labs': ['labB']})
    resp = _delete_share(s['client'], tid, 'bob', {'labs': ['labB']})
    assert resp.status_code == 403
    # The share is still there; bob changed nothing.
    assert [tuple(r) for r in _shares(s['data_root'], tid)] == [('lab', 'labB')]


def test_share_and_revoke_mixed_labs_and_chair_group(lab_universe):
    # Share with a lab and the chair group, then revoke just the lab.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])
    _post_share(s['client'], tid, 'alice',
                {'labs': ['labB'], 'chair_group': True})
    rows = {tuple(r) for r in _shares(s['data_root'], tid)}
    assert rows == {('lab', 'labB'), ('chair_group', 'cgX')}

    _delete_share(s['client'], tid, 'alice', {'labs': ['labB']})
    # Chair-group share remains, so visibility summarises as chair_group.
    assert _get(s['client'], f'/api/tracks/{tid}', 'alice') \
        .get_json()['visibility'] == 'chair_group'
    assert [tuple(r) for r in _shares(s['data_root'], tid)] == \
        [('chair_group', 'cgX')]


# --- GET /api/tracks/<id>/export (task 6.3; Req 5.1, 5.2, 5.3, 5.4) --------
#
# Export is a VIEW action, so it gates on access.can_view_track (which requires
# track_status='ready') rather than the owner-lab metadata exception. These
# tests store a real file under tracks/<id>/ so send_file has bytes to stream.


def _store_track_file(data_root, track_id, filename, content):
    """Write a real stored file for a track and return its path (str)."""
    track_dir = data_root / 'tracks' / track_id
    track_dir.mkdir(parents=True, exist_ok=True)
    path = track_dir / filename
    path.write_bytes(content)
    return str(path)


def _insert_track_with_file(data_root, genome_id, *, filename, content,
                            owner_lab='labA', visibility='lab_only',
                            track_status='ready', track_type='bigwig',
                            name='A Track'):
    """Insert a tracks row backed by a real stored file; return (id, path)."""
    track_id = uuid.uuid4().hex
    path = _store_track_file(data_root, track_id, filename, content)
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO tracks '
            '(id, genome_id, name, type, owner_lab, visibility, track_status, '
            ' error, file_path, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
            (track_id, genome_id, name, track_type, owner_lab, visibility,
             track_status, None, path, 'alice'))
        conn.commit()
    finally:
        conn.close()
    return track_id, path


def test_export_bigwig_returns_bytes_named_by_track(lab_universe):
    # Req 5.1, 5.4: a viewable BigWig track returns its stored bytes as an
    # attachment named by the track's name with a .bw extension.
    s = lab_universe
    content = b'fake bigwig payload'
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='coverage.bw',
        content=content, track_type='bigwig', name='My Coverage')

    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'alice')
    assert resp.status_code == 200
    assert resp.data == content
    disposition = resp.headers['Content-Disposition']
    assert 'attachment' in disposition
    assert 'My Coverage.bw' in disposition


def test_export_feature_track_returns_uploaded_file(lab_universe):
    # Req 5.2, 5.4: a feature track comes back in its uploaded format, named by
    # the track name with the original extension.
    s = lab_universe
    content = b'Chr1\t10\t20\tfeat\n'
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='features.bed',
        content=content, track_type='bed', name='Gene Models')

    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'alice')
    assert resp.status_code == 200
    assert resp.data == content
    assert 'Gene Models.bed' in resp.headers['Content-Disposition']


def test_export_gff_feature_track_keeps_extension(lab_universe):
    s = lab_universe
    content = b'Chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n'
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='genes.gff3',
        content=content, track_type='gff', name='Genes')

    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'alice')
    assert resp.status_code == 200
    assert resp.data == content
    assert 'Genes.gff3' in resp.headers['Content-Disposition']


def test_export_non_viewer_denied_403_no_bytes(lab_universe):
    # Req 5.3: a member who may not view the track is denied 403 with no bytes.
    s = lab_universe
    content = b'secret bigwig'
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='coverage.bw',
        content=content, track_type='bigwig', name='Private')

    # bob (labB) is not in the Owner_Lab and the track is lab_only.
    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'bob')
    assert resp.status_code == 403
    assert content not in resp.data
    assert 'error' in resp.get_json()


def test_export_converting_track_denied_even_for_owner(lab_universe):
    # Export is a view action: a converting track is viewable by no one
    # (can_view_track requires 'ready'), so even the owner lab is denied 403.
    s = lab_universe
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='coverage.bw',
        content=b'partial', track_type='bigwig', name='Pending',
        track_status='converting')

    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'alice')
    assert resp.status_code == 403


def test_export_shared_lab_member_can_download(lab_universe):
    # A member of a lab the track is shared with may export it (Req 5.1 via 7.3).
    s = lab_universe
    content = b'shared bigwig'
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='coverage.bw',
        content=content, track_type='bigwig', name='Shared')
    _post_share(s['client'], tid, 'alice', {'labs': ['labB']})

    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'bob')
    assert resp.status_code == 200
    assert resp.data == content


def test_export_unknown_track_404(lab_universe):
    s = lab_universe
    resp = _get(s['client'], '/api/tracks/does-not-exist/export', 'alice')
    assert resp.status_code == 404


def test_export_unauthenticated_401(lab_universe):
    s = lab_universe
    tid, _ = _insert_track_with_file(
        s['data_root'], s['genome_id'], filename='coverage.bw',
        content=b'x', track_type='bigwig')
    resp = s['client'].get(f'/api/tracks/{tid}/export')
    assert resp.status_code == 401


def test_export_missing_file_returns_404(lab_universe):
    # A ready track whose stored file is gone fails clearly rather than 500.
    s = lab_universe
    tid = _insert_track(s['data_root'], s['genome_id'])  # file_path='/tmp/x'
    conn = db.connect(s['data_root'] / 'genome.db')
    try:
        conn.execute('UPDATE tracks SET file_path = ? WHERE id = ?',
                     (str(s['data_root'] / 'tracks' / tid / 'gone.bw'), tid))
        conn.commit()
    finally:
        conn.close()
    resp = _get(s['client'], f'/api/tracks/{tid}/export', 'alice')
    assert resp.status_code == 404
