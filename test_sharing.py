"""
Sharing between labs: a shared track takes its genome with it.

labA owns a genome, seen by labA alone, and a track on it. labX (carol) and
labY (dave) see neither until something is shared with them.
"""

import uuid

import pytest

import access
import app as app_module
import db
import identity
import labs_reader


@pytest.fixture
def world(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True

    genome, suite_genome = uuid.uuid4().hex, uuid.uuid4().hex
    track, track2, suite_track = uuid.uuid4().hex, uuid.uuid4().hex, uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    conn.executemany(
        'INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, '
        'created_by) VALUES (?, ?, ?, ?, ?, ?)',
        [(genome, 'Genome A', 'A1', 'labA', 'owner_lab_only', 'alice'),
         (suite_genome, 'Suite', 'S1', 'suite', 'suite', 'admin')])
    conn.executemany(
        'INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, '
        'file_path, created_by) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
        [(track, genome, 'Peaks', 'bed', 'labA', 'lab_only', 'ready', None, 'alice'),
         (track2, genome, 'Peaks 2', 'bed', 'labA', 'lab_only', 'ready', None, 'alice'),
         (suite_track, suite_genome, 'Mine on suite', 'bed', 'labA', 'lab_only', 'ready', None, 'alice')])
    conn.commit()
    conn.close()

    members = {'alice': 'labA', 'carol': 'labX', 'dave': 'labY'}
    lab_group = {'labA': 'cgA', 'labX': 'cgX', 'labY': 'cgX'}
    groups = {'cgA': ['labA'], 'cgX': ['labX', 'labY']}
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: [
        {'lab_id': members[u], 'lab': members[u], 'role': 'member'}] if u in members else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: lab_group.get(lab))
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: groups.get(g, []))
    monkeypatch.setattr(labs_reader, 'chair_group_exists', lambda g: g in groups)
    monkeypatch.setattr(labs_reader, 'all_labs', lambda: [
        {'id': lab, 'name': 'Lab ' + lab[-1], 'chair_group_id': g} for lab, g in lab_group.items()])
    monkeypatch.setattr(labs_reader, 'all_chair_groups', lambda: [
        {'id': g, 'name': 'Group ' + g[-1]} for g in groups])
    monkeypatch.setattr(labs_reader, 'lab_name', lambda lab: 'Lab ' + lab[-1])
    monkeypatch.setattr(labs_reader, 'chair_group_name', lambda g: 'Group ' + g[-1])
    client = application.test_client()

    def call(method, path, user='alice', body=None):
        return client.open(path, method=method, json=body, headers={'X-Suite-User': user})

    return {'call': call, 'root': data_root, 'genome': genome, 'track': track,
            'track2': track2, 'suite_track': suite_track, 'suite_genome': suite_genome}


def _sees_genome(w, user, genome=None):
    return w['call']('GET', f"/api/genomes/{genome or w['genome']}", user=user).status_code == 200


def _sees_track(w, user, track=None):
    library = w['call']('GET', f"/api/genomes/{w['genome']}/tracks", user=user)
    return library.status_code == 200 and any(
        t['id'] == (track or w['track']) for t in library.get_json())


def _lab_shares(w):
    conn = db.connect(w['root'] / 'genome.db')
    try:
        return [(r['lab_id'], r['auto']) for r in conn.execute(
            'SELECT lab_id, auto FROM genome_lab_shares WHERE genome_id = ? ORDER BY lab_id',
            (w['genome'],))]
    finally:
        conn.close()


def test_sharing_a_track_shares_its_genome_with_that_lab(world):
    assert not _sees_genome(world, 'carol')
    resp = world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    assert resp.status_code == 200, resp.get_json()
    assert resp.get_json()['genome_shared_with'] == [
        {'type': 'lab', 'id': 'labX', 'name': 'Lab X'}]
    assert _lab_shares(world) == [('labX', 1)]
    assert _sees_genome(world, 'carol') and _sees_track(world, 'carol')
    # Only the lab the track went to; and only that track, not its neighbour.
    assert not _sees_genome(world, 'dave')
    assert not _sees_track(world, 'carol', world['track2'])


def test_unsharing_the_track_takes_the_genome_back(world):
    world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    world['call']('POST', f"/api/tracks/{world['track2']}/share", body={'labs': ['labX']})
    world['call']('DELETE', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    # Peaks 2 still needs it.
    assert _lab_shares(world) == [('labX', 1)]
    world['call']('DELETE', f"/api/tracks/{world['track2']}/share", body={'labs': ['labX']})
    assert _lab_shares(world) == []
    assert not _sees_genome(world, 'carol')


def test_a_purposeful_genome_share_outlasts_the_track(world):
    world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    # The owner then shares the genome on purpose: the automatic share becomes theirs.
    world['call']('POST', f"/api/genomes/{world['genome']}/share", body={'labs': ['labX']})
    assert _lab_shares(world) == [('labX', 0)]
    world['call']('DELETE', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    assert _lab_shares(world) == [('labX', 0)]
    assert _sees_genome(world, 'carol')


def test_unsharing_the_genome_keeps_it_for_a_shared_track(world):
    world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    resp = world['call']('DELETE', f"/api/genomes/{world['genome']}/share", body={'labs': ['labX']})
    assert resp.status_code == 200
    assert [k['id'] for k in resp.get_json()['kept']] == ['labX']
    assert _sees_genome(world, 'carol')


def test_a_lab_that_already_sees_the_genome_gets_nothing_extra(world):
    world['call']('POST', f"/api/genomes/{world['genome']}/share", body={'chair_group_id': 'cgX'})
    resp = world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labY']})
    assert resp.get_json()['genome_shared_with'] == []
    assert _lab_shares(world) == []
    assert _sees_track(world, 'dave')


def test_a_track_on_a_suite_genome_needs_no_genome_share(world):
    resp = world['call']('POST', f"/api/tracks/{world['suite_track']}/share", body={'labs': ['labX']})
    assert resp.status_code == 200
    assert resp.get_json()['genome_shared_with'] == []
    assert _sees_genome(world, 'carol', world['suite_genome'])


def test_who_sees_what_is_readable(world):
    world['call']('POST', f"/api/tracks/{world['track']}/share", body={'labs': ['labX']})
    g = world['call']('GET', f"/api/genomes/{world['genome']}/shares").get_json()
    assert g['can_change'] is True
    assert g['owner_lab'] == {'type': 'lab', 'id': 'labA', 'name': 'Lab A'}
    assert g['shares'] == [{'type': 'lab', 'id': 'labX', 'name': 'Lab X', 'auto': True}]
    t = world['call']('GET', f"/api/tracks/{world['track']}/shares").get_json()
    assert t['shares'] == [{'type': 'lab', 'id': 'labX', 'name': 'Lab X'}]
    assert t['owner_chair_group']['id'] == 'cgA'
    # carol may see both now, but not change them.
    assert world['call']('GET', f"/api/genomes/{world['genome']}/shares", user='carol').get_json()['can_change'] is False
    labs = world['call']('GET', '/api/labs').get_json()
    assert [l['id'] for l in labs['labs']] == ['labA', 'labX', 'labY']
    assert 'registration_code' not in labs['labs'][0]
    # dave sees neither.
    assert world['call']('GET', f"/api/tracks/{world['track']}/shares", user='dave').status_code == 404
