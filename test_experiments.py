"""
A track's experiment (to sort and group the track list by), and finding Seed
Bank lines, lots and plants to tag a track with.
"""

import sqlite3
import uuid

import pytest

import app as app_module
import db
import identity
import labs_reader
import tags
import tracks


@pytest.fixture
def seeded(tmp_path, monkeypatch):
    """Alice is in Lab A, which owns two tracks: one its own, one shared with everyone. Bob is in Lab B."""
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root
        SEEDS_URL = '/seeds'

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    genome_id, own, suite = uuid.uuid4().hex, uuid.uuid4().hex, uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                     "VALUES (?, 'G', 'asm', 'suite', 'suite', 'alice')", (genome_id,))
        for tid, name, visibility in ((own, 'Rep1 coverage', 'lab_only'), (suite, 'Peaks', 'suite')):
            conn.execute('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, created_by) '
                         "VALUES (?, ?, ?, 'bigwig', 'labA', ?, 'ready', 'alice')", (tid, genome_id, name, visibility))
        conn.commit()
    finally:
        conn.close()
    labs = {'alice': [{'lab_id': 'labA', 'lab': 'Lab A', 'role': 'member'}],
            'bob': [{'lab_id': 'labB', 'lab': 'Lab B', 'role': 'member'}]}
    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: labs.get(u, []))
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: 'plant' if lab in ('labA', 'labS') else None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: ['labA', 'labS'] if g == 'plant' else [])
    monkeypatch.setattr(labs_reader, 'lab_name', lambda lab: {'labA': 'Lab A', 'labS': 'Smith Lab', 'labB': 'Lab B'}.get(lab))
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    return {'client': application.test_client(), 'app': application, 'genome_id': genome_id,
            'own': own, 'suite': suite, 'data_root': data_root}


def who(name):
    return {'X-Suite-User': name}


def test_the_owner_lab_sets_a_tracks_experiment(seeded):
    c = seeded['client']
    r = c.patch(f"/api/tracks/{seeded['own']}", json={'experiment': '  Heat   stress RNA-seq '}, headers=who('alice'))
    assert r.status_code == 200, r.get_json()
    assert r.get_json()['experiment'] == 'Heat stress RNA-seq'
    listed = {t['id']: t for t in c.get(f"/api/genomes/{seeded['genome_id']}/tracks", headers=who('alice')).get_json()}
    assert listed[seeded['own']]['experiment'] == 'Heat stress RNA-seq'
    assert listed[seeded['own']]['can_change'] is True and listed[seeded['own']]['created_at']
    # Blank takes it out of any experiment.
    assert c.patch(f"/api/tracks/{seeded['own']}", json={'experiment': ' '}, headers=who('alice')).get_json()['experiment'] is None


def test_only_the_owner_lab_changes_it(seeded):
    c = seeded['client']
    # Bob sees the track shared with everyone, but it is not his lab's.
    assert c.patch(f"/api/tracks/{seeded['suite']}", json={'experiment': 'Mine'}, headers=who('bob')).status_code == 403
    listed = {t['id']: t for t in c.get(f"/api/genomes/{seeded['genome_id']}/tracks", headers=who('bob')).get_json()}
    assert listed[seeded['suite']]['can_change'] is False
    # Lab A's own track he does not see at all.
    assert c.patch(f"/api/tracks/{seeded['own']}", json={'experiment': 'Mine'}, headers=who('bob')).status_code == 404
    assert c.patch(f"/api/tracks/{seeded['own']}", json={}, headers=who('alice')).status_code == 400


def test_an_upload_names_its_experiment(seeded):
    with seeded['app'].test_request_context(method='POST', data={'experiment': 'Drought 2026'}):
        conn = db.connect(seeded['data_root'] / 'genome.db')
        try:
            tracks._insert_track(conn, 'new1', seeded['genome_id'], 'Drought rep1', 'bigwig', 'labA',
                                 'alice', 'ready', None)
            conn.commit()
            assert conn.execute("SELECT experiment FROM tracks WHERE id = 'new1'").fetchone()[0] == 'Drought 2026'
        finally:
            conn.close()
    assert tracks.clean_experiment('x' * 300) == 'x' * tracks.EXPERIMENT_MAX
    assert tracks.clean_experiment(None) is None


def test_the_seed_bank_search_shows_only_what_the_member_sees(seeded, tmp_path, monkeypatch):
    seeds = tmp_path / 'seeds.db'
    conn = sqlite3.connect(seeds)
    conn.executescript("""
        CREATE TABLE lines (id TEXT, code TEXT, lab TEXT, name TEXT, deleted_at TEXT);
        CREATE TABLE lots (id TEXT, code TEXT, lab TEXT, line_id TEXT, deleted_at TEXT);
        CREATE TABLE plants (id TEXT, code TEXT, lab TEXT, lot_id TEXT, line_id TEXT, number INTEGER, deleted_at TEXT);
        INSERT INTO lines VALUES ('l1', 'GL-L-0001', 'Lab A', 'Col-0 pUBQ:GFP', NULL),
                                 ('l2', 'SL-L-0001', 'Smith Lab', 'Col-0 smith', NULL),
                                 ('l3', 'BL-L-0001', 'Lab B', 'Col-0 bob', NULL),
                                 ('l4', 'GL-L-0002', 'Lab A', 'Col-0 gone', '2026-01-01');
        INSERT INTO lots VALUES ('s1', 'GL-S-0001', 'Lab A', 'l1', NULL);
        INSERT INTO plants VALUES ('p1', 'GL-P-0001', 'Lab A', 's1', 'l1', 3, NULL);
    """)
    conn.commit()
    conn.close()
    monkeypatch.setenv('SEEDS_DB', str(seeds))
    found = seeded['client'].get('/api/seed-stocks?q=Col-0', headers=who('alice')).get_json()['stocks']
    codes = {s['code']: s for s in found}
    # Lab A's and its chair group mate Smith Lab's; not Lab B's, nor a deleted line.
    assert set(codes) == {'GL-L-0001', 'SL-L-0001', 'GL-S-0001', 'GL-P-0001'}
    assert codes['GL-P-0001']['kind'] == 'plant' and '#3' in codes['GL-P-0001']['label']
    assert {s['code'] for s in seeded['client'].get('/api/seed-stocks?q=Col-0', headers=who('bob')).get_json()['stocks']} \
        == {'BL-L-0001'}


def test_a_tags_seed_bank_link_is_the_seed_banks_own_address(seeded):
    with seeded['app'].app_context():
        assert tags._seed_stock_link('plant', 'GL-P-0001') == '/seeds/#plant=GL-P-0001'
        assert tags._seed_stock_link('lot', 'GL-S-0001') == '/seeds/#lot=GL-S-0001'


def test_the_notebook_finds_tracks_across_genomes(seeded):
    c = seeded['client']
    c.patch(f"/api/tracks/{seeded['own']}", json={'experiment': 'Heat stress'}, headers=who('alice'))
    found = c.get('/api/tracks?q=heat', headers=who('alice')).get_json()['tracks']
    assert [t['id'] for t in found] == [seeded['own']]
    assert found[0]['genome'] == 'G' and found[0]['genome_id'] == seeded['genome_id']
    # Every word must match: the track's name, experiment or genome.
    assert {t['id'] for t in c.get('/api/tracks?q=g', headers=who('alice')).get_json()['tracks']} \
        == {seeded['own'], seeded['suite']}
    # Bob finds only what he sees: the track shared with everyone.
    assert [t['id'] for t in c.get('/api/tracks', headers=who('bob')).get_json()['tracks']] == [seeded['suite']]


def test_a_tracks_metadata_names_its_genome(seeded):
    # What a notebook link (?track=) opens the track's genome by.
    body = seeded['client'].get(f"/api/tracks/{seeded['own']}", headers=who('alice')).get_json()
    assert body['genome_id'] == seeded['genome_id']
