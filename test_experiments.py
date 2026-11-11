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


def test_someone_in_no_lab_finds_their_own_seed_records(seeded):
    # The Seed Bank keeps their records under user:<name>, and shows them those.
    assert tags._seed_labs_seen('carol') == {'user:carol'}
    assert 'lab a' in tags._seed_labs_seen('alice')


def test_a_track_names_the_paper_and_geo_data_it_comes_from(seeded):
    c = seeded['client']
    r = c.patch(f"/api/tracks/{seeded['own']}", headers=who('alice'), json={
        'citation': 'Stroud et al. (2013)  Cell', 'doi': 'https://doi.org/10.1016/j.cell.2012.10.054',
        'geo': 'gse39901, gsm980986 GSE39901'})
    assert r.status_code == 200, r.get_json()
    body = r.get_json()
    assert body['citation'] == 'Stroud et al. (2013) Cell'
    assert body['doi'] == '10.1016/j.cell.2012.10.054'          # kept without its address
    assert body['geo'] == 'GSE39901 GSM980986'                  # upper-case, each once
    listed = {t['id']: t for t in c.get(f"/api/genomes/{seeded['genome_id']}/tracks", headers=who('alice')).get_json()}
    assert listed[seeded['own']]['geo'] == 'GSE39901 GSM980986'
    # The notebook's search finds it by its GEO series.
    assert [t['id'] for t in c.get('/api/tracks?q=gse39901', headers=who('alice')).get_json()['tracks']] == [seeded['own']]
    # Only what is sent changes; the experiment is left alone.
    c.patch(f"/api/tracks/{seeded['own']}", headers=who('alice'), json={'experiment': 'Methylomes'})
    again = c.get(f"/api/tracks/{seeded['own']}", headers=who('alice')).get_json()
    assert again['experiment'] == 'Methylomes' and again['doi'] == '10.1016/j.cell.2012.10.054'


def test_a_malformed_doi_or_geo_accession_is_refused(seeded):
    c = seeded['client']
    for body in ({'doi': 'not a doi'}, {'geo': 'GSE39901 SRR12345'}):
        r = c.patch(f"/api/tracks/{seeded['own']}", headers=who('alice'), json=body)
        assert r.status_code == 400 and r.get_json()['error']
    assert c.get(f"/api/tracks/{seeded['own']}", headers=who('alice')).get_json()['geo'] is None
    # An upload with a malformed one is refused before its file is read.
    import io
    r = c.post(f"/api/genomes/{seeded['genome_id']}/tracks", headers=who('alice'), content_type='multipart/form-data',
               data={'name': 'x', 'geo': 'GEO123', 'file': (io.BytesIO(b'Chr1\t1\t2\n'), 'x.bed')})
    assert r.status_code == 400 and 'GEO' in r.get_json()['error']


def test_the_owner_lab_removes_a_track_and_everything_of_it(seeded):
    c, root = seeded['client'], seeded['data_root']
    folder = root / 'tracks' / seeded['own']
    folder.mkdir(parents=True)
    (folder / 'cov.bw').write_bytes(b'x')
    conn = db.connect(root / 'genome.db')
    try:
        conn.execute("UPDATE tracks SET file_path = ? WHERE id = ?", (str(folder / 'cov.bw'), seeded['own']))
        conn.execute("INSERT INTO seed_stock_tags (id, target_type, target_id, stock_kind, stock_id, created_by) "
                     "VALUES ('t1', 'track', ?, 'line', 'GL-L-0001', 'alice')", (seeded['own'],))
        conn.execute("INSERT INTO track_shares (track_id, target_type, target_id, granted_by) "
                     "VALUES (?, 'lab', 'labB', 'alice')", (seeded['own'],))
        conn.commit()
    finally:
        conn.close()
    # Bob sees it now (shared with his lab) but may not remove it.
    assert c.delete(f"/api/tracks/{seeded['own']}", headers=who('bob')).status_code == 403
    r = c.delete(f"/api/tracks/{seeded['own']}", headers=who('alice'))
    assert r.status_code == 200 and r.get_json()['removed'] == seeded['own']
    assert not folder.exists()
    conn = db.connect(root / 'genome.db')
    try:
        for sql in ('SELECT 1 FROM tracks WHERE id = ?', 'SELECT 1 FROM track_shares WHERE track_id = ?',
                    "SELECT 1 FROM seed_stock_tags WHERE target_id = ?"):
            assert conn.execute(sql, (seeded['own'],)).fetchone() is None, sql
    finally:
        conn.close()
    assert c.get(f"/api/tracks/{seeded['own']}", headers=who('alice')).status_code == 404
    assert c.delete(f"/api/tracks/{seeded['own']}", headers=who('alice')).status_code == 404


def test_removing_never_reaches_outside_the_data_root(tmp_path):
    outside = tmp_path / 'elsewhere.bw'
    outside.write_bytes(b'x')
    root = tmp_path / 'data'
    (root / 'tracks').mkdir(parents=True)
    tracks.remove_track_files(root, '../..', str(outside), [str(tmp_path / 'x' / 'in.bam')])
    assert outside.exists() and root.exists()


def test_the_seed_bank_asks_what_is_tagged_with_one_of_its_records(seeded):
    c = seeded['client']
    conn = db.connect(seeded['data_root'] / 'genome.db')
    try:
        for tid, target_type, target in (('t1', 'track', seeded['own']), ('t2', 'track', seeded['suite']),
                                         ('t3', 'genome', seeded['genome_id'])):
            conn.execute("INSERT INTO seed_stock_tags (id, target_type, target_id, stock_kind, stock_id, created_by) "
                         "VALUES (?, ?, ?, 'line', 'GL-L-0001', 'alice')", (tid, target_type, target))
        conn.commit()
    finally:
        conn.close()
    seen = c.get('/api/seed-stocks/GL-L-0001/tagged', headers=who('alice')).get_json()
    assert {t['id'] for t in seen['tracks']} == {seeded['own'], seeded['suite']}
    assert seen['tracks'][0]['genome'] == 'G' and [g['id'] for g in seen['genomes']] == [seeded['genome_id']]
    # Bob sees only the track shared with everyone (and the genome everyone sees).
    bob = c.get('/api/seed-stocks/GL-L-0001/tagged', headers=who('bob')).get_json()
    assert [t['id'] for t in bob['tracks']] == [seeded['suite']]
    assert c.get('/api/seed-stocks/GL-L-9999/tagged', headers=who('alice')).get_json() == {'tracks': [], 'genomes': []}
