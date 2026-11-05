"""
Tests for plasmid and seed-stock tagging (task 14.1).

These exercise ``POST/GET/DELETE /api/tracks/<id>/tags`` and
``/api/genomes/<id>/tags`` end to end against a seeded temporary database, with
the external-tool validators (``plasmid_exists`` / ``seed_stock_lookup``)
monkeypatched so the tests run without a real MoClo/Editor/Seed Bank store —
exercising the endpoint logic (view guard, missing-id rejection naming the id,
duplicate rejection retaining the existing row, listing with links, delete).

``register_tags`` is wired onto the test app here because the app factory does
not call it yet (the orchestrator adds that line); the registrar is a single
call, mirroring how the module will be wired.
"""

import uuid

import pytest

import app as app_module
import db
import labs_reader
import tags


@pytest.fixture
def seeded(tmp_path, monkeypatch):
    """A built app + client with one viewable genome and one viewable track."""
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root
        MOCLO_URL = '/'
        EDITOR_URL = '/editor'
        SEEDS_URL = '/seeds'

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    # The app factory wires register_tags, so the endpoints are already present.

    genome_id = uuid.uuid4().hex
    track_id = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO reference_genomes '
            '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            (genome_id, 'Test Genome', 'asm1', 'labA', 'owner_lab_only', 'alice'))
        conn.execute(
            'INSERT INTO tracks '
            '(id, genome_id, name, type, owner_lab, visibility, track_status, '
            ' created_by) VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
            (track_id, genome_id, 'Coverage', 'bigwig', 'labA', 'lab_only',
             'ready', 'alice'))
        conn.commit()
    finally:
        conn.close()

    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A',
                                    'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])

    client = application.test_client()

    def post(kind, external, body_key, target='track', target_id=None):
        tid = target_id or (track_id if target == 'track' else genome_id)
        return client.post(
            f'/api/{target}s/{tid}/tags',
            json={'kind': kind, body_key: external},
            headers={'X-Suite-User': 'alice'})

    return {'client': client, 'genome_id': genome_id, 'track_id': track_id,
            'data_root': data_root, 'post': post}


def _plasmid_rows(data_root):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute('SELECT * FROM plasmid_tags').fetchall()
    finally:
        conn.close()


def _stock_rows(data_root):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute('SELECT * FROM seed_stock_tags').fetchall()
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Add a valid plasmid tag (Req 11.1)
# ---------------------------------------------------------------------------

def test_add_plasmid_tag_records_row_with_link(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'plasmid_exists', lambda pid: pid == 'pGFP')
    resp = seeded['post']('plasmid', 'pGFP', 'plasmid_id')
    assert resp.status_code == 201
    body = resp.get_json()
    assert body['plasmid_id'] == 'pGFP'
    assert body['kind'] == 'plasmid'
    assert 'pGFP' in body['link']  # navigable MoClo/Editor link (11.6)
    rows = _plasmid_rows(seeded['data_root'])
    assert len(rows) == 1
    assert rows[0]['plasmid_id'] == 'pGFP'


# ---------------------------------------------------------------------------
# Missing id rejected naming it, target unchanged (Req 11.3)
# ---------------------------------------------------------------------------

def test_missing_plasmid_rejected_naming_it(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'plasmid_exists', lambda pid: False)
    resp = seeded['post']('plasmid', 'pGhost', 'plasmid_id')
    assert resp.status_code == 404
    assert 'pGhost' in resp.get_json()['error']
    assert _plasmid_rows(seeded['data_root']) == []


def test_missing_seed_stock_rejected_naming_it(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'seed_stock_lookup', lambda sid: (False, None))
    resp = seeded['post']('seed_stock', 'L-ghost', 'stock_id')
    assert resp.status_code == 404
    assert 'L-ghost' in resp.get_json()['error']
    assert _stock_rows(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Duplicate rejected, existing retained (Req 11.4)
# ---------------------------------------------------------------------------

def test_duplicate_plasmid_rejected_retains_existing(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'plasmid_exists', lambda pid: True)
    first = seeded['post']('plasmid', 'pGFP', 'plasmid_id')
    assert first.status_code == 201
    second = seeded['post']('plasmid', 'pGFP', 'plasmid_id')
    assert second.status_code == 409
    # Exactly one row retained.
    assert len(_plasmid_rows(seeded['data_root'])) == 1


# ---------------------------------------------------------------------------
# Seed line tag records the Family_Tree key (Req 11.5) + list (11.6)
# ---------------------------------------------------------------------------

def test_seed_line_tag_records_kind_and_lists_with_link(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'seed_stock_lookup', lambda sid: (True, 'line'))
    resp = seeded['post']('seed_stock', 'Line-42', 'stock_id')
    assert resp.status_code == 201
    body = resp.get_json()
    assert body['stock_kind'] == 'line'
    assert body['stock_id'] == 'Line-42'

    rows = _stock_rows(seeded['data_root'])
    assert len(rows) == 1 and rows[0]['stock_id'] == 'Line-42'

    listed = seeded['client'].get(
        f"/api/tracks/{seeded['track_id']}/tags",
        headers={'X-Suite-User': 'alice'}).get_json()
    assert len(listed['seed_stock_tags']) == 1
    tag = listed['seed_stock_tags'][0]
    assert 'Line-42' in tag['link'] and 'lines' in tag['link']


# ---------------------------------------------------------------------------
# Delete removes the row (Req 11.7, 11.8)
# ---------------------------------------------------------------------------

def test_delete_removes_tag(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'seed_stock_lookup', lambda sid: (True, 'plant'))
    created = seeded['post']('seed_stock', 'Plant-7', 'stock_id').get_json()
    tag_id = created['id']
    resp = seeded['client'].delete(
        f"/api/tracks/{seeded['track_id']}/tags/{tag_id}",
        headers={'X-Suite-User': 'alice'})
    assert resp.status_code == 204
    assert _stock_rows(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Tagging on a genome target works too (Req 11.1, same endpoints)
# ---------------------------------------------------------------------------

def test_tag_on_genome_target(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'plasmid_exists', lambda pid: True)
    resp = seeded['post']('plasmid', 'pRED', 'plasmid_id', target='genome')
    assert resp.status_code == 201
    rows = _plasmid_rows(seeded['data_root'])
    assert len(rows) == 1 and rows[0]['target_type'] == 'genome'


# ---------------------------------------------------------------------------
# Unviewable target is a 404 (view guard)
# ---------------------------------------------------------------------------

def test_tag_on_unviewable_target_404(seeded, monkeypatch):
    monkeypatch.setattr(tags, 'plasmid_exists', lambda pid: True)
    # bob is in no lab, so the lab_only track is not viewable.
    resp = seeded['client'].post(
        f"/api/tracks/{seeded['track_id']}/tags",
        json={'kind': 'plasmid', 'plasmid_id': 'pGFP'},
        headers={'X-Suite-User': 'bob'})
    assert resp.status_code == 404
    assert _plasmid_rows(seeded['data_root']) == []
