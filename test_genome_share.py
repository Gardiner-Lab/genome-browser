"""
Tests for the genome sharing endpoints (task 5.3).

These exercise ``POST`` and ``DELETE /api/genomes/<id>/share`` end to end
against a seeded temporary database: an Owner_Lab member may share a genome
with a whole chair group (which then lets that chair group's labs view it), a
non-owner is refused with a 403 and nothing changed, an individual-lab target
and a non-existent chair group are both rejected, a revoke removes the share,
and every mutation touches only the one genome. Labs and chair groups are
monkeypatched on ``labs_reader`` so the tests run without a ``suite.db``.
"""

import uuid

import pytest

import access
import app as app_module
import db
import labs_reader


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def seeded(tmp_path, monkeypatch):
    """A built app + client over a temp data root with two genomes.

    genome_a is owned by labA (alice's lab); genome_b is owned by labA too but
    exists only to prove a mutation touches one genome. The lab/chair-group
    graph: labA is in chair group cgOwner; labX is in chair group cgOther
    (the share target); carol is in labX only; bob is in labZ (no shared
    group), so bob cannot change sharing and cannot view until shared.
    """
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True

    genome_a = uuid.uuid4().hex
    genome_b = uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.executemany(
            'INSERT INTO reference_genomes '
            '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
            'VALUES (?, ?, ?, ?, ?, ?)',
            [(genome_a, 'Genome A', 'asmA', 'labA', 'owner_lab_only', 'alice'),
             (genome_b, 'Genome B', 'asmB', 'labA', 'owner_lab_only', 'alice')])
        conn.commit()
    finally:
        conn.close()

    members = {
        'alice': [{'lab_id': 'labA', 'lab': 'Lab A', 'role': 'member'}],
        'carol': [{'lab_id': 'labX', 'lab': 'Lab X', 'role': 'member'}],
        'bob': [{'lab_id': 'labZ', 'lab': 'Lab Z', 'role': 'member'}],
    }
    lab_group = {'labA': 'cgOwner', 'labX': 'cgOther', 'labZ': None}
    group_labs = {'cgOwner': ['labA'], 'cgOther': ['labX']}

    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: members.get(u, []))
    monkeypatch.setattr(labs_reader, 'chair_group_of',
                        lambda lab: lab_group.get(lab))
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group',
                        lambda g: group_labs.get(g, []))
    monkeypatch.setattr(labs_reader, 'chair_group_exists',
                        lambda g: g in group_labs)

    client = application.test_client()

    def share(genome_id, body, user='alice', method='POST'):
        return client.open(
            f'/api/genomes/{genome_id}/share', method=method, json=body,
            headers={'X-Suite-User': user})

    return {'client': client, 'data_root': data_root, 'share': share,
            'genome_a': genome_a, 'genome_b': genome_b}


def _shares(data_root, genome_id=None):
    conn = db.connect(data_root / 'genome.db')
    try:
        if genome_id is None:
            return conn.execute('SELECT * FROM genome_shares').fetchall()
        return conn.execute(
            'SELECT * FROM genome_shares WHERE genome_id = ?',
            (genome_id,)).fetchall()
    finally:
        conn.close()


def _can_view(data_root, member, genome_id):
    """Whether ``member`` may view the genome, resolved over the temp DB.

    Uses ``access.can_view_genome`` (the single visibility resolver the share
    endpoint feeds) against the test's own database, so the check reflects the
    effect of a share without depending on an unrelated read endpoint.
    """
    conn = db.connect(data_root / 'genome.db')
    try:
        genome = conn.execute(
            'SELECT id, name, assembly_id, owner_lab, genome_visibility '
            'FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()
        return access.can_view_genome(member, genome, conn)
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Owner shares a chair group -> row added, that chair group's labs can view
# (Req 12.3)
# ---------------------------------------------------------------------------


def test_owner_shares_chair_group(seeded):
    resp = seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    assert resp.status_code == 200
    assert resp.get_json()['shared'] is True

    rows = _shares(seeded['data_root'], seeded['genome_a'])
    assert len(rows) == 1
    assert rows[0]['chair_group_id'] == 'cgOther'
    assert rows[0]['granted_by'] == 'alice'

    # carol (labX in cgOther) can now view the genome; before sharing she
    # could not (owner_lab_only default, and labX is not labA).
    assert _can_view(seeded['data_root'], 'carol', seeded['genome_a'])


def test_share_is_idempotent(seeded):
    seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    resp = seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    assert resp.status_code == 200
    assert len(_shares(seeded['data_root'], seeded['genome_a'])) == 1


def test_share_accepts_target_type_pair(seeded):
    resp = seeded['share'](
        seeded['genome_a'],
        {'target_type': 'chair_group', 'target_id': 'cgOther'})
    assert resp.status_code == 200
    assert len(_shares(seeded['data_root'], seeded['genome_a'])) == 1


# ---------------------------------------------------------------------------
# Non-owner -> 403, nothing changed (Req 12.6)
# ---------------------------------------------------------------------------


def test_non_owner_forbidden(seeded):
    resp = seeded['share'](
        seeded['genome_a'], {'chair_group_id': 'cgOther'}, user='carol')
    # carol can view genome_a only after it is shared; before sharing she gets
    # a 404 (not visible). Share it first as the owner, then carol (a viewer
    # but not an owner-lab member) must still be refused with 403.
    assert resp.status_code == 404

    seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    resp = seeded['share'](
        seeded['genome_a'], {'chair_group_id': 'cgOwner'}, user='carol')
    assert resp.status_code == 403
    assert 'permitted' in resp.get_json()['error'].lower()

    # The owner's single share is untouched; carol added nothing.
    rows = _shares(seeded['data_root'], seeded['genome_a'])
    assert [r['chair_group_id'] for r in rows] == ['cgOther']


# ---------------------------------------------------------------------------
# Individual-lab target rejected (Req 12.4)
# ---------------------------------------------------------------------------


def test_individual_lab_target_rejected(seeded):
    resp = seeded['share'](
        seeded['genome_a'], {'target_type': 'lab', 'target_id': 'labX'})
    assert resp.status_code == 400
    assert 'chair group' in resp.get_json()['error'].lower()
    assert _shares(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Non-existent chair group rejected (Req 12.8)
# ---------------------------------------------------------------------------


def test_nonexistent_chair_group_rejected(seeded):
    resp = seeded['share'](
        seeded['genome_a'], {'chair_group_id': 'cgGhost'})
    assert resp.status_code == 400
    assert 'invalid' in resp.get_json()['error'].lower()
    assert _shares(seeded['data_root']) == []


def test_missing_target_rejected(seeded):
    resp = seeded['share'](seeded['genome_a'], {})
    assert resp.status_code == 400
    assert _shares(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Revoke removes the share; default viewers keep it (Req 12.7)
# ---------------------------------------------------------------------------


def test_revoke_removes_share(seeded):
    seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    assert len(_shares(seeded['data_root'], seeded['genome_a'])) == 1

    resp = seeded['share'](
        seeded['genome_a'], {'chair_group_id': 'cgOther'}, method='DELETE')
    assert resp.status_code == 200
    assert resp.get_json()['shared'] is False
    assert _shares(seeded['data_root'], seeded['genome_a']) == []

    # carol (cgOther) loses view immediately; alice (owner lab) keeps it.
    assert not _can_view(seeded['data_root'], 'carol', seeded['genome_a'])
    assert _can_view(seeded['data_root'], 'alice', seeded['genome_a'])


def test_revoke_unshared_is_harmless(seeded):
    resp = seeded['share'](
        seeded['genome_a'], {'chair_group_id': 'cgOther'}, method='DELETE')
    assert resp.status_code == 200
    assert _shares(seeded['data_root']) == []


# ---------------------------------------------------------------------------
# Locality: only the named genome changes (Req 12.5)
# ---------------------------------------------------------------------------


def test_only_named_genome_changes(seeded):
    seeded['share'](seeded['genome_a'], {'chair_group_id': 'cgOther'})
    # genome_b, owned by the same lab, is untouched.
    assert _shares(seeded['data_root'], seeded['genome_b']) == []
    all_rows = _shares(seeded['data_root'])
    assert {r['genome_id'] for r in all_rows} == {seeded['genome_a']}


# ---------------------------------------------------------------------------
# Missing / unviewable genome -> 404; unauthenticated -> 401
# ---------------------------------------------------------------------------


def test_unknown_genome_404(seeded):
    resp = seeded['share']('does-not-exist', {'chair_group_id': 'cgOther'})
    assert resp.status_code == 404


def test_unauthenticated_401(seeded):
    resp = seeded['client'].post(
        f"/api/genomes/{seeded['genome_a']}/share",
        json={'chair_group_id': 'cgOther'})
    assert resp.status_code == 401
