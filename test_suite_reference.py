"""
Suite reference genomes: an administrator's genome that belongs to the whole
suite, not a lab, and that every signed-in member sees — with the tracks an
administrator adds to it.
"""

import io

import pytest

import app as app_module
import identity
import labs_reader

FASTA = b'>Chr1\nACGTACGTAC\n>ChrM\nACGT\n'


@pytest.fixture
def suite(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True

    # 'admin' is an administrator in no lab; alice is in labA, bob in labB.
    labs = {'alice': 'labA', 'bob': 'labB'}
    monkeypatch.setattr(identity, 'is_admin', lambda u: u == 'admin')
    monkeypatch.setattr(
        labs_reader, 'labs_of',
        lambda u: [{'lab_id': labs[u], 'lab': labs[u], 'role': 'member'}]
        if u in labs else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    client = application.test_client()

    def as_user(user):
        return {'X-Suite-User': user}

    def add_genome(user, suite_reference='true', assembly='TAIR10'):
        return client.post('/api/genomes', data={
            'name': 'Arabidopsis thaliana (TAIR10 / Araport11)',
            'assembly_id': assembly,
            'suite_reference': suite_reference,
            'file': (io.BytesIO(FASTA), 'TAIR10_chr_all.fas'),
        }, headers=as_user(user), content_type='multipart/form-data')

    def add_track(user, genome_id, body=b'Chr1\t.\tgene\t1\t8\t.\t+\t.\tID=g1\n'):
        return client.post(f'/api/genomes/{genome_id}/tracks', data={
            'name': 'Araport11 genes',
            'file': (io.BytesIO(body), 'genes.gff3'),
        }, headers=as_user(user), content_type='multipart/form-data')

    return {'client': client, 'as': as_user, 'add_genome': add_genome,
            'add_track': add_track}


def test_admin_in_no_lab_adds_a_suite_reference_every_member_sees(suite):
    resp = suite['add_genome']('admin')
    assert resp.status_code == 201, resp.get_json()
    body = resp.get_json()
    assert body['owner_lab'] == 'suite'
    assert body['genome_visibility'] == 'suite'

    for user in ('alice', 'bob', 'admin'):
        listed = suite['client'].get('/api/genomes',
                                     headers=suite['as'](user)).get_json()
        assert [g['id'] for g in listed] == [body['id']], user
        detail = suite['client'].get(f"/api/genomes/{body['id']}",
                                     headers=suite['as'](user))
        assert detail.status_code == 200


def test_only_an_administrator_adds_a_suite_reference(suite):
    resp = suite['add_genome']('alice')
    assert resp.status_code == 403
    listed = suite['client'].get('/api/genomes',
                                 headers=suite['as']('alice')).get_json()
    assert listed == []


def test_without_the_flag_an_upload_is_still_the_lab_s(suite):
    resp = suite['add_genome']('alice', suite_reference='')
    assert resp.status_code == 201
    assert resp.get_json()['owner_lab'] == 'labA'
    assert suite['client'].get('/api/genomes',
                               headers=suite['as']('bob')).get_json() == []


def test_admin_track_on_a_suite_reference_is_seen_by_every_member(suite):
    genome_id = suite['add_genome']('admin').get_json()['id']
    resp = suite['add_track']('admin', genome_id)
    assert resp.status_code == 201, resp.get_json()
    track = resp.get_json()
    assert track['owner_lab'] == 'suite'
    assert track['visibility'] == 'suite'
    for user in ('alice', 'bob'):
        library = suite['client'].get(f'/api/genomes/{genome_id}/tracks',
                                      headers=suite['as'](user)).get_json()
        assert [t['id'] for t in library] == [track['id']], user
        region = suite['client'].get(
            f"/api/tracks/{track['id']}/region?chrom=Chr1&start=0&end=10",
            headers=suite['as'](user))
        assert region.status_code == 200, user


def test_a_member_s_track_on_a_suite_reference_stays_their_lab_s(suite):
    genome_id = suite['add_genome']('admin').get_json()['id']
    resp = suite['add_track']('alice', genome_id)
    assert resp.status_code == 201
    track = resp.get_json()
    assert track['owner_lab'] == 'labA'
    assert track['visibility'] == 'lab_only'
    bob_sees = suite['client'].get(f'/api/genomes/{genome_id}/tracks',
                                   headers=suite['as']('bob')).get_json()
    assert bob_sees == []


def test_suite_reference_cannot_be_reshared_by_a_member(suite):
    genome_id = suite['add_genome']('admin').get_json()['id']
    resp = suite['client'].post(f'/api/genomes/{genome_id}/share',
                                json={'chair_group': 'g1'},
                                headers=suite['as']('alice'))
    assert resp.status_code in (400, 403)
