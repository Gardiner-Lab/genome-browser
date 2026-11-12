"""
The files the viewer reads: an indexed reference and indexed track data, served
by range request to members who may see them, and to no one else.

These need samtools, bgzip and tabix, which are in the image; where they are
not installed (a bare test runner) the tests skip.
"""

import gzip
import io
import shutil
import subprocess

import pytest

import app as app_module
import identity
import labs_reader
import serving

NEEDS_TOOLS = pytest.mark.skipif(
    not all(shutil.which(t) for t in ('samtools', 'bgzip', 'tabix', 'sort')),
    reason='samtools/bgzip/tabix are in the image, not on this runner')

# Uneven line lengths within Chr1: samtools cannot index it as it is.
FASTA = b'>Chr1 first\nACGTACGTAC\nACG\nTTTTT\n>ChrM\nGGGGCCCCAA\n'
GFF = (b'##gff-version 3\n'
       b'Chr1\t.\tgene\t9\t12\t.\t+\t.\tID=g2;Note=a\x1eb\n'
       b'Chr1\t.\tgene\t1\t5\t.\t+\t.\tID=g1\n'
       b'###\n'
       b'ChrM\t.\tgene\t2\t4\t.\t-\t.\tID=m1\n')


@pytest.fixture
def served(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    labs = {'alice': 'labA', 'bob': 'labB'}
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    monkeypatch.setattr(
        labs_reader, 'labs_of',
        lambda u: [{'lab_id': labs[u], 'lab': labs[u], 'role': 'member'}]
        if u in labs else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    client = application.test_client()

    genome = client.post('/api/genomes', data={
        'name': 'Test', 'assembly_id': 'T1',
        'file': (io.BytesIO(FASTA), 'ref.fa')},
        headers={'X-Suite-User': 'alice'},
        content_type='multipart/form-data').get_json()
    track = client.post(f"/api/genomes/{genome['id']}/tracks", data={
        'name': 'Genes', 'file': (io.BytesIO(GFF), 'genes.gff3')},
        headers={'X-Suite-User': 'alice'},
        content_type='multipart/form-data').get_json()
    return {'client': client, 'genome': genome['id'], 'track': track['id'],
            'data_root': data_root}


def _get(served, path, user='alice', **headers):
    headers['X-Suite-User'] = user
    return served['client'].get(path, headers=headers)


@NEEDS_TOOLS
def test_reference_is_indexed_and_read_by_range(served):
    fai = _get(served, f"/api/genomes/{served['genome']}/reference.fa.fai")
    assert fai.status_code == 200
    rows = [line.split('\t') for line in fai.get_data(as_text=True).splitlines()]
    assert [(r[0], int(r[1])) for r in rows] == [('Chr1', 18), ('ChrM', 10)]

    # The uneven FASTA was re-wrapped so it could be indexed; the sequence is
    # the same and a range of it reads back exactly.
    offset, = (int(r[2]) for r in rows if r[0] == 'Chr1')
    part = _get(served, f"/api/genomes/{served['genome']}/reference.fa",
                Range=f'bytes={offset}-{offset + 17}')
    assert part.status_code == 206
    assert part.get_data(as_text=True) == 'ACGTACGTACACGTTTTT'


@NEEDS_TOOLS
def test_feature_track_is_sorted_compressed_and_tabix_indexed(served, tmp_path):
    base = f"/api/tracks/{served['track']}/features.gff3.gz"
    data = _get(served, base)
    index = _get(served, base + '.tbi')
    assert data.status_code == 200 and index.status_code == 200
    text = gzip.decompress(data.get_data()).decode('utf-8')
    lines = text.splitlines()
    assert lines[0] == '##gff-version 3'
    assert '###' not in lines
    body = [l.split('\t')[3] for l in lines if not l.startswith('#')]
    assert body == ['1', '9', '2']  # Chr1 sorted by start, then ChrM
    assert 'Note=a b' in text       # the separator stayed in its line, as a space

    # tabix can query what was served.
    (tmp_path / 'f.gff3.gz').write_bytes(data.get_data())
    (tmp_path / 'f.gff3.gz.tbi').write_bytes(index.get_data())
    hits = subprocess.run(['tabix', str(tmp_path / 'f.gff3.gz'), 'Chr1:8-20'],
                          capture_output=True, text=True, check=True).stdout
    assert 'ID=g2' in hits and 'ID=g1' not in hits


@NEEDS_TOOLS
def test_wrong_extension_and_other_labs_get_nothing(served):
    assert _get(served, f"/api/tracks/{served['track']}/features.bed.gz").status_code == 404
    assert _get(served, f"/api/tracks/{served['track']}/coverage.bw").status_code == 404
    for path in (f"/api/genomes/{served['genome']}/reference.fa",
                 f"/api/genomes/{served['genome']}/reference.fa.fai",
                 f"/api/tracks/{served['track']}/features.gff3.gz",
                 f"/api/tracks/{served['track']}/features.gff3.gz.tbi"):
        assert _get(served, path, user='bob').status_code == 404, path


@NEEDS_TOOLS
def test_a_region_comes_with_its_sequence_and_annotations(served):
    base = f"/api/genomes/{served['genome']}/sequence"
    # Chr1 is ACGTACGTAC ACG TTTTT; bases 0..12 hold the two genes (g1 1-5, g2 9-12).
    plus = _get(served, f"{base}?chrom=Chr1&start=0&end=12&tracks={served['track']}")
    assert plus.status_code == 200, plus.get_json()
    body = plus.get_json()
    assert body['sequence'] == 'ACGTACGTACAC'
    feats = sorted((f['name'], f['start'], f['end'], f['strand']) for f in body['features'])
    assert feats == [('g1', 0, 5, 1), ('g2', 8, 12, 1)]

    # The minus strand: the reverse complement, and the features turned with it.
    minus = _get(served, f"{base}?chrom=Chr1&start=0&end=12&strand=-1&tracks={served['track']}").get_json()
    assert minus['sequence'] == 'GTGTACGTACGT'
    feats = sorted((f['name'], f['start'], f['end'], f['strand']) for f in minus['features'])
    assert feats == [('g1', 7, 12, -1), ('g2', 0, 4, -1)]


@NEEDS_TOOLS
def test_a_region_is_checked(served):
    base = f"/api/genomes/{served['genome']}/sequence"
    assert _get(served, f'{base}?chrom=Chr9&start=0&end=5').status_code == 400
    assert _get(served, f'{base}?chrom=Chr1&start=5&end=5').status_code == 400
    assert _get(served, f'{base}?chrom=Chr1&start=x&end=5').status_code == 400
    # Clamped to the chromosome, and nothing for another lab.
    assert _get(served, f'{base}?chrom=ChrM&start=0&end=999').get_json()['sequence'] == 'GGGGCCCCAA'
    assert _get(served, f'{base}?chrom=ChrM&start=0&end=5', user='bob').status_code == 404


@NEEDS_TOOLS
def test_a_region_downloads_as_fasta_of_any_length(served):
    base = f"/api/genomes/{served['genome']}/sequence.fasta"
    plus = _get(served, f'{base}?chrom=Chr1&start=0&end=12&name=AT1G01010.1')
    assert plus.status_code == 200
    assert 'attachment' in plus.headers['Content-Disposition']
    assert plus.get_data(as_text=True) == '>AT1G01010.1 T1 Chr1:1-12\nACGTACGTACAC\n'
    minus = _get(served, f'{base}?chrom=Chr1&start=0&end=12&strand=-1').get_data(as_text=True)
    assert minus == '>T1 Chr1:1-12(-)\nGTGTACGTACGT\n'
    # The whole chromosome, wrapped at 60, past the 2 Mb the JSON route allows.
    whole = _get(served, f'{base}?chrom=Chr1&start=0&end=99999999').get_data(as_text=True)
    assert whole.splitlines()[1:] == ['ACGTACGTACACGTTTTT']
    assert _get(served, f'{base}?chrom=Chr1&start=0&end=5', user='bob').status_code == 404
    assert _get(served, f'{base}?chrom=Nope&start=0&end=5').status_code == 400


def test_a_track_for_another_genome_is_offered_there(served):
    client = served['client']
    other = client.post('/api/genomes', data={
        'name': 'Other', 'assembly_id': 'O1',
        'file': (io.BytesIO(b'>contig_1\n' + b'A' * 50 + b'\n'), 'other.fa')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data').get_json()
    bed = b'contig_1\t1\t20\tx\n'
    res = client.post(f"/api/genomes/{served['genome']}/tracks", data={
        'name': 'Peaks', 'file': (io.BytesIO(bed), 'peaks.bed')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data')
    assert res.status_code == 400
    assert [g['id'] for g in res.get_json()['fits_genomes']] == [other['id']]
    # Bob cannot see either genome, so is offered nothing of alice's.
    lib = client.get(f"/api/genomes/{other['id']}/tracks",
                     headers={'X-Suite-User': 'alice'}).get_json()
    assert lib == []


def test_a_partly_fitting_feature_file_goes_where_it_all_fits(served):
    client = served['client']
    other = client.post('/api/genomes', data={
        'name': 'Bigger', 'assembly_id': 'B1',
        'file': (io.BytesIO(b'>Chr1\n' + b'A' * 40 + b'\n>contig_9\n' + b'C' * 40 + b'\n'), 'b.fa')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data').get_json()
    bed = b'Chr1\t1\t10\tx\ncontig_9\t1\t10\ty\n'
    res = client.post(f"/api/genomes/{served['genome']}/tracks", data={
        'name': 'Mixed', 'file': (io.BytesIO(bed), 'mixed.bed')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data')
    assert res.status_code == 400
    assert [g['id'] for g in res.get_json()['fits_genomes']] == [other['id']]


def test_a_gzipped_fasta_is_unpacked(served):
    res = served['client'].post('/api/genomes', data={
        'name': 'Packed', 'assembly_id': 'P1',
        'file': (io.BytesIO(gzip.compress(b'>c1\nACGT\nAC\n')), 'p.fa.gz')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data')
    assert res.status_code == 201, res.get_json()
    assert res.get_json()['chromosomes'] == [{'name': 'c1', 'length_bp': 6}]


def test_rewrap_keeps_the_sequence(tmp_path):
    src, dst = tmp_path / 'a.fa', tmp_path / 'b.fa'
    src.write_bytes(b'>x desc\nAC\nGTA\n\nCG\n>y\nT\n')
    serving._rewrap_fasta(src, dst, width=3)
    assert dst.read_text() == '>x desc\nACG\nTAC\nG\n>y\nT\n'


def test_split_keeps_separator_characters_inside_lines(tmp_path):
    src = tmp_path / 'f.gff'
    src.write_bytes(b'#h\r\nChr1\t.\tg\t1\t2\t.\t+\t.\tNote=x\x1ey\x0cz\r\ntrack name=x\n')
    head, body = tmp_path / 'h', tmp_path / 'b'
    serving._split_header_and_body(src, head, body, 'gff')
    assert head.read_text() == '#h\n'
    assert body.read_text(encoding='utf-8') == \
        'Chr1\t.\tg\t1\t2\t.\t+\t.\tNote=x y z\n'


def test_an_added_genome_is_prepared_for_checking_primers(served):
    """Its BLAST database is asked for on upload, for the worker to make in the background."""
    import primer_check
    res = served['client'].post('/api/genomes', data={
        'name': 'Fresh', 'assembly_id': 'F1',
        'file': (io.BytesIO(b'>c1\n' + b'ACGT' * 30 + b'\n'), 'f.fa')},
        headers={'X-Suite-User': 'alice'}, content_type='multipart/form-data')
    assert res.status_code == 201, res.get_json()
    assert primer_check.state(served['data_root'], res.get_json()['id'])['state'] == 'building'
