"""Primers checked against a genome with BLAST (primer_check.py)."""
import random
import shutil
import uuid

import pytest

import app as app_module
import db
import identity
import labs_reader
import primer_check
import serving


def _hit(chrom='Chr1', pident=100.0, length=25, qstart=1, qend=25, sstart=101, send=125, bitscore=50.0, evalue=1e-5):
    return {'chrom': chrom, 'pident': pident, 'length': length, 'qstart': qstart, 'qend': qend,
            'sstart': sstart, 'send': send, 'bitscore': bitscore, 'evalue': evalue}


PRIMER = {'id': '7', 'sequence': 'A' * 25}


def test_a_primer_with_one_perfect_site_passes():
    v = primer_check.verdict(PRIMER, [_hit(), _hit(pident=90.0, length=15, bitscore=20.0, sstart=900, send=886)], 0.75, 80)
    assert v['ok'] and v['site'] == {'chrom': 'Chr1', 'start': 101, 'end': 125, 'strand': '+', 'identity': 100.0, 'length': 25}
    # The next best is shown, on the minus strand, though too short to matter (15 of 25 nt).
    assert v['other']['strand'] == '−' and v['other']['start'] == 886 and v['off_target'] == 0


def test_a_second_good_match_turns_a_primer_down():
    near = _hit(chrom='Chr3', pident=92.0, length=24, bitscore=40.0)
    v = primer_check.verdict(PRIMER, [_hit(), near], 0.75, 80)
    assert not v['ok'] and v['off_target'] == 1 and v['other']['chrom'] == 'Chr3'
    # The cutoffs are the caller's: a match over 24 of 25 nt passes if it must cover all of it.
    assert primer_check.verdict(PRIMER, [_hit(), near], 1.0, 80)['ok']
    # Two perfect sites: the second is off-target.
    assert not primer_check.verdict(PRIMER, [_hit(), _hit(chrom='Chr2')], 0.75, 80)['ok']
    # Any hit counts, not only the second best (the script looked at the second only).
    third = _hit(chrom='Chr4', pident=85.0, length=20, bitscore=30.0)
    weak = _hit(chrom='Chr5', pident=100.0, length=12, bitscore=35.0)
    assert not primer_check.verdict(PRIMER, [_hit(), weak, third], 0.75, 80)['ok']


def test_a_primer_not_in_the_genome_passes_without_a_site():
    v = primer_check.verdict(PRIMER, [], 0.75, 80)
    assert v['ok'] and v['site'] is None and v['other'] is None and v['hits'] == 0
    # Not found whole, but a partial match: no site, and judged by that match.
    v = primer_check.verdict(PRIMER, [_hit(pident=95.0, length=22, qend=22)], 0.75, 80)
    assert not v['ok'] and v['site'] is None


def test_blast_output_is_read_by_query():
    out = 'q0\tChr1 CHROMOSOME dumped\t100.000\t25\t1\t25\t101\t125\t1e-05\t50.1\nq1\tChr2\t88.0\t20\t3\t22\t50\t31\t0.1\t30\nbad line\n'
    parsed = primer_check.parse(out)
    assert parsed['q0'][0]['chrom'] == 'Chr1' and parsed['q1'][0]['send'] == 31 and set(parsed) == {'q0', 'q1'}


@pytest.mark.parametrize('primers, message', [
    ([], 'list'), ('ACGT', 'list'), ([{'id': 1, 'sequence': 'ACGT'}], 'nt'),
    ([{'id': 1, 'sequence': 'ACGTYACGTACGTACGTAC'}], 'nt'),
    ([{'id': 1, 'sequence': 'ACGT' * 5}, {'id': 1, 'sequence': 'ACGT' * 5}], 'twice'),
    ([{'id': i, 'sequence': 'ACGT' * 5} for i in range(primer_check.MAX_PRIMERS + 1)], 'at most'),
])
def test_primers_are_checked(primers, message):
    with pytest.raises(ValueError, match=message):
        primer_check.clean_primers(primers)


# --- the endpoint -------------------------------------------------------------------------------

def _genome_sequence():
    rng = random.Random(3)
    unique = ''.join(rng.choices('ACGT', k=3000))
    repeat = ''.join(rng.choices('ACGT', k=40))
    # The repeat is in the genome twice, once on each strand.
    rc = repeat[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    return unique[:1000] + repeat + unique[1000:2000] + rc + unique[2000:], unique, repeat


@pytest.fixture
def genome(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    genome_id = uuid.uuid4().hex
    seq, unique, repeat = _genome_sequence()
    fasta = serving.reference_path(data_root, genome_id)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta.write_text('>Chr1 test chromosome\n' + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n')
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                     "VALUES (?, 'Testcress', 'asm', 'labA', 'lab_only', 'alice')", (genome_id,))
        conn.execute('INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)', (genome_id, 'Chr1', len(seq)))
        conn.commit()
    finally:
        conn.close()
    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A', 'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    return {'client': application.test_client(), 'id': genome_id, 'data_root': data_root,
            'unique': unique, 'repeat': repeat}


def _post(genome, primers, user='alice', **extra):
    return genome['client'].post(f"/api/genomes/{genome['id']}/primer-check", headers={'X-Suite-User': user},
                                 json={'primers': primers, **extra})


def test_only_those_who_may_see_the_genome_check_against_it(genome):
    assert _post(genome, [{'id': 'a', 'sequence': 'ACGT' * 6}], user='mallory').status_code == 404
    assert _post(genome, [{'id': 'a', 'sequence': 'ACGU' * 6}]).status_code == 400
    assert _post(genome, [{'id': 'a', 'sequence': 'ACGT' * 6}], max_match=2).status_code == 400


def test_the_first_check_asks_the_worker_for_a_database(genome, monkeypatch):
    r = _post(genome, [{'id': 'a', 'sequence': 'ACGT' * 6}])
    assert r.status_code == 202 and 'Testcress' in r.get_json()['message']
    assert primer_check.state(genome['data_root'], genome['id'])['state'] == 'building'
    # A database that could not be made says why, rather than waiting for ever.
    monkeypatch.setattr(primer_check, 'MAKEBLASTDB', 'no-such-makeblastdb')
    assert primer_check.build_next(genome['data_root'])
    assert not primer_check.build_next(genome['data_root'])          # nothing more asked for
    r = _post(genome, [{'id': 'a', 'sequence': 'ACGT' * 6}])
    assert r.status_code == 503 and 'could not be prepared' in r.get_json()['error']


@pytest.mark.skipif(not shutil.which('blastn') or not shutil.which('makeblastdb'),
                    reason='BLAST+ is in the image, not on this machine (the image CI job runs this)')
def test_primers_are_checked_with_blast(genome):
    assert _post(genome, [{'id': 'a', 'sequence': 'ACGT' * 6}]).status_code == 202
    assert primer_check.build_next(genome['data_root'])
    unique, repeat = genome['unique'], genome['repeat']
    r = _post(genome, [{'id': 'once', 'sequence': unique[500:526]},
                       {'id': 'twice', 'sequence': repeat[5:31]},
                       {'id': 'nowhere', 'sequence': 'ACGTTGCA' * 3}])
    assert r.status_code == 200, r.get_json()
    v = {p['id']: p for p in r.get_json()['primers']}
    assert v['once']['ok'] and v['once']['site'] == {'chrom': 'Chr1', 'start': 501, 'end': 526, 'strand': '+',
                                                      'identity': 100.0, 'length': 26}
    # Twice, once on each strand: whichever BLAST puts first is its site, the other off-target.
    assert not v['twice']['ok'] and v['twice']['off_target'] == 1
    assert {v['twice']['site']['strand'], v['twice']['other']['strand']} == {'+', '−'}
    assert v['nowhere']['ok'] and v['nowhere']['site'] is None


def test_the_worker_asks_at_start_for_genomes_with_no_database(genome, monkeypatch):
    conn = db.connect(genome['data_root'] / 'genome.db')
    try:
        assert primer_check.request_missing(conn, genome['data_root']) == 1
        assert primer_check.state(genome['data_root'], genome['id'])['state'] == 'building'
        # Made, it is not asked for again; a sequence changed since, it is.
        ready = primer_check.db_dir(genome['data_root'], genome['id'])
        ready.mkdir(parents=True)
        (primer_check._wanted(genome['data_root']) / genome['id']).unlink()
        (ready / 'ready').touch()
        assert primer_check.request_missing(conn, genome['data_root']) == 0
        import os, time
        later = time.time() + 60
        os.utime(primer_check._reference(genome['data_root'], genome['id']), (later, later))
        assert primer_check.request_missing(conn, genome['data_root']) == 1
    finally:
        conn.close()
