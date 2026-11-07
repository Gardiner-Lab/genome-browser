"""The genome the browser ships with is added once, as the suite's, and never twice."""

import gzip

import pytest

import access
import bundled
import db

FASTA = '>Chr1 test\n' + 'ACGT' * 500 + '\n>ChrM\n' + 'GATTACA' * 50 + '\n'
GFF = '##gff-version 3\nChr1\tAraport11\tgene\t10\t200\t.\t+\t.\tID=AT1G01010;Name=NAC001\n'


@pytest.fixture
def setup(tmp_path, monkeypatch):
    source = tmp_path / 'bundled'
    source.mkdir()
    (source / bundled.FASTA).write_bytes(gzip.compress(FASTA.encode()))
    (source / bundled.GFF).write_bytes(gzip.compress(GFF.encode()))
    root = tmp_path / 'data'
    for sub in ('genomes', 'tracks', 'scratch'):
        (root / sub).mkdir(parents=True)
    db.init_db(root / 'genome.db')
    monkeypatch.delenv('GENOME_BUNDLED_ARAPORT11', raising=False)
    conn = db.connect(root / 'genome.db')
    yield conn, root, source
    conn.close()


def test_it_is_added_once_as_the_suites(setup):
    conn, root, source = setup
    genome_id = bundled.seed(conn, root, source)
    assert genome_id
    genome = conn.execute('SELECT * FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()
    assert genome['owner_lab'] == access.SUITE_OWNER and genome['genome_visibility'] == access.SUITE_VISIBILITY
    assert dict(conn.execute('SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?', (genome_id,)).fetchall()) \
        == {'Chr1': 2000, 'ChrM': 350}
    track = conn.execute('SELECT * FROM tracks WHERE genome_id = ?', (genome_id,)).fetchone()
    assert track['visibility'] == access.SUITE_VISIBILITY and track['track_status'] == 'ready'
    assert (root / 'genomes' / genome_id / 'reference.fasta').exists()
    # Not again, even after an administrator deletes it.
    conn.execute('DELETE FROM tracks')
    conn.execute('DELETE FROM chrom_sizes')
    conn.execute('DELETE FROM reference_genomes')
    conn.commit()
    assert bundled.seed(conn, root, source) is None


def test_a_suite_tair10_uploaded_by_hand_is_kept(setup):
    conn, root, source = setup
    conn.execute("INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) "
                 "VALUES ('mine', 'Mine', 'TAIR10', 'suite', 'suite', 'admin')")
    conn.commit()
    assert bundled.seed(conn, root, source) is None
    assert conn.execute('SELECT count(*) FROM reference_genomes').fetchone()[0] == 1


def test_it_can_be_turned_off(setup, monkeypatch):
    conn, root, source = setup
    monkeypatch.setenv('GENOME_BUNDLED_ARAPORT11', '0')
    assert bundled.seed(conn, root, source) is None
