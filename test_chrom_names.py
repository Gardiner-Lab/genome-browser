"""
Tests for the chromosome-name check (``chrom_names``) and its use on the
worker's BAM/wig path.
"""

import uuid

import pytest

import chrom_names
import conversion
import db
import readers
from chrom_names import match_names, name_key


# ---------------------------------------------------------------------------
# Name keys and matching
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('spellings, key', [
    (['1', 'Chr1', 'chr1', 'CHR1', 'chr01', 'Chromosome 1', 'chromosome_1',
      'chrom1', 'Chr_1'], '1'),
    (['M', 'MT', 'Mt', 'ChrM', 'chrMT', 'mito', 'mitochondria',
      'Mitochondrion'], 'M'),
    (['C', 'Pt', 'Cp', 'ChrC', 'chrPt', 'chloroplast', 'Plastid'], 'C'),
    (['scaffold_12', 'Scaffold_12'], 'scaffold_12'),
])
def test_name_key_folds_styles(spellings, key):
    assert {name_key(s) for s in spellings} == {key}


def test_distinct_chromosomes_keep_distinct_keys():
    assert name_key('Chr1') != name_key('Chr10')
    assert name_key('ChrM') != name_key('ChrC')
    assert name_key('2L') != name_key('2R')


def test_tair10_fasta_against_araport11_gff():
    # TAIR10's FASTA as the genome, Araport11's GFF3 as the track.
    genome = ['1', '2', '3', '4', '5', 'mitochondria', 'chloroplast']
    track = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'ChrC', 'ChrM']
    match = match_names(track, genome)
    assert match.unmatched == []
    assert match.renames == {
        'Chr1': '1', 'Chr2': '2', 'Chr3': '3', 'Chr4': '4', 'Chr5': '5',
        'ChrC': 'chloroplast', 'ChrM': 'mitochondria'}


def test_exact_names_are_not_renamed():
    match = match_names(['Chr1', 'Chr2'], ['Chr1', 'Chr2'])
    assert match.renames == {} and match.unmatched == []
    assert match.report() is None


def test_ambiguous_genome_names_are_never_a_rename_target():
    # The genome has both 'chr1' and 'Chr1': '1' could mean either.
    match = match_names(['1', 'chr1'], ['chr1', 'Chr1'])
    assert match.renames == {}
    assert match.unmatched == ['1']


def test_unmatched_and_all_unmatched():
    partial = match_names(['Chr1', 'scaffold_9'], ['1', '2'])
    assert partial.renames == {'Chr1': '1'}
    assert partial.unmatched == ['scaffold_9']
    assert not partial.all_unmatched
    none = match_names(['contig_1'], ['1', '2'])
    assert none.all_unmatched
    report = partial.report()
    assert report['renamed'] == {'Chr1': '1'}
    assert report['not_in_genome'] == ['scaffold_9']
    assert len(report['messages']) == 2


def test_separator_characters_inside_an_attribute_do_not_split_a_line():
    # Araport11's AT5G64510 note carries an ASCII record separator (0x1E);
    # str.splitlines broke the line there and read its tail as chromosome "c".
    gff = ('Chr5\tAraport11\tgene\t10\t90\t.\t+\t.\tID=AT5G64510;'
           'Note=tunicamy\x1ec ER stress-inducible protein\x0cx y\r\n'
           'Chr5\tAraport11\tgene\t95\t99\t.\t+\t.\tID=g2\n')
    assert readers.scan_gff_intervals(gff).chroms == ['Chr5']
    import regions
    feats = regions._scan_gff_features(gff)
    assert [f[0] for f in feats] == ['Chr5', 'Chr5']
    assert readers.split_lines('a\r\nb\x1ec\n') == ['a', 'b\x1ec']


def test_rename_header_uses_genome_names_and_drops_unmatched():
    header = readers.HeaderResult(chrom_extents={'chr1': 90, 'x': 5})
    match = match_names(header.chroms, ['Chr1'])
    renamed = chrom_names.rename_header(header, match)
    assert renamed.chrom_extents == {'Chr1': 90}


# ---------------------------------------------------------------------------
# Rewriting files
# ---------------------------------------------------------------------------


def test_rename_feature_lines():
    renames = {'1': 'Chr1'}
    assert chrom_names.rename_feature_line(
        '1\t.\tgene\t1\t9\t.\t+\t.\tID=1\n', 'gff', renames) == \
        'Chr1\t.\tgene\t1\t9\t.\t+\t.\tID=1\n'
    assert chrom_names.rename_feature_line(
        '##sequence-region 1 1 100\n', 'gff', renames) == \
        '##sequence-region Chr1 1 100\n'
    assert chrom_names.rename_feature_line(
        '# 1 is a comment\n', 'gff', renames) == '# 1 is a comment\n'
    assert chrom_names.rename_feature_line(
        '1 5 9 name\r\n', 'bed', renames) == 'Chr1 5 9 name\r\n'


def test_rename_wig_leaves_step_data_alone():
    wig = ['track type=wiggle_0\n',
           'variableStep chrom=1 span=1\n',
           '1 0.5\n',
           '2 0.7\n',
           'fixedStep chrom=2 start=1 step=1\n',
           '1\n',
           '1\t10\t20\t3.0\n']
    out = list(chrom_names.rename_wig_lines(wig, {'1': 'Chr1', '2': 'Chr2'}))
    assert out == ['track type=wiggle_0\n',
                   'variableStep chrom=Chr1 span=1\n',
                   '1 0.5\n',
                   '2 0.7\n',
                   'fixedStep chrom=Chr2 start=1 step=1\n',
                   '1\n',
                   'Chr1\t10\t20\t3.0\n']


def test_rewrite_bigwig_renames_header_and_keeps_values(tmp_path):
    pyBigWig = pytest.importorskip('pyBigWig')
    src, dst = tmp_path / 'in.bw', tmp_path / 'out.bw'
    bw = pyBigWig.open(str(src), 'w')
    bw.addHeader([('1', 1000), ('Pt', 200)])
    bw.addEntries(['1', '1'], [0, 100], ends=[50, 150], values=[1.5, 2.5])
    bw.addEntries(['Pt'], [10], ends=[20], values=[7.0])
    bw.close()

    chrom_names.rewrite_bigwig(str(src), str(dst), {'1': 'Chr1', 'Pt': 'ChrC'})

    out = pyBigWig.open(str(dst))
    try:
        assert out.chroms() == {'Chr1': 1000, 'ChrC': 200}
        assert out.intervals('Chr1') == ((0, 50, 1.5), (100, 150, 2.5))
        assert out.intervals('ChrC') == ((10, 20, 7.0),)
    finally:
        out.close()


def test_rename_sam_header():
    header = '@HD\tVN:1.6\n@SQ\tSN:1\tLN:100\n@SQ\tSN:Pt\tLN:50\n@PG\tID:x\n'
    assert chrom_names.rename_sam_header(header, {'1': 'Chr1', 'Pt': 'ChrC'}) == \
        '@HD\tVN:1.6\n@SQ\tSN:Chr1\tLN:100\n@SQ\tSN:ChrC\tLN:50\n@PG\tID:x\n'


# ---------------------------------------------------------------------------
# The worker path: a wig in another style is renamed before conversion
# ---------------------------------------------------------------------------


def _seed_converting_track(data_root, chroms):
    db.init_db(data_root / 'genome.db')
    genome_id, track_id = uuid.uuid4().hex, uuid.uuid4().hex
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute(
            'INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, '
            'genome_visibility, created_by) VALUES (?, ?, ?, ?, ?, ?)',
            (genome_id, 'G', 'asm', 'labA', 'owner_lab_only', 'alice'))
        conn.executemany(
            'INSERT INTO chrom_sizes (genome_id, chrom, length_bp) '
            'VALUES (?, ?, ?)', [(genome_id, c, n) for c, n in chroms])
        conn.execute(
            'INSERT INTO tracks (id, genome_id, name, type, owner_lab, '
            'track_status, created_by) VALUES (?, ?, ?, ?, ?, ?, ?)',
            (track_id, genome_id, 'cov', 'bigwig', 'labA', 'converting',
             'alice'))
        conn.commit()
    finally:
        conn.close()
    return genome_id, track_id


def _track_row(data_root, track_id):
    conn = db.connect(data_root / 'genome.db')
    try:
        return conn.execute('SELECT * FROM tracks WHERE id = ?',
                            (track_id,)).fetchone()
    finally:
        conn.close()


def test_worker_renames_wig_before_conversion(tmp_path):
    data_root = tmp_path / 'genome'
    genome_id, track_id = _seed_converting_track(data_root, [('Chr1', 1000)])
    scratch = data_root / 'scratch' / 'job1'
    scratch.mkdir(parents=True)
    wig = scratch / 'input.wig'
    wig.write_text('variableStep chrom=1\n1 0.5\n10 0.7\n', encoding='utf-8')

    seen = {}

    def fake_wig_to_bigwig(input_wig, sizes, out):
        seen['wig'] = open(input_wig, encoding='utf-8').read()
        open(out, 'wb').write(b'bw')

    tools = conversion.Toolchain(wig_to_bigwig=fake_wig_to_bigwig)
    outcome = conversion.convert_job(
        track_id=track_id, genome_id=genome_id, kind='wig', input_path=wig,
        data_root=data_root, samtools_threads=1, toolchain=tools)
    assert outcome.status == 'ready', outcome.error
    assert seen['wig'] == 'variableStep chrom=Chr1\n1 0.5\n10 0.7\n'
    import json
    check = json.loads(_track_row(data_root, track_id)['chrom_check'])
    assert check['renamed'] == {'1': 'Chr1'}


def test_worker_reheaders_bam_in_another_style(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome'
    genome_id, track_id = _seed_converting_track(data_root, [('Chr1', 1000)])
    scratch = data_root / 'scratch' / 'job2'
    scratch.mkdir(parents=True)
    bam = scratch / 'input.bam'
    bam.write_bytes(b'bam')
    monkeypatch.setattr(readers, '_samtools_header_reader',
                        lambda path: '@SQ\tSN:chr1\tLN:1000\n')
    monkeypatch.setattr(
        readers, 'read_bam_header',
        lambda path, reader=None: readers.parse_bam_header_text(
            '@SQ\tSN:chr1\tLN:1000\n'))

    calls = []

    def fake_reheader(src, dst, renames):
        calls.append(renames)
        open(dst, 'wb').write(b'reheadered')

    def fake_coverage(sorted_bam, out, threads):
        open(out, 'wb').write(b'bw')

    tools = conversion.Toolchain(
        samtools_sort=lambda i, o, t, m: open(o, 'wb').write(
            open(i, 'rb').read()),
        samtools_index=lambda b: None,
        bamcoverage=fake_coverage,
        samtools_reheader=fake_reheader)
    monkeypatch.setattr(conversion.shutil, 'which', lambda name: '/bin/x')
    outcome = conversion.convert_job(
        track_id=track_id, genome_id=genome_id, kind='bam', input_path=bam,
        data_root=data_root, samtools_threads=1, toolchain=tools)
    assert outcome.status == 'ready', outcome.error
    assert calls == [{'chr1': 'Chr1'}]


def test_worker_refuses_unmatched_names_listing_the_genome(tmp_path):
    data_root = tmp_path / 'genome'
    genome_id, track_id = _seed_converting_track(
        data_root, [('Chr1', 1000), ('Chr2', 500)])
    scratch = data_root / 'scratch' / 'job3'
    scratch.mkdir(parents=True)
    wig = scratch / 'input.wig'
    wig.write_text('variableStep chrom=contig_7\n1 0.5\n', encoding='utf-8')
    outcome = conversion.convert_job(
        track_id=track_id, genome_id=genome_id, kind='wig', input_path=wig,
        data_root=data_root, samtools_threads=1,
        toolchain=conversion.Toolchain())
    assert outcome.status == 'failed'
    assert 'contig_7' in outcome.error
    assert 'Chr1, Chr2' in outcome.error


def test_existing_database_gains_the_column(tmp_path):
    path = tmp_path / 'old.db'
    conn = db.connect(path)
    conn.execute('CREATE TABLE tracks (id TEXT PRIMARY KEY, genome_id TEXT, '
                 'name TEXT, type TEXT, owner_lab TEXT, visibility TEXT, '
                 'track_status TEXT, error TEXT, file_path TEXT, '
                 'created_by TEXT, created_at TEXT)')
    conn.commit()
    conn.close()
    db.init_db(path)
    conn = db.connect(path)
    try:
        cols = {row[1] for row in conn.execute('PRAGMA table_info(tracks)')}
    finally:
        conn.close()
    assert 'chrom_check' in cols
