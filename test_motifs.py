"""
Marking a sequence on the genome: finding a guide or a primer on both strands,
and keeping the hits as a feature track of their own.
"""

import uuid

import pytest

import motifs


def write_fasta(path, chroms):
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        for name, seq in chroms.items():
            handle.write(f'>{name} some description\n')
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + '\n')
    return path


def test_a_sequence_is_found_on_both_strands(tmp_path):
    guide = 'GACGTTAGCACCTGTAGGAT'
    rc = motifs.reverse_complement(guide)
    # Once as written, once as its reverse complement, and padding that is neither.
    seq = 'A' * 100 + guide + 'C' * 50 + rc + 'T' * 80
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    hits, capped = motifs.search(fasta, [motifs.Query('sgRNA1', guide)], {'Chr1': len(seq)})
    assert not capped
    assert [(h.start, h.strand) for h in hits] == [(100, '+'), (170, '-')]
    assert {h.chrom for h in hits} == {'Chr1'} and {h.name for h in hits} == {'sgRNA1'}
    # Its coordinates are the bases themselves.
    assert seq[hits[0].start:hits[0].end] == guide
    assert motifs.reverse_complement(seq[hits[1].start:hits[1].end]) == guide


def test_a_match_across_the_readers_chunk_boundary_is_still_found(tmp_path, monkeypatch):
    monkeypatch.setattr(motifs, 'CHUNK', 1000)
    primer = 'TTCGCAAGTGCCATAGGTCA'
    # Straddling the boundary: 10 bases before it, 10 after.
    seq = 'A' * 990 + primer + 'G' * 500
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    hits, _ = motifs.search(fasta, [motifs.Query('P1', primer)], {'Chr1': len(seq)})
    assert [(h.start, h.end) for h in hits] == [(990, 1010)], 'the hit over the boundary is missed or doubled'


def test_iupac_codes_match_any_of_their_bases(tmp_path):
    seq = 'TTTT' + 'GGTCTCA' + 'AAAA' + 'GGTCTCG' + 'CCCC'
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    hits, _ = motifs.search(fasta, [motifs.Query('BsaI', 'GGTCTCNNNNN')], {'Chr1': len(seq)})
    assert len(hits) == 2 and all(h.name == 'BsaI' for h in hits)


def test_overlapping_matches_are_all_found(tmp_path):
    # A sequence that repeats head to tail is at every position, not every other one.
    seq = 'CA' * 30
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    hits, _ = motifs.search(fasta, [motifs.Query('r', 'CACACACA')], {'Chr1': len(seq)})
    assert [h.start for h in hits][:3] == [0, 2, 4]


def test_a_palindrome_is_reported_once(tmp_path):
    site = 'GAATTCGAATTC'                      # its own reverse complement
    assert motifs.reverse_complement(site) == site
    seq = 'TTTGCATTAG' + site + 'ACCTGAGTTA'
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    hits, _ = motifs.search(fasta, [motifs.Query('EcoRI pair', site)], {'Chr1': len(seq)})
    assert len(hits) == 1 and hits[0].strand == '+', 'a palindrome is not two features'


def test_several_named_sequences_at_once_and_their_numbering(tmp_path):
    one, two = 'GACGTTAGCACCTGTAGGAT', 'TTCGCAAGTGCCATAGGTCA'
    seq = one + 'A' * 50 + two + 'C' * 50 + one
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': seq})
    queries = motifs.parse_queries(f'sgRNA1\t{one}\nsgRNA2 {two}', 'Guides')
    assert [q.name for q in queries] == ['sgRNA1', 'sgRNA2']
    hits = motifs.number_hits(motifs.search(fasta, queries, {'Chr1': len(seq)})[0])
    # The one found twice is numbered; the one found once keeps its plain name.
    assert sorted(h.name for h in hits) == ['sgRNA1.1', 'sgRNA1.2', 'sgRNA2']
    bed = motifs.bed_text(hits).splitlines()
    assert len(bed) == 3 and bed[0].split('\t')[0] == 'Chr1' and bed[0].split('\t')[5] in '+-'


def test_what_a_member_types_is_read_forgivingly():
    one = 'gacgttagcacctgtaggat'
    queries = motifs.parse_queries(f'# a comment\n\n>sgRNA1\n{one}\n  GGTCTCAAAAAA  \n', 'Guides')
    assert queries[0] == motifs.Query('sgRNA1', one.upper())
    assert queries[1].name == 'Guides', 'an unnamed sequence takes the track name'
    # A FASTA header, spaces and dashes in the sequence, and lower case are all fine.
    assert motifs.parse_queries('P1 ggtc-tca.aaaa', 'T')[0].sequence == 'GGTCTCAAAAA'


@pytest.mark.parametrize('text, says', [
    ('ACGTACGTX ', 'not a base'),
    ('ACGTA', 'at least'),
    ('', 'at least one sequence'),
    ('A' * 1001, 'at most'),
])
def test_a_sequence_that_cannot_be_searched_says_why(text, says):
    with pytest.raises(ValueError) as raised:
        motifs.parse_queries(text, 'T')
    assert says in str(raised.value)


def test_one_sequence_everywhere_does_not_crowd_out_a_guide(tmp_path, monkeypatch):
    """An enzyme site is in a genome thousands of times; a guide typed beside it is found all the same."""
    monkeypatch.setattr(motifs, 'MAX_PER_SEQUENCE', 10)
    guide = 'GACGTTAGCACCTGTAGGAT'
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': 'ACGTACGTAC' * 500 + guide})
    hits, capped = motifs.search(fasta, [motifs.Query('site', 'ACGTACGT'), motifs.Query('sgRNA1', guide)],
                                 {'Chr1': 5020})
    assert capped == ['site'], 'the common site should be the one that stops'
    assert sum(1 for h in hits if h.name == 'site') == 10
    assert [h.name for h in hits if h.name == 'sgRNA1'] == ['sgRNA1'], 'the guide was crowded out'


def test_a_chromosome_the_genome_does_not_list_is_skipped(tmp_path):
    guide = 'GACGTTAGCACCTGTAGGAT'
    fasta = write_fasta(tmp_path / 'g.fa', {'Chr1': 'A' * 40 + guide, 'scaffold_9': 'T' * 40 + guide})
    hits, _ = motifs.search(fasta, [motifs.Query('g', guide)], {'Chr1': 60})
    assert [h.chrom for h in hits] == ['Chr1']


# ---------------------------------------------------------------------------
# The endpoint: a search becomes a track like any other
# ---------------------------------------------------------------------------

import sqlite3                                                  # noqa: E402

import app as app_module                                        # noqa: E402
import db                                                       # noqa: E402
import identity                                                 # noqa: E402
import labs_reader                                              # noqa: E402
import serving                                                  # noqa: E402

GUIDE = 'GACGTTAGCACCTGTAGGAT'


@pytest.fixture
def genome(tmp_path, monkeypatch):
    """A tiny one-chromosome genome of Lab A's, with the guide in it twice."""
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    genome_id = uuid.uuid4().hex
    seq = 'ACGTTGCAAT' * 20 + GUIDE + 'TTGCACGATT' * 20 + motifs.reverse_complement(GUIDE) + 'GGCATTACGA' * 20
    fasta = serving.reference_path(data_root, genome_id)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(fasta, {'Chr1': seq})
    conn = db.connect(data_root / 'genome.db')
    try:
        conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                     "VALUES (?, 'G', 'asm', 'labA', 'lab_only', 'alice')", (genome_id,))
        conn.execute('INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
                     (genome_id, 'Chr1', len(seq)))
        conn.commit()
    finally:
        conn.close()
    monkeypatch.setattr(labs_reader, 'labs_of',
                        lambda u: [{'lab_id': 'labA', 'lab': 'Lab A', 'role': 'member'}] if u == 'alice' else [])
    monkeypatch.setattr(labs_reader, 'current_lab', lambda u: 'labA' if u == 'alice' else None)
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    return {'client': application.test_client(), 'id': genome_id, 'data_root': data_root}


def test_marking_a_sequence_makes_a_feature_track(genome):
    c = genome['client']
    r = c.post(f"/api/genomes/{genome['id']}/tracks/from-sequence", headers={'X-Suite-User': 'alice'},
               json={'name': 'CRISPR guides', 'sequences': f'sgRNA1 {GUIDE}', 'experiment': 'Editing'})
    assert r.status_code == 201, r.get_json()
    body = r.get_json()
    assert body['found'] == 2 and body['per_sequence'] == {'sgRNA1': 2} and body['type'] == 'bed'

    # It is an ordinary track of the lab's from here on, and it is on the genome's list.
    listed = {t['name']: t for t in c.get(f"/api/genomes/{genome['id']}/tracks",
                                          headers={'X-Suite-User': 'alice'}).get_json()}
    assert listed['CRISPR guides']['experiment'] == 'Editing'
    assert listed['CRISPR guides']['type'] == 'bed' and listed['CRISPR guides']['can_change'] is True

    # Its file is BED: one line per place, on the strand it reads on.
    conn = db.connect(genome['data_root'] / 'genome.db')
    try:
        path = conn.execute('SELECT file_path FROM tracks WHERE id = ?', (body['id'],)).fetchone()[0]
    finally:
        conn.close()
    lines = [l.split('\t') for l in open(path, encoding='utf-8').read().splitlines()]
    assert [l[0] for l in lines] == ['Chr1', 'Chr1']
    assert sorted(l[3] for l in lines) == ['sgRNA1.1', 'sgRNA1.2']
    assert sorted(l[5] for l in lines) == ['+', '-'], 'both strands are marked'


def test_a_sequence_that_is_not_there_makes_no_track(genome):
    c = genome['client']
    r = c.post(f"/api/genomes/{genome['id']}/tracks/from-sequence", headers={'X-Suite-User': 'alice'},
               json={'name': 'nothing', 'sequences': 'CCCCCCCCGGGGGGGGCCCCCCCC'})
    assert r.status_code == 400 and 'Not found on this genome' in r.get_json()['error']
    assert c.get(f"/api/genomes/{genome['id']}/tracks", headers={'X-Suite-User': 'alice'}).get_json() == []


def test_marking_a_sequence_needs_a_name_a_sequence_and_a_lab(genome):
    c = genome['client']
    url = f"/api/genomes/{genome['id']}/tracks/from-sequence"
    assert c.post(url, headers={'X-Suite-User': 'alice'}, json={'sequences': GUIDE}).status_code == 400
    assert c.post(url, headers={'X-Suite-User': 'alice'}, json={'name': 'x', 'sequences': ''}).status_code == 400
    bad = c.post(url, headers={'X-Suite-User': 'alice'}, json={'name': 'x', 'sequences': 'ACGTZZZZACGT'})
    assert bad.status_code == 400 and 'not a base' in bad.get_json()['error']
    # Someone who cannot see the genome is told it does not exist.
    assert c.post(url, headers={'X-Suite-User': 'bob'}, json={'name': 'x', 'sequences': GUIDE}).status_code == 404
