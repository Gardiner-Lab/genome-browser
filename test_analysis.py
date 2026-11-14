"""The Analyse panel: features as regions, the job queue, who sees what, and (with deepTools) each analysis."""

import json
import shutil
import uuid

import pytest

import analysis
import app as app_module
import db
import identity
import labs_reader

GFF = ('##gff-version 3\n'
       'Chr1\tx\tgene\t1001\t3000\t.\t+\t.\tID=g1;Name=GENE1\n'
       'Chr1\tx\tgene\t6001\t9000\t.\t-\t.\tID=g2;Name=GENE2\n'
       'Chr1\tx\ttransposable_element\t12001\t14000\t.\t+\t.\tID=t1;Name=TE1\n'
       'Chr2\tx\tgene\t2001\t4000\t.\t+\t.\tID=g3;Name=GENE3\n')
HAS_DEEPTOOLS = shutil.which('computeMatrix') is not None


def _bigwig(path, high):
    """A small coverage track: high over the genes, low elsewhere."""
    import pyBigWig
    bw = pyBigWig.open(str(path), 'w')
    bw.addHeader([('Chr1', 20000), ('Chr2', 10000)])
    for chrom, length in (('Chr1', 20000), ('Chr2', 10000)):
        starts = list(range(0, length, 100))
        values = [high if (1000 <= s < 3000 or 6000 <= s < 9000 or (chrom == 'Chr2' and 2000 <= s < 4000)) else 1.0
                  for s in starts]
        bw.addEntries([chrom] * len(starts), starts, ends=[s + 100 for s in starts], values=values)
    bw.close()


@pytest.fixture
def world(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    (data_root / 'tracks').mkdir(parents=True, exist_ok=True)
    gff = data_root / 'tracks' / 'genes.gff'
    gff.write_text(GFF)
    genome = uuid.uuid4().hex
    ids = {k: uuid.uuid4().hex for k in ('genes', 'wt', 'mut', 'theirs')}
    conn = db.connect(data_root / 'genome.db')
    conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                 "VALUES (?, 'Plant', 'P1', 'labA', 'owner_lab_only', 'alice')", (genome,))
    conn.executemany('INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
                     [(genome, 'Chr1', 20000), (genome, 'Chr2', 10000)])
    rows = [(ids['genes'], 'Genes', 'gff', 'labA', str(gff))]
    try:
        for key, high in (('wt', 10.0), ('mut', 4.0), ('theirs', 7.0)):
            path = data_root / 'tracks' / f'{key}.bw'
            _bigwig(path, high)
            rows.append((ids[key], key.upper(), 'bigwig', 'labB' if key == 'theirs' else 'labA', str(path)))
    except ImportError:
        pass
    conn.executemany('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, file_path, '
                     "created_by) VALUES (?, ?, ?, ?, ?, 'lab_only', 'ready', ?, 'alice')",
                     [(i, genome, n, t, lab, p) for i, n, t, lab, p in rows])
    conn.commit()
    conn.close()
    members = {'alice': 'labA', 'bob': 'labB'}
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: [
        {'lab_id': members[u], 'lab': members[u], 'role': 'member'}] if u in members else [])
    monkeypatch.setattr(labs_reader, 'current_lab', lambda u: members.get(u))
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    monkeypatch.setattr(labs_reader, 'demo_lab_ids', lambda: set())
    client = application.test_client()

    def call(method, path, user='alice', body=None):
        return client.open(path, method=method, json=body, headers={'X-Suite-User': user})
    return {'call': call, 'genome': genome, 'ids': ids, 'root': data_root, 'gff': gff}


def _run_queue(world):
    conn = db.connect(world['root'] / 'genome.db')
    try:
        while True:
            job = analysis.claim_next(conn)
            if job is None:
                return
            analysis.run_job(conn, world['root'], job)
    finally:
        conn.close()


def test_features_as_regions(world, tmp_path):
    out = tmp_path / 'r.bed'
    assert analysis.region_bed(str(world['gff']), 'gff', 'gene', None, out) == 3
    assert out.read_text().splitlines()[1] == 'Chr1\t6000\t9000\tGENE2\t0\t-'
    assert analysis.region_bed(str(world['gff']), 'gff', None, ['te1', 'gene3'], out) == 2
    types = dict(analysis.feature_types(str(world['gff']), 'gff'))
    assert types == {'gene': 3, 'transposable_element': 1}


def test_who_may_start_and_see_an_analysis(world):
    w, ids = world, world['ids']
    path = f"/api/genomes/{w['genome']}/analyses"
    regions = [{'track': ids['genes'], 'type': 'gene'}]
    assert w['call']('POST', path, 'bob', {'kind': 'genes', 'tracks': [ids['theirs']], 'regions': regions}).status_code == 404
    bad = w['call']('POST', path, 'alice', {'kind': 'genes', 'tracks': [ids['theirs']], 'regions': regions})
    assert bad.status_code == 400 and 'not yours' in bad.get_json()['error']
    assert w['call']('POST', path, 'alice', {'kind': 'nonsense'}).status_code == 400
    assert w['call']('POST', path, 'alice', {'kind': 'genes', 'tracks': [ids['genes']], 'regions': regions}).status_code == 400
    assert w['call']('POST', path, 'alice', {'kind': 'compare', 'a': ids['wt'], 'b': ids['wt']}).status_code == 400
    ok = w['call']('POST', path, 'alice', {'kind': 'genes', 'tracks': [ids['wt']], 'regions': regions})
    assert ok.status_code == 202
    job = ok.get_json()['id']
    assert w['call']('GET', f'/api/analyses/{job}').get_json()['state'] == 'queued'
    assert w['call']('GET', f'/api/analyses/{job}', 'bob').status_code == 404      # only the one who asked
    types = w['call']('GET', f"/api/tracks/{ids['genes']}/feature-types").get_json()['types']
    assert types[0] == {'type': 'gene', 'count': 3}


def test_a_failed_analysis_says_why(world):
    w, ids = world, world['ids']
    w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice',
              {'kind': 'genes', 'tracks': [ids['wt']], 'regions': [{'track': ids['genes'], 'type': 'gene', 'names': ['NOPE']}]})
    _run_queue(w)
    conn = db.connect(w['root'] / 'genome.db')
    row = conn.execute('SELECT state, error FROM analysis_jobs').fetchone()
    conn.close()
    assert row['state'] == 'failed' and 'were found' in row['error']


def test_profile_counts_missing_coverage_as_zero(monkeypatch, tmp_path):
    # Without it, bins a coverage track has no entries for are NaN: black in the heatmap.
    runs = []
    monkeypatch.setattr(analysis, 'region_bed', lambda *a, **k: 5)
    monkeypatch.setattr(analysis, '_run', lambda args, cwd, what: runs.append(args))
    monkeypatch.setattr(analysis, '_features_drawn', lambda matrix: [5])
    monkeypatch.setattr(analysis, '_off_chromosome_ends', lambda *a: [])
    analysis.run_profile({'regions': [{'track': 'g', 'type': 'gene'}], 'heatmap': True},
                         [{'file_path': 'wt.bw', 'name': 'WT'}], {'g': {'file_path': 'g.gff', 'type': 'gff', 'name': 'Genes'}},
                         tmp_path)
    compute = next(a for a in runs if a[0] == 'computeMatrix')
    assert '--missingDataAsZero' in compute
    assert [a[0] for a in runs] == ['computeMatrix', 'plotProfile', 'plotHeatmap']


def test_profile_choices_reach_compute_matrix(monkeypatch, tmp_path):
    def run_with(**choices):
        runs = []
        monkeypatch.setattr(analysis, 'region_bed', lambda *a, **k: 5)
        monkeypatch.setattr(analysis, '_run', lambda args, cwd, what: runs.append(args))
        monkeypatch.setattr(analysis, '_features_drawn', lambda matrix: [5])
        monkeypatch.setattr(analysis, '_off_chromosome_ends', lambda *a: [])
        out = analysis.run_profile({'regions': [{'track': 'g', 'type': 'gene'}], **choices},
                                   [{'file_path': 'wt.bw', 'name': 'WT'}],
                                   {'g': {'file_path': 'g.gff', 'type': 'gff', 'name': 'Genes'}}, tmp_path)
        return next(a for a in runs if a[0] == 'computeMatrix'), out
    args, out = run_with()                                   # as before they could choose: both on
    assert '--skipZeros' in args and '--missingDataAsZero' in args
    assert out['skip_zeros'] is True and out['missing_as_zero'] is True
    args, out = run_with(skip_zeros=False, missing_as_zero=False)
    assert '--skipZeros' not in args and '--missingDataAsZero' not in args
    assert out['skip_zeros'] is False and out['missing_as_zero'] is False
    args, _ = run_with(skip_zeros=False, missing_as_zero=True)
    assert '--skipZeros' not in args and '--missingDataAsZero' in args


def test_features_whose_window_runs_off_a_chromosome_end(tmp_path):
    pyBigWig = pytest.importorskip('pyBigWig')
    import gzip
    bw = pyBigWig.open(str(tmp_path / 't.bw'), 'w')
    bw.addHeader([('Chr1', 10_000)])
    bw.addEntries(['Chr1'], [0], ends=[10_000], values=[1.0])
    bw.close()
    with gzip.open(tmp_path / 'matrix.gz', 'wt') as fh:
        fh.write('@{"group_boundaries": [0, 4]}\n')
        for row in (('Chr1', 500, 1500, 'NEAR_START', '+'), ('Chr1', 4000, 6000, 'MIDDLE', '+'),
                    ('Chr1', 9000, 9800, 'MINUS_NEAR_END', '-'), ('Chr1', 300, 2500, 'MINUS_NEAR_START', '-')):
            fh.write('\t'.join(map(str, row[:4])) + '\t0\t' + row[4] + '\t1\t1\n')
    tracks = [{'file_path': str(tmp_path / 't.bw')}]

    def off(mode, before=2000, after=1000):
        return {r[0]: (r[5], r[6]) for r in analysis._off_chromosome_ends(tmp_path / 'matrix.gz', tracks, mode, before, after)}

    # Chr1 is 10 kb; 2 kb before and 1 kb after each feature, read along its own strand.
    # Scaled: + spans [start - 2000, end + 1000]; - spans [start - 1000, end + 2000].
    assert off('scale') == {'NEAR_START': (1500, 0), 'MINUS_NEAR_END': (0, 1800), 'MINUS_NEAR_START': (700, 0)}
    # Around the start (TSS): a - feature starts at its end coordinate, [end - 1000, end + 2000].
    assert off('start') == {'NEAR_START': (1500, 0), 'MINUS_NEAR_END': (0, 1800)}
    # Around the end (TES): + at its end, [end - 2000, end + 1000]; - at its start, [start - 1000, start + 2000].
    assert off('end') == {'NEAR_START': (500, 0), 'MINUS_NEAR_END': (0, 1000), 'MINUS_NEAR_START': (700, 0)}
    # Short flanks: nothing runs off.
    assert off('scale', before=200, after=200) == {}


def test_profile_label_counts_the_features_drawn(monkeypatch, tmp_path):
    # computeMatrix leaves out features without signal and shorter than a bin: say so, not the count found.
    import gzip
    runs = []
    monkeypatch.setattr(analysis, 'region_bed', lambda *a, **k: 40)

    def fake_run(args, cwd, what):
        runs.append(args)
        if args[0] == 'computeMatrix':
            with gzip.open(tmp_path / 'matrix.gz', 'wt') as fh:
                fh.write('@' + json.dumps({'group_boundaries': [0, 31, 71]}) + '\n')
    monkeypatch.setattr(analysis, '_run', fake_run)
    monkeypatch.setattr(analysis, '_off_chromosome_ends', lambda *a: [])
    out = analysis.run_profile({'regions': [{'track': 'g', 'type': 'gene'}, {'track': 'g', 'type': 'transposable_element'}]},
                               [{'file_path': 'wt.bw', 'name': 'WT'}],
                               {'g': {'file_path': 'g.gff', 'type': 'gff', 'name': 'Genes'}}, tmp_path)
    assert out['regions'] == ['gene (31 of 40)', 'transposable_element (40)']
    assert out['left_out'] == {'gene': 9}
    plotted = next(a for a in runs if a[0] == 'plotProfile')
    assert plotted[plotted.index('--regionsLabel') + 1:][:2] == ['gene (31 of 40)', 'transposable_element (40)']


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_every_analysis_with_deeptools(world):
    w, ids = world, world['ids']
    path = f"/api/genomes/{w['genome']}/analyses"
    regions = [{'track': ids['genes'], 'type': 'gene'}, {'track': ids['genes'], 'type': 'transposable_element'}]
    asks = {
        'profile': {'kind': 'profile', 'tracks': [ids['wt'], ids['mut']], 'regions': regions,
                    'heatmap': True, 'before': 1000, 'after': 1000, 'body': 1000, 'bin': 100},
        'genes': {'kind': 'genes', 'tracks': [ids['wt'], ids['mut']], 'regions': regions[:1]},
        'correlation': {'kind': 'correlation', 'tracks': [ids['wt'], ids['mut']], 'bin': 1000},
        'compare': {'kind': 'compare', 'a': ids['wt'], 'b': ids['mut'], 'operation': 'subtract', 'bin': 100,
                    'name': 'WT minus MUT'},
    }
    jobs = {}
    for kind, ask in asks.items():          # at most three waiting at once, per person
        jobs[kind] = w['call']('POST', path, 'alice', ask)
        assert jobs[kind].status_code == 202, jobs[kind].get_json()
        _run_queue(w)
    done = {k: w['call']('GET', f"/api/analyses/{r.get_json()['id']}").get_json() for k, r in jobs.items()}
    for kind, job in done.items():
        assert job['state'] == 'done', (kind, job.get('error'))
    assert {'profile.svg', 'heatmap.png', 'matrix.gz'} <= set(done['profile']['result']['files'])
    svg = w['call']('GET', f"/api/analyses/{done['profile']['id']}/files/profile.svg")
    assert svg.status_code == 200 and b'<svg' in svg.data and 'sandbox' in svg.headers['Content-Security-Policy']
    assert w['call']('GET', f"/api/analyses/{done['profile']['id']}/files/profile.svg", 'bob').status_code == 404
    genes = done['genes']['result']
    assert genes['count'] == 3 and genes['preview'][0]['values'][0] == pytest.approx(10.0)
    table = w['call']('GET', f"/api/analyses/{done['genes']['id']}/files/signal-per-feature.csv").get_data(as_text=True)
    assert table.splitlines()[0] == 'name,chrom,start,end,strand,WT,MUT,WT minus MUT'
    assert 'GENE1,Chr1,1001,3000,+,10,4,6' in table
    assert {'correlation.svg', 'pca.svg'} <= set(done['correlation']['result']['files'])
    # The comparison is a new coverage track of alice's lab, shown in the genome's track list.
    new = done['compare']['result']['track_id']
    listed = {t['id']: t for t in w['call']('GET', f"/api/genomes/{w['genome']}/tracks").get_json()}
    assert listed[new]['name'] == 'WT minus MUT' and listed[new]['type'] == 'bigwig'
    means = w['call']('POST', f"/api/genomes/{w['genome']}/region-means", 'alice',
                      {'tracks': [ids['wt'], new], 'chrom': 'Chr1', 'start': 1000, 'end': 3000}).get_json()
    assert [round(t['mean'], 3) for t in means['tracks']] == [10.0, 6.0]


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_the_recipe_makes_the_figure_again(world, tmp_path):
    """recipe.zip holds what is needed to run it again: the commands, by the tracks' download names, and the regions."""
    import io
    import subprocess
    import zipfile
    w, ids = world, world['ids']
    res = w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice',
                    {'kind': 'genes', 'tracks': [ids['wt'], ids['mut']], 'regions': [{'track': ids['genes'], 'type': 'gene'}]})
    _run_queue(w)
    job = res.get_json()['id']
    recipe = w['call']('GET', f'/api/analyses/{job}/recipe').get_json()
    assert recipe['inputs'] == ['MUT.bw', 'WT.bw'] and recipe['regions'] == ['regions.bed']
    assert all(str(w['root']) not in c for c in recipe['commands'])       # no server paths
    blob = w['call']('GET', f'/api/analyses/{job}/files/recipe.zip').data
    z = zipfile.ZipFile(io.BytesIO(blob))
    assert {'README.md', 'run.sh', 'regions.bed', 'name_features.py'} <= set(z.namelist())
    # Run it, in a folder with the tracks under their download names: the same table comes out.
    z.extractall(tmp_path)
    for key, name in (('wt', 'WT.bw'), ('mut', 'MUT.bw')):
        shutil.copy(w['root'] / 'tracks' / f'{key}.bw', tmp_path / name)
    subprocess.run(['sh', 'run.sh'], cwd=tmp_path, check=True, capture_output=True)
    again = (tmp_path / 'signal-per-feature.csv').read_text().splitlines()
    assert again[0] == 'name,chrom,start,end,strand,WT,MUT' and again[1].startswith('GENE1,Chr1,1001,3000,+,10')
    assert w['call']('GET', f'/api/analyses/{job}/recipe', 'bob').status_code == 404
