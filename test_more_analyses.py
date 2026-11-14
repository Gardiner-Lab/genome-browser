"""Fingerprints, profiles centred on a middle or a summit, k-means clusters as lists and tracks, and
methylation against expression."""
import csv
import io
import shutil
import subprocess
import uuid
import zipfile

import numpy as np
import pytest

pyBigWig = pytest.importorskip('pyBigWig')

import analysis
import db
import plots
from test_analysis import HAS_DEEPTOOLS, _run_queue, world  # noqa: F401  (the same genome, tracks and members)


# --- the pieces --------------------------------------------------------------------------------

def test_fingerprint_numbers():
    even = plots.fingerprint_metrics([5.0] * 1000)
    assert even['zero_bins'] == 0 and even['top1'] == pytest.approx(0.01) and even['auc'] == pytest.approx(0.5, abs=0.01)
    peaks = plots.fingerprint_metrics([0.0] * 990 + [100.0] * 10)
    assert peaks['zero_bins'] == pytest.approx(0.99) and peaks['top1'] == pytest.approx(1.0)
    assert peaks['elbow'] == pytest.approx(0.99) and peaks['auc'] < 0.01
    assert plots.fingerprint_metrics([0.0] * 10)['top1'] is None          # no signal: nothing to say
    x, y = plots.fingerprint_curve([3, float('nan'), 1])
    assert list(y) == pytest.approx([0, 0.25, 1.0])


def test_expression_groups_and_their_names():
    expr = [0, 5, 1, None, 0, 9, 2, 7, 3, 4, 8, 6, 10]
    group = plots.expression_groups(expr, 5)
    assert group[0] == group[4] == 0 and group[3] is None
    # Ten expressed, five groups: two in each, lowest first.
    by = {g: sorted(e for e, x in zip(expr, group) if x == g) for g in range(1, 6)}
    assert by == {1: [1, 2], 2: [3, 4], 3: [5, 6], 4: [7, 8], 5: [9, 10]}
    assert plots.group_labels(5) == ['no reads', 'lowest 20%', '20–40%', '40–60%', '60–80%', 'highest 20%']
    assert plots.group_labels(4)[1] == 'lowest 25%' and plots.group_labels(4)[-1] == 'highest 25%'


def test_summits_of_a_narrowpeak_file(tmp_path):
    peaks = tmp_path / 'peaks.narrowPeak'
    peaks.write_text('Chr1\t100\t300\tp1\t500\t.\t9\t-1\t-1\t50\n'      # summit 50 bp in
                     'Chr1\t1000\t1200\tp2\t500\t.\t9\t-1\t-1\t-1\n'    # none called: the middle
                     'Chr2\t10\t20\n'                                    # plain BED: the middle
                     'track name=x\n')
    out = tmp_path / 's.bed'
    assert analysis.summit_bed(str(peaks), None, out) == (3, 1)
    rows = [line.split('\t')[:4] for line in out.read_text().splitlines()]
    assert rows == [['Chr1', '150', '151', 'p1'], ['Chr1', '1100', '1101', 'p2'], ['Chr2', '15', '16', 'Chr2:11']]
    assert analysis.summit_bed(str(peaks), ['P2'], out) == (1, 0)


def test_the_stretch_before_each_gene_by_its_strand(tmp_path):
    genes = tmp_path / 'g.bed'
    genes.write_text('Chr1\t5000\t6000\tA\t0\t+\nChr1\t5000\t6000\tB\t0\t-\nChr1\t300\t900\tC\t0\t+\n'
                     'Chr1\t9500\t9900\tD\t0\t-\n')
    out = tmp_path / 'u.bed'
    analysis.upstream_bed(genes, out, 1000, {'Chr1': 10000})
    assert [line.split('\t')[:4] for line in out.read_text().splitlines()] == [
        ['Chr1', '4000', '5000', 'A'], ['Chr1', '6000', '7000', 'B'],
        ['Chr1', '0', '300', 'C'], ['Chr1', '9900', '10000', 'D']]


def test_clusters_become_lists(tmp_path):
    (tmp_path / 'heatmap-regions.bed').write_text(
        '#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStart\tdeepTools_group\n'
        'Chr1\t10\t20\tG2\t0\t-\t10\t20\t0\t1\t10\t0\tcluster_2\n'
        'Chr1\t0\t5\tG1\t0\t+\t0\t5\t0\t1\t5\t0\tcluster_1\n'
        'Chr2\t0\t5\tG3\t0\t+\t0\t5\t0\t1\t5\t0\tcluster_1\n')
    found = analysis._clusters(tmp_path)
    assert [(c['cluster'], c['count'], c['first']) for c in found] == [(1, 2, ['G1', 'G3']), (2, 1, ['G2'])]
    assert (tmp_path / 'cluster-1.txt').read_text() == 'G1\nG3\n'
    assert analysis.cluster_bed(tmp_path, 2) == 'Chr1\t10\t20\tG2\t0\t-\n'


# --- the whole way, with deepTools ------------------------------------------------------------------

N_GENES = 36


def _genes_world(w):
    """36 genes of 300 bp on Chr1, every 500 bp; an RNA-seq track whose level rises gene by gene (a third
    of them unread), a CG track methylated less the more a gene is read, and peaks with summits."""
    root = w['root'] / 'tracks'
    gff = root / 'many.gff'
    with open(gff, 'w') as fh:
        fh.write('##gff-version 3\n')
        for i in range(N_GENES):
            s = 1000 + i * 500
            fh.write(f'Chr1\tx\tgene\t{s + 1}\t{s + 300}\t.\t{"+-"[i % 2]}\t.\tID=g{i};Name=G{i}\n')
    rna = root / 'rna.bw'
    bw = pyBigWig.open(str(rna), 'w')
    bw.addHeader([('Chr1', 20000), ('Chr2', 10000)])
    starts = [1000 + i * 500 for i in range(N_GENES) if i % 3]
    bw.addEntries(['Chr1'] * len(starts), starts, ends=[s + 300 for s in starts],
                  values=[float((s - 1000) // 500) for s in starts])
    bw.close()
    cg = root / 'cg.bw'
    bw = pyBigWig.open(str(cg), 'w')
    bw.addHeader([('Chr1', 20000), ('Chr2', 10000)])
    cs, vs = [], []
    for i in range(N_GENES):
        s = 1000 + i * 500
        level = 0.9 if i % 3 == 0 else 0.9 - i / N_GENES
        for c in range(s, s + 300, 20):
            cs += [c, c + 1]
            vs += [level, -level]
    bw.addEntries(['Chr1'] * len(cs), cs, ends=[c + 1 for c in cs], values=vs)
    bw.close()
    peaks = root / 'peaks.narrowPeak'
    peaks.write_text(''.join(f'Chr1\t{1000 + i * 500}\t{1300 + i * 500}\tpeak{i}\t100\t.\t5\t-1\t-1\t{100 + i}\n'
                             for i in range(N_GENES)))
    ids = {k: uuid.uuid4().hex for k in ('many', 'rna', 'cg', 'peaks')}
    conn = db.connect(w['root'] / 'genome.db')
    conn.executemany('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, file_path, '
                     "created_by) VALUES (?, ?, ?, ?, 'labA', 'lab_only', 'ready', ?, 'alice')",
                     [(ids['many'], w['genome'], 'Many genes', 'gff', str(gff)),
                      (ids['rna'], w['genome'], 'RNA', 'bigwig', str(rna)),
                      (ids['cg'], w['genome'], 'CG', 'bigwig', str(cg)),
                      (ids['peaks'], w['genome'], 'Peaks', 'bed', str(peaks))])
    conn.commit()
    conn.close()
    return ids


def _done(w, ask):
    res = w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice', ask)
    assert res.status_code == 202, res.get_json()
    _run_queue(w)
    job = w['call']('GET', f"/api/analyses/{res.get_json()['id']}").get_json()
    assert job['state'] == 'done', job.get('error')
    return job


def _file(w, job, name, user='alice'):
    return w['call']('GET', f"/api/analyses/{job['id']}/files/{name}", user)


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_fingerprint_of_coverage_tracks(world):
    w = world
    ids = _genes_world(w)
    job = _done(w, {'kind': 'fingerprint', 'tracks': [w['ids']['wt'], ids['rna'], ids['cg']], 'bin': 500})
    r = job['result']
    assert r['skipped_methylation'] == ['CG'] and [m['track'] for m in r['metrics']] == ['WT', 'RNA']
    wt, rna = r['metrics']
    assert wt['zero_bins'] == 0 and rna['zero_bins'] > 0.5          # RNA only on two thirds of the genes
    assert rna['top1'] > wt['top1']
    assert b'<svg' in _file(w, job, 'fingerprint.svg').data
    table = _file(w, job, 'fingerprint-metrics.tab').get_data(as_text=True).splitlines()
    assert table[0].startswith('track\tbins\t') and table[1].startswith('WT\t60\t')
    alone = w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice', {'kind': 'fingerprint', 'tracks': [ids['cg']]})
    _run_queue(w)
    failed = w['call']('GET', f"/api/analyses/{alone.get_json()['id']}").get_json()
    assert failed['state'] == 'failed' and 'methylation' in failed['error']


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_profiles_on_summits_and_middles_and_clusters_kept_as_tracks(world):
    w = world
    ids = _genes_world(w)
    peaks = _done(w, {'kind': 'profile', 'tracks': [ids['rna']], 'mode': 'summit', 'before': 1000, 'after': 1000,
                      'bin': 50, 'regions': [{'track': ids['peaks'], 'type': 'region'}]})
    assert peaks['result']['summits'] == {'peaks': N_GENES, 'with_summit': N_GENES}
    middle = _done(w, {'kind': 'profile', 'tracks': [ids['rna']], 'mode': 'center', 'before': 1000, 'after': 1000,
                       'bin': 50, 'regions': [{'track': ids['many'], 'type': 'gene'}]})
    assert middle['result']['mode'] == 'center'
    recipe = w['call']('GET', f"/api/analyses/{middle['id']}/recipe").get_json()
    assert '--referencePoint center' in ' '.join(recipe['commands'][0].replace('\\', ' ').split())

    job = _done(w, {'kind': 'profile', 'tracks': [ids['rna'], ids['cg']], 'heatmap': True, 'kmeans': 2,
                    'before': 500, 'after': 500, 'body': 500, 'bin': 50,
                    'regions': [{'track': ids['many'], 'type': 'gene'}]})
    r = job['result']
    assert [c['cluster'] for c in r['clusters']] == [1, 2]
    assert sum(c['count'] for c in r['clusters']) == N_GENES
    assert {'clusters.csv', 'cluster-1.txt', 'cluster-2.txt'} <= set(r['files'])
    listed = _file(w, job, 'cluster-1.txt').get_data(as_text=True).split()
    assert len(listed) == r['clusters'][0]['count'] and all(n.startswith('G') for n in listed)
    rows = list(csv.DictReader(io.StringIO(_file(w, job, 'clusters.csv').get_data(as_text=True))))
    assert {row['cluster'] for row in rows} == {'1', '2'} and rows[0]['chrom'] == 'Chr1'
    assert _file(w, job, 'cluster-1.txt', 'bob').status_code == 404

    path = f"/api/analyses/{job['id']}/clusters/1/track"
    assert w['call']('POST', path, 'bob', {}).status_code == 404
    assert w['call']('POST', f"/api/analyses/{job['id']}/clusters/9/track", 'alice', {}).status_code == 404
    made = w['call']('POST', path, 'alice', {})
    assert made.status_code == 201, made.get_json()
    track = made.get_json()
    assert track['type'] == 'bed' and track['features'] == r['clusters'][0]['count'] and 'cluster 1 of 2' in track['name']
    tracks = {t['id']: t for t in w['call']('GET', f"/api/genomes/{w['genome']}/tracks").get_json()}
    assert tracks[track['id']]['owner_lab'] == 'labA'
    # And it can be measured again, as features.
    again = _done(w, {'kind': 'genes', 'tracks': [ids['rna']], 'regions': [{'track': track['id'], 'type': 'region'}]})
    assert again['result']['count'] == r['clusters'][0]['count']


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_methylation_against_expression(world, tmp_path):
    w = world
    ids = _genes_world(w)
    regions = [{'track': ids['many'], 'type': 'gene'}]
    bad = w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice',
                    {'kind': 'methexpr', 'tracks': [ids['cg']], 'regions': regions})
    assert bad.status_code == 400 and 'expression' in bad.get_json()['error']
    job = _done(w, {'kind': 'methexpr', 'tracks': [ids['cg']], 'expression': ids['rna'], 'regions': regions,
                    'groups': 4})
    r = job['result']
    assert r['count'] == N_GENES and r['expression'] == 'RNA' and r['methylation'] == ['CG']
    assert [g['count'] for g in r['groups']] == [12, 6, 6, 6, 6]
    cg = r['tracks'][0]
    assert cg['rho'] < -0.9 and cg['rho_expressed'] < -0.9 and cg['n'] == N_GENES
    means = cg['mean_by_group']
    assert means[0] == pytest.approx(0.9) and means[1] > means[2] > means[3] > means[4]
    assert {'methexpr.svg', 'groups-profile.svg', 'methylation-vs-expression.csv', 'recipe.zip'} <= set(r['files'])
    table = _file(w, job, 'methylation-vs-expression.csv').get_data(as_text=True).splitlines()
    assert table[0] == 'name,chrom,start,end,strand,expression,CG,expression_group'
    assert table[1] == 'G0,Chr1,1001,1300,+,0,0.9,no reads'
    assert table[-1].startswith('G35,Chr1,18501,18800,-,35,') and table[-1].endswith(',highest 25%')

    # The promoter: the stretch before each gene holds no cytosines here, so no level at all.
    up = _done(w, {'kind': 'methexpr', 'tracks': [ids['cg']], 'expression': ids['rna'], 'regions': regions,
                   'region': 'up1000', 'group_profile': False})
    assert up['result']['where'] == 'over the 1 kb before each gene' and 'groups-profile.svg' not in up['result']['files']
    assert up['result']['tracks'][0]['n'] < N_GENES

    # The recipe makes the same table again.
    z = zipfile.ZipFile(io.BytesIO(_file(w, job, 'recipe.zip').data))
    assert {'plots.py', 'bigwig_means.py', 'regions.bed', 'run.sh'} <= set(z.namelist())
    z.extractall(tmp_path)
    shutil.copy(w['root'] / 'tracks' / 'rna.bw', tmp_path / 'RNA.bw')
    shutil.copy(w['root'] / 'tracks' / 'cg.bw', tmp_path / 'CG.bw')
    subprocess.run(['sh', 'run.sh'], cwd=tmp_path, check=True, capture_output=True)
    assert (tmp_path / 'methylation-vs-expression.csv').read_text().splitlines() == table
