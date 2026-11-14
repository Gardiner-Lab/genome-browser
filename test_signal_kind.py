"""Methylation tracks are found from their files and measured as a level, whatever their strand signs."""
import gzip
import json
import shutil
import uuid

import numpy as np
import pytest

pyBigWig = pytest.importorskip('pyBigWig')

import analysis
import bigwig_means
import db
import signal_kind
from test_analysis import HAS_DEEPTOOLS, _run_queue, world  # noqa: F401  (the same genome, tracks and members)

GENES = {'Chr1': [(1000, 3000), (6000, 9000)], 'Chr2': [(2000, 4000)]}


def _cg_track(path, inside=0.8, outside=0.1, step=250, signed=True, scale=1.0):
    """CG calls as GEO gives them: each CG a pair of cytosines side by side, the minus strand's with a
    minus sign; ``inside`` methylated over the genes, ``outside`` elsewhere, a pair every ``step``."""
    bw = pyBigWig.open(str(path), 'w')
    bw.addHeader([('Chr1', 20000), ('Chr2', 10000)])
    for chrom, length in (('Chr1', 20000), ('Chr2', 10000)):
        starts, values = [], []
        for s in range(0, length - 1, step):
            v = (inside if any(a <= s < b for a, b in GENES.get(chrom, [])) else outside) * scale
            starts += [s, s + 1]
            values += [v, -v if signed else v]
        bw.addEntries([chrom] * len(starts), starts, ends=[s + 1 for s in starts], values=values)
    bw.close()
    return path


def _dense(path, chain=True):
    """Coverage at single-base resolution: every base its own entry, joined end to end."""
    rng = np.random.default_rng(1)
    bw = pyBigWig.open(str(path), 'w')
    bw.addHeader([('Chr1', 20000)])
    starts = list(range(0, 20000)) if chain else list(range(0, 20000, 7))
    bw.addEntries(['Chr1'] * len(starts), starts, ends=[s + 1 for s in starts],
                  values=rng.uniform(0, 60, len(starts)).tolist())
    bw.close()
    return path


def test_methylation_is_told_from_coverage(tmp_path):
    assert signal_kind.detect(_cg_track(tmp_path / 'cg.bw')) == 'methylation'
    assert signal_kind.detect(_cg_track(tmp_path / 'cg-plus.bw', signed=False)) == 'methylation'
    assert signal_kind.detect(_cg_track(tmp_path / 'cg-percent.bw', scale=100, signed=False)) == 'methylation'
    # Single-base coverage joined into long runs is coverage, whatever its values.
    assert signal_kind.detect(_dense(tmp_path / 'dense.bw')) == 'coverage'
    import test_analysis
    test_analysis._bigwig(tmp_path / 'rna.bw', 10.0)
    assert signal_kind.detect(tmp_path / 'rna.bw') == 'coverage'
    # Values past 100 are never a methylation level.
    assert signal_kind.detect(_cg_track(tmp_path / 'big.bw', scale=1000)) == 'coverage'


def test_the_level_file_drops_the_sign_once(tmp_path):
    cg = _cg_track(tmp_path / 'cg.bw')
    level = signal_kind.level_file(str(cg))
    assert level.endswith('cg.level.bw')
    bw = pyBigWig.open(level)
    try:
        assert bw.intervals('Chr1', 1000, 1002) == ((1000, 1001, pytest.approx(0.8)), (1001, 1002, pytest.approx(0.8)))
        assert bw.header()['minVal'] >= 0
    finally:
        bw.close()
    mtime = (tmp_path / 'cg.level.bw').stat().st_mtime_ns
    assert signal_kind.level_file(str(cg)) == level and (tmp_path / 'cg.level.bw').stat().st_mtime_ns == mtime
    plus = _cg_track(tmp_path / 'plus.bw', signed=False)
    assert signal_kind.level_file(str(plus)) == str(plus)          # nothing to drop: the track itself


def test_a_level_is_the_mean_over_the_cytosines(tmp_path):
    bw = pyBigWig.open(signal_kind.level_file(str(_cg_track(tmp_path / 'cg.bw'))))
    try:
        assert bigwig_means.mean(bw, 'Chr1', 1000, 3000, level=True) == pytest.approx(0.8)
        assert bigwig_means.mean(bw, 'Chr1', 1000, 3000) == pytest.approx(0.8 * 16 / 2000)   # per base: diluted
        assert bigwig_means.mean(bw, 'Chr1', 1002, 1200, level=True) is None                 # no cytosine there
        assert bigwig_means.bin_means(bw, 'Chr1', 0, 1000, 100, level=True) == pytest.approx(
            [0.1, np.nan, 0.1, np.nan, np.nan, 0.1, np.nan, 0.1, np.nan, np.nan], nan_ok=True)
        # The last bin may be shorter.
        assert bigwig_means.bin_means(bw, 'Chr1', 0, 1050, 300, level=True) == pytest.approx(
            [0.1, 0.1, 0.1, 0.8], nan_ok=True)
    finally:
        bw.close()


def test_fill_sets_missing_values_of_some_samples_to_zero(tmp_path):
    header = {'sample_boundaries': [0, 2, 4], 'group_boundaries': [0, 1]}
    m = tmp_path / 'matrix.gz'
    with gzip.open(m, 'wt') as fh:
        fh.write('@' + json.dumps(header) + '\n')
        fh.write('Chr1\t0\t10\tG\t.\t+\tnan\t1\tnan\t2\n')
    bigwig_means.main(['fill', str(m), '--samples', '1', '-o', str(m)])
    with gzip.open(m, 'rt') as fh:
        lines = fh.read().splitlines()
    assert json.loads(lines[0][1:]) == header
    assert lines[1] == 'Chr1\t0\t10\tG\t.\t+\t0\t1\tnan\t2'


def _add_cg(world, name='CG', signal=None):
    path = _cg_track(world['root'] / 'tracks' / f'{name}.bw')
    track = uuid.uuid4().hex
    conn = db.connect(world['root'] / 'genome.db')
    conn.execute('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, file_path, '
                 "created_by, signal) VALUES (?, ?, ?, 'bigwig', 'labA', 'lab_only', 'ready', ?, 'alice', ?)",
                 (track, world['genome'], name, str(path), signal))
    conn.commit()
    conn.close()
    return track


def test_the_track_list_says_which_and_its_lab_can_change_it(world):
    w = world
    cg = _add_cg(w)
    listed = {t['id']: t for t in w['call']('GET', f"/api/genomes/{w['genome']}/tracks").get_json()}
    assert listed[cg]['signal'] == 'methylation' and listed[w['ids']['wt']]['signal'] == 'coverage'
    assert listed[w['ids']['genes']]['signal'] is None
    conn = db.connect(w['root'] / 'genome.db')
    assert tuple(conn.execute('SELECT signal, signal_set_by FROM tracks WHERE id = ?', (cg,)).fetchone()) == \
        ('methylation', 'detected')
    conn.close()
    res = w['call']('PATCH', f'/api/tracks/{cg}', 'alice', {'signal': 'coverage'})
    assert res.status_code == 200 and res.get_json()['signal'] == 'coverage'
    assert w['call']('PATCH', f'/api/tracks/{cg}', 'alice', {'signal': 'reads'}).status_code == 400
    assert w['call']('PATCH', f"/api/tracks/{w['ids']['genes']}", 'alice', {'signal': 'coverage'}).status_code == 400
    assert w['call']('PATCH', f"/api/tracks/{w['ids']['theirs']}", 'alice', {'signal': 'methylation'}).status_code in (403, 404)
    listed = {t['id']: t for t in w['call']('GET', f"/api/genomes/{w['genome']}/tracks").get_json()}
    assert listed[cg]['signal'] == 'coverage'          # kept: not found again from the file


def test_the_bar_chart_gives_a_methylation_level(world):
    w = world
    cg = _add_cg(w)
    means = w['call']('POST', f"/api/genomes/{w['genome']}/region-means", 'alice',
                      {'tracks': [cg, w['ids']['wt']], 'chrom': 'Chr1', 'start': 1000, 'end': 3000}).get_json()
    got = {t['id']: t for t in means['tracks']}
    assert got[cg]['mean'] == pytest.approx(0.8) and got[cg]['max'] == pytest.approx(0.8)
    assert got[cg]['signal'] == 'methylation' and got[w['ids']['wt']]['mean'] == pytest.approx(10.0)


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_every_analysis_measures_methylation_as_a_level(world, tmp_path):
    w, ids = world, world['ids']
    cg = _add_cg(w)
    path = f"/api/genomes/{w['genome']}/analyses"
    genes = [{'track': ids['genes'], 'type': 'gene'}]
    asks = {
        'profile': {'kind': 'profile', 'tracks': [ids['wt'], cg], 'regions': genes, 'heatmap': True,
                    'before': 1000, 'after': 1000, 'body': 1000, 'bin': 100},
        'genes': {'kind': 'genes', 'tracks': [ids['wt'], cg], 'regions': genes},
        'correlation': {'kind': 'correlation', 'tracks': [ids['wt'], cg], 'bin': 1000},
    }
    done = {}
    for kind, ask in asks.items():
        res = w['call']('POST', path, 'alice', ask)
        assert res.status_code == 202, res.get_json()
        _run_queue(w)
        done[kind] = w['call']('GET', f"/api/analyses/{res.get_json()['id']}").get_json()
        assert done[kind]['state'] == 'done', (kind, done[kind].get('error'))
        assert done[kind]['result']['methylation'] == ['CG'] and done[kind]['result']['unsigned'] == ['CG']

    # The table: the genes' level, not the strands cancelled to 0.
    table = w['call']('GET', f"/api/analyses/{done['genes']['id']}/files/signal-per-feature.csv").get_data(as_text=True)
    assert 'GENE1,Chr1,1001,3000,+,10,0.8,9.2' in table

    # The profile: WT's stretches without reads are zero, as chosen; CG's steps without a cytosine left out.
    blob = w['call']('GET', f"/api/analyses/{done['profile']['id']}/files/matrix.gz").data
    lines = gzip.decompress(blob).decode().splitlines()
    bounds = json.loads(lines[0][1:])['sample_boundaries']
    rows = np.array([[float(v) for v in line.split('\t')[6:]] for line in lines[1:]])
    wt, cgm = rows[:, bounds[0]:bounds[1]], rows[:, bounds[1]:bounds[2]]
    assert not np.isnan(wt).any()
    assert np.isnan(cgm).any()
    body = cgm[:, 10:20]                                 # the scaled gene bodies
    assert np.nanmean(body) == pytest.approx(0.8)
    assert done['profile']['result']['coverage_tracks'] == 1

    # The recipe makes the copy without signs and measures from it, by download names.
    recipe = w['call']('GET', f"/api/analyses/{done['genes']['id']}/recipe").get_json()
    joined = '\n'.join(recipe['commands']).replace(' \\\n    ', ' ')     # unwrapped
    assert 'bigwig_means.py unsigned CG.bw -o CG.level.bw' in joined
    assert '--as coverage level' in joined and str(w['root']) not in joined
    assert recipe['inputs'] == ['CG.bw', 'WT.bw']


@pytest.mark.skipif(not HAS_DEEPTOOLS, reason='deepTools is not installed here (it is in the image)')
def test_the_methylation_recipe_runs_again(world, tmp_path):
    import io
    import subprocess
    import zipfile
    w, ids = world, world['ids']
    cg = _add_cg(w)
    res = w['call']('POST', f"/api/genomes/{w['genome']}/analyses", 'alice',
                    {'kind': 'genes', 'tracks': [ids['wt'], cg], 'regions': [{'track': ids['genes'], 'type': 'gene'}]})
    _run_queue(w)
    blob = w['call']('GET', f"/api/analyses/{res.get_json()['id']}/files/recipe.zip").data
    zipfile.ZipFile(io.BytesIO(blob)).extractall(tmp_path)
    shutil.copy(w['root'] / 'tracks' / 'wt.bw', tmp_path / 'WT.bw')
    shutil.copy(w['root'] / 'tracks' / 'CG.bw', tmp_path / 'CG.bw')
    subprocess.run(['sh', 'run.sh'], cwd=tmp_path, check=True, capture_output=True)
    again = (tmp_path / 'signal-per-feature.csv').read_text().splitlines()
    assert again[1].startswith('GENE1,Chr1,1001,3000,+,10,0.8')


def test_a_profile_of_methylation_alone_leaves_missing_out(monkeypatch, tmp_path):
    runs = []
    cg = _cg_track(tmp_path / 'cg.bw')
    monkeypatch.setattr(analysis, 'region_bed', lambda *a, **k: 5)
    monkeypatch.setattr(analysis, '_run', lambda args, cwd, what, **step: runs.append(args))
    monkeypatch.setattr(analysis, '_features_drawn', lambda matrix: [5])
    monkeypatch.setattr(analysis, '_off_chromosome_ends', lambda *a: [])
    analysis._RECIPE.clear()
    out = analysis.run_profile({'regions': [{'track': 'g', 'type': 'gene'}], 'heatmap': True},
                               [{'file_path': str(cg), 'name': 'CG', 'signal': 'methylation'}],
                               {'g': {'file_path': 'g.gff', 'type': 'gff', 'name': 'Genes'}}, tmp_path)
    compute = next(a for a in runs if a[0] == 'computeMatrix')
    assert '--missingDataAsZero' not in compute and str(tmp_path / 'cg.level.bw') in compute
    heat = next(a for a in runs if a[0] == 'plotHeatmap')
    assert heat[heat.index('--missingDataColor') + 1] == '0.85'
    profile = next(a for a in runs if a[0] == 'plotProfile')
    assert profile[profile.index('--yAxisLabel') + 1] == 'methylation level'
    assert out['coverage_tracks'] == 0 and out['methylation'] == ['CG']
