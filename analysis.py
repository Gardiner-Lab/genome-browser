"""
Analyses of the coverage tracks on a genome, drawn as figures: the Analyse panel.

    profile      the average signal of chosen tracks across a set of features (genes,
                 transposons, the regions of a BED track): a metagene profile, and with
                 it, if asked, a heatmap of every feature, sorted or clustered
    compare      one coverage track against another, A minus B or their (log2) ratio,
                 kept as a new coverage track of the member's lab
    correlation  how alike tracks are along the genome (replicates, mutants):
                 a correlation heatmap and a PCA
    genes        each feature's mean signal in each track, as a table (CSV)

They are drawn by deepTools (installed in the image; Ramirez et al. 2016), which take
from seconds to minutes on a whole genome, so each is a job: the web app checks the
member may see the genome and every track, and queues it (``analysis_jobs``); the
genome worker (worker.py) runs the queue, one at a time, between conversions; the page
asks after it until it is done, and fetches the files it made.

Two quick ones answer at once, without a job: the tracks' means over a region (for a bar
chart of the selection) and nothing else; the synteny dot plot is drawn by the page itself
from the links track.

A job's files are kept in ``analysis/<id>/`` for KEEP_DAYS, for the member who asked.
"""

import csv
import json
import math
import os
import re
import shutil
import sqlite3
import subprocess
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import feature_index
import regions

KINDS = ('profile', 'compare', 'correlation', 'genes', 'motif')
MAX_TRACKS = 12
MAX_REGION_SETS = 4
KEEP_DAYS = 7
COLORMAPS = ('viridis', 'magma', 'RdBu_r', 'Blues', 'Reds')
TIMEOUT_SECONDS = 1800
THREADS = int(os.environ.get('ANALYSIS_THREADS', '2'))
OPERATIONS = {'subtract': 'A minus B', 'log2': 'log2(A / B)', 'ratio': 'A / B', 'add': 'A plus B', 'mean': 'mean of A and B'}

SCHEMA = """
CREATE TABLE IF NOT EXISTS analysis_jobs (
  id          TEXT PRIMARY KEY,
  genome_id   TEXT NOT NULL,
  member      TEXT NOT NULL,
  owner_lab   TEXT,
  kind        TEXT NOT NULL,
  params      TEXT NOT NULL,
  state       TEXT NOT NULL DEFAULT 'queued',
  error       TEXT,
  result      TEXT,
  created_at  TEXT NOT NULL DEFAULT (datetime('now')),
  started_at  TEXT,
  finished_at TEXT
);
CREATE INDEX IF NOT EXISTS analysis_jobs_queue ON analysis_jobs (state, created_at);
"""


class AnalysisError(Exception):
    """Something the member can act on, as a sentence."""


# The commands the job running now has run, in order (the worker runs one job at a time).
_RECIPE: List[List[str]] = []


def ensure_schema(conn) -> None:
    conn.executescript(SCHEMA)


def job_dir(data_root: Path, job_id: str) -> Path:
    return Path(data_root) / 'analysis' / job_id


# --- the queue -----------------------------------------------------------------------

def enqueue(conn, genome_id: str, member: str, owner_lab: Optional[str], kind: str, params: dict) -> str:
    ensure_schema(conn)
    job_id = uuid.uuid4().hex
    conn.execute('INSERT INTO analysis_jobs (id, genome_id, member, owner_lab, kind, params) VALUES (?, ?, ?, ?, ?, ?)',
                 (job_id, genome_id, member, owner_lab, kind, json.dumps(params)))
    conn.commit()
    return job_id


def claim_next(conn) -> Optional[sqlite3.Row]:
    ensure_schema(conn)
    conn.execute('BEGIN IMMEDIATE')
    try:
        row = conn.execute("SELECT * FROM analysis_jobs WHERE state = 'queued' ORDER BY created_at LIMIT 1").fetchone()
        if row is not None:
            conn.execute("UPDATE analysis_jobs SET state = 'running', started_at = datetime('now') WHERE id = ?",
                         (row['id'],))
        conn.execute('COMMIT')
    except Exception:
        conn.execute('ROLLBACK')
        raise
    return row


def recover(conn) -> int:
    """A job left running by a restart is run again."""
    ensure_schema(conn)
    n = conn.execute("UPDATE analysis_jobs SET state = 'queued', started_at = NULL WHERE state = 'running'").rowcount
    conn.commit()
    return n


def tidy(conn, data_root: Path) -> None:
    """Jobs older than KEEP_DAYS go, with their files."""
    ensure_schema(conn)
    old = [r[0] for r in conn.execute(
        f"SELECT id FROM analysis_jobs WHERE created_at < datetime('now', '-{KEEP_DAYS} days')")]
    for job_id in old:
        shutil.rmtree(job_dir(data_root, job_id), ignore_errors=True)
        conn.execute('DELETE FROM analysis_jobs WHERE id = ?', (job_id,))
    conn.commit()


# --- inputs ---------------------------------------------------------------------------

def feature_types(path: str, kind: str) -> List[Tuple[str, int]]:
    idx = feature_index.ensure(path, regions.feature_parser(kind))
    conn = sqlite3.connect(f'file:{idx}?mode=ro', uri=True)
    try:
        return [(t or 'feature', n) for t, n in conn.execute(
            'SELECT type, count(*) FROM features GROUP BY type ORDER BY 2 DESC')]
    finally:
        conn.close()


def region_bed(path: str, kind: str, ftype: Optional[str], names: Optional[List[str]], out: Path,
               chroms: Optional[set] = None) -> int:
    """The features of one type (or all, for BED) as a 6-column BED; returns how many."""
    idx = feature_index.ensure(path, regions.feature_parser(kind))
    conn = sqlite3.connect(f'file:{idx}?mode=ro', uri=True)
    wanted = {n.strip().casefold() for n in (names or []) if n.strip()}
    n = 0
    try:
        sql = 'SELECT chrom, start, "end", name, strand FROM features'
        args: list = []
        if ftype:
            sql += ' WHERE type = ?'
            args.append(ftype)
        sql += ' ORDER BY chrom, start'
        with open(out, 'w', encoding='utf-8') as fh:
            for chrom, start, end, name, strand in conn.execute(sql, args):
                if wanted and (name or '').casefold() not in wanted:
                    continue
                if chroms is not None and chrom not in chroms:
                    continue
                if end <= start:
                    continue
                label = re.sub(r'\s+', '_', name or f'{chrom}:{start + 1}')
                fh.write(f"{chrom}\t{start}\t{end}\t{label}\t0\t{'-' if strand == -1 else '+'}\n")
                n += 1
    finally:
        conn.close()
    return n


def _label(name: str) -> str:
    """A track name as deepTools shows it in a legend: no characters that break its arguments."""
    return re.sub(r'\s+', ' ', name).strip()[:60] or 'track'


def _run(args: List[str], cwd: Path, what: str) -> None:
    _RECIPE.append(list(args))
    env = dict(os.environ, MPLBACKEND='Agg', MPLCONFIGDIR=str(cwd / '.mpl'))
    try:
        done = subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS, env=env)
    except subprocess.TimeoutExpired:
        raise AnalysisError(f'{what} took too long: try fewer tracks or features.')
    if done.returncode != 0:
        detail = (done.stderr or done.stdout or '').strip().splitlines()[-2:]
        raise AnalysisError(f'{what} failed' + (': ' + ' '.join(detail) if detail else '.'))


# --- the analyses ----------------------------------------------------------------------

def run_profile(params: dict, tracks: List[dict], region_tracks: Dict[str, dict], work: Path) -> dict:
    sets, labels = [], []
    for i, rs in enumerate(params.get('regions') or []):
        rt = region_tracks.get(rs.get('track'))
        if rt is None:
            continue
        bed = work / f'regions-{i + 1}.bed'
        n = region_bed(rt['file_path'], rt['type'], rs.get('type') or None, rs.get('names'), bed)
        if n:
            sets.append(bed)
            labels.append(f"{rs.get('type') or rt['name']} ({n:,})")
    if not sets:
        raise AnalysisError('None of the chosen features were found: choose another type, or check the names.')
    mode = params.get('mode') if params.get('mode') in ('scale', 'start', 'end') else 'scale'
    before, after = (int(_clamp(params.get(k), 0, 20000, 2000)) for k in ('before', 'after'))
    body = int(_clamp(params.get('body'), 500, 20000, 3000))
    binsize = int(_clamp(params.get('bin'), 10, 5000, 100))
    # A coverage track has no entries where nothing was read, so without --missingDataAsZero those
    # bins are NaN: plotHeatmap paints them black and plotProfile leaves them out of the mean.
    args = ['computeMatrix', 'scale-regions' if mode == 'scale' else 'reference-point',
            '-S', *[t['file_path'] for t in tracks], '-R', *[str(s) for s in sets],
            '--samplesLabel', *[_label(t['name']) for t in tracks],
            '-b', str(before), '-a', str(after), '-bs', str(binsize), '-p', str(THREADS), '--skipZeros',
            '--missingDataAsZero', '-o', 'matrix.gz', '--quiet']
    if mode == 'scale':
        args += ['-m', str(body)]
    else:
        args += ['--referencePoint', 'TSS' if mode == 'start' else 'TES']
    _run(args, work, 'Computing the profile')
    # Region-set labels as chosen (computeMatrix names them after the files).
    title = (params.get('title') or '').strip()[:100]
    start_label, end_label = ('start', 'end')
    common = ['--regionsLabel', *labels, '--startLabel', start_label, '--endLabel', end_label,
              '--refPointLabel', 'start' if mode == 'start' else 'end']
    profile_args = ['plotProfile', '-m', 'matrix.gz', '-o', 'profile.svg', '--perGroup', *common,
                    '--plotTitle', title or 'Mean signal', '--yAxisLabel', 'mean signal']
    _run(profile_args, work, 'Drawing the profile')
    files = ['profile.svg', 'matrix.gz']
    if params.get('heatmap'):
        k = int(_clamp(params.get('kmeans'), 0, 8, 0))
        colormap = str(params.get('colormap') or 'viridis')
        if colormap not in COLORMAPS:
            colormap = 'viridis'
        heat_args = ['plotHeatmap', '-m', 'matrix.gz', '-o', 'heatmap.png', '--dpi', '110', *common,
                     '--colorMap', colormap, '--heatmapHeight', '18', '--plotTitle', title or 'Signal per feature']
        if k >= 2:
            heat_args += ['--kmeans', str(k)]
        _run(heat_args, work, 'Drawing the heatmap')
        files.insert(1, 'heatmap.png')
    return {'files': files, 'regions': labels}


def run_compare(params: dict, a: dict, b: dict, work: Path) -> dict:
    op = params.get('operation') if params.get('operation') in OPERATIONS else 'subtract'
    binsize = int(_clamp(params.get('bin'), 10, 10000, 100))
    args = ['bigwigCompare', '-b1', a['file_path'], '-b2', b['file_path'], '--operation', op,
            '-bs', str(binsize), '-p', str(THREADS), '-o', 'compare.bw']
    if op in ('log2', 'ratio'):
        args += ['--pseudocount', str(_clamp(params.get('pseudocount'), 0.0001, 1000, 1))]
    _run(args, work, 'Comparing the tracks')
    return {'files': ['compare.bw'], 'operation': OPERATIONS[op]}


def run_correlation(params: dict, tracks: List[dict], work: Path) -> dict:
    binsize = int(_clamp(params.get('bin'), 1000, 1_000_000, 10000))
    method = params.get('method') if params.get('method') in ('pearson', 'spearman') else 'pearson'
    labels = [_label(t['name']) for t in tracks]
    _run(['multiBigwigSummary', 'bins', '-b', *[t['file_path'] for t in tracks], '--labels', *labels,
          '-bs', str(binsize), '-p', str(THREADS), '-o', 'summary.npz', '--outRawCounts', 'bins.tab'],
         work, 'Summarising the tracks')
    _run(['plotCorrelation', '-in', 'summary.npz', '-c', method, '-p', 'heatmap', '--plotNumbers',
          '--skipZeros', '--removeOutliers', '--colorMap', 'RdYlBu_r', '-o', 'correlation.svg',
          '--plotTitle', f'{method.capitalize()} correlation, {binsize:,} bp bins',
          '--outFileCorMatrix', 'correlation.tab'], work, 'Drawing the correlation')
    _run(['plotPCA', '-in', 'summary.npz', '-o', 'pca.svg', '--plotTitle', 'Principal components'],
         work, 'Drawing the PCA')
    return {'files': ['correlation.svg', 'pca.svg', 'correlation.tab', 'bins.tab']}


def run_genes(params: dict, tracks: List[dict], region_tracks: Dict[str, dict], work: Path) -> dict:
    rs = (params.get('regions') or [{}])[0]
    rt = region_tracks.get(rs.get('track'))
    if rt is None:
        raise AnalysisError('Choose the features to measure.')
    bed = work / 'regions.bed'
    n = region_bed(rt['file_path'], rt['type'], rs.get('type') or None, rs.get('names'), bed)
    if not n:
        raise AnalysisError('None of the chosen features were found.')
    labels = [_label(t['name']) for t in tracks]
    _run(['multiBigwigSummary', 'BED-file', '-b', *[t['file_path'] for t in tracks], '--BED', str(bed),
          '--labels', *labels, '-p', str(THREADS), '-o', 'genes.npz', '--outRawCounts', 'raw.tab'],
         work, 'Measuring each feature')
    # deepTools keeps coordinates only: the names come back from the BED, by position.
    named: Dict[Tuple[str, int, int], Tuple[str, str]] = {}
    with open(bed, encoding='utf-8') as fh:
        for line in fh:
            c, s, e, name, _, strand = line.rstrip('\n').split('\t')
            named.setdefault((c, int(s), int(e)), (name, strand))
    rows = []
    with open(work / 'raw.tab', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            c, s, e = f[0], int(f[1]), int(f[2])
            values = [None if v in ('nan', '') else float(v) for v in f[3:]]
            name, strand = named.get((c, s, e), (f'{c}:{s + 1}', '.'))
            rows.append((name, c, s + 1, e, strand, values))
    rows.sort(key=lambda r: (r[1], r[2]))
    two = len(tracks) == 2
    with open(work / 'signal-per-feature.csv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'chrom', 'start', 'end', 'strand'] + [t['name'] for t in tracks]
                   + ([f"{tracks[0]['name']} minus {tracks[1]['name']}"] if two else []))
        for name, c, s, e, strand, values in rows:
            extra = []
            if two:
                extra = ['' if None in values else f'{values[0] - values[1]:.6g}']
            w.writerow([name, c, s, e, strand] + ['' if v is None else f'{v:.6g}' for v in values] + extra)
    # The first rows to show: by the largest change between two tracks, else by the first track.
    def key(r):
        v = r[5]
        if two:
            return -abs(v[0] - v[1]) if None not in v else 0
        return -(v[0] or 0)
    top = sorted(rows, key=key)[:50]
    preview = [{'name': r[0], 'chrom': r[1], 'start': r[2], 'end': r[3],
                'values': [None if v is None else round(v, 6) for v in r[5]]} for r in top]
    return {'files': ['signal-per-feature.csv'], 'count': len(rows), 'preview': preview,
            'columns': [t['name'] for t in tracks]}


def _clamp(value, low, high, default):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(v):
        return default
    return max(low, min(high, v))


# --- running a job (the worker) -------------------------------------------------------

def run_job(conn, data_root: Path, job) -> None:
    """Run one claimed job and record what came of it."""
    params = json.loads(job['params'])
    _RECIPE.clear()
    work = job_dir(data_root, job['id'])
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    started = time.time()
    try:
        ids = list(params.get('tracks') or []) + [r.get('track') for r in params.get('regions') or []]
        ids = [i for i in ids + [params.get('a'), params.get('b')] if i]
        rows = {}
        if ids:
            rows = {r['id']: dict(r) for r in conn.execute(
                f'SELECT id, name, type, file_path, track_status FROM tracks WHERE id IN ({",".join("?" * len(ids))})',
                ids)}
        coverage = [rows[i] for i in params.get('tracks') or [] if i in rows and rows[i]['type'] == 'bigwig'
                    and rows[i]['track_status'] == 'ready']
        region_tracks = {i: rows[i] for i in rows if rows[i]['type'] in ('gff', 'bed')}
        kind = job['kind']
        if kind == 'profile':
            if not coverage:
                raise AnalysisError('Choose at least one coverage track.')
            result = run_profile(params, coverage, region_tracks, work)
        elif kind == 'correlation':
            if len(coverage) < 2:
                raise AnalysisError('Choose at least two coverage tracks to compare.')
            result = run_correlation(params, coverage, work)
        elif kind == 'genes':
            if not coverage:
                raise AnalysisError('Choose at least one coverage track.')
            result = run_genes(params, coverage, region_tracks, work)
        elif kind == 'motif':
            # Marking a sequence on a big genome: too slow to hold a request open, so it waits here
            # with the rest. It makes its own BED track (tracks.py), as 'compare' makes a BigWig one.
            import tracks as tracks_module
            result = tracks_module.run_motif_job(conn, data_root, job, params)
            result['files'] = []
        elif kind == 'compare':
            a, b = rows.get(params.get('a')), rows.get(params.get('b'))
            if not a or not b or a['type'] != 'bigwig' or b['type'] != 'bigwig':
                raise AnalysisError('Choose two coverage tracks.')
            result = run_compare(params, a, b, work)
            result['track_id'] = _add_track(conn, data_root, job, params, a, b, work / 'compare.bw', result['operation'])
            result['files'] = []
        else:
            raise AnalysisError(f'Unknown analysis {kind}.')
        result['seconds'] = round(time.time() - started)
        if kind not in ('compare', 'motif'):
            write_recipe(work, job, params, rows, list(_RECIPE))
            result['files'].append('recipe.zip')
        conn.execute("UPDATE analysis_jobs SET state = 'done', result = ?, finished_at = datetime('now') WHERE id = ?",
                     (json.dumps(result), job['id']))
    except Exception as exc:   # every failure is the member's to see, never the worker's end
        message = str(exc) if isinstance(exc, AnalysisError) else f'The analysis failed: {exc}'
        conn.execute("UPDATE analysis_jobs SET state = 'failed', error = ?, finished_at = datetime('now') WHERE id = ?",
                     (message[:1000], job['id']))
    finally:
        shutil.rmtree(work / '.mpl', ignore_errors=True)
        conn.commit()


def _deeptools_version() -> str:
    try:
        out = subprocess.run(['computeMatrix', '--version'], capture_output=True, text=True, timeout=30)
        return (out.stdout or out.stderr).strip().split()[-1]
    except Exception:
        return '3.5.1'


def write_recipe(work: Path, job, params: dict, rows: Dict[str, dict], commands: List[List[str]]) -> None:
    """
    recipe.zip: everything needed to make the figure again outside the suite: README.md (what it is, the
    inputs, the software), run.sh (the commands as run, with the tracks named by their download names),
    the region files used, and for a table, the script that adds the features' names.
    """
    import zipfile
    from tracks import _export_filename
    names = {}
    for row in rows.values():
        names[row['file_path']] = _export_filename(row['name'], row['type'], row['file_path'])
    regions = sorted(p.name for p in work.glob('regions*.bed'))

    def portable(arg: str) -> str:
        if arg in names:
            return names[arg]
        if arg.startswith(str(work)):
            return Path(arg).name
        return arg

    def shell(args: List[str]) -> str:
        parts = []
        for a in (portable(x) for x in args):
            parts.append(a if re.fullmatch(r'[A-Za-z0-9_./:=+-]+', a) else "'" + a.replace("'", "'\\''") + "'")
        out, line = [], ''
        for p in parts:
            if len(line) + len(p) > 90 and line:
                out.append(line + ' \\')
                line = '    ' + p
            else:
                line = (line + ' ' + p) if line else p
        out.append(line)
        return '\n'.join(out)

    version = _deeptools_version()
    inputs = sorted(set(names[p] for p in names if p in {a for c in commands for a in c}))
    run = ['#!/bin/sh', '# Makes the figure again: run in a folder holding the files below.', 'set -e', '']
    run += [shell(c) + '\n' for c in commands]
    if job['kind'] == 'genes':
        run.append('python3 name_features.py   # adds the features\u2019 names: signal-per-feature.csv')
    readme = [f"# {job['kind'].capitalize()}, from the Trellis Tools Genome Browser", '',
              f"Made {time.strftime('%Y-%m-%d %H:%M')} UTC by {job['member']}.", '',
              '## Inputs', '',
              'Download each track from the Genome Browser (Export, in its settings) and keep its name:', '']
    readme += [f'- `{n}`' for n in inputs]
    if regions:
        readme += ['', 'The features used, as BED files (included here):', ''] + [f'- `{r}`' for r in regions]
    readme += ['', '## Software', '',
               f'- deepTools {version} (`pip install deeptools=={version}`, or `conda install -c bioconda deeptools={version}`)',
               '- Python 3 with matplotlib (comes with deepTools)', '',
               '## Run', '', 'In a folder with the inputs and these files:', '', '```sh', 'sh run.sh', '```', '',
               'The commands, as they were run:', '', '```sh'] + [shell(c) for c in commands] + ['```', '',
               'Cite deepTools: Ramirez F et al. (2016) deepTools2: a next generation web server for deep-sequencing '
               'data analysis. Nucleic Acids Research 44: W160-W165.', '']
    with zipfile.ZipFile(work / 'recipe.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('README.md', '\n'.join(readme))
        z.writestr('run.sh', '\n'.join(run) + '\n')
        for r in regions:
            z.write(work / r, r)
        if job['kind'] == 'genes':
            z.writestr('name_features.py', NAME_FEATURES)
    (work / 'recipe.json').write_text(json.dumps({'commands': [shell(c) for c in commands], 'inputs': inputs,
                                                  'regions': regions, 'deeptools': version}))


# The table's last step, as a script of its own for the recipe (run_genes does the same).
NAME_FEATURES = '''"""Adds each feature's name (from regions.bed) to deepTools' raw.tab: signal-per-feature.csv."""
import csv

named = {}
for line in open('regions.bed'):
    c, s, e, name, _, strand = line.rstrip('\\n').split('\\t')
    named.setdefault((c, int(s), int(e)), (name, strand))
rows, header = [], None
for line in open('raw.tab'):
    f = line.rstrip('\\n').split('\\t')
    if line.startswith('#'):
        header = [h.strip("'") for h in f[3:]]
        continue
    name, strand = named.get((f[0], int(f[1]), int(f[2])), (f[0] + ':' + str(int(f[1]) + 1), '.'))
    rows.append([name, f[0], int(f[1]) + 1, int(f[2]), strand] + f[3:])
rows.sort(key=lambda r: (r[1], r[2]))
with open('signal-per-feature.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['name', 'chrom', 'start', 'end', 'strand'] + (header or []))
    w.writerows(rows)
'''


def _add_track(conn, data_root: Path, job, params: dict, a: dict, b: dict, made: Path, operation: str) -> str:
    """The comparison, kept as a coverage track of the member's lab, as an upload is."""
    track_id = uuid.uuid4().hex
    folder = Path(data_root) / 'tracks' / track_id
    folder.mkdir(parents=True, exist_ok=True)
    name = (params.get('name') or '').strip()[:200] or \
        operation.replace('A', a['name'], 1).replace('B', b['name'], 1)
    stored = folder / (re.sub(r'[^A-Za-z0-9_.-]+', '_', name)[:80] + '.bw')
    shutil.move(str(made), str(stored))
    conn.execute('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, file_path, '
                 "created_by) VALUES (?, ?, ?, 'bigwig', ?, 'lab_only', 'ready', ?, ?)",
                 (track_id, job['genome_id'], name, job['owner_lab'], str(stored), job['member']))
    return track_id


# --- the endpoints -----------------------------------------------------------------------

def register_analysis(app) -> None:
    from flask import abort, jsonify, request, send_file

    import access
    import db
    import identity
    import labs_reader

    def _db():
        return db.connect(Path(app.config['DATA_ROOT']) / 'genome.db')

    def _viewable(conn, member, genome_id, ids):
        rows = {}
        for i in ids:
            row = conn.execute('SELECT id, genome_id, name, type, owner_lab, visibility, track_status, file_path '
                               'FROM tracks WHERE id = ?', (str(i),)).fetchone()
            if row is None or row['genome_id'] != genome_id or not access.can_view_track(member, row, conn):
                raise AnalysisError('A track chosen is not on this genome, or not yours to see.')
            rows[row['id']] = row
        return rows

    @app.route('/api/genomes/<genome_id>/analyses', methods=['POST'])
    def start_analysis(genome_id):
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401
        body = request.get_json(silent=True) or {}
        kind = body.get('kind')
        if kind not in KINDS:
            return jsonify({'error': 'Choose what to draw.'}), 400
        conn = _db()
        try:
            genome = conn.execute('SELECT id, owner_lab, genome_visibility FROM reference_genomes WHERE id = ?',
                                  (genome_id,)).fetchone()
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            tracks = [str(t) for t in (body.get('tracks') or []) if t][:MAX_TRACKS]
            sets = [r for r in (body.get('regions') or []) if isinstance(r, dict)][:MAX_REGION_SETS]
            params = {k: body.get(k) for k in ('mode', 'before', 'after', 'body', 'bin', 'heatmap', 'kmeans',
                                                'colormap', 'title', 'method', 'operation', 'pseudocount', 'name')}
            params['tracks'] = tracks
            params['regions'] = [{'track': str(r.get('track') or ''), 'type': str(r.get('type') or '')[:80] or None,
                                  'names': [str(n)[:80] for n in (r.get('names') or [])][:5000]} for r in sets]
            params['a'], params['b'] = str(body.get('a') or ''), str(body.get('b') or '')
            ids = tracks + [r['track'] for r in params['regions']] + [x for x in (params['a'], params['b']) if x]
            try:
                rows = _viewable(conn, member, genome_id, ids)
            except AnalysisError as exc:
                return jsonify({'error': str(exc)}), 400
            owner_lab = None
            if kind == 'compare':
                if not params['a'] or not params['b'] or params['a'] == params['b']:
                    return jsonify({'error': 'Choose two different coverage tracks.'}), 400
                owner_lab = labs_reader.current_lab(member)
                if not owner_lab:
                    return jsonify({'error': 'You must belong to a lab to keep a new track.'}), 400
            elif not tracks:
                return jsonify({'error': 'Choose at least one coverage track.'}), 400
            if any(rows[t]['type'] != 'bigwig' for t in tracks):
                return jsonify({'error': 'Only coverage (BigWig) tracks can be measured.'}), 400
            if kind in ('profile', 'genes') and not params['regions']:
                return jsonify({'error': 'Choose the features: a gene annotation or BED track, and a type.'}), 400
            busy = conn.execute("SELECT count(*) FROM analysis_jobs WHERE member = ? AND state IN ('queued', 'running')",
                                (member,)).fetchone()[0] if conn.execute(
                "SELECT 1 FROM sqlite_master WHERE name = 'analysis_jobs'").fetchone() else 0
            if busy >= 3:
                return jsonify({'error': 'You have three analyses waiting already: let them finish first.'}), 429
            job_id = enqueue(conn, genome_id, member, owner_lab, kind, params)
            ahead = conn.execute("SELECT count(*) FROM analysis_jobs WHERE state IN ('queued', 'running') AND id <> ?",
                                 (job_id,)).fetchone()[0]
        finally:
            conn.close()
        return jsonify({'id': job_id, 'state': 'queued', 'ahead': ahead}), 202

    def _own_job(job_id):
        member = identity.require_member()
        conn = _db()
        try:
            ensure_schema(conn)
            job = conn.execute('SELECT * FROM analysis_jobs WHERE id = ?', (job_id,)).fetchone()
            ahead = conn.execute("SELECT count(*) FROM analysis_jobs WHERE state IN ('queued', 'running') "
                                 'AND created_at < ?', (job['created_at'],)).fetchone()[0] if job else 0
        finally:
            conn.close()
        if job is None or job['member'] != member:
            abort(404)
        return job, ahead

    @app.route('/api/analyses/<job_id>')
    def analysis_status(job_id):
        job, ahead = _own_job(job_id)
        out = {'id': job['id'], 'kind': job['kind'], 'state': job['state'], 'error': job['error'],
               'created_at': job['created_at'], 'started_at': job['started_at'], 'finished_at': job['finished_at'],
               'ahead': ahead if job['state'] == 'queued' else 0}
        if job['result']:
            out['result'] = json.loads(job['result'])
        return jsonify(out)

    @app.route('/api/analyses/<job_id>/files/<name>')
    def analysis_file(job_id, name):
        job, _ = _own_job(job_id)
        result = json.loads(job['result'] or '{}')
        if name not in (result.get('files') or []):
            abort(404)
        path = job_dir(Path(app.config['DATA_ROOT']), job_id) / name
        if not path.exists():
            abort(404)
        types = {'.svg': 'image/svg+xml', '.png': 'image/png', '.csv': 'text/csv', '.tab': 'text/plain',
                 '.gz': 'application/gzip'}
        res = send_file(str(path), mimetype=types.get(path.suffix, 'application/octet-stream'),
                        as_attachment=path.suffix not in ('.svg', '.png'), download_name=name, max_age=0)
        res.headers['X-Content-Type-Options'] = 'nosniff'
        res.headers['Content-Security-Policy'] = "default-src 'none'; style-src 'unsafe-inline'; sandbox"
        return res

    @app.route('/api/analyses/<job_id>/recipe')
    def analysis_recipe(job_id):
        """The commands that made it, and its inputs: for the code page in the notebook."""
        job, _ = _own_job(job_id)
        path = job_dir(Path(app.config['DATA_ROOT']), job_id) / 'recipe.json'
        if job['state'] != 'done' or not path.exists():
            abort(404)
        return jsonify(json.loads(path.read_text()))

    @app.route('/api/tracks/<track_id>/feature-types')
    def track_feature_types(track_id):
        member = identity.require_member()
        conn = _db()
        try:
            row = conn.execute('SELECT id, genome_id, name, type, owner_lab, visibility, track_status, file_path '
                               'FROM tracks WHERE id = ?', (track_id,)).fetchone()
            if row is None or row['type'] not in ('gff', 'bed') or not access.can_view_track(member, row, conn):
                abort(404)
        finally:
            conn.close()
        return jsonify({'types': [{'type': t, 'count': n} for t, n in feature_types(row['file_path'], row['type'])]})

    @app.route('/api/genomes/<genome_id>/region-means', methods=['POST'])
    def region_means(genome_id):
        """Each coverage track's mean (and maximum) over one region: for a bar chart of the selection."""
        import pyBigWig
        import chrom_names
        member = identity.require_member()
        body = request.get_json(silent=True) or {}
        try:
            chrom, start, end = str(body['chrom']), int(body['start']), int(body['end'])
        except (KeyError, TypeError, ValueError):
            return jsonify({'error': 'Give the region: chrom, start and end.'}), 400
        if end <= start:
            return jsonify({'error': 'The region is empty.'}), 400
        conn = _db()
        try:
            genome = conn.execute('SELECT id, owner_lab, genome_visibility FROM reference_genomes WHERE id = ?',
                                  (genome_id,)).fetchone()
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            try:
                rows = _viewable(conn, member, genome_id, [str(t) for t in (body.get('tracks') or [])][:40])
            except AnalysisError as exc:
                return jsonify({'error': str(exc)}), 400
        finally:
            conn.close()
        out = []
        for tid in [str(t) for t in (body.get('tracks') or [])][:40]:
            row = rows.get(tid)
            if row is None or row['type'] != 'bigwig' or not row['file_path']:
                continue
            bw = pyBigWig.open(row['file_path'])
            try:
                theirs = bw.chroms() or {}
                name = chrom if chrom in theirs else chrom_names.match_names([chrom], list(theirs)).renames.get(chrom, chrom)
                if name not in theirs:
                    mean = peak = None
                else:
                    e = min(end, theirs[name])
                    mean = (bw.stats(name, start, e, type='mean') or [None])[0]
                    peak = (bw.stats(name, start, e, type='max') or [None])[0]
            finally:
                bw.close()
            out.append({'id': tid, 'name': row['name'], 'mean': mean, 'max': peak})
        return jsonify({'chrom': chrom, 'start': start, 'end': end, 'tracks': out})
