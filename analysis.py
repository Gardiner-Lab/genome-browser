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
import sys
import time
import uuid
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import feature_index
import regions
import signal_kind

KINDS = ('profile', 'compare', 'correlation', 'genes', 'motif', 'fingerprint', 'methexpr')
PROFILE_MODES = ('scale', 'start', 'end', 'center', 'summit')
# Where a methylation-against-expression plot measures methylation: the gene, or the stretch before it.
METHEXPR_REGIONS = {'body': 'over each gene', 'up1000': 'over the 1 kb before each gene',
                    'up2000': 'over the 2 kb before each gene'}
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

# Where the job running now has got to: run_job points this at a function that records a step and a
# fraction done, which the dashboard's My jobs shows as a progress bar.
_PROGRESS: List[Optional[Callable[[str, float], None]]] = [None]


def _progress(step: str, fraction: float) -> None:
    if _PROGRESS[0] is not None:
        _PROGRESS[0](step, max(0.0, min(1.0, fraction)))


def _progress_recorder(conn, job_id: str) -> Callable[[str, float], None]:
    """Records progress in the job's row, at most about once a second unless the step changes."""
    last = {'step': None, 'at': 0.0}

    def record(step: str, fraction: float) -> None:
        now = time.time()
        if step == last['step'] and now - last['at'] < 1.0:
            return
        last.update(step=step, at=now)
        try:
            conn.execute('UPDATE analysis_jobs SET step = ?, progress = ? WHERE id = ?',
                         (step, round(fraction, 4), job_id))
            conn.commit()
        except sqlite3.Error:
            pass        # a progress bar is never worth failing the job over
    return record


def ensure_schema(conn) -> None:
    conn.executescript(SCHEMA)
    # Added in 0.13.0: where a running job has got to, for the progress bars on the suite's dashboard.
    have = {r[1] for r in conn.execute('PRAGMA table_info(analysis_jobs)')}
    for column, kind in (('step', 'TEXT'), ('progress', 'REAL')):
        if column not in have:
            conn.execute(f'ALTER TABLE analysis_jobs ADD COLUMN {column} {kind}')


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
            conn.execute("UPDATE analysis_jobs SET state = 'running', started_at = datetime('now'), "
                         "step = NULL, progress = 0 WHERE id = ?", (row['id'],))
        conn.execute('COMMIT')
    except Exception:
        conn.execute('ROLLBACK')
        raise
    return row


def recover(conn) -> int:
    """A job left running by a restart is run again."""
    ensure_schema(conn)
    n = conn.execute("UPDATE analysis_jobs SET state = 'queued', started_at = NULL, step = NULL, progress = NULL "
                     "WHERE state = 'running'").rowcount
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


def summit_bed(path: str, names: Optional[List[str]], out: Path) -> Tuple[int, int]:
    """
    A BED track's peaks as their summits, one base each: a narrowPeak file's 10th column is the summit's
    distance from the peak's start (-1 when not called). A peak without one is put at its middle.
    Returns how many peaks, and how many of them had a summit.
    """
    import gzip
    wanted = {n.strip().casefold() for n in (names or []) if n.strip()}
    opener = gzip.open if open(path, 'rb').read(2) == b'\x1f\x8b' else open
    n = with_summit = 0
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh, open(out, 'w', encoding='utf-8') as w:
        for line in fh:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            f = line.rstrip('\r\n').split('\t') if '\t' in line else line.split()
            try:
                chrom, start, end = f[0], int(f[1]), int(f[2])
            except (IndexError, ValueError):
                continue
            if end <= start:
                continue
            name = f[3] if len(f) > 3 and f[3] not in ('', '.') else f'{chrom}:{start + 1}'
            if wanted and name.casefold() not in wanted:
                continue
            at = (start + end) // 2
            if len(f) >= 10:
                try:
                    offset = int(f[9])
                except ValueError:
                    offset = -1
                if 0 <= offset < end - start:
                    at = start + offset
                    with_summit += 1
            label = re.sub(r'\s+', '_', name)
            strand = f[5] if len(f) > 5 and f[5] in ('+', '-') else '+'
            w.write(f'{chrom}\t{at}\t{at + 1}\t{label}\t0\t{strand}\n')
            n += 1
    return n, with_summit


def _label(name: str) -> str:
    """A track name as deepTools shows it in a legend: no characters that break its arguments."""
    return re.sub(r'\s+', ' ', name).strip()[:60] or 'track'


def _run(args: List[str], cwd: Path, what: str, at: Optional[float] = None,
         until: Optional[float] = None, expect: float = 0.0, recipe: Optional[List[str]] = None) -> None:
    """
    Run one command of the analysis.

    For the progress bar: ``at`` is how far through the job it starts. deepTools says nothing of how
    far it has got, so for a long step, ``expect`` (its likely seconds) moves the bar from ``at``
    towards ``until`` as time passes, never past 95% of the way: it only reaches ``until`` when the
    step is really done and the next one begins.
    """
    if at is not None:
        _progress(what, at)
    _RECIPE.append(list(recipe or args))
    env = dict(os.environ, MPLBACKEND='Agg', MPLCONFIGDIR=str(cwd / '.mpl'))
    started = time.time()
    proc = subprocess.Popen(args, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    while True:
        try:
            out, err = proc.communicate(timeout=1.0)
            break
        except subprocess.TimeoutExpired:
            elapsed = time.time() - started
            if elapsed > TIMEOUT_SECONDS:
                proc.kill()
                proc.communicate()
                raise AnalysisError(f'{what} took too long: try fewer tracks or features.')
            if at is not None and until is not None and expect > 0:
                _progress(what, at + (until - at) * min(0.95, elapsed / expect))
    if proc.returncode != 0:
        detail = (err or out or '').strip().splitlines()[-2:]
        raise AnalysisError(f'{what} failed' + (': ' + ' '.join(detail) if detail else '.'))


# --- the analyses ----------------------------------------------------------------------

def _is_level(t: dict) -> bool:
    """Whether a coverage track is measured as a methylation level (signal_kind.py)."""
    return t.get('signal') == signal_kind.METHYLATION


def _source(t: dict) -> str:
    """
    The file to measure a track from. A methylation track that gives minus-strand cytosines a minus
    sign is measured from a copy with the sign dropped, made once beside it; the recipe makes it too.
    """
    if not _is_level(t) or not signal_kind.is_signed(t['file_path']):
        return t['file_path']
    path = signal_kind.level_file(t['file_path'])
    step = ['python3', 'bigwig_means.py', 'unsigned', t['file_path'], '-o', path]
    if step not in _RECIPE:
        _RECIPE.append(step)
    return path


def _measured_as(tracks: List[dict], out: dict) -> dict:
    """Says in a result which tracks were measured as a methylation level, and which had a sign dropped."""
    level = [t['name'] for t in tracks if _is_level(t)]
    if level:
        out['methylation'] = level
        unsigned = [t['name'] for t in tracks if _is_level(t) and signal_kind.is_signed(t['file_path'])]
        if unsigned:
            out['unsigned'] = unsigned
    return out


def _as(tracks: List[dict]) -> List[str]:
    """bigwig_means.py's --as, when any track is measured as a level (the default is coverage for all)."""
    if not any(_is_level(t) for t in tracks):
        return []
    return ['--as', *['level' if _is_level(t) else 'coverage' for t in tracks]]


def run_profile(params: dict, tracks: List[dict], region_tracks: Dict[str, dict], work: Path) -> dict:
    sets, names, found = [], [], []
    mode = params.get('mode') if params.get('mode') in PROFILE_MODES else 'scale'
    summits = [0, 0]
    for i, rs in enumerate(params.get('regions') or []):
        rt = region_tracks.get(rs.get('track'))
        if rt is None:
            continue
        bed = work / f'regions-{i + 1}.bed'
        if mode == 'summit' and rt['type'] == 'bed':
            # Peaks: each one's summit (a narrowPeak file's 10th column), or its middle without one.
            n, with_summit = summit_bed(rt['file_path'], rs.get('names'), bed)
            summits[0] += n
            summits[1] += with_summit
        else:
            n = region_bed(rt['file_path'], rt['type'], rs.get('type') or None, rs.get('names'), bed)
        if n:
            sets.append(bed)
            names.append(rs.get('type') or rt['name'])
            found.append(n)
    if not sets:
        raise AnalysisError('None of the chosen features were found: choose another type, or check the names.')
    before, after = (int(_clamp(params.get(k), 0, 20000, 2000)) for k in ('before', 'after'))
    body = int(_clamp(params.get('body'), 500, 20000, 3000))
    binsize = int(_clamp(params.get('bin'), 10, 5000, 100))
    # The member's two choices, each on unless turned off (as every profile was before they could choose).
    # --skipZeros leaves out features with no signal anywhere in the window. --missingDataAsZero counts
    # stretches a track has no data for as zero; without it they are NaN, which plotHeatmap paints
    # black and plotProfile leaves out of the mean.
    skip_zeros = params.get('skip_zeros') is not False
    missing_as_zero = params.get('missing_as_zero') is not False
    # A methylation track is measured as its level (signal_kind.py): without its sign, and with the
    # stretches between cytosines left out whatever the choice above, which is then for the coverage
    # tracks alone. deepTools' flag is for every track at once, so with both kinds the coverage
    # tracks' missing values are set to zero after (bigwig_means.py fill).
    levels = [_is_level(t) for t in tracks]
    paths = [_source(t) for t in tracks]
    fill_after = missing_as_zero and any(levels) and not all(levels)
    args = ['computeMatrix', 'scale-regions' if mode == 'scale' else 'reference-point',
            '-S', *paths, '-R', *[str(s) for s in sets],
            '--samplesLabel', *[_label(t['name']) for t in tracks],
            '-b', str(before), '-a', str(after), '-bs', str(binsize), '-p', str(THREADS),
            *(['--skipZeros'] if skip_zeros else []),
            *(['--missingDataAsZero'] if missing_as_zero and not any(levels) else []),
            '-o', 'matrix.gz', '--quiet']
    if mode == 'scale':
        args += ['-m', str(body)]
    else:
        # A summit is a one-base feature (summit_bed), or a gene's middle where it has none: its centre.
        args += ['--referencePoint', {'start': 'TSS', 'end': 'TES'}.get(mode, 'center')]
    # How long computeMatrix is likely to take, for the progress bar only: about 0.05 ms for each
    # feature, track and step (30,000 genes, one track, 70 steps took 95-118 s in the worker).
    steps = (before + after + (body if mode == 'scale' else 0)) / binsize
    expect = sum(found) * len(tracks) * max(steps, 1) * 5e-5
    _run(args, work, 'Computing the profile', at=0.0, until=0.85 if params.get('heatmap') else 0.9, expect=expect)
    if fill_after:
        coverage = [str(i + 1) for i, lv in enumerate(levels) if not lv]
        _measure(['fill', str(work / 'matrix.gz'), '--samples', *coverage, '-o', str(work / 'matrix.gz')],
                 'Counting stretches without reads as zero', span=(0.85, 0.85))
    # Label each set with the features actually drawn: computeMatrix leaves out those with no signal
    # at all (--skipZeros) and those shorter than a bin, so the count found would overstate it.
    drawn = _features_drawn(work / 'matrix.gz')
    if len(drawn) != len(sets):
        raise AnalysisError('None of the features in one of the sets has any signal here: '
                            'choose another type, or check the names.')
    labels, left_out = [], {}
    for name, n, kept in zip(names, found, drawn):
        labels.append(f'{name} ({kept:,})' if kept == n else f'{name} ({kept:,} of {n:,})')
        if kept < n:
            left_out[name] = n - kept
    # Region-set labels as chosen (computeMatrix names them after the files).
    title = (params.get('title') or '').strip()[:100]
    start_label, end_label = ('start', 'end')
    common = ['--regionsLabel', *labels, '--startLabel', start_label, '--endLabel', end_label,
              '--refPointLabel', {'start': 'start', 'end': 'end', 'center': 'middle'}.get(mode, 'summit')]
    y_label = 'methylation level' if all(levels) else 'mean signal'
    profile_args = ['plotProfile', '-m', 'matrix.gz', '-o', 'profile.svg', '--perGroup', *common,
                    '--plotTitle', title or ('Methylation' if all(levels) else 'Mean signal'), '--yAxisLabel', y_label]
    _run(profile_args, work, 'Drawing the profile', at=0.85 if params.get('heatmap') else 0.9)
    files = ['profile.svg', 'matrix.gz']
    if params.get('heatmap'):
        k = int(_clamp(params.get('kmeans'), 0, 8, 0))
        colormap = str(params.get('colormap') or 'viridis')
        if colormap not in COLORMAPS:
            colormap = 'viridis'
        heat_matrix = 'matrix.gz'
        if any(levels):
            # A methylation track has steps with no cytosine in nearly every row, and plotHeatmap paints a
            # pixel holding any of them as missing: the whole heatmap came out grey. For the heatmap only,
            # each empty step takes its value from the measured steps beside it in the same row. The
            # profiles over its columns are still drawn from the matrix as measured (--profileMatrix).
            heat_matrix = 'heatmap-matrix.gz'
            _measure(['interpolate', str(work / 'matrix.gz'), '--samples',
                      *[str(i + 1) for i, lv in enumerate(levels) if lv], '-o', str(work / heat_matrix)],
                     'Drawing the heatmap', span=(0.9, 0.9))
        # With k-means the rows are regrouped into clusters, which deepTools names itself: our names for
        # the feature sets no longer apply, and plotHeatmap refuses them ("length new labels").
        heat_common = common[1 + len(labels):] if k >= 2 else common
        heat_args = ['plotHeatmap', '-m', heat_matrix, '-o', 'heatmap.png', '--dpi', '110', *heat_common,
                     '--colorMap', colormap, '--heatmapHeight', '18', '--plotTitle', title or 'Signal per feature']
        if k >= 2:
            # Each feature's cluster, as plotHeatmap drew them: for gene lists (_clusters).
            heat_args += ['--kmeans', str(k), '--outFileSortedRegions', 'heatmap-regions.bed']
        if any(levels):
            # A row with no cytosine measured anywhere: light grey, not deepTools' black, which reads as
            # a value at the bottom of the scale.
            heat_args += ['--missingDataColor', '0.85', '--profileMatrix', 'matrix.gz']
        if not missing_as_zero and levels.count(False):
            # Coverage left with gaps (as chosen): each pixel one feature's value, not a blend that any
            # gap turns to the missing colour.
            heat_args += ['--interpolationMethod', 'nearest']
        # Through bigwig_means.py, which turns the tracks' names over their columns 45 degrees, so long
        # names don't run into each other.
        heat_args = ['heatmap', *heat_args[1:]]
        _run([sys.executable, str(Path(__file__).with_name('bigwig_means.py')), *heat_args], work,
             'Drawing the heatmap', at=0.9, recipe=['python3', 'bigwig_means.py', *heat_args])
        files.insert(1, 'heatmap.png')
        if k >= 2:
            clusters = _clusters(work)
            files += ['clusters.csv'] + [c['file'] for c in clusters]
    out = _measured_as(tracks, {'files': files, 'regions': labels, 'skip_zeros': skip_zeros,
                                'missing_as_zero': missing_as_zero, 'coverage_tracks': levels.count(False)})
    out['mode'] = mode
    if mode == 'summit':
        out['summits'] = {'peaks': summits[0], 'with_summit': summits[1]}
    if params.get('heatmap') and int(_clamp(params.get('kmeans'), 0, 8, 0)) >= 2:
        out['clusters'] = [{k: c[k] for k in ('cluster', 'count', 'file', 'first')} for c in clusters]
    if params.get('heatmap') and any(levels):
        out['heatmap_filled'] = [t['name'] for t in tracks if _is_level(t)]
    if left_out:
        out['left_out'] = left_out
    # Features whose window runs past a chromosome end have bins with no genome under them: those count
    # as zero (or are left out, without --missingDataAsZero) like any stretch without reads. Usually a
    # handful; they are listed so the member can judge whether it matters for theirs.
    off = _off_chromosome_ends(work / 'matrix.gz', tracks, mode, before, after)
    if off:
        with open(work / 'off-chromosome-ends.csv', 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['name', 'chrom', 'start', 'end', 'strand', 'bases_past_start', 'bases_past_end'])
            w.writerows(off)
        files.append('off-chromosome-ends.csv')
        out['off_end'] = {'count': len(off), 'features': [r[0] for r in off[:20]]}
    return out


def _off_chromosome_ends(matrix: Path, tracks: List[dict], mode: str, before: int, after: int) -> List[list]:
    """The drawn features whose window (the flanks, by strand) runs past either end of its chromosome."""
    import gzip
    import pyBigWig
    sizes: Dict[str, int] = {}
    for t in tracks:
        bw = pyBigWig.open(t['file_path'])
        try:
            for c, n in (bw.chroms() or {}).items():
                sizes[c] = min(n, sizes.get(c, n))
        finally:
            bw.close()
    out = []
    with gzip.open(matrix, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('@'):
                continue
            f = line.split('	', 6)
            chrom, s, e, name, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]
            if chrom not in sizes:
                continue
            minus = strand == '-'
            if mode == 'scale':
                left, right = (s - after, e + before) if minus else (s - before, e + after)
            elif mode in ('center', 'summit'):
                at = (s + e) // 2
                left, right = (at - after, at + before) if minus else (at - before, at + after)
            else:
                # The reference point: a feature's start (TSS) or end (TES), by its strand.
                at = (e if minus else s) if mode == 'start' else (s if minus else e)
                left, right = (at - after, at + before) if minus else (at - before, at + after)
            past_start, past_end = max(0, -left), max(0, right - sizes[chrom])
            if past_start or past_end:
                out.append([name, chrom, s + 1, e, strand, past_start, past_end])
    return out


def _clusters(work: Path) -> List[dict]:
    """
    The heatmap's k-means clusters as lists: clusters.csv (every feature, its cluster), and for each
    cluster cluster-N.txt, its features' names one to a line, to paste into a GO-enrichment tool.
    From plotHeatmap's --outFileSortedRegions, whose 13th column is the cluster ('cluster_1', ...).
    """
    found: Dict[str, List[list]] = {}
    with open(work / 'heatmap-regions.bed', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            group = f[12] if len(f) > 12 else 'cluster_1'
            found.setdefault(group, []).append([f[3], f[0], int(f[1]) + 1, int(f[2]), f[5] if len(f) > 5 else '.'])
    order = sorted(found, key=lambda g: int(re.sub(r'\D', '', g) or 0))
    out = []
    with open(work / 'clusters.csv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'chrom', 'start', 'end', 'strand', 'cluster'])
        for n, group in enumerate(order, 1):
            rows = found[group]
            w.writerows(r + [n] for r in rows)
            name = f'cluster-{n}.txt'
            # A feature in several sets (a gene in two lists) once only.
            names = list(dict.fromkeys(r[0] for r in rows))
            (work / name).write_text('\n'.join(names) + '\n', encoding='utf-8')
            out.append({'cluster': n, 'count': len(rows), 'file': name, 'first': names[:5], 'rows': rows})
    return out


def cluster_bed(work: Path, cluster: int) -> str:
    """One cluster's features as a BED file's text (for a track of them), from clusters.csv."""
    lines = []
    with open(work / 'clusters.csv', newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if int(row['cluster']) == cluster:
                lines.append(f"{row['chrom']}\t{int(row['start']) - 1}\t{row['end']}\t{row['name']}\t0\t"
                             f"{row['strand'] if row['strand'] in ('+', '-') else '.'}\n")
    return ''.join(lines)


def _plot(args: List[str], work: Path, what: str, at: float) -> None:
    """Run a plots.py figure (in its own process, as matplotlib wants), recorded for the recipe."""
    _run([sys.executable, str(Path(__file__).with_name('plots.py')), *args], work, what, at=at,
         recipe=['python3', 'plots.py', *args])


def run_fingerprint(params: dict, tracks: List[dict], work: Path) -> dict:
    # A methylation level is not a count of reads: how evenly it is spread says nothing of enrichment.
    skipped = [t['name'] for t in tracks if _is_level(t)]
    tracks = [t for t in tracks if not _is_level(t)]
    if not tracks:
        raise AnalysisError('A fingerprint is of coverage tracks (reads): the tracks chosen are all methylation tracks.')
    binsize = int(_clamp(params.get('bin'), 100, 100_000, 1000))
    labels = [_label(t['name']) for t in tracks]
    _measure(['bins', *[t['file_path'] for t in tracks], '--labels', *labels, '-bs', str(binsize),
              '-o', str(work / 'fingerprint.npz')], 'Summing each track in bins', span=(0.0, 0.85))
    title = (params.get('title') or '').strip()[:100]
    _plot(['fingerprint', 'fingerprint.npz', '-o', 'fingerprint.svg', '--metrics', 'fingerprint-metrics.tab',
           '--summary', 'fingerprint.json', '--title', title or f'Fingerprint, {binsize:,} bp bins'],
          work, 'Drawing the fingerprint', at=0.85)
    found = json.loads((work / 'fingerprint.json').read_text())
    out = {'files': ['fingerprint.svg', 'fingerprint-metrics.tab'], 'bin': binsize, 'metrics': found}
    if skipped:
        out['skipped_methylation'] = skipped
    return out


def upstream_bed(src: Path, out: Path, length: int, sizes: Dict[str, int]) -> None:
    """
    For each feature of a BED file, the ``length`` bases before it by its strand (its promoter, for a
    gene), cut at the chromosome's ends; a row for every feature, in order, even when nothing is left.
    """
    with open(src, encoding='utf-8') as fh, open(out, 'w', encoding='utf-8', newline='\n') as w:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            s, e = int(f[1]), int(f[2])
            if f[5] == '-':
                a, b = e, min(e + length, sizes.get(f[0], e + length))
            else:
                a, b = max(0, s - length), s
            w.write('\t'.join([f[0], str(a), str(max(a, b))] + f[3:]) + '\n')


def run_methexpr(params: dict, expression: dict, tracks: List[dict], region_tracks: Dict[str, dict],
                 work: Path) -> dict:
    not_level = [t['name'] for t in tracks if not _is_level(t)]
    if not_level:
        raise AnalysisError('Not methylation tracks: ' + ', '.join(not_level) + '. Tick methylation tracks only '
                            '(a track’s Details says which it is), and choose the RNA-seq track as the expression.')
    if _is_level(expression):
        raise AnalysisError('The expression track is a methylation track: choose an RNA-seq coverage track.')
    rs = (params.get('regions') or [{}])[0]
    rt = region_tracks.get(rs.get('track'))
    if rt is None:
        raise AnalysisError('Choose the features: the genes.')
    bed = work / 'regions.bed'
    n = region_bed(rt['file_path'], rt['type'], rs.get('type') or None, rs.get('names'), bed)
    if n < 10:
        raise AnalysisError('Fewer than ten of the chosen features were found: choose the genes.')
    where = params.get('region') if params.get('region') in METHEXPR_REGIONS else 'body'
    groups = int(_clamp(params.get('groups'), 2, 10, 5))
    meth_bed = bed
    if where != 'body':
        import pyBigWig
        bw = pyBigWig.open(expression['file_path'])
        try:
            sizes = dict(bw.chroms() or {})
        finally:
            bw.close()
        meth_bed = work / 'regions-upstream.bed'
        upstream_bed(bed, meth_bed, int(where[2:]), sizes)
    _measure(['features', str(bed), expression['file_path'], '--labels', _label(expression['name']),
              '-o', str(work / 'expression.tab')], 'Measuring expression', span=(0.0, 0.25))
    levels = [_source(t) for t in tracks]
    _measure(['features', str(meth_bed), *levels, '--labels', *[_label(t['name']) for t in tracks],
              '--as', *['level'] * len(tracks), '-o', str(work / 'methylation.tab')],
             'Measuring methylation', span=(0.25, 0.55))
    title = (params.get('title') or '').strip()[:100]
    _plot(['methexpr', '--genes', 'regions.bed', '--expression', 'expression.tab', '--methylation', 'methylation.tab',
           '-o', 'methexpr.svg', '--table', 'methylation-vs-expression.csv', '--groups', str(groups),
           '--group-beds', 'group', '--where', METHEXPR_REGIONS[where], '--summary', 'methexpr.json',
           *(['--title', title] if title else [])], work, 'Drawing methylation against expression', at=0.55)
    found = json.loads((work / 'methexpr.json').read_text())
    files = ['methexpr.svg', 'methylation-vs-expression.csv']
    if params.get('group_profile') is not False and found['beds']:
        # The classic figure: the level along the genes of each expression group, one line a group.
        _run(['computeMatrix', 'scale-regions', '-S', *levels, '-R', *[b['file'] for b in found['beds']],
              '--samplesLabel', *[_label(t['name']) for t in tracks], '-b', '2000', '-a', '2000', '-m', '3000',
              '-bs', '100', '-p', str(THREADS), '--skipZeros', '-o', 'groups-matrix.gz', '--quiet'],
             work, 'Computing the profile of each group', at=0.65, until=0.95, expect=n * len(tracks) * 70 * 5e-5)
        _run(['plotProfile', '-m', 'groups-matrix.gz', '-o', 'groups-profile.svg',
              '--regionsLabel', *[b['label'] for b in found['beds']], '--startLabel', 'start', '--endLabel', 'end',
              '--yAxisLabel', 'methylation level', '--plotTitle', 'Methylation along the genes of each expression group',
              '--colors', *_group_colours(found['beds'])], work, 'Drawing the profile of each group', at=0.95)
        files.insert(1, 'groups-profile.svg')
        files.append('groups-matrix.gz')
    return _measured_as(tracks, {'files': files, 'count': n, 'region': where, 'where': METHEXPR_REGIONS[where],
                                 'expression': expression['name'], 'groups': found['groups'],
                                 'no_expression_value': found['no_expression_value'], 'tracks': found['tracks']})


def _group_colours(beds: List[dict]) -> List[str]:
    """Grey for the no-reads group, then viridis from lowest to highest expression, as the boxes."""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import cm, colors
    expressed = [b for b in beds if not b['file'].endswith('-0.bed')]
    out = []
    for b in beds:
        if b['file'].endswith('-0.bed'):
            out.append('#bbbbbb')
        else:
            i = expressed.index(b)
            out.append(colors.to_hex(cm.viridis(0.15 + 0.75 * i / max(1, len(expressed) - 1))))
    return out


def _features_drawn(matrix: Path) -> List[int]:
    """How many features of each set computeMatrix kept, from its header's group boundaries."""
    import gzip
    with gzip.open(matrix, 'rt', encoding='utf-8') as fh:
        header = fh.readline()
    bounds = json.loads(header[1:])['group_boundaries'] if header.startswith('@') else [0]
    return [b - a for a, b in zip(bounds, bounds[1:])]


def run_compare(params: dict, a: dict, b: dict, work: Path) -> dict:
    op = params.get('operation') if params.get('operation') in OPERATIONS else 'subtract'
    binsize = int(_clamp(params.get('bin'), 10, 10000, 100))
    args = ['bigwigCompare', '-b1', _source(a), '-b2', _source(b), '--operation', op,
            '-bs', str(binsize), '-p', str(THREADS), '-o', 'compare.bw']
    if _is_level(a) or _is_level(b):
        # Between cytosines a methylation track has no data, which is not zero: leave those steps empty.
        args.append('--skipNonCoveredRegions')
    if op in ('log2', 'ratio'):
        args += ['--pseudocount', str(_clamp(params.get('pseudocount'), 0.0001, 1000, 1))]
    _run(args, work, 'Comparing the tracks', at=0.0)
    return {'files': ['compare.bw'], 'operation': OPERATIONS[op]}


def _measure(args: List[str], what: str, span: Tuple[float, float] = (0.0, 1.0)) -> None:
    """Run bigwig_means.py's step in this process, and record it for the recipe as the command it is.

    ``span`` is the part of the job it takes, for the progress bar, which it moves as it goes."""
    import bigwig_means
    _RECIPE.append(['python3', 'bigwig_means.py', *args])
    _progress(what, span[0])
    try:
        # Paths under work are written as bare names in the recipe (write_recipe).
        bigwig_means.main(args, progress=lambda f: _progress(what, span[0] + (span[1] - span[0]) * f))
    except SystemExit as exc:
        raise AnalysisError(f'{what} failed' + (f': {exc.code}' if isinstance(exc.code, str) else '.'))


def run_correlation(params: dict, tracks: List[dict], work: Path) -> dict:
    binsize = int(_clamp(params.get('bin'), 1000, 1_000_000, 10000))
    method = params.get('method') if params.get('method') in ('pearson', 'spearman') else 'pearson'
    labels = [_label(t['name']) for t in tracks]
    # Not multiBigwigSummary: it averages over covered bases only and leaves a bin with no reads NaN,
    # which plotCorrelation drops, so bins where only one replicate has reads left the correlation.
    _measure(['bins', *[_source(t) for t in tracks], '--labels', *labels, '-bs', str(binsize),
              *_as(tracks), '-o', str(work / 'summary.npz'), '--outRawCounts', str(work / 'bins.tab')],
             'Summarising the tracks', span=(0.0, 0.7))
    _run(['plotCorrelation', '-in', 'summary.npz', '-c', method, '-p', 'heatmap', '--plotNumbers',
          '--skipZeros', '--removeOutliers', '--colorMap', 'RdYlBu_r', '-o', 'correlation.svg',
          '--plotTitle', f'{method.capitalize()} correlation, {binsize:,} bp bins',
          '--outFileCorMatrix', 'correlation.tab'], work, 'Drawing the correlation', at=0.7)
    _run(['plotPCA', '-in', 'summary.npz', '-o', 'pca.svg', '--plotTitle', 'Principal components'],
         work, 'Drawing the PCA', at=0.85)
    return _measured_as(tracks, {'files': ['correlation.svg', 'pca.svg', 'correlation.tab', 'bins.tab']})


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
    # Not multiBigwigSummary: it averages over covered bases only (a gene with 10 over half its
    # length read 10, not 5) and gives a feature with no reads no value at all.
    _measure(['features', str(bed), *[_source(t) for t in tracks], '--labels', *labels, *_as(tracks),
              '-o', str(work / 'raw.tab')], 'Measuring each feature', span=(0.0, 0.95))
    _progress('Writing the table', 0.95)
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
    return _measured_as(tracks, {'files': ['signal-per-feature.csv'], 'count': len(rows), 'preview': preview,
                                 'columns': [t['name'] for t in tracks]})


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
    _PROGRESS[0] = _progress_recorder(conn, job['id'])
    try:
        ids = list(params.get('tracks') or []) + [r.get('track') for r in params.get('regions') or []]
        ids = [i for i in ids + [params.get('a'), params.get('b'), params.get('expression')] if i]
        rows = {}
        if ids:
            rows = {r['id']: dict(r) for r in conn.execute(
                f'SELECT id, name, type, file_path, track_status, signal FROM tracks '
                f'WHERE id IN ({",".join("?" * len(ids))})', ids)}
            for r in rows.values():
                if r['type'] == 'bigwig' and r['track_status'] == 'ready':
                    r['signal'] = signal_kind.ensure(conn, r)
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
        elif kind == 'fingerprint':
            if not coverage:
                raise AnalysisError('Choose at least one coverage track.')
            result = run_fingerprint(params, coverage, work)
        elif kind == 'methexpr':
            expression = rows.get(params.get('expression'))
            if not expression or expression['type'] != 'bigwig' or expression['track_status'] != 'ready':
                raise AnalysisError('Choose the expression track: an RNA-seq coverage track.')
            if not coverage:
                raise AnalysisError('Tick at least one methylation track.')
            result = run_methexpr(params, expression, coverage, region_tracks, work)
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
        conn.execute("UPDATE analysis_jobs SET state = 'done', result = ?, step = NULL, progress = 1, "
                     "finished_at = datetime('now') WHERE id = ?",
                     (json.dumps(result), job['id']))
    except Exception as exc:   # every failure is the member's to see, never the worker's end
        message = str(exc) if isinstance(exc, AnalysisError) else f'The analysis failed: {exc}'
        conn.execute("UPDATE analysis_jobs SET state = 'failed', error = ?, finished_at = datetime('now') WHERE id = ?",
                     (message[:1000], job['id']))
    finally:
        _PROGRESS[0] = None
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
    names, made = {}, {}
    for row in rows.values():
        names[row['file_path']] = _export_filename(row['name'], row['type'], row['file_path'])
        if row['type'] == 'bigwig' and row['file_path']:
            # A methylation track's copy without signs is made by the recipe (bigwig_means.py unsigned).
            stem = names[row['file_path']]
            made[str(signal_kind.level_path(row['file_path']))] = \
                (stem[:-3] if stem.endswith('.bw') else stem) + '.level.bw'
    regions = sorted(p.name for p in work.glob('regions*.bed'))

    def portable(arg: str) -> str:
        if arg in names:
            return names[arg]
        if arg in made:
            return made[arg]
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
    uses_means = any(c[:2] == ['python3', 'bigwig_means.py'] for c in commands)
    uses_plots = any(c[:2] == ['python3', 'plots.py'] for c in commands)
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
               '- Python 3 with matplotlib (comes with deepTools)']
    if uses_means:
        readme += ['- `bigwig_means.py` (included): each track’s mean over a feature or bin, with the bases a track '
                   'has no data for counted as zero (for a methylation track, `--as level`: left out); and a '
                   'methylation track without its strand signs (`unsigned`); and the heatmap, drawn by deepTools’ '
                   'plotHeatmap with the tracks’ names turned 45 degrees (`heatmap`). It needs pyBigWig and numpy, which '
                   'come with deepTools.']
    if uses_plots:
        readme += ['- `plots.py` (included): the figures deepTools has no tool for with coverage (BigWig) tracks: '
                   'the fingerprint, and methylation against expression. It needs numpy, scipy and matplotlib, which '
                   'come with deepTools.']
    readme += ['',
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
        if uses_means:
            z.write(Path(__file__).with_name('bigwig_means.py'), 'bigwig_means.py')
        if uses_plots:
            z.write(Path(__file__).with_name('plots.py'), 'plots.py')
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
                                                'colormap', 'title', 'method', 'operation', 'pseudocount', 'name',
                                                'skip_zeros', 'missing_as_zero', 'region', 'groups',
                                                'group_profile')}
            params['tracks'] = tracks
            params['regions'] = [{'track': str(r.get('track') or ''), 'type': str(r.get('type') or '')[:80] or None,
                                  'names': [str(n)[:80] for n in (r.get('names') or [])][:5000]} for r in sets]
            params['a'], params['b'] = str(body.get('a') or ''), str(body.get('b') or '')
            params['expression'] = str(body.get('expression') or '') if kind == 'methexpr' else ''
            ids = tracks + [r['track'] for r in params['regions']] + \
                [x for x in (params['a'], params['b'], params['expression']) if x]
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
            if kind == 'methexpr':
                if not params['expression'] or rows[params['expression']]['type'] != 'bigwig':
                    return jsonify({'error': 'Choose the expression track: an RNA-seq coverage track.'}), 400
                if params['expression'] in tracks:
                    return jsonify({'error': 'The expression track is ticked as a methylation track too.'}), 400
            if kind in ('profile', 'genes', 'methexpr') and not params['regions']:
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
        out = {'id': job['id'], 'kind': job['kind'], 'genome_id': job['genome_id'], 'state': job['state'],
               'error': job['error'], 'step': job['step'], 'progress': job['progress'],
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
                 '.txt': 'text/plain', '.gz': 'application/gzip'}
        res = send_file(str(path), mimetype=types.get(path.suffix, 'application/octet-stream'),
                        as_attachment=path.suffix not in ('.svg', '.png'), download_name=name, max_age=0)
        res.headers['X-Content-Type-Options'] = 'nosniff'
        res.headers['Content-Security-Policy'] = "default-src 'none'; style-src 'unsafe-inline'; sandbox"
        return res

    @app.route('/api/analyses/<job_id>/clusters/<int:cluster>/track', methods=['POST'])
    def analysis_cluster_track(job_id, cluster):
        """One k-means cluster of a heatmap, kept as a BED track of the member's lab: to view, or to draw
        a profile or a table of it."""
        import tracks as tracks_module
        job, _ = _own_job(job_id)
        result = json.loads(job['result'] or '{}')
        if job['state'] != 'done' or cluster not in [c['cluster'] for c in result.get('clusters') or []]:
            abort(404)
        member = identity.require_member()
        owner_lab = labs_reader.current_lab(member)
        if not owner_lab:
            return jsonify({'error': 'You must belong to a lab to keep a new track.'}), 400
        text = cluster_bed(job_dir(Path(app.config['DATA_ROOT']), job_id), cluster)
        if not text:
            return jsonify({'error': 'That cluster has no features.'}), 400
        body = request.get_json(silent=True) or {}
        params = json.loads(job['params'] or '{}')
        name = (str(body.get('name') or '').strip()[:200] or
                f"{(params.get('title') or '').strip() or 'Heatmap'}, cluster {cluster} of {len(result['clusters'])}")
        conn = _db()
        try:
            made = tracks_module.add_bed_track(conn, Path(app.config['DATA_ROOT']), job['genome_id'], owner_lab,
                                               member, name, text)
        finally:
            conn.close()
        return jsonify(made), 201

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
        import bigwig_means
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
        conn = _db()
        try:
            kinds = {tid: signal_kind.ensure(conn, conn.execute('SELECT * FROM tracks WHERE id = ?', (tid,)).fetchone())
                     for tid, row in rows.items() if row['type'] == 'bigwig' and row['file_path']}
        finally:
            conn.close()
        for tid in [str(t) for t in (body.get('tracks') or [])][:40]:
            row = rows.get(tid)
            if row is None or row['type'] != 'bigwig' or not row['file_path']:
                continue
            # A methylation track: its level over the cytosines in the region, without strand signs.
            level = kinds.get(tid) == signal_kind.METHYLATION
            bw = pyBigWig.open(signal_kind.level_file(row['file_path']) if level else row['file_path'])
            try:
                theirs = bw.chroms() or {}
                name = chrom if chrom in theirs else chrom_names.match_names([chrom], list(theirs)).renames.get(chrom, chrom)
                if name not in theirs:
                    mean = peak = None
                else:
                    # Bases without reads count as zero (bigwig_means.py): pyBigWig's own mean skips them.
                    e = min(end, theirs[name])
                    mean = bigwig_means.mean(bw, name, start, e, level)
                    peak = (bw.stats(name, start, e, type='max', exact=True) or [None])[0]
                    if peak is None and not level:
                        peak = 0.0
            finally:
                bw.close()
            out.append({'id': tid, 'name': row['name'], 'mean': mean, 'max': peak,
                        'signal': kinds.get(tid, signal_kind.COVERAGE)})
        return jsonify({'chrom': chrom, 'start': start, 'end': end, 'tracks': out})
