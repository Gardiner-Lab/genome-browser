"""
My jobs: the member's analyses and uploads being prepared, with how far each has got.

The Genome Browser runs its slow work in the genome worker, one job at a time: analyses (profiles,
replicate checks, tables, comparisons, marking sequences) in ``analysis_jobs``, and uploads being
turned into coverage tracks (BAM, WIG) in ``conversion_jobs``. This module shows a member theirs:

  GET /api/jobs              their jobs from the last KEEP_DAYS, newest first, as JSON
  GET /jobs                  a page of them with progress bars, which refreshes itself while any run;
                             ``?embed=1`` leaves out the header, for the suite's dashboard
  GET /suite-manifest.json   the dashboard's "My jobs" tab (README of Trellis Tools, "The dashboard")

How far a job has got: an analysis records its step and a fraction done as it goes (analysis.py,
``_progress``); marking sequences has an estimate of how long it takes, so its bar runs against
that; preparing an upload has no measure, so its bar only says it is working.
"""
import json
from pathlib import Path
from typing import Dict, List, Optional

from flask import Response, current_app, jsonify, request

import analysis
import db
import identity

KINDS = {
    'profile': 'Profile around features, and heatmap',
    'correlation': 'Replicate check',
    'genes': 'Signal per feature (table)',
    'compare': 'Comparing two tracks',
    'motif': 'Marking sequences on the genome',
    'bam': 'Preparing an upload (BAM to coverage)',
    'wig': 'Preparing an upload (WIG to BigWig)',
}
LIMIT = 50
# A marked sequence's bar never claims to be done before it is.
MOTIF_CAP = 0.95


def register_jobs(app) -> None:

    def _db():
        return db.connect(Path(app.config['DATA_ROOT']) / 'genome.db')

    @app.route('/api/jobs')
    def my_jobs():
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401
        conn = _db()
        try:
            analysis.ensure_schema(conn)
            jobs = _analysis_jobs(conn, member) + _conversion_jobs(conn, member)
        finally:
            conn.close()
        mount = _mount()
        for job in jobs:
            job['open_url'] = _open_url(mount, job)
        # Running first, then waiting, then the rest; newest first within each.
        order = {'running': 0, 'queued': 1}
        jobs.sort(key=lambda j: j['created_at'] or '', reverse=True)
        jobs.sort(key=lambda j: order.get(j['state'], 2))
        return jsonify({'jobs': jobs[:LIMIT], 'keep_days': analysis.KEEP_DAYS,
                        'busy': any(j['state'] in ('queued', 'running') for j in jobs)})

    @app.route('/jobs')
    def jobs_page():
        static = _static()
        return Response(_page_html(static, _mount(), request.args.get('embed') == '1'), mimetype='text/html')

    @app.route('/suite-manifest.json')
    def suite_manifest():
        mount = _mount()
        return jsonify({'tabs': [{'id': 'genome-jobs', 'title': 'My jobs', 'icon': '⏳', 'section': 'me',
                                  'src': f'{mount}/jobs?embed=1'}]})


# --- the jobs -----------------------------------------------------------------------------

def _analysis_jobs(conn, member: str) -> List[dict]:
    rows = conn.execute(
        'SELECT j.*, g.name AS genome_name, '
        "  (julianday('now') - julianday(j.started_at)) * 86400 AS running_for, "
        "  (julianday(coalesce(j.finished_at, 'now')) - julianday(coalesce(j.started_at, j.created_at))) * 86400 "
        '    AS took '
        'FROM analysis_jobs j LEFT JOIN reference_genomes g ON g.id = j.genome_id '
        "WHERE j.member = ? AND j.created_at >= datetime('now', ?) "
        'ORDER BY j.created_at DESC LIMIT ?',
        (member, f'-{analysis.KEEP_DAYS} days', LIMIT)).fetchall()
    names = _track_names(conn, rows)
    out = []
    for r in rows:
        params = _json(r['params'])
        result = _json(r['result'])
        job = {'id': r['id'], 'source': 'analysis', 'kind': r['kind'],
               'what': KINDS.get(r['kind'], r['kind']), 'about': _about(r['kind'], params, result, names),
               'genome': {'id': r['genome_id'], 'name': r['genome_name'] or 'a genome no longer here'},
               'state': r['state'], 'step': None, 'progress': None, 'ahead': 0, 'error': r['error'],
               'created_at': r['created_at'], 'started_at': r['started_at'], 'finished_at': r['finished_at'],
               'seconds': round(r['took']) if r['took'] is not None else None}
        if r['state'] == 'queued':
            job['progress'] = 0.0
            job['ahead'] = conn.execute("SELECT count(*) FROM analysis_jobs WHERE state IN ('queued', 'running') "
                                        'AND created_at < ?', (r['created_at'],)).fetchone()[0]
        elif r['state'] == 'running':
            job['step'] = r['step']
            if r['kind'] == 'motif':
                estimate = float(params.get('estimate_seconds') or 0)
                if estimate > 0 and r['running_for'] is not None:
                    job['progress'] = min(MOTIF_CAP, r['running_for'] / estimate)
                    job['estimate_seconds'] = round(estimate)
                job['step'] = 'Looking through the genome'
            elif r['progress'] is not None:
                job['progress'] = float(r['progress'])
        elif r['state'] == 'done':
            job['progress'] = 1.0
            job['track_id'] = result.get('track_id') or (result.get('id') if r['kind'] == 'motif' else None)
            job['files'] = result.get('files') or []
        out.append(job)
    return out


def _conversion_jobs(conn, member: str) -> List[dict]:
    rows = conn.execute(
        'SELECT c.id, c.track_id, c.kind, c.state, c.error, c.created_at, c.started_at, c.finished_at, '
        '  t.name AS track_name, t.genome_id, g.name AS genome_name, '
        "  (julianday(coalesce(c.finished_at, 'now')) - julianday(coalesce(c.started_at, c.created_at))) * 86400 "
        '    AS took '
        'FROM conversion_jobs c JOIN tracks t ON t.id = c.track_id '
        'LEFT JOIN reference_genomes g ON g.id = t.genome_id '
        "WHERE t.created_by = ? AND c.created_at >= datetime('now', ?) "
        'ORDER BY c.created_at DESC LIMIT ?',
        (member, f'-{analysis.KEEP_DAYS} days', LIMIT)).fetchall()
    out = []
    for r in rows:
        job = {'id': r['id'], 'source': 'upload', 'kind': r['kind'], 'what': KINDS.get(r['kind'], 'Preparing an upload'),
               'about': r['track_name'], 'genome': {'id': r['genome_id'], 'name': r['genome_name'] or 'a genome'},
               'state': r['state'], 'step': None, 'progress': None, 'ahead': 0, 'error': r['error'],
               'created_at': r['created_at'], 'started_at': r['started_at'], 'finished_at': r['finished_at'],
               'seconds': round(r['took']) if r['took'] is not None else None, 'track_id': r['track_id']}
        if r['state'] == 'queued':
            job['progress'] = 0.0
            job['ahead'] = conn.execute("SELECT count(*) FROM conversion_jobs WHERE state IN ('queued', 'running') "
                                        'AND created_at < ?', (r['created_at'],)).fetchone()[0]
        elif r['state'] == 'running':
            job['step'] = 'Turning the file into a coverage track'   # no measure of how far: the bar only moves
        elif r['state'] == 'done':
            job['progress'] = 1.0
        out.append(job)
    return out


def _track_names(conn, rows) -> Dict[str, str]:
    ids = set()
    for r in rows:
        p = _json(r['params'])
        ids.update(str(i) for i in (p.get('tracks') or []))
        ids.update(str(p[k]) for k in ('a', 'b') if p.get(k))
    if not ids:
        return {}
    marks = ','.join('?' * len(ids))
    return {r['id']: r['name'] for r in conn.execute(f'SELECT id, name FROM tracks WHERE id IN ({marks})', list(ids))}


def _about(kind: str, params: dict, result: dict, names: Dict[str, str]) -> str:
    """One line on what the job is of: its tracks, and the features or name."""
    tracks = [names.get(str(t), 'a track no longer here') for t in (params.get('tracks') or [])]
    listed = ', '.join(tracks[:3]) + (f' and {len(tracks) - 3} more' if len(tracks) > 3 else '')
    if kind == 'compare':
        if params.get('name'):
            return str(params['name'])
        return f"{names.get(str(params.get('a')), 'A')} against {names.get(str(params.get('b')), 'B')}"
    if kind == 'motif':
        return str(params.get('name') or 'sequence')
    if kind in ('profile', 'genes'):
        over = ', '.join(result.get('regions') or []) or \
            ', '.join(r.get('type') or 'features' for r in params.get('regions') or [])
        return f'{listed} over {over}' if over else listed
    return listed


def _open_url(mount: str, job: dict) -> Optional[str]:
    """Where to see a finished job's result: the browser on its genome, showing it."""
    if job['state'] != 'done':
        return None
    if job.get('track_id'):
        return f"{mount}/?track={job['track_id']}"
    if job['source'] == 'analysis':
        return f"{mount}/?job={job['id']}"
    return None


def _json(text) -> dict:
    try:
        value = json.loads(text) if text else {}
    except ValueError:
        return {}
    return value if isinstance(value, dict) else {}


def _static() -> str:
    prefix = str(current_app.config.get('URL_PREFIX', '/genome')).rstrip('/')
    return f'{prefix}/static' if prefix else '/static'


def _mount() -> str:
    return str(current_app.config.get('URL_PREFIX', '/genome')).rstrip('/')


# --- the page -----------------------------------------------------------------------------

def _page_html(static: str, mount: str, embed: bool) -> str:
    header = '' if embed else f"""
  <header class="tt-bar">
    <a class="tt-home" href="/" title="All of Trellis Tools"><img src="{static}/trellis/trellis.svg" alt=""><span>Trellis Tools</span></a>
    <a class="tt-module" href="{mount}/"><img src="{static}/trellis/genome.svg" alt="">Genome Browser</a><span class="tt-sub">My jobs</span>
    <nav class="tt-links"><a href="{mount}/jobs">My jobs</a><a href="{mount}/help">Help</a><a href="{mount}/about">About</a></nav>
  </header>"""
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>My jobs · Genome Browser</title>
  <link rel="icon" href="{static}/trellis/genome.svg" type="image/svg+xml">
  <link rel="stylesheet" href="{static}/trellis/trellis.css">
  <link rel="stylesheet" href="{static}/jobs.css">
</head>
<body class="tt{' embedded' if embed else ''}" data-mount="{mount}">{header}
  <main class="jobs">
    <div class="jobs-head">
      <h1>My jobs</h1>
      <p class="jobs-hint">Analyses and uploads you started in the Genome Browser: only yours, and nobody else sees
        them here. The genome worker runs one job at a time for everyone, in the order they were asked for, so a
        waiting job says how many are ahead of it (other people's too, as a number only). Finished ones are kept for
        <span id="keep-days">7</span> days.</p>
    </div>
    <p id="jobs-error" class="jobs-error" role="alert" hidden></p>
    <ul id="jobs-list" class="jobs-list" aria-live="polite"></ul>
    <p id="jobs-empty" class="jobs-empty" hidden>Nothing yet. Start an analysis from <b>Analyse</b> in the
      <a href="{mount}/" target="_top">Genome Browser</a>, and it shows here while it runs.</p>
  </main>
  <script src="{static}/jobs.js" defer></script>
</body>
</html>
"""
