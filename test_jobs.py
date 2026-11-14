"""My jobs: each member sees their own jobs, with how far each has got and where to open it."""
import json
import sqlite3
import uuid

import pytest

import analysis
import db
from test_analysis import world  # noqa: F401  (the same genome, tracks and members)


def _job(world, member, kind, state, params=None, result=None, **extra):
    conn = db.connect(world['root'] / 'genome.db')
    analysis.ensure_schema(conn)
    job_id = uuid.uuid4().hex
    conn.execute('INSERT INTO analysis_jobs (id, genome_id, member, owner_lab, kind, params, state, result, created_at, '
                 'started_at, finished_at, step, progress) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                 (job_id, world['genome'], member, 'labA', kind, json.dumps(params or {}), state,
                  json.dumps(result) if result is not None else None,
                  extra.get('created_at', "2099-01-01 00:00:00"), extra.get('started_at'), extra.get('finished_at'),
                  extra.get('step'), extra.get('progress')))
    conn.commit()
    conn.close()
    return job_id


def _now(conn, offset):
    return conn.execute("SELECT datetime('now', ?)", (offset,)).fetchone()[0]


def test_my_jobs_with_progress_and_where_to_open_them(world):
    conn = db.connect(world['root'] / 'genome.db')
    t = {k: _now(conn, k) for k in ('-50 seconds', '-40 seconds', '-30 seconds', '-20 seconds', '-10 seconds', '-5 seconds')}
    conn.close()
    ids = world['ids']
    done = _job(world, 'alice', 'profile', 'done', {'tracks': [ids['wt']], 'regions': [{'type': 'gene'}]},
                {'files': ['profile.svg'], 'regions': ['gene (3)']}, created_at=t['-50 seconds'],
                started_at=t['-50 seconds'], finished_at=t['-40 seconds'])
    made = _job(world, 'alice', 'compare', 'done', {'a': ids['wt'], 'b': ids['mut']}, {'track_id': 'newtrack', 'files': []},
                created_at=t['-40 seconds'], started_at=t['-40 seconds'], finished_at=t['-30 seconds'])
    running = _job(world, 'alice', 'correlation', 'running', {'tracks': [ids['wt'], ids['mut']]},
                   created_at=t['-30 seconds'], started_at=t['-20 seconds'], step='Summarising the tracks', progress=0.42)
    waiting = _job(world, 'alice', 'genes', 'queued', {'tracks': [ids['wt']]}, created_at=t['-10 seconds'])
    failed = _job(world, 'alice', 'profile', 'failed', {'tracks': [ids['wt']]}, created_at=t['-5 seconds'])
    _job(world, 'bob', 'profile', 'queued', {}, created_at=t['-20 seconds'])     # ahead of alice's, but not hers to see

    res = world['call']('GET', '/api/jobs', 'alice')
    assert res.status_code == 200
    body = res.get_json()
    assert body['busy'] is True
    jobs = {j['id']: j for j in body['jobs']}
    assert set(jobs) == {done, made, running, waiting, failed}
    # Running first, then waiting, then the rest, newest first.
    assert [j['id'] for j in body['jobs']] == [running, waiting, failed, made, done]

    r = jobs[running]
    assert r['what'] == 'Replicate check' and r['about'] == 'WT, MUT' and r['genome']['name'] == 'Plant'
    assert r['progress'] == pytest.approx(0.42) and r['step'] == 'Summarising the tracks' and r['open_url'] is None
    w = jobs[waiting]
    assert w['progress'] == 0 and w['ahead'] == 2          # the running one, and bob's queued before it
    d = jobs[done]
    assert d['progress'] == 1 and d['open_url'] == f'/genome/?job={done}' and d['about'] == 'WT over gene (3)'
    assert d['seconds'] == 10
    assert jobs[made]['open_url'] == '/genome/?track=newtrack' and jobs[made]['about'] == 'WT against MUT'
    assert jobs[failed]['open_url'] is None

    bob = world['call']('GET', '/api/jobs', 'bob').get_json()
    assert [j['kind'] for j in bob['jobs']] == ['profile'] and bob['jobs'][0]['ahead'] == 1
    assert world['call']('GET', '/api/jobs', '').status_code == 401


def test_marking_sequences_runs_against_its_estimate(world):
    conn = db.connect(world['root'] / 'genome.db')
    started = _now(conn, '-30 seconds')
    conn.close()
    _job(world, 'alice', 'motif', 'running', {'name': 'my primers', 'estimate_seconds': 60},
         created_at=started, started_at=started)
    job = world['call']('GET', '/api/jobs', 'alice').get_json()['jobs'][0]
    assert job['about'] == 'my primers' and job['estimate_seconds'] == 60
    assert job['progress'] == pytest.approx(0.5, abs=0.05)
    # Long past its estimate, it still does not claim to be done.
    conn = db.connect(world['root'] / 'genome.db')
    conn.execute("UPDATE analysis_jobs SET started_at = datetime('now', '-1 hour')")
    conn.commit()
    conn.close()
    assert world['call']('GET', '/api/jobs', 'alice').get_json()['jobs'][0]['progress'] == pytest.approx(0.95)


def test_uploads_being_prepared(world):
    conn = db.connect(world['root'] / 'genome.db')
    conn.execute('INSERT INTO conversion_jobs (id, track_id, kind, input_path, state, started_at) '
                 "VALUES ('c1', ?, 'bam', '/x.bam', 'running', datetime('now'))", (world['ids']['wt'],))
    conn.commit()
    conn.close()
    jobs = world['call']('GET', '/api/jobs', 'alice').get_json()['jobs']
    assert len(jobs) == 1
    j = jobs[0]
    assert j['source'] == 'upload' and j['about'] == 'WT' and j['progress'] is None   # no measure: the bar only moves
    assert world['call']('GET', '/api/jobs', 'bob').get_json()['jobs'] == []


def test_the_page_and_the_dashboard_tab(world):
    page = world['call']('GET', '/jobs', 'alice')
    assert page.status_code == 200 and b'My jobs' in page.data and b'tt-bar' in page.data
    embedded = world['call']('GET', '/jobs?embed=1', 'alice')
    assert b'tt-bar' not in embedded.data and b'class="tt embedded"' in embedded.data
    tab = world['call']('GET', '/suite-manifest.json').get_json()['tabs'][0]
    assert tab == {'id': 'genome-jobs', 'title': 'My jobs', 'icon': '⏳', 'section': 'me', 'src': '/genome/jobs?embed=1'}


def test_an_older_database_gains_the_progress_columns(tmp_path):
    conn = sqlite3.connect(tmp_path / 'old.db')
    conn.execute('CREATE TABLE analysis_jobs (id TEXT PRIMARY KEY, genome_id TEXT NOT NULL, member TEXT NOT NULL, '
                 "owner_lab TEXT, kind TEXT NOT NULL, params TEXT NOT NULL, state TEXT NOT NULL DEFAULT 'queued', "
                 "error TEXT, result TEXT, created_at TEXT NOT NULL DEFAULT (datetime('now')), started_at TEXT, "
                 'finished_at TEXT)')
    conn.execute("INSERT INTO analysis_jobs (id, genome_id, member, kind, params) VALUES ('j', 'g', 'm', 'profile', '{}')")
    analysis.ensure_schema(conn)
    analysis.ensure_schema(conn)                   # and again, as every request does
    cols = {r[1] for r in conn.execute('PRAGMA table_info(analysis_jobs)')}
    assert {'step', 'progress'} <= cols
    assert conn.execute("SELECT step, progress FROM analysis_jobs WHERE id = 'j'").fetchone() == (None, None)


def test_progress_is_recorded_as_a_job_goes(world, monkeypatch):
    seen = []
    conn = db.connect(world['root'] / 'genome.db')
    analysis.ensure_schema(conn)
    job_id = _job(world, 'alice', 'compare', 'queued', {'a': world['ids']['wt'], 'b': world['ids']['mut']})
    job = analysis.claim_next(conn)
    assert job['id'] == job_id
    assert conn.execute('SELECT progress FROM analysis_jobs WHERE id = ?', (job_id,)).fetchone()[0] == 0

    def fake_compare(params, a, b, work):
        analysis._progress('Comparing the tracks', 0.0)
        seen.append(conn.execute('SELECT step, progress FROM analysis_jobs WHERE id = ?', (job_id,)).fetchone())
        (work / 'compare.bw').write_bytes(b'')
        return {'files': ['compare.bw'], 'operation': 'A minus B'}
    monkeypatch.setattr(analysis, 'run_compare', fake_compare)
    monkeypatch.setattr(analysis, '_add_track', lambda *a, **k: 'newtrack')
    analysis.run_job(conn, world['root'], job)
    assert tuple(seen[0]) == ('Comparing the tracks', 0.0)
    row = conn.execute('SELECT state, step, progress FROM analysis_jobs WHERE id = ?', (job_id,)).fetchone()
    assert tuple(row) == ('done', None, 1.0)
    assert analysis._PROGRESS[0] is None          # nothing left listening once the job is over
    conn.close()
