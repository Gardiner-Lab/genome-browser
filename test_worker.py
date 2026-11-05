"""
Tests for the Conversion_Worker loop and startup recovery (task 9.2, ``worker.py``).

These drive the two testable entry points directly against a temporary
``genome.db`` with ``conversion.convert_job`` and ``conversion.has_free_disk_floor``
stubbed, so no real toolchain or disk is touched:

* :func:`worker.recover` re-enqueues a stranded ``running`` job (Req 13.18/13.20)
  and rescues an orphaned ``converting`` track — re-enqueue when its scratch
  input still exists, fail with a descriptive error when it is gone (Req 13.18).
* :func:`worker.run_once` claims the oldest queued job, calls ``convert_job``, and
  marks the job ``done``/``failed`` from the outcome (Req 13.12/13.16); defers
  without claiming when the free-disk floor is below the floor (Req 13.13); and
  returns ``False`` on an empty queue so ``main`` sleeps (Req 13.17).
"""

from pathlib import Path

import pytest

import conversion
import conversion_queue
import db
import worker


@pytest.fixture
def conn(tmp_path):
    """An initialised temp genome.db with a genome and a couple of tracks."""
    db_path = tmp_path / 'genome.db'
    db.init_db(db_path)
    connection = db.connect(db_path)
    connection.execute(
        'INSERT INTO reference_genomes '
        '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
        "VALUES ('g1', 'G', 'asm', 'labA', 'owner_lab_only', 'alice')")
    connection.executemany(
        'INSERT INTO tracks (id, genome_id, name, type, owner_lab, '
        'track_status, created_by) VALUES (?, ?, ?, ?, ?, ?, ?)',
        [(f't{i}', 'g1', f'Track {i}', 'bigwig', 'labA', 'converting', 'alice')
         for i in range(1, 5)])
    connection.commit()
    yield connection
    connection.close()


@pytest.fixture
def config(tmp_path):
    """A WorkerConfig pointing at the temp data root, floor disabled by default."""
    data_root = tmp_path
    (data_root / 'scratch').mkdir(exist_ok=True)
    return worker.WorkerConfig(
        data_root=data_root, samtools_threads=2, min_free_gb=0)


def _states(connection):
    return {r['id']: r['state']
            for r in connection.execute(
                'SELECT id, state FROM conversion_jobs').fetchall()}


def _track_status(connection, track_id):
    row = connection.execute(
        'SELECT track_status, error FROM tracks WHERE id = ?',
        (track_id,)).fetchone()
    return row['track_status'], row['error']


def _only_converting(connection, *track_ids):
    """
    Leave only ``track_ids`` in 'converting'; set the rest 'ready'.

    The fixture seeds several converting tracks so FKs hold; a recover test
    wants only its subject track(s) to be the orphan(s) under examination, so
    the other converting tracks are moved to a terminal state first.
    """
    keep = set(track_ids)
    for row in connection.execute(
            "SELECT id FROM tracks WHERE track_status = 'converting'").fetchall():
        if row['id'] not in keep:
            connection.execute(
                "UPDATE tracks SET track_status = 'ready' WHERE id = ?",
                (row['id'],))
    connection.commit()


# --- load_config reads the same env vars app.py's Config reads ------------


def test_load_config_reads_env(monkeypatch):
    monkeypatch.setenv('GENOME_DATA_ROOT', '/data/genome')
    monkeypatch.setenv('SAMTOOLS_THREADS', '4')
    monkeypatch.setenv('CONVERT_MIN_FREE_GB', '25')
    cfg = worker.load_config()
    assert cfg.data_root == Path('/data/genome')
    assert cfg.samtools_threads == 4
    assert cfg.min_free_gb == 25
    assert cfg.scratch_root == Path('/data/genome') / 'scratch'


def test_load_config_defaults(monkeypatch):
    for var in ('GENOME_DATA_ROOT', 'SAMTOOLS_THREADS', 'CONVERT_MIN_FREE_GB'):
        monkeypatch.delenv(var, raising=False)
    cfg = worker.load_config()
    assert cfg.data_root == Path('/genome')
    assert cfg.samtools_threads == 2
    assert cfg.min_free_gb == 50


# --- recover: re-enqueue a running job (Req 13.18, 13.20) -----------------


def test_recover_requeues_running_job(conn):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j/input.bam')
    conn.execute("UPDATE conversion_jobs SET state = 'running' WHERE id = ?",
                 (jid,))
    conn.commit()

    summary = worker.recover(conn, Path('/unused'))

    assert summary['running_requeued'] == 1
    assert _states(conn)[jid] == 'queued'


# --- recover: orphaned converting track with input present -> re-enqueue ---


def test_recover_reenqueues_orphan_track_when_input_exists(conn, tmp_path):
    # t2 is converting but has no job row at all (its durable job was lost),
    # yet an earlier job record names a scratch input that still exists.
    _only_converting(conn, 't2')
    scratch = tmp_path / 'scratch' / 'job2'
    scratch.mkdir(parents=True)
    input_path = scratch / 'input.bam'
    input_path.write_bytes(b'BAMDATA')
    # A terminal (done) job record carrying the kind + input path, but nothing
    # queued/running -> t2 is an orphan the queue alone would never touch.
    conn.execute(
        'INSERT INTO conversion_jobs '
        '(id, track_id, kind, input_path, state) '
        "VALUES ('oldjob', 't2', 'bam', ?, 'done')", (str(input_path),))
    conn.commit()

    summary = worker.recover(conn, tmp_path)

    assert summary['re_enqueued'] == 1
    assert summary['failed'] == 0
    # A fresh queued job for t2 now exists, pointing at the surviving input.
    row = conn.execute(
        "SELECT kind, input_path, state FROM conversion_jobs "
        "WHERE track_id = 't2' AND state = 'queued'").fetchone()
    assert row is not None
    assert row['kind'] == 'bam'
    assert row['input_path'] == str(input_path)
    # The track is still converting (it will be converted), not failed.
    assert _track_status(conn, 't2')[0] == 'converting'


# --- recover: orphaned converting track with input gone -> fail -----------


def test_recover_fails_orphan_track_when_input_missing(conn, tmp_path):
    # t3 is converting with a terminal job whose scratch input no longer exists.
    _only_converting(conn, 't3')
    conn.execute(
        'INSERT INTO conversion_jobs '
        '(id, track_id, kind, input_path, state) '
        "VALUES ('gonejob', 't3', 'bam', ?, 'failed')",
        (str(tmp_path / 'scratch' / 'gone' / 'input.bam'),))
    conn.commit()

    summary = worker.recover(conn, tmp_path)

    assert summary['failed'] == 1
    status, error = _track_status(conn, 't3')
    assert status == 'failed'
    assert error  # a descriptive, non-empty reason
    assert 'resume' in error.lower() or 'scratch' in error.lower()
    # No new queued job was created for t3.
    assert conn.execute(
        "SELECT COUNT(*) AS c FROM conversion_jobs "
        "WHERE track_id = 't3' AND state = 'queued'").fetchone()['c'] == 0


def test_recover_leaves_track_with_active_job_untouched(conn, tmp_path):
    # t4 is converting and HAS a queued job -> not an orphan, must be left alone.
    _only_converting(conn, 't4')
    conversion_queue.enqueue(conn, 't4', 'wig', '/s/j4/input.wig')

    summary = worker.recover(conn, tmp_path)

    assert summary['re_enqueued'] == 0
    assert summary['failed'] == 0
    assert _track_status(conn, 't4')[0] == 'converting'


# --- run_once: claim + convert + mark_done on success (Req 13.12, 13.16) ---


def test_run_once_processes_job_and_marks_done(conn, config, monkeypatch):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j1/input.bam')
    calls = {}

    def fake_convert_job(*, track_id, genome_id, kind, input_path, data_root,
                         samtools_threads):
        calls.update(dict(track_id=track_id, genome_id=genome_id, kind=kind,
                          input_path=input_path,
                          samtools_threads=samtools_threads))
        return conversion.ConversionOutcome(status='ready',
                                            file_path='/genome/tracks/t1/t1.bw')

    monkeypatch.setattr(conversion, 'convert_job', fake_convert_job)
    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: True)

    processed = worker.run_once(conn, config.data_root, config)

    assert processed is True
    # convert_job was called with the job's fields and the track's genome_id.
    assert calls['track_id'] == 't1'
    assert calls['genome_id'] == 'g1'
    assert calls['kind'] == 'bam'
    assert calls['input_path'] == '/s/j1/input.bam'
    assert calls['samtools_threads'] == 2
    # The job was closed as done.
    row = conn.execute(
        'SELECT state, error, finished_at FROM conversion_jobs WHERE id = ?',
        (jid,)).fetchone()
    assert row['state'] == 'done'
    assert row['error'] is None
    assert row['finished_at'] is not None


# --- run_once: marks failed from a failed outcome (Req 13.5, 13.20) -------


def test_run_once_marks_failed_from_outcome(conn, config, monkeypatch):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j1/input.bam')

    monkeypatch.setattr(
        conversion, 'convert_job',
        lambda **kw: conversion.ConversionOutcome(
            status='failed', error='chrom Chr9 not found'))
    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: True)

    processed = worker.run_once(conn, config.data_root, config)

    assert processed is True
    row = conn.execute(
        'SELECT state, error FROM conversion_jobs WHERE id = ?',
        (jid,)).fetchone()
    assert row['state'] == 'failed'
    assert row['error'] == 'chrom Chr9 not found'


# --- run_once: claims the oldest queued job first -------------------------


def test_run_once_claims_oldest_first(conn, config, monkeypatch):
    for i, ts in enumerate(
            ['2024-01-01 00:00:01', '2024-01-01 00:00:02'], start=1):
        conn.execute(
            'INSERT INTO conversion_jobs '
            '(id, track_id, kind, input_path, state, created_at) '
            "VALUES (?, ?, 'bam', ?, 'queued', ?)",
            (f'j{i}', f't{i}', f'/s/j{i}/input.bam', ts))
    conn.commit()

    seen = []
    monkeypatch.setattr(
        conversion, 'convert_job',
        lambda **kw: seen.append(kw['track_id'])
        or conversion.ConversionOutcome(status='ready', file_path='x'))
    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: True)

    worker.run_once(conn, config.data_root, config)
    assert seen == ['t1']  # the oldest
    assert _states(conn)['j1'] == 'done'
    assert _states(conn)['j2'] == 'queued'


# --- run_once: defers (no claim) when below the free-disk floor (Req 13.13)


def test_run_once_defers_below_disk_floor(conn, config, monkeypatch):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j1/input.bam')
    called = {'convert': False}

    def must_not_run(**kw):
        called['convert'] = True
        raise AssertionError('convert_job must not be called below the floor')

    monkeypatch.setattr(conversion, 'convert_job', must_not_run)
    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: False)

    processed = worker.run_once(conn, config.data_root, config)

    assert processed is False
    assert called['convert'] is False
    # The job is left queued (not claimed / not running).
    assert _states(conn)[jid] == 'queued'


# --- run_once: empty queue returns False (main would then sleep, Req 13.17)


def test_run_once_returns_false_on_empty_queue(conn, config, monkeypatch):
    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: True)
    monkeypatch.setattr(
        conversion, 'convert_job',
        lambda **kw: pytest.fail('convert_job called with empty queue'))

    assert worker.run_once(conn, config.data_root, config) is False


# --- run_once: a deleted track closes the job failed, no convert ----------


def test_run_once_fails_job_when_track_missing(conn, config, monkeypatch):
    # A job whose track row is absent (belt-and-braces: FK cascade usually
    # removes the job too). run_once must close it failed without calling convert.
    conn.execute(
        'INSERT INTO conversion_jobs '
        '(id, track_id, kind, input_path, state) '
        "VALUES ('orphanjob', 't1', 'bam', '/s/jx/input.bam', 'queued')")
    conn.commit()
    # Drop the track out from under the job without cascading to the job row,
    # leaving the job pointing at a missing track.
    conn.execute('PRAGMA foreign_keys = OFF')
    conn.execute("DELETE FROM tracks WHERE id = 't1'")
    conn.commit()
    conn.execute('PRAGMA foreign_keys = ON')

    monkeypatch.setattr(conversion, 'has_free_disk_floor',
                        lambda *a, **k: True)
    monkeypatch.setattr(
        conversion, 'convert_job',
        lambda **kw: pytest.fail('convert_job called for a missing track'))

    processed = worker.run_once(conn, config.data_root, config)

    assert processed is True
    row = conn.execute(
        "SELECT state, error FROM conversion_jobs WHERE id = 'orphanjob'"
    ).fetchone()
    assert row['state'] == 'failed'
    assert row['error']
