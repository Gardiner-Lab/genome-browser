"""
Tests for the durable conversion queue (task 8.1, ``conversion_queue.py``).

These exercise the queue logic directly against a temporary ``genome.db``: a
committed ``enqueue`` row survives a reconnect (durability, Req 13.3/13.20),
``claim_next_job`` takes the oldest queued row and flips exactly one to
``running`` — including under a simulated double-claim where only one caller may
win (concurrency 1, Req 13.12) — and ``recover_running`` re-enqueues a stranded
``running`` row on startup (Req 13.18/13.20).
"""

import threading
import time

import pytest

import conversion_queue
import db


@pytest.fixture
def conn(tmp_path):
    """An initialised temp genome.db with one genome + track to key jobs to."""
    db_path = tmp_path / 'genome.db'
    db.init_db(db_path)
    connection = db.connect(db_path)
    # A genome and a converting track so conversion_jobs' FK to tracks holds.
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


def _states(connection):
    rows = connection.execute(
        'SELECT id, state FROM conversion_jobs').fetchall()
    return {r['id']: r['state'] for r in rows}


# --- enqueue: durable queued row (Req 13.3, 13.20) -------------------------


def test_enqueue_commits_queued_row_that_survives_reconnect(conn, tmp_path):
    job_id = conversion_queue.enqueue(conn, 't1', 'bam', '/scratch/x/input.bam')

    # A fresh connection sees the row: enqueue committed, so an accepted job is
    # durable the moment the response would return.
    other = db.connect(tmp_path / 'genome.db')
    try:
        row = other.execute(
            'SELECT id, track_id, kind, input_path, state, attempts '
            'FROM conversion_jobs WHERE id = ?', (job_id,)).fetchone()
    finally:
        other.close()
    assert row is not None
    assert row['track_id'] == 't1'
    assert row['kind'] == 'bam'
    assert row['input_path'] == '/scratch/x/input.bam'
    assert row['state'] == 'queued'
    assert row['attempts'] == 0


# --- claim_next_job: oldest first, stamps running (Req 13.12) --------------


def test_claim_returns_none_on_empty_queue(conn):
    assert conversion_queue.claim_next_job(conn) is None


def test_claim_takes_oldest_queued_row(conn):
    # Insert three jobs with explicit increasing created_at so order is defined.
    for i, ts in enumerate(
            ['2024-01-01 00:00:01', '2024-01-01 00:00:02',
             '2024-01-01 00:00:03'], start=1):
        conn.execute(
            'INSERT INTO conversion_jobs '
            '(id, track_id, kind, input_path, state, created_at) '
            "VALUES (?, ?, 'bam', ?, 'queued', ?)",
            (f'j{i}', f't{i}', f'/s/j{i}/input.bam', ts))
    conn.commit()

    claimed = conversion_queue.claim_next_job(conn)
    assert claimed is not None
    assert claimed['id'] == 'j1'            # the oldest
    assert claimed['state'] == 'running'
    assert claimed['attempts'] == 1
    assert claimed['started_at'] is not None

    # The next claim takes the second-oldest; j1 stays running.
    second = conversion_queue.claim_next_job(conn)
    assert second['id'] == 'j2'
    states = _states(conn)
    assert states == {'j1': 'running', 'j2': 'running', 'j3': 'queued'}


def test_claim_drains_queue_then_returns_none(conn):
    ids = []
    for i in range(1, 4):
        jid = conversion_queue.enqueue(conn, f't{i}', 'wig', f'/s/j{i}/input.wig')
        ids.append(jid)
        time.sleep(0.001)
    claimed = []
    while True:
        job = conversion_queue.claim_next_job(conn)
        if job is None:
            break
        claimed.append(job['id'])
    # Every enqueued job was claimed exactly once, none left queued.
    assert sorted(claimed) == sorted(ids)
    assert all(s == 'running' for s in _states(conn).values())


# --- concurrency 1: a double-claim gives the row to exactly one caller -----


def test_simulated_double_claim_claims_each_row_once(tmp_path):
    """
    Two threads racing claim_next_job over the same queue must never both get
    the same job: each queued row is claimed by exactly one thread (Req 13.12).
    """
    db_path = tmp_path / 'genome.db'
    db.init_db(db_path)
    setup = db.connect(db_path)
    setup.execute(
        'INSERT INTO reference_genomes '
        '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
        "VALUES ('g1', 'G', 'asm', 'labA', 'owner_lab_only', 'alice')")
    n_jobs = 20
    for i in range(n_jobs):
        setup.execute(
            'INSERT INTO tracks (id, genome_id, name, type, owner_lab, '
            'track_status, created_by) VALUES (?, ?, ?, ?, ?, ?, ?)',
            (f't{i}', 'g1', f'T{i}', 'bigwig', 'labA', 'converting', 'alice'))
        setup.execute(
            'INSERT INTO conversion_jobs '
            '(id, track_id, kind, input_path, state, created_at) '
            "VALUES (?, ?, 'bam', ?, 'queued', ?)",
            (f'j{i}', f't{i}', f'/s/j{i}/input.bam',
             f'2024-01-01 00:00:{i:02d}'))
    setup.commit()
    setup.close()

    claimed_by = {'A': [], 'B': []}
    barrier = threading.Barrier(2)

    def worker(label):
        connection = db.connect(db_path)
        try:
            barrier.wait()
            while True:
                job = conversion_queue.claim_next_job(connection)
                if job is None:
                    # Queue may just be momentarily contended; retry briefly,
                    # then stop when truly drained.
                    remaining = connection.execute(
                        "SELECT COUNT(*) AS c FROM conversion_jobs "
                        "WHERE state = 'queued'").fetchone()['c']
                    if remaining == 0:
                        break
                    continue
                claimed_by[label].append(job['id'])
        finally:
            connection.close()

    threads = [threading.Thread(target=worker, args=(lbl,))
               for lbl in ('A', 'B')]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)

    all_claimed = claimed_by['A'] + claimed_by['B']
    # No job claimed twice, and every job claimed exactly once.
    assert len(all_claimed) == len(set(all_claimed)), 'a job was claimed twice'
    assert set(all_claimed) == {f'j{i}' for i in range(n_jobs)}

    verify = db.connect(db_path)
    try:
        running = verify.execute(
            "SELECT COUNT(*) AS c FROM conversion_jobs "
            "WHERE state = 'running'").fetchone()['c']
        # attempts incremented exactly once per row (claimed once).
        bad = verify.execute(
            'SELECT COUNT(*) AS c FROM conversion_jobs '
            'WHERE attempts != 1').fetchone()['c']
    finally:
        verify.close()
    assert running == n_jobs
    assert bad == 0


# --- mark_done / mark_failed close a job (Req 13.5, 13.20) -----------------


def test_mark_done_sets_done_and_finished_at(conn):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j/input.bam')
    conversion_queue.claim_next_job(conn)
    conversion_queue.mark_done(conn, jid)
    row = conn.execute(
        'SELECT state, error, finished_at FROM conversion_jobs WHERE id = ?',
        (jid,)).fetchone()
    assert row['state'] == 'done'
    assert row['error'] is None
    assert row['finished_at'] is not None


def test_mark_failed_records_reason(conn):
    jid = conversion_queue.enqueue(conn, 't1', 'bam', '/s/j/input.bam')
    conversion_queue.claim_next_job(conn)
    conversion_queue.mark_failed(conn, jid, 'chrom Chr9 not found')
    row = conn.execute(
        'SELECT state, error, finished_at FROM conversion_jobs WHERE id = ?',
        (jid,)).fetchone()
    assert row['state'] == 'failed'
    assert row['error'] == 'chrom Chr9 not found'
    assert row['finished_at'] is not None


# --- recover_running re-enqueues stranded jobs on startup (Req 13.18, 13.20)


def test_recover_running_requeues_only_running_rows(conn):
    q = conversion_queue.enqueue(conn, 't1', 'bam', '/s/q/input.bam')
    r = conversion_queue.enqueue(conn, 't2', 'bam', '/s/r/input.bam')
    d = conversion_queue.enqueue(conn, 't3', 'wig', '/s/d/input.wig')
    # Claim r (-> running) and mark d done; q stays queued.
    conn.execute("UPDATE conversion_jobs SET state = 'running' WHERE id = ?", (r,))
    conn.execute("UPDATE conversion_jobs SET state = 'done' WHERE id = ?", (d,))
    conn.commit()

    reset = conversion_queue.recover_running(conn)
    assert reset == 1
    states = _states(conn)
    assert states[q] == 'queued'   # untouched
    assert states[r] == 'queued'   # recovered back to queued
    assert states[d] == 'done'     # untouched
    # The recovered row is claimable again.
    row = conn.execute(
        'SELECT started_at FROM conversion_jobs WHERE id = ?', (r,)).fetchone()
    assert row['started_at'] is None
