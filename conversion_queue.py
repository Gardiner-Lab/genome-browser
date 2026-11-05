"""
The durable conversion queue over the ``conversion_jobs`` table.

A BAM/wig upload is never converted inside the request. The web app streams the
original into ``scratch/<job_id>/`` and commits a ``conversion_jobs`` row in the
``queued`` state *before* answering the member (:func:`enqueue`); the single
Conversion_Worker (task 9.2) later claims the oldest queued row
(:func:`claim_next_job`), calls ``conversion.convert_job``, and closes the job
out (:func:`mark_done` / :func:`mark_failed`). Because the row is committed to
the WAL database before the HTTP response returns, an accepted job survives a
crash or restart of either container — nothing lives only in memory or in the
member's session (Requirements 13.3, 13.17, 13.20).

The module is named ``conversion_queue`` rather than ``queue`` on purpose: the
module runs with its own directory on ``sys.path`` (flat imports like ``import
db``), so a file named ``queue.py`` here would shadow the Python standard
library's ``queue`` module that ``requests``/``urllib3`` (pulled in by
``notebook.py``) import, breaking the app at startup. ``conversion_queue`` keeps
the durable-queue logic in its own module without that collision.

Two properties make this a correct queue rather than just a table:

* **Serialized claim, concurrency 1.** :func:`claim_next_job` is a single
  ``UPDATE ... WHERE id = (SELECT ... ORDER BY created_at LIMIT 1)`` run inside
  an ``IMMEDIATE`` transaction. SQLite grants one writer at a time, so even if
  two workers (or a racing double-claim) run it at once, exactly one row flips
  to ``running`` and the loser's ``UPDATE`` matches nothing — the two can never
  grab the same job (Requirement 13.12).

* **Startup recovery.** :func:`recover_running` flips any ``running`` row left
  behind by a crash back to ``queued`` so the worker picks it up again, so no
  Track is stranded in ``converting`` with no worker on it (Requirements 13.18,
  13.20). The worker calls it once on startup before its first claim.

All writes go through the module's own WAL database via
``db.connect(data_root / 'genome.db')``, the one coordination point the app and
worker share, exactly like the rest of the module.
"""

import sqlite3
import uuid
from typing import Optional


def enqueue(conn: sqlite3.Connection, track_id: str, kind: str,
            input_path: str) -> str:
    """
    Commit one ``queued`` conversion job for a just-written scratch input.

    Called from the upload path once the original BAM/wig has been fully
    streamed to ``input_path`` under ``scratch/<job_id>/`` and the ``converting``
    track row exists. A fresh uuid job id is minted, a ``conversion_jobs`` row is
    inserted in the ``queued`` state keyed to ``track_id``, and the connection is
    **committed** before returning so the accepted job is durable the moment the
    HTTP response is sent — it survives a crash or restart with no worker running
    yet (Requirements 13.3, 13.20).

    ``kind`` is ``'bam'`` or ``'wig'``; ``input_path`` is the absolute path to
    the streamed original the worker will convert. Returns the job id.
    """
    job_id = uuid.uuid4().hex
    conn.execute(
        'INSERT INTO conversion_jobs (id, track_id, kind, input_path, state) '
        "VALUES (?, ?, ?, ?, 'queued')",
        (job_id, track_id, kind, str(input_path)))
    conn.commit()
    return job_id


def claim_next_job(conn: sqlite3.Connection) -> Optional[sqlite3.Row]:
    """
    Claim the oldest ``queued`` job, flipping exactly one row to ``running``.

    The claim is a single ``UPDATE`` whose ``WHERE`` names the one oldest queued
    row (``ORDER BY created_at, id LIMIT 1``), run inside an ``IMMEDIATE``
    transaction so SQLite takes the write lock up front. The row is set
    ``running``, its ``attempts`` incremented, and ``started_at`` stamped. Under
    a racing double-claim only one caller's ``UPDATE`` matches the row (the other
    runs against a row that is no longer ``queued`` and changes nothing), so a
    job can never be claimed twice — the concurrency-1 guarantee (Requirement
    13.12).

    Returns the claimed job row (id, track_id, kind, input_path, attempts, ...)
    or ``None`` when the queue holds no ``queued`` job. On ``None`` the
    transaction is rolled back so no lock is held.
    """
    # Begin an immediate write transaction: SQLite acquires the RESERVED lock
    # now rather than lazily, so two concurrent claims serialize here and only
    # one proceeds to mutate the chosen row. The isolation level is switched to
    # manual for this one transaction and restored in the finally, so the
    # connection is left exactly as it was found.
    previous_isolation = conn.isolation_level
    conn.isolation_level = None  # manage BEGIN/COMMIT/ROLLBACK explicitly
    try:
        conn.execute('BEGIN IMMEDIATE')
        try:
            # Pick the one oldest queued row inside the write lock, then flip
            # exactly that row by id. Holding the RESERVED lock makes the
            # select-then-update atomic against any other claimer: a racing
            # claim blocks on BEGIN IMMEDIATE until this COMMIT, so it can never
            # read the same queued row (it is already 'running' by then). The id
            # is captured here rather than re-derived from a post-UPDATE read,
            # which is ambiguous when several rows share the one-second
            # datetime('now') stamp.
            target = conn.execute(
                "SELECT id FROM conversion_jobs WHERE state = 'queued' "
                'ORDER BY created_at, id LIMIT 1').fetchone()
            if target is None:
                # Nothing queued: release the lock without writing.
                conn.execute('ROLLBACK')
                return None

            job_id = target['id']
            conn.execute(
                "UPDATE conversion_jobs "
                "SET state = 'running', "
                "    attempts = attempts + 1, "
                "    started_at = datetime('now') "
                "WHERE id = ?",
                (job_id,))
            row = conn.execute(
                'SELECT id, track_id, kind, input_path, state, attempts, '
                'error, created_at, started_at, finished_at '
                'FROM conversion_jobs WHERE id = ?', (job_id,)).fetchone()
            conn.execute('COMMIT')
            return row
        except sqlite3.Error:
            conn.execute('ROLLBACK')
            raise
    finally:
        conn.isolation_level = previous_isolation


def mark_done(conn: sqlite3.Connection, job_id: str) -> None:
    """
    Close a successful job: ``state='done'`` with ``finished_at`` stamped.

    The worker calls this after ``conversion.convert_job`` has driven the track
    row to ``ready`` and stored the BigWig. The job row is retained (not
    deleted) so its terminal state and timing stay visible independently of any
    session (Requirement 13.20).
    """
    conn.execute(
        "UPDATE conversion_jobs SET state = 'done', error = NULL, "
        "finished_at = datetime('now') WHERE id = ?",
        (job_id,))
    conn.commit()


def mark_failed(conn: sqlite3.Connection, job_id: str, error: str) -> None:
    """
    Close a failed job: ``state='failed'`` with the reason and ``finished_at``.

    Called after ``conversion.convert_job`` has set the track ``failed`` with
    its descriptive error; the same reason is recorded on the job row so the
    failure is legible from the queue as well as the track (Requirements 13.5,
    13.20).
    """
    conn.execute(
        "UPDATE conversion_jobs SET state = 'failed', error = ?, "
        "finished_at = datetime('now') WHERE id = ?",
        (error, job_id))
    conn.commit()


def recover_running(conn: sqlite3.Connection) -> int:
    """
    Re-enqueue every ``running`` job on startup, returning how many were reset.

    A ``running`` row at startup means the worker died mid-job (a clean finish
    leaves ``done``/``failed``). Flipping it back to ``queued`` — clearing
    ``started_at`` so it re-stamps on the next claim, while leaving ``attempts``
    as the record of how many tries it has had — lets the worker pick it up
    again, so no accepted conversion is silently lost and no Track is left in
    ``converting`` with no worker processing it (Requirements 13.18, 13.20). The
    worker calls this once before its first :func:`claim_next_job`.
    """
    cur = conn.execute(
        "UPDATE conversion_jobs SET state = 'queued', started_at = NULL "
        "WHERE state = 'running'")
    conn.commit()
    return cur.rowcount
