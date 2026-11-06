"""
The single serialized Conversion_Worker loop (the ``genome-worker`` container).

This module is the entry point compose runs as ``python worker.py``. It is the
one process that actually converts an uploaded BAM/wig into a BigWig: the web
app (``genome-browser``) only accepts the upload, writes the original into
``scratch/<job_id>/`` and commits a ``queued`` ``conversion_jobs`` row; this
worker, in its own container with its own CPU/memory limits, drains that durable
queue one job at a time. Running conversion in a separate process keeps gunicorn
responsive and makes serialization trivial — one container, one loop, concurrency
1 (Requirements 13.12, 13.16).

The worker owns *only* queue polling and the per-job hand-off; the toolchain run,
the track status writes, and the scratch cleanup all live in ``conversion.py``,
and the claim transaction lives in ``conversion_queue.py``. The worker reads
nothing but the durable SQLite database: there is no in-memory or session state,
so an accepted job keeps running and completes whether or not the uploading
member's browser is still open, and survives a restart of either container
(Requirements 13.17, 13.20).

Startup recovery (:func:`recover`, run once before the first claim) repairs the
two ways a crash can strand work:

* A ``conversion_jobs`` row left ``running`` means the worker died mid-job; it is
  flipped back to ``queued`` so the loop picks it up again
  (``conversion_queue.recover_running``, Requirements 13.18, 13.20).
* A track left ``converting`` with **no** active (``queued``/``running``) job is
  an orphan — its durable job was lost. If its scratch input still exists the
  worker re-enqueues it; if the input is gone it cannot be resumed, so the track
  is set ``failed`` with a descriptive error. Either way no track is left
  ``converting`` forever with no worker on it (Requirement 13.18).

The main loop (:func:`run_once` one step, :func:`main` forever) is deliberately
split so the step is testable without looping forever: a test drives
:func:`run_once` and :func:`recover` directly, stubbing ``conversion.convert_job``
and ``conversion.has_free_disk_floor``. Each step honours the free-disk floor
gate before claiming (leave the job queued and wait when below floor, Requirement
13.13), claims the oldest queued job, converts it, and records the job's terminal
state with ``mark_done``/``mark_failed``.

Run standalone (not via ``create_app``), so configuration is read straight from
the same environment variables ``app.py``'s ``Config`` reads — ``GENOME_DATA_ROOT``
(default ``/genome``), ``SAMTOOLS_THREADS`` (default 2), ``CONVERT_MIN_FREE_GB``
(default 50) — rather than importing the Flask app.
"""

import os
import time
from dataclasses import dataclass
from pathlib import Path

import conversion
import conversion_queue
import db

# How long to wait before polling again when the queue is empty or the
# free-disk floor is below the safety floor. Short enough that an accepted job
# starts promptly, long enough that an idle worker is not a busy-loop.
POLL_INTERVAL_SECONDS = 2.0


@dataclass
class WorkerConfig:
    """
    The worker's settings, read once from the environment.

    Mirrors the subset of ``app.py``'s ``Config`` the worker needs. Read here
    directly (not via ``create_app``) because the worker runs standalone as its
    own container command; the defaults match ``Config`` so the worker behaves
    identically to the web side on the shared data root.
    """

    data_root: Path
    samtools_threads: int
    min_free_gb: int

    @property
    def scratch_root(self) -> Path:
        """The scratch area the free-disk floor is measured against."""
        return self.data_root / 'scratch'


def load_config() -> WorkerConfig:
    """Read the worker's configuration from the environment (same vars as app.py)."""
    return WorkerConfig(
        data_root=Path(os.environ.get('GENOME_DATA_ROOT', '/genome')),
        samtools_threads=int(os.environ.get('SAMTOOLS_THREADS', '2')),
        min_free_gb=int(os.environ.get('CONVERT_MIN_FREE_GB', '50')),
    )


# ---------------------------------------------------------------------------
# Startup recovery (Requirements 13.18, 13.20)
# ---------------------------------------------------------------------------


def recover(conn, data_root) -> dict:
    """
    Repair work stranded by a crash, once, before the loop's first claim.

    Two independent repairs, so no accepted conversion is silently lost and no
    track is left ``converting`` with no worker on it (Requirement 13.18):

    1. **Re-enqueue running jobs.** ``conversion_queue.recover_running`` flips
       every ``running`` ``conversion_jobs`` row back to ``queued`` (a clean
       finish would have left ``done``/``failed``), so the worker re-claims it.
    2. **Rescue orphaned converting tracks.** A track still ``converting`` whose
       job row is gone (no ``queued`` or ``running`` job keyed to it) cannot be
       picked up by the queue alone. If its scratch input still exists the track
       is re-enqueued (a fresh ``queued`` job); if the input is gone the job
       cannot be resumed, so the track is set ``failed`` with a descriptive
       error. Step 1 runs first, so a track that *does* have a recovered running
       job is no longer an orphan here.

    Returns a small summary (``{'running_requeued', 're_enqueued', 'failed'}``)
    for logging. All writes go through the shared WAL database.
    """
    running_requeued = conversion_queue.recover_running(conn)

    data_root = Path(data_root)
    re_enqueued = 0
    failed = 0

    # Tracks still 'converting' that have no active (queued/running) job: their
    # durable job was lost, so the queue alone will never touch them again.
    orphans = conn.execute(
        "SELECT t.id, t.genome_id FROM tracks t "
        "WHERE t.track_status = 'converting' "
        "AND NOT EXISTS ("
        "    SELECT 1 FROM conversion_jobs j "
        "    WHERE j.track_id = t.id AND j.state IN ('queued', 'running'))"
    ).fetchall()

    for track in orphans:
        track_id = track['id']
        # The most recent job for this track tells us the kind and where its
        # scratch input was written; without it we cannot rebuild a job.
        last_job = conn.execute(
            'SELECT kind, input_path FROM conversion_jobs '
            'WHERE track_id = ? ORDER BY created_at DESC, id DESC LIMIT 1',
            (track_id,)).fetchone()

        if last_job is not None and last_job['input_path'] \
                and Path(last_job['input_path']).exists():
            # The uploaded original is still in scratch: re-enqueue so the loop
            # converts it (Requirement 13.18, resume/re-enqueue branch).
            conversion_queue.enqueue(
                conn, track_id, last_job['kind'], last_job['input_path'])
            re_enqueued += 1
        else:
            # No job record, or the scratch input is gone: the conversion cannot
            # be resumed, so fail the track with a descriptive reason rather than
            # leave it converting forever (Requirement 13.18, cannot-resume
            # branch).
            conn.execute(
                "UPDATE tracks SET track_status = 'failed', "
                "error = ?, file_path = NULL WHERE id = ?",
                ('Conversion could not be resumed after a restart: the '
                 'uploaded input is no longer available in the scratch area.',
                 track_id))
            failed += 1

    conn.commit()
    return {'running_requeued': running_requeued,
            're_enqueued': re_enqueued, 'failed': failed}


# ---------------------------------------------------------------------------
# One poll-claim-process step (Requirements 13.12, 13.13, 13.16, 13.17)
# ---------------------------------------------------------------------------


def run_once(conn, data_root, config: WorkerConfig) -> bool:
    """
    Do one poll-claim-process step; return whether a job was processed.

    The step, reading only the durable queue (no session/in-memory state, so it
    runs independently of the uploader's browser, Requirement 13.17):

    1. **Free-disk floor gate.** If ``conversion.has_free_disk_floor`` reports
       free disk in the scratch area is below ``CONVERT_MIN_FREE_GB``, do *not*
       claim — leave jobs queued and return ``False`` so the caller waits. The
       worker never starts heavy work it may not be able to finish (Requirement
       13.13).
    2. **Claim the oldest queued job.** ``conversion_queue.claim_next_job`` flips
       exactly one ``queued`` row to ``running`` in a serialized transaction;
       with one worker this is the concurrency-1 guarantee (Requirements 13.12,
       13.16). ``None`` means the queue is empty — return ``False``.
    3. **Convert and close the job.** Look up the track's ``genome_id`` (the job
       row carries only ``track_id``/``kind``/``input_path``) and call
       ``conversion.convert_job``, which runs the toolchain, drives the track row
       ``ready``/``failed``, and cleans scratch. Record the job's terminal state
       from the returned :class:`conversion.ConversionOutcome` with
       ``mark_done``/``mark_failed``. Return ``True``.

    Returns ``True`` when a job was processed, ``False`` when the step did
    nothing (floor not met, or queue empty) so the caller knows to sleep.
    """
    data_root = Path(data_root)

    # 1. Honour the free-disk floor before claiming anything.
    if not conversion.has_free_disk_floor(config.scratch_root,
                                          config.min_free_gb):
        return False

    # 2. Claim the oldest queued job (serialized; concurrency 1).
    job = conversion_queue.claim_next_job(conn)
    if job is None:
        return False

    job_id = job['id']
    track_id = job['track_id']

    # 3. Resolve the genome for the track (the job row does not carry it), then
    #    convert and close the job out from the outcome.
    track_row = conn.execute(
        'SELECT genome_id FROM tracks WHERE id = ?', (track_id,)).fetchone()
    if track_row is None:
        # The track was deleted out from under the job (cascade would normally
        # remove the job too; this is the belt-and-braces case). Close the job
        # as failed rather than call convert_job with no genome.
        conversion_queue.mark_failed(
            conn, job_id,
            'Track no longer exists; conversion cannot proceed.')
        return True

    outcome = conversion.convert_job(
        track_id=track_id,
        genome_id=track_row['genome_id'],
        kind=job['kind'],
        input_path=job['input_path'],
        data_root=data_root,
        samtools_threads=config.samtools_threads,
    )

    if outcome.succeeded:
        conversion_queue.mark_done(conn, job_id)
    else:
        conversion_queue.mark_failed(
            conn, job_id,
            outcome.error or 'Conversion failed for an unknown reason.')
    return True


# ---------------------------------------------------------------------------
# The forever loop (the container entry point)
# ---------------------------------------------------------------------------


def main() -> None:
    """
    The ``genome-worker`` container entry point: recover, then loop forever.

    Reads configuration from the environment, opens the module's own WAL
    database, runs :func:`recover` once to repair anything a crash stranded, then
    loops :func:`run_once`. When a step does nothing (queue empty, or the
    free-disk floor is not met) it sleeps a short poll interval before trying
    again; when it processes a job it loops straight on to the next without
    sleeping, so a backlog drains as fast as one job at a time allows.
    """
    config = load_config()
    # The web app normally creates the schema and data directories on
    # startup, but the worker is a separate container that may come up first
    # on a fresh genome-data volume. Create them here too (idempotent) so the
    # worker never queries a database before its tables exist.
    for sub in ('', 'genomes', 'tracks', 'scratch'):
        (config.data_root / sub).mkdir(parents=True, exist_ok=True)
    db.init_db(config.data_root / 'genome.db')
    conn = db.connect(config.data_root / 'genome.db')
    try:
        summary = recover(conn, config.data_root)
        print(
            f'[genome-worker] startup recovery: '
            f"{summary['running_requeued']} running job(s) re-queued, "
            f"{summary['re_enqueued']} orphaned track(s) re-enqueued, "
            f"{summary['failed']} track(s) failed (input gone).",
            flush=True)

        while True:
            processed = run_once(conn, config.data_root, config)
            if not processed:
                time.sleep(POLL_INTERVAL_SECONDS)
    finally:
        conn.close()


if __name__ == '__main__':
    main()
