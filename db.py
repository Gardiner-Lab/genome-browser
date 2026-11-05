"""
The Genome Browser's own SQLite database, and read-only openers for the
shared stores.

The module keeps its data in one SQLite database on the genome-data volume
(``/genome/genome.db``), opened in WAL mode. WAL matters here because two
processes share the file: the web app (gunicorn) accepts uploads and writes
job rows, and the separate conversion worker reads and updates them. WAL lets
readers and a writer work at once without blocking each other, which is what
makes the durable job queue the coordination point between the two containers
rather than a message broker.

The shared stores — ``suite.db`` for labs and chair groups, ``moclo.db`` for
the account directory — belong to other services and are mounted read-only.
They are always opened ``mode=ro`` so this module can never write to them, and
so a crash here cannot corrupt another tool's data. These openers are used by
``identity.py`` and ``labs_reader.py``; everything that writes goes through
``connect()`` against this module's own database.
"""

import os
import sqlite3
from pathlib import Path
from typing import Optional

# The module's own database on the genome-data volume.
GENOME_DB = Path(os.environ.get('GENOME_DATA_ROOT', '/genome')) / 'genome.db'

# The schema lives beside this module; init_db applies it if present. Task 1.2
# creates schema.sql with the full data model; until then init_db is a no-op
# when the file is absent, so the skeleton comes up on an empty volume.
SCHEMA_FILE = Path(__file__).with_name('schema.sql')


def _db_path(path=None) -> Path:
    return Path(path) if path is not None else GENOME_DB


def connect(path=None) -> sqlite3.Connection:
    """
    A read-write connection to the module's own database, in WAL mode.

    Foreign keys are enforced (the schema relies on ``ON DELETE CASCADE``), and
    a generous busy timeout lets the app and the worker wait briefly for each
    other's write lock rather than failing outright.
    """
    db = _db_path(path)
    db.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute('PRAGMA journal_mode = WAL')
    conn.execute('PRAGMA foreign_keys = ON')
    conn.execute('PRAGMA busy_timeout = 30000')
    return conn


def open_readonly(path) -> Optional[sqlite3.Connection]:
    """
    A read-only connection to a shared store, or None when it is not there.

    Used for ``suite.db`` and ``moclo.db``, which belong to other services and
    are mounted read-only; ``mode=ro`` makes a write impossible and keeps this
    module from ever touching another tool's data.
    """
    db = Path(path)
    if not db.exists():
        return None
    try:
        conn = sqlite3.connect(f'file:{db}?mode=ro', uri=True, timeout=5)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error:
        return None


def open_suite() -> Optional[sqlite3.Connection]:
    """Read-only ``suite.db`` (labs and chair groups), or None when absent."""
    return open_readonly(os.environ.get('SUITE_DATABASE_PATH', '/suite/suite.db'))


def open_moclo() -> Optional[sqlite3.Connection]:
    """Read-only ``moclo.db`` (the account directory), or None when absent."""
    return open_readonly(os.environ.get('MOCLO_DB', '/moclo/moclo.db'))


def init_db(path=None) -> None:
    """
    Create the module's tables if they are not there yet.

    Applies ``schema.sql`` (idempotent ``CREATE TABLE IF NOT EXISTS`` ...)
    against the module's own database so a fresh volume comes up ready. Task
    1.2 writes schema.sql with the full data model; while it is absent this
    just opens the database so WAL mode is set and the file exists.
    """
    conn = connect(path)
    try:
        if SCHEMA_FILE.exists():
            conn.executescript(SCHEMA_FILE.read_text())
            conn.commit()
    finally:
        conn.close()


def database_reachable(path=None) -> bool:
    """
    Whether the module's database can be opened and queried, for /health.

    Runs a trivial query so a missing file, a permissions problem, or a locked
    or corrupt database all read as unreachable (health returns 503).
    """
    try:
        conn = connect(path)
        try:
            conn.execute('SELECT 1')
            return True
        finally:
            conn.close()
    except sqlite3.Error:
        return False
