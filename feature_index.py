"""
An index beside each stored BED/GFF track, so a region read is a lookup.

The region feed used to read and parse the whole feature file for every window
a viewer asked for. That is fine for a few hundred features and impossible for
a genome annotation: Araport11's GFF3 is 289 MB and ~930,000 features, read
again on every pan and zoom. Instead, each feature file gets a small SQLite file
beside it (``<file>.idx.sqlite``) holding every feature as ``(chrom, start,
end, name, type, strand)`` in 0-based half-open coordinates, indexed by
``(chrom, start)``.

An overlap query needs features with ``start < window.end`` and ``end >
window.start``. The index serves the first directly; for the second it also
bounds ``start >= window.start - longest`` — the longest feature on that
chromosome, recorded at build time — so the scan stays within the window.

The index is built when a track is stored and, for a track stored before
indexes existed (or before this layout), the first time it is read. It is
written to a temporary name and moved into place, so a reader never sees a
half-built index and two concurrent first reads just build it twice. It is
rebuilt whenever the feature file is newer than it.
"""

import os
import sqlite3
import uuid
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

# (chrom, start, end, name, type, strand); strand is 1, -1 or 0.
Feature = Tuple[str, int, int, str, str, int]

# The layout of the index; an index of another version is rebuilt.
VERSION = 2

_BATCH = 50_000


def index_path(file_path) -> Path:
    return Path(str(file_path) + '.idx.sqlite')


def is_current(file_path) -> bool:
    idx = index_path(file_path)
    try:
        if idx.stat().st_mtime < Path(file_path).stat().st_mtime:
            return False
        conn = sqlite3.connect(f'file:{idx}?mode=ro', uri=True)
        try:
            return conn.execute('PRAGMA user_version').fetchone()[0] == VERSION
        finally:
            conn.close()
    except (OSError, sqlite3.Error):
        return False


def build(file_path, parse: Callable[[Iterable[str]], Iterable[Feature]]) -> Path:
    """
    Build (or rebuild) the index for a feature file and return its path.

    ``parse`` turns the file's lines into features; it is the region feed's own
    BED or GFF parser, so the index holds exactly what a full read would.
    """
    target = index_path(file_path)
    tmp = target.with_name(f'{target.name}.{uuid.uuid4().hex}.tmp')
    conn = sqlite3.connect(str(tmp))
    try:
        conn.execute('PRAGMA journal_mode = OFF')
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('CREATE TABLE features (chrom TEXT NOT NULL, '
                     'start INTEGER NOT NULL, "end" INTEGER NOT NULL, '
                     'name TEXT NOT NULL, type TEXT NOT NULL, '
                     'strand INTEGER NOT NULL)')
        longest: Dict[str, int] = {}
        batch: List[Feature] = []
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            for feature in parse(handle):
                chrom, start, end = feature[0], feature[1], feature[2]
                span = end - start
                if span > longest.get(chrom, -1):
                    longest[chrom] = span
                batch.append(feature)
                if len(batch) >= _BATCH:
                    conn.executemany(
                        'INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)', batch)
                    batch.clear()
        if batch:
            conn.executemany('INSERT INTO features VALUES (?, ?, ?, ?, ?, ?)',
                             batch)
        conn.execute('CREATE INDEX features_by_start ON features (chrom, start)')
        conn.execute('CREATE TABLE longest (chrom TEXT PRIMARY KEY, '
                     'span INTEGER NOT NULL)')
        conn.executemany('INSERT INTO longest VALUES (?, ?)', longest.items())
        conn.execute(f'PRAGMA user_version = {VERSION}')
        conn.commit()
    except BaseException:
        conn.close()
        tmp.unlink(missing_ok=True)
        raise
    conn.close()
    os.replace(tmp, target)
    return target


def ensure(file_path, parse) -> Path:
    """The index for a feature file, building it first when missing or stale."""
    if is_current(file_path):
        return index_path(file_path)
    return build(file_path, parse)


def _query(file_path, parse, chrom: str, start: int, end: int,
           limit: Optional[int] = None):
    idx = ensure(file_path, parse)
    conn = sqlite3.connect(f'file:{idx}?mode=ro', uri=True)
    try:
        row = conn.execute('SELECT span FROM longest WHERE chrom = ?',
                           (chrom,)).fetchone()
        if row is None:
            return []
        sql = ('SELECT chrom, start, "end", name, type, strand FROM features '
               'WHERE chrom = ? AND start < ? AND start >= ? AND "end" > ? '
               'ORDER BY start, "end"')
        args = [chrom, end, start - int(row[0]), start]
        if limit is not None:
            sql += ' LIMIT ?'
            args.append(limit)
        return conn.execute(sql, args).fetchall()
    finally:
        conn.close()


def overlapping(file_path, parse, chrom: str, start: int,
                end: int) -> Optional[List[Dict[str, object]]]:
    """
    The features on ``chrom`` overlapping ``[start, end)``, sorted by start, as
    ``{chrom, start, end, name}`` (the region feed's shape).

    Returns None when the index cannot be built or read (an unreadable file),
    so the caller can fall back to answering with no features.
    """
    try:
        rows = _query(file_path, parse, chrom, start, end)
    except (OSError, sqlite3.Error):
        return None
    return [{'chrom': c, 'start': s, 'end': e, 'name': n}
            for c, s, e, n, _, _ in rows]


def details(file_path, parse, chrom: str, start: int, end: int,
            limit: int = 5000) -> Optional[List[Feature]]:
    """Like :func:`overlapping`, with each feature's type and strand, as tuples."""
    try:
        return [tuple(r) for r in _query(file_path, parse, chrom, start, end,
                                         limit)]
    except (OSError, sqlite3.Error):
        return None
