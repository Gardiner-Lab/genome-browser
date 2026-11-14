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

# Finding a feature by its name (the search box). Built with the index; an index built before it had
# one gets it on the first search (a few seconds for a genome annotation), not a rebuild.
_BY_NAME = 'CREATE INDEX IF NOT EXISTS features_by_name ON features (name COLLATE NOCASE)'


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
        conn.execute(_BY_NAME)
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


# A GFF feature's other names, as the search box may be given them: Araport11's symbol= (FLC) and
# full_name= (FLOWERING LOCUS C), GFF3's Alias=, and gene_name= from Ensembl-style files.
_ALIAS_KEYS = ('symbol', 'Alias', 'gene_name', 'full_name')


def _aliases(file_path) -> Iterable[Feature]:
    """
    Each top-level GFF feature's other names, as (alias, chrom, start, end, name, type, strand): a gene
    or a transposon, not its transcripts and exons, which share them. 0-based starts, as the index.
    """
    from urllib.parse import unquote
    with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) < 9 or 'Parent=' in f[8]:
                continue
            attrs = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            name = unquote(attrs.get('Name') or attrs.get('ID') or '').strip()
            found = set()
            for key in _ALIAS_KEYS:
                value = attrs.get(key)
                if not value:
                    continue
                for part in ([value] if key == 'full_name' else value.split(',')):
                    alias = unquote(part).strip()
                    if alias and alias.casefold() != name.casefold():
                        found.add(alias)
            if not found:
                continue
            try:
                start, end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            strand = 1 if f[6] == '+' else -1 if f[6] == '-' else 0
            for alias in found:
                yield (alias, f[0], start, end, name, f[2], strand)


def prepare_search(file_path, parse, track_type: str) -> None:
    """
    Make a feature file's index ready for the search box: a name index, and for a GFF file a table of
    its features' other names (aliases). Built once, on the first search or when the worker starts
    (seconds for Araport11); an index rebuilt from a newer file gets them again.
    """
    idx = ensure(file_path, parse)
    conn = sqlite3.connect(str(idx), timeout=60)
    try:
        conn.execute(_BY_NAME)
        conn.commit()
        if track_type != 'gff':
            return
        conn.execute('BEGIN IMMEDIATE')       # one builder at a time; the others find it made
        if conn.execute("SELECT 1 FROM sqlite_master WHERE name = 'aliases'").fetchone():
            conn.rollback()
            return
        conn.execute('CREATE TABLE aliases (alias TEXT NOT NULL, chrom TEXT NOT NULL, start INTEGER NOT NULL, '
                     '"end" INTEGER NOT NULL, name TEXT NOT NULL, type TEXT NOT NULL, strand INTEGER NOT NULL)')
        batch = []
        for row in _aliases(file_path):
            batch.append(row)
            if len(batch) >= _BATCH:
                conn.executemany('INSERT INTO aliases VALUES (?, ?, ?, ?, ?, ?, ?)', batch)
                batch.clear()
        if batch:
            conn.executemany('INSERT INTO aliases VALUES (?, ?, ?, ?, ?, ?, ?)', batch)
        conn.execute('CREATE INDEX aliases_by_alias ON aliases (alias COLLATE NOCASE)')
        conn.commit()
    finally:
        conn.close()


def named(file_path, parse, name: str, limit: int = 20, track_type: str = 'gff') -> List[Feature]:
    """
    The features called ``name``, ignoring case: by their own name (a gene's AGI code; its transcripts
    and exons are named apart, ``AT3G01835.1``), or, in a GFF file, by a symbol, alias or full name
    (``FLC``, ``FLOWERING LOCUS C``). Genes first, then the longest. An unreadable file finds nothing.
    """
    try:
        try:
            prepare_search(file_path, parse, track_type)
        except sqlite3.Error:
            pass                          # read-only, or busy: the search scans the names instead
        conn = sqlite3.connect(f'file:{index_path(file_path)}?mode=ro', uri=True)
        try:
            rows = conn.execute(
                'SELECT chrom, start, "end", name, type, strand FROM features WHERE name = ? COLLATE NOCASE '
                "ORDER BY type <> 'gene', \"end\" - start DESC LIMIT ?", (name, limit)).fetchall()
            if conn.execute("SELECT 1 FROM sqlite_master WHERE name = 'aliases'").fetchone():
                rows += conn.execute(
                    'SELECT chrom, start, "end", name, type, strand FROM aliases WHERE alias = ? COLLATE NOCASE '
                    "ORDER BY type <> 'gene', \"end\" - start DESC LIMIT ?", (name, limit)).fetchall()
        finally:
            conn.close()
    except (OSError, sqlite3.Error):
        return []
    out = list(dict.fromkeys(tuple(r) for r in rows))
    out.sort(key=lambda r: (r[4] != 'gene', -(r[2] - r[1])))
    return out[:limit]


def details(file_path, parse, chrom: str, start: int, end: int,
            limit: int = 5000) -> Optional[List[Feature]]:
    """Like :func:`overlapping`, with each feature's type and strand, as tuples."""
    try:
        return [tuple(r) for r in _query(file_path, parse, chrom, start, end,
                                         limit)]
    except (OSError, sqlite3.Error):
        return None
