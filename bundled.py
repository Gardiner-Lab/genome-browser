"""
The genome the browser ships with: Arabidopsis thaliana, TAIR10 with the Araport11 annotation.

The image carries the two files compressed (BUNDLED_DIR, /opt/bundled; see the Dockerfile):

    TAIR10_chr_all.fasta.gz                        the TAIR10 chromosomes (Chr1..Chr5, ChrC, ChrM)
    Araport11_GFF3_genes_transposons.20250813.gff.gz   Araport11 genes and transposons

On its first start on a genome-data volume, the worker (worker.py) adds them as the suite's
reference genome, seen by every lab (and by the teaching demo), with Araport11 as a suite track
on it: exactly what an administrator gets by uploading them with "suite reference". After that it
does nothing: an installation that has a suite TAIR10 already (uploaded by hand, as before this
shipped) keeps its own, and one whose administrator deleted it is not given it back
(the 'bundled.araport11' mark in the settings table).

Set GENOME_BUNDLED_ARAPORT11=0 to not add it. TAIR's data is CC BY 4.0: see NOTICE.
"""

import gzip
import os
import shutil
import sqlite3
import uuid
from pathlib import Path
from typing import Optional

import access
import feature_index
import readers
import regions
import serving

FASTA = 'TAIR10_chr_all.fasta.gz'
GFF = 'Araport11_GFF3_genes_transposons.20250813.gff.gz'
GENOME_NAME = 'Arabidopsis thaliana (TAIR10 / Araport11)'
ASSEMBLY = 'TAIR10'
TRACK_NAME = 'Araport11 genes and transposons (2025-08-13)'
MARK = 'bundled.araport11'
BY = 'system'


def bundled_dir() -> Path:
    return Path(os.environ.get('BUNDLED_DIR', '/opt/bundled'))


def _wanted() -> bool:
    return os.environ.get('GENOME_BUNDLED_ARAPORT11', '1').strip().lower() not in ('0', 'false', 'no', 'off')


def _marked(conn) -> bool:
    conn.execute('CREATE TABLE IF NOT EXISTS settings (key TEXT PRIMARY KEY, value TEXT)')
    return conn.execute('SELECT 1 FROM settings WHERE key = ?', (MARK,)).fetchone() is not None


def _mark(conn, value: str) -> None:
    conn.execute('INSERT OR REPLACE INTO settings (key, value) VALUES (?, ?)', (MARK, value))
    conn.commit()


def _gunzip(source: Path, dest: Path) -> None:
    with gzip.open(source, 'rb') as packed, open(dest, 'wb') as out:
        shutil.copyfileobj(packed, out, 1 << 20)


def seed(conn: sqlite3.Connection, data_root: Path, source: Optional[Path] = None) -> Optional[str]:
    """
    Add the bundled genome and its track, once. Returns the new genome's id, or None when it was
    not added (not wanted, not in the image, added before, or a suite TAIR10 is there already).
    """
    source = source or bundled_dir()
    if not _wanted() or _marked(conn):
        return None
    if not (source / FASTA).exists() or not (source / GFF).exists():
        return None
    existing = conn.execute('SELECT id FROM reference_genomes WHERE owner_lab = ? AND assembly_id = ?',
                            (access.SUITE_OWNER, ASSEMBLY)).fetchone()
    if existing is not None:
        _mark(conn, 'present')         # uploaded by hand already: keep that one
        return None

    genome_id = str(uuid.uuid4())
    genome_dir = data_root / 'genomes' / genome_id
    track_id = uuid.uuid4().hex
    track_dir = data_root / 'tracks' / track_id
    try:
        genome_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = genome_dir / 'reference.fasta'
        _gunzip(source / FASTA, fasta_path)
        result = readers.parse_fasta(fasta_path.read_text(encoding='utf-8', errors='replace'))
        if not result.valid:
            raise ValueError(result.reason)

        track_dir.mkdir(parents=True, exist_ok=True)
        gff_path = track_dir / GFF[:-len('.gz')]
        _gunzip(source / GFF, gff_path)

        conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                     'VALUES (?, ?, ?, ?, ?, ?)',
                     (genome_id, GENOME_NAME, ASSEMBLY, access.SUITE_OWNER, access.SUITE_VISIBILITY, BY))
        conn.executemany('INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
                         [(genome_id, chrom, length) for chrom, length in result.contigs])
        conn.execute('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, '
                     'file_path, created_by) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
                     (track_id, genome_id, TRACK_NAME, 'gff', access.SUITE_OWNER, access.SUITE_VISIBILITY,
                      'ready', str(gff_path), BY))
        conn.execute('CREATE TABLE IF NOT EXISTS settings (key TEXT PRIMARY KEY, value TEXT)')
        conn.execute('INSERT OR REPLACE INTO settings (key, value) VALUES (?, ?)', (MARK, genome_id))
        conn.commit()
    except Exception:
        conn.rollback()
        shutil.rmtree(genome_dir, ignore_errors=True)
        shutil.rmtree(track_dir, ignore_errors=True)
        raise

    # Indexes for the viewer, as an upload makes them; if one fails, the first view makes it.
    serving.prepare_reference(fasta_path)
    try:
        feature_index.build(gff_path, regions.feature_parser('gff'))
    except (OSError, sqlite3.Error):
        pass
    serving.prepare_feature_track(gff_path, 'gff')
    return genome_id
