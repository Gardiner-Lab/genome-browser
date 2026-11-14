"""
What a coverage track's values are, so that each analysis measures it the right way.

Most coverage tracks (BigWig) hold read coverage: a value on every base read, and nothing where
nothing was read, which means zero there. A methylation track is different: one value per cytosine,
the fraction (or percent) of reads methylated there, and nothing between cytosines, which means
"not measured", not zero. Many published ones (GEO especially) also keep both strands in one file
and give minus-strand cytosines a minus sign, so averaging them cancels the strands against each
other: a CG pair methylated at 0.9 reads +0.9 and -0.9, and averages 0.

So for a track whose ``signal`` is ``methylation``, every analysis measures its *level*:

  * values are taken without their sign (``level_file`` writes a copy with the sign dropped, for
    deepTools, which reads files);
  * stretches without data are left out of a mean, never counted as zero.

``signal`` is found from the file the first time it is needed (``detect``) and kept in the
``tracks`` row; the track's own lab can change it in the track's Details.
"""
import os
from pathlib import Path
from typing import Optional

COVERAGE = 'coverage'
METHYLATION = 'methylation'
KINDS = (COVERAGE, METHYLATION)
LABELS = {COVERAGE: 'Coverage (reads)', METHYLATION: 'Methylation (a level per cytosine)'}

# How much of a track to look at: the start of its longest chromosome.
_SAMPLE_BASES = 2_000_000


def detect(path: str) -> str:
    """
    ``methylation`` when the file looks like per-cytosine calls, else ``coverage``.

    Per-cytosine calls are entries one base long, standing alone or in short chains (a CG pair is
    two, side by side), with values within 0-1 or 0-100, signed or not. Read coverage is long runs
    of joined entries, or single bases joined into long chains, with values that run into the
    hundreds or thousands. On uploaded tracks: CG calls from GEO were 100% single bases in chains of
    1.6, within -1 to 1; RNA-seq coverage, chains of 40,000 entries and values up to 75,000.
    """
    import pyBigWig
    bw = pyBigWig.open(str(path))
    try:
        chroms = bw.chroms() or {}
        if not chroms:
            return COVERAGE
        h = bw.header()
        if max(abs(h.get('minVal') or 0), abs(h.get('maxVal') or 0)) > 100:
            return COVERAGE
        entries = []
        for chrom, size in sorted(chroms.items(), key=lambda kv: -kv[1]):
            entries += bw.intervals(chrom, 0, min(size, _SAMPLE_BASES)) or []
            if len(entries) >= 2000:
                break
    finally:
        bw.close()
    if len(entries) < 100:
        return COVERAGE
    single = sum(1 for s, e, _ in entries if e - s == 1) / len(entries)
    breaks = 1 + sum(1 for a, b in zip(entries, entries[1:]) if b[0] != a[1])
    chain = len(entries) / breaks
    return METHYLATION if single >= 0.8 and chain <= 4 else COVERAGE


def ensure(conn, row) -> str:
    """The track's ``signal``, found from its file and kept the first time it is asked for."""
    keys = row.keys() if hasattr(row, 'keys') else row
    kind = row['signal'] if 'signal' in keys else None
    if kind in KINDS:
        return kind
    if row['type'] != 'bigwig' or not row['file_path'] or not Path(row['file_path']).exists():
        return COVERAGE
    try:
        kind = detect(row['file_path'])
    except Exception:       # an unreadable file is the analysis's to report, not this
        return COVERAGE
    conn.execute("UPDATE tracks SET signal = ?, signal_set_by = 'detected' WHERE id = ? AND signal IS NULL",
                 (kind, row['id']))
    conn.commit()
    return kind


def check_all(conn, say=print) -> dict:
    """
    Every ready coverage track not yet looked at: found coverage or methylation, and for a methylation
    track that gives minus-strand cytosines a minus sign, its copy without signs made now rather than
    at its first analysis. The worker runs it on starting, so tracks added before an update are
    checked then. Returns the counts; ``say`` is told of each methylation track and any trouble.
    """
    counts = {COVERAGE: 0, METHYLATION: 0, 'failed': 0}
    rows = conn.execute("SELECT * FROM tracks WHERE type = 'bigwig' AND track_status = 'ready' "
                        'AND signal IS NULL').fetchall()
    for row in rows:
        try:
            if not row['file_path'] or not Path(row['file_path']).exists():
                raise FileNotFoundError(row['file_path'])
            kind = detect(row['file_path'])
            conn.execute("UPDATE tracks SET signal = ?, signal_set_by = 'detected' WHERE id = ? AND signal IS NULL",
                         (kind, row['id']))
            conn.commit()
            counts[kind] += 1
            if kind == METHYLATION:
                signed = is_signed(row['file_path'])
                if signed:
                    level_file(row['file_path'])
                say(f"methylation track: {row['name']}" + (' (strand signs: copy without them made)' if signed else ''))
        except Exception as exc:
            counts['failed'] += 1
            say(f"could not check the track {row['name']}: {exc}")
    return counts


def is_signed(path: str) -> bool:
    """Whether any value is below zero (minus-strand cytosines marked with a minus sign)."""
    # Not the header's minVal: pyBigWig gives it as a whole number, so -0.8 reads 0. Each
    # chromosome's minimum, from the file's summaries, is quick and keeps the fraction.
    import pyBigWig
    bw = pyBigWig.open(str(path))
    try:
        for chrom in bw.chroms() or {}:
            low = (bw.stats(chrom, type='min') or [None])[0]
            if low is not None and low < 0:
                return True
        return False
    finally:
        bw.close()


def level_path(path: str) -> Path:
    p = Path(path)
    return p.with_name(p.stem + '.level.bw')


def level_file(path: str) -> str:
    """
    The file to measure a methylation track's level from: the track itself when no value is below
    zero, else a copy beside it with every value's sign dropped, made once and kept.
    """
    if not is_signed(path):
        return str(path)
    out = level_path(path)
    if out.exists() and out.stat().st_mtime >= Path(path).stat().st_mtime:
        return str(out)
    import bigwig_means
    tmp = out.with_name(out.name + f'.{os.getpid()}.tmp')
    bigwig_means.unsigned(path, str(tmp))
    os.replace(tmp, out)
    return str(out)
