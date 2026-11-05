"""
The windowed region endpoint: ``GET /api/tracks/<id>/region``.

This is the data feed behind the interactive browser. Given a track and a
visible window, it returns just the data for that window — summarized coverage
bins for a BigWig track (Requirement 8.6) or the features overlapping the range
for a feature track (Requirement 8.7) — so the client (igv.js and the SVG
export path) draws only what is on screen rather than downloading a whole track.

It is a thin HTTP layer over two pure pieces it does not own:

- :mod:`locus` does the coordinate math. The requested ``chrom``/``start``/
  ``end`` are turned into a :class:`locus.Locus` and clamped to the chromosome
  bounds (Requirement 8.3); an unknown chromosome is rejected naming it
  (Requirement 8.2), exactly as a locus string typed into the browser would be.

- :mod:`readers` already knows the file formats; the BED/GFF overlap scan reuses
  its interval readers. The one heavy read, a windowed BigWig summary, goes
  through pyBigWig, imported *lazily* inside :func:`_bigwig_bins` so this module
  (and the tests that stub the bin reader) import without the native library
  present.

Access control is never reinvented here: the track must exist and be viewable
by the signed-in member (:func:`access.can_view_track`, which also requires
``track_status='ready'``), or the endpoint answers 404 — the same "learn
nothing about tracks out of reach" guard the read endpoints use. Identity is
proxy-trusted via :func:`identity.require_member`, and the module relies on the
shared ``identity.Unauthenticated`` handler registered in ``genomes.py`` for the
401 case rather than re-registering one.

Wired by :func:`register_regions`, mirroring ``register_genomes`` /
``register_tracks`` / ``register_figures`` so the app factory stays a thin
wiring layer and this module owns its own route.
"""

from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from flask import Flask, current_app, jsonify, request

import access
import db
import identity
import locus as locus_math
import readers

# The default and maximum number of coverage bins a BigWig region is summarized
# into. The client asks for roughly one bin per screen pixel; the ceiling keeps
# a pathological ``bins=`` from asking pyBigWig for a huge summary.
DEFAULT_BINS = 700
MAX_BINS = 5000


def register_regions(app: Flask) -> None:
    """
    Wire ``GET /api/tracks/<id>/region`` onto the application.

    Called once from ``create_app`` with a single line, mirroring the other
    ``register_*`` registrars. Keeping the route here lets this module own the
    region feed without touching the factory body or the sibling endpoint
    modules.
    """

    @app.route('/api/tracks/<track_id>/region')
    def track_region(track_id):
        """
        Return the data for a track over a window: coverage bins or features.

        Requires a signed-in member and a viewable, ready track (else 404). The
        query carries ``chrom`` (required) and optional ``start``/``end``
        (defaulting to the whole chromosome) and ``bins`` (BigWig only). The
        window is clamped to the chromosome via :mod:`locus`; an unknown
        chromosome is rejected naming it (Requirements 8.2, 8.3). A BigWig track
        returns summarized coverage bins over the window (Requirement 8.6); a
        feature track returns the features overlapping it (Requirement 8.7).
        """
        member = identity.require_member()

        chrom = (request.args.get('chrom') or '').strip()
        if not chrom:
            return jsonify({'error': 'A chrom query parameter is required.'}), 400

        conn = db.connect(_data_root() / 'genome.db')
        try:
            track = _load_track(conn, track_id)
            if track is None or not access.can_view_track(member, track, conn):
                # Same 404 whether the track is absent, not viewable, or not yet
                # ready, so a member learns nothing about tracks out of reach.
                return jsonify({'error': 'Track not found'}), 404

            chrom_sizes = _chrom_sizes(conn, track['genome_id'])
            try:
                window = _requested_window(chrom, request.args, chrom_sizes)
            except locus_math.UnknownChrom as exc:
                return jsonify({'error': str(exc), 'chrom': exc.chrom}), 400
            except locus_math.LocusParseError as exc:
                return jsonify({'error': str(exc)}), 400

            bins = _requested_bins(request.args)
            return _region_body(track, window, bins)
        finally:
            conn.close()


# ---------------------------------------------------------------------------
# Window + bin request parsing
# ---------------------------------------------------------------------------


def _requested_window(chrom: str, args, chrom_sizes: Dict[str, int]) -> locus_math.Locus:
    """
    Build the clamped :class:`locus.Locus` for the requested region (8.2, 8.3).

    ``chrom`` must be in ``chrom_sizes`` or :class:`locus.UnknownChrom` is
    raised (naming it). ``start``/``end`` are optional 0-based half-open
    coordinates; absent, they default to the whole chromosome. The window is
    clamped to the chromosome bounds, so an out-of-range request snaps to the
    edge rather than erroring (Requirement 8.3).
    """
    if chrom not in chrom_sizes:
        raise locus_math.UnknownChrom(chrom)
    length = int(chrom_sizes[chrom])

    start = _int_arg(args.get('start'), default=0)
    end = _int_arg(args.get('end'), default=length)
    return locus_math.clamp(
        locus_math.Locus(chrom=chrom, start=start, end=end), chrom_sizes)


def _requested_bins(args) -> int:
    """The bin count for a BigWig summary, bounded to ``1..MAX_BINS``."""
    bins = _int_arg(args.get('bins'), default=DEFAULT_BINS)
    return max(1, min(MAX_BINS, bins))


def _int_arg(value, default: int) -> int:
    """Parse an integer query arg, falling back to ``default`` on bad input."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Building the response per track kind
# ---------------------------------------------------------------------------


def _region_body(track, window: locus_math.Locus, bins: int):
    """
    Assemble the JSON region response for the track's kind.

    A BigWig track (``type='bigwig'``) carries summarized coverage bins
    (Requirement 8.6); a BED/GFF feature track carries the overlapping features
    (Requirement 8.7). The window (already clamped) is echoed back so the client
    knows the actual range it received after clamping.
    """
    body = {
        'track_id': track['id'],
        'type': track['type'],
        'chrom': window.chrom,
        'start': window.start,
        'end': window.end,
    }
    file_path = track['file_path'] if _has(track, 'file_path') else None

    if track['type'] == 'bigwig':
        values = _bigwig_bins(file_path, window, bins)
        body['bins'] = bins
        body['values'] = values
        return jsonify(body)

    features = _feature_overlaps(track['type'], file_path, window)
    body['features'] = features
    return jsonify(body)


def _bigwig_bins(file_path: Optional[str], window: locus_math.Locus,
                 bins: int,
                 reader: Optional[Callable[..., List[Optional[float]]]] = None
                 ) -> List[Optional[float]]:
    """
    Summarize a BigWig into ``bins`` mean-coverage values over the window (8.6).

    pyBigWig's ``stats(chrom, start, end, type='mean', nBins=bins)`` returns one
    summary value per bin across ``[start, end)`` — exactly the per-pixel
    coverage the browser plots. pyBigWig is imported lazily so this module loads
    without the native library, and the actual read is isolated behind
    ``reader`` so a test can inject a stub returning a small values list instead
    of opening a real BigWig. A missing file or an empty window yields an empty
    list; a bin with no data comes back as ``None`` (JSON ``null``).
    """
    if not file_path or window.span <= 0:
        return []
    read = reader if reader is not None else _pybigwig_stats
    try:
        return read(file_path, window.chrom, window.start, window.end, bins)
    except Exception as exc:  # pragma: no cover - defensive around the native lib
        current_app.logger.warning('BigWig region read failed: %s', exc)
        return []


def _pybigwig_stats(path: str, chrom: str, start: int, end: int,
                    bins: int) -> List[Optional[float]]:
    """
    The default BigWig bin reader: pyBigWig ``stats`` (imported lazily).

    Kept private and injectable so :func:`_bigwig_bins` can be tested without
    pyBigWig or a real file. Returns one mean value per bin over ``[start, end)``,
    with ``None`` for bins the file has no data for.
    """
    import pyBigWig  # imported lazily; native lib only needed for a real read

    bw = None
    try:
        bw = pyBigWig.open(path)
        values = bw.stats(chrom, start, end, type='mean', nBins=bins)
    finally:
        if bw is not None:
            bw.close()
    return [None if v is None else float(v) for v in (values or [])]


def _feature_overlaps(track_type: str, file_path: Optional[str],
                      window: locus_math.Locus) -> List[Dict[str, object]]:
    """
    The features of a BED/GFF track overlapping the window (Requirement 8.7).

    The uploaded feature file is read with the matching :mod:`readers` scanner
    and every interval on the window's chromosome that overlaps ``[start, end)``
    is returned as ``{chrom, start, end, name}``, sorted by start. Overlap is
    half-open: a feature overlaps when ``feature.start < window.end`` and
    ``feature.end > window.start``. A missing or unreadable file yields an empty
    list, so the lane renders empty rather than erroring.
    """
    if not file_path:
        return []
    path = Path(file_path)
    if not path.exists():
        return []
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return []

    intervals = _scan_features(track_type, text)
    out: List[Dict[str, object]] = []
    for chrom, start, end, name in intervals:
        if chrom != window.chrom:
            continue
        if start < window.end and end > window.start:
            out.append({'chrom': chrom, 'start': start, 'end': end,
                        'name': name})
    out.sort(key=lambda f: (f['start'], f['end']))
    return out


def _scan_features(track_type: str, text: str
                   ) -> List[Tuple[str, int, int, str]]:
    """
    Parse a feature file into ``(chrom, start, end, name)`` tuples.

    BED is 0-based half-open, GFF/GFF3 is 1-based inclusive; both are normalised
    here to 0-based half-open coordinates so overlap testing against the window
    is uniform. Lines that do not parse into coordinates are skipped. This reads
    the body directly rather than through ``readers``' chromosome-extent scan,
    because the region feed needs each individual interval, not just the furthest
    coordinate per chromosome.
    """
    if track_type == 'bed':
        return _scan_bed_features(text)
    if track_type == 'gff':
        return _scan_gff_features(text)
    return []


def _scan_bed_features(text: str) -> List[Tuple[str, int, int, str]]:
    """Parse BED ``chrom start end [name ...]`` (0-based half-open) rows."""
    features: List[Tuple[str, int, int, str]] = []
    for raw in text.splitlines():
        line = raw.rstrip('\n')
        if _is_comment(line):
            continue
        fields = line.split('\t') if '\t' in line else line.split()
        if len(fields) < 3:
            continue
        start = _int_or_none(fields[1])
        end = _int_or_none(fields[2])
        if start is None or end is None:
            continue
        name = fields[3] if len(fields) >= 4 else ''
        features.append((fields[0], start, end, name))
    return features


def _scan_gff_features(text: str) -> List[Tuple[str, int, int, str]]:
    """
    Parse GFF/GFF3 features (1-based inclusive) into 0-based half-open tuples.

    Columns are ``seqid source type start end score strand phase attributes``;
    the 1-based inclusive ``start``/``end`` become a 0-based half-open window
    (``start-1`` .. ``end``). A ``Name=``/``ID=`` attribute, when present, is
    used as the feature name.
    """
    features: List[Tuple[str, int, int, str]] = []
    for raw in text.splitlines():
        line = raw.rstrip('\n')
        if _is_comment(line):
            continue
        fields = line.split('\t') if '\t' in line else line.split()
        if len(fields) < 5:
            continue
        start = _int_or_none(fields[3])
        end = _int_or_none(fields[4])
        if start is None or end is None:
            continue
        name = _gff_name(fields[8]) if len(fields) >= 9 else ''
        features.append((fields[0], start - 1, end, name))
    return features


def _gff_name(attributes: str) -> str:
    """Pull a human name from a GFF attributes column (``Name=`` then ``ID=``)."""
    parts = {}
    for token in attributes.replace(' ', '').split(';'):
        if '=' in token:
            key, _, value = token.partition('=')
            parts[key.lower()] = value
    return parts.get('name') or parts.get('id') or ''


def _is_comment(line: str) -> bool:
    stripped = line.strip()
    return (stripped == ''
            or stripped.startswith('#')
            or stripped.startswith('track ')
            or stripped.startswith('browser '))


def _int_or_none(value: str) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Database reads
# ---------------------------------------------------------------------------


def _load_track(conn, track_id: str):
    """The tracks row for this id, or None when there is none."""
    return conn.execute(
        'SELECT id, genome_id, name, type, owner_lab, visibility, '
        'track_status, file_path FROM tracks WHERE id = ?',
        (track_id,)).fetchone()


def _chrom_sizes(conn, genome_id: str) -> Dict[str, int]:
    """The genome's ``chrom -> length_bp`` map, for clamping the window."""
    rows = conn.execute(
        'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?',
        (genome_id,)).fetchall()
    return {row['chrom']: int(row['length_bp']) for row in rows}


def _data_root() -> Path:
    """The configured data root (the genome-data volume)."""
    return Path(current_app.config['DATA_ROOT'])


def _has(row, key: str) -> bool:
    try:
        return key in row.keys()
    except AttributeError:
        return False
