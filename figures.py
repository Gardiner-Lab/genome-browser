"""
The SVG figure-export endpoint: ``POST /api/figure.svg``.

This is the thin HTTP layer over the pure SVG builder in :mod:`render`. It does
the three things the builder deliberately does not: find out who is asking
(:mod:`identity`), read the member's persisted ``Browser_View`` and the stored
track data from the module's database and the genome-data volume, and apply
access control so a member only ever sees tracks they may view
(:func:`access.can_view_track`). It then assembles a :class:`render.ViewSpec`
of the *enabled* tracks in Browser_View order and returns
``render.build_svg(view)`` as an ``image/svg+xml`` document.

Kept in its own module with a :func:`register_figures` registrar (mirroring
``app.register_health``) so it can be wired with a single line in the app
factory, without this module and the concurrently-developed genome/view
endpoints having to edit the factory body at the same time.

Requirements: 6.4 (disabled tracks excluded), 9.1 (enabled tracks top-to-bottom
in order + the locus), 9.2/9.4 (each track in its Rendering_Options colour),
9.3 (empty-view SVG with the scale, no track), 9.5 (invalid colour rejected,
previous retained, error indicated), 9.6 (label verbatim ≤200 chars), 9.7
(coordinate scale spanning the locus), 9.8 (well-formed document).
"""

import json
import re
from typing import Dict, List, Optional, Tuple

from flask import Response, jsonify, request

import access
import db
import identity
import render


def register_figures(app) -> None:
    """
    Register ``POST /api/figure.svg`` on the app (called by the app factory).

    One registrar, mirroring ``register_health``, so the app factory wires this
    with a single line and this module owns its own route.
    """

    @app.route('/api/figure.svg', methods=['POST'])
    def figure_svg():
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401

        payload = request.get_json(silent=True) or {}
        genome_id = (payload.get('genome_id') or '').strip()
        if not genome_id:
            return jsonify({'error': 'A genome_id is required.'}), 400

        conn = db.connect(app.config['DATA_ROOT'] / 'genome.db')
        try:
            genome = _genome_row(conn, genome_id)
            if genome is None:
                return jsonify({'error': 'Genome not found.'}), 404
            if not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'You may not view this genome.'}), 403

            view = _build_view(conn, member, genome, payload)
        finally:
            conn.close()

        # Report any colour rejections as a header so the SVG stays a pure
        # document; the previous colour is already retained in the view
        # (Requirement 9.5).
        svg = render.build_svg(view.spec)
        resp = Response(svg, mimetype='image/svg+xml')
        resp.headers['Content-Disposition'] = 'attachment; filename="figure.svg"'
        if view.colour_errors:
            resp.headers['X-Colour-Warning'] = ' '.join(view.colour_errors)[:400]
        return resp


# ---------------------------------------------------------------------------
# Assembling a ViewSpec from the persisted Browser_View + stored data
# ---------------------------------------------------------------------------


class _BuiltView:
    """A rendered :class:`render.ViewSpec` plus any colour-rejection messages."""

    def __init__(self, spec: render.ViewSpec, colour_errors: List[str]):
        self.spec = spec
        self.colour_errors = colour_errors


def _build_view(conn, member: str, genome, payload: Dict) -> "_BuiltView":
    """
    Turn the persisted Browser_View (and request overrides) into a ViewSpec.

    The enabled-track order and the per-track Rendering_Options come from the
    member's stored ``browser_views`` row, with the request body allowed to
    override the locus/enabled/options for an export of exactly what the member
    is looking at. Disabled tracks never enter the list (Requirement 6.4); a
    track the member may not view is dropped too. Each track's colour is run
    through :func:`render.apply_colour`, so an invalid colour is rejected and
    the previous colour retained (Requirement 9.5).
    """
    stored = _stored_view(conn, member, genome['id'])

    locus_str = payload.get('locus') or stored.get('locus') or ''
    enabled = payload.get('enabled')
    if enabled is None:
        enabled = stored.get('enabled', [])
    options = payload.get('options')
    if options is None:
        options = stored.get('options', {})
    prev_options = stored.get('options', {})

    chrom_sizes = _chrom_sizes(conn, genome['id'])
    locus = _parse_locus(locus_str, chrom_sizes)

    colour_errors: List[str] = []
    tracks: List[render.TrackSpec] = []
    for track_id in enabled:
        row = _track_row(conn, track_id, genome['id'])
        if row is None or not access.can_view_track(member, row, conn):
            continue  # unknown, wrong genome, not-ready, or not viewable
        opt_in = options.get(track_id, {}) if isinstance(options, dict) else {}
        prev = prev_options.get(track_id, {}) if isinstance(prev_options, dict) else {}
        render_opts, err = _render_options(opt_in, prev)
        if err:
            colour_errors.append(f'{row["name"]}: {err}')
        tracks.append(render.TrackSpec(
            name=row['name'],
            kind='bigwig' if row['type'] == 'bigwig' else 'feature',
            options=render_opts,
            signal=[],      # region data read in task 11.1; empty lane here
            features=[],
        ))

    label = _figure_label(payload.get('label'))
    spec = render.ViewSpec(
        locus=locus,
        tracks=tracks,
        label=label,
        genome_name=genome['name'] if _has(genome, 'name') else '',
    )
    return _BuiltView(spec, colour_errors)


def _render_options(opt_in: Dict, prev: Dict) -> Tuple[render.RenderOptions, Optional[str]]:
    """
    Build :class:`render.RenderOptions`, validating the colour (Requirement 9.5).

    The requested colour (if any) is validated against the previous colour via
    :func:`render.apply_colour`: valid replaces it, invalid is rejected and the
    previous colour kept, with the returned error surfaced to the caller.
    Height and signal range are coerced to numbers, falling back to defaults.
    """
    previous_colour = str(prev.get('colour') or render.DEFAULT_TRACK_COLOUR)
    requested = opt_in.get('colour') if isinstance(opt_in, dict) else None
    colour, err = render.apply_colour(previous_colour, requested)

    height = _as_int(opt_in.get('height'), default=60)
    signal_min = _as_float(opt_in.get('signal_min'))
    signal_max = _as_float(opt_in.get('signal_max'))
    return render.RenderOptions(
        colour=colour, height=height,
        signal_min=signal_min, signal_max=signal_max), err


def _figure_label(raw) -> Optional[render.FigureLabel]:
    """Build a FigureLabel from the request, or None when none was supplied."""
    if not raw:
        return None
    if isinstance(raw, str):
        text = raw
        return render.FigureLabel(text=text) if text else None
    if isinstance(raw, dict):
        text = raw.get('text') or ''
        if not text:
            return None
        return render.FigureLabel(
            text=text,
            x=_as_float(raw.get('x')) if raw.get('x') is not None else 10.0,
            y=_as_float(raw.get('y')) if raw.get('y') is not None else 20.0)
    return None


# ---------------------------------------------------------------------------
# Database reads (all tolerant of a schema still being filled in by siblings)
# ---------------------------------------------------------------------------


def _genome_row(conn, genome_id: str):
    try:
        return conn.execute(
            'SELECT * FROM reference_genomes WHERE id = ?',
            (genome_id,)).fetchone()
    except Exception:
        return None


def _track_row(conn, track_id: str, genome_id: str):
    try:
        return conn.execute(
            'SELECT * FROM tracks WHERE id = ? AND genome_id = ?',
            (track_id, genome_id)).fetchone()
    except Exception:
        return None


def _chrom_sizes(conn, genome_id: str) -> Dict[str, int]:
    try:
        rows = conn.execute(
            'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?',
            (genome_id,)).fetchall()
    except Exception:
        return {}
    return {r['chrom']: int(r['length_bp']) for r in rows}


def _stored_view(conn, member: str, genome_id: str) -> Dict:
    """
    The member's persisted Browser_View for this genome, defaults when absent.

    Returns ``locus`` (string), ``enabled`` (ordered track ids) and ``options``
    (track id -> option dict). A missing row, missing table (sibling task still
    building it), or malformed JSON degrades to empty defaults rather than
    failing the export.
    """
    out = {'locus': '', 'enabled': [], 'options': {}}
    try:
        row = conn.execute(
            'SELECT locus, enabled_json, options_json FROM browser_views '
            'WHERE username = ? AND genome_id = ?',
            (member, genome_id)).fetchone()
    except Exception:
        return out
    if row is None:
        return out
    out['locus'] = row['locus'] or ''
    out['enabled'] = _json_list(row['enabled_json'])
    out['options'] = _json_dict(row['options_json'])
    return out


# ---------------------------------------------------------------------------
# Locus parsing (lenient; the authoritative parse/clamp lives in task 11.1)
# ---------------------------------------------------------------------------

_LOCUS_RE = re.compile(r'^\s*([^:]+?)\s*(?::\s*([\d,]+)\s*-\s*([\d,]+)\s*)?$')


def _parse_locus(locus_str: str, chrom_sizes: Dict[str, int]) -> render.Locus:
    """
    Parse ``Chr1:1,000-5,000`` (or bare ``Chr1``) into a :class:`render.Locus`.

    Coordinates are clamped to ``0..length`` when the chromosome length is
    known, and a bare chromosome spans its whole length (or ``0..0`` when the
    length is unknown). A string that will not parse yields a harmless
    zero-width locus so the builder still produces a well-formed SVG with a
    scale (Requirement 9.3/9.8). This is intentionally lenient; the canonical
    locus math is task 11.1's ``genome/locus.py``.
    """
    match = _LOCUS_RE.match(locus_str or '')
    if not match:
        chrom = (locus_str or '').strip() or 'unknown'
        return render.Locus(chrom=chrom, start=0, end=chrom_sizes.get(chrom, 0))

    chrom = match.group(1).strip()
    length = chrom_sizes.get(chrom)
    if match.group(2) is None:
        end = length if length is not None else 0
        return render.Locus(chrom=chrom, start=0, end=end)

    start = int(match.group(2).replace(',', ''))
    end = int(match.group(3).replace(',', ''))
    if end < start:
        start, end = end, start
    if length is not None:
        start = max(0, min(start, length))
        end = max(0, min(end, length))
    else:
        start = max(0, start)
        end = max(start, end)
    return render.Locus(chrom=chrom, start=start, end=end)


# ---------------------------------------------------------------------------
# Small coercions
# ---------------------------------------------------------------------------


def _has(row, key: str) -> bool:
    try:
        return key in row.keys()
    except Exception:
        return False


def _json_list(text) -> List[str]:
    try:
        value = json.loads(text) if text else []
    except (TypeError, ValueError):
        return []
    return [str(x) for x in value] if isinstance(value, list) else []


def _json_dict(text) -> Dict:
    try:
        value = json.loads(text) if text else {}
    except (TypeError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


def _as_int(value, default: int) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _as_float(value) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
