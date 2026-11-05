"""
Server-side SVG figure generation for the Genome Browser.

For publication export the module renders an SVG on the server rather than
leaning on igv.js's client-side export, because the server already holds the
stored BigWig/feature data and the member's ``Rendering_Options`` and so can
produce deterministic, repeatable output with exact ``Figure_Label`` placement,
strict colour validation, and a guaranteed well-formed document (design
"igv.js for interactive view, server-side SVG for export").

This module is split into two halves on purpose:

* **Pure SVG building** — :func:`build_svg` takes a plain :class:`ViewSpec`
  (the visible :class:`Locus`, the ordered list of *enabled* :class:`TrackSpec`
  with their :class:`RenderOptions`, and an optional :class:`FigureLabel`) and
  returns an SVG string. It touches no database, no filesystem, and no HTTP, so
  the property tests in task 12.2 drive it directly with constructed views.
  Every piece of member text is XML-escaped and every attribute is emitted from
  validated values, so the output is always well-formed XML — including the
  empty-view case (Requirements 9.3, 9.8).

* **Colour validation** — :func:`validate_colour` is the one place a colour is
  judged valid. An invalid colour is *rejected*: the caller keeps the track's
  previous colour and surfaces an error, and the SVG never contains an unsafe
  colour string (Requirements 9.4, 9.5).

The requirement mapping:

- Enabled tracks are drawn top-to-bottom in Browser_View order; disabled tracks
  are excluded (the caller passes only enabled tracks, in order) — 6.4, 9.1.
- Each track is drawn with its current Rendering_Options (colour, height,
  signal range) — 9.2, 9.4.
- The empty view still yields an SVG with the locus and its coordinate scale
  and no track drawn — 9.3.
- An invalid colour is rejected and the previous colour retained, with an error
  indication — 9.5.
- A Figure_Label is drawn verbatim at its position, capped at 200 characters —
  9.6.
- A coordinate scale spans the locus start and end — 9.7.
- The document is well-formed and parses under a conforming XML/SVG parser —
  9.8.
"""

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from xml.sax.saxutils import escape, quoteattr


# ---------------------------------------------------------------------------
# Colour validation (Requirements 9.4, 9.5)
# ---------------------------------------------------------------------------
#
# A colour is accepted only if it is one the SVG spec itself accepts, so that a
# validated colour can be written straight into a `fill`/`stroke` attribute
# without any further escaping and the document stays well-formed. We accept:
#   * #RGB / #RRGGBB / #RGBA / #RRGGBBAA hex,
#   * rgb()/rgba() with integer or percentage channels,
#   * the SVG/CSS named colours.
# Anything else (a bare word that is not a known name, a string with angle
# brackets or quotes, an empty value) is rejected: the track keeps its previous
# colour and the caller reports the colour invalid.

_HEX_RE = re.compile(r'^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{4}|'
                     r'[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$')
_RGB_RE = re.compile(
    r'^rgba?\(\s*'
    r'(\d{1,3}%?|\d{1,3}\.\d+%)\s*,\s*'
    r'(\d{1,3}%?|\d{1,3}\.\d+%)\s*,\s*'
    r'(\d{1,3}%?|\d{1,3}\.\d+%)'
    r'(?:\s*,\s*(0|1|0?\.\d+|\d{1,3}%))?'
    r'\s*\)$')

# The SVG 1.1 / CSS named colours. Kept as a frozenset so membership is exact
# and case-insensitive (we lower-case before testing).
_NAMED_COLOURS = frozenset({
    'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige',
    'bisque', 'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown',
    'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral',
    'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan',
    'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey', 'darkkhaki',
    'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred',
    'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray',
    'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue',
    'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'floralwhite',
    'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod',
    'gray', 'grey', 'green', 'greenyellow', 'honeydew', 'hotpink', 'indianred',
    'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen',
    'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan',
    'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey',
    'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue',
    'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow',
    'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine',
    'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen',
    'mediumslateblue', 'mediumspringgreen', 'mediumturquoise',
    'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 'moccasin',
    'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange',
    'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise',
    'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum',
    'powderblue', 'purple', 'rebeccapurple', 'red', 'rosybrown', 'royalblue',
    'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna',
    'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow',
    'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'tomato',
    'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow',
    'yellowgreen', 'transparent',
})

# The colour a track falls back to when none has ever been set.
DEFAULT_TRACK_COLOUR = '#3b5b92'

# The hard cap on a Figure_Label, in characters (Requirement 9.6).
MAX_LABEL_CHARS = 200


def validate_colour(value: str) -> bool:
    """
    Whether ``value`` is a colour the figure may safely draw (Req 9.4, 9.5).

    Accepts hex (``#RGB``/``#RGBA``/``#RRGGBB``/``#RRGGBBAA``), ``rgb()``/
    ``rgba()`` with integer or percentage channels, and the SVG/CSS named
    colours (case-insensitive). Everything else — an unknown word, a value with
    angle brackets or quotes, ``None`` or an empty string — is invalid, so the
    caller rejects it and keeps the track's previous colour (Requirement 9.5).
    Because only values accepted here ever reach an SVG attribute, the document
    cannot be made malformed through a colour.
    """
    if not isinstance(value, str):
        return False
    colour = value.strip()
    if not colour:
        return False
    if _HEX_RE.match(colour):
        return True
    if _RGB_RE.match(colour.replace(' ', '')) or _RGB_RE.match(colour):
        return True
    return colour.lower() in _NAMED_COLOURS


def apply_colour(previous: str, requested: Optional[str]) -> Tuple[str, Optional[str]]:
    """
    Resolve a requested colour against the track's previous colour (Req 9.5).

    Returns ``(colour, error)``. A valid ``requested`` colour becomes the new
    colour with no error. An invalid one is rejected: the ``previous`` colour is
    returned unchanged alongside an error message stating the colour is invalid,
    so a caller applying a colour change never loses the prior value. A
    ``requested`` of ``None`` means "no change" and returns the previous colour.
    """
    prev = previous if (isinstance(previous, str) and previous) else DEFAULT_TRACK_COLOUR
    if requested is None:
        return prev, None
    if validate_colour(requested):
        return requested.strip(), None
    return prev, (f'The colour {requested!r} is not a valid colour; '
                  f'the track keeps its previous colour.')


# ---------------------------------------------------------------------------
# The view spec the pure builder renders (no DB / HTTP / filesystem)
# ---------------------------------------------------------------------------


@dataclass
class Locus:
    """
    The visible region a figure is drawn for.

    ``chrom`` names the chromosome/contig; ``start``/``end`` are its base-pair
    bounds (0-based start, end exclusive, ``start <= end``). The coordinate
    scale spans exactly ``start``..``end`` (Requirement 9.7). A zero-width locus
    (``start == end``) is tolerated so the builder never divides by zero.
    """

    chrom: str
    start: int
    end: int

    @property
    def span(self) -> int:
        """Width of the locus in base pairs (never negative)."""
        return max(0, self.end - self.start)


@dataclass
class RenderOptions:
    """
    A track's display settings (Requirement 9.2).

    ``colour`` is the already-validated colour to draw the track in (the caller
    validates with :func:`validate_colour` / :func:`apply_colour` before
    building). ``height`` is the track's drawn height in pixels. For a BigWig
    coverage track ``signal_min``/``signal_max`` bound the plotted signal range;
    they are ignored for feature tracks.
    """

    colour: str = DEFAULT_TRACK_COLOUR
    height: int = 60
    signal_min: Optional[float] = None
    signal_max: Optional[float] = None


@dataclass
class TrackSpec:
    """
    One enabled track to draw, with the data already read for the locus.

    ``kind`` is ``'bigwig'`` (drawn as a coverage plot, Requirement 8.6) or a
    feature kind (``'bed'``/``'gff'``, drawn as positioned boxes, Requirement
    8.7). ``signal`` is a list of ``(start, end, value)`` coverage bins for a
    BigWig track; ``features`` is a list of ``(start, end)`` intervals for a
    feature track. Coordinates are in the same base-pair space as the locus.
    The builder only *draws* what it is given, so the endpoint decides what data
    to read; an empty list renders an empty (but present) track lane.
    """

    name: str
    kind: str = 'bigwig'
    options: RenderOptions = field(default_factory=RenderOptions)
    signal: List[Tuple[float, float, float]] = field(default_factory=list)
    features: List[Tuple[int, int]] = field(default_factory=list)


@dataclass
class FigureLabel:
    """
    Member-supplied text drawn on the figure (Requirement 9.6).

    ``text`` is drawn verbatim, capped at 200 characters by the builder;
    ``x``/``y`` place it in the SVG user space at the member's chosen position.
    """

    text: str
    x: float = 10.0
    y: float = 20.0


@dataclass
class ViewSpec:
    """
    Everything the pure builder needs to render one figure.

    ``tracks`` is already filtered to the *enabled* tracks and already ordered
    in Browser_View order, so the builder draws them top-to-bottom as given and
    never sees a disabled track (Requirements 6.4, 9.1). ``genome_name`` and
    ``label`` are optional annotations.
    """

    locus: Locus
    tracks: List[TrackSpec] = field(default_factory=list)
    label: Optional[FigureLabel] = None
    genome_name: str = ''


# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------

_WIDTH = 900              # drawing width in px
_MARGIN_LEFT = 60
_MARGIN_RIGHT = 20
_PLOT_LEFT = _MARGIN_LEFT
_PLOT_RIGHT = _WIDTH - _MARGIN_RIGHT
_PLOT_WIDTH = _PLOT_RIGHT - _PLOT_LEFT
_HEADER_H = 40            # top band: locus title
_SCALE_H = 34            # coordinate-scale band
_TRACK_GAP = 10
_LABEL_BAND_H = 28        # bottom band reserved for the figure label


# ---------------------------------------------------------------------------
# Pure SVG building (Requirements 9.1, 9.2, 9.3, 9.6, 9.7, 9.8)
# ---------------------------------------------------------------------------


def build_svg(view: ViewSpec) -> str:
    """
    Render a :class:`ViewSpec` into a well-formed SVG document string.

    Draws, top to bottom: a header naming the genome and locus, the coordinate
    scale spanning the locus start/end (Requirement 9.7), then each enabled
    track in the order given — a BigWig as a coverage plot, a feature track as
    positioned boxes — using that track's Rendering_Options (Requirements 9.1,
    9.2). With no tracks the scale is still drawn and no track lane appears
    (the empty-view case, Requirement 9.3). A Figure_Label, if present, is drawn
    verbatim up to 200 characters at its position (Requirement 9.6).

    Every text node is XML-escaped and every attribute value comes from a
    validated number or a validated colour, so the result always parses under a
    conforming XML/SVG parser (Requirement 9.8), the empty view included.
    """
    locus = view.locus
    tracks = list(view.tracks)

    # Vertical layout: header, scale, each track lane, then the label band.
    top = _HEADER_H + _SCALE_H
    lanes: List[Tuple[TrackSpec, int, int]] = []  # (track, y_top, height)
    y = top
    for track in tracks:
        height = _clamp_int(track.options.height, lo=20, hi=600, default=60)
        lanes.append((track, y, height))
        y += height + _TRACK_GAP
    total_height = y + _LABEL_BAND_H

    body: List[str] = []
    body.append(_header_svg(view, locus))
    body.append(_scale_svg(locus))
    for track, lane_top, height in lanes:
        body.append(_track_svg(track, locus, lane_top, height))
    if view.label is not None:
        body.append(_label_svg(view.label))

    parts: List[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{_WIDTH}" height="{total_height}" '
        f'viewBox="0 0 {_WIDTH} {total_height}" '
        f'font-family="sans-serif">')
    # A white background so the figure is usable on any viewer.
    parts.append(
        f'<rect x="0" y="0" width="{_WIDTH}" height="{total_height}" '
        f'fill="#ffffff"/>')
    parts.extend(body)
    parts.append('</svg>')
    return '\n'.join(parts)


def _header_svg(view: ViewSpec, locus: Locus) -> str:
    """The top band: the genome name (when given) and the locus as text."""
    locus_text = f'{locus.chrom}:{locus.start}-{locus.end}'
    title = locus_text
    if view.genome_name:
        title = f'{view.genome_name}  {locus_text}'
    return (f'<text x="{_PLOT_LEFT}" y="24" font-size="16" '
            f'fill="#222222">{escape(title)}</text>')


def _scale_svg(locus: Locus) -> str:
    """
    The coordinate scale spanning the locus start and end (Requirement 9.7).

    A horizontal axis with ticks and labels at the start, end, and a few evenly
    spaced interior positions, so the figure always carries a scale — in the
    empty-view case too (Requirement 9.3).
    """
    y_axis = _HEADER_H + _SCALE_H - 12
    pieces: List[str] = [f'<g stroke="#444444" fill="#444444">']
    pieces.append(
        f'<line x1="{_PLOT_LEFT}" y1="{y_axis}" x2="{_PLOT_RIGHT}" '
        f'y2="{y_axis}" stroke-width="1"/>')

    ticks = _scale_ticks(locus)
    for coord in ticks:
        x = _coord_to_x(coord, locus)
        pieces.append(
            f'<line x1="{x:.2f}" y1="{y_axis}" x2="{x:.2f}" '
            f'y2="{y_axis + 6}" stroke-width="1"/>')
        pieces.append(
            f'<text x="{x:.2f}" y="{y_axis - 4}" font-size="10" '
            f'text-anchor="middle" stroke="none">{escape(str(coord))}</text>')
    pieces.append('</g>')
    return ''.join(pieces)


def _scale_ticks(locus: Locus) -> List[int]:
    """Start, end, and a few interior coordinates for the scale labels."""
    start, end = locus.start, locus.end
    if end <= start:
        return [start]
    divisions = 4
    coords = [start + round((end - start) * i / divisions)
              for i in range(divisions + 1)]
    # De-duplicate while keeping order (a tiny locus can collapse ticks).
    seen = set()
    ordered: List[int] = []
    for c in coords:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered


def _track_svg(track: TrackSpec, locus: Locus, lane_top: int, height: int) -> str:
    """
    One track lane: its name, a frame, and its data drawn in its colour.

    A BigWig track is a coverage plot scaled to its signal range; a feature
    track is a row of boxes positioned by start/end. The colour is the track's
    validated Rendering_Options colour (Requirements 9.2, 9.4).
    """
    colour = track.options.colour if validate_colour(track.options.colour) \
        else DEFAULT_TRACK_COLOUR
    pieces: List[str] = []
    pieces.append(f'<g>')
    # Lane frame.
    pieces.append(
        f'<rect x="{_PLOT_LEFT}" y="{lane_top}" width="{_PLOT_WIDTH}" '
        f'height="{height}" fill="none" stroke="#dddddd" stroke-width="1"/>')
    # A small colour swatch in the track's colour, next to its name. This makes
    # each enabled track carry its chosen colour in the figure even when the
    # drawn data is empty (region data is read elsewhere), so the colour the
    # member set is always shown (Requirements 9.2, 9.4).
    pieces.append(
        f'<rect x="{_PLOT_LEFT}" y="{lane_top - 11}" width="9" height="9" '
        f'fill={quoteattr(colour)}/>')
    # Track name.
    pieces.append(
        f'<text x="{_PLOT_LEFT + 13}" y="{lane_top - 2}" font-size="11" '
        f'fill="#333333">{escape(track.name)}</text>')

    if track.kind == 'bigwig':
        pieces.append(_coverage_svg(track, locus, lane_top, height, colour))
    else:
        pieces.append(_features_svg(track, locus, lane_top, height, colour))
    pieces.append('</g>')
    return ''.join(pieces)


def _coverage_svg(track: TrackSpec, locus: Locus, lane_top: int,
                  height: int, colour: str) -> str:
    """A BigWig coverage plot: filled bars per bin, scaled to the signal range."""
    bins = track.signal or []
    if not bins:
        return ''
    opt = track.options
    values = [v for _, _, v in bins if v is not None]
    lo = opt.signal_min if opt.signal_min is not None else (min(values) if values else 0.0)
    hi = opt.signal_max if opt.signal_max is not None else (max(values) if values else 1.0)
    if hi <= lo:
        hi = lo + 1.0
    base_y = lane_top + height
    rects: List[str] = []
    for start, end, value in bins:
        if value is None:
            continue
        x1 = _coord_to_x(start, locus)
        x2 = _coord_to_x(end, locus)
        w = max(0.5, x2 - x1)
        frac = (value - lo) / (hi - lo)
        frac = min(1.0, max(0.0, frac))
        bar_h = frac * height
        rects.append(
            f'<rect x="{x1:.2f}" y="{base_y - bar_h:.2f}" width="{w:.2f}" '
            f'height="{bar_h:.2f}" fill={quoteattr(colour)}/>')
    return ''.join(rects)


def _features_svg(track: TrackSpec, locus: Locus, lane_top: int,
                  height: int, colour: str) -> str:
    """A feature track: one box per interval, positioned by start/end."""
    features = track.features or []
    if not features:
        return ''
    box_h = max(6, min(height - 8, 16))
    y = lane_top + (height - box_h) / 2
    boxes: List[str] = []
    for start, end in features:
        x1 = _coord_to_x(start, locus)
        x2 = _coord_to_x(end, locus)
        w = max(1.0, x2 - x1)
        boxes.append(
            f'<rect x="{x1:.2f}" y="{y:.2f}" width="{w:.2f}" '
            f'height="{box_h}" rx="2" fill={quoteattr(colour)}/>')
    return ''.join(boxes)


def _label_svg(label: FigureLabel) -> str:
    """
    Draw the Figure_Label verbatim, capped at 200 characters (Requirement 9.6).

    The member's exact text is used, truncated to at most
    :data:`MAX_LABEL_CHARS` characters, XML-escaped, and placed at the label's
    position. Non-finite coordinates fall back to sensible defaults so the
    attribute is always a valid number and the document stays well-formed.
    """
    text = (label.text or '')[:MAX_LABEL_CHARS]
    x = _finite(label.x, default=10.0)
    y = _finite(label.y, default=20.0)
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-size="13" '
            f'fill="#000000">{escape(text)}</text>')


# ---------------------------------------------------------------------------
# Coordinate / numeric helpers
# ---------------------------------------------------------------------------


def _coord_to_x(coord: float, locus: Locus) -> float:
    """Map a base-pair coordinate to an x pixel inside the plot area."""
    span = locus.span
    if span <= 0:
        return float(_PLOT_LEFT)
    frac = (coord - locus.start) / span
    frac = min(1.0, max(0.0, frac))
    return _PLOT_LEFT + frac * _PLOT_WIDTH


def _clamp_int(value, lo: int, hi: int, default: int) -> int:
    """A bounded integer for a layout dimension, tolerating bad input."""
    try:
        n = int(value)
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, n))


def _finite(value, default: float) -> float:
    """A finite float for an SVG coordinate, falling back on bad input."""
    try:
        n = float(value)
    except (TypeError, ValueError):
        return default
    if n != n or n in (float('inf'), float('-inf')):  # NaN / inf
        return default
    return n
