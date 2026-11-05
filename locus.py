"""
Locus parsing and navigation math for the Genome Browser (Requirement 8).

This module is the single place the browser's coordinate arithmetic lives, and
it is deliberately *pure*: it knows about a genome only as a ``chrom_sizes``
mapping (``chrom -> length_bp``), and it touches no database, filesystem, or
HTTP layer. Everything an IGV-style navigation needs — parse a locus string,
clamp a region to a chromosome's bounds, centre on a locus, zoom while keeping
the centre fixed, pan while keeping the span fixed — is a plain function over a
small :class:`Locus` value, so the property tests in task 11.3 can drive it
directly and the region endpoint (``genome/regions.py``) and the Browser_View
persistence endpoint (task 11.2) can reuse exactly the same rules.

The one invariant everything here maintains is ``0 <= start <= end <= length``
(Requirement 8.3): a parsed or navigated locus is always a valid half-open
``[start, end)`` window that fits inside its chromosome. Two things break that
invariant from the outside and are handled explicitly:

- **An unknown chromosome** (Requirement 8.2). A chromosome name that is not in
  ``chrom_sizes`` cannot be clamped — there is no length to clamp to — so
  parsing/clamping raises :class:`UnknownChrom`, whose message names the
  offending chromosome. The endpoint turns that into a rejection naming the
  unknown chromosome/contig.

- **An out-of-range coordinate** (Requirement 8.3). A coordinate past the
  chromosome's end, or a negative start, is *clamped* to the chromosome bounds
  rather than rejected, so navigating off the end of a chromosome snaps to the
  end instead of erroring.

Coordinates are 0-based, half-open (``[start, end)``), the same convention the
renderer (:mod:`render`) and BigWig/BED use, so a whole chromosome of length
``L`` is the locus ``0..L`` and its span is ``L``. A locus string written the
way a member types it — ``Chr1:1,000-5,000`` with thousands separators, or a
bare ``Chr1`` — is 1-based inclusive in the *display* form, so ``parse`` reads
the typed start as 1-based and stores it 0-based; a bare chromosome spans the
whole chromosome.
"""

import re
from dataclasses import dataclass
from typing import Mapping, Optional


class UnknownChrom(ValueError):
    """
    A locus names a chromosome/contig not in the reference genome (Req 8.2).

    Raised by :func:`parse` and :func:`clamp` when the chromosome is absent from
    ``chrom_sizes``; :attr:`chrom` is the offending name so the caller can name
    it in the rejection message the requirement asks for. It subclasses
    ``ValueError`` so a caller that only wants "bad locus" can catch that.
    """

    def __init__(self, chrom: str):
        self.chrom = chrom
        super().__init__(
            f'Unknown chromosome/contig "{chrom}"; it is not in this '
            f'reference genome.')


class LocusParseError(ValueError):
    """
    A locus string is not shaped like ``Chrom`` or ``Chrom:start-end`` (Req 8.1).

    Distinct from :class:`UnknownChrom`: the string itself is malformed (empty,
    no chromosome, non-numeric coordinates) rather than naming a chromosome that
    merely does not exist.
    """


@dataclass(frozen=True)
class Locus:
    """
    A half-open genomic window ``[start, end)`` on one chromosome/contig.

    ``chrom`` is the chromosome/contig name (matched case-sensitively against
    ``chrom_sizes``, as the alignment constraint requires elsewhere).
    ``start``/``end`` are 0-based base-pair coordinates with ``start <= end``.
    A :class:`Locus` produced by this module always satisfies
    ``0 <= start <= end <= length`` for its chromosome (Requirement 8.3); the
    dataclass itself does not enforce that (it is a plain value), the functions
    that build one do.
    """

    chrom: str
    start: int
    end: int

    @property
    def span(self) -> int:
        """The window width in base pairs (``end - start``), never negative."""
        return max(0, self.end - self.start)

    @property
    def centre(self) -> float:
        """The midpoint coordinate, used to keep zoom centred (Req 8.4)."""
        return (self.start + self.end) / 2.0

    def format(self) -> str:
        """
        Render back to the display form ``Chrom:start-end`` (1-based inclusive).

        The inverse of :func:`parse`: the stored 0-based half-open window is
        shown 1-based inclusive the way a member reads a locus, with thousands
        separators for legibility. A zero-width window shows a single base.
        """
        display_start = self.start + 1
        display_end = max(self.start + 1, self.end)
        return f'{self.chrom}:{display_start:,}-{display_end:,}'


# A locus string: a chromosome name, optionally followed by ``:start-end`` with
# optional thousands separators and surrounding whitespace. The chromosome part
# is everything up to the first colon, so names with dots/underscores work; a
# colon only introduces the coordinate range.
_LOCUS_RE = re.compile(
    r'^\s*(?P<chrom>[^:\s][^:]*?)\s*'
    r'(?::\s*(?P<start>[\d,]+)\s*-\s*(?P<end>[\d,]+)\s*)?$')


def parse(locus_str: str, chrom_sizes: Mapping[str, int]) -> Locus:
    """
    Parse a member-typed locus into a clamped :class:`Locus` (Req 8.1, 8.2, 8.3).

    Accepts ``Chr1:1,000-5,000`` (1-based inclusive display coordinates, commas
    ignored) or a bare ``Chr1`` (the whole chromosome, ``0..length``). The
    chromosome must exist in ``chrom_sizes`` or :class:`UnknownChrom` is raised
    naming it (Requirement 8.2). A coordinate range outside the chromosome is
    clamped to ``0..length`` rather than rejected (Requirement 8.3); a reversed
    range (``end < start``) is normalised. The returned locus always satisfies
    ``0 <= start <= end <= length``.

    ``chrom_sizes`` maps each chromosome/contig name to its length in base
    pairs; it is read, never mutated.
    """
    match = _LOCUS_RE.match(locus_str or '')
    if not match:
        raise LocusParseError(
            f'"{locus_str}" is not a valid locus; expected a chromosome name '
            f'optionally followed by ":start-end".')

    chrom = match.group('chrom').strip()
    if chrom not in chrom_sizes:
        raise UnknownChrom(chrom)
    length = int(chrom_sizes[chrom])

    if match.group('start') is None:
        # A bare chromosome spans the whole chromosome, 0..length.
        return Locus(chrom=chrom, start=0, end=length)

    # Display coordinates are 1-based inclusive; store 0-based half-open.
    display_start = int(match.group('start').replace(',', ''))
    display_end = int(match.group('end').replace(',', ''))
    start = display_start - 1
    end = display_end
    if end < start:
        start, end = end, start
    return clamp(Locus(chrom=chrom, start=start, end=end), chrom_sizes)


def clamp(locus: Locus, chrom_sizes: Mapping[str, int]) -> Locus:
    """
    Clamp a locus to ``0 <= start <= end <= length`` for its chromosome (8.3).

    The chromosome must exist in ``chrom_sizes`` or :class:`UnknownChrom` is
    raised naming it (Requirement 8.2) — there is no length to clamp an unknown
    chromosome to. A start below 0 snaps to 0 and an end past the chromosome
    length snaps to the length, so navigating off either edge lands on the edge
    rather than erroring. A reversed window is normalised so ``start <= end``.
    """
    if locus.chrom not in chrom_sizes:
        raise UnknownChrom(locus.chrom)
    length = int(chrom_sizes[locus.chrom])

    start, end = locus.start, locus.end
    if end < start:
        start, end = end, start
    start = max(0, min(start, length))
    end = max(0, min(end, length))
    if start > end:            # both past the same edge; collapse to it
        start = end
    return Locus(chrom=locus.chrom, start=start, end=end)


def center_on(locus: Locus, chrom_sizes: Mapping[str, int],
              span: Optional[int] = None) -> Locus:
    """
    Centre the view on ``locus`` (Requirement 8.1), keeping a chosen span.

    Used when a member enters a locus and the browser displays the view centred
    on it. The returned window has the given ``span`` (defaulting to the locus's
    own span) centred on the locus's midpoint, then clamped to the chromosome
    bounds (Requirement 8.3). Clamping may shift the window inward near an edge;
    the span is preserved where the chromosome is wide enough to hold it and
    shrinks only when the whole chromosome is narrower than the requested span.
    """
    width = locus.span if span is None else max(0, int(span))
    centre = locus.centre
    start = int(round(centre - width / 2.0))
    end = start + width
    return _place_window(locus.chrom, start, end, width, chrom_sizes)


def zoom(locus: Locus, factor: float, chrom_sizes: Mapping[str, int]) -> Locus:
    """
    Zoom in/out while keeping the view centred (Requirement 8.4).

    ``factor`` scales the span: ``factor < 1`` zooms in (a smaller span),
    ``factor > 1`` zooms out (a larger span); ``factor == 1`` is a no-op. The
    centre of the window is preserved — the new window is the scaled span placed
    symmetrically around the current centre — then clamped to the chromosome
    bounds (Requirement 8.3). A zoom-out that would exceed the chromosome is
    clamped to the whole chromosome. ``factor`` must be positive.

    Preserving the centre is exact when the clamp does not bite; near a
    chromosome edge the clamped window keeps its span and butts against the edge,
    which is the IGV behaviour (you cannot scroll past the end, but the zoom
    level is kept).
    """
    if factor <= 0:
        raise ValueError('zoom factor must be positive')
    if locus.chrom not in chrom_sizes:
        raise UnknownChrom(locus.chrom)

    new_span = max(1, int(round(locus.span * factor)))
    centre = locus.centre
    start = int(round(centre - new_span / 2.0))
    end = start + new_span
    return _place_window(locus.chrom, start, end, new_span, chrom_sizes)


def pan(locus: Locus, delta: int, chrom_sizes: Mapping[str, int]) -> Locus:
    """
    Pan left/right by ``delta`` base pairs, keeping the span (Requirement 8.5).

    A positive ``delta`` shifts the window right (toward higher coordinates), a
    negative ``delta`` left. The span is preserved exactly: the window slides as
    a unit and, when it hits a chromosome edge, stops against the edge with its
    span intact rather than shrinking (Requirement 8.3 clamp). Only when the
    span is wider than the whole chromosome does the window shrink to the
    chromosome.
    """
    if locus.chrom not in chrom_sizes:
        raise UnknownChrom(locus.chrom)
    span = locus.span
    start = locus.start + int(delta)
    end = start + span
    return _place_window(locus.chrom, start, end, span, chrom_sizes)


def pan_direction(locus: Locus, direction: str, chrom_sizes: Mapping[str, int],
                  fraction: float = 0.5) -> Locus:
    """
    Pan by a direction word (``'left'``/``'right'``), a span fraction each step.

    A convenience over :func:`pan` for the browser's pan buttons: it shifts by
    ``fraction`` of the current span (default half a screen, the IGV default) in
    the named direction, preserving the span (Requirement 8.5). An unrecognised
    direction is a no-op (returns the locus clamped).
    """
    step = max(1, int(round(locus.span * fraction)))
    key = (direction or '').strip().lower()
    if key in ('left', 'l', '<', 'back'):
        return pan(locus, -step, chrom_sizes)
    if key in ('right', 'r', '>', 'forward'):
        return pan(locus, step, chrom_sizes)
    return clamp(locus, chrom_sizes)


def _place_window(chrom: str, start: int, end: int, span: int,
                  chrom_sizes: Mapping[str, int]) -> Locus:
    """
    Place a fixed-span window inside the chromosome, sliding it off no edge.

    Shared by centre/zoom/pan: given a desired ``[start, end)`` of width
    ``span``, slide the window so it fits inside ``0..length`` while keeping its
    span, and only shrink the span when the chromosome is narrower than the
    window. This is what makes zoom keep its level and pan keep its span at the
    chromosome edges (Requirements 8.3, 8.4, 8.5).
    """
    if chrom not in chrom_sizes:
        raise UnknownChrom(chrom)
    length = int(chrom_sizes[chrom])

    if span >= length:
        # The window is at least the whole chromosome; show all of it.
        return Locus(chrom=chrom, start=0, end=length)

    if start < 0:
        start, end = 0, span
    elif end > length:
        end, start = length, length - span
    return Locus(chrom=chrom, start=start, end=end)
