"""
The BigWig alignment constraint for the Genome Browser.

A track aligns to a reference genome only when every chromosome/contig name it
references is an *exact, case-sensitive* match of a name in the genome and every
coordinate it reaches on a matched chromosome lies within that chromosome's
length. This module is the one place that decision is made, so the direct-add
path (a BigWig upload, checked in the request) and the Conversion_Worker
validation path (a BAM or wig input, checked before heavy conversion) agree by
construction rather than by two copies of the rule staying in step
(Requirements 3.3, 3.4, 3.5, 13.7, 13.8).

The module is pure logic. It takes a track's parsed header (a
``readers.HeaderResult`` — its ``chroms`` names and its per-chromosome
``chrom_extents``) and the reference genome's ``chrom_sizes`` mapping, and
returns an :class:`AlignmentResult` saying whether the track is accepted and,
when it is not, exactly which chromosome names did not match and which
coordinates ran past a chromosome's end. It touches no database and no
filesystem: the **persist-nothing contract** lives in the caller, which must
store no track row and no file when the result is a rejection (Requirements
3.4, 3.5, 13.8). Returning the decision rather than acting on it is what lets
both callers honour that contract the same way.
"""

from dataclasses import dataclass, field
from typing import Dict, List

from readers import HeaderResult


@dataclass
class AlignmentResult:
    """
    The outcome of checking a track header against a reference genome.

    ``accepted`` is the single gate a caller branches on. When it is false the
    caller rejects the track and persists nothing (Requirements 3.4, 3.5,
    13.8); the lists below say why, so the caller can surface a message that
    names every offender:

    * ``unmatched_chroms`` — chromosome/contig names in the track that do not
      exactly (case-sensitively) match any genome chromosome, in the order the
      header mentions them (Requirement 3.4).
    * ``out_of_range`` — one entry per matched chromosome whose coordinates run
      past the genome's length, each a :class:`CoordinateViolation` naming the
      chromosome, the offending coordinate, and the allowed length
      (Requirement 3.5).

    ``messages`` is the ready-to-surface, human-readable rendering of those two
    lists, one line per offender; it is empty when the track is accepted.
    """

    accepted: bool
    unmatched_chroms: List[str] = field(default_factory=list)
    out_of_range: List["CoordinateViolation"] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)

    @property
    def error_message(self) -> str:
        """The messages joined into one string, or '' when accepted."""
        return ' '.join(self.messages)


@dataclass(frozen=True)
class CoordinateViolation:
    """
    A coordinate on a matched chromosome that exceeds the genome's length.

    ``chrom`` is the matched chromosome name, ``coordinate`` is the furthest
    coordinate the track reaches there, and ``length`` is that chromosome's
    length in the reference genome — the bound the coordinate exceeded
    (Requirement 3.5).
    """

    chrom: str
    coordinate: int
    length: int


def check_alignment(
    header: HeaderResult,
    chrom_sizes: Dict[str, int],
) -> AlignmentResult:
    """
    Check a track's header against a reference genome's ``chrom_sizes``.

    The track is accepted *if and only if* every chromosome/contig name the
    header mentions is an exact case-sensitive match of a name in
    ``chrom_sizes`` **and** every coordinate the header reaches on a matched
    chromosome is within that chromosome's length (Requirements 3.3, 3.5). The
    same call serves both the direct-add path and the worker validation path,
    so the two always agree (Requirements 13.7, 13.8).

    ``header`` is a :class:`readers.HeaderResult`: ``header.chroms`` is every
    name the file references and ``header.chrom_extents`` is the furthest
    coordinate it reaches per chromosome (a bare BigWig header, with no
    intervals, has its extent equal to its declared length, so a well-formed
    BigWig never trips the range check against its own genome). ``chrom_sizes``
    maps each genome chromosome name to its length in base pairs.

    On any violation the result lists *every* unmatched name (Requirement 3.4)
    and *every* out-of-range coordinate (Requirement 3.5) — the check does not
    stop at the first offender, so the member sees the full picture in one
    message. This function only decides; the caller persists nothing on a
    rejection (the persist-nothing contract, Requirements 3.4, 3.5, 13.8).

    An unparsable header (``header.reason`` set, e.g. a BAM with no ``@SQ``
    lines or a BigWig that would not open) is itself a rejection: the track's
    chromosomes could not be read, so it cannot be shown to align, and the
    reason is surfaced as the message.
    """
    # A header the reader could not parse cannot be shown to align. Reject it
    # with the reader's reason rather than silently treating an empty header as
    # a match-everything pass.
    if header.reason:
        return AlignmentResult(accepted=False, messages=[header.reason])

    unmatched: List[str] = []
    violations: List[CoordinateViolation] = []

    for chrom in header.chroms:
        if chrom not in chrom_sizes:
            # Exact, case-sensitive: `in` on the dict is case-sensitive, so a
            # name differing only in case (e.g. "chr1" vs "Chr1") is unmatched.
            unmatched.append(chrom)
            continue
        # Matched name: check the furthest coordinate the track reaches here
        # against the genome's length for this chromosome.
        length = chrom_sizes[chrom]
        extent = header.chrom_extents.get(chrom)
        if extent is not None and extent > length:
            violations.append(
                CoordinateViolation(chrom=chrom, coordinate=extent,
                                    length=length))

    if not unmatched and not violations:
        return AlignmentResult(accepted=True)

    return AlignmentResult(
        accepted=False,
        unmatched_chroms=unmatched,
        out_of_range=violations,
        messages=_build_messages(unmatched, violations),
    )


def _build_messages(
    unmatched: List[str],
    violations: List["CoordinateViolation"],
) -> List[str]:
    """
    Render the offenders into member-facing messages, one line per concern.

    The unmatched-name line names every unmatched chromosome/contig
    (Requirement 3.4); each out-of-range line identifies the chromosome and the
    out-of-range coordinate against its length (Requirement 3.5).
    """
    messages: List[str] = []
    if unmatched:
        names = ', '.join(unmatched)
        plural = 's' if len(unmatched) > 1 else ''
        messages.append(
            f'Chromosome/contig name{plural} not found in the reference '
            f'genome (names must match exactly, including case): {names}.')
    for violation in violations:
        messages.append(
            f'Coordinate {violation.coordinate} on chromosome '
            f'"{violation.chrom}" exceeds its length of {violation.length} '
            f'in the reference genome.')
    return messages
