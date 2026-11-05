"""
Format parsing for the Genome Browser's uploads.

This module knows how to read the file formats the module accepts, and nothing
about the database, sharing, or HTTP layers. Everything here answers one of two
questions: *what is this file* (format detection + the supported-formats list
that upload rejection quotes), and *what chromosomes/contigs and intervals does
it describe* (so the upload and alignment paths can check it against a reference
genome's ``chrom_sizes``).

Three groups of readers live here:

* **FASTA** — custom reference-genome upload. ``parse_fasta`` validates the
  file per Requirement 2.5 (at least one ``>`` header followed by one or more
  lines of nucleotide sequence) and returns each contig's name and length in
  base pairs, or a clear reason the file is not valid FASTA.

* **Track headers** — BigWig, BED, GFF/GFF3, BAM, and wig. These feed the
  alignment constraint (task 4.2): the chromosome/contig names a track mentions
  and, where the format carries them, the per-chromosome coordinate extents.

* **Format detection** — ``detect_format`` maps an uploaded filename to one of
  the supported formats or flags it unsupported, and ``SUPPORTED_FORMATS`` /
  ``supported_formats_message`` provide the list that a rejection message
  quotes (Requirement 4.6).

The one heavy, external-tool read — a BAM header — is factored behind
``read_bam_header`` with an injectable reader so property tests can stub it with
a small in-memory fixture instead of shelling out to ``samtools``. pyBigWig is
imported lazily for the same reason: the parsing logic a test exercises does not
drag in the native library unless a real BigWig is read.
"""

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Supported formats and detection (Requirement 4.6)
# ---------------------------------------------------------------------------
#
# Exactly the formats Requirement 4 names: BigWig stored directly, BED and
# GFF/GFF3 stored directly as feature tracks, and BAM and wig/wiggle accepted
# as conversion inputs. FASTA is the reference-genome upload format, kept
# separate from the track formats because it is never a track.

# Format key -> the file extensions that map to it. Detection is by extension
# because the module routes on format before reading the (possibly huge) body,
# and an uploaded BAM must never be sniffed into memory first.
TRACK_FORMAT_EXTENSIONS: Dict[str, Tuple[str, ...]] = {
    'bigwig': ('.bw', '.bigwig'),
    'bed': ('.bed',),
    'gff': ('.gff', '.gff3'),
    'bam': ('.bam',),
    'wig': ('.wig', '.wiggle', '.bedgraph', '.bg'),
}

# How each supported track format is handled, in the order a rejection message
# lists them. Used both for routing (stored directly vs converted) and to build
# the human-readable supported-formats list.
SUPPORTED_FORMATS: Tuple[Dict[str, object], ...] = (
    {'format': 'bigwig', 'label': 'BigWig', 'handling': 'stored directly',
     'track_type': 'bigwig'},
    {'format': 'bed', 'label': 'BED', 'handling': 'stored directly',
     'track_type': 'bed'},
    {'format': 'gff', 'label': 'GFF/GFF3', 'handling': 'stored directly',
     'track_type': 'gff'},
    {'format': 'bam', 'label': 'BAM', 'handling': 'converted to BigWig',
     'track_type': 'bigwig'},
    {'format': 'wig', 'label': 'wig/wiggle', 'handling': 'converted to BigWig',
     'track_type': 'bigwig'},
)

# Formats stored as-is vs routed through the conversion worker (task 8/9).
DIRECT_STORE_FORMATS = frozenset({'bigwig', 'bed', 'gff'})
CONVERSION_FORMATS = frozenset({'bam', 'wig'})

# The reference-genome upload format (never a track).
FASTA_EXTENSIONS = ('.fa', '.fasta', '.fna', '.fsa')


def detect_format(filename: str) -> Optional[str]:
    """
    Map an uploaded track filename to a supported format key, or None.

    Detection is by extension (case-insensitive), including a stripped ``.gz``
    so ``reads.bed.gz`` reads as BED. A None result means the format is not one
    the module accepts; the caller rejects it with ``supported_formats_message``
    (Requirement 4.6). FASTA is deliberately not returned here — it is the
    reference-genome upload format, matched separately by ``is_fasta_filename``.
    """
    ext = _extension(filename)
    for fmt, extensions in TRACK_FORMAT_EXTENSIONS.items():
        if ext in extensions:
            return fmt
    return None


def is_fasta_filename(filename: str) -> bool:
    """Whether a filename looks like FASTA (reference-genome upload)."""
    return _extension(filename) in FASTA_EXTENSIONS


def track_type_for_format(fmt: str) -> Optional[str]:
    """
    The stored ``tracks.type`` a format produces, or None for an unknown format.

    BAM and wig both end up as a BigWig track once converted, so they map to
    ``bigwig`` too; the caller tells store-now from convert with
    ``DIRECT_STORE_FORMATS`` / ``CONVERSION_FORMATS``.
    """
    for spec in SUPPORTED_FORMATS:
        if spec['format'] == fmt:
            return str(spec['track_type'])
    return None


def supported_formats_message() -> str:
    """
    The supported-formats list an unsupported-format rejection quotes (4.6).

    Phrased so the message tells the member both which formats are accepted and
    how each is handled, since BAM/wig are accepted but converted rather than
    stored.
    """
    parts = [f"{spec['label']} ({spec['handling']})" for spec in SUPPORTED_FORMATS]
    return 'Supported track formats: ' + ', '.join(parts) + '.'


def _extension(filename: str) -> str:
    """Lower-cased final extension, ignoring a trailing ``.gz``."""
    name = (filename or '').strip().lower()
    if name.endswith('.gz'):
        name = name[:-len('.gz')]
    _, ext = os.path.splitext(name)
    return ext


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------


@dataclass
class FastaResult:
    """
    The outcome of validating and parsing a FASTA upload.

    ``valid`` is the single gate the upload path checks; on an invalid file
    ``reason`` carries the member-facing explanation (Requirement 2.5) and
    ``contigs`` is empty. On a valid file ``contigs`` lists ``(name, length_bp)``
    in file order, de-duplicated check left to the caller against existing
    assemblies.
    """

    valid: bool
    reason: str = ''
    contigs: List[Tuple[str, int]] = field(default_factory=list)

    @property
    def chrom_sizes(self) -> Dict[str, int]:
        """Contig name -> length, convenient for the ``chrom_sizes`` write."""
        return {name: length for name, length in self.contigs}


@dataclass
class HeaderResult:
    """
    The chromosomes/contigs a track file references, for the alignment check.

    ``chrom_lengths`` maps each chromosome name to the declared length where the
    format carries one (BigWig, BAM, wig fixed/variableStep ``chrom=`` with a
    span). ``chrom_extents`` maps each chromosome to the largest coordinate the
    file actually reaches, which the coordinate-in-range check (task 4.2) tests
    against the reference length; a format that declares a length but no
    intervals (a bare BigWig header) leaves the extent equal to the length.
    """

    chrom_lengths: Dict[str, int] = field(default_factory=dict)
    chrom_extents: Dict[str, int] = field(default_factory=dict)
    reason: str = ''

    @property
    def chroms(self) -> List[str]:
        """Every chromosome/contig name the file mentions, in a stable order."""
        names = list(self.chrom_lengths)
        for name in self.chrom_extents:
            if name not in self.chrom_lengths:
                names.append(name)
        return names


# ---------------------------------------------------------------------------
# FASTA (Requirements 2.3, 2.5)
# ---------------------------------------------------------------------------

# A nucleotide sequence line: IUPAC bases (ACGT + ambiguity codes), the gap
# character, and case-insensitive, nothing else. ``U`` allows RNA; ``N`` and the
# ambiguity codes allow soft-masked / unknown bases real assemblies contain.
_NUCLEOTIDE_RE = re.compile(r'^[ACGTUNRYSWKMBDHV.-]+$', re.IGNORECASE)


def parse_fasta(text: str) -> FastaResult:
    """
    Validate a FASTA file and return its contig names and lengths (2.3, 2.5).

    Valid FASTA, per Requirement 2.5, is at least one header line starting with
    ``>`` followed by one or more lines of nucleotide sequence characters. This
    reader enforces exactly that: a file with no header, a header with no
    sequence under it, or a sequence line with a non-nucleotide character is
    invalid, and ``reason`` names which. Each contig's length is the total count
    of sequence characters under its header (gaps and newlines excluded from the
    line, gap characters themselves counted as bases of the assembly).

    ``text`` is the decoded file body. Callers reading a real upload decode it
    (FASTA is ASCII) before handing it here, which keeps this reader pure and
    trivially testable.
    """
    header: Optional[str] = None
    length = 0
    seen_sequence = False
    contigs: List[Tuple[str, int]] = []

    def flush() -> Optional[str]:
        """Close the current contig; return a reason if it had no sequence."""
        nonlocal header, length, seen_sequence
        if header is not None:
            if length == 0:
                return (f'Header "{header}" has no sequence lines; a FASTA '
                        f'header must be followed by nucleotide sequence.')
            contigs.append((header, length))
        header, length, seen_sequence = None, 0, False
        return None

    for raw in text.splitlines():
        line = raw.strip()
        if line == '':
            continue  # blank lines separate records / pad a file harmlessly
        if line.startswith('>'):
            reason = flush()
            if reason:
                return FastaResult(valid=False, reason=reason)
            header = line[1:].strip().split()[0] if line[1:].strip() else ''
            if header == '':
                return FastaResult(
                    valid=False,
                    reason='A FASTA header line (">") has no sequence name.')
            continue
        if header is None:
            return FastaResult(
                valid=False,
                reason='File does not start with a FASTA header line (">"); '
                       'it is not valid FASTA.')
        if not _NUCLEOTIDE_RE.match(line):
            bad = _first_bad_char(line)
            return FastaResult(
                valid=False,
                reason=(f'Sequence for "{header}" contains a non-nucleotide '
                        f'character {bad!r}; FASTA sequence must be nucleotide '
                        f'characters only.'))
        length += _sequence_length(line)
        seen_sequence = True

    reason = flush()
    if reason:
        return FastaResult(valid=False, reason=reason)

    if not contigs:
        return FastaResult(
            valid=False,
            reason='File contains no FASTA records (no ">" header followed by '
                   'sequence); it is not valid FASTA.')
    return FastaResult(valid=True, contigs=contigs)


def _sequence_length(line: str) -> int:
    """Count sequence characters on a line (every char; the line is validated)."""
    return len(line)


def _first_bad_char(line: str) -> str:
    for ch in line:
        if not _NUCLEOTIDE_RE.match(ch):
            return ch
    return ''


# ---------------------------------------------------------------------------
# BigWig header (Requirement 3.3; alignment check feeds task 4.2)
# ---------------------------------------------------------------------------


def read_bigwig_header(path: str) -> HeaderResult:
    """
    Read a BigWig's chromosome names and lengths with pyBigWig.

    pyBigWig's ``chroms()`` returns the header's ``name -> length`` map, which is
    exactly the chromosome set and per-chromosome extent the alignment check
    needs: a BigWig's declared chromosome length is the furthest coordinate its
    signal can reach, so the length doubles as the extent. pyBigWig is imported
    here, not at module load, so tests that never touch a real BigWig do not need
    the native library present.
    """
    try:
        import pyBigWig  # imported lazily; native lib only needed for real reads
    except ImportError as exc:  # pragma: no cover - environment-dependent
        return HeaderResult(reason=f'pyBigWig is not available: {exc}')

    bw = None
    try:
        bw = pyBigWig.open(path)
        chroms = bw.chroms() or {}
    except (RuntimeError, OSError) as exc:
        return HeaderResult(reason=f'Could not read BigWig header: {exc}')
    finally:
        if bw is not None:
            bw.close()

    lengths = {str(name): int(length) for name, length in chroms.items()}
    return HeaderResult(chrom_lengths=lengths, chrom_extents=dict(lengths))


# ---------------------------------------------------------------------------
# BED / GFF / GFF3 interval scan (feeds the alignment check)
# ---------------------------------------------------------------------------


def scan_bed_intervals(text: str) -> HeaderResult:
    """
    Scan a BED file for the chromosomes and the furthest coordinate each reaches.

    BED is ``chrom<TAB>start<TAB>end`` with 0-based, half-open coordinates; the
    chromosome's extent is the largest ``end`` seen for it. Track/browser/comment
    lines and blank lines are skipped. A line that names a chromosome but cannot
    be parsed into coordinates still registers the chromosome (so a name mismatch
    is still caught) but contributes no extent.
    """
    extents: Dict[str, int] = {}
    for raw in text.splitlines():
        line = raw.rstrip('\n')
        if _is_feature_comment(line):
            continue
        fields = line.split('\t') if '\t' in line else line.split()
        if not fields:
            continue
        chrom = fields[0]
        extents.setdefault(chrom, 0)
        end = _int_or_none(fields[2]) if len(fields) >= 3 else None
        if end is not None:
            extents[chrom] = max(extents[chrom], end)
    return HeaderResult(chrom_extents=extents)


def scan_gff_intervals(text: str) -> HeaderResult:
    """
    Scan a GFF/GFF3 file for chromosomes and their furthest coordinate.

    GFF3 columns are ``seqid source type start end ...`` with 1-based, inclusive
    coordinates, so a feature's extent is its ``end``. ``##sequence-region``
    pragmas, when present, declare a seqid's end coordinate up front and are
    honoured as an extent; all other ``##``/``#`` lines and blanks are skipped.
    """
    extents: Dict[str, int] = {}
    for raw in text.splitlines():
        line = raw.rstrip('\n')
        if line.startswith('##sequence-region'):
            parts = line.split()
            if len(parts) >= 4:
                end = _int_or_none(parts[3])
                if end is not None:
                    extents[parts[1]] = max(extents.get(parts[1], 0), end)
            continue
        if _is_feature_comment(line):
            continue
        fields = line.split('\t') if '\t' in line else line.split()
        if not fields:
            continue
        seqid = fields[0]
        extents.setdefault(seqid, 0)
        end = _int_or_none(fields[4]) if len(fields) >= 5 else None
        if end is not None:
            extents[seqid] = max(extents[seqid], end)
    return HeaderResult(chrom_extents=extents)


def _is_feature_comment(line: str) -> bool:
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
# wig header (feeds the worker's validation path)
# ---------------------------------------------------------------------------


def scan_wig_chroms(text: str) -> HeaderResult:
    """
    Read the chromosome declarations from a wig file (fixed/variable/bedGraph).

    ``fixedStep``/``variableStep`` declare ``chrom=`` on a header line, and the
    furthest coordinate a chromosome reaches is tracked so the worker can check
    it against the reference length: fixedStep advances by ``step`` from ``start``
    across its data lines, variableStep reads a position per line, and bedGraph
    lines carry ``chrom start end value`` directly. Enough to catch a chromosome
    name mismatch or an out-of-range coordinate before heavy conversion.
    """
    extents: Dict[str, int] = {}
    mode: Optional[str] = None
    chrom: Optional[str] = None
    step = 1
    span = 1
    position = 0

    for raw in text.splitlines():
        line = raw.strip()
        if line == '' or line.startswith('#') or line.startswith('track '):
            continue
        lower = line.lower()
        if lower.startswith('fixedstep'):
            attrs = _wig_attrs(line)
            mode, chrom = 'fixed', attrs.get('chrom')
            step = _int_or_none(attrs.get('step', '')) or 1
            span = _int_or_none(attrs.get('span', '')) or 1
            position = _int_or_none(attrs.get('start', '')) or 1
            if chrom:
                extents.setdefault(chrom, 0)
            continue
        if lower.startswith('variablestep'):
            attrs = _wig_attrs(line)
            mode, chrom = 'variable', attrs.get('chrom')
            span = _int_or_none(attrs.get('span', '')) or 1
            if chrom:
                extents.setdefault(chrom, 0)
            continue
        fields = line.split()
        if mode == 'fixed' and chrom:
            extents[chrom] = max(extents[chrom], position + span - 1)
            position += step
            continue
        if mode == 'variable' and chrom and fields:
            pos = _int_or_none(fields[0])
            if pos is not None:
                extents[chrom] = max(extents[chrom], pos + span - 1)
            continue
        # bedGraph: chrom start end value
        if len(fields) >= 3:
            bg_chrom = fields[0]
            end = _int_or_none(fields[2])
            extents.setdefault(bg_chrom, 0)
            if end is not None:
                extents[bg_chrom] = max(extents[bg_chrom], end)
    return HeaderResult(chrom_extents=extents)


def _wig_attrs(line: str) -> Dict[str, str]:
    """Parse ``key=value`` tokens from a wig declaration line."""
    attrs: Dict[str, str] = {}
    for token in line.split()[1:]:
        if '=' in token:
            key, _, value = token.partition('=')
            attrs[key.strip().lower()] = value.strip()
    return attrs


# ---------------------------------------------------------------------------
# BAM header (external tool; factored so tests can stub it)
# ---------------------------------------------------------------------------


def _samtools_header_reader(path: str) -> str:
    """
    Default BAM-header reader: ``samtools view -H`` (or pysam if present).

    Kept private and injectable so the parsing in ``read_bam_header`` can be
    tested against a captured header without samtools or pysam installed. Prefer
    pysam when it is importable (no subprocess), otherwise shell out to
    ``samtools view -H`` and return its text.
    """
    try:
        import pysam  # type: ignore
        return str(pysam.view('-H', path))
    except ImportError:
        pass
    if shutil.which('samtools') is None:
        raise RuntimeError('samtools is not available to read the BAM header')
    completed = subprocess.run(
        ['samtools', 'view', '-H', path],
        capture_output=True, text=True, check=True)
    return completed.stdout


def parse_bam_header_text(header_text: str) -> HeaderResult:
    """
    Parse a SAM/BAM header's ``@SQ`` lines into chromosome names and lengths.

    Each reference sequence is an ``@SQ`` line with ``SN:<name>`` and
    ``LN:<length>`` tab-separated fields; a BAM's declared sequence length is the
    furthest coordinate a read can map to, so it serves as both the declared
    length and the extent for the alignment check. Pure text in, so property
    tests feed it a captured header directly.
    """
    lengths: Dict[str, int] = {}
    for raw in header_text.splitlines():
        if not raw.startswith('@SQ'):
            continue
        name: Optional[str] = None
        length: Optional[int] = None
        for field_token in raw.split('\t')[1:]:
            if field_token.startswith('SN:'):
                name = field_token[3:]
            elif field_token.startswith('LN:'):
                length = _int_or_none(field_token[3:])
        if name is not None and length is not None:
            lengths[name] = length
    if not lengths:
        return HeaderResult(reason='BAM header has no @SQ reference sequences.')
    return HeaderResult(chrom_lengths=lengths, chrom_extents=dict(lengths))


def read_bam_header(
    path: str,
    reader: Callable[[str], str] = _samtools_header_reader,
) -> HeaderResult:
    """
    Read a BAM's chromosome names and lengths for the worker validation path.

    The heavy, external-tool read is isolated in ``reader`` (defaulting to
    ``samtools``/pysam) so a test injects a stub returning a captured header and
    exercises the ``@SQ`` parsing without the tool installed. A failure to run
    the reader comes back as a ``reason`` rather than an exception, so the worker
    can fail the job with a descriptive error.
    """
    try:
        header_text = reader(path)
    except (RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        return HeaderResult(reason=f'Could not read BAM header: {exc}')
    return parse_bam_header_text(header_text)


# ---------------------------------------------------------------------------
# Convenience: header read dispatched by format (worker/upload validation)
# ---------------------------------------------------------------------------


def read_track_header(fmt: str, path: str) -> HeaderResult:
    """
    Read a track file's chromosome header for the format, dispatching by key.

    A single entry point for the upload (BigWig direct-add) and worker (BAM/wig)
    validation paths: it reads the file at ``path`` as the given format and
    returns the chromosomes and extents the alignment check needs. Text formats
    are read from disk here; the BAM path uses the default (stubbable) reader.
    """
    if fmt == 'bigwig':
        return read_bigwig_header(path)
    if fmt == 'bam':
        return read_bam_header(path)
    if fmt in ('bed',):
        return scan_bed_intervals(_read_text(path))
    if fmt in ('gff',):
        return scan_gff_intervals(_read_text(path))
    if fmt in ('wig',):
        return scan_wig_chroms(_read_text(path))
    return HeaderResult(reason=f'Unsupported format for header read: {fmt!r}')


def _read_text(path: str) -> str:
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        return handle.read()
