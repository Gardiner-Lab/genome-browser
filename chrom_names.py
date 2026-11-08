"""
Chromosome-name reconciliation: the check every track upload runs against its
reference genome, and the renames it can make on its own.

The same chromosome is written many ways. TAIR10's FASTA calls Arabidopsis
chromosome 1 ``1`` and its mitochondrion ``mitochondria``; Araport11's GFF3
calls them ``Chr1`` and ``ChrM``; UCSC writes ``chr1`` and Ensembl ``1`` and
``Mt``. A track whose names differ from its genome's only in that style is the
same data and should simply be drawn, so this module works out, for each name a
track uses, which of the genome's chromosomes it means:

* an **exact** match is kept as it is;
* a name that differs only in style is **renamed** to the genome's name — the
  case, a ``chr``/``chrom``/``chromosome`` prefix and the separator after it,
  leading zeros on a number, and the common spellings of the organelles
  (``M``, ``MT``, ``mito``, ``mitochondria``, ``mitochondrion`` are one, and
  ``C``, ``Pt``, ``Cp``, ``chloroplast``, ``plastid`` are the other);
* anything else is **unmatched** and left for the caller to refuse or report.

A rename is made only when it is unambiguous: if two of the genome's own names
reduce to the same key (a genome with both ``chr1`` and ``Chr1``), neither is a
rename target, and a track name that reduces to that key stays unmatched.

The module is pure logic, like ``align``: it decides, and the callers rewrite
files and choose what to refuse. :func:`rename_header` applies a decision to a
parsed header so ``align.check_alignment`` can then check the coordinates
against the genome's lengths under the genome's own names.
"""

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional

from readers import HeaderResult

# Spellings of the two organelle genomes, after the chr/chromosome prefix has
# been removed and the name lower-cased.
_MITO = {'m', 'mt', 'mito', 'mitochondria', 'mitochondrion', 'mitochondrial'}
_PLASTID = {'c', 'pt', 'cp', 'chloroplast', 'plastid', 'chloroplastic'}

_PREFIX_RE = re.compile(r'^(chromosome|chrom|chr)[\s_.\-]?', re.IGNORECASE)


def name_key(name: str) -> str:
    """
    The style-free key two spellings of one chromosome share.

    ``Chr1``, ``chr01``, ``chromosome_1`` and ``1`` all give ``'1'``; ``ChrM``,
    ``MT`` and ``mitochondria`` give ``'M'``; ``ChrC``, ``Pt`` and
    ``chloroplast`` give ``'C'``. Anything else is lower-cased with the prefix
    removed, so ``scaffold_12`` and ``Scaffold_12`` share a key too.
    """
    stripped = name.strip()
    rest = _PREFIX_RE.sub('', stripped)
    if rest == '':
        # The name *was* the prefix (a contig literally called "chr"): keep it.
        rest = stripped
    lowered = rest.lower()
    if lowered in _MITO:
        return 'M'
    if lowered in _PLASTID:
        return 'C'
    if lowered.isdigit():
        return str(int(lowered))
    return lowered


@dataclass
class NameMatch:
    """
    How a track's chromosome names map onto a genome's.

    ``renames`` maps each track name that needs changing to the genome name it
    means; ``unmatched`` lists the track names that are not in the genome under
    any spelling, in the order the track first used them; ``matched`` is every
    genome chromosome the track reaches (exactly or through a rename).
    """

    renames: Dict[str, str] = field(default_factory=dict)
    unmatched: List[str] = field(default_factory=list)
    matched: List[str] = field(default_factory=list)

    @property
    def all_unmatched(self) -> bool:
        """True when the track uses names and none of them are in the genome."""
        return bool(self.unmatched) and not self.matched

    def target(self, name: str) -> str:
        """The genome's name for a track name (itself when kept or unmatched)."""
        return self.renames.get(name, name)

    def report(self) -> Optional[dict]:
        """
        The member-facing account of the check, or None when nothing changed.

        ``renamed`` is the ``{track name: genome name}`` map, ``not_in_genome``
        the unmatched names, and ``messages`` both as sentences.
        """
        if not self.renames and not self.unmatched:
            return None
        return {
            'renamed': dict(self.renames),
            'not_in_genome': list(self.unmatched),
            'messages': self.messages(),
        }

    def messages(self) -> List[str]:
        out: List[str] = []
        if self.renames:
            pairs = ', '.join(f'{src} -> {dst}'
                              for src, dst in self.renames.items())
            out.append('Chromosome names were renamed to match the reference '
                       f'genome: {pairs}.')
        if self.unmatched:
            names = ', '.join(self.unmatched)
            plural = 's are' if len(self.unmatched) > 1 else ' is'
            out.append(f'Chromosome name{plural} not in the reference genome '
                       f'under any spelling: {names}.')
        return out


def match_names(track_names: Iterable[str],
                genome_names: Iterable[str]) -> NameMatch:
    """
    Decide, for every name a track uses, which genome chromosome it means.

    Exact (case-sensitive) matches win; then a name whose :func:`name_key` is
    the key of exactly one genome chromosome is renamed to it; everything else
    is unmatched. Order follows the track's first use of each name.
    """
    genome_list = list(genome_names)
    genome_set = set(genome_list)
    by_key: Dict[str, List[str]] = {}
    for gname in genome_list:
        by_key.setdefault(name_key(gname), []).append(gname)

    result = NameMatch()
    seen = set()
    matched_seen = set()
    for name in track_names:
        if name in seen:
            continue
        seen.add(name)
        if name in genome_set:
            target = name
        else:
            candidates = by_key.get(name_key(name), [])
            target = candidates[0] if len(candidates) == 1 else None
            if target is None:
                result.unmatched.append(name)
                continue
            result.renames[name] = target
        if target not in matched_seen:
            matched_seen.add(target)
            result.matched.append(target)
    return result


def rename_header(header: HeaderResult, match: NameMatch) -> HeaderResult:
    """
    The header under the genome's names, for the coordinate check.

    Unmatched names are dropped, so ``align.check_alignment`` on the result
    checks only the chromosomes the track will actually be drawn on. Two track
    names merged into one genome chromosome keep the larger extent.
    """
    if header.reason:
        return header
    lengths: Dict[str, int] = {}
    for name, length in header.chrom_lengths.items():
        if name in match.unmatched:
            continue
        target = match.target(name)
        lengths[target] = max(lengths.get(target, 0), length)
    extents: Dict[str, int] = {}
    for name, extent in header.chrom_extents.items():
        if name in match.unmatched:
            continue
        target = match.target(name)
        extents[target] = max(extents.get(target, 0), extent)
    return HeaderResult(chrom_lengths=lengths, chrom_extents=extents)


# ---------------------------------------------------------------------------
# Rewriting text formats line by line
# ---------------------------------------------------------------------------


def rename_feature_line(line: str, fmt: str, renames: Dict[str, str]) -> str:
    """
    One BED/GFF line with its chromosome column renamed (comments kept as-is).

    The chromosome is the first column of a feature line in both formats; a
    GFF3 ``##sequence-region <seqid> <start> <end>`` pragma names one too.
    ``line`` keeps its own line ending.
    """
    if not renames:
        return line
    if fmt == 'gff' and line.startswith('##sequence-region'):
        return _SEQ_REGION_RE.sub(
            lambda m: m.group(1) + renames.get(m.group(2), m.group(2)), line,
            count=1)
    stripped = line.strip()
    if (stripped == '' or stripped.startswith('#')
            or stripped.startswith('track ') or stripped.startswith('browser ')):
        return line
    if fmt == 'bedpe':
        # Two chromosome columns: the first and the fourth.
        body = line.rstrip('\r\n')
        ending = line[len(body):]
        tabbed = '\t' in body
        fields = body.split('\t') if tabbed else body.split()
        for i in (0, 3):
            if len(fields) > i and fields[i] in renames:
                fields[i] = renames[fields[i]]
        return ('\t' if tabbed else ' ').join(fields) + ending
    sep = '\t' if '\t' in line else None
    if sep:
        first, _, rest = line.partition('\t')
        if first in renames:
            return renames[first] + '\t' + rest
        return line
    # Whitespace-separated BED: replace the leading token only.
    match = re.match(r'^(\S+)', line)
    if match and match.group(1) in renames:
        return renames[match.group(1)] + line[match.end(1):]
    return line


_SEQ_REGION_RE = re.compile(r'^(##sequence-region\s+)(\S+)')
_WIG_CHROM_RE = re.compile(r'(\bchrom=)(\S+)')


def rename_wig_lines(lines: Iterable[str],
                     renames: Dict[str, str]) -> Iterator[str]:
    """
    A wig/bedGraph file's lines with each chromosome renamed.

    ``fixedStep``/``variableStep`` declare ``chrom=<name>`` and are followed by
    data lines that are only numbers; those pass through untouched even when a
    number happens to be a chromosome name being renamed (``1`` -> ``Chr1``), so
    the reader keeps track of which kind of block it is in, as ``scan_wig_chroms``
    does. A bedGraph line (outside any step block, four columns) starts with
    the name.
    """
    in_step_block = False
    for line in lines:
        stripped = line.strip()
        lower = stripped.lower()
        if lower.startswith('fixedstep') or lower.startswith('variablestep'):
            in_step_block = True
            yield _WIG_CHROM_RE.sub(
                lambda m: m.group(1) + renames.get(m.group(2), m.group(2)),
                line)
            continue
        if (stripped == '' or stripped.startswith('#')
                or lower.startswith('track') or lower.startswith('browser')):
            yield line
            continue
        fields = stripped.split()
        if in_step_block and len(fields) < 3:
            yield line
            continue
        # A bedGraph line, which also ends any step block before it.
        in_step_block = False
        match = re.match(r'^(\s*)(\S+)', line)
        if match and match.group(2) in renames:
            line = (match.group(1) + renames[match.group(2)]
                    + line[match.end(2):])
        yield line


def rewrite_text_file(src, dst, fmt: str, renames: Dict[str, str]) -> None:
    """
    Copy a BED/GFF/wig file to ``dst`` with its chromosome names renamed.

    Streams line by line, so a large annotation is never held in memory twice.
    Line endings are kept as they were.
    """
    with open(src, 'r', encoding='utf-8', errors='replace', newline='') as inp, \
            open(dst, 'w', encoding='utf-8', newline='') as out:
        if fmt == 'wig':
            for line in rename_wig_lines(inp, renames):
                out.write(line)
        else:
            for line in inp:
                out.write(rename_feature_line(line, fmt, renames))


def rewrite_bigwig(src: str, dst: str, renames: Dict[str, str]) -> None:
    """
    Write a copy of a BigWig whose header uses the genome's chromosome names.

    A BigWig's names live in its header, so renaming means writing a new file:
    the same chromosomes in the same order under their new names, then every
    interval copied across. pyBigWig is imported here so the rest of the module
    (and its tests) never needs the native library.
    """
    import pyBigWig  # lazily: only real BigWig rewrites need it

    source = pyBigWig.open(src)
    target = pyBigWig.open(dst, 'w')
    try:
        chroms = source.chroms() or {}
        header = [(renames.get(name, name), int(length))
                  for name, length in chroms.items()]
        if len({name for name, _ in header}) != len(header):
            raise ValueError('two of the BigWig\'s chromosomes would get the '
                             'same name in the reference genome')
        target.addHeader(header)
        for name in chroms:
            intervals = source.intervals(name) or []
            if not intervals:
                continue
            new_name = renames.get(name, name)
            target.addEntries(
                [new_name] * len(intervals),
                [iv[0] for iv in intervals],
                ends=[iv[1] for iv in intervals],
                values=[float(iv[2]) for iv in intervals])
    finally:
        target.close()
        source.close()


def rename_sam_header(header_text: str, renames: Dict[str, str]) -> str:
    """A SAM header with each ``@SQ SN:`` renamed (for ``samtools reheader``)."""
    out = []
    # Split on newlines only (keeping them), never on str.splitlines' extras.
    pieces = header_text.split('\n')
    lines = [piece + '\n' for piece in pieces[:-1]]
    if pieces[-1]:
        lines.append(pieces[-1])
    for line in lines:
        if line.startswith('@SQ'):
            parts = line.rstrip('\r\n').split('\t')
            ending = line[len(line.rstrip('\r\n')):]
            for i, token in enumerate(parts):
                if token.startswith('SN:') and token[3:] in renames:
                    parts[i] = 'SN:' + renames[token[3:]]
            line = '\t'.join(parts) + ending
        out.append(line)
    return ''.join(out)
