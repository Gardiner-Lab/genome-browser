"""
Marking a sequence on the genome: a feature track made by searching for it.

A guide RNA, a primer, a restriction site or any other stretch that the gene
annotation does not name is found on both strands of the reference genome and
kept as an ordinary BED feature track. From then on it is an ordinary track:
shown beside the others, shared, exported, tagged and removed the same way.

What a member types is one sequence per line, each optionally named:

    sgRNA1   GACGTTAGCACCTGTAGGAT
    sgRNA2   TTCGCAAGTGCCATAGGTCA
    GGTCTC                                 (unnamed: named after the track)

IUPAC ambiguity codes are allowed (N, R, Y, W, S, K, M, B, D, H, V), so a
degenerate primer or a PAM-bearing guide can be searched as written. The search
is exact otherwise: it finds what is there, and does not score mismatches, so a
guide's off-targets with a base or two changed are not reported.

The reference is read as a stream, in overlapping chunks, so a chromosome never
has to fit in memory and a match that straddles a chunk boundary is still found.
"""

import re
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Tuple

# Bases each IUPAC code may stand for (ambiguity codes expand to a character class).
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'R': 'Y', 'Y': 'R',
              'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D',
              'V': 'B', 'N': 'N'}

MIN_LENGTH = 8          # shorter than this is in every genome many times over
MAX_LENGTH = 1000
MAX_QUERIES = 50
MAX_PER_SEQUENCE = 500  # one sequence's own limit, so a common site cannot crowd out a guide
MAX_HITS = 5000         # the whole track
CHUNK = 4_000_000       # bases held at once while reading the reference


class Query(NamedTuple):
    name: str
    sequence: str


class Hit(NamedTuple):
    chrom: str
    start: int          # 0-based, as BED counts
    end: int
    name: str
    strand: str


def reverse_complement(seq: str) -> str:
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def parse_queries(text: str, fallback_name: str) -> List[Query]:
    """
    The lines a member typed, as named sequences.

    A line is ``name<space or tab>SEQUENCE`` or a bare sequence; blank lines and
    lines starting with ``#`` are skipped, and a FASTA ``>name`` header names the
    sequence under it. Raises ValueError, naming the line, when a sequence has a
    letter that is not a base or is too short or too long.
    """
    queries: List[Query] = []
    pending_name = ''
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('>'):
            pending_name = line[1:].strip()
            continue
        parts = line.replace('\t', ' ').split()
        if len(parts) >= 2:
            name, sequence = ' '.join(parts[:-1]), parts[-1]
        else:
            name, sequence = pending_name, parts[0]
        pending_name = ''
        sequence = sequence.upper().replace('-', '').replace('.', '')
        unknown = sorted(set(sequence) - set(IUPAC))
        if unknown:
            raise ValueError(f'Line {number} has {", ".join(unknown)}, which is not a base. '
                             'Use A, C, G, T and the IUPAC codes (N, R, Y, W, S, K, M, B, D, H, V).')
        if len(sequence) < MIN_LENGTH:
            raise ValueError(f'Line {number} is {len(sequence)} bases. A sequence of at least '
                             f'{MIN_LENGTH} is needed, or it matches all over the genome.')
        if len(sequence) > MAX_LENGTH:
            raise ValueError(f'Line {number} is {len(sequence):,} bases; at most {MAX_LENGTH:,} '
                             'can be searched for. Use a distinctive part of it.')
        queries.append(Query(name or fallback_name, sequence))
        if len(queries) > MAX_QUERIES:
            raise ValueError(f'At most {MAX_QUERIES} sequences at a time; make another track for the rest.')
    if not queries:
        raise ValueError('Type at least one sequence to look for.')
    return queries


def _pattern(sequence: str) -> 're.Pattern':
    """
    A regex matching ``sequence``, ambiguity codes included.

    Wrapped in a lookahead so overlapping matches are all found: a sequence that
    repeats head-to-tail is at every one of its positions, not every other one.
    """
    body = ''.join(b if len(IUPAC[b]) == 1 else f'[{IUPAC[b]}]' for b in sequence)
    return re.compile(f'(?=({body}))')


def _chunks(fasta: Path, overlap: int) -> Iterable[Tuple[str, str, int]]:
    """
    ``(chrom, bases, offset)`` over a FASTA, in chunks that overlap by ``overlap``
    bases so a match across a boundary is still whole in one of them.
    """
    chrom, held, offset = '', [], 0
    size = 0

    def flush(tail: bool):
        nonlocal held, size, offset
        if not chrom or not size:
            return None
        block = ''.join(held)
        out = (chrom, block, offset)
        if tail:
            held, size, offset = [], 0, 0
        else:
            keep = block[-overlap:] if overlap else ''
            offset += len(block) - len(keep)
            held, size = ([keep], len(keep)) if keep else ([], 0)
        return out

    with open(fasta, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('>'):
                done = flush(True)
                if done:
                    yield done
                chrom = line[1:].split()[0] if len(line) > 1 else ''
                held, size, offset = [], 0, 0
                continue
            if not chrom:
                continue
            held.append(line.strip().upper())
            size += len(held[-1])
            if size >= CHUNK:
                done = flush(False)
                if done:
                    yield done
    done = flush(True)
    if done:
        yield done


def search(fasta: Path, queries: List[Query], chroms: Dict[str, int]) -> Tuple[List[Hit], List[str]]:
    """
    Every place the sequences occur, on both strands: ``(hits, capped)``.

    A hit on the minus strand is where the reverse complement reads. A sequence
    that is its own reverse complement is reported once, on the plus strand.
    Chromosomes the genome does not list are skipped (the FASTA and the recorded
    sizes are the same file's, so this does not normally happen).

    Each sequence is allowed :data:`MAX_PER_SEQUENCE` hits of its own, and the
    search goes on for the others when one reaches it: otherwise a sequence that
    is all over the genome (an enzyme site, a short degenerate primer) would use
    up the whole track and a guide typed beside it would never be reached.
    ``capped`` names the sequences that stopped at their limit.
    """
    patterns = []
    for query in queries:
        rc = reverse_complement(query.sequence)
        patterns.append((query, '+', _pattern(query.sequence)))
        if rc != query.sequence:
            patterns.append((query, '-', _pattern(rc)))
    longest = max(len(q.sequence) for q in queries)

    found: Dict[Tuple[str, int, int, str], Hit] = {}
    counts: Dict[str, int] = {q.name: 0 for q in queries}
    capped: List[str] = []
    for chrom, block, offset in _chunks(fasta, longest - 1):
        if chroms and chrom not in chroms:
            continue
        if len(found) >= MAX_HITS:
            break
        for query, strand, pattern in patterns:
            if counts[query.name] >= MAX_PER_SEQUENCE:
                continue                # this one is full; the others go on
            for match in pattern.finditer(block):
                start = offset + match.start()
                key = (chrom, start, start + len(match.group(1)), query.name)
                if key in found:
                    continue            # the chunks' overlap, seen twice
                found[key] = Hit(chrom, start, start + len(match.group(1)), query.name, strand)
                counts[query.name] += 1
                if counts[query.name] >= MAX_PER_SEQUENCE:
                    if query.name not in capped:
                        capped.append(query.name)
                    break
            if len(found) >= MAX_HITS:
                break
    hits = sorted(found.values(), key=lambda h: (h.chrom, h.start, h.end, h.name))
    return hits, capped


def number_hits(hits: List[Hit]) -> List[Hit]:
    """Name each hit of a sequence found more than once ``name.1``, ``name.2``..."""
    totals: Dict[str, int] = {}
    for hit in hits:
        totals[hit.name] = totals.get(hit.name, 0) + 1
    seen: Dict[str, int] = {}
    out = []
    for hit in hits:
        if totals[hit.name] == 1:
            out.append(hit)
            continue
        seen[hit.name] = seen.get(hit.name, 0) + 1
        out.append(hit._replace(name=f'{hit.name}.{seen[hit.name]}'))
    return out


def bed_text(hits: List[Hit]) -> str:
    """The hits as BED: chrom, start, end, name, score, strand."""
    return ''.join(f'{h.chrom}\t{h.start}\t{h.end}\t{h.name}\t0\t{h.strand}\n' for h in hits)
