"""
Mean coverage over a region, counting the bases a track has no data for as zero.

A coverage track (BigWig) has no entries where nothing was read. pyBigWig's ``stats(type='mean')``,
and deepTools' ``multiBigwigSummary`` which uses it, average over the bases that do have entries:
a gene with 10 over half its length reads 10, not 5, and a stretch with no reads at all is NaN,
which ``plotCorrelation`` then drops, so the bins where replicates disagree most leave the
correlation. Here a region's mean is the signal summed over it divided by its length.

This file is also the script the recipes carry (it needs only pyBigWig, and numpy for ``bins``):

    python3 bigwig_means.py features regions.bed A.bw B.bw --labels A B -o raw.tab
    python3 bigwig_means.py bins A.bw B.bw --labels A B -bs 10000 -o summary.npz --outRawCounts bins.tab

Both write the tab-separated table ``multiBigwigSummary --outRawCounts`` does (0-based starts),
and ``bins`` the .npz it does, which ``plotCorrelation`` and ``plotPCA`` read.
"""
import argparse
import sys
from typing import List, Optional

# A big chromosome is read in pieces of this many bins, so memory stays small.
_CHUNK_BINS = 10_000


def mean(bw, chrom: str, start: int, end: int) -> Optional[float]:
    """The mean over [start, end), missing bases as zero; None when the track lacks the chromosome."""
    length = (bw.chroms() or {}).get(chrom)
    if length is None:
        return None
    end = min(end, length)
    if end <= start:
        return 0.0
    total = (bw.stats(chrom, start, end, type='sum', exact=True) or [None])[0]
    return (total or 0.0) / (end - start)


def bin_means(bw, chrom: str, start: int, end: int, size: int) -> List[float]:
    """Means over consecutive bins of ``size`` from start (the last one shorter), missing bases as zero."""
    out: List[float] = []
    length = (bw.chroms() or {}).get(chrom)
    if length is None:
        return out
    end = min(end, length)
    s = start
    while end - s >= size:
        n = min((end - s) // size, _CHUNK_BINS)
        sums = bw.stats(chrom, s, s + n * size, type='sum', nBins=n) or [None] * n
        out += [(v or 0.0) / size for v in sums]
        s += n * size
    if s < end:
        total = (bw.stats(chrom, s, end, type='sum') or [None])[0]
        out.append((total or 0.0) / (end - s))
    return out


def _header(labels: List[str]) -> str:
    return '\t'.join(["#'chr'", "'start'", "'end'"] + [f"'{x}'" for x in labels]) + '\n'


def _fmt(v: float) -> str:
    return f'{v:.6g}'


def features(bed: str, bigwigs: List[str], labels: List[str], out: str) -> int:
    """One row per line of the BED file: each track's mean over it. Returns the number of rows."""
    import pyBigWig
    handles = [pyBigWig.open(p) for p in bigwigs]
    n = 0
    try:
        with open(bed, encoding='utf-8') as fh, open(out, 'w', encoding='utf-8', newline='\n') as w:
            w.write(_header(labels))
            for line in fh:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                f = line.rstrip('\n').split('\t')
                chrom, s, e = f[0], int(f[1]), int(f[2])
                values = [mean(bw, chrom, s, e) for bw in handles]
                w.write('\t'.join([chrom, str(s), str(e)] + ['nan' if v is None else _fmt(v) for v in values]) + '\n')
                n += 1
    finally:
        for bw in handles:
            bw.close()
    return n


def bins(bigwigs: List[str], labels: List[str], size: int, npz: str, raw: Optional[str] = None) -> int:
    """The genome in bins of ``size``, on the chromosomes every track has; each track's mean per bin."""
    import numpy as np
    import pyBigWig
    handles = [pyBigWig.open(p) for p in bigwigs]
    try:
        sizes = [bw.chroms() or {} for bw in handles]
        chroms = [c for c in sizes[0] if all(c in s for s in sizes[1:])]
        if not chroms:
            raise SystemExit('The tracks have no chromosome name in common.')
        columns, rows = [[] for _ in handles], []
        for chrom in chroms:
            length = min(s[chrom] for s in sizes)
            per_track = [bin_means(bw, chrom, 0, length, size) for bw in handles]
            for col, values in zip(columns, per_track):
                col.extend(values)
            rows += [(chrom, s, min(s + size, length)) for s in range(0, length, size)]
        matrix = np.array(columns, dtype=float).T
    finally:
        for bw in handles:
            bw.close()
    np.savez_compressed(npz, matrix=matrix, labels=np.array(labels))
    if raw:
        with open(raw, 'w', encoding='utf-8', newline='\n') as w:
            w.write(_header(labels))
            for (chrom, s, e), values in zip(rows, matrix):
                w.write('\t'.join([chrom, str(s), str(e)] + [_fmt(v) for v in values]) + '\n')
    return len(rows)


def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0].strip())
    sub = p.add_subparsers(dest='what', required=True)
    f = sub.add_parser('features', help='each track’s mean over each line of a BED file')
    f.add_argument('bed')
    f.add_argument('bigwigs', nargs='+')
    f.add_argument('--labels', nargs='+', required=True)
    f.add_argument('-o', '--out', required=True)
    b = sub.add_parser('bins', help='each track’s mean in bins along the genome')
    b.add_argument('bigwigs', nargs='+')
    b.add_argument('--labels', nargs='+', required=True)
    b.add_argument('-bs', '--binSize', type=int, default=10000)
    b.add_argument('-o', '--out', required=True)
    b.add_argument('--outRawCounts')
    a = p.parse_args(argv)
    if len(a.labels) != len(a.bigwigs):
        p.error('give one label per BigWig file')
    if a.what == 'features':
        features(a.bed, a.bigwigs, a.labels, a.out)
    else:
        bins(a.bigwigs, a.labels, a.binSize, a.out, a.outRawCounts)


if __name__ == '__main__':
    main(sys.argv[1:])
