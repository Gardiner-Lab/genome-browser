"""
Mean coverage over a region, counting the bases a track has no data for as zero.

A coverage track (BigWig) has no entries where nothing was read. pyBigWig's ``stats(type='mean')``,
and deepTools' ``multiBigwigSummary`` which uses it, average over the bases that do have entries:
a gene with 10 over half its length reads 10, not 5, and a stretch with no reads at all is NaN,
which ``plotCorrelation`` then drops, so the bins where replicates disagree most leave the
correlation. Here a region's mean is the signal summed over it divided by its length.

Always read exactly (``exact=True``), never from the file's zoom levels: on real tracks the zoom-level
'sum' came back 0 for every bin, and the zoom-level 'mean' was some 10% off; exact costs about a
second per chromosome.

A methylation track (signal_kind.py) is measured as a *level* instead (``level=True``, or ``--as
level`` for its file): the mean over the bases that have data, the cytosines, and no value (NaN) where
there are none, since a stretch without cytosines is not unmethylated. pyBigWig's own exact 'mean'
is exactly that.

This file is also the script the recipes carry (it needs only pyBigWig, and numpy for ``bins``):

    python3 bigwig_means.py features regions.bed A.bw B.bw --labels A B -o raw.tab
    python3 bigwig_means.py bins A.bw B.bw --labels A B -bs 10000 -o summary.npz --outRawCounts bins.tab
    python3 bigwig_means.py features regions.bed RNA.bw CG.level.bw --labels RNA CG --as coverage level -o raw.tab
    python3 bigwig_means.py unsigned CG.bw -o CG.level.bw       # minus-strand cytosines' sign dropped
    python3 bigwig_means.py fill matrix.gz --samples 1 -o matrix.gz   # computeMatrix's NaN as 0, for sample 1
    python3 bigwig_means.py interpolate matrix.gz --samples 2 -o heatmap-matrix.gz   # for the heatmap only
    python3 bigwig_means.py heatmap -m matrix.gz -o heatmap.png ...   # plotHeatmap, column names at 45 degrees
    python3 bigwig_means.py heatmap -m heatmap-matrix.gz --profileMatrix matrix.gz ...   # profiles as measured

Both write the tab-separated table ``multiBigwigSummary --outRawCounts`` does (0-based starts),
and ``bins`` the .npz it does, which ``plotCorrelation`` and ``plotPCA`` read.
"""
import argparse
import sys
from typing import Callable, List, Optional

Progress = Optional[Callable[[float], None]]   # told the fraction done, as it goes

# A big chromosome is read in pieces of this many bins, so memory stays small.
_CHUNK_BINS = 10_000


def mean(bw, chrom: str, start: int, end: int, level: bool = False) -> Optional[float]:
    """The mean over [start, end), missing bases as zero; None when the track lacks the chromosome.

    ``level``: the mean over the bases with data only, and None where there are none."""
    length = (bw.chroms() or {}).get(chrom)
    if length is None:
        return None
    end = min(end, length)
    if end <= start:
        return None if level else 0.0
    if level:
        return (bw.stats(chrom, start, end, type='mean', exact=True) or [None])[0]
    total = (bw.stats(chrom, start, end, type='sum', exact=True) or [None])[0]
    return (total or 0.0) / (end - start)


def bin_means(bw, chrom: str, start: int, end: int, size: int, level: bool = False) -> List[float]:
    """Means over consecutive bins of ``size`` from start (the last one shorter), missing bases as zero.

    ``level``: each bin's mean over the bases with data only, and NaN for a bin with none."""
    out: List[float] = []
    length = (bw.chroms() or {}).get(chrom)
    if length is None:
        return out
    end = min(end, length)
    if level:
        nan = float('nan')
        s = start
        while s < end:
            n = min(-(-(end - s) // size), _CHUNK_BINS)
            stop = min(end, s + n * size)
            # Equal bins, but the last may be shorter: read the whole ones together, the last alone.
            whole = (stop - s) // size
            if whole:
                vals = bw.stats(chrom, s, s + whole * size, type='mean', nBins=whole, exact=True) or [None] * whole
                out += [nan if v is None else v for v in vals]
            if s + whole * size < stop:
                v = (bw.stats(chrom, s + whole * size, stop, type='mean', exact=True) or [None])[0]
                out.append(nan if v is None else v)
            s = stop
        return out
    s = start
    while end - s >= size:
        n = min((end - s) // size, _CHUNK_BINS)
        sums = bw.stats(chrom, s, s + n * size, type='sum', nBins=n, exact=True) or [None] * n
        out += [(v or 0.0) / size for v in sums]
        s += n * size
    if s < end:
        total = (bw.stats(chrom, s, end, type='sum', exact=True) or [None])[0]
        out.append((total or 0.0) / (end - s))
    return out


def _header(labels: List[str]) -> str:
    return '\t'.join(["#'chr'", "'start'", "'end'"] + [f"'{x}'" for x in labels]) + '\n'


def _fmt(v: float) -> str:
    return f'{v:.6g}'


def features(bed: str, bigwigs: List[str], labels: List[str], out: str, progress: Progress = None,
             levels: Optional[List[bool]] = None) -> int:
    """One row per line of the BED file: each track's mean over it. Returns the number of rows."""
    levels = levels or [False] * len(bigwigs)
    import pyBigWig
    with open(bed, encoding='utf-8') as fh:
        total = sum(1 for _ in fh) or 1
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
                values = [mean(bw, chrom, s, e, lv) for bw, lv in zip(handles, levels)]
                w.write('\t'.join([chrom, str(s), str(e)] + ['nan' if v is None else _fmt(v) for v in values]) + '\n')
                n += 1
                if progress and n % 500 == 0:
                    progress(n / total)
    finally:
        for bw in handles:
            bw.close()
    return n


def bins(bigwigs: List[str], labels: List[str], size: int, npz: str, raw: Optional[str] = None,
         progress: Progress = None, levels: Optional[List[bool]] = None) -> int:
    """The genome in bins of ``size``, on the chromosomes every track has; each track's mean per bin."""
    levels = levels or [False] * len(bigwigs)
    import numpy as np
    import pyBigWig
    handles = [pyBigWig.open(p) for p in bigwigs]
    try:
        sizes = [bw.chroms() or {} for bw in handles]
        chroms = [c for c in sizes[0] if all(c in s for s in sizes[1:])]
        if not chroms:
            raise SystemExit('The tracks have no chromosome name in common.')
        columns, rows = [[] for _ in handles], []
        lengths = {c: min(s[c] for s in sizes) for c in chroms}
        total, done = (sum(lengths.values()) * len(handles)) or 1, 0
        for chrom in chroms:
            length = lengths[chrom]
            per_track = []
            for bw, lv in zip(handles, levels):
                per_track.append(bin_means(bw, chrom, 0, length, size, lv))
                done += length
                if progress:
                    progress(done / total)
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


def unsigned(src: str, out: str) -> None:
    """A copy of a BigWig file with every value's sign dropped: a methylation track that gives
    minus-strand cytosines a minus sign, as its level on both strands."""
    import pyBigWig
    a = pyBigWig.open(str(src))
    b = pyBigWig.open(str(out), 'w')
    try:
        chroms = list((a.chroms() or {}).items())
        b.addHeader(chroms)
        for chrom, _ in chroms:
            entries = a.intervals(chrom) or []
            if entries:
                b.addEntries([chrom] * len(entries), [e[0] for e in entries],
                             ends=[e[1] for e in entries], values=[abs(e[2]) for e in entries])
    finally:
        a.close()
        b.close()


def fill(matrix: str, samples: List[int], out: str) -> None:
    """
    In a computeMatrix matrix (.gz), the missing values (NaN) of the given samples (1 for the first
    -S file) set to 0: what --missingDataAsZero does, but for some of the tracks only, so coverage
    tracks can count stretches without reads as zero while a methylation track leaves out the
    stretches without cytosines.
    """
    import gzip
    import json
    with gzip.open(matrix, 'rt', encoding='utf-8') as fh:
        lines = fh.read().split('\n')
    header = json.loads(lines[0][1:])
    bounds = header['sample_boundaries']
    cols = set()
    for i in samples:
        cols.update(range(bounds[i - 1], bounds[i]))
    body = []
    for line in lines[1:]:
        if not line:
            continue
        f = line.split('\t')
        values = f[6:]
        for c in cols:
            if values[c] in ('nan', 'NaN', ''):
                values[c] = '0'
        body.append('\t'.join(f[:6] + values))
    with gzip.open(out, 'wt', encoding='utf-8') as fh:
        fh.write(lines[0] + '\n' + ''.join(line + '\n' for line in body))


def interpolate(matrix: str, samples: List[int], out: str) -> None:
    """
    For drawing a heatmap of a methylation track: in each feature's row, a step with no cytosine
    measured (NaN) takes its value from the nearest measured steps on either side (linearly between
    them; the nearest one past the last). A row with nothing measured stays NaN.

    Without it a sparse track's heatmap is all one colour: nearly every row has some empty steps, and
    plotHeatmap, squeezing thousands of rows into each pixel, paints any pixel that touches one as
    missing. The mean profile is drawn from the matrix as measured, never from this.
    """
    import gzip
    import json
    import numpy as np
    with gzip.open(matrix, 'rt', encoding='utf-8') as fh:
        lines = fh.read().split('\n')
    bounds = json.loads(lines[0][1:])['sample_boundaries']
    spans = [(bounds[i - 1], bounds[i]) for i in samples]
    body = []
    for line in lines[1:]:
        if not line:
            continue
        f = line.split('\t')
        values = f[6:]
        for a, b in spans:
            seg = np.array([float(v) if v not in ('', 'nan', 'NaN') else np.nan for v in values[a:b]])
            ok = ~np.isnan(seg)
            if ok.any() and not ok.all():
                at = np.arange(len(seg))
                seg[~ok] = np.interp(at[~ok], at[ok], seg[ok])
                values[a:b] = [f'{v:.6g}' for v in seg]
        body.append('\t'.join(f[:6] + values))
    with gzip.open(out, 'wt', encoding='utf-8') as fh:
        fh.write(lines[0] + '\n' + ''.join(line + '\n' for line in body))


def heatmap(args: List[str]) -> None:
    """
    deepTools' plotHeatmap with the same arguments, but each column's name (the track's) tilted 45
    degrees, starting over the middle of its column, and the plot title moved above them. plotHeatmap
    writes them level and centred, so two long track names run into each other, and it has no option
    to turn them.

    ``--profileMatrix FILE`` (not plotHeatmap's own): the profiles drawn over the heatmap's columns
    are averaged from that matrix instead of the one the heatmap shows. With ``-m heatmap-matrix.gz``
    (``interpolate``) and ``--profileMatrix matrix.gz`` the heatmap has its gaps filled while the
    profiles are as measured. Each feature's row is found in it by its place, name and strand, so the
    profiles follow the heatmap's sorting and clusters.
    """
    profile_matrix = None
    if '--profileMatrix' in args:
        at = args.index('--profileMatrix')
        profile_matrix = args[at + 1]
        args = args[:at] + args[at + 2:]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.axes
    import matplotlib.pyplot as plt
    from deeptools import plotHeatmap

    set_title = matplotlib.axes.Axes.set_title

    def tilted(ax, label, *a, **kw):
        kw.update(rotation=45, rotation_mode='anchor', ha='left', va='bottom')
        text = set_title(ax, label, *a, **kw)
        text.set_x(0.5)
        return text

    savefig = plt.savefig

    def save(*a, **kw):
        fig = plt.gcf()
        top = fig._suptitle
        if top is not None and top.get_text():
            renderer = fig.canvas.get_renderer()
            heights = [ax.title.get_window_extent(renderer).y1 for ax in fig.axes if ax.title.get_text()]
            if heights:
                y = fig.transFigure.inverted().transform((0, max(heights)))[1]
                top.set_y(y + 0.01)
                top.set_va('bottom')
        return savefig(*a, **kw)

    add_profile = plotHeatmap.addProfilePlot

    def measured_profile(hm, *a, **kw):
        import numpy as np
        from deeptools import heatmapper
        source = heatmapper.heatmapper()
        source.read_matrix_file(profile_matrix)
        rows: dict = {}
        for region, row in zip(source.matrix.regions, source.matrix.matrix):
            rows.setdefault(_region_key(region), []).append(row)
        shown = hm.matrix.matrix
        try:
            hm.matrix.matrix = np.ma.vstack([rows[_region_key(r)].pop(0) for r in hm.matrix.regions])
        except (KeyError, IndexError):
            raise SystemExit(f'{profile_matrix} does not hold the same features as the heatmap')
        try:
            return add_profile(hm, *a, **kw)
        finally:
            hm.matrix.matrix = shown

    matplotlib.axes.Axes.set_title = tilted
    plt.savefig = save
    if profile_matrix:
        plotHeatmap.addProfilePlot = measured_profile
    try:
        plotHeatmap.main(args)
    finally:
        matplotlib.axes.Axes.set_title = set_title
        plt.savefig = savefig
        plotHeatmap.addProfilePlot = add_profile


def _region_key(region) -> tuple:
    """A feature in a deepTools matrix: chromosome, its blocks, name and strand."""
    chrom, blocks, name, _group, strand = region[:5]
    return chrom, tuple(tuple(b) for b in blocks), name, strand


def main(argv: Optional[List[str]] = None, progress: Progress = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    if argv[:1] == ['heatmap']:
        # Its arguments are plotHeatmap's own, passed on as they are.
        heatmap(argv[1:])
        return
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0].strip())
    sub = p.add_subparsers(dest='what', required=True)
    f = sub.add_parser('features', help='each track’s mean over each line of a BED file')
    f.add_argument('bed')
    f.add_argument('bigwigs', nargs='+')
    f.add_argument('--labels', nargs='+', required=True)
    f.add_argument('-o', '--out', required=True)
    f.add_argument('--as', dest='measure', nargs='+', choices=('coverage', 'level'),
                   help='per file: coverage (missing bases are zero, the default) or level (left out)')
    b = sub.add_parser('bins', help='each track’s mean in bins along the genome')
    b.add_argument('bigwigs', nargs='+')
    b.add_argument('--labels', nargs='+', required=True)
    b.add_argument('-bs', '--binSize', type=int, default=10000)
    b.add_argument('-o', '--out', required=True)
    b.add_argument('--outRawCounts')
    b.add_argument('--as', dest='measure', nargs='+', choices=('coverage', 'level'),
                   help='per file: coverage (missing bases are zero, the default) or level (left out)')
    u = sub.add_parser('unsigned', help='a copy of a BigWig file with every value made positive')
    u.add_argument('bigwig')
    u.add_argument('-o', '--out', required=True)
    m = sub.add_parser('fill', help='missing values of some samples of a computeMatrix matrix set to 0')
    m.add_argument('matrix')
    m.add_argument('--samples', nargs='+', type=int, required=True, help='1 for the first -S file')
    m.add_argument('-o', '--out', required=True)
    i = sub.add_parser('interpolate', help='for a heatmap: empty steps of some samples filled from their neighbours')
    i.add_argument('matrix')
    i.add_argument('--samples', nargs='+', type=int, required=True, help='1 for the first -S file')
    i.add_argument('-o', '--out', required=True)
    a = p.parse_args(argv)
    if a.what == 'interpolate':
        interpolate(a.matrix, a.samples, a.out)
        return
    if a.what == 'unsigned':
        unsigned(a.bigwig, a.out)
        return
    if a.what == 'fill':
        fill(a.matrix, a.samples, a.out)
        return
    if len(a.labels) != len(a.bigwigs):
        p.error('give one label per BigWig file')
    if a.measure and len(a.measure) != len(a.bigwigs):
        p.error('give one --as per BigWig file')
    levels = [m == 'level' for m in a.measure] if a.measure else None
    if a.what == 'features':
        features(a.bed, a.bigwigs, a.labels, a.out, progress, levels)
    else:
        bins(a.bigwigs, a.labels, a.binSize, a.out, a.outRawCounts, progress, levels)


if __name__ == '__main__':
    main(sys.argv[1:])
