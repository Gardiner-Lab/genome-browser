"""
Figures the Genome Browser draws itself, where deepTools has no tool for a coverage (BigWig) track.

This file is also a script the recipes carry (it needs numpy, scipy and matplotlib, which come with
deepTools):

    python3 plots.py fingerprint summary.npz -o fingerprint.svg --metrics fingerprint-metrics.tab
    python3 plots.py methexpr --genes regions.bed --expression expression.tab --methylation methylation.tab \\
        -o methexpr.svg --table methylation-vs-expression.csv --groups 5 --group-beds group --ranked ranked
    python3 plots.py genecorr --genes regions.bed --expression expression.tab --methylation methylation.tab \\
        -o genecorr.svg --table per-gene-correlation.csv

``fingerprint``: deepTools' plotFingerprint reads BAM files only. This draws the same curve from the
bins ``bigwig_means.py bins`` makes of each track: the bins sorted by their signal, lowest first, and
the running share of the track's whole signal. A track spread evenly over the genome follows the
diagonal; one whose signal sits in a few places stays low and rises steeply at the right.

``methexpr``: each feature's (gene's) expression, from an RNA-seq coverage track's mean over it, and
each methylation track's level over it (or over the stretch before it), from two
``bigwig_means.py features`` tables of the same BED file. The features are put in expression groups:
those with no reads, then the rest in equal-sized groups by expression. For each methylation track it
draws the level per group (boxes), and the level against expression (each feature a point, with
Spearman's rank correlation). It writes the table, and each group as a BED file, for a profile of
each group (computeMatrix, then ``plots.py profile``: plotProfile with its legend below the panels). With ``--ranked``, a file per methylation track ranks the genes by how
much more methylated than expressed they are (percentile ranks): methylated, silent genes first.

``genecorr``: several samples, each a methylation track and an RNA-seq track of the same material
(WT and mutants, tissues, time points), in the same order in both tables. For each gene, Spearman's
correlation between its methylation level and its expression across the samples, with an exact
p-value (all orderings of the samples, up to nine) and Benjamini-Hochberg q-values; the genes ranked
from the most negative to the most positive.
"""
import argparse
import csv
import json
import sys
from typing import List, Optional


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'font.family': 'sans-serif'})
    return plt


# --- deepTools' plotProfile, with its legend below ---------------------------------------------------

def profile(args: List[str]) -> None:
    """
    deepTools' plotProfile with the same arguments, but the legend moved out of the panels into a row
    below them: with several lines (an expression group each) it sat on top of the curves, and
    plotProfile can only put it inside.
    """
    plt = _plt()
    from deeptools import plotProfile
    savefig = plt.savefig

    def save(*a, **kw):
        fig = plt.gcf()
        handles, labels = [], []
        for ax in fig.axes:
            legend = ax.get_legend()
            if legend is not None:
                if not handles:
                    handles, labels = ax.get_legend_handles_labels()
                legend.remove()
        if handles:
            fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.0),
                       ncol=min(len(labels), 3), frameon=False, fontsize=9)
        kw['bbox_inches'] = 'tight'
        return savefig(*a, **kw)

    plt.savefig = save
    try:
        plotProfile.main(args)
    finally:
        plt.savefig = savefig


# --- fingerprint ------------------------------------------------------------------------------

def fingerprint_curve(values):
    """The curve for one track's bins: x the share of bins (sorted, lowest first), y the share of signal."""
    import numpy as np
    v = np.sort(np.nan_to_num(np.asarray(values, dtype=float), nan=0.0))
    v = np.clip(v, 0, None)
    n = len(v)
    x = np.arange(1, n + 1) / n
    total = v.sum()
    y = np.cumsum(v) / total if total > 0 else np.zeros(n)
    return x, y


def fingerprint_metrics(values) -> dict:
    """
    What the curve says, in numbers. ``zero_bins``: the share of bins with no signal (deepTools'
    X-intercept). ``top1``: the share of all signal in the 1% of bins with most. ``elbow``: where the
    curve is furthest below the diagonal, as the share of bins left of it (deepTools' Elbow Point).
    ``auc``: the area under the curve, 0.5 for signal spread evenly and less the more it is gathered.
    """
    import numpy as np
    v = np.clip(np.nan_to_num(np.asarray(values, dtype=float), nan=0.0), 0, None)
    n = len(v)
    if n == 0 or v.sum() <= 0:
        return {'bins': n, 'zero_bins': 1.0 if n else None, 'top1': None, 'elbow': None, 'auc': None}
    x, y = fingerprint_curve(v)
    top = max(1, int(round(n * 0.01)))
    return {'bins': int(n),
            'zero_bins': float((v == 0).mean()),
            'top1': float(np.sort(v)[-top:].sum() / v.sum()),
            'elbow': float(x[int(np.argmax(x - y))]),
            'auc': float(np.trapz(np.concatenate([[0.0], y]), np.concatenate([[0.0], x])))}


def fingerprint(npz: str, out: str, metrics: Optional[str] = None, title: str = '') -> List[dict]:
    import numpy as np
    plt = _plt()
    data = np.load(npz)
    matrix, labels = data['matrix'], [str(x) for x in data['labels']]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.plot([0, 1], [0, 1], ls='--', color='#999', lw=1, label='signal spread evenly')
    found = []
    for i, label in enumerate(labels):
        x, y = fingerprint_curve(matrix[:, i])
        # Every bin is too many points for an SVG: a thousand steps along the curve draw the same line.
        keep = np.unique(np.linspace(0, len(x) - 1, min(len(x), 1000)).astype(int))
        ax.plot(np.concatenate([[0], x[keep]]), np.concatenate([[0], y[keep]]), lw=1.6, label=label)
        found.append({'track': label, **fingerprint_metrics(matrix[:, i])})
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('bins, sorted by signal (share of the genome)')
    ax.set_ylabel('share of the track’s whole signal')
    ax.set_title(title or 'Fingerprint')
    ax.legend(loc='upper left', fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    if metrics:
        with open(metrics, 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['track', 'bins', 'share_of_bins_without_signal', 'share_of_signal_in_top_1pct_bins',
                        'elbow_point', 'area_under_curve'])
            for m in found:
                w.writerow([m['track'], m['bins']] + ['' if m[k] is None else f'{m[k]:.4f}'
                                                       for k in ('zero_bins', 'top1', 'elbow', 'auc')])
    return found


# --- methylation against expression ---------------------------------------------------------------

def _table(path: str):
    """A bigwig_means.py features table: its labels, and each row's values (None for no data)."""
    labels, rows = [], []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                labels = [h.strip("'") for h in f[3:]]
                continue
            rows.append([None if v in ('nan', '') else float(v) for v in f[3:]])
    return labels, rows


def group_labels(groups: int) -> List[str]:
    """The expression groups' names: no reads, then equal shares from lowest to highest."""
    step = 100 / groups
    out = ['no reads']
    for g in range(groups):
        lo, hi = round(g * step), round((g + 1) * step)
        out.append('lowest {}%'.format(hi) if g == 0 else 'highest {}%'.format(100 - lo) if g == groups - 1
                   else '{}–{}%'.format(lo, hi))
    return out


def expression_groups(expression: List[Optional[float]], groups: int) -> List[Optional[int]]:
    """
    Each feature's group: 0 for no reads (an expression of 0), then 1 (lowest) to ``groups`` (highest),
    the expressed features split into equal-sized groups by expression (ties kept in order). None for a
    feature the track has no value for (its chromosome is not in the track).
    """
    out: List[Optional[int]] = [None] * len(expression)
    expressed = sorted((v, i) for i, v in enumerate(expression) if v is not None and v > 0)
    for i, v in enumerate(expression):
        if v is not None and v <= 0:
            out[i] = 0
    n = len(expressed)
    for rank, (_, i) in enumerate(expressed):
        out[i] = 1 + rank * groups // n
    return out


def _spearman(x, y):
    from scipy.stats import spearmanr
    if len(x) < 3:
        return None
    rho = spearmanr(x, y).correlation
    return None if rho != rho else float(rho)


def percentiles(values: List[Optional[float]]) -> List[Optional[float]]:
    """Each value's rank among those given (ties share the mean rank), from 0 (lowest) to 100 (highest)."""
    from scipy.stats import rankdata
    have = [i for i, v in enumerate(values) if v is not None]
    out: List[Optional[float]] = [None] * len(values)
    if not have:
        return out
    ranks = rankdata([values[i] for i in have])
    span = max(1, len(have) - 1)
    for i, r in zip(have, ranks):
        out[i] = (r - 1) / span * 100
    return out


def write_ranked(path: str, bed, expr, level, group, names) -> int:
    """
    The genes of one methylation track, ranked by how much more methylated than expressed they are:
    its methylation's percentile among the genes minus its expression's. At the top, methylated genes
    with little or no expression; at the foot, unmethylated genes that are highly expressed; genes as
    methylated as they are expressed (either both high or both low) in the middle.
    """
    e_pct, m_pct = percentiles(expr), percentiles(level)
    rows = [(m - e, i) for i, (e, m) in enumerate(zip(e_pct, m_pct)) if e is not None and m is not None]
    rows.sort(key=lambda r: (-r[0], r[1]))
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['rank', 'name', 'chrom', 'start', 'end', 'strand', 'expression', 'expression_percentile',
                    'methylation', 'methylation_percentile', 'methylated_minus_expressed', 'expression_group'])
        for rank, (score, i) in enumerate(rows, 1):
            f = bed[i]
            w.writerow([rank, f[3], f[0], int(f[1]) + 1, f[2], f[5] if len(f) > 5 else '.', f'{expr[i]:.6g}',
                        f'{e_pct[i]:.1f}', f'{level[i]:.6g}', f'{m_pct[i]:.1f}', f'{score:.1f}',
                        '' if group[i] is None else names[group[i]]])
    return len(rows)


def methexpr(genes: str, expression: str, methylation: str, out: str, table: str, groups: int = 5,
             group_beds: Optional[str] = None, title: str = '', where: str = 'over each gene',
             ranked: Optional[str] = None) -> dict:
    import numpy as np
    plt = _plt()
    _, expr_rows = _table(expression)
    meth_labels, meth_rows = _table(methylation)
    bed = [line.rstrip('\n').split('\t') for line in open(genes, encoding='utf-8') if line.strip()]
    if not (len(bed) == len(expr_rows) == len(meth_rows)):
        raise SystemExit('The tables do not have a row for each feature of the BED file.')
    expr = [r[0] for r in expr_rows]
    group = expression_groups(expr, groups)
    names = group_labels(groups)
    counts = [sum(1 for g in group if g == k) for k in range(groups + 1)]

    with open(table, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'chrom', 'start', 'end', 'strand', 'expression'] + meth_labels + ['expression_group'])
        for f, e, m, g in zip(bed, expr, meth_rows, group):
            w.writerow([f[3], f[0], int(f[1]) + 1, f[2], f[5] if len(f) > 5 else '.',
                        '' if e is None else f'{e:.6g}'] + ['' if v is None else f'{v:.6g}' for v in m] +
                       ['' if g is None else names[g]])
    beds = []
    if group_beds:
        for k in range(groups + 1):
            if not counts[k]:
                continue
            path = f'{group_beds}-{k}.bed'
            with open(path, 'w', encoding='utf-8', newline='\n') as fh:
                fh.writelines('\t'.join(f) + '\n' for f, g in zip(bed, group) if g == k)
            beds.append({'file': path, 'label': f'{names[k]} ({counts[k]:,})', 'count': counts[k]})

    ranked_files = []
    if ranked:
        for t, label in enumerate(meth_labels):
            path = f'{ranked}-{t + 1}.csv'
            n = write_ranked(path, bed, expr, [r[t] for r in meth_rows], group, names)
            ranked_files.append({'track': label, 'file': path, 'genes': n})

    tracks = []
    fig, axes = plt.subplots(len(meth_labels), 2, figsize=(11, 3.9 * len(meth_labels)), squeeze=False,
                             gridspec_kw={'width_ratios': [1.25, 1]})
    for t, label in enumerate(meth_labels):
        level = [r[t] for r in meth_rows]
        box, ax = axes[t]
        per_group = [[lv for lv, g in zip(level, group) if g == k and lv is not None] for k in range(groups + 1)]
        shown = [k for k in range(groups + 1) if per_group[k]]
        if shown:
            parts = box.boxplot([per_group[k] for k in shown], positions=range(len(shown)), widths=0.6,
                                showfliers=False, patch_artist=True, medianprops={'color': '#222'})
            for k, patch in zip(shown, parts['boxes']):
                patch.set_facecolor('#bbbbbb' if k == 0 else plt.cm.viridis(0.15 + 0.75 * (k - 1) / max(1, groups - 1)))
            box.plot(range(len(shown)), [float(np.mean(per_group[k])) for k in shown], 'o', color='#d62728',
                     ms=4, label='mean')
            box.set_xticks(range(len(shown)))
            box.set_xticklabels([f'{names[k]}\n(n={len(per_group[k]):,})' for k in shown], fontsize=8)
            box.legend(loc='upper right', fontsize=8, frameon=False)
        box.set_ylabel(f'{label}\nmethylation level')
        box.set_xlabel('expression group')
        box.set_title(f'{label}: methylation {where}, by expression', fontsize=10)

        pairs = [(e, lv) for e, lv in zip(expr, level) if e is not None and lv is not None]
        rho = _spearman([p[0] for p in pairs], [p[1] for p in pairs]) if pairs else None
        expressed = [(e, lv) for e, lv in pairs if e > 0]
        rho_expressed = _spearman([p[0] for p in expressed], [p[1] for p in expressed]) if expressed else None
        if pairs:
            xs = np.log2(np.array([p[0] for p in pairs]) + 1)
            ys = np.array([p[1] for p in pairs])
            hb = ax.hexbin(xs, ys, gridsize=45, bins='log', cmap='viridis', mincnt=1)
            fig.colorbar(hb, ax=ax, label='features (log scale)')
        ax.set_xlabel('expression, log2(mean coverage + 1)')
        ax.set_ylabel('methylation level')
        ax.set_title('Spearman ρ = ' + ('n/a' if rho is None else f'{rho:.2f}') + f' (n={len(pairs):,})'
                     + ('' if rho_expressed is None else f'; expressed only {rho_expressed:.2f}'), fontsize=10)
        tracks.append({'track': label, 'rho': rho, 'rho_expressed': rho_expressed, 'n': len(pairs),
                       'mean_by_group': [None if not per_group[k] else float(np.mean(per_group[k]))
                                         for k in range(groups + 1)]})
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return {'groups': [{'label': names[k], 'count': counts[k]} for k in range(groups + 1)],
            'no_expression_value': sum(1 for g in group if g is None),
            'tracks': tracks, 'beds': beds, 'ranked': ranked_files}


# --- per gene, across samples ---------------------------------------------------------------------

_NULL: dict = {}


def _null_rho(n: int):
    """Spearman's rho of every ordering of n samples against 1..n, without ties: its exact null."""
    import itertools
    import numpy as np
    if n not in _NULL:
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int16)
        d = perms - np.arange(n)
        _NULL[n] = np.sort(np.abs(1 - 6 * (d * d).sum(axis=1) / (n * (n * n - 1))))
    return _NULL[n]


def gene_correlation(meth: List[Optional[float]], expr: List[Optional[float]]):
    """
    One gene's Spearman correlation across the samples where it has both values: (rho, p, n). p is
    two-sided: exact over every ordering of the samples for up to nine, else from scipy's t
    approximation. None for rho where either value does not vary (it cannot rise or fall together).
    """
    import numpy as np
    from scipy.stats import spearmanr
    pairs = [(m, e) for m, e in zip(meth, expr) if m is not None and e is not None]
    n = len(pairs)
    if n < 3:
        return None, None, n
    m = [p[0] for p in pairs]
    e = [p[1] for p in pairs]
    if len(set(m)) < 2 or len(set(e)) < 2:
        return None, None, n
    res = spearmanr(m, e)
    rho = float(res.correlation)
    if rho != rho:
        return None, None, n
    if n <= 9:
        null = _null_rho(n)
        # Ties make rho take values between the null's; a hair of slack keeps equal ones equal.
        p = float((len(null) - np.searchsorted(null, abs(rho) - 1e-9)) / len(null))
    else:
        p = float(res.pvalue)
    return rho, p, n


def bh(pvalues: List[Optional[float]]) -> List[Optional[float]]:
    """Benjamini-Hochberg q-values over the p-values given (None left out)."""
    have = sorted((p, i) for i, p in enumerate(pvalues) if p is not None)
    q: List[Optional[float]] = [None] * len(pvalues)
    m, running = len(have), 1.0
    for k in range(m - 1, -1, -1):
        p, i = have[k]
        running = min(running, p * m / (k + 1))
        q[i] = running
    return q


def genecorr(genes: str, expression: str, methylation: str, out: str, table: str, title: str = '',
             where: str = 'over each gene') -> dict:
    import math
    import numpy as np
    plt = _plt()
    labels, expr_rows = _table(expression)
    meth_labels, meth_rows = _table(methylation)
    bed = [line.rstrip('\n').split('\t') for line in open(genes, encoding='utf-8') if line.strip()]
    if not (len(bed) == len(expr_rows) == len(meth_rows)):
        raise SystemExit('The tables do not have a row for each feature of the BED file.')
    if len(labels) != len(meth_labels):
        raise SystemExit('Give as many expression tracks as methylation tracks, sample by sample.')
    samples = meth_labels
    found = [gene_correlation(m, e) for m, e in zip(meth_rows, expr_rows)]
    q = bh([p for _, p, _ in found])
    order = sorted(range(len(bed)), key=lambda i: (found[i][0] is None, found[i][0] if found[i][0] is not None else 0,
                                                  found[i][1] if found[i][1] is not None else 1))
    with open(table, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['rank', 'name', 'chrom', 'start', 'end', 'strand', 'samples', 'spearman_rho', 'p', 'q']
                   + [f'methylation {s}' for s in samples] + [f'expression {s}' for s in samples])
        for rank, i in enumerate(order, 1):
            f, (rho, p, n) = bed[i], found[i]
            w.writerow([rank if rho is not None else '', f[3], f[0], int(f[1]) + 1, f[2], f[5] if len(f) > 5 else '.', n,
                        '' if rho is None else f'{rho:.4f}', '' if p is None else f'{p:.4g}',
                        '' if q[i] is None else f'{q[i]:.4g}']
                       + ['' if v is None else f'{v:.6g}' for v in meth_rows[i]]
                       + ['' if v is None else f'{v:.6g}' for v in expr_rows[i]])

    rhos = np.array([r for r, _, _ in found if r is not None])
    tested = len(rhos)
    best_p = min((p for _, p, _ in found if p is not None), default=None)
    negative = [i for i in order if found[i][0] is not None and found[i][0] < 0][:3]
    positive = [i for i in reversed(order) if found[i][0] is not None and found[i][0] > 0][:3]
    shown = negative[:2] + positive[:2]          # four small plots: the two most negative, the two most positive
    fig = plt.figure(figsize=(12, 7.5))
    grid = fig.add_gridspec(2, 4)
    ax = fig.add_subplot(grid[:, 0:2] if not shown else grid[0, 0:2])
    if tested:
        ax.hist(rhos, bins=np.linspace(-1, 1, 2 * len(samples) + 1 if len(samples) <= 9 else 41),
                color='#3b528b', edgecolor='white')
    ax.set_xlim(-1, 1)
    ax.set_xlabel('Spearman ρ across the samples (methylation against expression)')
    ax.set_ylabel('genes')
    ax.set_title(f'{tested:,} genes that vary in both, over {len(samples)} samples', fontsize=10)
    side = fig.add_subplot(grid[1, 0:2]) if shown else None
    if side is not None:
        side.axis('off')
        lines = ['Most negative: ' + ', '.join(bed[i][3] for i in negative),
                 'Most positive: ' + ', '.join(bed[i][3] for i in positive)]
        if best_p is not None:
            lines.append(f'Smallest p possible with {len(samples)} samples: {2 / math.factorial(len(samples)):.3g}'
                         if len(samples) <= 9 else '')
        side.text(0, 0.9, '\n'.join(x for x in lines if x), va='top', fontsize=9, wrap=True)
    cells = [grid[r, c] for r in range(2) for c in range(2, 4)]
    for cell, i in zip(cells + [None] * 2, shown[:4]):
        if cell is None:
            break
        a = fig.add_subplot(cell)
        m, e = meth_rows[i], expr_rows[i]
        xs = [math.log2(v + 1) if v is not None else None for v in e]
        for s, x, y in zip(samples, xs, m):
            if x is not None and y is not None:
                a.plot(x, y, 'o', color='#21918c')
                a.annotate(s, (x, y), fontsize=7, xytext=(3, 3), textcoords='offset points')
        rho, p, n = found[i]
        a.set_title(f'{bed[i][3]}: ρ {rho:.2f}, p {p:.2g}', fontsize=9)
        a.set_xlabel('log2(expression + 1)', fontsize=8)
        a.set_ylabel(f'methylation {where}', fontsize=8)
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return {'samples': samples, 'genes': len(bed), 'tested': tested,
            'significant': sum(1 for v in q if v is not None and v < 0.05),
            'smallest_p_possible': 2 / math.factorial(len(samples)) if len(samples) <= 9 else None,
            'most_negative': [{'name': bed[i][3], 'rho': found[i][0], 'p': found[i][1], 'q': q[i]} for i in order[:10]
                              if found[i][0] is not None],
            'most_positive': [{'name': bed[i][3], 'rho': found[i][0], 'p': found[i][1], 'q': q[i]}
                              for i in reversed(order) if found[i][0] is not None][:10]}


def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0].strip())
    sub = p.add_subparsers(dest='what', required=True)
    f = sub.add_parser('fingerprint', help='the fingerprint of the tracks in a bigwig_means.py bins .npz')
    f.add_argument('npz')
    f.add_argument('-o', '--out', required=True)
    f.add_argument('--metrics')
    f.add_argument('--title', default='')
    f.add_argument('--summary', help='also write what was found as JSON')
    m = sub.add_parser('methexpr', help='methylation against expression, per feature')
    m.add_argument('--genes', required=True, help='the BED file both tables were measured over (or made from)')
    m.add_argument('--expression', required=True, help='bigwig_means.py features table of the expression track')
    m.add_argument('--methylation', required=True, help='bigwig_means.py features table of the methylation tracks')
    m.add_argument('-o', '--out', required=True)
    m.add_argument('--table', required=True)
    m.add_argument('--groups', type=int, default=5)
    m.add_argument('--group-beds', help='write each expression group as PREFIX-N.bed')
    m.add_argument('--title', default='')
    m.add_argument('--where', default='over each gene')
    m.add_argument('--summary', help='also write what was found as JSON')
    m.add_argument('--ranked', help='write each methylation track\'s genes ranked as PREFIX-N.csv')
    g = sub.add_parser('genecorr', help='per gene, methylation against expression across samples')
    g.add_argument('--genes', required=True)
    g.add_argument('--expression', required=True, help='bigwig_means.py features table, a column per sample')
    g.add_argument('--methylation', required=True, help='the same, the methylation tracks in the same order')
    g.add_argument('-o', '--out', required=True)
    g.add_argument('--table', required=True)
    g.add_argument('--title', default='')
    g.add_argument('--where', default='over each gene')
    g.add_argument('--summary', help='also write what was found as JSON')
    if argv is None:
        argv = sys.argv[1:]
    if argv[:1] == ['profile']:
        # Its arguments are plotProfile's own, passed on as they are.
        profile(argv[1:])
        return
    a = p.parse_args(argv)
    if a.what == 'fingerprint':
        found = fingerprint(a.npz, a.out, a.metrics, a.title)
    elif a.what == 'genecorr':
        found = genecorr(a.genes, a.expression, a.methylation, a.out, a.table, a.title, a.where)
    else:
        found = methexpr(a.genes, a.expression, a.methylation, a.out, a.table, max(2, min(10, a.groups)),
                         a.group_beds, a.title, a.where, a.ranked)
    if a.summary:
        with open(a.summary, 'w', encoding='utf-8') as fh:
            json.dump(found, fh)


if __name__ == '__main__':
    main(sys.argv[1:])
