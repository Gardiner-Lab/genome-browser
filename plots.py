"""
Figures the Genome Browser draws itself, where deepTools has no tool for a coverage (BigWig) track.

This file is also a script the recipes carry (it needs numpy, scipy and matplotlib, which come with
deepTools):

    python3 plots.py fingerprint summary.npz -o fingerprint.svg --metrics fingerprint-metrics.tab
    python3 plots.py methexpr --genes regions.bed --expression expression.tab --methylation methylation.tab \\
        -o methexpr.svg --table methylation-vs-expression.csv --groups 5 --group-beds group

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
each group (computeMatrix).
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


def methexpr(genes: str, expression: str, methylation: str, out: str, table: str, groups: int = 5,
             group_beds: Optional[str] = None, title: str = '', where: str = 'over each gene') -> dict:
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
            'tracks': tracks, 'beds': beds}


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
    a = p.parse_args(argv)
    if a.what == 'fingerprint':
        found = fingerprint(a.npz, a.out, a.metrics, a.title)
    else:
        found = methexpr(a.genes, a.expression, a.methylation, a.out, a.table, max(2, min(10, a.groups)),
                         a.group_beds, a.title, a.where)
    if a.summary:
        with open(a.summary, 'w', encoding='utf-8') as fh:
            json.dump(found, fh)


if __name__ == '__main__':
    main(sys.argv[1:])
