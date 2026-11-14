"""
Circos plots of a genome and the tracks shown on it (https://circos.ca).

The page sends the genome on screen, the tracks it shows (each with its colour,
and how to draw it) and the plot's options; this turns them into Circos's own
files and runs Circos (Martin Krzywinski's; GPL; the Debian package in this
image) on them:

    karyotype.txt        the genome's chromosomes, the ideogram round the outside
    data/track-N.txt     each coverage or feature track, binned along the genome:
                           coverage (BigWig)  the mean signal in each step
                           features (GFF)     genes starting in each step (all features if it has no genes)
                           features (BED)     features starting in each step
                         drawn as a ring: histogram, line or heatmap
    data/links-N.txt     each links track (BEDPE): arcs across the circle
    data/labels.txt      gene names asked for, round the inside of the ideogram
    data/highlight.txt   the selection, shaded through the rings
    circos.conf          all of it, outermost first

What is drawn: the whole genome (small chromosomes left out or not), chosen
chromosomes, or one region. The picture comes back as SVG with a legend added,
and the files come back zipped, to change and run again with Circos itself.

Only tracks the member may view are drawn. Circos runs with a time limit and,
in each web worker, one at a time.
"""

import base64
import io
import math
import re
import shutil
import sqlite3
import subprocess
import tempfile
import threading
import zipfile
from html import escape
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import bigwig_means
import chrom_names
import feature_index
import regions

MAX_TRACKS = 12
MAX_LINKS = 25_000
MAX_LABELS = 100
RESOLUTIONS = {'coarse': 300, 'medium': 800, 'fine': 2000}
MAX_STEPS = 5000
STYLES = ('histogram', 'line', 'heatmap')
TIMEOUT_SECONDS = 150
# A chromosome smaller than this share of the genome (organelles, scaffolds) is left
# out when asked: it would be a sliver on the circle.
SMALL_SHARE = 0.01
DEFAULT_COLOURS = {'bigwig': (100, 100, 100), 'bed': (0, 0, 150), 'gff': (0, 0, 150), 'bedpe': (150, 60, 160)}

_one_at_a_time = threading.Lock()

Span = Tuple[str, int, int]          # a chromosome, and the stretch of it drawn (0-based, half-open)


class CircosError(Exception):
    """Something the member can act on, as a sentence."""


def available() -> bool:
    return shutil.which('circos') is not None


def rgb(colour: Optional[str], kind: str) -> Tuple[int, int, int]:
    """A track's colour, as the page has it (#rgb, #rrggbb, rgb(...)), as r, g, b."""
    text = (colour or '').strip()
    m = re.fullmatch(r'#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})(?:[0-9a-fA-F]{2})?', text)
    if m:
        h = m.group(1)
        if len(h) == 3:
            h = ''.join(c * 2 for c in h)
        return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    m = re.fullmatch(r'rgba?\(\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*(?:,[^)]*)?\)', text)
    if m:
        return tuple(min(255, int(v)) for v in m.groups())
    return DEFAULT_COLOURS.get(kind, (0, 0, 150))


def nice(value: float) -> int:
    """value rounded to 1, 2 or 5 times a power of ten (a step size or a tick spacing)."""
    if value <= 1:
        return 1
    power = 10 ** int(math.floor(math.log10(value)))
    for step in (1, 2, 5, 10):
        if value <= step * power:
            return step * power
    return 10 * power


def _ident(chrom: str) -> str:
    """A chromosome's name as Circos takes it: no spaces or punctuation."""
    return re.sub(r'[^A-Za-z0-9_.]', '_', chrom)


# --- what is drawn ------------------------------------------------------------------

_REGION = re.compile(r'^\s*([^:\s]+)\s*:\s*([\d,._]+)\s*[-–]\s*([\d,._]+)\s*$')


def parse_region(text: str, sizes: Dict[str, int]) -> Span:
    """'Chr1:1,000,001-5,000,000' (1-based, inclusive, as the browser writes it) as a span."""
    m = _REGION.match(text or '')
    if not m:
        raise CircosError('Write the region as Chr1:1,000,000-5,000,000.')
    chrom = m.group(1)
    if chrom not in sizes:
        match = chrom_names.match_names([chrom], list(sizes))
        chrom = match.renames.get(chrom, chrom)
    if chrom not in sizes:
        raise CircosError(f'The genome has no chromosome {m.group(1)}.')
    start, end = (int(re.sub(r'[,._]', '', v)) for v in (m.group(2), m.group(3)))
    start, end = max(0, min(start, end) - 1), min(sizes[chrom], max(start, end))
    if end - start < 1000:
        raise CircosError('Choose a region of at least 1 kb.')
    return chrom, start, end


def choose_spans(sizes: Dict[str, int], chroms: Optional[List[str]] = None, region: Optional[str] = None,
                 skip_small: bool = True) -> List[Span]:
    if region:
        return [parse_region(region, sizes)]
    if chroms:
        unknown = [c for c in chroms if c not in sizes]
        if unknown:
            raise CircosError(f'The genome has no chromosome {", ".join(unknown)}.')
        return [(c, 0, sizes[c]) for c in sizes if c in set(chroms)]
    total = sum(sizes.values()) or 1
    chosen = [(c, 0, n) for c, n in sizes.items() if not skip_small or n >= SMALL_SHARE * total]
    return chosen or [(c, 0, n) for c, n in sizes.items()]


# Kept for callers of the first version: one chromosome, or the whole genome.
def choose_chromosomes(sizes: Dict[str, int], chrom: Optional[str], skip_small: bool) -> List[Tuple[str, int]]:
    if chrom and chrom not in sizes:
        raise CircosError(f'The genome has no chromosome {chrom}.')
    return [(c, e) for c, _, e in choose_spans(sizes, [chrom] if chrom else None, None, skip_small)]


def step_size(spans: List[Span], resolution: str, step: Optional[int]) -> int:
    total = sum(e - s for _, s, e in spans)
    floor = max(100, math.ceil(total / MAX_STEPS))
    if step:
        return max(int(step), floor)
    return max(nice(total / RESOLUTIONS.get(resolution, RESOLUTIONS['medium'])), 1)


def _bins(span: Span, size: int):
    _, s, e = span
    start = s
    while start < e:
        yield start, min(start + size, e)
        start += size


def where_text(spans: List[Span], sizes: Dict[str, int], everything: bool) -> str:
    if everything:
        return 'whole genome'
    if len(spans) == 1 and (spans[0][1] > 0 or spans[0][2] < sizes.get(spans[0][0], 0)):
        c, s, e = spans[0]
        return f'{c}:{s + 1:,}-{e:,}'
    return ', '.join(c for c, _, _ in spans)


# --- each track, binned -----------------------------------------------------------

def bigwig_bins(path: str, spans: List[Span], size: int) -> List[Tuple[str, int, int, float]]:
    import pyBigWig
    bw = pyBigWig.open(path)
    try:
        theirs = bw.chroms() or {}
        match = chrom_names.match_names(list(theirs), [c for c, _, _ in spans])
        to_theirs = {v: k for k, v in match.renames.items()}
        out = []
        for chrom, s, e in spans:
            name = to_theirs.get(chrom, chrom)
            if name not in theirs:
                continue
            e = min(e, theirs[name])
            # Bases without reads count as zero (bigwig_means.py): pyBigWig's own mean skips them,
            # which drew a bin with reads over a tenth of it as high as one covered throughout.
            values = bigwig_means.bin_means(bw, name, s, e, size)
            for (bs, be), value in zip(_bins((chrom, s, e), size), values):
                if not math.isnan(value):
                    out.append((chrom, bs, be, float(value)))
        return out
    finally:
        bw.close()


def _index(path: str, kind: str):
    idx = feature_index.ensure(path, regions.feature_parser(kind))
    return sqlite3.connect(f'file:{idx}?mode=ro', uri=True)


def feature_bins(path: str, kind: str, spans: List[Span], size: int) -> Tuple[List, str]:
    """Features starting in each step, from the track's index; for GFF its genes, if it has any."""
    conn = _index(path, kind)
    try:
        theirs = [r[0] for r in conn.execute('SELECT DISTINCT chrom FROM features')]
        match = chrom_names.match_names(theirs, [c for c, _, _ in spans])
        to_theirs = {v: k for k, v in match.renames.items()}
        where, what = '', 'features'
        if kind == 'gff' and conn.execute("SELECT 1 FROM features WHERE type = 'gene' LIMIT 1").fetchone():
            where, what = " AND type = 'gene'", 'genes'
        out = []
        for chrom, s, e in spans:
            name = to_theirs.get(chrom, chrom)
            counts = dict(conn.execute(
                f'SELECT (start - ?) / ?, count(*) FROM features WHERE chrom = ? AND start >= ? AND start < ?{where} '
                'GROUP BY 1', (s, size, name, s, e)).fetchall())
            for i, (bs, be) in enumerate(_bins((chrom, s, e), size)):
                out.append((chrom, bs, be, float(counts.get(i, 0))))
        return out, what
    finally:
        conn.close()


def gene_labels(paths: List[Tuple[str, str]], names: List[str], spans: List[Span]) -> Tuple[List, List[str]]:
    """The genes named (by their name in the track, any case), found in the GFF tracks given, within the spans."""
    wanted = {n.strip().casefold(): n.strip() for n in names if n.strip()}
    found: Dict[str, Tuple[str, int, int, str]] = {}
    for path, kind in paths:
        if not wanted.keys() - found.keys():
            break
        conn = _index(path, kind)
        try:
            chroms = [r[0] for r in conn.execute('SELECT DISTINCT chrom FROM features')]
            renames = chrom_names.match_names(chroms, [c for c, _, _ in spans]).renames
            todo = sorted(wanted.keys() - found.keys())
            for i in range(0, len(todo), 400):
                chunk = todo[i:i + 400]
                for chrom, start, end, name, ftype in conn.execute(
                        f'SELECT chrom, start, end, name, type FROM features WHERE lower(name) IN '
                        f'({", ".join("?" * len(chunk))}) ORDER BY type = \'gene\' DESC', chunk):
                    key = (name or '').casefold()
                    if key in found:
                        continue
                    chrom = renames.get(chrom, chrom)
                    if any(c == chrom and start < e and end > s for c, s, e in spans):
                        found[key] = (chrom, start, end, name)
        finally:
            conn.close()
    missing = [wanted[k] for k in wanted if k not in found]
    return list(found.values()), missing


def links(path: str, spans: List[Span], limit: int = MAX_LINKS) -> Tuple[List, int]:
    """A BEDPE track's links with both ends drawn, strongest first (by its score column), at most limit."""
    chroms = {c for c, _, _ in spans}
    out = []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for n, line in enumerate(fh):
            if line.startswith(('#', 'track ', 'browser ')) or line.lower().startswith('chrom1'):
                continue
            f = line.rstrip('\r\n').split('\t') if '\t' in line else line.split()
            if len(f) < 6 or f[0] not in chroms or f[3] not in chroms:
                continue
            try:
                a, b = (f[0], int(f[1]), int(f[2])), (f[3], int(f[4]), int(f[5]))
            except ValueError:
                continue
            if not all(any(c == x[0] and x[1] < e and x[2] > s for c, s, e in spans) for x in (a, b)):
                continue
            try:
                score = float(f[7]) if len(f) > 7 else 0.0
            except ValueError:
                score = 0.0
            out.append((score, -n, a, b))
    total = len(out)
    out.sort(reverse=True)
    return [(a, b) for _, _, a, b in out[:limit]], total


# --- Circos's files -----------------------------------------------------------------

def _units(total: int) -> Tuple[int, str]:
    if total >= 10_000_000:
        return 1_000_000, 'Mb'
    if total >= 10_000:
        return 1_000, 'kb'
    return 1, 'bp'


def karyotype(chroms) -> str:
    """chroms: (name, length) pairs, or spans (whose chromosomes' lengths are given by their ends)."""
    lines = []
    for i, item in enumerate(chroms):
        chrom, length = item[0], item[-1]
        ident = _ident(chrom)
        lines.append(f'chr - {ident} {ident} 0 {length} {"vlgrey" if i % 2 else "lgrey"}')
    return '\n'.join(lines) + '\n'


def data_file(rows) -> str:
    return ''.join(f'{_ident(c)} {s} {max(s, e - 1)} {v:.6g}\n' for c, s, e, v in rows)


def links_file(pairs) -> str:
    return ''.join(f'{_ident(a[0])} {a[1]} {max(a[1], a[2] - 1)} {_ident(b[0])} {b[1]} {max(b[1], b[2] - 1)}\n'
                   for a, b in pairs)


def labels_file(found) -> str:
    return ''.join(f'{_ident(c)} {s} {max(s, e - 1)} {re.sub(r"[^A-Za-z0-9_.:-]", "_", name)}\n'
                   for c, s, e, name in found)


def _gradient(colour, steps: int = 9) -> List[Tuple[int, int, int]]:
    """From near white to the track's colour, for a heatmap."""
    return [tuple(round(255 + (c - 255) * (i + 1) / steps) for c in colour) for i in range(steps)]


def conf(spans, rings: List[dict], links_layers: Optional[List[dict]] = None, labels: Optional[str] = None,
         highlight: Optional[str] = None, sizes: Optional[Dict[str, int]] = None) -> str:
    """spans: (chrom, start, end), or (chrom, length) pairs for the whole of each chromosome."""
    spans = [s if len(s) == 3 else (s[0], 0, s[1]) for s in spans]
    sizes = sizes or {c: e for c, _, e in spans}
    links_layers = links_layers or []
    total = sum(e - s for _, s, e in spans)
    unit, _ = _units(total)
    spacing = nice(total / unit / 60)
    label_every = spacing * (5 if len(spans) > 1 else 2)
    whole = all(s == 0 and e >= sizes.get(c, e) for c, s, e in spans)
    shown = ';'.join(_ident(c) if (s == 0 and e >= sizes.get(c, e)) else f'{_ident(c)}:{s / unit:.6f}-{e / unit:.6f}'
                     for c, s, e in spans)

    colours = []
    blocks = []
    top = 0.97
    if labels:
        blocks.append(f"""<plot>
type          = text
file          = {labels}
r0            = 0.84r
r1            = 0.99r
label_size    = 24p
label_font    = default
color         = black
show_links    = yes
link_dims     = 2p,4p,6p,4p,2p
link_thickness = 1p
link_color    = grey
label_snuggle = yes
max_snuggle_distance = 2r
</plot>""")
        top = 0.82
    inner = 0.42 if links_layers else 0.25
    gap = 0.015
    weights = [max(1, min(3, int(r.get('height') or 2))) for r in rings]
    unit_h = min(0.06, (top - inner - gap * max(0, len(rings) - 1)) / max(1, sum(weights)))
    r1 = top
    for i, (ring, weight) in enumerate(zip(rings, weights)):
        r0 = r1 - unit_h * weight
        colour = ','.join(str(v) for v in ring['rgb'])
        style = ring.get('style') if ring.get('style') in STYLES else 'histogram'
        scale = ''.join(f'\n{k:<11} = {v:.6g}' for k, v in (('min', ring.get('min', 0)), ('max', ring.get('max')))
                        if v is not None)
        if style == 'heatmap':
            names = []
            for j, c in enumerate(_gradient(ring['rgb'])):
                names.append(f'ring{i + 1}_{j}')
                colours.append(f'{names[-1]} = {c[0]},{c[1]},{c[2]}')
            body = f'color       = {",".join(names)}\nstroke_thickness = 0'
        elif style == 'line':
            body = f'color       = {colour}\nthickness   = 3p'
        else:
            body = f'fill_color  = {colour}\ncolor       = {colour}\nthickness   = 1p\nextend_bin  = no'
        extras = '' if style == 'heatmap' else """
<backgrounds>
<background>
color = vvlgrey
</background>
</backgrounds>
<axes>
<axis>
color     = lgrey
thickness = 1
spacing   = 0.5r
</axis>
</axes>"""
        blocks.append(f"""<plot>
type        = {style}
file        = {ring['file']}
r0          = {r0:.3f}r
r1          = {r1:.3f}r
{body}{scale}{extras}
</plot>""")
        r1 = r0 - gap
    link_r = max(0.15, r1 - 0.01)
    link_blocks = [f"""<link>
file          = {layer['file']}
radius        = {link_r:.3f}r
bezier_radius = 0.1r
color         = {','.join(str(v) for v in layer['rgb'])},0.5
thickness     = 2p
</link>""" for layer in links_layers]
    # Over the rings (not the gene labels), see-through, so the bars still show beneath: a plot of
    # its own drawn last, since Circos draws its highlights layer under everything else.
    if highlight:
        blocks.append(f"""<plot>
type       = highlight
file       = {highlight}
r0         = {link_r:.3f}r
r1         = {top:.3f}r
fill_color = 255,190,0,0.3
stroke_thickness = 0
z          = 50
</plot>""")
    highlight_block = ''
    colour_block = ('<colors>\n' + '\n'.join(colours) + '\n</colors>\n') if colours else ''
    return f"""# Made by the Trellis Tools Genome Browser. Run: circos -conf circos.conf
karyotype = karyotype.txt
chromosomes_units = {unit}
chromosomes_display_default = {'yes' if whole else 'no'}
{'' if whole else f'chromosomes = {shown}'}
{colour_block}
<ideogram>
<spacing>
default = {0.006 if len(spans) > 1 else 0.0}r
</spacing>
radius           = 0.82r
thickness        = 40p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = default
label_radius     = 1r + 90p
label_size       = 48
label_parallel   = yes
</ideogram>

show_ticks       = yes
show_tick_labels = yes
<ticks>
radius     = 1r
color      = black
thickness  = 2p
multiplier = {1 / unit if unit > 1 else 1}
format     = {'%.1f' if spacing < 1 else '%d'}
<tick>
spacing    = {spacing}u
size       = 12p
</tick>
<tick>
spacing        = {label_every}u
size           = 22p
show_label     = yes
label_size     = 26p
label_offset   = 6p
</tick>
</ticks>
{highlight_block}
<plots>
{chr(10).join(blocks)}
</plots>
{'<links>' + chr(10) + chr(10).join(link_blocks) + chr(10) + '</links>' if link_blocks else ''}

<image>
<<include etc/image.conf>>
</image>
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
"""


def legend_svg(svg: str, title: str, rings: List[dict]) -> str:
    """The rings named, outermost first, in the top-left corner of Circos's picture."""
    rows = [f'<text x="40" y="70" font-family="Helvetica, Arial, sans-serif" font-size="44" '
            f'font-weight="bold" fill="#222">{escape(title)}</text>']
    for i, ring in enumerate(rings):
        y = 130 + i * 52
        colour = 'rgb({},{},{})'.format(*ring['rgb'])
        mark = (f'<line x1="40" y1="{y - 13}" x2="74" y2="{y - 13}" stroke="{colour}" stroke-width="6"/>'
                if ring.get('kind') == 'links' else f'<rect x="40" y="{y - 30}" width="34" height="34" fill="{colour}"/>')
        number = '' if ring.get('kind') == 'links' else f'{i + 1}. '
        rows.append(f'{mark}<text x="90" y="{y}" font-family="Helvetica, Arial, sans-serif" font-size="32" fill="#222">'
                    f'{number}{escape(ring["name"])} ({escape(ring["what"])})</text>')
    block = '<g id="trellis-legend">' + ''.join(rows) + '</g>'
    at = svg.rfind('</svg>')
    return svg[:at] + block + svg[at:] if at >= 0 else svg


def make(genome: dict, sizes: Dict[str, int], tracks: List[dict], chrom: Optional[str] = None,
         resolution: str = 'medium', skip_small: bool = True, chroms: Optional[List[str]] = None,
         region: Optional[str] = None, step: Optional[int] = None, shared_scale: bool = False,
         log_scale: bool = False, highlight: Optional[dict] = None, labels: Optional[List[str]] = None,
         title: Optional[str] = None) -> dict:
    """
    tracks: [{'id', 'name', 'type', 'file_path', 'color', 'style', 'height'}], viewable ones, in the order shown.
    Returns {'svg', 'bundle' (base64 zip), 'rings', 'bin_size', 'where', 'notes'}.
    """
    if not available():
        raise CircosError('Circos is not installed in this image.')
    if not tracks:
        raise CircosError('Show at least one track to draw.')
    if len(tracks) > MAX_TRACKS:
        raise CircosError(f'Circos draws at most {MAX_TRACKS} tracks at a time: hide some first.')
    if chrom and not chroms and not region:
        chroms = [chrom]
    spans = choose_spans(sizes, chroms, region, skip_small)
    everything = not chroms and not region
    size = step_size(spans, resolution, step)
    total = sum(e - s for _, s, e in spans)
    notes: List[str] = []

    work = Path(tempfile.mkdtemp(prefix='circos-'))
    try:
        (work / 'data').mkdir()
        shown_chroms = []
        for c, _, _ in spans:
            if c not in shown_chroms:
                shown_chroms.append(c)
        (work / 'karyotype.txt').write_text(karyotype([(c, sizes[c]) for c in shown_chroms]))
        rings, layers, legend = [], [], []
        for i, t in enumerate(tracks, 1):
            colour = rgb(t.get('color'), t['type'])
            if t['type'] == 'bedpe':
                pairs, found = links(t['file_path'], spans)
                name = f'data/links-{i}.txt'
                (work / name).write_text(links_file(pairs))
                if found > len(pairs):
                    notes.append(f"{t['name']}: the {len(pairs):,} strongest of {found:,} links are drawn.")
                layers.append({'file': name, 'rgb': colour})
                legend.append({'name': t['name'], 'what': f'{len(pairs):,} links', 'rgb': colour, 'kind': 'links'})
                continue
            if t['type'] == 'bigwig':
                rows, what = bigwig_bins(t['file_path'], spans, size), 'mean coverage'
                if log_scale:
                    rows = [(c, s, e, math.log10(1 + max(0.0, v))) for c, s, e, v in rows]
                    what = 'mean coverage, log10(1 + x)'
            elif t['type'] in ('bed', 'gff'):
                rows, what = feature_bins(t['file_path'], t['type'], spans, size)
                what += ' per step'
            else:
                continue
            name = f'data/track-{i}.txt'
            (work / name).write_text(data_file(rows))
            ring = {'name': t['name'], 'what': what, 'file': name, 'rgb': colour, 'kind': t['type'],
                    'style': t.get('style'), 'height': t.get('height'),
                    'top': max((r[3] for r in rows), default=0.0)}
            rings.append(ring)
            legend.append(ring)
        if not rings and not layers:
            raise CircosError('None of the tracks shown can be drawn by Circos.')
        if shared_scale:
            coverage = [r for r in rings if r['kind'] == 'bigwig']
            top = max((r['top'] for r in coverage), default=0)
            for r in coverage:
                r['max'] = top or 1
            if len(coverage) > 1:
                notes.append('Coverage rings share one scale.')

        label_file = None
        if labels:
            gffs = [(t['file_path'], 'gff') for t in tracks if t['type'] == 'gff']
            if not gffs:
                notes.append('Gene labels need a gene annotation (GFF) track shown.')
            else:
                found, missing = gene_labels(gffs, labels[:MAX_LABELS], spans)
                if found:
                    label_file = 'data/labels.txt'
                    (work / label_file).write_text(labels_file(found))
                if missing:
                    notes.append('Not found in the region drawn: ' + ', '.join(missing[:20])
                                 + (' and more' if len(missing) > 20 else '') + '.')
        highlight_file = None
        if highlight:
            try:
                hc, hs, he = str(highlight['chrom']), int(highlight['start']), int(highlight['end'])
            except (KeyError, TypeError, ValueError):
                hc = None
            if hc and any(c == hc and hs < e and he > s for c, s, e in spans):
                highlight_file = 'data/highlight.txt'
                (work / highlight_file).write_text(f'{_ident(hc)} {hs} {max(hs, he - 1)}\n')

        (work / 'circos.conf').write_text(conf(spans, rings, layers, label_file, highlight_file, sizes))
        # The configuration includes Circos's standard files as etc/..., as an ordinary Circos install
        # has them (so the files downloaded run anywhere); Debian keeps them in /etc/circos.
        if Path('/etc/circos/image.conf').exists():
            (work / 'etc').symlink_to('/etc/circos', target_is_directory=True)
        _, unit_name = _units(total)
        where = where_text(spans, sizes, everything)
        heading = (title or '').strip()[:120] or f"{genome['name']} ({genome['assembly_id']}), {where}"
        (work / 'README.txt').write_text(
            f"A Circos plot of {genome['name']} ({genome['assembly_id']}), {where}, from the Trellis Tools Genome Browser.\n"
            f'Steps of {size:,} bp. Ticks in {unit_name}.\n\nRings, outermost first:\n'
            + ''.join(f"  {i}. {r['name']}: {r['what']} ({r['file']}, drawn as {r.get('style') or 'histogram'})\n"
                      for i, r in enumerate(rings, 1))
            + ''.join(f"Links: {layer['file']}\n" for layer in layers)
            + ('Gene labels: data/labels.txt\n' if label_file else '')
            + ('Highlighted: data/highlight.txt\n' if highlight_file else '')
            + ''.join(f'Note: {n}\n' for n in notes)
            + '\nTo change it, edit circos.conf (and the data files) and run Circos (https://circos.ca):\n'
              '    circos -conf circos.conf\n')
        with _one_at_a_time:
            try:
                done = subprocess.run(['circos', '-conf', 'circos.conf', '-outputdir', str(work),
                                       '-outputfile', 'circos', '-silent'],
                                      cwd=work, capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
            except subprocess.TimeoutExpired:
                raise CircosError('Circos took too long: try less detail or fewer tracks.')
        svg_path = work / 'circos.svg'
        if done.returncode != 0 or not svg_path.exists():
            detail = (done.stderr or done.stdout or '').strip().splitlines()[-3:]
            raise CircosError('Circos could not draw it' + (': ' + ' '.join(detail) if detail else '.'))
        svg = legend_svg(svg_path.read_text(encoding='utf-8', errors='replace'), heading, legend)

        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as z:
            names = ['circos.conf', 'karyotype.txt', 'README.txt'] + [r['file'] for r in rings] + \
                    [layer['file'] for layer in layers] + [f for f in (label_file, highlight_file) if f]
            for name in names:
                z.write(work / name, f'circos/{name}')
            z.writestr('circos/circos.svg', svg)
            if (work / 'circos.png').exists():
                z.write(work / 'circos.png', 'circos/circos.png')
        return {'svg': svg, 'bundle': base64.b64encode(buf.getvalue()).decode('ascii'),
                'rings': [{'name': r['name'], 'what': r['what'], 'color': 'rgb({},{},{})'.format(*r['rgb'])}
                          for r in legend],
                'bin_size': size, 'where': where, 'notes': notes, 'title': heading}
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --- the endpoint ------------------------------------------------------------------

def register_circos(app) -> None:
    from flask import jsonify, request

    import access
    import db
    import identity

    @app.route('/api/genomes/<genome_id>/circos', methods=['POST'])
    def circos_plot(genome_id):
        """
        A Circos plot of this genome and the tracks named, in the order and colours given:
        {tracks: [{id, color, style, height}], chroms: [...] or region: 'Chr1:1-5,000,000' (else the whole
         genome), resolution: coarse|medium|fine, step: bp, skip_small, shared_scale, log_scale,
         highlight: {chrom, start, end}, labels: [gene names], title}.
        Answers {svg, bundle, rings, bin_size, where, notes, title}.
        """
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401
        body = request.get_json(silent=True) or {}
        asked = body.get('tracks') if isinstance(body.get('tracks'), list) else []
        conn = db.connect(app.config['DATA_ROOT'] / 'genome.db')
        try:
            genome = conn.execute('SELECT id, name, assembly_id, owner_lab, genome_visibility FROM reference_genomes '
                                  'WHERE id = ?', (genome_id,)).fetchone()
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            sizes = dict(conn.execute('SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ? ORDER BY rowid',
                                      (genome_id,)).fetchall())
            tracks = []
            for item in asked[:MAX_TRACKS + 1]:
                if not isinstance(item, dict):
                    continue
                row = conn.execute('SELECT id, genome_id, name, type, owner_lab, visibility, track_status, file_path '
                                   'FROM tracks WHERE id = ?', (str(item.get('id') or ''),)).fetchone()
                if (row is None or row['genome_id'] != genome_id or not row['file_path']
                        or not access.can_view_track(member, row, conn)):
                    continue
                tracks.append({'id': row['id'], 'name': row['name'], 'type': row['type'],
                               'file_path': row['file_path'], 'color': str(item.get('color') or '')[:40],
                               'style': str(item.get('style') or 'histogram'), 'height': _int(item.get('height'), 2)})
        finally:
            conn.close()
        chroms = body.get('chroms') if isinstance(body.get('chroms'), list) else None
        highlight = body.get('highlight') if isinstance(body.get('highlight'), dict) else None
        labels = body.get('labels')
        if isinstance(labels, str):
            labels = re.split(r'[\s,;]+', labels)
        labels = [str(x)[:60] for x in labels if str(x).strip()][:MAX_LABELS] if isinstance(labels, list) else None
        try:
            out = make({'name': genome['name'], 'assembly_id': genome['assembly_id']}, sizes, tracks,
                       chrom=(str(body['chrom']) if body.get('chrom') else None),
                       resolution=str(body.get('resolution') or 'medium'),
                       skip_small=body.get('skip_small', True) is not False,
                       chroms=[str(c) for c in chroms][:200] if chroms else None,
                       region=str(body['region'])[:100] if body.get('region') else None,
                       step=_int(body.get('step'), None),
                       shared_scale=bool(body.get('shared_scale')), log_scale=bool(body.get('log_scale')),
                       highlight=highlight, labels=labels, title=str(body.get('title') or '')[:120])
        except CircosError as exc:
            return jsonify({'error': str(exc)}), 400
        return jsonify(out)


def _int(value, default):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default
