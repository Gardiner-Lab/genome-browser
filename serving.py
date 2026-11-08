"""
The files the in-browser viewer reads: a genome's sequence and each track's
data, indexed so igv.js fetches only the window on screen.

igv.js reads everything by HTTP range request. A reference is a FASTA plus its
``.fai`` index; a BigWig carries its own index; a BED or GFF is read through a
coordinate-sorted, bgzip-compressed copy and its tabix ``.tbi`` index, which is
also what lets igv.js put a GFF's exons, UTRs and CDS back together into gene
models. None of these are what was uploaded, so this module makes them:

* ``reference.fasta.fai`` with ``samtools faidx``. A FASTA whose lines are not
  of one length within a sequence cannot be indexed as it is, so it is rewritten
  once with 60-base lines (the sequence itself unchanged) and indexed then.
* ``<track file>.bgz`` and ``.bgz.tbi`` with ``sort``, ``bgzip`` and ``tabix``:
  comment lines first, then every feature sorted by chromosome and start.

Each is made on upload when the tools are there, and otherwise on the first
request for it; it is written under a temporary name and moved into place, and
remade when its source is newer. Every route checks the same visibility as the
rest of the module (``access.can_view_genome`` / ``access.can_view_track``) and
answers 404 for what the member may not see.

Routes (relative to the stripped ``/genome`` mount):

- ``GET /api/genomes/<id>/reference.fa`` and ``reference.fa.fai``
- ``GET /api/genomes/<id>/sequence.fasta?chrom=&start=&end=&strand=&name=`` — the
  same as a FASTA download, any length
- ``GET /api/genomes/<id>/sequence?chrom=&start=&end=&strand=&tracks=`` — a
  stretch of sequence and the annotations in it, for the Sequence Editor and BLAST
- ``GET /api/tracks/<id>/coverage.bw`` — a BigWig track's file
- ``GET /api/tracks/<id>/features.<bed|gff3>.gz`` and ``….gz.tbi`` — a feature
  track's sorted, compressed copy and its index
"""

import os
import re
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Optional

from flask import Flask, Response, abort, current_app, jsonify, request, send_file

import access
import db
import feature_index
import identity
import regions


class ToolMissing(RuntimeError):
    """A tool this module needs (samtools, bgzip, tabix, sort) is not installed."""


def register_serving(app: Flask) -> None:

    @app.route('/api/genomes/<genome_id>/reference.fa')
    def reference_fasta(genome_id):
        path = _viewable_reference(genome_id)
        # Make sure the index exists before the sequence is read: igv.js asks
        # for the index first, but a FASTA rewrite must never happen under a
        # reader that already has the old line layout.
        _ensure_or_503(lambda: ensure_fasta_index(path))
        return _send(path, 'text/plain')

    @app.route('/api/genomes/<genome_id>/reference.fa.fai')
    def reference_index(genome_id):
        path = _viewable_reference(genome_id)
        index = _ensure_or_503(lambda: ensure_fasta_index(path))
        return _send(index, 'text/plain')

    @app.route('/api/genomes/<genome_id>/sequence')
    def sequence_of_region(genome_id):
        """
        A stretch of a genome's sequence, with the annotations inside it.

        Query: ``chrom``, ``start`` (0-based) and ``end`` (exclusive), up to
        :data:`MAX_REGION` bases; ``strand=-1`` returns the reverse complement,
        as a feature on the minus strand reads; ``tracks`` is a comma-separated
        list of feature tracks (of this genome, visible to the member) whose
        features inside the region come back too, on the region's own
        coordinates and strand, in the Sequence Editor's shape. This is what the
        browser sends to the Sequence Editor and to BLAST.
        """
        member = identity.require_member()
        path = _viewable_reference(genome_id)
        chrom = (request.args.get('chrom') or '').strip()
        try:
            start = int(request.args.get('start', ''))
            end = int(request.args.get('end', ''))
        except ValueError:
            return jsonify({'error': 'start and end must be whole numbers.'}), 400
        strand = -1 if request.args.get('strand') in ('-1', '-') else 1

        conn = db.connect(_data_root() / 'genome.db')
        try:
            row = conn.execute(
                'SELECT length_bp FROM chrom_sizes WHERE genome_id = ? AND chrom = ?',
                (genome_id, chrom)).fetchone()
            if row is None:
                return jsonify({'error': f'This genome has no sequence called "{chrom}".'}), 400
            length = int(row['length_bp'])
            start, end = max(0, start), min(end, length)
            if end <= start:
                return jsonify({'error': 'The region is empty.'}), 400
            if end - start > MAX_REGION:
                return jsonify({
                    'error': f'The region is {end - start:,} bp; at most '
                             f'{MAX_REGION:,} bp can be sent at once.'}), 400
            tracks = _region_tracks(conn, member, genome_id,
                                    request.args.get('tracks', ''))
        finally:
            conn.close()

        seq = _ensure_or_503(lambda: read_region(path, chrom, start, end))
        if strand == -1:
            seq = reverse_complement(seq)
        features = []
        for track in tracks:
            found = feature_index.details(
                track['file_path'], regions.feature_parser(track['type']),
                chrom, start, end) or []
            features.extend(_on_region(found, start, end, strand))
        return jsonify({'chrom': chrom, 'start': start, 'end': end,
                        'strand': strand, 'sequence': seq,
                        'features': features[:MAX_FEATURES]})

    @app.route('/api/genomes/<genome_id>/sequence.fasta')
    def sequence_as_fasta(genome_id):
        """
        A stretch of sequence as a FASTA file to download, of any length.

        The same query as ``/sequence`` (``chrom``, ``start`` 0-based, ``end``
        exclusive, ``strand``, plus an optional ``name`` for the header line),
        but streamed from ``samtools faidx`` rather than built in memory, so a
        region too long for the Sequence Editor, or a whole chromosome, can
        still be taken away.
        """
        path = _viewable_reference(genome_id)
        chrom = (request.args.get('chrom') or '').strip()
        try:
            start = int(request.args.get('start', ''))
            end = int(request.args.get('end', ''))
        except ValueError:
            return jsonify({'error': 'start and end must be whole numbers.'}), 400
        strand = -1 if request.args.get('strand') in ('-1', '-') else 1
        conn = db.connect(_data_root() / 'genome.db')
        try:
            genome = conn.execute('SELECT assembly_id FROM reference_genomes WHERE id = ?',
                                  (genome_id,)).fetchone()
            row = conn.execute(
                'SELECT length_bp FROM chrom_sizes WHERE genome_id = ? AND chrom = ?',
                (genome_id, chrom)).fetchone()
        finally:
            conn.close()
        if row is None:
            return jsonify({'error': f'This genome has no sequence called "{chrom}".'}), 400
        start, end = max(0, start), min(end, int(row['length_bp']))
        if end <= start:
            return jsonify({'error': 'The region is empty.'}), 400

        where = f'{chrom}:{start + 1}-{end}' + ('(-)' if strand == -1 else '')
        name = re.sub(r'[^\w.:()+-]+', '_', (request.args.get('name') or '').strip())[:80]
        title = ' '.join(p for p in (name, genome['assembly_id'], where) if p)
        _ensure_or_503(lambda: ensure_fasta_index(path))
        samtools = _ensure_or_503(lambda: _need('samtools'))
        args = [samtools, 'faidx', '-n', '60']
        if strand == -1:
            args.append('-i')
        args += [str(path), f'{chrom}:{start + 1}-{end}']

        def body():
            proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            try:
                first = True
                for line in proc.stdout:
                    if first:
                        # samtools names the record itself; this header says
                        # what it is and where it came from.
                        first = False
                        yield f'>{title}\n'.encode('utf-8')
                        continue
                    yield line
            finally:
                proc.stdout.close()
                proc.wait()

        filename = re.sub(r'[^\w.-]+', '_', title).strip('_') + '.fasta'
        return Response(body(), mimetype='text/plain', headers={
            'Content-Disposition': f'attachment; filename="{filename}"',
            'Cache-Control': 'private, no-store'})

    @app.route('/api/tracks/<track_id>/coverage.bw')
    def track_bigwig(track_id):
        track = _viewable_track(track_id)
        if track['type'] != 'bigwig':
            abort(404)
        return _send(Path(track['file_path']), 'application/octet-stream')

    @app.route('/api/tracks/<track_id>/links.bedpe')
    def track_links(track_id):
        """A links track (BEDPE), whole: the viewer draws it as arcs."""
        track = _viewable_track(track_id)
        if track['type'] != 'bedpe':
            abort(404)
        return _send(Path(track['file_path']), 'text/plain')

    @app.route('/api/tracks/<track_id>/features.<ext>.gz')
    def track_features(track_id, ext):
        track = _viewable_feature_track(track_id, ext)
        data, _ = _ensure_or_503(
            lambda: ensure_tabix(Path(track['file_path']), track['type']))
        return _send(data, 'application/octet-stream')

    @app.route('/api/tracks/<track_id>/features.<ext>.gz.tbi')
    def track_features_index(track_id, ext):
        track = _viewable_feature_track(track_id, ext)
        _, index = _ensure_or_503(
            lambda: ensure_tabix(Path(track['file_path']), track['type']))
        return _send(index, 'application/octet-stream')


# ---------------------------------------------------------------------------
# Access
# ---------------------------------------------------------------------------


def _data_root() -> Path:
    return Path(current_app.config['DATA_ROOT'])


def _viewable_reference(genome_id: str) -> Path:
    member = identity.require_member()
    conn = db.connect(_data_root() / 'genome.db')
    try:
        genome = conn.execute(
            'SELECT id, owner_lab, genome_visibility FROM reference_genomes '
            'WHERE id = ?', (genome_id,)).fetchone()
        if genome is None or not access.can_view_genome(member, genome, conn):
            abort(404)
    finally:
        conn.close()
    path = reference_path(_data_root(), genome_id)
    if not path.exists():
        abort(404)
    return path


def _viewable_track(track_id: str):
    member = identity.require_member()
    conn = db.connect(_data_root() / 'genome.db')
    try:
        track = conn.execute(
            'SELECT id, genome_id, type, owner_lab, visibility, track_status, '
            'file_path FROM tracks WHERE id = ?', (track_id,)).fetchone()
        if track is None or not access.can_view_track(member, track, conn):
            abort(404)
    finally:
        conn.close()
    if not track['file_path'] or not Path(track['file_path']).exists():
        abort(404)
    return track


def _viewable_feature_track(track_id: str, ext: str):
    track = _viewable_track(track_id)
    if track['type'] not in ('bed', 'gff') or ext != feature_extension(track['type']):
        abort(404)
    return track


def _send(path: Path, mimetype: str):
    # conditional=True answers Range requests with 206, which is how igv.js reads.
    response = send_file(str(path), mimetype=mimetype, conditional=True,
                         max_age=0)
    response.headers['Cache-Control'] = 'private, no-cache'
    return response


def _ensure_or_503(make):
    try:
        return make()
    except ToolMissing as exc:
        current_app.logger.error('Cannot index for the viewer: %s', exc)
        abort(503)
    except (OSError, subprocess.CalledProcessError) as exc:
        current_app.logger.error('Indexing for the viewer failed: %s', exc)
        abort(500)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def reference_path(data_root: Path, genome_id: str) -> Path:
    return Path(data_root) / 'genomes' / genome_id / 'reference.fasta'


def feature_extension(track_type: str) -> str:
    """The format a feature track is served as: ``bed`` or ``gff3``."""
    return 'bed' if track_type == 'bed' else 'gff3'


def _newer_or_same(made: Path, source: Path) -> bool:
    try:
        return made.stat().st_mtime >= source.stat().st_mtime
    except OSError:
        return False


def _need(tool: str) -> str:
    found = shutil.which(tool)
    if found is None:
        raise ToolMissing(f'{tool} is not installed')
    return found


# ---------------------------------------------------------------------------
# The reference: samtools faidx
# ---------------------------------------------------------------------------


def ensure_fasta_index(fasta: Path) -> Path:
    """The ``.fai`` beside a FASTA, made (and the FASTA re-wrapped) when needed."""
    fasta = Path(fasta)
    index = Path(str(fasta) + '.fai')
    if _newer_or_same(index, fasta):
        return index
    samtools = _need('samtools')
    tmp_dir = Path(tempfile.mkdtemp(prefix='faidx-', dir=str(fasta.parent)))
    try:
        # Index a hard link in a private directory, so the .fai is written
        # there and only moved beside the FASTA once it is complete.
        staged = tmp_dir / fasta.name
        os.link(fasta, staged)
        result = subprocess.run([samtools, 'faidx', str(staged)],
                                capture_output=True, text=True)
        if result.returncode != 0:
            # Uneven line lengths: re-wrap the sequence and index that.
            staged.unlink()
            _rewrap_fasta(fasta, staged)
            subprocess.run([samtools, 'faidx', str(staged)],
                           capture_output=True, text=True, check=True)
            os.replace(staged, fasta)
        os.replace(Path(str(staged) + '.fai'), index)
        # The index must not look older than a FASTA just re-wrapped.
        os.utime(index)
        return index
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _rewrap_fasta(source: Path, target: Path, width: int = 60) -> None:
    """Copy a FASTA with every sequence wrapped at ``width`` bases."""
    with open(source, 'r', encoding='ascii', errors='replace') as inp, \
            open(target, 'w', encoding='ascii', newline='\n') as out:
        pending = ''
        for raw in inp:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if pending:
                    out.write(pending + '\n')
                    pending = ''
                out.write(line + '\n')
                continue
            pending += line
            while len(pending) >= width:
                out.write(pending[:width] + '\n')
                pending = pending[width:]
        if pending:
            out.write(pending + '\n')


# ---------------------------------------------------------------------------
# Feature tracks: sort | bgzip, tabix
# ---------------------------------------------------------------------------


def ensure_tabix(source: Path, track_type: str):
    """``(data, index)``: the sorted bgzip copy of a BED/GFF and its .tbi."""
    source = Path(source)
    data = Path(str(source) + '.bgz')
    index = Path(str(data) + '.tbi')
    if _newer_or_same(data, source) and _newer_or_same(index, data):
        return data, index
    sort, bgzip, tabix = _need('sort'), _need('bgzip'), _need('tabix')
    tmp_dir = Path(tempfile.mkdtemp(prefix='tabix-', dir=str(source.parent)))
    try:
        header = tmp_dir / 'header'
        body = tmp_dir / 'body'
        _split_header_and_body(source, header, body, track_type)

        sorted_body = tmp_dir / 'sorted'
        start_key = '-k2,2n' if track_type == 'bed' else '-k4,4n'
        env = dict(os.environ, LC_ALL='C')
        subprocess.run([sort, '-t', '\t', '-k1,1', start_key, '-S', '256M',
                        '-T', str(tmp_dir), '-o', str(sorted_body), str(body)],
                       check=True, capture_output=True, env=env)

        staged = tmp_dir / 'features.bgz'
        with open(staged, 'wb') as out:
            proc = subprocess.Popen([bgzip, '-c'], stdin=subprocess.PIPE,
                                    stdout=out)
            for part in (header, sorted_body):
                with open(part, 'rb') as handle:
                    shutil.copyfileobj(handle, proc.stdin, 1024 * 1024)
            proc.stdin.close()
            if proc.wait() != 0:
                raise subprocess.CalledProcessError(proc.returncode, 'bgzip')
        preset = 'bed' if track_type == 'bed' else 'gff'
        indexed = subprocess.run([tabix, '-f', '-p', preset, str(staged)],
                                 capture_output=True, text=True)
        if indexed.returncode != 0:
            raise OSError(f'tabix could not index the track: '
                          f'{indexed.stderr.strip()}')
        os.replace(staged, data)
        os.replace(Path(str(staged) + '.tbi'), index)
        return data, index
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# Control characters other than tab, and the Unicode line separators, become
# spaces in the viewer's copy: htslib takes a control character near the start
# of a file for binary data and will not index it, and a line separator could
# split the line in a reader. They only ever occur in free-text attributes
# (Araport11's note on AT5G64510 has a 0x1E); the stored upload keeps them.
_CONTROL_TO_SPACE = {c: ' ' for c in list(range(0x00, 0x09)) + list(range(0x0a, 0x20))
                     + [0x7f, 0x85, 0x2028, 0x2029]}


def _split_header_and_body(source: Path, header: Path, body: Path,
                           track_type: str) -> None:
    """
    Comment lines to ``header``, feature lines to ``body``.

    ``track``/``browser`` lines are dropped (tabix cannot skip them), as are
    GFF3's ``###`` separators and an embedded ``##FASTA`` section, which mean
    nothing once the features are sorted. Lines end at a newline only, so a
    stray separator character inside an attribute stays in its line, as a
    space (``_CONTROL_TO_SPACE``).
    """
    with open(source, 'r', encoding='utf-8', errors='replace',
              newline='') as inp, \
            open(header, 'w', encoding='utf-8', newline='\n') as head, \
            open(body, 'w', encoding='utf-8', newline='\n') as rows:
        in_fasta = False
        for raw in _newline_lines(inp):
            line = raw.rstrip('\r')
            stripped = line.strip()
            if in_fasta or not stripped:
                continue
            if stripped.startswith('##FASTA'):
                in_fasta = True
                continue
            if stripped.startswith('#'):
                if stripped != '###':
                    head.write(line + '\n')
                continue
            if stripped.startswith('track ') or stripped.startswith('browser '):
                continue
            if '\t' not in line:
                # Whitespace-separated BED: tabix wants tabs.
                line = '\t'.join(line.split())
            rows.write(line.translate(_CONTROL_TO_SPACE) + '\n')


def _newline_lines(handle):
    """A text file's lines split on ``\\n`` only (``newline=''`` keeps ``\\r``)."""
    pending = ''
    while True:
        chunk = handle.read(1024 * 1024)
        if not chunk:
            break
        pending += chunk
        *lines, pending = pending.split('\n')
        yield from lines
    if pending:
        yield pending


# ---------------------------------------------------------------------------
# A region's sequence and its annotations
# ---------------------------------------------------------------------------

# The Sequence Editor holds sequences up to 2 Mb; NCBI's BLAST form takes less,
# which it says itself.
MAX_REGION = 2_000_000
MAX_FEATURES = 5000

_COMPLEMENT = str.maketrans('ACGTRYKMBDHVNacgtrykmbdhvn', 'TGCAYRMKVHDBNtgcayrmkvhdbn')

# Feature types that only repeat their children's extent, left out of what is
# sent: a gene's mRNA and exons say it better, and ``chromosome`` would cover
# everything.
_SKIP_TYPES = {'chromosome', 'region'}


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_region(fasta: Path, chrom: str, start: int, end: int) -> str:
    """The bases of ``[start, end)`` on ``chrom`` (0-based), via the .fai index."""
    ensure_fasta_index(fasta)
    samtools = _need('samtools')
    result = subprocess.run(
        [samtools, 'faidx', str(fasta), f'{chrom}:{start + 1}-{end}'],
        capture_output=True, text=True, check=True)
    return ''.join(line.strip() for line in result.stdout.splitlines()[1:])


def _region_tracks(conn, member: str, genome_id: str, raw: str):
    """The requested feature tracks that are this genome's and the member's to see."""
    out = []
    for track_id in [t.strip() for t in raw.split(',') if t.strip()][:20]:
        track = conn.execute(
            'SELECT id, genome_id, type, owner_lab, visibility, track_status, '
            'file_path FROM tracks WHERE id = ?', (track_id,)).fetchone()
        if (track is None or track['genome_id'] != genome_id
                or track['type'] not in ('bed', 'gff') or not track['file_path']
                or not access.can_view_track(member, track, conn)):
            continue
        out.append(track)
    return out


def _on_region(found, start: int, end: int, strand: int):
    """Features clipped to the region, on its coordinates (and strand)."""
    length = end - start
    out = []
    for _chrom, f_start, f_end, name, f_type, f_strand in found:
        if f_type in _SKIP_TYPES:
            continue
        s, e = max(f_start, start) - start, min(f_end, end) - start
        if e <= s:
            continue
        fs = f_strand or 1
        if strand == -1:
            s, e, fs = length - e, length - s, -fs
        out.append({'name': name or f_type, 'type': _editor_type(f_type),
                    'start': s, 'end': e, 'strand': fs,
                    'note': f_type if _editor_type(f_type) != f_type else ''})
    return out


def _editor_type(gff_type: str) -> str:
    """A GFF feature type as the GenBank-style type the Sequence Editor knows."""
    return {'five_prime_UTR': "5'UTR", 'three_prime_UTR': "3'UTR",
            'transposable_element': 'mobile_element',
            'transposon_fragment': 'mobile_element'}.get(gff_type, gff_type)


def prepare_feature_track(path: Path, track_type: str) -> Optional[str]:
    """
    Make a feature track's viewer copy now; the reason it could not, or None.

    Called at upload so the first view is quick. A missing tool or a failure is
    not an upload error: the first request tries again.
    """
    try:
        ensure_tabix(Path(path), track_type)
        return None
    except (ToolMissing, OSError, subprocess.CalledProcessError) as exc:
        return str(exc)


def prepare_reference(fasta: Path) -> Optional[str]:
    """Index a reference now (on upload); the reason it could not, or None."""
    try:
        ensure_fasta_index(Path(fasta))
        return None
    except (ToolMissing, OSError, subprocess.CalledProcessError) as exc:
        return str(exc)
