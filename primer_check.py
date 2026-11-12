"""
Primers checked against a reference genome with BLAST: does each bind the genome at one place only?

For the Sequence Editor's bisulfite primers, and any other primer: the editor's page sends the
primers (as their sites read in the genome, before conversion) to
``POST /api/genomes/<id>/primer-check``. The test is the one in ``design_bisulfite_primers.py`` by
Colette L. Picard (2016), which BLASTed each primer against TAIR10: a primer is turned down when a
hit other than its own site, a perfect match over its whole length, covers more than
``max_match`` of its length at more than ``min_identity`` percent identity. The script judged the
second line of BLAST's output, which lists every hit on the best-matching chromosome first, so in
effect it searched only the primer's own chromosome; here every hit on every chromosome counts.
That is far stricter at the script's cutoffs (75% of the primer, 80% identity), under which most
primers in AT-rich plant DNA have a near-match somewhere, so the defaults are 90% and 90%. A
primer found nowhere passes, with no site, as the script let it (it may be from another genome,
or a plasmid).

BLAST needs a database made from the genome's sequence (makeblastdb), a minute for Arabidopsis
and longer for a large genome, so it is made once per genome by the worker, the first time that
genome is checked, and kept under ``<data>/blastdb/<genome id>/``. Until it is ready the endpoint
answers 202, and the page asks again.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

BLASTN = 'blastn'
MAKEBLASTDB = 'makeblastdb'
MAX_PRIMERS = 200
PRIMER_LENGTH = (10, 60)
BLAST_TIMEOUT = 170          # under gunicorn's 180 s
BUILD_TIMEOUT = 6 * 3600
FIELDS = ('qseqid', 'stitle', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')


class NotReady(Exception):
    """The genome's BLAST database is being made: ask again in a little while."""


def _root(data_root: Path) -> Path:
    return Path(data_root) / 'blastdb'


def db_dir(data_root: Path, genome_id: str) -> Path:
    return _root(data_root) / genome_id


def _wanted(data_root: Path) -> Path:
    return _root(data_root) / 'wanted'


def _reference(data_root: Path, genome_id: str) -> Path:
    return Path(data_root) / 'genomes' / genome_id / 'reference.fasta'


def state(data_root: Path, genome_id: str) -> Dict[str, Any]:
    """``ready``, ``building``, ``failed`` (with why) or ``none`` for this genome's database."""
    folder = db_dir(data_root, genome_id)
    ready = folder / 'ready'
    reference = _reference(data_root, genome_id)
    if ready.exists() and (not reference.exists() or ready.stat().st_mtime >= reference.stat().st_mtime):
        return {'state': 'ready'}
    failed = folder.with_name(f'{genome_id}.failed')
    if failed.exists():
        return {'state': 'failed', 'error': failed.read_text(encoding='utf-8', errors='replace').strip()}
    if (_wanted(data_root) / genome_id).exists():
        return {'state': 'building'}
    return {'state': 'none'}


def request_build(data_root: Path, genome_id: str) -> None:
    """Ask the worker to make this genome's database (again, if the last try failed)."""
    wanted = _wanted(data_root)
    wanted.mkdir(parents=True, exist_ok=True)
    db_dir(data_root, genome_id).with_name(f'{genome_id}.failed').unlink(missing_ok=True)
    (wanted / genome_id).touch()


def build_next(data_root: Path) -> bool:
    """For the worker: make one database that was asked for. Gives whether there was one to make."""
    wanted = _wanted(data_root)
    try:
        asked = sorted(wanted.iterdir(), key=lambda p: p.stat().st_mtime)
    except FileNotFoundError:
        return False
    if not asked:
        return False
    marker = asked[0]
    genome_id = marker.name
    reference = _reference(data_root, genome_id)
    folder = db_dir(data_root, genome_id)
    failed = folder.with_name(f'{genome_id}.failed')
    work = folder.with_name(f'{genome_id}.building')
    print(f'[genome-worker] making the BLAST database of genome {genome_id}.', flush=True)
    try:
        if not reference.exists():
            raise RuntimeError("the genome's sequence is not there")
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        done = subprocess.run([MAKEBLASTDB, '-in', str(reference), '-dbtype', 'nucl', '-out', str(work / 'reference'),
                               '-title', genome_id], capture_output=True, text=True, timeout=BUILD_TIMEOUT)
        if done.returncode != 0:
            raise RuntimeError((done.stderr or done.stdout or 'makeblastdb failed').strip()[-500:])
        (work / 'ready').touch()
        shutil.rmtree(folder, ignore_errors=True)
        os.replace(work, folder)
        print(f'[genome-worker] BLAST database of genome {genome_id} ready.', flush=True)
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        shutil.rmtree(work, ignore_errors=True)
        failed.write_text(str(exc), encoding='utf-8')
        print(f'[genome-worker] could not make the BLAST database of genome {genome_id}: {exc}', flush=True)
    finally:
        marker.unlink(missing_ok=True)
    return True


def clean_primers(raw: Any) -> List[Dict[str, str]]:
    """``[{id, sequence}]``, each an unambiguous DNA sequence of a primer's length."""
    if not isinstance(raw, list) or not raw:
        raise ValueError('Send the primers to check as a list of {id, sequence}.')
    if len(raw) > MAX_PRIMERS:
        raise ValueError(f'Check at most {MAX_PRIMERS} primers at a time.')
    out, seen = [], set()
    for i, p in enumerate(raw):
        if not isinstance(p, dict):
            raise ValueError('Each primer is {id, sequence}.')
        pid = str(p.get('id', i)).strip()[:40] or str(i)
        seq = ''.join(str(p.get('sequence') or '').split()).upper()
        if pid in seen:
            raise ValueError(f'Primer id {pid} is given twice.')
        if not PRIMER_LENGTH[0] <= len(seq) <= PRIMER_LENGTH[1] or set(seq) - set('ACGT'):
            raise ValueError(f'Primer {pid}: {PRIMER_LENGTH[0]}–{PRIMER_LENGTH[1]} nt of A, C, G and T, as the site '
                             'reads in the genome.')
        seen.add(pid)
        out.append({'id': pid, 'sequence': seq})
    return out


def _hit(row: Dict[str, Any]) -> Dict[str, Any]:
    start, end = sorted((row['sstart'], row['send']))
    return {'chrom': row['chrom'], 'start': start, 'end': end, 'strand': '+' if row['sstart'] <= row['send'] else '−',
            'identity': round(row['pident'], 1), 'length': row['length']}


def verdict(primer: Dict[str, str], hits: List[Dict[str, Any]], max_match: float, min_identity: float) -> Dict[str, Any]:
    """Whether one primer binds once only, from its BLAST hits."""
    n = len(primer['sequence'])
    hits = sorted(hits, key=lambda h: (-h['bitscore'], h['evalue']))
    site = next((h for h in hits if h['pident'] == 100.0 and h['length'] == n and h['qstart'] == 1 and h['qend'] == n), None)
    others = [h for h in hits if h is not site]
    offending = [h for h in others if h['length'] > max_match * n and h['pident'] > min_identity]
    other = (offending or others or [None])[0]
    return {'id': primer['id'], 'ok': not offending, 'hits': len(hits), 'off_target': len(offending),
            'site': _hit(site) if site else None, 'other': _hit(other) if other else None}


def parse(output: str) -> Dict[str, List[Dict[str, Any]]]:
    by_query: Dict[str, List[Dict[str, Any]]] = {}
    for line in output.splitlines():
        cols = line.split('\t')
        if len(cols) != len(FIELDS):
            continue
        row = dict(zip(FIELDS, cols))
        by_query.setdefault(row['qseqid'], []).append({
            'chrom': (row['stitle'].split() or [''])[0], 'pident': float(row['pident']), 'length': int(row['length']),
            'qstart': int(row['qstart']), 'qend': int(row['qend']), 'sstart': int(row['sstart']),
            'send': int(row['send']), 'evalue': float(row['evalue']), 'bitscore': float(row['bitscore'])})
    return by_query


def check(data_root: Path, genome_id: str, primers: List[Dict[str, str]], max_match: float = 0.9,
          min_identity: float = 90.0) -> List[Dict[str, Any]]:
    """A verdict for each primer; :class:`NotReady` while the genome's database is being made."""
    now = state(data_root, genome_id)
    if now['state'] == 'failed':
        raise RuntimeError(f"This genome could not be prepared for searching: {now['error']}")
    if now['state'] != 'ready':
        if now['state'] == 'none':
            request_build(data_root, genome_id)
        raise NotReady()
    # BLAST names queries by the first word of their FASTA line: number them, and map back.
    query = ''.join(f'>q{i}\n{p["sequence"]}\n' for i, p in enumerate(primers))
    done = subprocess.run(
        [BLASTN, '-task', 'blastn', '-db', str(db_dir(data_root, genome_id) / 'reference'),
         '-outfmt', '6 ' + ' '.join(FIELDS), '-max_target_seqs', '100', '-num_threads', '2'],
        input=query, capture_output=True, text=True, timeout=BLAST_TIMEOUT)
    if done.returncode != 0:
        raise RuntimeError(f"BLAST could not search this genome: {(done.stderr or '').strip()[-300:]}")
    hits = parse(done.stdout)
    return [verdict(p, hits.get(f'q{i}', []), max_match, min_identity) for i, p in enumerate(primers)]


def register_primer_check(app) -> None:
    from flask import jsonify, request

    import access
    import db
    import identity

    @app.route('/api/genomes/<genome_id>/primer-check', methods=['POST'])
    def primer_check(genome_id):
        """
        ``{"primers": [{"id": "7", "sequence": "TTGCC..."}], "max_match": 0.9, "min_identity": 90}``:
        ``{"primers": [{id, ok, hits, off_target, site, other}]}``, or 202 while the genome is being
        prepared for searching (the first time it is checked).
        """
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401
        body = request.get_json(silent=True) or {}
        try:
            primers = clean_primers(body.get('primers'))
            max_match = float(body.get('max_match', 0.9))
            min_identity = float(body.get('min_identity', 90))
            if not 0 < max_match <= 1 or not 0 <= min_identity <= 100:
                raise ValueError('max_match is a fraction of the primer (0–1), min_identity a percentage.')
        except (TypeError, ValueError) as exc:
            return jsonify({'error': str(exc)}), 400
        data_root = Path(app.config['DATA_ROOT'])
        conn = db.connect(data_root / 'genome.db')
        try:
            genome = conn.execute('SELECT id, name, assembly_id, owner_lab, genome_visibility '
                                  'FROM reference_genomes WHERE id = ?', (genome_id,)).fetchone()
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            name = genome['name']
        finally:
            conn.close()
        if not _reference(data_root, genome_id).exists():
            return jsonify({'error': "This genome's sequence is not ready yet."}), 503
        try:
            verdicts = check(data_root, genome_id, primers, max_match, min_identity)
        except NotReady:
            return jsonify({'state': 'building', 'message': f'Preparing {name} for searching: the first time, a minute '
                            'or so for Arabidopsis, longer for a large genome…'}), 202
        except (RuntimeError, OSError, subprocess.SubprocessError) as exc:
            return jsonify({'error': str(exc)}), 503
        return jsonify({'genome': name, 'primers': verdicts})
