"""
Sharing between labs: who sees a genome or a track, and keeping a shared
track's genome in sight.

A genome is seen by its owner lab's chair group (or the owner lab alone), by
the chair groups and labs it has been shared with, or, as a suite reference, by
everyone. A track is seen by its owner lab and the labs and chair group it has
been shared with (``access.py`` resolves both). A track is only any use to
someone who can also see its genome, so whenever a track's sharing changes this
module brings its genome's sharing into line with it:

* every lab, or chair group, a track on the genome is shared with, and that
  cannot otherwise see the genome, is given the genome too, as a share marked
  ``auto`` (*through a shared track*);
* such an automatic share is removed again once no shared track needs it;
* shares the genome's owner lab made on purpose (``auto = 0``) are never
  touched here, and sharing on purpose something already shared automatically
  makes it a purposeful share.

The automatic share is granted by whoever shared the track, even when the
genome belongs to another lab: sharing a track you are allowed to share would
otherwise show the other lab nothing.

Routes (relative to the ``/genome`` mount), for the page's Share dialog:

- ``GET /api/labs``: every lab and chair group, to choose from
- ``GET /api/genomes/<id>/shares``: who sees a genome, and why
- ``GET /api/tracks/<id>/shares``: who sees a track
"""

from pathlib import Path
from typing import Dict, List, Set

from flask import Flask, current_app, jsonify

import access
import db
import identity
import labs_reader


def register_sharing(app: Flask) -> None:

    @app.route('/api/labs')
    def list_labs():
        member = identity.require_member()
        groups = labs_reader.all_chair_groups()
        names = {g['id']: g['name'] for g in groups}
        labs = [dict(lab, chair_group=names.get(lab['chair_group_id']))
                for lab in labs_reader.all_labs()]
        if labs_reader.in_demo(member):
            # The teaching demo is shown its own labs and chair group, not the real ones.
            demo = labs_reader.demo_lab_ids()
            labs = [lab for lab in labs if lab['id'] in demo]
            mine = {lab['chair_group_id'] for lab in labs}
            groups = [g for g in groups if g['id'] in mine]
        return jsonify({'labs': labs, 'chair_groups': groups})

    @app.route('/api/genomes/<genome_id>/shares')
    def genome_shares(genome_id):
        member = identity.require_member()
        conn = db.connect(_db_path())
        try:
            genome = conn.execute(
                'SELECT id, name, owner_lab, genome_visibility FROM reference_genomes '
                'WHERE id = ?', (genome_id,)).fetchone()
            if genome is None or not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'Reference genome not found'}), 404
            return jsonify(genome_share_body(conn, member, genome))
        finally:
            conn.close()

    @app.route('/api/tracks/<track_id>/shares')
    def track_shares(track_id):
        member = identity.require_member()
        conn = db.connect(_db_path())
        try:
            track = conn.execute(
                'SELECT id, genome_id, name, owner_lab, visibility, track_status '
                'FROM tracks WHERE id = ?', (track_id,)).fetchone()
            if track is None or not (access.can_view_track(member, track, conn)
                                     or access.is_owner_lab_member(member, track['owner_lab'])):
                return jsonify({'error': 'Track not found'}), 404
            return jsonify(track_share_body(conn, member, track))
        finally:
            conn.close()


def _db_path() -> Path:
    return Path(current_app.config['DATA_ROOT']) / 'genome.db'


def _named(kind: str, ident: str) -> Dict:
    name = (labs_reader.lab_name(ident) if kind == 'lab'
            else labs_reader.chair_group_name(ident))
    return {'type': kind, 'id': ident, 'name': name or ident}


def genome_share_body(conn, member: str, genome) -> Dict:
    """Who sees a genome: its default viewers, and each share with whether it is automatic."""
    suite = genome['genome_visibility'] == access.SUITE_VISIBILITY
    owner = genome['owner_lab']
    default = []
    if suite:
        default.append({'type': 'everyone', 'name': 'Everyone in the suite'})
    else:
        default.append(dict(_named('lab', owner), role='owner'))
        group = labs_reader.chair_group_of(owner)
        if group and genome['genome_visibility'] == 'chair_group_default':
            # Removable: taking it away narrows the genome to its owner lab (genomes.py).
            default.append(dict(_named('chair_group', group), role='default'))
    shares = []
    for row in conn.execute('SELECT chair_group_id, auto FROM genome_shares WHERE genome_id = ? '
                            'ORDER BY granted_at', (genome['id'],)):
        shares.append(dict(_named('chair_group', row['chair_group_id']), auto=bool(row['auto'])))
    for row in conn.execute('SELECT lab_id, auto FROM genome_lab_shares WHERE genome_id = ? '
                            'ORDER BY granted_at', (genome['id'],)):
        shares.append(dict(_named('lab', row['lab_id']), auto=bool(row['auto'])))
    return {
        'id': genome['id'], 'name': genome['name'], 'suite': suite,
        'owner_lab': _named('lab', owner) if not suite else None,
        'can_change': (not suite) and access.is_owner_lab_member(member, owner),
        'default': default, 'shares': shares,
        # The owner lab's chair group, which sees it by default: removable, and addable again.
        'owner_chair_group': (_named('chair_group', labs_reader.chair_group_of(owner))
                              if not suite and labs_reader.chair_group_of(owner) else None),
    }


def track_share_body(conn, member: str, track) -> Dict:
    """Who sees a track: its owner lab and the labs and chair group it is shared with."""
    owner = track['owner_lab']
    suite = owner == access.SUITE_OWNER
    shares = [_named(row['target_type'], row['target_id']) for row in conn.execute(
        'SELECT target_type, target_id FROM track_shares WHERE track_id = ? ORDER BY granted_at',
        (track['id'],))]
    group = None if suite else labs_reader.chair_group_of(owner)
    return {
        'id': track['id'], 'name': track['name'], 'genome_id': track['genome_id'],
        'suite': suite,
        'owner_lab': _named('lab', owner) if not suite else None,
        'owner_chair_group': _named('chair_group', group) if group else None,
        'can_change': (not suite) and access.is_owner_lab_member(member, owner),
        'shares': shares,
    }


# ---------------------------------------------------------------------------
# Keeping a shared track's genome in sight
# ---------------------------------------------------------------------------


def sync_genome_access(conn, genome_id: str, member: str) -> List[Dict]:
    """
    Bring a genome's automatic shares into line with its tracks' sharing.

    Returns the shares this call added, as ``{type, id, name}``, so the caller
    can tell the member that the genome was shared too. The caller commits.
    """
    genome = conn.execute(
        'SELECT id, owner_lab, genome_visibility FROM reference_genomes WHERE id = ?',
        (genome_id,)).fetchone()
    if genome is None:
        return []
    if genome['genome_visibility'] == access.SUITE_VISIBILITY:
        # Everyone sees a suite reference already: nothing automatic is needed.
        conn.execute('DELETE FROM genome_shares WHERE genome_id = ? AND auto = 1', (genome_id,))
        conn.execute('DELETE FROM genome_lab_shares WHERE genome_id = ? AND auto = 1', (genome_id,))
        return []

    seen = _viewers_without_auto(conn, genome)
    need_labs: Set[str] = set()
    need_groups: Set[str] = set()
    for row in conn.execute(
            'SELECT s.target_type, s.target_id FROM track_shares s '
            'JOIN tracks t ON t.id = s.track_id WHERE t.genome_id = ?', (genome_id,)):
        if row['target_type'] == 'lab':
            if row['target_id'] not in seen:
                need_labs.add(row['target_id'])
        elif row['target_type'] == 'chair_group':
            labs = set(labs_reader.labs_in_chair_group(row['target_id']))
            if labs and not labs <= seen:
                need_groups.add(row['target_id'])

    added: List[Dict] = []
    have_labs = {r['lab_id']: r['auto'] for r in conn.execute(
        'SELECT lab_id, auto FROM genome_lab_shares WHERE genome_id = ?', (genome_id,))}
    have_groups = {r['chair_group_id']: r['auto'] for r in conn.execute(
        'SELECT chair_group_id, auto FROM genome_shares WHERE genome_id = ?', (genome_id,))}
    for lab in sorted(need_labs - set(have_labs)):
        conn.execute('INSERT INTO genome_lab_shares (genome_id, lab_id, granted_by, auto) '
                     'VALUES (?, ?, ?, 1)', (genome_id, lab, member))
        added.append(_named('lab', lab))
    for group in sorted(need_groups - set(have_groups)):
        conn.execute('INSERT INTO genome_shares (genome_id, chair_group_id, granted_by, auto) '
                     'VALUES (?, ?, ?, 1)', (genome_id, group, member))
        added.append(_named('chair_group', group))
    for lab, auto in have_labs.items():
        if auto and lab not in need_labs:
            conn.execute('DELETE FROM genome_lab_shares WHERE genome_id = ? AND lab_id = ? AND auto = 1',
                         (genome_id, lab))
    for group, auto in have_groups.items():
        if auto and group not in need_groups:
            conn.execute('DELETE FROM genome_shares WHERE genome_id = ? AND chair_group_id = ? AND auto = 1',
                         (genome_id, group))
    return added


def _viewers_without_auto(conn, genome) -> Set[str]:
    """The labs that see a genome by its default and its purposeful shares alone."""
    labs: Set[str] = set()
    owner = genome['owner_lab']
    labs.add(owner)
    if genome['genome_visibility'] == 'chair_group_default':
        group = labs_reader.chair_group_of(owner)
        if group:
            labs.update(labs_reader.labs_in_chair_group(group))
    for row in conn.execute('SELECT chair_group_id FROM genome_shares WHERE genome_id = ? AND auto = 0',
                            (genome['id'],)):
        labs.update(labs_reader.labs_in_chair_group(row['chair_group_id']))
    for row in conn.execute('SELECT lab_id FROM genome_lab_shares WHERE genome_id = ? AND auto = 0',
                            (genome['id'],)):
        labs.add(row['lab_id'])
    return labs
