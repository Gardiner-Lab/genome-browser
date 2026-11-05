"""
Who may view a genome or a track, and who may change their sharing.

The Genome Browser has two deliberately different sharing schemes, and this
module is the single place both are resolved, over the module's own database
(``genome_shares``/``track_shares``) and the suite's lab/chair-group reader
(``labs_reader``). Endpoints never reinvent the rules; they call these helpers.

- **Genome visibility (Requirement 12).** A reference genome defaults to its
  Owner_Lab's *chair group*: every lab in that chair group may view it. WHERE
  the Owner_Lab belongs to no chair group, the default is the Owner_Lab only.
  That default is recorded as ``genome_visibility`` —
  ``'chair_group_default'`` or ``'owner_lab_only'``. On top of the default a
  member of the Owner_Lab may share the genome with *whole chair groups*,
  recorded in ``genome_shares``; the labs of every shared chair group may then
  view it too. Sharing a genome with an individual lab is never allowed (that
  rejection lives in the sharing endpoint, task 5.3).

- **Track visibility (Requirement 7).** A track defaults to ``lab_only`` — only
  the Owner_Lab's members may view it. It can be shared with one or more
  *specific labs* and/or the Owner_Lab's *chair group*, recorded as rows in
  ``track_shares`` (``target_type`` ``'lab'`` or ``'chair_group'``). A member
  may view a track when their lab is the Owner_Lab, a shared lab, or a lab in a
  shared chair group — **and only when ``track_status='ready'``** (Requirement
  13.6): a ``converting`` or ``failed`` track is viewable by no one.

Both resolvers answer in terms of the member's *labs* (a member may be in
more than one), so a member sees a genome or track when any of their labs is
in the viewable set.

- **Owner-lab guard (Requirements 7.6, 12.6).** Every sharing mutation — the
  genome and track share/unshare endpoints (tasks 5.3, 6.2) — calls
  :func:`require_owner_lab` first. It succeeds only when the member belongs to
  the Owner_Lab, and otherwise raises :class:`PermissionDenied` so the mutation
  changes nothing and the endpoint returns a 403 with a clear message.

Every lab/chair-group lookup goes through ``labs_reader``, which reads
``suite.db`` read-only and degrades to the empty answer when the store is
absent; this module adds no lab records of its own.
"""

from typing import Iterable, List, Mapping, Optional, Set

import db
import labs_reader


class PermissionDenied(Exception):
    """
    A member tried to change sharing they do not own.

    Raised by :func:`require_owner_lab` so a sharing endpoint can let it become
    an HTTP 403 while the mutation leaves every share row unchanged
    (Requirements 7.6, 12.6).
    """


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _get(row, key, default=None):
    """
    Read ``key`` from a genome/track record that may be a sqlite Row or a dict.

    Endpoints hand these resolvers whatever they already have — a ``sqlite3.Row``
    straight from a query or a plain dict in a test — so both are accepted.
    """
    if row is None:
        return default
    if isinstance(row, Mapping):
        return row.get(key, default)
    try:
        keys = row.keys()
    except AttributeError:
        return getattr(row, key, default)
    return row[key] if key in keys else default


def _member_lab_ids(username: str) -> Set[str]:
    """Every lab id this member belongs to (a member may be in more than one)."""
    return {m['lab_id'] for m in labs_reader.labs_of(username)}


def _labs_of_chair_groups(group_ids: Iterable[str]) -> Set[str]:
    """The union of every lab id in any of the given chair groups."""
    labs: Set[str] = set()
    for group_id in group_ids:
        if group_id:
            labs.update(labs_reader.labs_in_chair_group(group_id))
    return labs


# ---------------------------------------------------------------------------
# Genome visibility (Requirement 12)
# ---------------------------------------------------------------------------

def genome_viewer_labs(genome, conn=None) -> Set[str]:
    """
    Every lab id permitted to view this genome.

    The default viewers come from ``genome_visibility``: when it is
    ``'chair_group_default'`` and the Owner_Lab is in a chair group, that is
    every lab in the Owner_Lab's chair group; otherwise (``'owner_lab_only'``,
    or a ``chair_group_default`` genome whose Owner_Lab is in no chair group)
    it is the Owner_Lab alone. Added to that are the labs of every chair group
    the genome has been shared with (``genome_shares``). (Requirements 12.1,
    12.2, 12.3.)
    """
    owner_lab = _get(genome, 'owner_lab')
    visibility = _get(genome, 'genome_visibility', 'chair_group_default')

    labs: Set[str] = set()
    if owner_lab:
        labs.add(owner_lab)
        if visibility == 'chair_group_default':
            owner_group = labs_reader.chair_group_of(owner_lab)
            if owner_group:
                labs.update(labs_reader.labs_in_chair_group(owner_group))

    genome_id = _get(genome, 'id')
    if genome_id:
        labs.update(_labs_of_chair_groups(_shared_chair_groups(genome_id, conn)))
    return labs


def _shared_chair_groups(genome_id: str, conn=None) -> List[str]:
    """The chair group ids a genome has been shared with (``genome_shares``)."""
    owns = conn is None
    conn = conn or db.connect()
    try:
        rows = conn.execute(
            'SELECT chair_group_id FROM genome_shares WHERE genome_id = ?',
            (genome_id,)).fetchall()
    finally:
        if owns:
            conn.close()
    return [r['chair_group_id'] for r in rows]


def can_view_genome(username: str, genome, conn=None) -> bool:
    """
    Whether this member may view the genome.

    True when any of the member's labs is in the genome's viewer set (its
    default viewers plus any shared chair groups). (Requirements 12.1, 12.2,
    12.3.)
    """
    username = (username or '').strip()
    if not username or genome is None:
        return False
    return bool(_member_lab_ids(username) & genome_viewer_labs(genome, conn))


# ---------------------------------------------------------------------------
# Track visibility (Requirement 7, gated by 13.6)
# ---------------------------------------------------------------------------

def track_viewer_labs(track, conn=None) -> Set[str]:
    """
    Every lab id permitted to view this track, ignoring Track_Status.

    Always the Owner_Lab (``lab_only`` default), plus every lab and every lab
    of every chair group the track has been shared with (``track_shares``).
    :func:`can_view_track` additionally requires ``track_status='ready'``; this
    resolver answers the lab side only so sharing logic can reuse it.
    (Requirements 7.2, 7.3, 7.4.)
    """
    labs: Set[str] = set()
    owner_lab = _get(track, 'owner_lab')
    if owner_lab:
        labs.add(owner_lab)

    track_id = _get(track, 'id')
    if track_id:
        shared_labs, shared_groups = _track_share_targets(track_id, conn)
        labs.update(shared_labs)
        labs.update(_labs_of_chair_groups(shared_groups))
    return labs


def _track_share_targets(track_id: str, conn=None):
    """The (lab ids, chair group ids) a track has been shared with."""
    owns = conn is None
    conn = conn or db.connect()
    try:
        rows = conn.execute(
            'SELECT target_type, target_id FROM track_shares WHERE track_id = ?',
            (track_id,)).fetchall()
    finally:
        if owns:
            conn.close()
    labs = [r['target_id'] for r in rows if r['target_type'] == 'lab']
    groups = [r['target_id'] for r in rows if r['target_type'] == 'chair_group']
    return labs, groups


def can_view_track(username: str, track, conn=None) -> bool:
    """
    Whether this member may view the track.

    A track is viewable only while its ``track_status`` is ``'ready'``: a
    ``converting`` or ``failed`` track is viewable by no one (Requirement
    13.6). When ready, the member may view it if any of their labs is the
    Owner_Lab, a shared lab, or a lab in a shared chair group (Requirements
    7.1, 7.2, 7.3, 7.4).
    """
    username = (username or '').strip()
    if not username or track is None:
        return False
    if _get(track, 'track_status', 'ready') != 'ready':
        return False
    return bool(_member_lab_ids(username) & track_viewer_labs(track, conn))


# ---------------------------------------------------------------------------
# Owner-lab guard (Requirements 7.6, 12.6)
# ---------------------------------------------------------------------------

def is_owner_lab_member(username: str, owner_lab: str) -> bool:
    """Whether this member belongs to the Owner_Lab — the one lab that may share."""
    username = (username or '').strip()
    if not username or not owner_lab:
        return False
    return owner_lab in _member_lab_ids(username)


def require_owner_lab(username: str, owner_lab: str) -> None:
    """
    Guard every sharing mutation: the member must belong to the Owner_Lab.

    Returns quietly when the member is in the Owner_Lab; otherwise raises
    :class:`PermissionDenied` so the caller changes no share rows and returns a
    403 with a message that the member may not change this sharing.
    (Requirements 7.6, 12.6.)
    """
    if not is_owner_lab_member(username, owner_lab):
        raise PermissionDenied(
            'You are not permitted to change this sharing; '
            'only members of the owner lab may.')
