"""
Each member's labs and chair groups, read from the suite's own store.

The Genome Browser keeps no lab records of its own. It reads ``suite.db`` —
the sign-in service's database, where ``portal/labs.py`` keeps labs, chair
groups and who is in which — to resolve each member's ``Owner_Lab`` and
``Chair_Group`` and to enforce visibility. The store is always opened
``mode=ro`` through ``db.open_suite``, so this module can never write to it.

The tables mirror ``portal/labs.py``:

- ``labs(id, name, chair_group_id)`` — a lab, and the chair group it is in (if any).
- ``chair_groups(id, name)`` — a group of labs that work together.
- ``lab_members(username, lab_id, role, lab_admin, joined_at)`` — everyone's
  labs, with their role in each. The member's *first* lab (earliest joined) is
  their current ``Owner_Lab``, the lab anything they create belongs to.

Every function tolerates an absent or unreadable store by returning the empty
answer (no labs, no membership), so the module degrades rather than failing
when ``suite.db`` is not mounted — as in development and tests.
"""

import sqlite3
from typing import List, Optional

import db


def _query(sql: str, params: tuple = ()) -> List[sqlite3.Row]:
    """Run a read-only query against ``suite.db``; [] when it is absent or unreadable."""
    conn = db.open_suite()
    if conn is None:
        return []
    try:
        return conn.execute(sql, params).fetchall()
    except sqlite3.Error:
        return []
    finally:
        conn.close()


def labs_of(username: str) -> List[dict]:
    """
    This member's labs, first joined first: ``[{'lab_id', 'lab', 'role'}]``.

    Mirrors ``labs.py``'s ordering (by ``joined_at`` then ``rowid``) so the
    first entry is the member's current lab. A lab lead who also administers
    the lab (``lab_admin``) reads as role ``lab_admin``, as in ``labs.py``.
    """
    username = (username or '').strip()
    if not username:
        return []
    rows = _query(
        'SELECT lab_members.lab_id, lab_members.role, lab_members.lab_admin, '
        'labs.name AS lab FROM lab_members JOIN labs ON labs.id = lab_members.lab_id '
        'WHERE lab_members.username = ? ORDER BY lab_members.joined_at, lab_members.rowid',
        (username,))
    return [{'lab_id': r['lab_id'], 'lab': r['lab'], 'role': _role(r)} for r in rows]


def _role(row: sqlite3.Row) -> str:
    """A lab lead with ``lab_admin`` set reads as 'lab_admin', else the stored role."""
    return 'lab_admin' if row['role'] == 'lab_lead' and row['lab_admin'] else row['role']


def current_lab(username: str) -> Optional[str]:
    """
    The member's first/current lab id — the ``Owner_Lab`` for anything they
    create — or ``None`` when they are in no lab.
    """
    labs = labs_of(username)
    return labs[0]['lab_id'] if labs else None


def role_in(username: str, lab_id: str) -> Optional[str]:
    """Their role in this lab, or ``None`` when they are not in it."""
    return next((m['role'] for m in labs_of(username) if m['lab_id'] == lab_id), None)


def chair_group_of(lab_id: str) -> Optional[str]:
    """
    The chair group this lab is in, or ``None`` when it is in none (or no such lab).
    """
    if not lab_id:
        return None
    rows = _query('SELECT chair_group_id FROM labs WHERE id = ?', (lab_id,))
    return rows[0]['chair_group_id'] if rows else None


def labs_in_chair_group(group_id: str) -> List[str]:
    """Every lab id in this chair group (empty when the group is empty or unknown)."""
    if not group_id:
        return []
    rows = _query('SELECT id FROM labs WHERE chair_group_id = ? ORDER BY name COLLATE NOCASE',
                  (group_id,))
    return [r['id'] for r in rows]


def chair_group_exists(group_id: str) -> bool:
    """Whether a chair group with this id exists in ``suite.db``."""
    if not group_id:
        return False
    return bool(_query('SELECT 1 FROM chair_groups WHERE id = ?', (group_id,)))
