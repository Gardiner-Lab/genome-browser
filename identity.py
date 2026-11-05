"""
Who is making the request, and who is an administrator.

The Genome Browser keeps no accounts of its own. Behind the proxy it trusts
the ``X-Suite-User`` header, which nginx sets after the portal's ``/api/check``
confirms the session; the header names the signed-in member. Administrator
status lives in the shared account directory — MoClo's ``users.is_admin``
column — exactly as ``portal/directory.py`` reads it, so every tool agrees on
who an administrator is and this module never keeps a copy.

Both reads are read-only: ``is_admin`` opens ``moclo.db`` ``mode=ro`` through
``db.open_moclo``, so this module can never write to another tool's store.
"""

import sqlite3
from typing import Optional

from flask import request

import db

# The header the proxy sets to the signed-in member (see the design's
# "Identity is proxy-trusted" convention). It is the only source of identity.
SUITE_USER_HEADER = 'X-Suite-User'


class Unauthenticated(Exception):
    """
    No signed-in member on the request.

    Raised by :func:`require_member` so an endpoint can let it turn into a 401;
    behind the proxy this should not happen for a protected route, because
    nginx only passes a request on once ``/api/check`` has set the header.
    """


def current_member() -> Optional[str]:
    """
    The signed-in member from ``X-Suite-User``, or ``None`` when it is absent.

    Reads the header the proxy sets; surrounding whitespace is stripped and an
    empty header reads as absent. Callers that require a member use
    :func:`require_member`, which turns the absent case into 401 behaviour.
    """
    member = (request.headers.get(SUITE_USER_HEADER) or '').strip()
    return member or None


def require_member() -> str:
    """
    The signed-in member, or :class:`Unauthenticated` when the header is absent.

    The 401 behaviour the design calls for: an endpoint that needs a member
    lets this raise, and the app turns it into an HTTP 401. The proxy normally
    guarantees the header, so this guards the unproxied/misconfigured case.
    """
    member = current_member()
    if member is None:
        raise Unauthenticated('No signed-in member on the request')
    return member


def is_admin(username: str) -> bool:
    """
    Whether this member is a suite administrator, from MoClo's ``users.is_admin``.

    Reads the shared directory read-only, mirroring ``portal/directory.py``: an
    older MoClo schema without the ``is_admin`` column, an unreachable store, a
    missing account, or an empty username all read as not an administrator.
    """
    username = (username or '').strip()
    if not username:
        return False
    conn = db.open_moclo()
    if conn is None:
        return False
    try:
        row = conn.execute('SELECT * FROM users WHERE username = ?', (username,)).fetchone()
    except sqlite3.Error:
        return False
    finally:
        conn.close()
    if row is None:
        return False
    # Older MoClo releases have no is_admin column; absent reads as not an admin.
    keys = row.keys()
    return bool(row['is_admin']) if 'is_admin' in keys else False
