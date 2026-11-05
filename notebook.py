"""
The direct-to-notebook SVG-export endpoint: ``POST /api/figure/notebook``.

Where :mod:`figures` lets a member *download* the SVG of their current view,
this module sends that same SVG *straight into the member's Lab Notebook* so
they can record a figure with an experiment without a download/re-upload round
trip (Requirement 10). It is deliberately thin: it reuses the exact
view-assembly and pure SVG builder that ``POST /api/figure.svg`` uses
(:func:`figures._build_view` + :func:`render.build_svg`), so the figure a
member files in their notebook is byte-for-byte the one the download endpoint
would have produced. The only new behaviour here is the HTTP call out to the
Lab Notebook and the failure handling around it.

What the endpoint does (Requirement 10):

- **10.1** It POSTs the SVG to ``NOTEBOOK_URL`` *as the current member*,
  forwarding the ``X-Suite-User`` header so the Lab Notebook creates the entry
  under the member's own account rather than the service's.
- **10.4** The payload carries the ``Reference_Genome`` name and the visible
  ``Locus`` alongside the SVG, so the notebook entry records *what* the figure
  shows.
- **10.2** When the Lab Notebook accepts the figure, the endpoint returns a
  reference to the created entry (its id and, when the notebook returns one,
  its URL) so the member can jump straight to it.
- **10.3** When the Lab Notebook is unreachable or rejects the request, the
  endpoint returns a descriptive error *and* keeps the figure downloadable: the
  SVG itself is returned in the error body (plus the ``/api/figure.svg`` route
  the member can re-request), so a notebook outage never loses the figure.

The outbound HTTP call is isolated behind :func:`post_to_notebook` so the
property/integration tests (task 13.2) can stub a success, a rejection, and an
unreachable notebook without a live service.

EXTERNAL HANDOFF
----------------
There is no Lab Notebook SVG-ingest HTTP API in this repository today (the only
``notebook`` artefacts are UI glyphs, not a service with an entry-creation
endpoint). This module therefore implements the *client* side against an
assumed ingest contract and records the server side as an external handoff:

    The Lab Notebook must expose an SVG-ingest endpoint at
    ``POST {NOTEBOOK_URL}/api/entries`` that

      * trusts the forwarded ``X-Suite-User`` header as the entry's author
        (same proxy-identity convention the rest of the suite uses),
      * accepts a JSON body ``{genome, locus, filename, content_type, svg}``
        where ``svg`` is the figure document as text, and
      * returns JSON identifying the created entry, e.g.
        ``{"id": "...", "url": "/notebook/entries/..."}`` (any of ``id`` /
        ``entry_id`` and ``url`` / ``href`` are accepted here).

If the Notebook team settles on a different path or envelope, only
:func:`post_to_notebook` and :func:`_entry_reference` need to change; the
endpoint and the failure contract (10.3) stay as they are.

Kept in its own module with a :func:`register_notebook` registrar (mirroring
``figures.register_figures``) so the app factory wires it with a single line.

Requirements: 10.1 (SVG sent as the member via X-Suite-User), 10.2 (reference
to the created entry on success), 10.3 (descriptive error + downloadable SVG on
failure), 10.4 (genome name + locus included in the payload).
"""

from typing import Dict, Optional, Tuple

import requests
from flask import Response, jsonify, request

import access
import db
import figures
import identity
import render

# How the module reaches the Lab Notebook's SVG-ingest endpoint. See the
# EXTERNAL HANDOFF note above: the client assumes an ``/api/entries`` ingest
# path under NOTEBOOK_URL. A short timeout keeps a slow/hung notebook from
# wedging a request; a timeout is treated as "unreachable" (Requirement 10.3).
_NOTEBOOK_INGEST_PATH = '/api/entries'
_NOTEBOOK_TIMEOUT_SECONDS = 10


def register_notebook(app) -> None:
    """
    Register ``POST /api/figure/notebook`` on the app (called by the factory).

    One registrar, mirroring ``register_figures``, so the app factory wires
    this with a single line and this module owns its own route.
    """

    @app.route('/api/figure/notebook', methods=['POST'])
    def figure_to_notebook():
        try:
            member = identity.require_member()
        except identity.Unauthenticated:
            return jsonify({'error': 'Sign-in required.'}), 401

        payload = request.get_json(silent=True) or {}
        genome_id = (payload.get('genome_id') or '').strip()
        if not genome_id:
            return jsonify({'error': 'A genome_id is required.'}), 400

        # Build exactly the SVG the download endpoint would produce, reusing
        # figures' view assembly so the notebook gets the same figure.
        conn = db.connect(app.config['DATA_ROOT'] / 'genome.db')
        try:
            genome = figures._genome_row(conn, genome_id)
            if genome is None:
                return jsonify({'error': 'Genome not found.'}), 404
            if not access.can_view_genome(member, genome, conn):
                return jsonify({'error': 'You may not view this genome.'}), 403

            view = figures._build_view(conn, member, genome, payload)
            genome_name = genome['name'] if figures._has(genome, 'name') else ''
        finally:
            conn.close()

        svg = render.build_svg(view.spec)
        locus_text = _locus_text(view.spec.locus)

        notebook_url = app.config.get('NOTEBOOK_URL', '')
        try:
            reference = post_to_notebook(
                notebook_url, member, svg, genome_name, locus_text)
        except NotebookUnavailable as exc:
            # 10.3: the notebook could not take the figure. Return a
            # descriptive error AND keep the figure downloadable — the SVG
            # itself travels in the body so nothing is lost, and the
            # download route is named so the member can re-request it.
            body = {
                'error': ('Could not send the figure to your Lab Notebook: '
                          f'{exc}. The figure is still available to download.'),
                'notebook_delivered': False,
                'svg': svg,
                'download_url': '/api/figure.svg',
                'genome': genome_name,
                'locus': locus_text,
            }
            resp = jsonify(body)
            resp.status_code = 502
            return resp

        # 10.2: the notebook accepted the figure — return a reference to the
        # entry it created so the member can open it.
        body = {
            'notebook_delivered': True,
            'entry': reference,
            'genome': genome_name,
            'locus': locus_text,
        }
        if view.colour_errors:
            body['colour_warning'] = ' '.join(view.colour_errors)[:400]
        return jsonify(body), 201


# ---------------------------------------------------------------------------
# The isolated HTTP call to the Lab Notebook (stubbed by tests, task 13.2)
# ---------------------------------------------------------------------------


class NotebookUnavailable(Exception):
    """
    The Lab Notebook could not accept the figure.

    Raised for any failure to deliver — the notebook unreachable, a timeout, or
    a non-success HTTP status — so the endpoint can turn it into the
    descriptive-error-plus-downloadable-SVG response (Requirement 10.3). The
    message is suitable to show the member.
    """


def post_to_notebook(notebook_url: str, member: str, svg: str,
                     genome_name: str, locus_text: str) -> Dict:
    """
    POST the SVG figure to the Lab Notebook as ``member``; return the entry ref.

    The one place the module talks to the Lab Notebook, kept small and pure of
    Flask so tests can stub it. It forwards the member's identity in the
    ``X-Suite-User`` header (Requirement 10.1) so the notebook files the entry
    under the member's own account, and includes the ``Reference_Genome`` name
    and the visible ``Locus`` in the JSON body (Requirement 10.4) so the entry
    records what the figure shows.

    Returns a reference to the created entry (``{'id', 'url'}``) on success
    (Requirement 10.2). Raises :class:`NotebookUnavailable` with a member-facing
    message when the notebook is unreachable, times out, or answers with a
    non-success status (Requirement 10.3) — this function never returns on
    failure, so the caller always has a clean success/failure split.

    See the module's EXTERNAL HANDOFF note for the assumed ingest contract.
    """
    endpoint = _join_url(notebook_url, _NOTEBOOK_INGEST_PATH)
    headers = {
        identity.SUITE_USER_HEADER: member,   # 10.1: act as the member
        'Accept': 'application/json',
    }
    body = {
        'genome': genome_name,                # 10.4: what the figure shows
        'locus': locus_text,                  # 10.4
        'filename': 'figure.svg',
        'content_type': 'image/svg+xml',
        'svg': svg,
    }

    try:
        resp = requests.post(
            endpoint, json=body, headers=headers,
            timeout=_NOTEBOOK_TIMEOUT_SECONDS)
    except requests.RequestException as exc:
        # Connection refused, DNS failure, timeout, etc. — unreachable.
        raise NotebookUnavailable('the Lab Notebook could not be reached') from exc

    if not (200 <= resp.status_code < 300):
        raise NotebookUnavailable(
            f'the Lab Notebook rejected the figure (status {resp.status_code})')

    return _entry_reference(resp)


def _entry_reference(resp: "requests.Response") -> Dict:
    """
    Pull a created-entry reference out of the notebook's success response.

    Tolerant of the exact envelope the notebook settles on (see EXTERNAL
    HANDOFF): accepts ``id`` or ``entry_id`` for the identifier and ``url`` or
    ``href`` for a link, and always returns at least an (possibly empty) id so
    the caller has a stable shape.
    """
    try:
        data = resp.json()
    except ValueError:
        data = {}
    if not isinstance(data, dict):
        data = {}
    entry_id = data.get('id') or data.get('entry_id') or ''
    url = data.get('url') or data.get('href') or ''
    reference: Dict[str, str] = {'id': str(entry_id)}
    if url:
        reference['url'] = str(url)
    return reference


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _locus_text(locus: "render.Locus") -> str:
    """The visible locus as ``Chrom:start-end`` text for the notebook payload."""
    return f'{locus.chrom}:{locus.start}-{locus.end}'


def _join_url(base: str, path: str) -> str:
    """Join an ingest path onto NOTEBOOK_URL without doubling the slash."""
    base = (base or '').rstrip('/')
    path = '/' + path.lstrip('/')
    return f'{base}{path}'
