"""
The browser page and the vendored igv.js it runs on.

The page was written to boot igv.js from ``static/igv.min.js``, but the file
itself was missing for the tool's first five releases: every visitor got the
"did not load" note instead of a browser, and nothing failed to say so. These
tests hold the vendoring in place - the bundle is there, the page asks for it,
the asset route serves it, and the MIT licence it has to travel with is beside
it - so the same gap cannot reopen quietly.

    GENOME_DATA_ROOT="$(mktemp -d)" python -m pytest test_browser_page.py
"""

from pathlib import Path

import pytest

import app as app_module

HERE = Path(__file__).resolve().parent
STATIC = HERE / 'static'


@pytest.fixture
def client(tmp_path):
    """The app over a temporary data root; the page needs no database."""

    class Config(app_module.Config):
        DATA_ROOT = tmp_path / 'genome-data'

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    return application.test_client()


def test_the_viewer_is_vendored():
    """igv.min.js is in static/, and is the bundle rather than a placeholder."""
    bundle = STATIC / 'igv.min.js'
    assert bundle.is_file(), 'igv.min.js is not vendored into static/'
    head = bundle.read_text(encoding='utf-8', errors='replace')[:400]
    assert 'igv' in head and 'function' in head, 'igv.min.js does not look like the bundle'
    # The real bundle is well over a megabyte; a stub or an error page is not.
    assert bundle.stat().st_size > 500_000, 'igv.min.js is too small to be the bundle'


def test_the_viewer_licence_travels_with_it():
    """
    MIT asks that the notice goes wherever the software goes, so it ships here.
    """
    licence = STATIC / 'igv.js-LICENSE.txt'
    assert licence.is_file(), 'the igv.js licence is not beside the bundle'
    text = licence.read_text(encoding='utf-8')
    assert 'MIT License' in text
    assert 'Broad Institute' in text and 'Regents of the University of California' in text
    assert 'shall be included in' in text
    notice = (HERE / 'NOTICE').read_text(encoding='utf-8')
    assert 'igv.js' in notice, 'NOTICE does not record igv.js'


def test_the_page_boots_the_vendored_viewer(client):
    """The page asks for the vendored script under its own mount, not a CDN."""
    page = client.get('/')
    assert page.status_code == 200
    html = page.get_data(as_text=True)
    assert '/genome/static/igv.min.js' in html
    assert '/genome/static/browser.js' in html
    script = client.get('/static/browser.js').get_data(as_text=True)
    assert 'igv.createBrowser' in script
    assert 'igv-container' in html
    # Nothing is fetched from the internet: an isolated network must be enough.
    for host in ('cdn.jsdelivr.net', 'unpkg.com', 'cdnjs.cloudflare.com', 'http://', 'https://'):
        assert host not in html.replace('https://doi.org', ''), f'the page reaches out to {host}'


def test_the_asset_route_serves_the_viewer(client):
    """
    Flask's own static handler is off, so views.py serves the directory itself.
    """
    served = client.get('/static/igv.min.js')
    assert served.status_code == 200
    assert len(served.get_data()) > 500_000
    assert client.get('/static/igv.js-LICENSE.txt').status_code == 200
    # And still refuses to climb out of it. (An absolute path makes the URL
    # `/static//etc/passwd`, which Flask redirects to its normalised form
    # before routing, so follow the redirect and judge where it lands.)
    for crafted in ('../app.py', '..%2fapp.py', '/etc/passwd'):
        answer = client.get(f'/static/{crafted}', follow_redirects=True)
        assert answer.status_code != 200, f'{crafted} was served'


# ---------------------------------------------------------------------------
# Help and About, and the guides that point into them
# ---------------------------------------------------------------------------


def test_help_and_about_are_linked_and_served(client):
    page = client.get('/').get_data(as_text=True)
    assert 'href="/genome/help"' in page and 'href="/genome/about"' in page

    help_page = client.get('/help')
    assert help_page.status_code == 200
    text = help_page.get_data(as_text=True)
    assert '<span class="tt-sub">Help</span>' in text
    assert '/genome/static/trellis/modules.js' in text

    about = client.get('/about')
    assert about.status_code == 200
    assert f'Genome Browser {app_module.APP_VERSION}' in about.get_data(as_text=True)


def test_every_help_link_lands_on_a_section(client):
    import re
    help_text = client.get('/help').get_data(as_text=True)
    sections = set(re.findall(r'<section id="([^"]+)"', help_text))
    linked = set(re.findall(r'/genome/help#([\w-]+)', client.get('/').get_data(as_text=True)))
    linked |= set(re.findall(r'href="#([\w-]+)"', help_text))
    assert linked, 'the page links into the help'
    assert linked <= sections, f'links to missing help sections: {sorted(linked - sections)}'


def test_the_guides_are_on_the_page(client):
    page = client.get('/').get_data(as_text=True)
    for marker in ('id="intro-guide"', 'id="controls"', 'id="controls-toggle"',
                   'Shift</kbd> + drag', 'Double-click a feature'):
        assert marker in page, marker
