# Static assets

This directory holds the module's client-side assets, served by the browser
page at `/genome/` through `views.py`'s own `/static/<asset>` route, since
`app.py` turns Flask's built-in static handler off.

`igv.min.js` is **igv.js 3.8.10, vendored here unchanged** rather than loaded
from a CDN, so the suite has no outbound dependency at run time and works on an
isolated network. Its licence sits beside it as `igv.js-LICENSE.txt`, and
`../NOTICE` records it along with the libraries the bundle carries inside it.

    source  https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js
    sha256  87b991ae2a720978172e269dd5127c495e56df5af31144dff19c83124e2f5744

To move to another version: take `dist/igv.min.js` from that igv.js release,
replace the file, and update the version and checksum here, in `../NOTICE` and
in `../README.md`. No Python names a version, so nothing else changes.
