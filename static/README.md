# Static assets

This directory holds the module's client-side assets, served by the browser
page at `/genome/`.

The IGV browser (igv.js) is **vendored** here rather than loaded from a CDN, so
the suite has no outbound dependency at runtime and works on an isolated
network. The vendored `igv.min.js` (and any assets it needs) is added when the
browser page is built (task 11.2); this placeholder keeps the directory present
in version control until then.
