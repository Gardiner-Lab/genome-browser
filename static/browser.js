/*
 * The Genome Browser page: choose a genome, its tracks and where to start, then
 * browse them in igv.js; add genomes and tracks; select a stretch of sequence
 * and send it to the Sequence Editor or to NCBI BLAST. Everything is read
 * through this tool's own API, under the same mount (data-prefix on <body>,
 * normally /genome), with the signed-in session the proxy checks. Nothing is
 * fetched from the internet; only BLAST, when a member asks for it, sends the
 * selected sequence to NCBI.
 *
 *   /genome/                  the landing page: choose what to show
 *   /genome/?genome=<id>      the browser, opened on the member's saved view of it
 */
(function () {
    'use strict';

    var PREFIX = document.body.getAttribute('data-prefix') || '';
    var API = PREFIX + '/api';
    // The Sequence Editor, beside this tool in the suite (EDITOR_URL); unset, its button is hidden.
    var EDITOR = (document.body.getAttribute('data-editor') || '').replace(/\/$/, '');

    var el = function (id) { return document.getElementById(id); };
    var state = {
        genomes: [],
        genome: null,          // the chosen genome (from the list)
        detail: null,          // its chromosomes
        tracks: [],            // its track library
        selected: new Set(),   // chosen track ids
        view: null,            // the saved view for the chosen genome
        browser: null,         // the igv browser, while the view is open
        saveTimer: null,
        me: null,              // the signed-in member, from the suite's session
        pending: [],           // conversions this page started: {id, name, genomeId}
        selection: null,       // {chrom, start, end, strand, name}: 0-based, end exclusive
        lastClick: { t: 0, x: 0, y: 0 },
        quietUntil: 0          // igv's click popup is skipped until then (a double-click)
    };

    // ------------------------------------------------------------------
    // Small helpers
    // ------------------------------------------------------------------

    function getJSON(url) {
        return fetch(url, { credentials: 'same-origin', headers: { 'Accept': 'application/json' } })
            .then(function (r) {
                if (r.status === 401) { signIn(); throw new Error('Not signed in'); }
                return r.json().then(function (body) {
                    if (!r.ok) { throw new Error((body && body.error) || ('HTTP ' + r.status)); }
                    return body;
                });
            });
    }

    function postJSON(url, body) {
        return fetch(url, {
            method: 'POST', credentials: 'same-origin',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify(body)
        }).then(function (r) {
            return r.json().then(function (data) {
                if (!r.ok) { throw new Error((data && data.error) || ('HTTP ' + r.status)); }
                return data;
            });
        });
    }

    function signIn() {
        window.location.href = '/account/login?next=' + encodeURIComponent(location.pathname + location.search);
    }

    function make(tag, attrs, children) {
        var node = document.createElement(tag);
        Object.keys(attrs || {}).forEach(function (k) {
            if (k === 'text') { node.textContent = attrs[k]; }
            else if (k === 'class') { node.className = attrs[k]; }
            else if (k.slice(0, 2) === 'on') { node.addEventListener(k.slice(2), attrs[k]); }
            else { node.setAttribute(k, attrs[k]); }
        });
        (children || []).forEach(function (c) { if (c) { node.appendChild(c); } });
        return node;
    }

    // The browser's current locus as "Chr1:3,600-3,722": igv.js reports the
    // window's end as a fraction of a base when zoomed, which is not a locus
    // anyone would type, so the coordinates are rounded.
    function currentLocus(browser) {
        var loci = browser.currentLoci() || '';
        var text = String(Array.isArray(loci) ? (loci[0] || '') : loci).split(' ')[0];
        var m = /^(.+):([\d.,]+)-([\d.,]+)$/.exec(text);
        if (!m) { return text; }
        var start = Math.max(1, Math.round(Number(m[2].replace(/,/g, ''))));
        var end = Math.max(start, Math.round(Number(m[3].replace(/,/g, ''))));
        return m[1] + ':' + start.toLocaleString('en-US') + '-' + end.toLocaleString('en-US');
    }

    function formatBp(bp) {
        if (bp >= 1e9) { return (bp / 1e9).toFixed(2) + ' Gb'; }
        if (bp >= 1e6) { return (bp / 1e6).toFixed(1) + ' Mb'; }
        if (bp >= 1e3) { return (bp / 1e3).toFixed(1) + ' kb'; }
        return bp + ' bp';
    }

    function trackKind(type) {
        if (type === 'gff') { return 'Gene annotation (GFF)'; }
        if (type === 'bed') { return 'Intervals (BED)'; }
        if (type === 'bigwig') { return 'Coverage (BigWig)'; }
        return type;
    }

    // A starting window: the saved locus, else the first 100 kb of the first
    // chromosome, so gene annotations are drawn straight away.
    function defaultLocus() {
        if (state.view && state.view.locus) { return state.view.locus; }
        var chroms = (state.detail && state.detail.chromosomes) || [];
        if (!chroms.length) { return ''; }
        var first = chroms[0];
        return first.name + ':1-' + Math.min(first.length_bp, 100000).toLocaleString('en-US');
    }

    // ------------------------------------------------------------------
    // Who is signed in (the suite's session), for the header
    // ------------------------------------------------------------------

    function loadUser() {
        fetch('/account/api/session', { credentials: 'same-origin' })
            .then(function (r) { return r.ok ? r.json() : null; })
            .then(function (s) {
                if (s && s.user) {
                    state.me = s.user;
                    el('user').textContent = s.user.username;
                }
            })
            .catch(function () { /* the header just stays without a name */ });
    }

    // ------------------------------------------------------------------
    // The landing page
    // ------------------------------------------------------------------

    function showChoose(genomeId, keepSelection) {
        closeBrowser();
        el('view').hidden = true;
        el('choose').hidden = false;
        el('crumb').hidden = true;
        document.title = 'Genome Browser · Trellis Tools';
        el('genome-list').replaceChildren(make('p', { class: 'gb-loading', text: 'Loading genomes…' }));
        getJSON(API + '/genomes').then(function (genomes) {
            state.genomes = genomes;
            renderGenomes();
            var pick = genomeId || (state.genome && state.genome.id) || (genomes.length === 1 ? genomes[0].id : null);
            if (pick) { chooseGenome(pick, keepSelection); }
        }).catch(function (err) {
            el('genome-list').replaceChildren(make('p', { class: 'gb-error', text: 'Could not load the genomes: ' + err.message }));
        });
    }

    function renderGenomes() {
        var list = el('genome-list');
        if (!state.genomes.length) {
            list.replaceChildren(make('div', { class: 'gb-empty' }, [
                make('strong', { text: 'No genomes yet.' }),
                document.createTextNode(' A reference genome is added by uploading its FASTA. Ask your lab lead or an administrator to add the one you need.')
            ]));
            el('step-tracks').hidden = true;
            el('step-locus').hidden = true;
            el('open-row').hidden = true;
            return;
        }
        list.replaceChildren.apply(list, state.genomes.map(function (g) {
            var meta = g.assembly_id + ' · ' + g.chromosome_count + ' sequence' +
                (g.chromosome_count === 1 ? '' : 's') + ' · ' + formatBp(g.length_bp);
            return make('button', {
                type: 'button', class: 'gb-genome', 'data-id': g.id,
                'aria-pressed': state.genome && state.genome.id === g.id ? 'true' : 'false',
                onclick: function () { chooseGenome(g.id, false); }
            }, [
                make('span', { class: 'gb-genome-name', text: g.name }),
                make('span', { class: 'gb-genome-meta', text: meta }),
                make('span', {
                    class: 'gb-badge' + (g.suite ? ' everyone' : ''),
                    text: g.suite ? 'Everyone' : 'Your lab'
                })
            ]);
        }));
    }

    function chooseGenome(id, keepSelection) {
        var genome = state.genomes.find(function (g) { return g.id === id; });
        if (!genome) { return; }
        var changed = !state.genome || state.genome.id !== id;
        state.genome = genome;
        Array.prototype.forEach.call(document.querySelectorAll('.gb-genome'), function (b) {
            b.setAttribute('aria-pressed', b.getAttribute('data-id') === id ? 'true' : 'false');
        });
        ['step-tracks', 'step-locus', 'open-row'].forEach(function (s) { el(s).hidden = false; });
        el('track-list').replaceChildren(make('li', {}, [make('p', { class: 'gb-loading', text: 'Loading tracks…' })]));
        el('locus-error').hidden = true;

        Promise.all([
            getJSON(API + '/genomes/' + encodeURIComponent(id)),
            getJSON(API + '/genomes/' + encodeURIComponent(id) + '/tracks'),
            getJSON(API + '/view?genome_id=' + encodeURIComponent(id)).catch(function () { return null; })
        ]).then(function (results) {
            state.detail = results[0];
            state.tracks = results[1];
            state.view = results[2];
            if (changed || !keepSelection) {
                var saved = (state.view && state.view.enabled) || [];
                var available = new Set(state.tracks.map(function (t) { return t.id; }));
                var keep = saved.filter(function (t) { return available.has(t); });
                // Never opened before: start with every track if there are a few.
                if (!state.view || !state.view.updated_at) {
                    keep = state.tracks.length <= 6 ? state.tracks.map(function (t) { return t.id; }) : [];
                }
                // A track just added here starts ticked.
                state.pending.forEach(function (p) {
                    if (p.genomeId === id && p.status === 'done' && available.has(p.id)) { keep.push(p.id); }
                });
                state.selected = new Set(keep);
                el('locus').value = defaultLocus();
            }
            renderTracks();
            renderChroms();
            renderPending();
            updateOpenHint();
        }).catch(function (err) {
            el('track-list').replaceChildren(make('li', {}, [make('p', { class: 'gb-error', text: 'Could not load this genome: ' + err.message })]));
        });
    }

    function renderTracks() {
        var list = el('track-list');
        var rows = [make('li', { class: 'gb-sequence-row' }, [
            make('span', { class: 'gb-dot', text: '✓', 'aria-hidden': 'true' }),
            make('span', { text: 'Reference sequence and ruler: always shown' })
        ])];
        if (!state.tracks.length) {
            rows.push(make('li', { class: 'gb-sequence-row' }, [
                make('span', { class: 'gb-dot', text: '' }),
                make('span', { text: 'No tracks for this genome yet. Tracks are BigWig, BED or GFF files added to it.' })
            ]));
        }
        state.tracks.forEach(function (t) {
            var box = make('input', { type: 'checkbox', value: t.id });
            box.checked = state.selected.has(t.id);
            box.addEventListener('change', function () {
                if (box.checked) { state.selected.add(t.id); } else { state.selected.delete(t.id); }
                updateOpenHint();
            });
            var kind = trackKind(t.type) + (t.visibility === 'suite' ? ' · shared with everyone' : '');
            rows.push(make('li', {}, [make('label', { class: 'gb-track' }, [
                box,
                make('span', { class: 'gb-track-text' }, [
                    make('span', { class: 'gb-track-name', text: t.name }),
                    make('span', { class: 'gb-track-kind', text: kind })
                ])
            ])]));
        });
        list.replaceChildren.apply(list, rows);
        el('track-actions').hidden = state.tracks.length < 2;
    }

    function setAllTracks(on) {
        state.selected = new Set(on ? state.tracks.map(function (t) { return t.id; }) : []);
        renderTracks();
        updateOpenHint();
    }

    function renderChroms() {
        var chroms = (state.detail && state.detail.chromosomes) || [];
        var shown = chroms.slice(0, 40);
        var box = el('chrom-list');
        box.replaceChildren.apply(box, shown.map(function (c) {
            return make('button', {
                type: 'button', class: 'gb-chip', title: formatBp(c.length_bp), text: c.name,
                onclick: function () { el('locus').value = c.name; el('locus-error').hidden = true; }
            });
        }));
        if (chroms.length > shown.length) {
            box.appendChild(make('span', { class: 'gb-hint', text: '+ ' + (chroms.length - shown.length) + ' more: type a name' }));
        }
    }

    function updateOpenHint() {
        var n = state.selected.size;
        el('open-hint').textContent = n === 0
            ? 'The sequence alone; you can add tracks in the browser later.'
            : n + ' track' + (n === 1 ? '' : 's') + ' selected.';
    }

    // Checks a typed locus against the genome's chromosomes before opening, so a
    // typo is caught here rather than as an empty browser.
    function checkLocus(text) {
        var chroms = (state.detail && state.detail.chromosomes) || [];
        if (!text) { return null; }
        var name = text.split(':')[0].trim();
        if (chroms.some(function (c) { return c.name === name; })) { return null; }
        var lower = name.toLowerCase();
        var near = chroms.find(function (c) { return c.name.toLowerCase() === lower; });
        return 'There is no sequence called "' + name + '" in this genome' +
            (near ? ' (did you mean ' + near.name + '?)' : '') + '. Pick one below.';
    }

    function openFromChoice() {
        if (!state.genome) { return; }
        var locus = el('locus').value.trim();
        var problem = checkLocus(locus);
        if (problem) {
            el('locus-error').textContent = problem;
            el('locus-error').hidden = false;
            el('locus').focus();
            return;
        }
        var order = state.tracks.map(function (t) { return t.id; }).filter(function (id) { return state.selected.has(id); });
        history.pushState({ genome: state.genome.id }, '', PREFIX + '/?genome=' + encodeURIComponent(state.genome.id));
        // Save first, so the view reopens the same way and the URL is enough;
        // each track's settings (height, colour...) are kept as they were.
        var options = (state.view && state.view.options) || {};
        postJSON(API + '/view', { genome_id: state.genome.id, locus: locus || null, enabled: order, options: options })
            .catch(function () { return null; })
            .then(function () { openBrowser(state.genome.id, { locus: locus, enabled: order, options: options }); });
    }

    // ------------------------------------------------------------------
    // The browser
    // ------------------------------------------------------------------

    // A track's igv.js config, with the settings saved for it in this view
    // (height, colour, gene display, coverage scale; see trackSettings()).
    function trackConfig(t, saved) {
        var base = API + '/tracks/' + encodeURIComponent(t.id);
        var opt = saved || {};
        var config;
        if (t.type === 'bigwig') {
            config = { name: t.name, type: 'wig', format: 'bigwig', url: base + '/coverage.bw', trellisId: t.id };
            if (opt.scale === 'fixed' && isFinite(opt.min) && isFinite(opt.max)) {
                config.autoscale = false;
                config.min = Number(opt.min);
                config.max = Number(opt.max);
            } else if (opt.scale === 'group') {
                config.autoscaleGroup = SHARED_SCALE;
            }
        } else {
            var ext = t.type === 'bed' ? 'bed' : 'gff3';
            config = {
                name: t.name, type: 'annotation', format: ext,
                url: base + '/features.' + ext + '.gz', indexURL: base + '/features.' + ext + '.gz.tbi',
                indexed: true, displayMode: 'EXPANDED', visibilityWindow: 2000000, trellisId: t.id
            };
            if (DISPLAY_MODES.indexOf(opt.displayMode) >= 0) { config.displayMode = opt.displayMode; }
        }
        if (isFinite(opt.height) && opt.height >= 15) { config.height = Math.min(1000, Number(opt.height)); }
        if (typeof opt.color === 'string' && /^(#[0-9a-f]{3,8}|rgba?\([\d.,\s]+\))$/i.test(opt.color)) { config.color = opt.color; }
        return config;
    }

    var DISPLAY_MODES = ['EXPANDED', 'SQUISHED', 'COLLAPSED'];
    // The one autoscale group "Shared" puts tracks in, so it survives a reload.
    var SHARED_SCALE = 'trellis-shared-scale';

    function openBrowser(genomeId, chosen) {
        closeBrowser();
        el('choose').hidden = true;
        el('view').hidden = false;
        // The first time, the controls are shown until dismissed.
        if (state.firstControls) { showControls(true); }
        el('igv-container').replaceChildren();
        el('saved').textContent = '';

        var genomeP = state.genomes.length ? Promise.resolve(state.genomes) : getJSON(API + '/genomes');
        Promise.all([
            genomeP,
            getJSON(API + '/genomes/' + encodeURIComponent(genomeId) + '/tracks'),
            chosen ? Promise.resolve(chosen) : getJSON(API + '/view?genome_id=' + encodeURIComponent(genomeId))
        ]).then(function (results) {
            state.genomes = results[0];
            var genome = state.genomes.find(function (g) { return g.id === genomeId; });
            if (!genome) { throw new Error('This genome does not exist, or you cannot see it.'); }
            state.genome = genome;
            state.tracks = results[1];
            var view = results[2] || {};
            // Never opened: go to the landing page with it chosen instead.
            if (!chosen && !view.updated_at) {
                history.replaceState({}, '', PREFIX + '/');
                showChoose(genomeId, false);
                return;
            }
            var byId = {};
            state.tracks.forEach(function (t) { byId[t.id] = t; });
            var enabled = (view.enabled || []).filter(function (id) { return byId[id]; });
            state.selected = new Set(enabled);

            el('view-title').textContent = genome.name;
            el('view-hint').textContent = genome.assembly_id + ' · ' + enabled.length + ' track' + (enabled.length === 1 ? '' : 's');
            // The genome's name is in the toolbar under the header; in the
            // header too, it would push the suite's links onto a second row.
            document.title = genome.name + ' · Genome Browser';

            if (typeof igv === 'undefined' || !igv || typeof igv.createBrowser !== 'function') {
                el('vendor-note').hidden = false;
                return;
            }
            var ref = API + '/genomes/' + encodeURIComponent(genome.id) + '/reference.fa';
            var options = {
                // igv.js labels the genome by its id, so give it the assembly's.
                reference: { id: genome.assembly_id, name: genome.name, fastaURL: ref, indexURL: ref + '.fai' },
                tracks: enabled.map(function (id) { return trackConfig(byId[id], (view.options || {})[id]); }),
                showCenterGuide: false,
                showSampleNames: false,
                // Several tracks at once are changed in Track settings instead.
                showMultiSelectButton: false
            };
            if (view.locus) { options.locus = view.locus; }
            return igv.createBrowser(el('igv-container'), options).then(function (browser) {
                state.browser = browser;
                window.__genomeBrowser = browser;
                browser.on('locuschange', scheduleSave);
                browser.on('trackremoved', function () { scheduleSave(); renderTrackPanel(); });
                browser.on('trackorderchanged', function () { scheduleSave(); renderTrackPanel(); });
                browser.on('trackheightchange', scheduleSave);
                // A region drawn with igv's region tool becomes the selection.
                browser.on('roiadded', function (rois) {
                    var r = rois && rois[0];
                    if (!r || state.placingRoi) { return; }
                    setSelection({ chrom: r.chr, start: r.start, end: r.end, strand: 1, name: '' }, false);
                });
                // A double-click selects the feature instead of zooming; igv's
                // details popup for its first click is skipped.
                browser.on('trackclick', function () {
                    if (Date.now() < state.quietUntil) { return false; }
                    return undefined;
                });
            });
        }).catch(function (err) {
            el('igv-container').replaceChildren(make('div', { class: 'gb-empty' }, [
                make('strong', { text: 'The browser could not open. ' }),
                document.createTextNode(err && err.message ? err.message : String(err))
            ]));
        });
    }

    function closeBrowser() {
        setSelection(null);
        setSelectMode(false);
        showTrackPanel(false);
        state.panelChecked = new Set();
        if (state.saveTimer) { clearTimeout(state.saveTimer); state.saveTimer = null; }
        if (state.browser && typeof igv !== 'undefined' && igv.removeBrowser) {
            try { igv.removeBrowser(state.browser); } catch (e) { /* already gone */ }
        }
        state.browser = null;
    }

    // Keeps the member's view (where they are, which tracks, in what order) so
    // the genome reopens the same way.
    function scheduleSave() {
        if (state.saveTimer) { clearTimeout(state.saveTimer); }
        state.saveTimer = setTimeout(saveView, 800);
    }

    function saveView() {
        var browser = state.browser;
        if (!browser || !state.genome) { return; }
        var locus = currentLocus(browser);
        var tracks = browser.findTracks(function (t) { return t.config && t.config.trellisId; });
        var enabled = tracks.map(function (t) { return t.config.trellisId; });
        var options = {};
        tracks.forEach(function (t) { options[t.config.trellisId] = trackSettings(t); });
        postJSON(API + '/view', { genome_id: state.genome.id, locus: locus || null, enabled: enabled, options: options })
            .then(function () {
                el('saved').textContent = 'View saved';
                state.selected = new Set(enabled);
                el('view-hint').textContent = state.genome.assembly_id + ' · ' + enabled.length + ' track' + (enabled.length === 1 ? '' : 's');
            })
            .catch(function () { el('saved').textContent = 'Not saved'; });
    }

    // ------------------------------------------------------------------
    // Selecting sequence: double-click a feature, or draw a region
    // ------------------------------------------------------------------

    // The feature under a click, from igv's own hit-testing of the track the
    // click landed on (so in expanded mode it is the row's transcript). Of
    // several, the shortest is the most specific.
    function featureAt(event) {
        var browser = state.browser;
        if (!browser || !browser.trackViews) { return null; }
        for (var i = 0; i < browser.trackViews.length; i++) {
            var view = browser.trackViews[i];
            var track = view.track;
            if (!track || track.type !== 'annotation' || typeof track.clickedFeatures !== 'function') { continue; }
            var ports = view.viewports || [];
            for (var j = 0; j < ports.length; j++) {
                var port = ports[j];
                // By position, not by the click's target: igv lays overlays
                // (labels, regions) over its tracks, outside their elements.
                if (!port.viewportElement) { continue; }
                var box = port.viewportElement.getBoundingClientRect();
                if (event.clientX < box.left || event.clientX > box.right ||
                    event.clientY < box.top || event.clientY > box.bottom) { continue; }
                try {
                    var click = port.createClickState(event);
                    var found = click ? (track.clickedFeatures(click) || []) : [];
                    found = found.filter(function (f) { return f && f.chr && f.end > f.start; });
                    if (!found.length) { return null; }
                    found.sort(function (a, b) { return (a.end - a.start) - (b.end - b.start); });
                    return found[0];
                } catch (e) {
                    return null;
                }
            }
        }
        return null;
    }

    // Runs before igv's own click handling (capture phase): the second click of
    // a double-click on a feature selects it, and igv never sees that click, so
    // it does not zoom in.
    // The Sequence Editor is for working on a gene or a construct, not a
    // chromosome arm: longer selections are offered as a FASTA download.
    var EDITOR_MAX = 30000;
    // What the sequence API sends at once (for BLAST and copying); the FASTA
    // download has no limit.
    var INLINE_MAX = 2000000;

    function onClickCapture(event) {
        // The click that ends a drag-selection is not igv's to act on.
        if (state.swallowClick) {
            state.swallowClick = false;
            event.stopPropagation();
            event.preventDefault();
            return;
        }
        var now = Date.now();
        var last = state.lastClick;
        var isDouble = now - last.t < 450 && Math.abs(event.clientX - last.x) < 6 && Math.abs(event.clientY - last.y) < 6;
        state.lastClick = { t: isDouble ? 0 : now, x: event.clientX, y: event.clientY };
        if (!isDouble) { return; }
        var feature = featureAt(event);
        if (!feature) { return; }
        event.stopPropagation();
        event.preventDefault();
        state.quietUntil = now + 1200;
        var strand = feature.strand === '-' || feature.strand === -1 ? -1 : 1;
        setSelection({
            chrom: feature.chr, start: feature.start, end: feature.end, strand: strand,
            name: feature.name || feature.id || feature.type || ''
        }, true);
    }

    // ------------------------------------------------------------------
    // Dragging out any stretch: Shift+drag, or a plain drag in Select mode
    // ------------------------------------------------------------------

    // The tracks' data area on screen, and the window of the genome it shows.
    function trackArea() {
        var browser = state.browser;
        var frames = browser && browser.referenceFrameList;
        if (!frames || frames.length !== 1 || !browser.trackViews) { return null; }
        var top = Infinity, bottom = -Infinity, left = null, right = null;
        browser.trackViews.forEach(function (view) {
            var port = view.viewports && view.viewports[0];
            if (!port || !port.viewportElement || !view.track || view.track.type === 'ideogram') { return; }
            var box = port.viewportElement.getBoundingClientRect();
            if (!box.width || !box.height) { return; }
            top = Math.min(top, box.top);
            bottom = Math.max(bottom, box.bottom);
            left = box.left;
            right = box.right;
        });
        if (left === null) { return null; }
        return { frame: frames[0], left: left, right: right, top: top, bottom: bottom };
    }

    function bpAt(x, area) {
        var clamped = Math.min(Math.max(x, area.left), area.right);
        return area.frame.start + (clamped - area.left) * area.frame.bpPerPixel;
    }

    function setSelectMode(on) {
        state.selectMode = on;
        el('select-mode').setAttribute('aria-pressed', on ? 'true' : 'false');
        document.body.classList.toggle('gb-selecting', on);
    }

    function onMouseDownCapture(event) {
        if (event.button !== 0 || !(event.shiftKey || state.selectMode)) { return; }
        var area = trackArea();
        if (!area || event.clientX < area.left || event.clientX > area.right ||
            event.clientY < area.top || event.clientY > area.bottom) { return; }
        // Ours, not igv's: it would otherwise start moving along the genome.
        event.stopPropagation();
        event.preventDefault();
        var band = el('drag-band');
        var from = event.clientX;
        var label = make('span', {});
        band.replaceChildren(label);
        band.style.top = area.top + 'px';
        band.style.height = (area.bottom - area.top) + 'px';
        function draw(to) {
            var a = Math.max(Math.min(from, to), area.left), b = Math.min(Math.max(from, to), area.right);
            band.style.left = a + 'px';
            band.style.width = Math.max(1, b - a) + 'px';
            var bp = Math.round(Math.abs(bpAt(to, area) - bpAt(from, area)));
            label.textContent = bp.toLocaleString('en-US') + ' bp';
        }
        draw(from);
        band.hidden = false;
        function move(e) { draw(e.clientX); }
        function up(e) {
            window.removeEventListener('mousemove', move, true);
            window.removeEventListener('mouseup', up, true);
            band.hidden = true;
            // A click without a drag stays a click, so double-clicking a
            // feature still selects it in Select mode.
            if (Math.abs(e.clientX - from) < 3) { return; }
            state.swallowClick = true;
            setTimeout(function () { state.swallowClick = false; }, 300);
            var a = bpAt(Math.min(from, e.clientX), area), b = bpAt(Math.max(from, e.clientX), area);
            var start = Math.max(0, Math.floor(a)), end = Math.ceil(b);
            if (end <= start) { return; }
            setSelection({ chrom: area.frame.chr, start: start, end: end, strand: 1, name: '' }, true);
        }
        window.addEventListener('mousemove', move, true);
        window.addEventListener('mouseup', up, true);
    }

    // "Chr1:3,631-5,899" (1-based, inclusive, as the browser shows it).
    function selectTyped(event) {
        event.preventDefault();
        var text = el('select-locus').value.trim().replace(/,/g, '');
        var m = /^([^:\s]+):(\d+)\s*(?:-|\.\.)\s*(\d+)$/.exec(text);
        if (!m) { selectionError('Type a region as chromosome:start-end, for example Chr1:3,631-5,899.'); return; }
        var chrom = m[1], a = Number(m[2]), b = Number(m[3]);
        var start = Math.min(a, b) - 1, end = Math.max(a, b);
        var genome = state.browser && state.browser.genome;
        var chr = genome && genome.getChromosome ? genome.getChromosome(chrom) : null;
        if (genome && genome.getChromosome && !chr) { selectionError('This genome has no sequence called "' + chrom + '".'); return; }
        if (chr && chr.bpLength) { end = Math.min(end, chr.bpLength); }
        if (end <= start) { selectionError('That region is empty.'); return; }
        var name = chr && chr.name ? chr.name : chrom;
        setSelection({ chrom: name, start: start, end: end, strand: 1, name: '' }, true);
        // Show it, with a little room on either side.
        var pad = Math.max(20, Math.round((end - start) * 0.1));
        if (state.browser && state.browser.search) {
            state.browser.search(name + ':' + Math.max(1, start + 1 - pad) + '-' + (end + pad));
        }
    }

    // ------------------------------------------------------------------
    // The selection itself
    // ------------------------------------------------------------------

    function setSelection(sel, highlight) {
        state.selection = sel;
        el('selection-error').hidden = true;
        el('selection-note').hidden = true;
        var box = el('selection');
        var text = el('selection-text');
        if (!sel) {
            box.classList.remove('has');
            el('selection-actions').hidden = true;
            text.replaceChildren(
                make('strong', { text: 'Nothing selected.' }),
                make('span', { text: 'Shift-drag across the tracks (or turn on Select region) to highlight any stretch; double-click a feature to select exactly it.' }));
            if (state.browser && state.browser.clearROIs) { try { state.browser.clearROIs(); } catch (e) { /* gone */ } }
            return;
        }
        var length = sel.end - sel.start;
        box.classList.add('has');
        el('selection-actions').hidden = false;
        var where = sel.chrom + ':' + (sel.start + 1).toLocaleString('en-US') + '-' + sel.end.toLocaleString('en-US');
        text.replaceChildren(
            make('span', { class: 'gb-sel-label', text: 'Selected' }),
            make('strong', { text: (sel.name ? sel.name + '  ' : '') + where }),
            make('span', { text: length.toLocaleString('en-US') + ' bp' + (sel.strand === -1 ? ', minus strand (sent as it reads)' : '') }));

        // What each action can take.
        var editor = el('send-editor');
        editor.hidden = !EDITOR;
        editor.disabled = length > EDITOR_MAX;
        editor.title = editor.disabled ? 'The Sequence Editor takes up to 30 kb' : '';
        el('send-blast').disabled = el('copy-seq').disabled = length > INLINE_MAX;
        var note = '';
        if (length > INLINE_MAX) {
            note = 'This is ' + formatBp(length) + ': too long to open, BLAST or copy here. Download it as FASTA.';
        } else if (EDITOR && length > EDITOR_MAX) {
            note = 'This is ' + formatBp(length) + ', more than the Sequence Editor takes (30 kb). Download it as FASTA, or select less.';
        }
        el('selection-note').textContent = note;
        el('selection-note').hidden = !note;

        if (highlight && state.browser && state.browser.loadROI) {
            state.placingRoi = true;
            try { state.browser.clearROIs(); } catch (e) { /* none */ }
            // igv.js takes a set of regions, each with its features.
            Promise.resolve(state.browser.loadROI([{
                name: 'Selection', color: 'rgba(43, 128, 64, 0.28)',
                features: [{ chr: sel.chrom, start: sel.start, end: sel.end }]
            }]))
                .catch(function () { /* the selection stands without the highlight */ })
                .then(function () { state.placingRoi = false; });
        }
    }

    function selectionError(message) {
        el('selection-error').textContent = message;
        el('selection-error').hidden = false;
    }

    // The selection's sequence, and the annotations of the feature tracks on
    // screen inside it, on its own coordinates.
    function fetchSelection() {
        var sel = state.selection;
        var featureTracks = state.browser
            ? state.browser.findTracks(function (t) { return t.config && t.config.trellisId && t.type === 'annotation'; })
                .map(function (t) { return t.config.trellisId; })
            : [];
        var query = new URLSearchParams({
            chrom: sel.chrom, start: String(sel.start), end: String(sel.end),
            strand: String(sel.strand), tracks: featureTracks.join(',')
        });
        return getJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/sequence?' + query.toString());
    }

    function selectionTitle(sel) {
        var where = sel.chrom + ':' + (sel.start + 1) + '-' + sel.end + (sel.strand === -1 ? '(-)' : '');
        return (sel.name ? sel.name + ' ' : '') + state.genome.assembly_id + ' ' + where;
    }

    function sendToEditor() {
        var sel = state.selection;
        if (!sel || !EDITOR || sel.end - sel.start > EDITOR_MAX) { return; }
        // Opened now, while the click still counts, so no popup blocker stops it.
        var tab = window.open('', '_blank');
        if (tab) { tab.document.title = 'Opening in the Sequence Editor…'; tab.document.body.textContent = 'Opening in the Sequence Editor…'; }
        fetchSelection().then(function (region) {
            var title = selectionTitle(sel);
            return postJSON(EDITOR + '/api/sequences', {
                document: {
                    name: sel.name || title, sequence: region.sequence, circular: false,
                    description: title + ' from ' + state.genome.name + ' (Genome Browser)',
                    features: region.features
                },
                origin: 'Genome Browser: ' + title
            });
        }).then(function (saved) {
            // The editor's page, as the other tools link it (its root redirects there without the query).
            var url = EDITOR + '/editor?open=' + encodeURIComponent(saved.sequence.id);
            if (tab) { tab.location.href = url; } else { window.location.href = url; }
        }).catch(function (err) {
            if (tab) { tab.close(); }
            selectionError('Could not open it in the Sequence Editor: ' + err.message);
        });
    }

    function sendToBlast() {
        var sel = state.selection;
        if (!sel || sel.end - sel.start > INLINE_MAX) { return; }
        var name = 'blast-' + Date.now();
        var tab = window.open('', name);
        fetchSelection().then(function (region) {
            var title = selectionTitle(sel);
            var form = make('form', { method: 'POST', action: 'https://blast.ncbi.nlm.nih.gov/Blast.cgi', target: name, hidden: 'hidden' });
            var fields = { PAGE_TYPE: 'BlastSearch', PROGRAM: 'blastn', BLAST_PROGRAMS: 'megaBlast',
                           LINK_LOC: 'blasthome', JOB_TITLE: title, QUERY: '>' + title + '\n' + region.sequence };
            Object.keys(fields).forEach(function (k) {
                var input = make('input', { type: 'hidden', name: k });
                input.value = fields[k];
                form.appendChild(input);
            });
            document.body.appendChild(form);
            form.submit();
            form.remove();
        }).catch(function (err) {
            if (tab) { tab.close(); }
            selectionError('Could not send it to BLAST: ' + err.message);
        });
    }

    function copySequence() {
        var sel = state.selection;
        if (!sel || sel.end - sel.start > INLINE_MAX) { return; }
        fetchSelection().then(function (region) {
            return navigator.clipboard.writeText(region.sequence).then(function () {
                el('saved').textContent = 'Copied ' + region.sequence.length.toLocaleString('en-US') + ' bp';
            });
        }).catch(function (err) { selectionError('Could not copy it: ' + err.message); });
    }

    // Any length: the server streams the FASTA, so this is just a link to it.
    function downloadFasta() {
        var sel = state.selection;
        if (!sel) { return; }
        var query = new URLSearchParams({
            chrom: sel.chrom, start: String(sel.start), end: String(sel.end),
            strand: String(sel.strand), name: sel.name || ''
        });
        var a = make('a', {
            href: API + '/genomes/' + encodeURIComponent(state.genome.id) + '/sequence.fasta?' + query.toString(),
            download: ''
        });
        document.body.appendChild(a);
        a.click();
        a.remove();
    }

    // ------------------------------------------------------------------
    // Adding genomes and tracks
    // ------------------------------------------------------------------

    function ownerNote() {
        var me = state.me;
        if (!me) { return ''; }
        return me.lab ? 'It will belong to ' + me.lab + '.' : 'You are in no lab yet, so only a suite reference can be added.';
    }

    // Sends a form, reporting progress (0-100, then null while the server
    // checks it); resolves {status, body}.
    function send(url, data, onProgress) {
        return new Promise(function (resolve, reject) {
            var xhr = new XMLHttpRequest();
            xhr.open('POST', url);
            xhr.withCredentials = true;
            xhr.upload.onprogress = function (e) {
                if (!e.lengthComputable || !onProgress) { return; }
                onProgress(e.loaded < e.total ? Math.round(100 * e.loaded / e.total) : null);
            };
            xhr.onload = function () {
                var body = {};
                try { body = JSON.parse(xhr.responseText); } catch (e) { /* not JSON */ }
                if (xhr.status === 401) { signIn(); }
                resolve({ status: xhr.status, body: body });
            };
            xhr.onerror = function () { reject(new Error('The upload did not reach the server.')); };
            xhr.send(data);
        });
    }

    // The genome dialog's single progress bar over send().
    function upload(url, data, form) {
        var bar = form.querySelector('.gb-progress');
        var meter = bar.querySelector('progress');
        var label = bar.querySelector('span');
        bar.hidden = false;
        meter.value = 0;
        label.textContent = 'Uploading…';
        return send(url, data, function (pct) {
            if (pct === null) { meter.value = 100; label.textContent = 'Checking the file…'; }
            else { meter.value = pct; label.textContent = pct + '%'; }
        }).finally(function () { bar.hidden = true; });
    }

    function formError(form, nodes) {
        var box = form.querySelector('.gb-error');
        box.replaceChildren.apply(box, nodes.map(function (n) { return typeof n === 'string' ? make('p', { text: n }) : n; }));
        box.hidden = false;
    }

    function busy(form, on) {
        Array.prototype.forEach.call(form.querySelectorAll('button, input'), function (b) { b.disabled = on; });
    }

    function openGenomeDialog() {
        var dialog = el('genome-dialog');
        var form = el('genome-form');
        form.reset();
        form.querySelector('.gb-error').hidden = true;
        el('genome-owner').textContent = ownerNote();
        el('suite-reference-row').hidden = !(state.me && state.me.is_admin);
        if (state.me && state.me.is_admin && !state.me.lab) { form.elements.suite_reference.checked = true; }
        dialog.showModal();
        form.elements.name.focus();
    }

    function submitGenome(event) {
        event.preventDefault();
        var form = el('genome-form');
        var f = form.elements;
        if (!f.name.value.trim() || !f.assembly_id.value.trim() || !f.file.files.length) {
            formError(form, ['A name, an assembly and a FASTA file are all needed.']);
            return;
        }
        var data = new FormData();
        data.append('name', f.name.value.trim());
        data.append('assembly_id', f.assembly_id.value.trim());
        if (f.suite_reference.checked) { data.append('suite_reference', 'true'); }
        data.append('file', f.file.files[0]);
        form.querySelector('.gb-error').hidden = true;
        busy(form, true);
        upload(API + '/genomes', data, form).then(function (res) {
            busy(form, false);
            if (res.status !== 201) {
                formError(form, [res.body.error || ('The upload failed (HTTP ' + res.status + ').')]);
                return;
            }
            el('genome-dialog').close();
            state.genome = null;
            showChoose(res.body.id, false);
        }).catch(function (err) { busy(form, false); formError(form, [err.message]); });
    }

    // ------------------------------------------------------------------
    // Adding tracks: one file or many, each its own row and its own upload
    // ------------------------------------------------------------------

    // The formats the server takes (readers.TRACK_FORMAT_EXTENSIONS), by
    // extension, ignoring a trailing .gz.
    var TRACK_FORMATS = {
        gff: 'Gene annotation (GFF)', gff3: 'Gene annotation (GFF)', bed: 'Intervals (BED)',
        bw: 'Coverage (BigWig)', bigwig: 'Coverage (BigWig)',
        bam: 'BAM, converted to coverage', wig: 'wig, converted to coverage', wiggle: 'wig, converted to coverage',
        bedgraph: 'bedGraph, converted to coverage', bg: 'bedGraph, converted to coverage'
    };

    function trackFormat(filename) {
        var name = filename.toLowerCase().replace(/\.gz$/, '');
        var dot = name.lastIndexOf('.');
        return dot < 0 ? null : (TRACK_FORMATS[name.slice(dot + 1)] || null);
    }

    // "leaf_rep1.bw" -> "leaf_rep1": a starting name, changed as the member likes.
    function defaultTrackName(filename) {
        return filename.replace(/\.gz$/i, '').replace(/\.[^.]+$/, '');
    }

    function openTrackDialog() {
        if (!state.genome) { return; }
        state.uploads = [];
        state.uploading = false;
        el('track-rows').replaceChildren();
        el('track-summary').hidden = true;
        el('track-form').querySelector('.gb-error').hidden = true;
        el('track-files').value = '';
        el('track-genome-name').textContent = state.genome.name;
        el('track-owner').textContent = state.genome.suite && state.me && state.me.is_admin
            ? 'As an administrator, your tracks on a suite genome are shared with everyone.'
            : ownerNote();
        updateTrackButtons();
        el('track-dialog').showModal();
        el('track-browse').focus();
    }

    function addTrackFiles(files) {
        Array.prototype.forEach.call(files, function (file) {
            if (state.uploads.some(function (u) { return u.file.name === file.name && u.file.size === file.size; })) { return; }
            var kind = trackFormat(file.name);
            var row = { file: file, state: kind ? 'ready' : 'skipped', kind: kind };
            row.name = make('input', {
                type: 'text', maxlength: '200', value: defaultTrackName(file.name),
                'aria-label': 'Track name for ' + file.name
            });
            row.status = make('span', { class: 'gb-upload-state' });
            row.meta = make('span', { class: 'gb-upload-meta' });
            row.bar = make('progress', { max: '100', value: '0', hidden: 'hidden' });
            row.msg = make('div', { class: 'gb-upload-msg' });
            row.remove = make('button', {
                type: 'button', class: 'gb-upload-remove', title: 'Leave this file out', 'aria-label': 'Leave out ' + file.name, text: '×',
                onclick: function () {
                    state.uploads = state.uploads.filter(function (u) { return u !== row; });
                    row.li.remove();
                    updateTrackButtons();
                }
            });
            row.li = make('li', { class: 'gb-upload' }, [
                row.name, make('span', {}, [row.status, row.remove]), row.meta, row.bar, row.msg
            ]);
            state.uploads.push(row);
            el('track-rows').appendChild(row.li);
            setRow(row, row.state);
        });
        updateTrackButtons();
    }

    function setRow(row, st, message) {
        row.state = st;
        row.li.setAttribute('data-state', st);
        row.status.textContent = {
            ready: 'Ready', skipped: 'Not a track', uploading: 'Uploading…', checking: 'Checking…',
            done: 'Added', converting: 'Added, converting', refused: 'Refused', failed: 'Failed'
        }[st] || st;
        row.meta.textContent = row.file.name + ' · ' + (row.kind || 'not a format the browser takes') + ' · ' + formatBytes(row.file.size);
        row.name.disabled = st !== 'ready';
        row.remove.hidden = st !== 'ready' && st !== 'skipped';
        row.bar.hidden = st !== 'uploading' && st !== 'checking';
        if (message !== undefined) {
            row.msg.replaceChildren.apply(row.msg, [].concat(message).map(function (n) {
                return typeof n === 'string' ? make('p', { text: n }) : n;
            }));
        }
    }

    function formatBytes(n) {
        if (n >= 1e9) { return (n / 1e9).toFixed(1) + ' GB'; }
        if (n >= 1e6) { return (n / 1e6).toFixed(1) + ' MB'; }
        if (n >= 1e3) { return Math.round(n / 1e3) + ' kB'; }
        return n + ' B';
    }

    function updateTrackButtons() {
        var ready = state.uploads.filter(function (u) { return u.state === 'ready'; }).length;
        var submit = el('track-submit');
        submit.disabled = state.uploading || ready === 0;
        submit.textContent = ready > 1 ? 'Upload ' + ready + ' tracks' : 'Upload';
        var finished = state.uploads.some(function (u) { return u.state !== 'ready' && u.state !== 'skipped'; });
        el('track-cancel').textContent = finished ? 'Close' : 'Cancel';
        el('track-cancel').disabled = state.uploading;
        el('track-browse').disabled = state.uploading;
    }

    // One file to one genome; resolves when its row says how it went.
    function uploadRow(row, genomeId) {
        var name = row.name.value.trim() || defaultTrackName(row.file.name);
        var data = new FormData();
        data.append('name', name);
        data.append('file', row.file);
        row.bar.value = 0;
        setRow(row, 'uploading', []);
        return send(API + '/genomes/' + encodeURIComponent(genomeId) + '/tracks', data, function (pct) {
            if (pct === null) { row.bar.removeAttribute('value'); setRow(row, 'checking'); }
            else { row.bar.value = pct; }
        }).then(function (res) {
            var b = res.body || {};
            if (res.status === 201 || res.status === 202) {
                var notes = (b.chrom_check && b.chrom_check.messages) || [];
                var where = genomeId === state.genome.id ? [] : ['Added to ' + genomeName(genomeId) + '.'];
                recordAdded(b, genomeId, res.status === 202);
                setRow(row, res.status === 202 ? 'converting' : 'done', where.concat(notes));
                return;
            }
            var nodes = [make('p', { text: (b.error || ('The upload failed (HTTP ' + res.status + ').')) + (b.detail ? ' ' + b.detail : '') })];
            if (b.fits_genomes && b.fits_genomes.length) {
                b.fits_genomes.forEach(function (g) {
                    nodes.push(make('button', {
                        type: 'button', class: 'tt-btn primary', text: 'Add it to ' + g.name + ' instead',
                        onclick: function () { runUploads([row], g.id); }
                    }));
                });
            } else if (b.genome_chroms && b.genome_chroms.length) {
                nodes.push(make('p', { text: 'This genome has ' + b.genome_chroms.slice(0, 10).join(', ') + (b.genome_chroms.length > 10 ? ', …' : '') + '.' }));
            }
            setRow(row, res.status === 400 || res.status === 409 || res.status === 413 ? 'refused' : 'failed', nodes);
        }).catch(function (err) {
            setRow(row, 'failed', [err.message]);
        });
    }

    function genomeName(id) {
        var g = state.genomes.find(function (x) { return x.id === id; });
        return g ? g.name : 'another genome';
    }

    // One after another: a large BAM should not compete with the rest for the
    // server, and the order on screen is the order they go.
    function runUploads(rows, genomeId) {
        if (state.uploading || !rows.length) { return; }
        state.uploading = true;
        updateTrackButtons();
        el('track-summary').hidden = true;
        var chain = Promise.resolve();
        rows.forEach(function (row) {
            chain = chain.then(function () {
                row.li.scrollIntoView({ block: 'nearest' });
                return uploadRow(row, genomeId);
            });
        });
        chain.then(function () {
            state.uploading = false;
            updateTrackButtons();
            summariseUploads();
            refreshAfterUploads();
        });
    }

    function submitTracks(event) {
        event.preventDefault();
        var ready = state.uploads.filter(function (u) { return u.state === 'ready'; });
        if (!ready.length) {
            formError(el('track-form'), ['Choose or drop at least one track file.']);
            return;
        }
        el('track-form').querySelector('.gb-error').hidden = true;
        runUploads(ready, state.genome.id);
    }

    function summariseUploads() {
        var count = function (st) { return state.uploads.filter(function (u) { return u.state === st; }).length; };
        var parts = [];
        if (count('done')) { parts.push(count('done') + ' added'); }
        if (count('converting')) { parts.push(count('converting') + ' added and converting'); }
        if (count('refused')) { parts.push(count('refused') + ' refused'); }
        if (count('failed')) { parts.push(count('failed') + ' failed'); }
        if (count('skipped')) { parts.push(count('skipped') + ' not a track format'); }
        var box = el('track-summary');
        box.textContent = parts.length ? parts.join(', ') + '.' + (count('refused') ? ' Each refused file says why.' : '') : '';
        box.hidden = !parts.length;
    }

    // Each added track: ticked on its genome, and listed under the tracks
    // (a conversion until it is ready).
    function recordAdded(body, genomeId, converting) {
        var notes = (body.chrom_check && body.chrom_check.messages) || [];
        state.pending.push({ id: body.id, name: body.name, genomeId: genomeId,
                             status: converting ? 'converting' : 'done', notes: notes });
        if (converting) { pollPending(); }
        if (state.genome && genomeId === state.genome.id) { state.selected.add(body.id); }
    }

    // Once a batch is through, the chosen genome's list is read again.
    function refreshAfterUploads() {
        if (state.genome) { chooseGenome(state.genome.id, true); }
        renderPending();
    }

    function closeTrackDialog() {
        if (state.uploading) { return; }
        el('track-dialog').close();
    }

    function renderPending() {
        var list = el('pending-list');
        var mine = state.pending.filter(function (p) { return state.genome && p.genomeId === state.genome.id; });
        // Tracks added cleanly are one line, however many there were; those
        // with something to say (renamed chromosomes, a conversion) have their own.
        var plain = mine.filter(function (p) { return p.status === 'done' && !p.notes.length; });
        var rest = mine.filter(function (p) { return plain.indexOf(p) < 0; });
        var rows = [];
        if (plain.length) {
            var names = plain.map(function (p) { return p.name; });
            var shown = names.slice(0, 4).join(', ') + (names.length > 4 ? ' and ' + (names.length - 4) + ' more' : '');
            rows.push(make('li', { class: 'done' }, [
                make('span', { text: (plain.length === 1 ? shown + ' was added.' : plain.length + ' tracks were added: ' + shown + '.') }),
                make('button', {
                    type: 'button', class: 'gb-link', text: 'Dismiss',
                    onclick: function () { state.pending = state.pending.filter(function (q) { return plain.indexOf(q) < 0; }); renderPending(); }
                })
            ]));
        }
        list.replaceChildren.apply(list, rows.concat(rest.map(function (p) {
            var text = p.status === 'converting' ? p.name + ': converting to coverage on the server…'
                : p.status === 'failed' ? p.name + ': the conversion failed. ' + (p.error || '')
                : p.name + ' was added.';
            var li = make('li', { class: p.status === 'converting' ? '' : p.status }, [make('span', { text: text })]);
            p.notes.forEach(function (n) { li.appendChild(make('span', { text: n })); });
            li.appendChild(make('button', {
                type: 'button', class: 'gb-link', text: 'Dismiss',
                onclick: function () { state.pending = state.pending.filter(function (q) { return q !== p; }); renderPending(); }
            }));
            return li;
        })));
    }

    function pollPending() {
        if (state.pollTimer) { return; }
        state.pollTimer = setInterval(function () {
            var waiting = state.pending.filter(function (p) { return p.status === 'converting'; });
            if (!waiting.length) { clearInterval(state.pollTimer); state.pollTimer = null; return; }
            waiting.forEach(function (p) {
                getJSON(API + '/tracks/' + encodeURIComponent(p.id) + '/status').then(function (s) {
                    if (s.track_status === 'converting') { return; }
                    p.status = s.track_status === 'ready' ? 'done' : 'failed';
                    p.error = s.error || '';
                    if (s.chrom_check && s.chrom_check.messages) { p.notes = s.chrom_check.messages; }
                    if (state.genome && state.genome.id === p.genomeId && !el('choose').hidden) { chooseGenome(p.genomeId, true); }
                    renderPending();
                }).catch(function () { /* try again next time */ });
            });
        }, 5000);
    }

    // ------------------------------------------------------------------
    // Track settings: height, colour, display and scale, several at once
    // ------------------------------------------------------------------

    // Colours that read on white and apart from each other, the suite's first.
    var PALETTE = ['#2b8040', '#1f5fa8', '#8b5a2b', '#b42318', '#7a4a8f', '#2a7480', '#c27c0e', '#555555', '#000000'];

    // What is kept for a track in the saved view: whatever it looks like now,
    // however it was changed (this panel or a track's own ⚙ menu).
    function trackSettings(t) {
        var s = { height: Math.round(t.height) };
        if (t.color) { s.color = t.color; }
        if (t.type === 'annotation' && DISPLAY_MODES.indexOf(t.displayMode) >= 0) { s.displayMode = t.displayMode; }
        if (t.type === 'wig') {
            if (t.autoscaleGroup) { s.scale = 'group'; }
            else if (t.autoscale === false && t.dataRange) { s.scale = 'fixed'; s.min = t.dataRange.min; s.max = t.dataRange.max; }
            else { s.scale = 'auto'; }
        }
        return s;
    }

    function panelTracks() {
        return state.browser
            ? state.browser.findTracks(function (t) { return t.config && t.config.trellisId; })
            : [];
    }

    function checkedTracks() {
        var ids = state.panelChecked || new Set();
        return panelTracks().filter(function (t) { return ids.has(t.config.trellisId); });
    }

    function showTrackPanel(on) {
        el('tracks-panel').hidden = !on;
        el('tracks-toggle').setAttribute('aria-expanded', on ? 'true' : 'false');
        if (on) { renderTrackPanel(); }
    }

    function renderTrackPanel() {
        if (el('tracks-panel').hidden) { return; }
        var tracks = panelTracks();
        var present = new Set(tracks.map(function (t) { return t.config.trellisId; }));
        state.panelChecked = new Set(Array.from(state.panelChecked || []).filter(function (id) { return present.has(id); }));
        var list = el('tp-tracks');
        if (!tracks.length) {
            list.replaceChildren(make('li', {}, [make('label', { text: 'No tracks are open. Add some with Change genome or tracks.' })]));
        } else {
            list.replaceChildren.apply(list, tracks.map(function (t) {
                var id = t.config.trellisId;
                var box = make('input', { type: 'checkbox' });
                box.checked = state.panelChecked.has(id);
                box.addEventListener('change', function () {
                    if (box.checked) { state.panelChecked.add(id); } else { state.panelChecked.delete(id); }
                    updateTrackControls();
                });
                var swatch = make('span', { class: 'gb-tp-swatch-mini' });
                swatch.style.background = t.color || (t.type === 'wig' ? 'rgb(150,150,150)' : 'rgb(0,0,150)');
                return make('li', {}, [make('label', {}, [
                    box, swatch,
                    make('span', { class: 'gb-tp-name', text: t.name }),
                    make('span', { class: 'gb-tp-kind', text: t.type === 'wig' ? 'coverage' : 'features' })
                ])]);
            }));
        }
        updateTrackControls();
    }

    function updateTrackControls() {
        var ticked = checkedTracks();
        el('tp-count').textContent = ticked.length ? ticked.length + ' ticked' : 'none ticked';
        el('tp-controls').classList.toggle('off', !ticked.length);
        el('tp-display-row').hidden = !ticked.some(function (t) { return t.type === 'annotation'; });
        el('tp-scale-row').hidden = !ticked.some(function (t) { return t.type === 'wig'; });
        el('tp-error').hidden = true;
    }

    // Applies a change to every ticked track, then redraws and keeps it.
    function eachTicked(change, onlyType) {
        var ticked = checkedTracks().filter(function (t) { return !onlyType || t.type === onlyType; });
        if (!ticked.length) { return; }
        ticked.forEach(function (t) {
            try { change(t); } catch (e) { /* one odd track does not stop the rest */ }
        });
        if (state.browser.updateViews) { state.browser.updateViews(); }
        renderTrackPanel();
        scheduleSave();
    }

    function setHeight(px) {
        var h = Math.round(Number(px));
        if (!isFinite(h) || h < 15 || h > 1000) {
            el('tp-error').textContent = 'A height is between 15 and 1000 pixels.';
            el('tp-error').hidden = false;
            return;
        }
        eachTicked(function (t) { t.trackView.setTrackHeight(h, true); });
    }

    function setColour(colour) {
        eachTicked(function (t) { t.color = colour; t.trackView.repaintViews(); });
    }

    function resetColour() {
        eachTicked(function (t) {
            t.color = t._initialColor || (t.constructor && t.constructor.defaultColor) || undefined;
            t.trackView.repaintViews();
        });
    }

    function setDisplay(mode) {
        eachTicked(function (t) {
            t.displayMode = mode;
            t.config.displayMode = mode;
            t.trackView.checkContentHeight();
            t.trackView.repaintViews();
        }, 'annotation');
    }

    function setScale(kind) {
        eachTicked(function (t) {
            if (kind === 'group') { t.autoscaleGroup = SHARED_SCALE; t.autoscale = true; }
            else { t.autoscaleGroup = undefined; t.autoscale = true; }
        }, 'wig');
    }

    function setRange(event) {
        event.preventDefault();
        var min = Number(el('tp-min').value), max = Number(el('tp-max').value);
        if (el('tp-min').value === '' || el('tp-max').value === '' || !isFinite(min) || !isFinite(max) || max <= min) {
            el('tp-error').textContent = 'Give a minimum and a larger maximum.';
            el('tp-error').hidden = false;
            return;
        }
        eachTicked(function (t) { t.setDataRange({ min: min, max: max }); }, 'wig');
    }

    function removeTicked() {
        var ticked = checkedTracks();
        if (!ticked.length) { return; }
        ticked.forEach(function (t) { state.browser.removeTrack(t); });
        state.panelChecked = new Set();
        renderTrackPanel();
        scheduleSave();
    }

    // ------------------------------------------------------------------
    // Figures: the view as SVG, to download or to a Lab Notebook page
    // ------------------------------------------------------------------

    // The Lab Notebook, beside this tool in the suite (NOTEBOOK_URL); unset, saving there is hidden.
    var NOTEBOOK = (document.body.getAttribute('data-notebook') || '').replace(/\/$/, '');

    // What is on screen, as igv.js draws it for export.
    function currentFigure() {
        var svg = state.browser.toSVG();
        var locus = currentLocus(state.browser);
        var tracks = state.browser.findTracks(function (t) { return t.config && t.config.trellisId; })
            .map(function (t) { return t.name; });
        var sel = state.selection;
        var base = ['genome-browser', state.genome.assembly_id, locus.replace(/,/g, '')].join('_')
            .replace(/[^\w.-]+/g, '_');
        return { svg: svg, locus: locus, tracks: tracks, selection: sel, filename: base };
    }

    // The SVG drawn onto a canvas, twice the size for a sharp page: the
    // notebook shows pictures inline but, so that no file can run script in
    // it, offers an SVG only as a download.
    function svgToPng(svg) {
        return new Promise(function (resolve, reject) {
            var url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' }));
            var img = new Image();
            img.onload = function () {
                var scale = 2;
                var canvas = document.createElement('canvas');
                canvas.width = Math.max(1, Math.round(img.naturalWidth * scale));
                canvas.height = Math.max(1, Math.round(img.naturalHeight * scale));
                var ctx = canvas.getContext('2d');
                ctx.fillStyle = '#fff';
                ctx.fillRect(0, 0, canvas.width, canvas.height);
                ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
                URL.revokeObjectURL(url);
                canvas.toBlob(function (blob) {
                    if (blob) { resolve(blob); } else { reject(new Error('The figure could not be turned into a picture.')); }
                }, 'image/png');
            };
            img.onerror = function () { URL.revokeObjectURL(url); reject(new Error('The figure could not be drawn.')); };
            img.src = url;
        });
    }

    function defaultFigureTitle(fig) {
        var what = fig.selection && fig.selection.name ? fig.selection.name + ', ' : '';
        return 'Genome Browser: ' + what + state.genome.assembly_id + ' ' + fig.locus;
    }

    function openFigureDialog() {
        if (!state.browser) { return; }
        var form = el('figure-form');
        form.reset();
        form.querySelector('.gb-error').hidden = true;
        el('figure-done').hidden = true;
        var fig;
        try { fig = currentFigure(); } catch (e) {
            selectionError('The view could not be drawn as a figure: ' + e.message);
            return;
        }
        state.figure = fig;
        var preview = el('figure-preview');
        if (preview.src) { URL.revokeObjectURL(preview.src); }
        preview.src = URL.createObjectURL(new Blob([fig.svg], { type: 'image/svg+xml' }));
        form.elements.title.value = defaultFigureTitle(fig);
        el('figure-save').hidden = !NOTEBOOK;
        el('figure-save').disabled = false;
        el('figure-where').textContent = NOTEBOOK
            ? 'Saving makes a page under today in your Lab Notebook, with the picture, a note of what it shows, and the SVG to download for papers and posters.'
            : '';
        el('figure-dialog').showModal();
        form.elements.title.focus();
        form.elements.title.select();
    }

    function downloadFigureSvg() {
        var fig = state.figure;
        if (!fig) { return; }
        var a = make('a', { href: URL.createObjectURL(new Blob([fig.svg], { type: 'image/svg+xml' })), download: fig.filename + '.svg' });
        document.body.appendChild(a);
        a.click();
        setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }

    function escapeHtml(text) {
        return String(text).replace(/[&<>"']/g, function (c) {
            return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
        });
    }

    // Today's day page, a new page under it, the two files attached to that
    // page, then its content naming them: four calls to the notebook's own API,
    // as the member, so the page is theirs.
    function saveFigureToNotebook(event) {
        event.preventDefault();
        var fig = state.figure;
        var form = el('figure-form');
        if (!fig || !NOTEBOOK) { return; }
        var title = form.elements.title.value.trim() || defaultFigureTitle(fig);
        var note = form.elements.note.value.trim();
        var NB = NOTEBOOK + '/api';
        form.querySelector('.gb-error').hidden = true;
        el('figure-save').disabled = true;
        el('figure-save').textContent = 'Saving…';
        var page;
        svgToPng(fig.svg).then(function (png) {
            return postJSON(NB + '/days', {}).then(function (res) {
                return postJSON(NB + '/notes', { parent_id: res.note.id, title: title, type: 'page' });
            }).then(function (res) {
                page = res.note;
                var files = new FormData();
                files.append('file', new File([png], fig.filename + '.png', { type: 'image/png' }));
                files.append('file', new File([fig.svg], fig.filename + '.svg', { type: 'image/svg+xml' }));
                return fetch(NB + '/notes/' + encodeURIComponent(page.id) + '/attachments', {
                    method: 'POST', credentials: 'same-origin', body: files
                }).then(function (r) {
                    return r.json().then(function (body) {
                        if (!r.ok) { throw new Error(body.error || ('the notebook said ' + r.status)); }
                        return body.attachments;
                    });
                });
            });
        }).then(function (added) {
            var png = added[0], svgFile = added[1];
            var where = state.genome.name + ' (' + state.genome.assembly_id + '), ' + fig.locus;
            var lines = ['<p><b>From the Genome Browser:</b> ' + escapeHtml(where) + '.</p>'];
            if (fig.tracks.length) {
                lines.push('<p>Tracks: ' + escapeHtml(fig.tracks.join('; ')) + '.</p>');
            }
            if (fig.selection) {
                var s = fig.selection;
                lines.push('<p>Selected: ' + escapeHtml((s.name ? s.name + ' ' : '') + s.chrom + ':' +
                    (s.start + 1).toLocaleString('en-US') + '-' + s.end.toLocaleString('en-US') +
                    ' (' + (s.end - s.start).toLocaleString('en-US') + ' bp' + (s.strand === -1 ? ', minus strand' : '') + ')') + '.</p>');
            }
            if (note) { lines.push('<p>' + escapeHtml(note).replace(/\n/g, '<br>') + '</p>'); }
            lines.push('<p><img data-attachment="' + escapeHtml(png.id) + '" alt="' + escapeHtml(title) + '"></p>');
            lines.push('<p>The figure as SVG, for papers and posters: <a data-attachment="' + escapeHtml(svgFile.id) + '">' +
                escapeHtml(svgFile.name || (fig.filename + '.svg')) + '</a></p>');
            return fetch(NB + '/notes/' + encodeURIComponent(page.id), {
                method: 'PUT', credentials: 'same-origin',
                headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
                body: JSON.stringify({ title: title, content: [{ type: 'text', html: lines.join('') }], base_version: page.version })
            }).then(function (r) {
                return r.json().then(function (body) {
                    if (!r.ok) { throw new Error(body.error || ('the notebook said ' + r.status)); }
                    return body;
                });
            });
        }).then(function () {
            var link = make('a', { href: NOTEBOOK + '/#note=' + encodeURIComponent(page.id), target: '_blank', rel: 'noopener', text: 'Open the page →' });
            var done = el('figure-done');
            done.replaceChildren(document.createTextNode('Saved to your Lab Notebook, under today: “' + title + '”. '), link);
            done.hidden = false;
            el('figure-save').textContent = 'Saved';
        }).catch(function (err) {
            formError(form, ['Could not save it to the Lab Notebook: ' + err.message +
                (page ? ' (a page was started under today; it may be incomplete)' : '') + '. You can still download the SVG.']);
            el('figure-save').disabled = false;
            el('figure-save').textContent = 'Save to Lab Notebook';
        });
    }

    // ------------------------------------------------------------------
    // Routing
    // ------------------------------------------------------------------

    function route() {
        var genome = new URLSearchParams(location.search).get('genome');
        if (genome) { openBrowser(genome, null); } else { showChoose(null, true); }
    }

    el('choose-form').addEventListener('submit', function (e) { e.preventDefault(); openFromChoice(); });
    el('tracks-all').addEventListener('click', function () { setAllTracks(true); });
    el('tracks-none').addEventListener('click', function () { setAllTracks(false); });
    el('locus').addEventListener('input', function () { el('locus-error').hidden = true; });
    el('change').addEventListener('click', function () {
        var id = state.genome && state.genome.id;
        // Back to the choices with what the browser had open: its tracks stay
        // ticked and it starts where the member was.
        if (state.browser) {
            el('locus').value = currentLocus(state.browser);
            saveView();
        }
        history.pushState({}, '', PREFIX + '/');
        showChoose(id, true);
    });
    window.addEventListener('popstate', route);

    el('add-genome').addEventListener('click', openGenomeDialog);
    el('add-track').addEventListener('click', openTrackDialog);
    el('genome-form').addEventListener('submit', submitGenome);
    el('track-form').addEventListener('submit', submitTracks);
    el('track-browse').addEventListener('click', function () { el('track-files').click(); });
    el('track-files').addEventListener('change', function () { addTrackFiles(el('track-files').files); el('track-files').value = ''; });
    el('track-cancel').addEventListener('click', closeTrackDialog);
    // Esc closes the dialog, but not while files are going up.
    el('track-dialog').addEventListener('cancel', function (e) { if (state.uploading) { e.preventDefault(); } });
    (function () {
        var drop = el('track-drop');
        ['dragenter', 'dragover'].forEach(function (type) {
            drop.addEventListener(type, function (e) { e.preventDefault(); if (!state.uploading) { drop.classList.add('over'); } });
        });
        ['dragleave', 'drop'].forEach(function (type) {
            drop.addEventListener(type, function () { drop.classList.remove('over'); });
        });
        drop.addEventListener('drop', function (e) {
            e.preventDefault();
            if (!state.uploading && e.dataTransfer && e.dataTransfer.files.length) { addTrackFiles(e.dataTransfer.files); }
        });
        // A file dropped beside the drop zone must not open in the tab instead.
        el('track-dialog').addEventListener('dragover', function (e) { e.preventDefault(); });
        el('track-dialog').addEventListener('drop', function (e) { e.preventDefault(); });
    })();
    Array.prototype.forEach.call(document.querySelectorAll('[data-close]'), function (b) {
        b.addEventListener('click', function () { b.closest('dialog').close(); });
    });

    // Guides: the landing page's intro and the browser's Controls panel show
    // themselves until the member says they have seen them (kept in this
    // browser only; nothing breaks when storage is unavailable).
    function seen(key) { try { return localStorage.getItem(key) === '1'; } catch (e) { return false; } }
    function markSeen(key) { try { localStorage.setItem(key, '1'); } catch (e) { /* private window */ } }
    el('intro-guide').hidden = seen('gb-intro-seen');
    el('intro-dismiss').addEventListener('click', function () {
        el('intro-guide').hidden = true;
        markSeen('gb-intro-seen');
    });
    function showControls(on) {
        el('controls').hidden = !on;
        el('controls-toggle').setAttribute('aria-expanded', on ? 'true' : 'false');
        if (!on) { markSeen('gb-controls-seen'); }
    }
    el('controls-toggle').addEventListener('click', function () { showControls(el('controls').hidden); });
    el('controls-close').addEventListener('click', function () { showControls(false); });
    state.firstControls = !seen('gb-controls-seen');

    el('tracks-toggle').addEventListener('click', function () { showTrackPanel(el('tracks-panel').hidden); });
    el('tracks-panel-close').addEventListener('click', function () { showTrackPanel(false); });
    el('tp-all').addEventListener('click', function () {
        state.panelChecked = new Set(panelTracks().map(function (t) { return t.config.trellisId; }));
        renderTrackPanel();
    });
    el('tp-none').addEventListener('click', function () { state.panelChecked = new Set(); renderTrackPanel(); });
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-height]'), function (b) {
        b.addEventListener('click', function () { setHeight(b.getAttribute('data-height')); });
    });
    el('tp-height-form').addEventListener('submit', function (e) { e.preventDefault(); setHeight(el('tp-height').value); });
    el('tp-swatches').replaceChildren.apply(el('tp-swatches'), PALETTE.map(function (c) {
        var b = make('button', { type: 'button', title: c, 'aria-label': 'Colour ' + c });
        b.style.background = c;
        b.addEventListener('click', function () { setColour(c); });
        return b;
    }));
    el('tp-color').addEventListener('change', function () { setColour(el('tp-color').value); });
    el('tp-color-reset').addEventListener('click', resetColour);
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-display]'), function (b) {
        b.addEventListener('click', function () { setDisplay(b.getAttribute('data-display')); });
    });
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-scale]'), function (b) {
        b.addEventListener('click', function () { setScale(b.getAttribute('data-scale')); });
    });
    el('tp-range-form').addEventListener('submit', setRange);
    el('tp-remove').addEventListener('click', removeTicked);

    el('figure-open').addEventListener('click', openFigureDialog);
    el('figure-form').addEventListener('submit', saveFigureToNotebook);
    el('figure-svg').addEventListener('click', downloadFigureSvg);
    el('figure-close').addEventListener('click', function () { el('figure-dialog').close(); });

    el('igv-container').addEventListener('click', onClickCapture, true);
    el('igv-container').addEventListener('mousedown', onMouseDownCapture, true);
    el('select-mode').addEventListener('click', function () { setSelectMode(!state.selectMode); });
    el('select-form').addEventListener('submit', selectTyped);
    // Escape leaves Select mode, then clears the selection.
    document.addEventListener('keydown', function (e) {
        if (e.key !== 'Escape' || el('view').hidden || document.querySelector('dialog[open]')) { return; }
        if (state.selectMode) { setSelectMode(false); } else if (state.selection) { setSelection(null); }
    });
    el('send-editor').addEventListener('click', sendToEditor);
    el('send-blast').addEventListener('click', sendToBlast);
    el('copy-seq').addEventListener('click', copySequence);
    el('save-fasta').addEventListener('click', downloadFasta);
    el('clear-selection').addEventListener('click', function () { setSelection(null); });

    loadUser();
    route();
})();
