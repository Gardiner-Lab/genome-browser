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
        if (type === 'bedpe') { return 'Links (BEDPE)'; }
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
        if (changed) { state.otherTracks = null; }
        state.genome = genome;
        Array.prototype.forEach.call(document.querySelectorAll('.gb-genome'), function (b) {
            b.setAttribute('aria-pressed', b.getAttribute('data-id') === id ? 'true' : 'false');
        });
        ['step-tracks', 'step-locus', 'open-row', 'share-genome'].forEach(function (s) { el(s).hidden = false; });
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
                // A track just added here starts ticked, as does one found by
                // searching from another genome.
                state.pending.forEach(function (p) {
                    if (p.genomeId === id && p.status === 'done' && available.has(p.id)) { keep.push(p.id); }
                });
                if (state.tickAfterLoad && available.has(state.tickAfterLoad)) { keep.push(state.tickAfterLoad); }
                state.tickAfterLoad = null;
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

    // ------------------------------------------------------------------
    // Finding a track: every word typed must appear in its name or kind
    // ------------------------------------------------------------------

    function searchWords(text) {
        return String(text || '').toLowerCase().split(/[\s,;]+/).filter(Boolean);
    }

    function trackText(t) {
        return (t.name + ' ' + trackKind(t.type) + ' ' + t.type + ' ' + (t.experiment || '')
            + ' ' + (t.citation || '') + ' ' + (t.doi || '') + ' ' + (t.geo || '')
            + (t.visibility === 'suite' ? ' shared everyone' : '')).toLowerCase();
    }

    function trackMatches(t, words) {
        if (!words.length) { return true; }
        var text = trackText(t);
        return words.every(function (w) { return text.indexOf(w) >= 0; });
    }

    // The name with each searched word marked, as nodes (never as HTML).
    function highlighted(name, words) {
        if (!words.length) { return [document.createTextNode(name)]; }
        var lower = name.toLowerCase();
        var marks = [];
        words.forEach(function (w) {
            var from = 0, at;
            while ((at = lower.indexOf(w, from)) >= 0) { marks.push([at, at + w.length]); from = at + w.length; }
        });
        if (!marks.length) { return [document.createTextNode(name)]; }
        marks.sort(function (a, b) { return a[0] - b[0]; });
        var merged = [];
        marks.forEach(function (m) {
            var last = merged[merged.length - 1];
            if (last && m[0] <= last[1]) { last[1] = Math.max(last[1], m[1]); } else { merged.push(m.slice()); }
        });
        var nodes = [], pos = 0;
        merged.forEach(function (m) {
            if (m[0] > pos) { nodes.push(document.createTextNode(name.slice(pos, m[0]))); }
            nodes.push(make('mark', { text: name.slice(m[0], m[1]) }));
            pos = m[1];
        });
        if (pos < name.length) { nodes.push(document.createTextNode(name.slice(pos))); }
        return nodes;
    }

    // ------------------------------------------------------------------
    // Sorting: by experiment (grouped under headings), name, or newest first.
    // The list's order is also the order tracks open in, top to bottom.
    // ------------------------------------------------------------------

    function trackSort() {
        try { return localStorage.getItem('gb-track-sort') || 'experiment'; } catch (e) { return 'experiment'; }
    }

    function byName(a, b) { return a.name.localeCompare(b.name, undefined, { sensitivity: 'base', numeric: true }); }

    function sortedTracks(list) {
        var mode = trackSort();
        var out = list.slice();
        if (mode === 'newest') {
            out.sort(function (a, b) { return String(b.created_at || '').localeCompare(String(a.created_at || '')) || byName(a, b); });
        } else if (mode === 'experiment') {
            // Named experiments first, A to Z; tracks in none last.
            out.sort(function (a, b) {
                var x = a.experiment || '', y = b.experiment || '';
                if (!x !== !y) { return x ? -1 : 1; }
                return x.localeCompare(y, undefined, { sensitivity: 'base', numeric: true }) || byName(a, b);
            });
        } else {
            out.sort(byName);
        }
        return out;
    }

    // Where a track comes from, for its row: its first GEO accession, else that it is published.
    function trackSource(t) {
        if (t.geo) { return ' · ' + t.geo.split(' ')[0] + (t.geo.indexOf(' ') > 0 ? ' …' : ''); }
        return t.citation || t.doi ? ' · published' : '';
    }

    function geoLink(acc) {
        return 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + encodeURIComponent(acc);
    }

    function experimentNames() {
        var seen = {};
        state.tracks.forEach(function (t) { if (t.experiment) { seen[t.experiment] = true; } });
        return Object.keys(seen).sort(function (a, b) { return a.localeCompare(b, undefined, { sensitivity: 'base' }); });
    }

    function fillExperimentNames() {
        var list = el('experiment-names');
        if (!list) { return; }
        list.replaceChildren.apply(list, experimentNames().map(function (n) { return make('option', { value: n }); }));
    }

    function experimentHeading(name, tracks) {
        return make('li', { class: 'gb-exp-head' }, [
            make('span', { class: 'gb-exp-name', text: name || 'No experiment' }),
            make('span', { class: 'gb-hint', text: tracks.length + ' track' + (tracks.length === 1 ? '' : 's') }),
            make('button', {
                type: 'button', class: 'gb-link', text: 'Tick these',
                title: 'Tick every track ' + (name ? 'in “' + name + '”' : 'in no experiment'),
                onclick: function () {
                    tracks.forEach(function (t) { state.selected.add(t.id); });
                    renderTracks();
                    updateOpenHint();
                }
            })
        ]);
    }

    function renderTracks() {
        var list = el('track-list');
        var words = searchWords(state.trackQuery);
        var shown = sortedTracks(state.tracks.filter(function (t) { return trackMatches(t, words); }));
        var grouping = trackSort() === 'experiment' && shown.some(function (t) { return t.experiment; });
        fillExperimentNames();
        var rows = [];
        if (!words.length) {
            rows.push(make('li', { class: 'gb-sequence-row' }, [
                make('span', { class: 'gb-dot', text: '✓', 'aria-hidden': 'true' }),
                make('span', { text: 'Reference sequence and ruler: always shown' })
            ]));
        }
        if (!state.tracks.length) {
            rows.push(make('li', { class: 'gb-sequence-row' }, [
                make('span', { class: 'gb-dot', text: '' }),
                make('span', { text: 'No tracks for this genome yet. Tracks are BigWig, BED or GFF files added to it.' })
            ]));
        } else if (!shown.length) {
            rows.push(make('li', { class: 'gb-no-match', text: 'No track on this genome matches “' + state.trackQuery.trim() + '”.' }));
        }
        var lastGroup = null;
        shown.forEach(function (t) {
            if (grouping && (t.experiment || '') !== lastGroup) {
                lastGroup = t.experiment || '';
                rows.push(experimentHeading(t.experiment, shown.filter(function (u) { return (u.experiment || '') === lastGroup; })));
            }
            var box = make('input', { type: 'checkbox', value: t.id });
            box.checked = state.selected.has(t.id);
            box.addEventListener('change', function () {
                if (box.checked) { state.selected.add(t.id); } else { state.selected.delete(t.id); }
                updateOpenHint();
            });
            var kind = trackKind(t.type) + (t.visibility === 'suite' ? ' · shared with everyone'
                : t.visibility === 'labs' ? ' · shared with other labs'
                : t.visibility === 'chair_group' ? ' · shared with the chair group' : '')
                + (!grouping && t.experiment ? ' · ' + t.experiment : '')
                + trackSource(t);
            rows.push(make('li', { class: 'gb-track-row' }, [
                make('label', { class: 'gb-track' }, [
                    box,
                    make('span', { class: 'gb-track-text' }, [
                        make('span', { class: 'gb-track-name' }, highlighted(t.name, words)),
                        make('span', { class: 'gb-track-kind', text: kind })
                    ])
                ]),
                make('button', {
                    type: 'button', class: 'gb-link gb-track-share', text: 'Details',
                    title: 'The experiment, plants and notebook pages of “' + t.name + '”',
                    onclick: function () { openDetails(t); }
                }),
                make('button', {
                    type: 'button', class: 'gb-link gb-track-share', text: 'Share',
                    title: 'See and change which labs see “' + t.name + '”',
                    onclick: function () { openShare('track', t); }
                })
            ]));
        });
        list.replaceChildren.apply(list, rows);
        el('track-actions').hidden = state.tracks.length < 2;
        el('track-search-count').textContent = words.length
            ? shown.length + ' of ' + state.tracks.length + ' track' + (state.tracks.length === 1 ? '' : 's')
            : '';
        el('tracks-shown').hidden = !words.length || !shown.length;
        renderElsewhere();
    }

    // Ticks every track the search shows (and keeps those already ticked).
    function tickShown() {
        var words = searchWords(state.trackQuery);
        state.tracks.forEach(function (t) { if (trackMatches(t, words)) { state.selected.add(t.id); } });
        renderTracks();
        updateOpenHint();
    }

    // The member's other genomes' track lists, read once when a search first
    // needs them, so a track on the wrong genome is still found.
    function otherTracks() {
        if (state.otherTracks) { return state.otherTracks; }
        var others = state.genomes.filter(function (g) { return !state.genome || g.id !== state.genome.id; });
        state.otherTracks = Promise.all(others.map(function (g) {
            return getJSON(API + '/genomes/' + encodeURIComponent(g.id) + '/tracks')
                .then(function (tracks) { return tracks.map(function (t) { t.genome = g; return t; }); })
                .catch(function () { return []; });
        })).then(function (lists) { return [].concat.apply([], lists); });
        return state.otherTracks;
    }

    function renderElsewhere() {
        var box = el('track-elsewhere');
        var words = searchWords(state.trackQuery);
        if (!words.length || state.genomes.length < 2) { box.hidden = true; return; }
        var query = state.trackQuery;
        otherTracks().then(function (all) {
            if (state.trackQuery !== query) { return; }   // typed on since
            var found = all.filter(function (t) {
                return (!state.genome || t.genome.id !== state.genome.id) && trackMatches(t, words);
            });
            box.hidden = !found.length;
            var list = el('track-elsewhere-list');
            list.replaceChildren.apply(list, found.slice(0, 20).map(function (t) {
                return make('li', {}, [
                    make('span', { class: 'gb-track-text' }, [
                        make('span', { class: 'gb-track-name' }, highlighted(t.name, words)),
                        make('span', { class: 'gb-track-kind', text: trackKind(t.type) + ' · on ' + t.genome.name })
                    ]),
                    make('button', {
                        type: 'button', class: 'tt-btn', text: 'Go to that genome',
                        onclick: function () {
                            state.tickAfterLoad = t.id;
                            chooseGenome(t.genome.id, false);
                        }
                    })
                ]);
            }));
            if (found.length > 20) {
                list.appendChild(make('li', { class: 'gb-no-match', text: (found.length - 20) + ' more: type more of the name to narrow it.' }));
            }
        });
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
        var order = sortedTracks(state.tracks).map(function (t) { return t.id; }).filter(function (id) { return state.selected.has(id); });
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
        } else if (t.type === 'bedpe') {
            // Links between two regions (synteny and the like), drawn as arcs between them.
            config = {
                name: t.name, type: 'interact', format: 'bedpe', url: base + '/links.bedpe',
                arcType: 'nested', arcOrientation: 'UP', height: 120, trellisId: t.id
            };
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
                // Several tracks at once are changed in Bulk track settings instead.
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
                fixIgvColourPicker();
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
    // igv's own "Set track color" picker (in its shadow root): its "More Colors" panel has its hex box turned
    // off and opens past the right edge of the window, so its Ok cannot be reached. Show the box (a colour as
    // #e41a1c, or rgb(...)) and keep the panel in view; and since igv tells nobody a colour changed, save the
    // view after a colour is chosen there, so it is there next time.
    function fixIgvColourPicker() {
        var root = el('igv-container').shadowRoot;
        if (!root || root.getElementById('trellis-colour-fix')) { return; }
        var style = document.createElement('style');
        style.id = 'trellis-colour-fix';
        style.textContent =
            '.picker_wrapper.no_editor .picker_editor { position: static !important; opacity: 1 !important; ' +
            'z-index: auto !important; width: 100%; order: 1; }' +
            '.picker_wrapper.no_editor .picker_editor input { font-family: monospace; }' +
            '.picker_wrapper.popup { left: auto !important; right: 0 !important; }';
        root.appendChild(style);
        var chosen = function (event) {
            var path = event.composedPath ? event.composedPath() : [];
            var inPicker = path.some(function (n) {
                return n.classList && (n.classList.contains('igv-ui-colorpicker-container') || n.classList.contains('picker_wrapper'));
            });
            if (!inPicker) { return; }
            if (event.type === 'keyup' && event.key !== 'Enter') { return; }
            setTimeout(function () { scheduleSave(); renderTrackPanel(); }, 300);
        };
        root.addEventListener('click', chosen, true);
        root.addEventListener('keyup', chosen, true);
        // Write the colour in the box as a hex code, the way people copy colours, when the panel opens.
        root.addEventListener('focusin', function (event) {
            var input = event.target;
            if (!input || !input.closest || !input.closest('.picker_editor')) { return; }
            var m = /rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(input.value);
            if (m) {
                input.value = '#' + [m[1], m[2], m[3]].map(function (v) { return ('0' + (+v).toString(16)).slice(-2); }).join('');
                input.select();
            }
        });
    }

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
        if (!isDouble) {
            clearIfOutside(event);
            return;
        }
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

    // A single click on the tracks outside the selected stretch clears it, as
    // clicking away from a selection does anywhere else. Not a click that
    // ended a drag along the genome, and not one inside the selection.
    function clearIfOutside(event) {
        var sel = state.selection;
        var down = state.downAt;
        if (!sel || !down) { return; }
        if (Math.abs(event.clientX - down.x) > 3 || Math.abs(event.clientY - down.y) > 3) { return; }
        var area = trackArea();
        if (!area || event.clientX < area.left || event.clientX > area.right ||
            event.clientY < area.top || event.clientY > area.bottom) { return; }
        var bp = bpAt(event.clientX, area);
        if (area.frame.chr === sel.chrom && bp >= sel.start && bp <= sel.end) { return; }
        setSelection(null);
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
                make('span', { text: 'Shift-drag across the tracks (or turn on Select region) to highlight any stretch; double-click a feature to select exactly it. Click outside a selection to clear it.' }));
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
        gff: 'Gene annotation (GFF)', gff3: 'Gene annotation (GFF)', bed: 'Intervals (BED)', bedpe: 'Links (BEDPE)',
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
        data.append('experiment', (el('track-experiment').value || '').trim());
        data.append('citation', (el('track-citation').value || '').trim());
        data.append('doi', (el('track-doi').value || '').trim());
        data.append('geo', (el('track-geo').value || '').trim());
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
        state.otherTracks = null;      // a search reads the other genomes' lists again
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

    function panelTrackMatches(t, words) {
        if (!words.length) { return true; }
        var text = (t.name + ' ' + (t.type === 'wig' ? 'coverage bigwig' : 'features gene annotation')).toLowerCase();
        return words.every(function (w) { return text.indexOf(w) >= 0; });
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
        var all = panelTracks();
        var present = new Set(all.map(function (t) { return t.config.trellisId; }));
        var words = searchWords(el('tp-search').value);
        var tracks = all.filter(function (t) { return panelTrackMatches(t, words); });
        state.panelChecked = new Set(Array.from(state.panelChecked || []).filter(function (id) { return present.has(id); }));
        var list = el('tp-tracks');
        if (!all.length) {
            list.replaceChildren(make('li', {}, [make('label', { text: 'No tracks are open. Add some with Change genome or tracks.' })]));
        } else if (!tracks.length) {
            list.replaceChildren(make('li', {}, [make('label', { text: 'No open track matches that.' })]));
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
                    make('span', { class: 'gb-tp-name' }, highlighted(t.name, words)),
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
        var shown = state.browser.findTracks(function (t) { return t.config && t.config.trellisId; });
        var tracks = shown.map(function (t) { return t.name; });
        var byId = {};
        (state.tracks || []).forEach(function (t) { byId[t.id] = t; });
        var recipe = { kind: 'view', tracks: shown.map(function (t) {
            var known = byId[t.config.trellisId] || {};
            var range = t.dataRange || {};
            return { id: t.config.trellisId, name: t.name, type: known.type, color: t.color || t.config.color,
                     height: t.height, displayMode: t.displayMode, autoscale: t.autoscale, min: range.min, max: range.max };
        }) };
        var sel = state.selection;
        var base = ['genome-browser', state.genome.assembly_id, locus.replace(/,/g, '')].join('_')
            .replace(/[^\w.-]+/g, '_');
        return { svg: svg, locus: locus, tracks: tracks, selection: sel, filename: base, recipe: recipe };
    }

    // The SVG drawn onto a canvas, twice the size for a sharp page: the
    // notebook shows pictures inline but, so that no file can run script in
    // it, offers an SVG only as a download.
    function svgToPng(svg) {
        return new Promise(function (resolve, reject) {
            var url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml' }));
            var img = new Image();
            img.onload = function () {
                var scale = img.naturalWidth > 1500 ? 1 : 2;
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
        if (fig.kind === 'circos') { return 'Circos plot: ' + (fig.title || state.genome.assembly_id + ', ' + fig.locus); }
        if (fig.kind === 'analysis') { return 'Genome Browser: ' + (fig.title || 'analysis') + ' (' + state.genome.assembly_id + ')'; }
        var what = fig.selection && fig.selection.name ? fig.selection.name + ', ' : '';
        return 'Genome Browser: ' + what + state.genome.assembly_id + ' ' + fig.locus;
    }

    function openFigureDialog(made) {
        if (!state.browser) { return; }
        var form = el('figure-form');
        form.reset();
        form.querySelector('.gb-error').hidden = true;
        el('figure-done').hidden = true;
        var fig = made && made.svg ? made : null;
        if (!fig) {
            try { fig = currentFigure(); } catch (e) {
                selectionError('The view could not be drawn as a figure: ' + e.message);
                return;
            }
        }
        el('figure-dialog-title').textContent = fig.kind === 'circos' ? 'Circos plot'
            : fig.kind === 'analysis' ? (fig.title || 'Analysis') : 'Save this view as a figure';
        el('figure-what').textContent = fig.kind === 'circos'
            ? 'The ' + fig.locus + ' with the tracks shown, as rings (outermost first, as in the legend), binned by ' +
              fig.binSize.toLocaleString('en-US') + ' bp. Circos’s own files can be downloaded to change it and run Circos again.'
            : 'Exactly what is on screen: the tracks, the position and any selection.';
        el('figure-circos-files').hidden = !fig.bundle;
        if (fig.kind === 'analysis') { el('figure-what').textContent = fig.about || ''; }
        var dl = el('figure-downloads');
        dl.replaceChildren.apply(dl, (fig.downloads || []).map(function (d) {
            return make('a', { class: 'gb-link', href: d.url, download: d.filename, text: d.label });
        }));
        if (fig.notes && fig.notes.length) {
            el('figure-what').textContent += ' ' + fig.notes.join(' ');
        }
        state.figure = fig;
        var preview = el('figure-preview');
        if (preview.src) { URL.revokeObjectURL(preview.src); }
        preview.src = URL.createObjectURL(new Blob([fig.svg], { type: 'image/svg+xml' }));
        form.elements.title.value = defaultFigureTitle(fig);
        el('figure-save').hidden = !NOTEBOOK;
        el('figure-save').disabled = false;
        el('figure-where').textContent = NOTEBOOK
            ? 'Saving makes a page under today in your Lab Notebook, with the picture, a note of what it shows, and the SVG ' +
              'to download for papers and posters; and under it a page with the code and files to make the figure again.'
            : '';
        el('figure-dialog').showModal();
        form.elements.title.focus();
        form.elements.title.select();
    }

    function downloadCircosFiles() {
        var fig = state.figure;
        if (!fig || !fig.bundle) { return; }
        var bytes = Uint8Array.from(atob(fig.bundle), function (c) { return c.charCodeAt(0); });
        var a = make('a', { href: URL.createObjectURL(new Blob([bytes], { type: 'application/zip' })), download: fig.filename + '-circos.zip' });
        document.body.appendChild(a);
        a.click();
        setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }

    // ------------------------------------------------------------------
    // Circos plots: the genome round a circle, each track shown as a ring
    // ------------------------------------------------------------------

    // The tracks on screen, top to bottom, with the colours they are drawn in.
    function shownTracks() {
        if (!state.browser) { return []; }
        return state.browser.findTracks(function (t) { return t.config && t.config.trellisId; }).map(function (t) {
            var type = t.type === 'interact' ? 'bedpe' : t.type === 'wig' ? 'bigwig' : 'features';
            return { id: t.config.trellisId, name: t.name, color: t.color || t.config.color || '', type: type };
        });
    }

    // What was chosen last time, per genome, so drawing again starts from it.
    function circosChoices() {
        state.circos = state.circos || {};
        var key = state.genome.id;
        state.circos[key] = state.circos[key] || { rings: {} };
        return state.circos[key];
    }

    function chromosomesOfGenome() {
        var g = state.browser && state.browser.genome;
        var names = (g && (g.chromosomeNames || g.wgChromosomeNames)) || [];
        return names.filter(function (n) { return n && n.toLowerCase() !== 'all'; });
    }

    function openCircosDialog() {
        if (!state.browser || !state.genome) { return; }
        var form = el('circos-form');
        var saved = circosChoices();
        form.querySelector('.gb-error').hidden = true;
        var tracks = shownTracks();

        var list = el('circos-tracks');
        list.replaceChildren.apply(list, tracks.length ? tracks.map(function (t, i) {
            var swatch = make('span', { class: 'gb-circos-swatch', style: 'background:' + (t.color || '#888') });
            var name = make('span', { class: 'gb-circos-name', text: (i + 1) + '. ' + t.name });
            if (t.type === 'bedpe') {
                return make('li', {}, [swatch, name, make('span', { class: 'gb-hint', text: 'links: arcs across the circle' })]);
            }
            var mine = saved.rings[t.id] || {};
            var style = make('select', { 'data-ring': t.id, 'data-what': 'style', 'aria-label': 'How to draw ' + t.name },
                [['histogram', 'Histogram'], ['line', 'Line'], ['heatmap', 'Heatmap']].map(function (o) {
                    return make('option', { value: o[0], text: o[1] });
                }));
            style.value = mine.style || 'histogram';
            var height = make('select', { 'data-ring': t.id, 'data-what': 'height', 'aria-label': 'How tall the ring for ' + t.name + ' is' },
                [['1', 'Thin'], ['2', 'Normal'], ['3', 'Tall']].map(function (o) { return make('option', { value: o[0], text: o[1] }); }));
            height.value = String(mine.height || 2);
            return make('li', {}, [swatch, name, style, height]);
        }) : [make('li', { class: 'gb-hint', text: 'No tracks are shown: show some first, then come back.' })]);

        // Where: the whole genome, chosen chromosomes, or a region (the one in view to start with).
        var chroms = chromosomesOfGenome();
        var frame = state.browser.referenceFrameList && state.browser.referenceFrameList[0];
        var inView = frame && frame.chr && frame.chr.toLowerCase() !== 'all' ? frame.chr : '';
        var picked = new Set(saved.chroms || (inView ? [inView] : chroms.slice(0, 1)));
        var boxes = el('circos-chroms');
        boxes.replaceChildren.apply(boxes, chroms.map(function (c) {
            var box = make('input', { type: 'checkbox', value: c });
            box.checked = picked.has(c);
            box.addEventListener('change', function () { form.elements.where.value = 'chroms'; });
            return make('label', {}, [box, document.createTextNode(' ' + c)]);
        }));
        form.elements.where.value = saved.where || 'genome';
        // A selection fills in the region (its exact stretch); else the last region, else the one in view.
        var picked0 = state.selection;
        form.elements.region.value = picked0
            ? picked0.chrom + ':' + (picked0.start + 1).toLocaleString('en-US') + '-' + picked0.end.toLocaleString('en-US')
            : saved.region || (inView ? currentLocus(state.browser) : '');
        ['title', 'resolution', 'step', 'labels'].forEach(function (k) {
            if (saved[k] !== undefined) { form.elements[k].value = saved[k]; }
        });
        ['skip_small', 'shared_scale', 'log_scale'].forEach(function (k) {
            if (saved[k] !== undefined) { form.elements[k].checked = saved[k]; }
        });
        var sel = state.selection;
        el('circos-highlight-row').hidden = !sel;
        // Off unless ticked: a selection fills in the region, but shading the rings is a choice.
        form.elements.highlight.checked = false;
        if (sel) {
            el('circos-highlight-text').textContent = (sel.name ? sel.name + ', ' : '') + sel.chrom + ':' +
                (sel.start + 1).toLocaleString('en-US') + '-' + sel.end.toLocaleString('en-US');
        }
        el('circos-draw').disabled = !tracks.length;
        el('circos-draw').textContent = 'Draw';
        el('circos-dialog').showModal();
    }

    function drawCircos(event) {
        event.preventDefault();
        var form = el('circos-form');
        var tracks = shownTracks();
        if (!tracks.length) { return; }
        form.querySelector('.gb-error').hidden = true;
        var saved = circosChoices();
        Array.prototype.forEach.call(form.querySelectorAll('[data-ring]'), function (s) {
            var id = s.getAttribute('data-ring');
            var r = saved.rings[id] = saved.rings[id] || {};
            r[s.getAttribute('data-what')] = s.getAttribute('data-what') === 'height' ? Number(s.value) : s.value;
        });
        var where = form.elements.where.value;
        var chroms = Array.prototype.filter.call(el('circos-chroms').querySelectorAll('input'), function (b) { return b.checked; })
            .map(function (b) { return b.value; });
        if (where === 'chroms' && !chroms.length) { formError(form, ['Tick at least one chromosome.']); return; }
        if (where === 'region' && !form.elements.region.value.trim()) {
            formError(form, ['Write the region, as Chr1:1,000,000-5,000,000.']);
            return;
        }
        var step = Number(form.elements.step.value) || null;
        Object.assign(saved, {
            where: where, chroms: chroms, region: form.elements.region.value.trim(), title: form.elements.title.value,
            resolution: form.elements.resolution.value, step: form.elements.step.value, labels: form.elements.labels.value,
            skip_small: form.elements.skip_small.checked, shared_scale: form.elements.shared_scale.checked,
            log_scale: form.elements.log_scale.checked, highlight: form.elements.highlight.checked
        });
        var sel = state.selection;
        el('circos-draw').disabled = true;
        el('circos-draw').textContent = 'Drawing… (up to a minute)';
        postJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/circos', {
            tracks: tracks.map(function (t) {
                var r = saved.rings[t.id] || {};
                return { id: t.id, color: t.color, style: r.style || 'histogram', height: r.height || 2 };
            }),
            chroms: where === 'chroms' ? chroms : null,
            region: where === 'region' ? form.elements.region.value.trim() : null,
            resolution: form.elements.resolution.value, step: step,
            skip_small: form.elements.skip_small.checked,
            shared_scale: form.elements.shared_scale.checked, log_scale: form.elements.log_scale.checked,
            highlight: sel && form.elements.highlight.checked ? { chrom: sel.chrom, start: sel.start, end: sel.end } : null,
            labels: form.elements.labels.value, title: form.elements.title.value.trim()
        }).then(function (out) {
            el('circos-dialog').close();
            var base = ['circos', state.genome.assembly_id, out.where].join('_').replace(/[^\w.-]+/g, '_').slice(0, 80);
            openFigureDialog({
                kind: 'circos', svg: out.svg, bundle: out.bundle, binSize: out.bin_size, notes: out.notes || [],
                locus: out.where, filename: base, selection: null, title: out.title,
                tracks: out.rings.map(function (r) { return r.name + ' (' + r.what + ')'; })
            });
        }).catch(function (err) {
            formError(form, ['Circos could not draw it: ' + err.message]);
        }).then(function () {
            el('circos-draw').disabled = false;
            el('circos-draw').textContent = 'Draw';
        });
    }

    // Typing in the region box chooses "a region".
    el('circos-form').elements.region.addEventListener('input', function () {
        el('circos-form').elements.where.value = 'region';
    });

    // ------------------------------------------------------------------
    // Analyse: profiles, heatmaps, differences, replicate checks, tables, bar charts, dot plots
    // ------------------------------------------------------------------

    var ANALYSES = {
        profile: { label: 'Profile around features, and heatmap',
            about: 'The average signal of the coverage tracks chosen across a set of features (genes, transposons, ' +
                   'the regions of a BED track), scaled from start to end, with the stretch before and after. Tick ' +
                   'Heatmap for every feature as a row as well. On a whole genome this takes a few minutes.',
            needs: ['tracks', 'features', 'profile'] },
        compare: { label: 'Difference between two tracks (a new track)',
            about: 'A new coverage track of one track against another (A minus B, or their log2 ratio), kept for your ' +
                   'lab: show it like any other track, in the browser and in Circos plots.',
            needs: ['compare'] },
        correlation: { label: 'Replicate check: correlation and PCA',
            about: 'How alike the coverage tracks chosen are along the whole genome, in bins: a correlation heatmap ' +
                   'and a principal component plot. Replicates should sit together, and mutants apart.',
            needs: ['tracks', 'correlation'] },
        genes: { label: 'Signal per feature (a table)',
            about: 'The mean signal of each coverage track over every feature of a type (for example each gene), ' +
                   'as a CSV table to sort and filter. With two tracks, their difference too.',
            needs: ['tracks', 'features'] },
        bars: { label: 'Tracks over the selection (bar chart)',
            about: 'The mean signal of each coverage track over the selection you made in the browser, as bars.',
            needs: ['tracks', 'selection'] },
        dotplot: { label: 'Synteny dot plot (links track)',
            about: 'A links (BEDPE) track drawn with the genome along both sides: each link a line from one region ' +
                   'to the other, so collinear blocks show as diagonals.',
            needs: ['links'] }
    };

    function genomeTracks(type) {
        return (state.tracks || []).filter(function (t) {
            return (Array.isArray(type) ? type.indexOf(t.type) >= 0 : t.type === type) && (!t.track_status || t.track_status === 'ready');
        });
    }

    function shownIds() {
        return new Set(shownTracks().map(function (t) { return t.id; }));
    }

    function openAnalyseDialog() {
        if (!state.genome) { return; }
        var form = el('analyse-form');
        form.querySelector('.gb-error').hidden = true;
        el('analyse-status').textContent = '';
        el('analyse-result').replaceChildren();
        var kind = el('analyse-kind');
        kind.replaceChildren.apply(kind, Object.keys(ANALYSES).map(function (k) {
            return make('option', { value: k, text: ANALYSES[k].label });
        }));
        kind.value = state.analyseKind || 'profile';

        // Coverage tracks, the shown ones ticked.
        var shown = shownIds();
        var coverage = genomeTracks('bigwig');
        var list = el('analyse-tracks');
        list.replaceChildren.apply(list, coverage.length ? coverage.map(function (t) {
            var box = make('input', { type: 'checkbox', value: t.id });
            box.checked = shown.has(t.id);
            return make('label', {}, [box, document.createTextNode(' ' + t.name)]);
        }) : [make('p', { class: 'gb-hint', text: 'This genome has no coverage tracks yet.' })]);

        // Feature tracks for profiles and tables, and links tracks for the dot plot.
        var features = genomeTracks(['gff', 'bed']);
        var featureSelect = el('analyse-feature-track');
        featureSelect.replaceChildren.apply(featureSelect, features.map(function (t) {
            return make('option', { value: t.id, text: t.name + (t.type === 'gff' ? ' (gene annotation)' : ' (BED)') });
        }));
        var links = genomeTracks('bedpe');
        var linkSelect = el('analyse-links-track');
        linkSelect.replaceChildren.apply(linkSelect, links.length ? links.map(function (t) {
            return make('option', { value: t.id, text: t.name });
        }) : [make('option', { value: '', text: 'No links (BEDPE) tracks on this genome' })]);
        ['analyse-a', 'analyse-b'].forEach(function (id, i) {
            var s = el(id);
            s.replaceChildren.apply(s, coverage.map(function (t) { return make('option', { value: t.id, text: t.name }); }));
            if (coverage[i]) { s.value = coverage[i].id; }
        });
        var sel = state.selection;
        el('analyse-selection-text').textContent = sel
            ? sel.chrom + ':' + (sel.start + 1).toLocaleString('en-US') + '-' + sel.end.toLocaleString('en-US') +
              ' (' + (sel.end - sel.start).toLocaleString('en-US') + ' bp)'
            : 'Nothing is selected: select a region in the browser first (Shift + drag, or type it).';
        showAnalyseFields();
        loadFeatureTypes();
        el('analyse-run').disabled = false;
        el('analyse-dialog').showModal();
    }

    function showAnalyseFields() {
        var kind = el('analyse-kind').value;
        state.analyseKind = kind;
        var needs = ANALYSES[kind].needs;
        el('analyse-about').textContent = ANALYSES[kind].about;
        Array.prototype.forEach.call(document.querySelectorAll('#analyse-form [data-for]'), function (part) {
            part.hidden = part.getAttribute('data-for').split(' ').every(function (n) { return needs.indexOf(n) < 0 && n !== kind; });
        });
        el('analyse-type-hint').textContent = kind === 'genes'
            ? 'Tick one type: each feature of it is a row of the table.'
            : 'Tick up to four: each is drawn as its own line (genes and transposons, for example).';
        el('analyse-run').textContent = kind === 'compare' ? 'Make the track' : 'Draw';
    }

    // The feature types of the chosen gene annotation or BED track, most common first.
    function loadFeatureTypes() {
        var id = el('analyse-feature-track').value;
        var box = el('analyse-types');
        if (!id) { box.replaceChildren(make('p', { class: 'gb-hint', text: 'This genome has no gene annotation or BED track.' })); return; }
        box.replaceChildren(make('p', { class: 'gb-hint', text: 'Loading the feature types…' }));
        getJSON(API + '/tracks/' + encodeURIComponent(id) + '/feature-types').then(function (res) {
            var preferred = ['gene', 'transposable_element', 'transposable_element_gene'];
            box.replaceChildren.apply(box, res.types.slice(0, 24).map(function (t, i) {
                var input = make('input', { type: 'checkbox', value: t.type });
                input.checked = preferred.indexOf(t.type) === 0 || (res.types.length === 1 && i === 0) ||
                    (!res.types.some(function (x) { return x.type === 'gene'; }) && i === 0);
                return make('label', {}, [input, document.createTextNode(' ' + t.type + ' (' + t.count.toLocaleString('en-US') + ')')]);
            }));
        }).catch(function (err) {
            box.replaceChildren(make('p', { class: 'gb-hint', text: 'Could not read its feature types: ' + err.message }));
        });
    }

    function checkedValues(container) {
        return Array.prototype.filter.call(container.querySelectorAll('input[type=checkbox]'), function (b) { return b.checked; })
            .map(function (b) { return b.value; });
    }

    function runAnalysis(event) {
        event.preventDefault();
        var form = el('analyse-form');
        form.querySelector('.gb-error').hidden = true;
        el('analyse-result').replaceChildren();
        var kind = el('analyse-kind').value;
        var tracks = checkedValues(el('analyse-tracks'));
        var types = checkedValues(el('analyse-types'));
        var names = form.elements.names.value.split(/[\s,;]+/).filter(Boolean);
        var featureTrack = el('analyse-feature-track').value;

        if (kind === 'bars') { return drawSelectionBars(tracks); }
        if (kind === 'dotplot') { return drawDotPlot(el('analyse-links-track').value); }

        var body = { kind: kind, tracks: tracks, title: form.elements.title.value.trim() };
        if (kind === 'profile' || kind === 'genes') {
            if (!featureTrack || !types.length) { formError(form, ['Choose the features: a track and at least one type.']); return; }
            if (kind === 'genes' && types.length > 1) { formError(form, ['For a table, tick one feature type.']); return; }
            if (types.length > 4) { formError(form, ['Tick at most four feature types.']); return; }
            body.regions = types.map(function (t) { return { track: featureTrack, type: t, names: names }; });
        }
        if (kind === 'profile') {
            Object.assign(body, { mode: form.elements.mode.value, before: Number(form.elements.flank.value),
                after: Number(form.elements.flank.value), body: Number(form.elements.body.value),
                bin: Number(form.elements.bin.value), heatmap: form.elements.heatmap.checked,
                kmeans: Number(form.elements.kmeans.value), colormap: form.elements.colormap.value });
        }
        if (kind === 'correlation') {
            if (tracks.length < 2) { formError(form, ['Tick at least two coverage tracks.']); return; }
            Object.assign(body, { method: form.elements.method.value, bin: Number(form.elements.corr_bin.value) });
        }
        if (kind === 'compare') {
            Object.assign(body, { a: el('analyse-a').value, b: el('analyse-b').value, operation: form.elements.operation.value,
                bin: Number(form.elements.compare_bin.value), pseudocount: Number(form.elements.pseudocount.value) || 1,
                name: form.elements.track_name.value.trim(), tracks: [] });
            if (body.a === body.b) { formError(form, ['Choose two different tracks.']); return; }
        }
        if (kind !== 'compare' && !tracks.length) { formError(form, ['Tick at least one coverage track.']); return; }

        el('analyse-run').disabled = true;
        el('analyse-status').textContent = 'Starting…';
        postJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/analyses', body).then(function (job) {
            watchAnalysis(job.id, Date.now());
        }).catch(function (err) {
            formError(form, ['Could not start it: ' + err.message]);
            el('analyse-run').disabled = false;
            el('analyse-status').textContent = '';
        });
    }

    // Asks after a job every few seconds until it is done, saying how long it has taken.
    function watchAnalysis(id, since) {
        var form = el('analyse-form');
        getJSON(API + '/analyses/' + encodeURIComponent(id)).then(function (job) {
            var secs = Math.round((Date.now() - since) / 1000);
            var took = secs < 90 ? secs + ' s' : Math.round(secs / 60) + ' min';
            if (job.state === 'queued' || job.state === 'running') {
                el('analyse-status').textContent = job.state === 'queued'
                    ? 'Waiting for ' + job.ahead + ' other job' + (job.ahead === 1 ? '' : 's') + ' to finish… (' + took + ')'
                    : 'Working… (' + took + '). You can close this and keep browsing; it carries on.';
                state.analyseTimer = setTimeout(function () { watchAnalysis(id, since); }, 3000);
                return;
            }
            el('analyse-run').disabled = false;
            if (job.state === 'failed') {
                el('analyse-status').textContent = '';
                formError(form, [job.error || 'It failed.']);
                return;
            }
            el('analyse-status').textContent = 'Done in ' + took + '.';
            showAnalysisResult(job);
        }).catch(function (err) {
            el('analyse-run').disabled = false;
            formError(form, ['Lost track of it: ' + err.message]);
        });
    }

    function analysisFile(job, name) {
        return API + '/analyses/' + encodeURIComponent(job.id) + '/files/' + encodeURIComponent(name);
    }

    function fetchText(url) {
        return fetch(url, { credentials: 'same-origin' }).then(function (r) {
            if (!r.ok) { throw new Error('HTTP ' + r.status); }
            return r.text();
        });
    }

    function fetchDataUrl(url) {
        return fetch(url, { credentials: 'same-origin' }).then(function (r) {
            if (!r.ok) { throw new Error('HTTP ' + r.status); }
            return r.blob();
        }).then(function (blob) {
            return new Promise(function (resolve) {
                var reader = new FileReader();
                reader.onload = function () { resolve(reader.result); };
                reader.readAsDataURL(blob);
            });
        });
    }

    function imageSize(src) {
        return new Promise(function (resolve, reject) {
            var img = new Image();
            img.onload = function () { resolve({ w: img.naturalWidth || 800, h: img.naturalHeight || 600 }); };
            img.onerror = function () { reject(new Error('A picture could not be read.')); };
            img.src = src;
        });
    }

    // Pictures side by side as one SVG, so they save and print as one figure.
    function composeFigure(sources, title) {
        return Promise.all(sources.map(function (src) {
            return (/^data:/.test(src) ? Promise.resolve(src) : fetchDataUrl(src)).then(function (data) {
                return imageSize(data).then(function (size) { return { data: data, size: size }; });
            });
        })).then(function (parts) {
            var height = Math.max.apply(null, parts.map(function (p) { return p.size.h; }));
            var top = title ? 50 : 0;
            var x = 0;
            var images = parts.map(function (p) {
                var out = '<image x="' + x + '" y="' + top + '" width="' + p.size.w + '" height="' + p.size.h + '" href="' + p.data + '"/>';
                x += p.size.w + 20;
                return out;
            });
            var width = x - 20;
            return '<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + (height + top) + '" viewBox="0 0 ' +
                width + ' ' + (height + top) + '"><rect width="100%" height="100%" fill="#fff"/>' +
                (title ? '<text x="10" y="32" font-family="Helvetica, Arial, sans-serif" font-size="24" font-weight="bold">' +
                    escapeHtml(title) + '</text>' : '') + images.join('') + '</svg>';
        });
    }

    function showAnalysisResult(job) {
        var r = job.result || {};
        var title = el('analyse-form').elements.title.value.trim();
        var base = ['genome-browser', state.genome.assembly_id, job.kind].join('_').replace(/[^\w.-]+/g, '_');
        var downloads = (r.files || []).filter(function (f) { return !/\.(svg|png)$/.test(f); }).map(function (f) {
            return { label: { 'matrix.gz': 'Matrix (deepTools)', 'correlation.tab': 'Correlation table', 'bins.tab': 'Values per bin',
                'signal-per-feature.csv': 'Table (CSV)' }[f] || f, url: analysisFile(job, f), filename: f };
        });
        if (job.kind === 'compare') {
            el('analyse-result').replaceChildren(make('div', { class: 'gb-notice', text: 'Made the new track, for your lab: adding it…' }));
            reloadTracksAndShow(r.track_id);
            return;
        }
        if (job.kind === 'genes') {
            var head = ['Name', 'Position'].concat(r.columns || []);
            var table = make('table', { class: 'gb-analyse-table' }, [
                make('thead', {}, [make('tr', {}, head.map(function (h) { return make('th', { text: h }); }))]),
                make('tbody', {}, (r.preview || []).map(function (row) {
                    return make('tr', {}, [make('td', { text: row.name }), make('td', { text: row.chrom + ':' + row.start.toLocaleString('en-US') })]
                        .concat(row.values.map(function (v) { return make('td', { class: 'num', text: v === null ? '' : String(v) }); })));
                }))]);
            el('analyse-result').replaceChildren(
                make('p', {}, [document.createTextNode((r.count || 0).toLocaleString('en-US') + ' features measured. ' +
                    ((r.columns || []).length === 2 ? 'The 50 with the largest difference between the two tracks:' : 'The 50 highest:') + ' '),
                    make('a', { href: downloads[0].url, download: 'signal-per-feature.csv', text: 'Download the whole table (CSV)' })]),
                make('div', { class: 'gb-analyse-table-wrap' }, [table]));
            return;
        }
        var pictures = (r.files || []).filter(function (f) { return /\.(svg|png)$/.test(f); });
        var heading = title || (job.kind === 'profile' ? 'Profile' + ((r.regions || []).length ? ': ' + r.regions.join(', ') : '')
            : 'Replicate check');
        composeFigure(pictures.map(function (f) { return analysisFile(job, f); }), heading).then(function (svg) {
            el('analyse-dialog').close();
            openFigureDialog({
                kind: 'analysis', svg: svg, locus: state.genome.assembly_id + ', ' + (job.kind === 'profile' ? 'profile' : 'replicate check'),
                filename: base, selection: null, title: heading, downloads: downloads,
                recipe: { kind: 'deeptools', jobId: job.id },
                about: job.kind === 'profile'
                    ? 'Drawn by deepTools from the tracks chosen, over ' + (r.regions || []).join(', ') + '.'
                    : 'Drawn by deepTools: correlation of the tracks along the genome, and their principal components.',
                tracks: checkedValues(el('analyse-tracks')).map(function (id) {
                    var t = state.tracks.find(function (x) { return x.id === id; });
                    return t ? t.name : id;
                })
            });
        }).catch(function (err) {
            formError(el('analyse-form'), ['The pictures could not be shown: ' + err.message]);
        });
    }

    // A new track from a comparison: list it, and show it in the view.
    function reloadTracksAndShow(trackId) {
        getJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/tracks').then(function (tracks) {
            state.tracks = tracks;
            if (trackId) { state.selected.add(trackId); }
            renderTracks();
            var t = tracks.find(function (x) { return x.id === trackId; });
            if (t && state.browser) {
                state.browser.loadTrack(trackConfig(t));
                scheduleSave();
            }
            var note = el('analyse-result').querySelector('.gb-notice');
            if (note && t) { note.textContent = 'Made the track “' + t.name + '”, for your lab, and showed it. It is in the track list too.'; }
        });
    }

    // Bars of each coverage track's mean over the selection, drawn here.
    function drawSelectionBars(tracks) {
        var form = el('analyse-form');
        var sel = state.selection;
        if (!sel) { formError(form, ['Select a region in the browser first.']); return; }
        if (!tracks.length) { formError(form, ['Tick at least one coverage track.']); return; }
        el('analyse-run').disabled = true;
        postJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/region-means',
            { tracks: tracks, chrom: sel.chrom, start: sel.start, end: sel.end }).then(function (res) {
            var colours = {};
            shownTracks().forEach(function (t) { colours[t.id] = t.color; });
            var rows = res.tracks.filter(function (t) { return t.mean !== null && t.mean !== undefined; });
            if (!rows.length) { throw new Error('none of the tracks has data there'); }
            var where = sel.chrom + ':' + (sel.start + 1).toLocaleString('en-US') + '-' + sel.end.toLocaleString('en-US');
            var title = form.elements.title.value.trim() || 'Mean signal over ' + (sel.name ? sel.name + ', ' : '') + where;
            var low = Math.min(0, Math.min.apply(null, rows.map(function (t) { return t.mean; })));
            var high = Math.max(0, Math.max.apply(null, rows.map(function (t) { return t.mean; })));
            var span = (high - low) || 1;
            var barW = 56, gap = 22, left = 190, top = 70, plotH = 320, labelH = 150;
            var width = left + rows.length * (barW + gap) + 20;
            var y = function (v) { return top + (high - v) / span * plotH; };
            var parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="' + Math.max(width, 420) + '" height="' + (top + plotH + labelH) +
                '" font-family="Helvetica, Arial, sans-serif"><rect width="100%" height="100%" fill="#fff"/>',
                '<text x="10" y="30" font-size="18" font-weight="bold">' + escapeHtml(title) + '</text>',
                '<line x1="' + left + '" y1="' + y(0) + '" x2="' + (width - 10) + '" y2="' + y(0) + '" stroke="#444"/>',
                '<line x1="' + left + '" y1="' + top + '" x2="' + left + '" y2="' + (top + plotH) + '" stroke="#444"/>'];
            [high, low, 0].filter(function (v, i, all) { return all.indexOf(v) === i; }).forEach(function (v) {
                parts.push('<text x="' + (left - 8) + '" y="' + (y(v) + 4) + '" font-size="12" text-anchor="end">' +
                    Number(v.toPrecision(3)) + '</text>');
            });
            rows.forEach(function (t, i) {
                var x = left + gap / 2 + i * (barW + gap);
                var y0 = y(Math.max(0, t.mean)), h = Math.abs(y(t.mean) - y(0));
                parts.push('<rect x="' + x + '" y="' + y0 + '" width="' + barW + '" height="' + Math.max(1, h) + '" fill="' +
                    escapeHtml(colours[t.id] || '#4a7cb0') + '"><title>' + escapeHtml(t.name + ': ' + t.mean) + '</title></rect>');
                parts.push('<text transform="translate(' + (x + barW / 2 + 4) + ',' + (top + plotH + 14) + ') rotate(-40)" ' +
                    'text-anchor="end" font-size="12">' +
                    escapeHtml(t.name.length > 28 ? t.name.slice(0, 27) + '…' : t.name) + '</text>');
            });
            parts.push('</svg>');
            el('analyse-dialog').close();
            openFigureDialog({ kind: 'analysis', svg: parts.join(''), locus: where, filename: 'selection-means_' + where.replace(/[^\w.-]+/g, '_'),
                selection: sel, title: title, about: 'The mean of each coverage track over the selection.',
                recipe: { kind: 'bars', chrom: sel.chrom, start: sel.start, end: sel.end,
                          tracks: rows.map(function (t) { return { id: t.id, color: colours[t.id] }; }) },
                tracks: rows.map(function (t) { return t.name + ': ' + Number(t.mean.toPrecision(4)); }) });
        }).catch(function (err) {
            formError(form, ['Could not draw it: ' + err.message]);
        }).then(function () { el('analyse-run').disabled = false; });
    }

    // The genome along both sides; each link a line from one region to the other (both ways round).
    // The genome's chromosomes and lengths: the landing page loads them; a page opened straight into the
    // browser (a link, a reload) has not, so they are fetched then.
    function genomeChromosomes() {
        if (state.detail && state.detail.id === state.genome.id && (state.detail.chromosomes || []).length) {
            return Promise.resolve(state.detail.chromosomes);
        }
        return getJSON(API + '/genomes/' + encodeURIComponent(state.genome.id)).then(function (detail) {
            state.detail = detail;
            return detail.chromosomes || [];
        });
    }

    function drawDotPlot(trackId) {
        var form = el('analyse-form');
        if (!trackId) { formError(form, ['Add a links (BEDPE) track to this genome first.']); return; }
        genomeChromosomes().then(function (all) { drawDotPlotOn(trackId, all); }).catch(function (err) {
            formError(form, ['The genome’s chromosomes could not be read: ' + err.message]);
        });
    }

    function drawDotPlotOn(trackId, all) {
        var form = el('analyse-form');
        var chroms = (all || []).slice();
        var total = chroms.reduce(function (s, c) { return s + c.length_bp; }, 0);
        if (!total) { formError(form, ['The genome’s chromosomes are not known.']); return; }
        var big = chroms.filter(function (c) { return c.length_bp >= total * 0.01; });
        if (big.length) { chroms = big; total = chroms.reduce(function (s, c) { return s + c.length_bp; }, 0); }
        var offset = {}, at = 0;
        chroms.forEach(function (c) { offset[c.name] = at; at += c.length_bp; });
        el('analyse-run').disabled = true;
        fetchText(API + '/tracks/' + encodeURIComponent(trackId) + '/links.bedpe').then(function (text) {
            var size = 760, pad = 70, scale = size / total;
            var t = state.tracks.find(function (x) { return x.id === trackId; }) || { name: 'links' };
            var colour = (shownTracks().find(function (x) { return x.id === trackId; }) || {}).color || '#7b3294';
            var lines = [], n = 0, skipped = 0;
            text.split('\n').some(function (line) {
                if (!line || line[0] === '#' || /^(track|browser|chrom1)/i.test(line)) { return false; }
                var f = line.split(line.indexOf('\t') >= 0 ? '\t' : /\s+/);
                if (f.length < 6 || !(f[0] in offset) || !(f[3] in offset)) { skipped++; return false; }
                var x1 = pad + (offset[f[0]] + Number(f[1])) * scale, x2 = pad + (offset[f[0]] + Number(f[2])) * scale;
                var y1 = pad + (offset[f[3]] + Number(f[4])) * scale, y2 = pad + (offset[f[3]] + Number(f[5])) * scale;
                if (![x1, x2, y1, y2].every(isFinite)) { skipped++; return false; }
                lines.push('<line x1="' + x1.toFixed(1) + '" y1="' + y1.toFixed(1) + '" x2="' + x2.toFixed(1) + '" y2="' + y2.toFixed(1) + '"/>');
                lines.push('<line x1="' + y1.toFixed(1) + '" y1="' + x1.toFixed(1) + '" x2="' + y2.toFixed(1) + '" y2="' + x2.toFixed(1) + '"/>');
                n++;
                return n >= 25000;
            });
            if (!n) { throw new Error('none of its links are on the chromosomes drawn'); }
            var grid = [], labels = [];
            chroms.forEach(function (c) {
                var p = pad + offset[c.name] * scale, mid = pad + (offset[c.name] + c.length_bp / 2) * scale;
                grid.push('<line x1="' + p + '" y1="' + pad + '" x2="' + p + '" y2="' + (pad + size) + '"/>' +
                    '<line x1="' + pad + '" y1="' + p + '" x2="' + (pad + size) + '" y2="' + p + '"/>');
                labels.push('<text x="' + mid + '" y="' + (pad - 10) + '" text-anchor="middle">' + escapeHtml(c.name) + '</text>' +
                    '<text x="' + (pad - 10) + '" y="' + (mid + 4) + '" text-anchor="end">' + escapeHtml(c.name) + '</text>');
            });
            var title = form.elements.title.value.trim() || 'Dot plot: ' + t.name;
            var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="' + (size + pad + 20) + '" height="' + (size + pad + 20) +
                '" font-family="Helvetica, Arial, sans-serif" font-size="13"><rect width="100%" height="100%" fill="#fff"/>' +
                '<text x="10" y="22" font-size="17" font-weight="bold">' + escapeHtml(title) + '</text>' +
                '<rect x="' + pad + '" y="' + pad + '" width="' + size + '" height="' + size + '" fill="none" stroke="#444"/>' +
                '<g stroke="#ccc" stroke-width="1">' + grid.join('') + '</g>' + labels.join('') +
                '<g stroke="' + escapeHtml(colour) + '" stroke-width="2" stroke-linecap="round">' + lines.join('') + '</g></svg>';
            el('analyse-dialog').close();
            openFigureDialog({ kind: 'analysis', svg: svg, locus: 'whole genome', filename: 'dotplot_' + t.name.replace(/[^\w.-]+/g, '_'),
                selection: null, title: title,
                recipe: { kind: 'dotplot', trackId: trackId, colour: colour, chroms: chroms },
                about: n.toLocaleString('en-US') + ' links drawn, each both ways round' +
                    (n >= 25000 ? ' (the first 25,000)' : '') + (skipped ? '; ' + skipped + ' not on the chromosomes drawn' : '') + '.',
                tracks: [t.name] });
        }).catch(function (err) {
            formError(form, ['Could not draw it: ' + err.message]);
        }).then(function () { el('analyse-run').disabled = false; });
    }

    el('analyse-kind').addEventListener('change', showAnalyseFields);
    el('analyse-feature-track').addEventListener('change', loadFeatureTypes);

    // ------------------------------------------------------------------
    // The code to make a figure again: a page under the figure's page in the notebook
    // ------------------------------------------------------------------

    function suiteUrl(path) { return location.origin + path; }

    function trackFile(id) {
        var t = (state.tracks || []).find(function (x) { return x.id === id; }) || { name: id, export_name: id };
        return { id: id, name: t.name, file: t.export_name || t.name, type: t.type,
                 url: suiteUrl(API + '/tracks/' + encodeURIComponent(id) + '/export') };
    }

    function hexColour(colour, fallback) {
        var c = String(colour || '').trim();
        if (/^#[0-9a-f]{6}$/i.test(c)) { return c; }
        if (/^#[0-9a-f]{3}$/i.test(c)) { return '#' + c.slice(1).split('').map(function (x) { return x + x; }).join(''); }
        var m = /rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i.exec(c);
        if (m) { return '#' + [m[1], m[2], m[3]].map(function (v) { return ('0' + Math.min(255, +v).toString(16)).slice(-2); }).join(''); }
        return fallback || '#4a7cb0';
    }

    function pyString(text) { return JSON.stringify(String(text)); }

    function codeBlock(text) { return '<pre><code>' + escapeHtml(text) + '</code></pre>'; }

    function inputsHtml(files) {
        if (!files.length) { return ''; }
        return '<h3>1. Get the inputs</h3><p>Download each into one folder, keeping the file name (signed in to the ' +
            'suite; or Export, in the track’s settings in the Genome Browser):</p><ul>' + files.map(function (f) {
                return '<li><a href="' + escapeHtml(f.url) + '">' + escapeHtml(f.file) + '</a>' +
                    (f.name && f.name !== f.file ? ': ' + escapeHtml(f.name) : '') + '</li>';
            }).join('') + '</ul>';
    }

    function textFile(name, text, type) {
        return new File([text], name, { type: type || 'text/plain' });
    }

    // Each kind of figure's recipe: {html, files: [File]}.
    function figureRecipe(fig) {
        var made = 'Made ' + new Date().toISOString().slice(0, 16).replace('T', ' ') + ' UTC' +
            (state.me && state.me.username ? ' by ' + state.me.username : '') + ', from ' +
            state.genome.name + ' (' + state.genome.assembly_id + ') in the Trellis Tools Genome Browser.';
        var r = fig.recipe || {};
        if (fig.kind === 'circos') { return Promise.resolve(circosRecipe(fig, made)); }
        if (r.kind === 'deeptools') { return deeptoolsRecipe(fig, r, made); }
        if (r.kind === 'bars') { return Promise.resolve(barsRecipe(fig, r, made)); }
        if (r.kind === 'dotplot') { return Promise.resolve(dotplotRecipe(fig, r, made)); }
        if (r.kind === 'view') { return Promise.resolve(viewRecipe(fig, r, made)); }
        return Promise.resolve(null);
    }

    function circosRecipe(fig, made) {
        var bytes = Uint8Array.from(atob(fig.bundle), function (c) { return c.charCodeAt(0); });
        var zip = new File([bytes], fig.filename + '-circos.zip', { type: 'application/zip' });
        var html = '<p>' + escapeHtml(made) + ' Drawn by Circos 0.69-9 (circos.ca). Everything it needs, the configuration ' +
            'and the data (each track already binned), is in the attached zip: no download from the suite.</p>' +
            '<h3>1. Install Circos 0.69</h3>' +
            codeBlock('# Debian or Ubuntu\nsudo apt install circos liblist-moreutils-perl\n\n# or with conda\nconda install -c bioconda circos=0.69.9') +
            '<h3>2. Unpack and run</h3>' +
            codeBlock('unzip ' + zip.name + '\ncd circos\n# A Debian install keeps Circos’s own files in /etc/circos; the configuration looks in etc/\n' +
                '[ -d etc ] || ln -s /etc/circos etc\ncircos -conf circos.conf      # writes circos.png and circos.svg') +
            '<p>The legend in the notebook’s picture is added by the Genome Browser (the rings are listed in README.txt). ' +
            'To change the plot, edit <code>circos.conf</code> (colours, ring heights, a heatmap instead of a histogram) or the ' +
            'data files, and run Circos again.</p>' +
            '<p>Cite Circos: Krzywinski M et al. (2009) Circos: an information aesthetic for comparative genomics. Genome Research 19: 1639–1645.</p>';
        return { html: html, files: [zip] };
    }

    function deeptoolsRecipe(fig, r, made) {
        return getJSON(API + '/analyses/' + encodeURIComponent(r.jobId) + '/recipe').then(function (rec) {
            return fetch(API + '/analyses/' + encodeURIComponent(r.jobId) + '/files/recipe.zip', { credentials: 'same-origin' })
                .then(function (res) {
                    if (!res.ok) { throw new Error('the recipe could not be fetched (HTTP ' + res.status + ')'); }
                    return res.blob();
                }).then(function (blob) {
                    var zip = new File([blob], fig.filename + '-recipe.zip', { type: 'application/zip' });
                    var byFile = {};
                    (state.tracks || []).forEach(function (t) { byFile[t.export_name] = t; });
                    var files = (rec.inputs || []).map(function (name) {
                        var t = byFile[name];
                        return t ? trackFile(t.id) : { file: name, name: '', url: '' };
                    });
                    var html = '<p>' + escapeHtml(made) + ' Drawn by deepTools ' + escapeHtml(rec.deeptools || '') + '.</p>' +
                        inputsHtml(files) +
                        (rec.regions && rec.regions.length ? '<p>The features used are in the attached zip, as BED files (' +
                            escapeHtml(rec.regions.join(', ')) + ').</p>' : '') +
                        '<h3>2. Install deepTools</h3>' +
                        codeBlock('pip install deeptools==' + (rec.deeptools || '3.5.1') + '\n# or: conda install -c bioconda deeptools=' + (rec.deeptools || '3.5.1')) +
                        '<h3>3. Run</h3><p>Unzip the attached recipe into the folder with the inputs and run <code>sh run.sh</code>. ' +
                        'The commands, as they were run:</p>' + codeBlock((rec.commands || []).join('\n\n')) +
                        '<p>Cite deepTools: Ramírez F et al. (2016) deepTools2: a next generation web server for deep-sequencing ' +
                        'data analysis. Nucleic Acids Research 44: W160–W165.</p>';
                    return { html: html, files: [zip] };
                });
        });
    }

    function barsRecipe(fig, r, made) {
        var files = r.tracks.map(function (t) { return Object.assign(trackFile(t.id), { colour: hexColour(t.color) }); });
        var script = [
            '"""' + fig.title + '.',
            made,
            'Needs: pip install pyBigWig matplotlib. Run in the folder with the BigWig files."""',
            'import matplotlib',
            'matplotlib.use("Agg")',
            'import matplotlib.pyplot as plt',
            'import pyBigWig',
            '',
            'TITLE = ' + pyString(fig.title),
            'CHROM, START, END = ' + pyString(r.chrom) + ', ' + r.start + ', ' + r.end + '   # 0-based, end not included',
            'TRACKS = [  # (file, label, colour)',
            files.map(function (f) { return '    (' + pyString(f.file) + ', ' + pyString(f.name) + ', ' + pyString(f.colour) + '),'; }).join('\n'),
            ']',
            '',
            'means = []',
            'for path, label, colour in TRACKS:',
            '    bw = pyBigWig.open(path)',
            '    means.append(bw.stats(CHROM, START, min(END, bw.chroms(CHROM)), type="mean")[0])',
            '    bw.close()',
            '',
            'shown = [(t, m) for t, m in zip(TRACKS, means) if m is not None]',
            'fig, ax = plt.subplots(figsize=(max(4, 0.9 * len(shown) + 1.5), 5))',
            'ax.bar(range(len(shown)), [m for _, m in shown], color=[t[2] for t, _ in shown])',
            'ax.set_xticks(range(len(shown)))',
            'ax.set_xticklabels([t[1] for t, _ in shown], rotation=40, ha="right")',
            'ax.axhline(0, color="#444", linewidth=0.8)',
            'ax.set_ylabel("mean signal")',
            'ax.set_title(TITLE)',
            'fig.tight_layout()',
            'fig.savefig("selection-means.svg")',
            'fig.savefig("selection-means.png", dpi=200)',
            'for (path, label, colour), m in zip(TRACKS, means):',
            '    print(f"{label}\\t{m}")',
            ''].join('\n');
        var name = fig.filename + '.py';
        var html = '<p>' + escapeHtml(made) + '</p>' + inputsHtml(files) +
            '<h3>2. Install</h3>' + codeBlock('pip install pyBigWig matplotlib') +
            '<h3>3. Run</h3><p>With the attached <code>' + escapeHtml(name) + '</code> in the same folder:</p>' +
            codeBlock('python3 ' + name + '      # writes selection-means.svg and .png, and prints each mean') +
            '<p>The script:</p>' + codeBlock(script);
        return { html: html, files: [textFile(name, script, 'text/x-python')] };
    }

    function dotplotRecipe(fig, r, made) {
        var file = trackFile(r.trackId);
        var script = [
            '"""' + fig.title + '.',
            made,
            'Needs: pip install matplotlib. Run in the folder with the links (BEDPE) file."""',
            'import matplotlib',
            'matplotlib.use("Agg")',
            'import matplotlib.pyplot as plt',
            'from matplotlib.collections import LineCollection',
            '',
            'TITLE = ' + pyString(fig.title),
            'LINKS = ' + pyString(file.file),
            'COLOUR = ' + pyString(hexColour(r.colour, '#7b3294')),
            'LIMIT = 25000   # links drawn at most, the first in the file',
            'CHROMS = [  # (name, length): the chromosomes drawn, in order',
            r.chroms.map(function (c) { return '    (' + pyString(c.name) + ', ' + c.length_bp + '),'; }).join('\n'),
            ']',
            '',
            'offset, at = {}, 0',
            'for name, length in CHROMS:',
            '    offset[name] = at',
            '    at += length',
            'segments = []',
            'for line in open(LINKS):',
            '    if not line.strip() or line.startswith(("#", "track", "browser")) or line.lower().startswith("chrom1"):',
            '        continue',
            '    f = line.rstrip("\\n").split("\\t") if "\\t" in line else line.split()',
            '    if len(f) < 6 or f[0] not in offset or f[3] not in offset:',
            '        continue',
            '    x1, x2 = offset[f[0]] + int(f[1]), offset[f[0]] + int(f[2])',
            '    y1, y2 = offset[f[3]] + int(f[4]), offset[f[3]] + int(f[5])',
            '    segments += [[(x1, y1), (x2, y2)], [(y1, x1), (y2, x2)]]   # both ways round',
            '    if len(segments) >= 2 * LIMIT:',
            '        break',
            '',
            'fig, ax = plt.subplots(figsize=(8, 8))',
            'ax.add_collection(LineCollection(segments, colors=COLOUR, linewidths=1.5))',
            'for name, length in CHROMS:',
            '    ax.axvline(offset[name], color="#ccc", linewidth=0.8)',
            '    ax.axhline(offset[name], color="#ccc", linewidth=0.8)',
            'ticks = [offset[n] + l / 2 for n, l in CHROMS]',
            'ax.set_xticks(ticks)',
            'ax.set_xticklabels([n for n, _ in CHROMS])',
            'ax.set_yticks(ticks)',
            'ax.set_yticklabels([n for n, _ in CHROMS])',
            'ax.set_xlim(0, at)',
            'ax.set_ylim(at, 0)',
            'ax.set_aspect("equal")',
            'ax.set_title(TITLE)',
            'fig.tight_layout()',
            'fig.savefig("dotplot.svg")',
            'fig.savefig("dotplot.png", dpi=200)',
            'print(len(segments) // 2, "links drawn")',
            ''].join('\n');
        var name = fig.filename + '.py';
        var html = '<p>' + escapeHtml(made) + '</p>' + inputsHtml([file]) +
            '<h3>2. Install</h3>' + codeBlock('pip install matplotlib') +
            '<h3>3. Run</h3><p>With the attached <code>' + escapeHtml(name) + '</code> in the same folder:</p>' +
            codeBlock('python3 ' + name + '      # writes dotplot.svg and dotplot.png') +
            '<p>The script:</p>' + codeBlock(script);
        return { html: html, files: [textFile(name, script, 'text/x-python')] };
    }

    // A browser view: the same igv.js, the same tracks and place, as a page to open locally.
    function viewRecipe(fig, r, made) {
        var genomeFile = state.genome.assembly_id.replace(/[^\w.-]+/g, '_') + '.fa';
        var files = [{ file: genomeFile, name: 'the reference sequence (FASTA)',
                       url: suiteUrl(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/reference.fa') }];
        var tracks = [], index = [];
        r.tracks.forEach(function (t) {
            var f = Object.assign(trackFile(t.id), {});
            files.push(f);
            var c = { name: t.name, color: t.color || undefined, height: t.height || undefined };
            if (t.type === 'bigwig') {
                Object.assign(c, { type: 'wig', format: 'bigwig', url: f.file });
                if (t.min !== undefined && t.max !== undefined && t.autoscale === false) { Object.assign(c, { autoscale: false, min: t.min, max: t.max }); }
            } else if (t.type === 'bedpe') {
                Object.assign(c, { type: 'interact', format: 'bedpe', url: f.file, arcType: 'nested' });
            } else {
                var ext = t.type === 'bed' ? 'bed' : 'gff3';
                Object.assign(c, { type: 'annotation', format: ext, url: f.file + '.gz', indexURL: f.file + '.gz.tbi',
                                   displayMode: t.displayMode || 'EXPANDED' });
                index.push(t.type === 'bed'
                    ? 'sort -k1,1 -k2,2n ' + JSON.stringify(f.file) + ' | bgzip > ' + JSON.stringify(f.file + '.gz') + ' && tabix -p bed ' + JSON.stringify(f.file + '.gz')
                    : '(grep "^#" ' + JSON.stringify(f.file) + '; grep -v "^#" ' + JSON.stringify(f.file) + ' | sort -k1,1 -k4,4n) | bgzip > ' +
                      JSON.stringify(f.file + '.gz') + ' && tabix -p gff ' + JSON.stringify(f.file + '.gz'));
            }
            tracks.push(c);
        });
        var session = { reference: { id: state.genome.assembly_id, name: state.genome.name, fastaURL: genomeFile, indexURL: genomeFile + '.fai' },
                        locus: fig.locus, tracks: tracks };
        if (fig.selection) {
            session.roi = [{ name: 'Selection', color: 'rgba(255,190,0,0.25)', features: [{ chr: fig.selection.chrom, start: fig.selection.start, end: fig.selection.end }] }];
        }
        var sessionText = JSON.stringify(session, null, 2);
        var page = ['<!doctype html>', '<meta charset="utf-8">', '<title>' + escapeHtml(fig.title || 'Genome Browser view') + '</title>',
            '<script src="https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js"></scr' + 'ipt>',
            '<div id="igv"></div>', '<p><button id="svg">Save as SVG</button></p>',
            '<script>',
            'fetch("session.json").then(function (r) { return r.json(); }).then(function (s) {',
            '  return igv.createBrowser(document.getElementById("igv"), s);',
            '}).then(function (browser) {',
            '  document.getElementById("svg").onclick = function () {',
            '    var a = document.createElement("a");',
            '    a.href = URL.createObjectURL(new Blob([browser.toSVG()], { type: "image/svg+xml" }));',
            '    a.download = "figure.svg"; a.click();',
            '  };',
            '});',
            '</scr' + 'ipt>', ''].join('\n');
        var html = '<p>' + escapeHtml(made) + ' Drawn by igv.js 3.8.10 (igv.org), as the Genome Browser draws it: ' +
            escapeHtml(fig.locus) + '.</p>' + inputsHtml(files) +
            '<h3>2. Index them</h3><p>With samtools and htslib (bgzip, tabix) installed:</p>' +
            codeBlock(['samtools faidx ' + genomeFile].concat(index).join('\n')) +
            '<h3>3. Open the view</h3><p>Put the attached <code>session.json</code> and <code>figure.html</code> in the same folder, then:</p>' +
            codeBlock('python3 -m http.server 8000\n# then open http://localhost:8000/figure.html, and press "Save as SVG"') +
            '<p>The view (<code>session.json</code>):</p>' + codeBlock(sessionText) +
            '<p>Cite igv.js: Robinson JT et al. (2023) igv.js: an embeddable JavaScript implementation of the Integrative ' +
            'Genomics Viewer (IGV). Bioinformatics 39: btac830.</p>';
        return { html: html, files: [textFile('session.json', sessionText, 'application/json'), textFile('figure.html', page, 'text/html')] };
    }

    // The recipe as a page under the figure's own page: its files attached, then its text.
    function saveRecipePage(NB, figurePage, title, fig) {
        return figureRecipe(fig).then(function (recipe) {
            if (!recipe) { return null; }
            var codePage;
            return postJSON(NB + '/notes', { parent_id: figurePage.id, title: 'How to make: ' + title, type: 'page' }).then(function (res) {
                codePage = res.note;
                var form = new FormData();
                recipe.files.forEach(function (f) { form.append('file', f); });
                return fetch(NB + '/notes/' + encodeURIComponent(codePage.id) + '/attachments', {
                    method: 'POST', credentials: 'same-origin', body: form
                }).then(function (r) {
                    return r.json().then(function (body) {
                        if (!r.ok) { throw new Error(body.error || ('the notebook said ' + r.status)); }
                        return body.attachments;
                    });
                });
            }).then(function (added) {
                var links = added.map(function (a, i) {
                    return '<a data-attachment="' + escapeHtml(a.id) + '">' + escapeHtml(a.name || recipe.files[i].name) + '</a>';
                }).join(', ');
                var html = '<h2>How to make this figure again</h2>' + recipe.html +
                    '<p><b>Attached:</b> ' + links + '.</p>';
                return fetch(NB + '/notes/' + encodeURIComponent(codePage.id), {
                    method: 'PUT', credentials: 'same-origin',
                    headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
                    body: JSON.stringify({ title: codePage.title, content: [{ type: 'text', html: html }], base_version: codePage.version })
                }).then(function (r) {
                    if (!r.ok) { return r.json().then(function (b) { throw new Error(b.error || ('the notebook said ' + r.status)); }); }
                    return codePage;
                });
            });
        });
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
            var lines = ['<p><b>' + (fig.kind === 'circos' ? 'A Circos plot from the Genome Browser:' : 'From the Genome Browser:') +
                '</b> ' + escapeHtml(where) + '.</p>'];
            if (fig.kind === 'analysis' && fig.about) {
                lines.push('<p>' + escapeHtml(fig.about) + (/deepTools/.test(fig.about)
                    ? ' Ram\u00edrez F et al. (2016) deepTools2: a next generation web server for deep-sequencing data analysis. ' +
                      'Nucleic Acids Research 44: W160\u2013W165.' : '') + '</p>');
            }
            if (fig.kind === 'circos') {
                lines.push('<p>Rings, outermost first, in bins of ' + fig.binSize.toLocaleString('en-US') +
                    ' bp. Drawn with Circos: Krzywinski M et al. (2009) Circos: an information aesthetic for comparative genomics. ' +
                    'Genome Research 19: 1639–1645.</p>');
            }
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
            // The code to make it again, as a page under it; the figure stands even if that fails.
            return saveRecipePage(NB, page, title, fig).then(function (codePage) {
                return codePage ? 'with the code to make it again on a page under it' : '';
            }, function (err) {
                return 'but the page with its code could not be made (' + err.message + ')';
            });
        }).then(function (codeNote) {
            var link = make('a', { href: NOTEBOOK + '/#note=' + encodeURIComponent(page.id), target: '_blank', rel: 'noopener', text: 'Open the page →' });
            var done = el('figure-done');
            done.replaceChildren(document.createTextNode('Saved to your Lab Notebook, under today: “' + title + '”' +
                (codeNote ? ', ' + codeNote : '') + '. '), link);
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
    // Sharing a genome or a track with other labs
    // ------------------------------------------------------------------

    function allLabs() {
        if (!state.labs) {
            state.labs = getJSON(API + '/labs').catch(function (err) { state.labs = null; throw err; });
        }
        return state.labs;
    }

    // ------------------------------------------------------------------
    // A track's details: its experiment, the Seed Bank plants it is tagged
    // with, and the Lab Notebook pages that link to it.
    // ------------------------------------------------------------------

    function sendJSON(method, url, body) {
        return fetch(url, {
            method: method, credentials: 'same-origin',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: body === undefined ? undefined : JSON.stringify(body)
        }).then(function (r) {
            return r.json().catch(function () { return {}; }).then(function (data) {
                if (!r.ok) { throw new Error((data && data.error) || ('HTTP ' + r.status)); }
                return data;
            });
        });
    }

    function detailError(text) {
        var box = el('detail-form').querySelector('.gb-error');
        box.textContent = text || '';
        box.hidden = !text;
    }

    function openDetails(t) {
        state.detailTrack = t;
        detailError('');
        el('detail-title').textContent = t.name;
        el('detail-intro').textContent = trackKind(t.type) + (t.can_change ? '' : ' · only its own lab changes these details');
        el('detail-experiment').value = t.experiment || '';
        el('detail-citation').value = t.citation || '';
        el('detail-doi').value = t.doi || '';
        el('detail-geo').value = t.geo || '';
        el('detail-saved').textContent = '';
        el('detail-edit').hidden = !t.can_change;
        renderFacts(t);
        el('detail-tag-add').hidden = !t.can_change;
        el('detail-remove').hidden = !t.can_change;
        el('detail-remove-confirm').hidden = true;
        el('detail-stock-search').value = '';
        el('detail-stock-results').replaceChildren();
        el('detail-tags').replaceChildren(make('li', { text: 'Loading…' }));
        el('detail-mentions').replaceChildren(make('li', { text: NOTEBOOK ? 'Loading…' : 'The Lab Notebook is not part of this suite.' }));
        fillExperimentNames();
        el('detail-dialog').showModal();
        loadDetailTags();
        loadMentions();
    }

    function loadDetailTags() {
        var t = state.detailTrack;
        getJSON(API + '/tracks/' + encodeURIComponent(t.id) + '/tags').then(function (body) {
            var tags = body.seed_stock_tags || [];
            var list = el('detail-tags');
            if (!tags.length) {
                list.replaceChildren(make('li', { class: 'gb-hint', text: 'Not tagged with a plant, line or seed lot yet.' }));
                return;
            }
            list.replaceChildren.apply(list, tags.map(function (tag) {
                var li = make('li', {}, [
                    make('a', { href: tag.link, target: '_blank', rel: 'noopener', text: tag.stock_id }),
                    make('span', { class: 'gb-share-why', text: tag.stock_kind === 'lot' ? 'seed lot' : tag.stock_kind })
                ]);
                if (t.can_change) {
                    li.appendChild(make('button', {
                        type: 'button', class: 'gb-link', text: 'Remove',
                        onclick: function () {
                            sendJSON('DELETE', API + '/tracks/' + encodeURIComponent(t.id) + '/tags/' + encodeURIComponent(tag.id))
                                .then(loadDetailTags).catch(function (err) { detailError(err.message); });
                        }
                    }));
                }
                return li;
            }));
        }).catch(function (err) { detailError('Could not read its tags: ' + err.message); });
    }

    function searchStocks() {
        var t = state.detailTrack;
        var q = el('detail-stock-search').value.trim();
        var list = el('detail-stock-results');
        if (!q) { list.replaceChildren(); return; }
        getJSON(API + '/seed-stocks?q=' + encodeURIComponent(q)).then(function (body) {
            if (el('detail-stock-search').value.trim() !== q) { return; }   // typed on since
            var stocks = body.stocks || [];
            if (!stocks.length) {
                list.replaceChildren(make('li', { class: 'gb-hint', text: 'Nothing in the Seed Bank you see matches “' + q + '”.' }));
                return;
            }
            list.replaceChildren.apply(list, stocks.map(function (s) {
                return make('li', {}, [
                    make('span', { class: 'gb-share-who', text: s.label }),
                    make('span', { class: 'gb-share-why', text: s.lab || '' }),
                    make('button', {
                        type: 'button', class: 'tt-btn', text: 'Tag',
                        onclick: function () {
                            postJSON(API + '/tracks/' + encodeURIComponent(t.id) + '/tags', { kind: 'seed_stock', stock_id: s.code })
                                .then(function () { el('detail-stock-search').value = ''; list.replaceChildren(); loadDetailTags(); })
                                .catch(function (err) { detailError(err.message); });
                        }
                    })
                ]);
            }));
        }).catch(function (err) { detailError('Could not search the Seed Bank: ' + err.message); });
    }

    function loadMentions() {
        if (!NOTEBOOK) { return; }
        var t = state.detailTrack;
        getJSON(NOTEBOOK + '/api/mentions?type=track&id=' + encodeURIComponent(t.id)).then(function (body) {
            var notes = body.notes || [];
            var list = el('detail-mentions');
            if (!notes.length) {
                list.replaceChildren(make('li', { class: 'gb-hint', text: 'No notebook page you can open links to it yet.' }));
                return;
            }
            list.replaceChildren.apply(list, notes.map(function (n) {
                return make('li', {}, [
                    make('a', { href: NOTEBOOK + '/#note=' + encodeURIComponent(n.id), target: '_blank', rel: 'noopener', text: n.title || 'Untitled' }),
                    make('span', { class: 'gb-share-why', text: n.owner ? 'by ' + n.owner : '' })
                ]);
            }));
        }).catch(function () {
            el('detail-mentions').replaceChildren(make('li', { class: 'gb-hint', text: 'The Lab Notebook could not be asked.' }));
        });
    }

    // What a member who cannot change the track sees: its experiment and where it comes from, linked.
    function renderFacts(t) {
        var facts = el('detail-facts');
        var rows = [];
        if (!t.can_change) {
            if (t.experiment) { rows.push(['Experiment', [document.createTextNode(t.experiment)]]); }
            if (t.citation) { rows.push(['Paper', [document.createTextNode(t.citation)]]); }
            if (t.doi) {
                rows.push(['DOI', [make('a', { href: 'https://doi.org/' + t.doi, target: '_blank', rel: 'noopener', text: t.doi })]]);
            }
            if (t.geo) {
                rows.push(['GEO', t.geo.split(' ').map(function (acc) {
                    return make('a', { href: geoLink(acc), target: '_blank', rel: 'noopener', text: acc });
                })]);
            }
        }
        facts.replaceChildren.apply(facts, [].concat.apply([], rows.map(function (r) {
            return [make('dt', { text: r[0] }), make('dd', {}, r[1])];
        })));
        facts.hidden = !rows.length;
    }

    var DETAIL_FIELDS = ['experiment', 'citation', 'doi', 'geo'];

    function saveDetails() {
        var t = state.detailTrack;
        var body = {};
        DETAIL_FIELDS.forEach(function (f) { body[f] = el('detail-' + f).value; });
        el('detail-saved').textContent = '';
        sendJSON('PATCH', API + '/tracks/' + encodeURIComponent(t.id), body)
            .then(function (saved) {
                DETAIL_FIELDS.forEach(function (f) {
                    t[f] = saved[f] || null;
                    el('detail-' + f).value = t[f] || '';
                });
                state.tracks.forEach(function (u) {
                    if (u.id === t.id) { DETAIL_FIELDS.forEach(function (f) { u[f] = t[f]; }); }
                });
                detailError('');
                el('detail-saved').textContent = 'Saved.';
                renderTracks();
            }).catch(function (err) { detailError(err.message); });
    }

    function removeTrack() {
        var t = state.detailTrack;
        el('detail-remove-yes').disabled = true;
        sendJSON('DELETE', API + '/tracks/' + encodeURIComponent(t.id)).then(function () {
            state.tracks = state.tracks.filter(function (u) { return u.id !== t.id; });
            state.selected.delete(t.id);
            el('detail-dialog').close();
            renderTracks();
            updateOpenHint();
        }).catch(function (err) {
            detailError('Could not remove it: ' + err.message);
        }).then(function () { el('detail-remove-yes').disabled = false; });
    }

    // ------------------------------------------------------------------
    // Marking a sequence on the genome: a guide, a primer or anything else
    // the annotation does not name, found on both strands and kept as a track.
    // ------------------------------------------------------------------

    function openMotifDialog() {
        if (!state.genome) { return; }
        el('motif-form').querySelector('.gb-error').hidden = true;
        el('motif-genome-name').textContent = state.genome.name;
        el('motif-owner').textContent = state.genome.suite && state.me && state.me.is_admin
            ? 'As an administrator, your tracks on a suite genome are shared with everyone.'
            : ownerNote();
        el('motif-summary').hidden = true;
        el('motif-go').disabled = false;
        fillExperimentNames();
        el('motif-dialog').showModal();
        el('motif-name').focus();
    }

    function submitMotif(e) {
        e.preventDefault();
        var form = el('motif-form');
        form.querySelector('.gb-error').hidden = true;
        var name = el('motif-name').value.trim();
        if (!name) { formError(form, ['Give the track a name.']); el('motif-name').focus(); return; }
        if (!el('motif-sequences').value.trim()) {
            formError(form, ['Type at least one sequence to look for.']);
            el('motif-sequences').focus();
            return;
        }
        el('motif-go').disabled = true;
        el('motif-summary').hidden = false;
        el('motif-summary').textContent = 'Looking through the genome…';
        postJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/tracks/from-sequence', {
            name: name,
            sequences: el('motif-sequences').value,
            experiment: el('motif-experiment').value,
            citation: el('motif-citation').value,
            doi: el('motif-doi').value,
            geo: el('motif-geo').value
        }).then(function (body) {
            el('motif-dialog').close();
            el('motif-sequences').value = '';
            el('motif-name').value = '';
            // The new track joins the list, ticked, as an upload's does.
            return getJSON(API + '/genomes/' + encodeURIComponent(state.genome.id) + '/tracks').then(function (tracks) {
                state.tracks = tracks;
                state.selected.add(body.id);
                renderTracks();
                updateOpenHint();
                var per = body.per_sequence || {};
                var parts = Object.keys(per).map(function (k) { return k + ': ' + per[k]; });
                el('track-search-count').textContent = 'Marked ' + body.found + ' place'
                    + (body.found === 1 ? '' : 's') + (parts.length > 1 ? ' (' + parts.join(', ') + ')' : '')
                    + ((body.capped || []).length
                        ? ' — ' + body.capped.join(', ') + ' is all over the genome, so only the first of it is marked'
                        : '');
            });
        }).catch(function (err) {
            el('motif-summary').hidden = true;
            el('motif-go').disabled = false;
            formError(form, [err.message]);
        });
    }

    // kind is 'genome' or 'track'; item is the genome or track from the lists.
    function openShare(kind, item) {
        state.share = { kind: kind, item: item };
        var form = el('share-form');
        form.querySelector('.gb-error').hidden = true;
        el('share-done').hidden = true;
        el('share-title').textContent = 'Share “' + item.name + '”';
        el('share-list').replaceChildren(make('li', { text: 'Loading…' }));
        el('share-add').hidden = true;
        el('share-intro').textContent = '';
        el('share-dialog').showModal();
        loadShare();
    }

    function sharePath() {
        var s = state.share;
        return API + (s.kind === 'genome' ? '/genomes/' : '/tracks/') + encodeURIComponent(s.item.id);
    }

    function loadShare() {
        var s = state.share;
        return Promise.all([getJSON(sharePath() + '/shares'), allLabs()]).then(function (res) {
            s.info = res[0];
            renderShare(res[0], res[1]);
        }).catch(function (err) {
            formError(el('share-form'), ['Could not read who sees it: ' + err.message]);
        });
    }

    function renderShare(info, labs) {
        var s = state.share;
        var list = el('share-list');
        var rows = [];
        function row(who, why, removable) {
            var li = make('li', {}, [make('span', { class: 'gb-share-who', text: who }), make('span', { class: 'gb-share-why', text: why })]);
            if (removable) {
                li.appendChild(make('button', {
                    type: 'button', class: 'gb-link', text: 'Remove',
                    onclick: function () { changeShare(removable, false); }
                }));
            }
            return li;
        }
        if (info.suite) {
            rows.push(row('Everyone in the suite', s.kind === 'genome' ? 'a suite reference genome' : 'an administrator\'s track on a suite genome', null));
        } else if (s.kind === 'genome') {
            info['default'].forEach(function (d) {
                // The owner lab always sees it; its chair group, by default, can be taken away (and added back).
                rows.push(row(d.name, d.role === 'owner' ? 'owner lab' : 'its chair group, by default',
                    info.can_change && d.role === 'default' ? { type: 'chair_group', id: d.id, name: d.name } : null));
            });
            info.shares.forEach(function (x) {
                var why = (x.type === 'chair_group' ? 'chair group, ' : '') +
                    (x.auto ? 'through a shared track' : 'shared');
                rows.push(row(x.name, why, info.can_change && !x.auto ? x : null));
            });
        } else {
            rows.push(row(info.owner_lab.name, 'owner lab', null));
            info.shares.forEach(function (x) {
                rows.push(row(x.name, x.type === 'chair_group' ? 'its chair group, shared' : 'shared', info.can_change ? x : null));
            });
        }
        list.replaceChildren.apply(list, rows);

        var intro = '';
        if (info.suite) {
            intro = 'Everyone in the suite sees this already.';
        } else if (!info.can_change) {
            intro = 'Only members of ' + info.owner_lab.name + ' can change who sees this.';
        } else if (s.kind === 'track') {
            intro = 'A lab you share this track with also gets its genome, if it cannot see it yet, so the track can be opened.';
        } else {
            intro = 'Labs and chair groups you share this genome with can open it and add their own tracks to it.';
        }
        el('share-intro').textContent = intro;
        el('share-add').hidden = info.suite || !info.can_change;
        if (el('share-add').hidden) { return; }

        // What it can still be shared with.
        var taken = new Set(info.shares.map(function (x) { return x.type + ':' + x.id; }));
        var owner = info.owner_lab && info.owner_lab.id;
        if (s.kind === 'genome') { info['default'].forEach(function (d) { taken.add(d.type + ':' + d.id); }); }
        var select = el('share-target');
        var options = [make('option', { value: '', text: 'Choose a lab or chair group…' })];
        var groupOptions = [];
        if (s.kind === 'genome') {
            var ownGroup = info.owner_chair_group && info.owner_chair_group.id;
            labs.chair_groups.forEach(function (g) {
                if (!taken.has('chair_group:' + g.id)) {
                    groupOptions.push(make('option', { value: 'chair_group:' + g.id,
                        text: g.name + (g.id === ownGroup ? ' (your chair group, every lab in it)' : ' (every lab in it)') }));
                }
            });
        } else if (info.owner_chair_group && !taken.has('chair_group:' + info.owner_chair_group.id)) {
            groupOptions.push(make('option', { value: 'chair_group:' + info.owner_chair_group.id, text: info.owner_chair_group.name + ' (your chair group, every lab in it)' }));
        }
        if (groupOptions.length) { options.push(make('optgroup', { label: 'Chair groups' }, groupOptions)); }
        var labOptions = labs.labs.filter(function (l) { return l.id !== owner && !taken.has('lab:' + l.id); })
            .map(function (l) { return make('option', { value: 'lab:' + l.id, text: l.name + (l.chair_group ? ' (' + l.chair_group + ')' : '') }); });
        if (labOptions.length) { options.push(make('optgroup', { label: 'Labs' }, labOptions)); }
        select.replaceChildren.apply(select, options);
        el('share-note').textContent = labOptions.length || groupOptions.length ? '' : 'It is shared with every lab already.';
    }

    function changeShare(target, grant) {
        var s = state.share;
        var body;
        if (s.kind === 'track') {
            body = target.type === 'chair_group' ? { chair_group: true } : { labs: [target.id] };
        } else {
            body = target.type === 'chair_group' ? { chair_group_id: target.id } : { labs: [target.id] };
        }
        var form = el('share-form');
        form.querySelector('.gb-error').hidden = true;
        el('share-done').hidden = true;
        return fetch(sharePath() + '/share', {
            method: grant ? 'POST' : 'DELETE', credentials: 'same-origin',
            headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
            body: JSON.stringify(body)
        }).then(function (r) {
            return r.json().then(function (b) {
                if (!r.ok) { throw new Error(b.error || ('HTTP ' + r.status)); }
                return b;
            });
        }).then(function (b) {
            var notes = [];
            var name = target.name || target.id;
            notes.push(grant ? 'Shared with ' + name + '.' : name + ' no longer sees it.');
            if (b.genome_shared_with && b.genome_shared_with.length) {
                notes.push('So that the track can be opened, the genome “' + (state.genome ? state.genome.name : '') +
                    '” is now shared with ' + b.genome_shared_with.map(function (g) { return g.name; }).join(', ') + ' too.');
            }
            if (b.kept && b.kept.length) {
                notes.push(b.kept.map(function (k) { return k.name; }).join(', ') +
                    ' still sees the genome, because a track on it is shared with them.');
            }
            el('share-done').textContent = notes.join(' ');
            el('share-done').hidden = false;
            state.otherTracks = null;
            return loadShare();
        }).catch(function (err) {
            formError(form, ['Could not change it: ' + err.message]);
        });
    }

    // ------------------------------------------------------------------
    // Routing
    // ------------------------------------------------------------------

    function route() {
        var params = new URLSearchParams(location.search);
        var genome = params.get('genome');
        var track = params.get('track');
        if (track && !genome) { openTrack(track); } else if (genome) { openBrowser(genome, null); } else { showChoose(null, true); }
    }

    // A link to one track (from a notebook page): its genome opens as the
    // member last had it, with that track on too.
    function openTrack(trackId) {
        getJSON(API + '/tracks/' + encodeURIComponent(trackId)).then(function (t) {
            if (!t.genome_id) { throw new Error('no genome'); }
            history.replaceState({ genome: t.genome_id }, '', PREFIX + '/?genome=' + encodeURIComponent(t.genome_id));
            return getJSON(API + '/view?genome_id=' + encodeURIComponent(t.genome_id)).catch(function () { return {}; })
                .then(function (view) {
                    view = view || {};
                    var enabled = (view.enabled || []).filter(function (id) { return id !== t.id; });
                    enabled.push(t.id);
                    openBrowser(t.genome_id, { locus: view.locus || '', enabled: enabled, options: view.options || {} });
                });
        }).catch(function () {
            history.replaceState({}, '', PREFIX + '/');
            showChoose(null, true);
            var note = el('choose-note');
            note.textContent = 'That track no longer exists, or it is not shared with you. Ask its lab to share it.';
            note.hidden = false;
        });
    }

    el('choose-form').addEventListener('submit', function (e) { e.preventDefault(); openFromChoice(); });
    el('track-sort').value = trackSort();
    el('track-sort').addEventListener('change', function () {
        try { localStorage.setItem('gb-track-sort', el('track-sort').value); } catch (e) { /* kept for this page only */ }
        renderTracks();
    });
    el('detail-close').addEventListener('click', function () { el('detail-dialog').close(); });
    el('detail-save').addEventListener('click', saveDetails);
    DETAIL_FIELDS.forEach(function (f) {
        el('detail-' + f).addEventListener('keydown', function (e) {
            if (e.key === 'Enter') { e.preventDefault(); saveDetails(); }
        });
        el('detail-' + f).addEventListener('input', function () { el('detail-saved').textContent = ''; });
    });
    el('detail-remove').addEventListener('click', function () {
        el('detail-remove-name').textContent = state.detailTrack.name;
        el('detail-remove-confirm').hidden = false;
        el('detail-remove').hidden = true;
        el('detail-remove-cancel').focus();
    });
    el('detail-remove-cancel').addEventListener('click', function () {
        el('detail-remove-confirm').hidden = true;
        el('detail-remove').hidden = false;
    });
    el('detail-remove-yes').addEventListener('click', removeTrack);
    var stockTimer = null;
    el('detail-stock-search').addEventListener('input', function () {
        clearTimeout(stockTimer);
        stockTimer = setTimeout(searchStocks, 250);
    });
    el('detail-stock-search').addEventListener('keydown', function (e) {
        if (e.key === 'Enter') { e.preventDefault(); searchStocks(); }
    });

    el('track-search').addEventListener('input', function () {
        state.trackQuery = el('track-search').value;
        renderTracks();
    });
    el('track-search').addEventListener('keydown', function (e) {
        if (e.key === 'Escape' && el('track-search').value) {
            e.preventDefault();
            el('track-search').value = '';
            state.trackQuery = '';
            renderTracks();
        }
    });
    el('tracks-shown').addEventListener('click', tickShown);
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

    el('share-genome').addEventListener('click', function () { if (state.genome) { openShare('genome', state.genome); } });
    el('share-go').addEventListener('click', function () {
        var value = el('share-target').value;
        if (!value) { formError(el('share-form'), ['Choose a lab or chair group first.']); return; }
        var cut = value.indexOf(':');
        var option = el('share-target').selectedOptions[0];
        changeShare({ type: value.slice(0, cut), id: value.slice(cut + 1), name: option ? option.textContent.replace(/ \(.*\)$/, '') : '' }, true);
    });
    el('share-close').addEventListener('click', function () {
        el('share-dialog').close();
        // Visibility labels in the list may have changed.
        if (state.genome && !el('choose').hidden) { chooseGenome(state.genome.id, true); }
    });

    el('add-genome').addEventListener('click', openGenomeDialog);
    el('add-track').addEventListener('click', openTrackDialog);
    el('mark-sequence').addEventListener('click', openMotifDialog);
    el('motif-cancel').addEventListener('click', function () { el('motif-dialog').close(); });
    el('motif-form').addEventListener('submit', submitMotif);
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
    // With a search typed, Select all ticks the tracks it shows (and keeps the rest).
    el('tp-all').addEventListener('click', function () {
        var words = searchWords(el('tp-search').value);
        panelTracks().forEach(function (t) {
            if (panelTrackMatches(t, words)) { state.panelChecked.add(t.config.trellisId); }
        });
        renderTrackPanel();
    });
    el('tp-search').addEventListener('input', renderTrackPanel);
    el('tp-none').addEventListener('click', function () { state.panelChecked = new Set(); renderTrackPanel(); });
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-height]'), function (b) {
        b.addEventListener('click', function () { setHeight(b.getAttribute('data-height')); });
    });
    el('tp-height-form').addEventListener('submit', function (e) { e.preventDefault(); setHeight(el('tp-height').value); });
    el('tp-swatches').replaceChildren.apply(el('tp-swatches'), PALETTE.map(function (c) {
        var b = make('button', { type: 'button', title: c, 'aria-label': 'Colour ' + c });
        b.style.background = c;
        b.addEventListener('click', function () { el('tp-hex').value = c; el('tp-color').value = c; setColour(c); });
        return b;
    }));
    el('tp-color').addEventListener('change', function () {
        el('tp-hex').value = el('tp-color').value;
        setColour(el('tp-color').value);
    });
    // A hex code: #1f5fa8, 1f5fa8 or the short #15a; anything else is refused with an example.
    el('tp-hex-form').addEventListener('submit', function (e) {
        e.preventDefault();
        var raw = el('tp-hex').value.trim().replace(/^#/, '');
        if (!/^([0-9a-f]{3}|[0-9a-f]{6})$/i.test(raw)) {
            el('tp-error').textContent = 'A hex colour is six (or three) of 0-9 and a-f, such as #1f5fa8.';
            el('tp-error').hidden = false;
            return;
        }
        if (raw.length === 3) { raw = raw.split('').map(function (c) { return c + c; }).join(''); }
        var hex = '#' + raw.toLowerCase();
        el('tp-hex').value = hex;
        el('tp-color').value = hex;
        setColour(hex);
    });
    el('tp-color-reset').addEventListener('click', resetColour);
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-display]'), function (b) {
        b.addEventListener('click', function () { setDisplay(b.getAttribute('data-display')); });
    });
    Array.prototype.forEach.call(document.querySelectorAll('#tp-controls [data-scale]'), function (b) {
        b.addEventListener('click', function () { setScale(b.getAttribute('data-scale')); });
    });
    el('tp-range-form').addEventListener('submit', setRange);
    el('tp-remove').addEventListener('click', removeTicked);

    el('figure-open').addEventListener('click', function () { openFigureDialog(); });
    el('circos-open').addEventListener('click', openCircosDialog);
    el('analyse-open').addEventListener('click', openAnalyseDialog);
    el('analyse-form').addEventListener('submit', runAnalysis);
    el('analyse-close').addEventListener('click', function () { el('analyse-dialog').close(); });
    el('circos-form').addEventListener('submit', drawCircos);
    el('circos-close').addEventListener('click', function () { el('circos-dialog').close(); });
    el('figure-circos-files').addEventListener('click', downloadCircosFiles);
    el('figure-form').addEventListener('submit', saveFigureToNotebook);
    el('figure-svg').addEventListener('click', downloadFigureSvg);
    el('figure-close').addEventListener('click', function () { el('figure-dialog').close(); });

    el('igv-container').addEventListener('click', onClickCapture, true);
    // Where each press on the tracks began, to tell a click from a drag.
    el('igv-container').addEventListener('mousedown', function (e) {
        state.downAt = { x: e.clientX, y: e.clientY };
    }, true);
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
