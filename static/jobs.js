/* My jobs (jobs.py): lists the member's jobs with progress bars, asking again every few seconds
   while any is waiting or running. Finished ones open in the Genome Browser. */
(function () {
    'use strict';

    var MOUNT = document.body.getAttribute('data-mount') || '/genome';
    var BUSY_EVERY = 3000, IDLE_EVERY = 20000;
    var timer = null;

    function el(id) { return document.getElementById(id); }

    function make(tag, attrs, children) {
        var node = document.createElement(tag);
        Object.keys(attrs || {}).forEach(function (k) {
            if (k === 'text') { node.textContent = attrs[k]; } else if (k === 'class') { node.className = attrs[k]; } else { node.setAttribute(k, attrs[k]); }
        });
        (children || []).forEach(function (c) { if (c) { node.appendChild(c); } });
        return node;
    }

    // SQLite keeps UTC as 'YYYY-MM-DD HH:MM:SS'.
    function when(stamp) { return stamp ? new Date(stamp.replace(' ', 'T') + 'Z') : null; }

    function duration(seconds) {
        if (seconds === null || seconds === undefined || isNaN(seconds)) { return ''; }
        seconds = Math.max(0, Math.round(seconds));
        if (seconds < 90) { return seconds + ' s'; }
        if (seconds < 5400) { return Math.round(seconds / 60) + ' min'; }
        return (seconds / 3600).toFixed(1) + ' h';
    }

    function ago(date) {
        if (!date) { return ''; }
        var s = (Date.now() - date.getTime()) / 1000;
        if (s < 60) { return 'just now'; }
        if (s < 3600) { return Math.round(s / 60) + ' min ago'; }
        if (s < 86400) { return Math.round(s / 3600) + ' h ago'; }
        return date.toLocaleDateString(undefined, { day: 'numeric', month: 'short' }) + ', ' +
            date.toLocaleTimeString(undefined, { hour: '2-digit', minute: '2-digit' });
    }

    var STATES = { queued: 'Waiting', running: 'Running', done: 'Done', failed: 'Failed' };

    function detail(job) {
        if (job.state === 'queued') {
            return job.ahead ? 'Waiting for ' + job.ahead + ' job' + (job.ahead === 1 ? '' : 's') + ' ahead of it.'
                : 'Next in line.';
        }
        if (job.state === 'running') {
            var parts = [];
            if (job.step) { parts.push(job.step + '…'); }
            if (typeof job.progress === 'number') { parts.push(Math.round(job.progress * 100) + '%'); }
            var started = when(job.started_at);
            if (started) { parts.push('running ' + duration((Date.now() - started.getTime()) / 1000)); }
            if (job.estimate_seconds) { parts.push('of about ' + duration(job.estimate_seconds)); }
            return parts.join(' · ');
        }
        if (job.state === 'done') { return 'Took ' + duration(job.seconds) + '.'; }
        return '';
    }

    function bar(job) {
        var known = typeof job.progress === 'number';
        var fill = make('span');
        if (known) { fill.style.width = (Math.max(0, Math.min(1, job.progress)) * 100).toFixed(1) + '%'; }
        if (job.state === 'failed') { fill.style.width = '100%'; }
        var outer = make('div', { class: 'job-bar' + (!known && job.state === 'running' ? ' unknown' : ''), role: 'progressbar',
            'aria-label': job.what, 'aria-valuemin': '0', 'aria-valuemax': '100' }, [fill]);
        if (known) { outer.setAttribute('aria-valuenow', String(Math.round(job.progress * 100))); }
        return outer;
    }

    function card(job) {
        var actions = [];
        if (job.open_url) {
            // In the dashboard this page is in a frame: open the browser in the whole window.
            actions.push(make('a', { class: 'job-open', href: job.open_url, target: '_top',
                text: job.source === 'analysis' && !job.track_id ? 'Open the result' : 'Open the track' }));
        }
        var started = when(job.created_at);
        return make('li', { class: 'job', 'data-state': job.state }, [
            make('div', { class: 'job-top' }, [
                make('span', { class: 'job-what', text: job.what }),
                make('span', { class: 'job-state', text: STATES[job.state] || job.state }),
                make('span', { class: 'job-when', text: started ? 'Asked ' + ago(started) : '' })
            ]),
            make('div', { class: 'job-about' }, [
                document.createTextNode(job.about || ''),
                make('span', { class: 'job-genome', text: (job.about ? ' · ' : '') + (job.genome && job.genome.name || '') })
            ]),
            job.state === 'done' ? null : bar(job),
            job.state === 'failed' && job.error ? make('div', { class: 'job-error', text: job.error }) : null,
            make('div', { class: 'job-detail' }, [
                make('span', { text: detail(job) }),
                actions.length ? make('span', { class: 'job-actions' }, actions) : null
            ])
        ]);
    }

    function render(data) {
        el('keep-days').textContent = String(data.keep_days || 7);
        var jobs = data.jobs || [];
        el('jobs-list').replaceChildren.apply(el('jobs-list'), jobs.map(card));
        el('jobs-empty').hidden = jobs.length > 0;
    }

    function load() {
        clearTimeout(timer);
        fetch(MOUNT + '/api/jobs', { credentials: 'same-origin', cache: 'no-store', headers: { Accept: 'application/json' } })
            .then(function (res) {
                if (res.status === 401) { throw new Error('Sign in again to see your jobs.'); }
                if (!res.ok) { throw new Error('Could not load your jobs (' + res.status + ').'); }
                return res.json();
            })
            .then(function (data) {
                el('jobs-error').hidden = true;
                render(data);
                schedule(data.busy ? BUSY_EVERY : IDLE_EVERY);
            })
            .catch(function (err) {
                el('jobs-error').textContent = err.message || String(err);
                el('jobs-error').hidden = false;
                schedule(IDLE_EVERY);
            });
    }

    function schedule(ms) {
        clearTimeout(timer);
        if (!document.hidden) { timer = setTimeout(load, ms); }
    }

    // Not asking while nobody is looking; asking at once on coming back.
    document.addEventListener('visibilitychange', function () {
        if (document.hidden) { clearTimeout(timer); } else { load(); }
    });
    load();
}());
