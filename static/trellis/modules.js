/*
 * The suite's other tools, as a row of icons in the top bar, so every tool reaches every other one.
 *
 * The list is the suite's own (landing/tools.json, served at /tools.json beside the landing page), so a tool
 * added there appears here without each tool being changed. Running on its own, outside the suite, that file
 * is not there and nothing is added.
 *
 * Include it after the header: <script src="…/static/trellis/modules.js" defer></script>
 */
(function () {
  'use strict';

  function suiteRoot() {
    // The "Trellis Tools" link in the bar points at the suite's own page: everything else sits beside it.
    const home = document.querySelector('.tt-bar .tt-home');
    return home ? new URL(home.getAttribute('href') || '/', window.location.href) : null;
  }

  function here(tool, root) {
    // The tool this page belongs to: its name in the bar, else the path we are on.
    const mine = (document.querySelector('.tt-bar .tt-module')?.textContent || '').trim().toLowerCase();
    if (mine && tool.name.toLowerCase() === mine) return true;
    const path = new URL(tool.href, root).pathname.replace(/\/$/, '');
    return !!path && window.location.pathname.replace(/\/$/, '').startsWith(path);
  }

  async function draw() {
    const bar = document.querySelector('.tt-bar');
    const root = suiteRoot();
    if (!bar || !root || document.querySelector('.tt-modules')) return;
    let tools;
    try {
      const res = await fetch(new URL('tools.json', root), { credentials: 'same-origin' });
      if (!res.ok) return;                       // on its own, or signed out: no switcher
      tools = (await res.json()).tools;
    } catch (err) {
      return;
    }
    if (!Array.isArray(tools) || !tools.length) return;
    const nav = document.createElement('nav');
    nav.className = 'tt-modules';
    nav.setAttribute('aria-label', 'The suite\'s tools');
    for (const tool of tools) {
      if (!tool || !tool.href || !tool.name || tool.status === 'planned') continue;
      const mine = here(tool, root);
      const a = document.createElement('a');
      a.className = `tt-module-link${mine ? ' here' : ''}`;
      a.href = new URL(tool.href, root).href;
      a.title = mine ? `${tool.name} (you are here)` : tool.name;
      if (mine) a.setAttribute('aria-current', 'page');
      if (tool.image) {
        const img = document.createElement('img');
        img.src = new URL(tool.image, root).href;
        img.alt = '';
        a.appendChild(img);
      }
      const label = document.createElement('span');
      label.textContent = tool.name;
      a.appendChild(label);
      nav.appendChild(a);
    }
    // After the tool's own name and any subtitle, before its tabs.
    const after = bar.querySelector('.tt-sub') || bar.querySelector('.tt-module');
    if (after) after.insertAdjacentElement('afterend', nav); else bar.appendChild(nav);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', draw);
  else draw();
}());
