/* OligoMet Studio tab (client side). Draws what R/studio.R sends:
     oms-studio-load   { mode, oligo, samples[], species[] }
     oms-studio-detail { id, sample, z0, env[], iso[], acquired[], theoretical[], score, has_spectrum }
     oms-studio-empty  { has_library }
   and reports clicks back through the Shiny inputs studio_pick, studio_demo,
   studio_mine and studio_nav. No data is made up here. */
(function () {
  'use strict';
  var S = { pl: null, sel: null, si: 0, filter: 'all', sort: { key: 'pct', dir: -1 }, detail: {}, hl: null,
            ions: { a: true, 'a-B': true, b: true, 'b-B': true, w: true, y: true }, root: null, view: null };

  var CLASSES = [
    { key: 'parent', label: 'Intact strand', color: 'var(--c-parent)' },
    { key: 'gal', label: 'Conjugate catabolites', color: 'var(--c-gal)' },
    { key: '3p', label: '3′ truncations', color: 'var(--c-3p)' },
    { key: '5p', label: '5′ truncations', color: 'var(--c-5p)' },
    { key: 'endo', label: 'Internal cleavage fragments', color: 'var(--c-endo)' },
    { key: 'ps', label: 'PS → PO desulfurisation', color: 'var(--c-ps)' },
    { key: 'other', label: 'Other', color: 'var(--c-other)' }
  ];
  var KIND_TEXT = {
    parent: 'Intact strand', conj_catab: 'Conjugate catabolite', exo_3p: '3′ end truncation', exo_5p: '5′ end truncation',
    endo_5frag: 'Internal cleavage, 5′ fragment', endo_5frag_p: 'Internal cleavage, 5′ fragment (3′-phosphate)',
    endo_3frag: 'Internal cleavage, 3′ fragment', endo_3frag_p: 'Internal cleavage, 3′ fragment (5′-phosphate)'
  };
  var FIVE = { a: 1, 'a-B': 1, b: 1, 'b-B': 1, c: 1, d: 1 };
  var LRE = /^(a-B|b-B|a|b|c|d|w|x|y|z)(\d+)\^(\d+)-$/;

  /* ---------- helpers ---------- */
  function q(s) { return S.root.querySelector(s); }
  function esc(t) { return String(t == null ? '' : t).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
  function fmt(x, d) { if (x == null || isNaN(x)) return '—'; d = d == null ? 2 : d; return Number(x).toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d }); }
  function sup(n) { return String(n).replace(/\d/g, function (d) { return '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]; }) + '⁻'; }
  function sv(tag, attrs) { var e = document.createElementNS('http://www.w3.org/2000/svg', tag); for (var k in attrs) e.setAttribute(k, attrs[k]); return e; }
  function svgText(parent, x, y, txt, attrs) { var t = sv('text', Object.assign({ x: x, y: y }, attrs || {})); t.textContent = txt; parent.appendChild(t); return t; }
  function clsOf(sp) { return CLASSES.filter(function (c) { return c.key === sp.cls; })[0] || CLASSES[CLASSES.length - 1]; }
  function cur() { return S.pl.species.filter(function (s) { return s.id === S.sel; })[0]; }
  function hasSamples() { return S.pl && S.pl.samples.length > 0; }
  function sampleId() { return hasSamples() ? S.pl.samples[S.si].id : ''; }
  function val(arr, i) { return arr && arr[i] != null ? arr[i] : null; }
  function toast(msg) { /* uses Shiny's own notification so it matches the app */
    if (window.Shiny && Shiny.notifications) Shiny.notifications.show({ html: esc(msg), type: 'message', duration: 3500 }); }
  function send(name, v) { if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue(name, v, { priority: 'event' }); }

  /* ---------- skeletons ---------- */
  function themeButton() { return '<button class="btn" data-act="theme" type="button" aria-label="Toggle light and dark theme">Theme</button>'; }
  function renderEmpty(msg) {
    S.view = 'empty';
    S.root.innerHTML = '<div class="top"><div class="brand"><div class="mark" aria-hidden="true">Om</div><div><h1>OligoMet Studio</h1><p>Metabolite workbench</p></div></div><div class="top-actions">' + themeButton() + '</div></div>' +
      '<div class="empty"><h2>' + esc(msg.title) + '</h2><p>' + esc(msg.body) + '</p><div class="row">' +
      '<button class="btn primary" data-act="nav" data-to="Library Generation" type="button">Open Library Generation</button>' +
      '<button class="btn" data-act="demo" type="button">Explore with demo data</button></div></div>';
    wire();
  }
  function renderShell() {
    S.view = 'studio';
    S.root.innerHTML =
      '<header class="top"><div class="brand"><div class="mark" aria-hidden="true">Om</div><div><h1>OligoMet Studio</h1><p>Sequence, species, MS/MS and signal in one view</p></div></div>' +
      '<div class="top-actions"><span class="pill">Oligo <b data-f="oname"></b></span>' + themeButton() + '</div></header>' +
      '<div class="notice" data-f="notice"></div>' +
      '<main class="studio">' +
      '<aside class="pane left" aria-label="Constructs and samples">' +
        '<div class="block"><div class="eyebrow">Constructs</div><div class="tree">' +
          '<button class="node active" type="button"><span class="nm"><span class="dot"></span><span data-f="nodeName"></span></span><small data-f="nodeSub"></small></button>' +
          '<button class="node" type="button" disabled><span class="nm"><span class="dot off"></span>Antisense strand</span><small>planned</small></button>' +
          '<button class="node" type="button" disabled><span class="nm"><span class="dot off"></span>Internal standard</span><small>not shown</small></button>' +
        '</div><p class="hint">The engine analyses one strand per run today. The tree is where duplex support would live.</p></div>' +
        '<div class="block"><div class="block-head"><div class="eyebrow">Samples</div><span class="hint" data-f="sampleHint"></span></div><div class="samples" data-f="samples" role="listbox" aria-label="Samples"></div></div>' +
        '<div class="block" data-f="actions"></div>' +
      '</aside>' +
      '<section class="pane center" aria-label="Analysis">' +
        '<div class="card"><div class="block-head"><div class="eyebrow">Sequence map, 5′ → 3′</div><span class="mono" data-f="mapNote" style="font-size:11px;color:var(--muted)"></span></div>' +
          '<div class="ribbon-row" data-f="ribbon"></div>' +
          '<div class="legend"><span><i style="background:var(--ome)"></i>2′-OMe</span><span><i style="background:var(--fl)"></i>2′-F</span><span><span class="star">*</span>phosphorothioate link</span>' +
          '<span><i style="background:var(--gal)"></i>conjugate</span><span><i style="background:var(--noise);opacity:.4"></i>not in this species</span><span><i style="box-shadow:0 0 0 2px var(--accent) inset;background:transparent"></i>covered by hovered ion</span></div></div>' +
        '<div class="card"><div class="block-head" style="margin-bottom:8px"><div class="eyebrow">Predicted species</div><span class="mono" data-f="tblCount" style="font-size:11px;color:var(--muted)"></span></div>' +
          '<div class="chips" data-f="chips" role="group" aria-label="Filter by class" style="margin-bottom:8px"></div>' +
          '<div class="tbl-wrap"><table><thead><tr><th><button data-sort="id" type="button">ID</button></th><th>Species</th><th class="num"><button data-sort="mass" type="button">Mass (Da)</button></th>' +
          '<th class="num" data-r>Δppm</th><th class="num" data-r><button data-sort="cov" type="button">MS/MS cov.</button></th><th class="num" data-r><button data-sort="pct" type="button">% of signal</button></th></tr></thead><tbody data-f="tbody"></tbody></table></div></div>' +
        '<div class="card"><div class="block-head" style="margin-bottom:6px;align-items:flex-start;flex-wrap:wrap"><div><div class="eyebrow">MS/MS mirror</div><div data-f="mirrorTitle" style="font-weight:600"></div></div><div class="stats" data-f="mirrorStats"></div></div>' +
          '<div class="toggles" data-f="toggles" role="group" aria-label="Ion series"></div>' +
          '<div class="plot-wrap" style="margin-top:6px"><svg class="plot" data-f="mirror" viewBox="0 0 760 330" role="img" aria-label="Acquired spectrum above, library spectrum below"></svg><div class="tip" data-f="tip"></div></div>' +
          '<div class="key" style="margin-top:4px"><span><i style="background:var(--accent)"></i>matched</span><span><i style="background:var(--noise)"></i>unmatched</span><span><i style="background:var(--theo)"></i>library ion (rule-based heuristic weight)</span></div></div>' +
      '</section>' +
      '<aside class="pane right" aria-label="Selected species">' +
        '<div class="card mcard"><div class="eyebrow" data-f="cardId"></div><h2 data-f="cardName"></h2><span class="kind" data-f="cardKind"></span><dl class="kv" data-f="cardKv"></dl>' +
          '<div style="margin-top:10px"><div class="eyebrow">Charge states [M−zH]<sup>z−</sup></div><div class="zrow" data-f="zrow"></div></div></div>' +
        '<div class="card" data-f="sparkCard"><div class="eyebrow">Signal across samples</div><div class="plot-wrap"><svg class="plot" data-f="spark" viewBox="0 0 300 120" role="img" aria-label="Percent of total signal across samples"></svg></div><p class="note" data-f="sparkNote"></p></div>' +
        '<div class="card"><div class="eyebrow">Isotope pattern at <span data-f="isoZ"></span></div><div class="plot-wrap"><svg class="plot" data-f="iso" viewBox="0 0 300 130" role="img" aria-label="Theoretical isotope pattern"></svg></div><p class="note" data-f="isoNote"></p></div>' +
        '<div class="card" data-f="compCard"><div class="eyebrow">Composition in <span data-f="compSample"></span></div><div class="stack" data-f="stack" style="margin-top:8px"></div><div class="comp" data-f="comp"></div></div>' +
      '</aside></main>' +
      '<p class="foot" data-f="foot"></p>';
    S.root.querySelectorAll('th button').forEach(function (b) {
      b.onclick = function () { var k = b.dataset.sort; S.sort = { key: k, dir: S.sort.key === k ? -S.sort.dir : (k === 'id' || k === 'mass' ? 1 : -1) }; renderTable(); };
    });
    wire();
  }
  function f(name) { return S.root.querySelector('[data-f="' + name + '"]'); }
  function wire() {
    S.root.querySelectorAll('[data-act]').forEach(function (b) {
      b.onclick = function () {
        var a = b.dataset.act;
        if (a === 'theme') { var n = S.root.getAttribute('data-theme') === 'dark' ? 'light' : 'dark'; S.root.setAttribute('data-theme', n); try { localStorage.setItem('oms-theme', n); } catch (e) {} }
        else if (a === 'demo') send('studio_demo', Date.now());
        else if (a === 'mine') send('studio_mine', Date.now());
        else if (a === 'nav') send('studio_nav', b.dataset.to);
      };
    });
  }

  /* ---------- left pane + notice ---------- */
  function renderFrame() {
    var pl = S.pl, o = pl.oligo;
    f('oname').textContent = o.name;
    f('nodeName').textContent = o.name;
    f('nodeSub').textContent = o.n + '-mer' + (o.conj3 ? ' + ' + o.conj3.split(' ·')[0] : '');
    var msg = {
      demo: '<b>Demo data.</b> Formulas, masses, charge-state m/z, fragment m/z and isotope patterns are computed by the OligoMetProfiler engine for the bundled Givosiran sense strand. Spectra, mass errors, retention times and signal percentages are simulated.',
      results: '<b>Your batch results.</b> Signal is the largest ' + esc(pl.signal_col) + ' per species and sample; percentages are shares of the summed signal of the species shown, with no response-factor correction.',
      library: '<b>Library only.</b> No batch results yet. This is the predicted library for the active oligo. Process samples in Batch Processing to fill in signal, mass errors and MS/MS.'
    }[pl.mode];
    f('notice').innerHTML = msg;
    f('sampleHint').textContent = hasSamples() ? '% parent' : '';
    var box = f('samples'); box.innerHTML = '';
    if (!hasSamples()) box.innerHTML = '<p class="hint">No samples loaded.</p>';
    var par = pl.species.filter(function (s) { return s.cls === 'parent'; })[0];
    pl.samples.forEach(function (s, i) {
      var b = document.createElement('button'); b.type = 'button'; b.className = 'sample' + (i === S.si ? ' active' : '');
      b.setAttribute('role', 'option'); b.setAttribute('aria-selected', i === S.si);
      b.innerHTML = '<span class="fn" title="' + esc(s.label) + '">' + esc(s.label) + '</span><span class="pr">' + (par && par.pct && par.pct[i] != null ? fmt(par.pct[i], 0) + '%' : '—') + '</span>';
      b.onclick = function () { S.si = i; renderAll(); requestDetail(); };
      box.appendChild(b);
    });
    var act = f('actions');
    act.innerHTML = '<button class="btn primary block" data-act="nav" data-to="Batch Processing" type="button">Open Batch Processing</button>' +
      '<button class="btn block" data-act="nav" data-to="Library Generation" type="button">Open Library Generation</button>' +
      (pl.mode === 'demo' ? '<button class="btn block" data-act="mine" type="button">Use my own data</button>' : '<button class="btn block" data-act="demo" type="button">Explore with demo data</button>');
    wire();
    S.root.querySelectorAll('[data-r]').forEach(function (e) { e.style.display = hasSamples() ? '' : 'none'; });
    f('sparkCard').style.display = hasSamples() ? '' : 'none';
    f('compCard').style.display = hasSamples() ? '' : 'none';
    f('foot').textContent = 'Library ion intensities in the mirror plot are the engine’s rule-based weights, not a calibrated intensity model. Percentages exclude any species the batch did not detect.';
  }

  /* ---------- ribbon ---------- */
  function renderRibbon() {
    var o = S.pl.oligo, m = cur(), row = f('ribbon'); row.innerHTML = '';
    if (o.conj5 || (m && m.conj5)) {
      var c5 = m ? m.conj5 : o.conj5; var e5 = document.createElement('span'); e5.className = 'conj c5' + (c5 ? '' : ' ghost'); e5.textContent = c5 || 'no 5′ cap'; row.appendChild(e5);
    }
    for (var i = 0; i < o.n; i++) {
      var pos = i + 1, gone = m && (pos < m.start || pos > m.end), inHl = S.hl && pos >= S.hl[0] && pos <= S.hl[1];
      var nt = document.createElement('div');
      nt.className = 'nt ' + o.sugars[i] + (gone ? ' gone' : '') + (inHl ? ' hl' : '');
      nt.title = 'Position ' + pos + ': ' + o.bases[i] + ', sugar ' + o.sugars[i];
      nt.innerHTML = '<div class="cell">' + esc(o.bases[i]) + '</div><div class="n">' + ([1, 5, 10, 15, 20, o.n].indexOf(pos) >= 0 ? pos : '') + '</div>';
      row.appendChild(nt);
      if (i < o.n - 1) {
        var l = o.linkages[i], lk = document.createElement('div');
        lk.className = 'lk' + (l === 's' ? ' s' : '') + (m && (pos < m.start || pos + 1 > m.end) ? ' gone' : '');
        lk.title = l === 's' ? 'phosphorothioate' : 'phosphodiester'; lk.textContent = l === 's' ? '*' : '·'; row.appendChild(lk);
      }
    }
    var c3 = m ? m.conj3 : o.conj3, e3 = document.createElement('span');
    e3.className = 'conj c3' + (/GalNAc|L96/.test(c3 || '') || (c3 && !/phosphate/.test(c3)) ? '' : ' ghost'); e3.textContent = c3 || 'no 3′ conjugate'; row.appendChild(e3);
    f('mapNote').textContent = m ? m.id + ' · positions ' + m.start + '–' + m.end + ' of ' + o.n : '';
  }

  /* ---------- table ---------- */
  function renderChips() {
    var box = f('chips'); box.innerHTML = '';
    var present = {}; S.pl.species.forEach(function (s) { present[s.cls] = 1; });
    [{ key: 'all', label: 'All' }].concat(CLASSES.filter(function (c) { return present[c.key]; })).forEach(function (c) {
      var b = document.createElement('button'); b.type = 'button'; b.className = 'chip'; b.textContent = c.label; b.setAttribute('aria-pressed', S.filter === c.key);
      b.onclick = function () { S.filter = c.key; renderChips(); renderTable(); }; box.appendChild(b);
    });
  }
  function renderTable() {
    var si = S.si, hs = hasSamples();
    var rows = S.pl.species.filter(function (s) { return S.filter === 'all' || s.cls === S.filter; });
    var k = S.sort.key, d = S.sort.dir; if (!hs && k !== 'id' && k !== 'mass') k = 'id';
    function v(s) { return k === 'id' ? s.id : k === 'mass' ? s.mass : k === 'cov' ? (val(s.cov, si) == null ? -1 : s.cov[si]) : (val(s.pct, si) == null ? -1 : s.pct[si]); }
    rows = rows.slice().sort(function (a, b) { var x = v(a), y = v(b); return (x > y ? 1 : x < y ? -1 : 0) * (k === 'id' || k === 'mass' ? (d === -1 && !hs ? 1 : d) : d); });
    var maxp = Math.max.apply(null, rows.map(function (s) { return hs && val(s.pct, si) != null ? s.pct[si] : 0; }).concat([1]));
    var tb = f('tbody'); tb.innerHTML = '';
    rows.forEach(function (s) {
      var tr = document.createElement('tr'); tr.tabIndex = 0; tr.className = s.id === S.sel ? 'sel' : ''; tr.setAttribute('aria-selected', s.id === S.sel);
      var p = hs ? val(s.pct, si) : null, ppm = hs ? val(s.ppm, si) : null, cv = hs ? val(s.cov, si) : null, col = clsOf(s).color;
      var html = '<td class="id">' + esc(s.id) + '</td><td class="name">' + esc(s.name) + '</td><td class="num mono">' + fmt(s.mass) + '</td>';
      if (hs) html += '<td class="num mono ' + (ppm == null ? 'dim' : Math.abs(ppm) > 5 ? 'bad' : 'ok') + '">' + (ppm == null ? '—' : (ppm >= 0 ? '+' : '') + fmt(ppm, 1)) + '</td>' +
        '<td class="num mono ' + (cv == null ? 'dim' : '') + '">' + (cv == null ? '—' : Math.round(100 * cv) + '%') + '</td>' +
        '<td class="num mono ' + (p == null ? 'dim' : '') + '">' + (p == null ? '—' : '<span class="bar" style="--c:' + col + ';width:' + Math.max(1, 46 * p / maxp) + 'px"></span>' + fmt(p, p < 1 ? 2 : 1)) + '</td>';
      tr.innerHTML = html;
      var pick = function () { S.sel = s.id; S.hl = null; renderAll(); requestDetail(); };
      tr.onclick = pick; tr.onkeydown = function (e) { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); } };
      tb.appendChild(tr);
    });
    f('tblCount').textContent = rows.length + ' of ' + S.pl.species.length + ' species';
  }

  /* ---------- detail ---------- */
  function dkey() { return S.sel + '|' + sampleId(); }
  function requestDetail() {
    if (!S.sel) return;
    if (S.detail[dkey()]) { renderDetail(); return; }
    renderDetail(); send('studio_pick', { id: S.sel, sample: sampleId(), nonce: Date.now() });
  }
  function curDetail() { return S.detail[dkey()] || null; }

  /* ---------- mirror ---------- */
  var MW = 760, ML = 34, MR = 10, MID = 154, TOPY = 14, BOTY = 296;
  function parseLbl(l) { var first = String(l || '').split(' / ')[0], m = LRE.exec(first); return m ? { ion: m[1], len: +m[2], z: +m[3], text: m[1] + m[2] + sup(m[3]) } : null; }
  function niceStep(span) { var raw = span / 8, c = [10, 20, 25, 50, 100, 200, 250, 500, 1000]; for (var i = 0; i < c.length; i++) if (c[i] >= raw) return c[i]; return 1000; }
  function renderMirror() {
    var svg = f('mirror'); svg.innerHTML = ''; var m = cur(), d = curDetail(), tip = f('tip');
    f('mirrorTitle').textContent = m ? m.id + ' · ' + m.name + (hasSamples() ? ' · ' + sampleId() : '') : '';
    if (!d) { svgText(svg, MW / 2, 160, 'Loading…', { 'text-anchor': 'middle' }); f('mirrorStats').innerHTML = ''; return; }
    var acq = (d.acquired || []).map(function (a) { var p = parseLbl(a[2]); return { mz: a[0], int: a[1], lab: a[2], ppm: a[3], p: p, ps: a[2] === 'PS', matched: !!a[2] }; })
      .filter(function (a) { return !a.matched || a.ps || (a.p && S.ions[a.p.ion] !== false); });
    var theo = (d.theoretical || []).map(function (t) { var p = parseLbl(t[2]); return { mz: t[0], int: t[1], lab: t[2], matched: t[3] === 1, p: p }; })
      .filter(function (t) { return t.p && S.ions[t.p.ion] !== false; });
    var all = acq.map(function (a) { return a.mz; }).concat(theo.map(function (t) { return t.mz; }));
    var lo = Math.min.apply(null, all.concat([1e9])), hi = Math.max.apply(null, all.concat([0]));
    if (!all.length) { svgText(svg, MW / 2, 160, 'No ions to show', { 'text-anchor': 'middle' }); return; }
    var step = niceStep(hi - lo), xmin = Math.max(0, Math.floor((lo - step * .3) / step) * step), xmax = Math.ceil((hi + step * .3) / step) * step;
    var X = function (mz) { return ML + (mz - xmin) / (xmax - xmin) * (MW - ML - MR); };
    var hTop = MID - TOPY - 4, hBot = BOTY - MID - 4;
    for (var t = Math.ceil(xmin / step) * step; t <= xmax; t += step) {
      var x = X(t); svg.appendChild(sv('line', { x1: x, x2: x, y1: TOPY, y2: BOTY, 'class': 'grid' })); svgText(svg, x, BOTY + 13, t, { 'text-anchor': 'middle' });
    }
    svg.appendChild(sv('line', { x1: ML, x2: MW - MR, y1: MID, y2: MID, 'class': 'ax' })); svg.appendChild(sv('line', { x1: ML, x2: MW - MR, y1: BOTY, y2: BOTY, 'class': 'ax' }));
    svgText(svg, ML - 5, TOPY + 8, '100%', { 'text-anchor': 'end' }); svgText(svg, ML - 5, MID - 2, '0', { 'text-anchor': 'end' }); svgText(svg, ML - 5, MID + 12, 'lib.', { 'text-anchor': 'end' });
    svgText(svg, MW - MR, BOTY + 26, 'm/z', { 'text-anchor': 'end' });
    theo.forEach(function (t2) {
      var x2 = X(t2.mz); svg.appendChild(sv('line', { x1: x2, x2: x2, y1: MID + 2, y2: MID + 2 + Math.max(2, hBot * .9 * t2.int / 100), stroke: t2.matched ? 'var(--accent)' : 'var(--theo)', 'stroke-width': t2.matched ? 1.6 : 1, opacity: t2.matched ? 1 : .55 }));
    });
    var labs = [], sites = {}, ppms = [];
    if (!d.has_spectrum) svgText(svg, (ML + MW - MR) / 2, TOPY + hTop / 2, 'No MS/MS spectrum captured for this species in this sample', { 'text-anchor': 'middle' });
    acq.slice().sort(function (a, b) { return a.matched - b.matched || a.int - b.int; }).forEach(function (a) {
      var x3 = X(a.mz), h = Math.max(1.5, a.int / 100 * hTop * .92), g = sv('g', { 'class': a.matched ? 'pk' : '' });
      g.appendChild(sv('line', { x1: x3, x2: x3, y1: MID - 1, y2: MID - 1 - h, stroke: a.matched ? 'var(--accent)' : 'var(--noise)', 'stroke-width': a.matched ? 1.6 : 1 }));
      if (a.matched) {
        g.appendChild(sv('line', { x1: x3, x2: x3, y1: MID - 1, y2: MID - 1 - Math.max(h, 12), stroke: 'transparent', 'stroke-width': 7 }));
        a.x = x3; a.h = h; a.text = a.ps ? 'PS' : a.p ? a.p.text : a.lab; labs.push(a);
        if (a.ppm != null) ppms.push(Math.abs(a.ppm));
        if (a.p && m) { var n = m.n, st = FIVE[a.p.ion] ? a.p.len : n - a.p.len; sites[st] = 1; }
        g.addEventListener('mouseenter', function () { showTip(a, x3, MID - 1 - h, svg, tip); setHl(a); });
        g.addEventListener('mouseleave', function () { tip.classList.remove('on'); setHl(null); });
      }
      svg.appendChild(g);
    });
    var placed = [];
    labs.sort(function (a, b) { return b.int - a.int; }).slice(0, 14).forEach(function (a) {
      var w = a.text.length * 6.2, y = MID - 1 - a.h - 4;
      if (placed.some(function (p) { return Math.abs(p.x - a.x) < (p.w + w) / 2 + 2 && Math.abs(p.y - y) < 11; })) return;
      placed.push({ x: a.x, y: y, w: w }); svgText(svg, a.x, y, a.text, { 'text-anchor': 'middle', 'class': 'lab' });
    });
    var nSites = m ? m.n - 1 : 0, cn = Object.keys(sites).length; ppms.sort(function (a, b) { return a - b; });
    f('mirrorStats').innerHTML = d.has_spectrum ? '<span>Coverage <b>' + (nSites ? Math.round(100 * cn / nSites) : 0) + '%</b> (' + cn + '/' + nSites + ' bonds)</span><span>Matched <b>' +
      labs.filter(function (a) { return !a.ps; }).length + '</b> ions</span><span>Median |Δ| <b>' + fmt(ppms.length ? ppms[Math.floor(ppms.length / 2)] : 0, 2) + '</b> ppm</span>' +
      (d.score && d.score.total_score != null ? '<span>Score <b>' + fmt(d.score.total_score, 0) + '</b>' + (d.score.confident ? ' · confident' : '') + '</span>' : '') : '<span class="dim">Library ions only</span>';
  }
  function showTip(a, x, y, svg, tip) {
    var r = svg.getBoundingClientRect();
    tip.innerHTML = esc(a.text) + '<br>m/z ' + fmt(a.mz, 4) + (a.ppm != null ? ' · Δ ' + (a.ppm >= 0 ? '+' : '') + fmt(a.ppm, 1) + ' ppm' : '') + (a.ps ? '<br>phosphorothioate marker' : '');
    tip.style.left = (x / MW * r.width) + 'px'; tip.style.top = (y / 330 * r.height) + 'px'; tip.classList.add('on');
  }
  function setHl(a) {
    var m = cur();
    if (!a || !a.p || !m) S.hl = null; else S.hl = FIVE[a.p.ion] ? [m.start, m.start + a.p.len - 1] : [m.end - a.p.len + 1, m.end];
    renderRibbon();
  }
  function renderToggles() {
    var box = f('toggles'); box.innerHTML = '<span class="eyebrow">Ion series</span>';
    [['a', 'a'], ['a-B', 'a–B'], ['b', 'b'], ['b-B', 'b–B'], ['w', 'w'], ['y', 'y']].forEach(function (p) {
      var l = document.createElement('label'); l.className = 'tg'; var i = document.createElement('input'); i.type = 'checkbox'; i.checked = S.ions[p[0]];
      i.onchange = function () { S.ions[p[0]] = i.checked; renderMirror(); }; l.appendChild(i); l.appendChild(document.createTextNode(p[1])); box.appendChild(l);
    });
  }

  /* ---------- right pane ---------- */
  function renderCard() {
    var m = cur(), si = S.si, hs = hasSamples(), d = curDetail(); if (!m) return;
    f('cardId').textContent = m.id + ' · ' + clsOf(m).label; f('cardName').textContent = m.name; f('cardKind').textContent = m.k_oxid > 0 ? 'Backbone desulfurisation' : (KIND_TEXT[m.kind] || m.kind);
    var z0 = hs ? val(m.z0, si) : null; if (z0 == null && d) z0 = d.z0;
    var rows = [['Formula', m.formula], ['Monoisotopic mass', fmt(m.mass, 4) + ' Da'], ['Length', m.n + ' nt']];
    if (hs) {
      var ppm = val(m.ppm, si), rt = val(m.rt, si), p = val(m.pct, si);
      rows.push(['Observed' + (z0 ? ', z=' + z0 + '−' : ''), fmt(val(m.obs, si), 4)], ['Δ vs theory', ppm == null ? '—' : (ppm >= 0 ? '+' : '') + fmt(ppm, 1) + ' ppm'],
        ['Retention time', rt == null ? '—' : fmt(rt, 2) + ' min'], ['% of signal', p == null ? 'not detected' : fmt(p, 2) + '%']);
    }
    f('cardKv').innerHTML = rows.map(function (r) { return '<dt>' + esc(r[0]) + '</dt><dd>' + esc(r[1]) + '</dd>'; }).join('');
    f('zrow').innerHTML = d ? d.env.map(function (e) { return '<span class="zc' + (e[0] === d.z0 ? ' best' : '') + '">' + e[0] + '− ' + fmt(e[1], 2) + '</span>'; }).join('') : '<span class="dim">Loading…</span>';
  }
  function renderSpark() {
    var m = cur(); if (!m || !hasSamples()) return; var svg = f('spark'); svg.innerHTML = '';
    var W = 300, H = 120, L = 30, R = 12, T = 12, B = 22, ss = S.pl.samples, vals = ss.map(function (s, i) { return val(m.pct, i); });
    var timed = ss.every(function (s) { return s.time != null; }), t0 = timed ? Math.min.apply(null, ss.map(function (s) { return s.time; })) : 0, t1 = timed ? Math.max.apply(null, ss.map(function (s) { return s.time; })) : ss.length - 1;
    var xs = function (i) { var t = timed ? ss[i].time : i; return L + (t1 === t0 ? .5 : (t - t0) / (t1 - t0)) * (W - L - R); };
    var vmax = Math.max.apply(null, vals.map(function (v) { return v == null ? 0 : v; }).concat([1])), ymax = vmax <= 10 ? Math.ceil(vmax / 2) * 2 : Math.ceil(vmax / 10) * 10;
    var Y = function (v) { return T + (1 - v / ymax) * (H - T - B); }, col = clsOf(m).color;
    [0, ymax / 2, ymax].forEach(function (v) { svg.appendChild(sv('line', { x1: L, x2: W - R, y1: Y(v), y2: Y(v), 'class': 'grid' })); svgText(svg, L - 4, Y(v) + 3, v, { 'text-anchor': 'end' }); });
    var step = Math.max(1, Math.ceil(ss.length / 6));
    ss.forEach(function (s, i) { if (i % step === 0 || i === ss.length - 1) svgText(svg, xs(i), H - 7, timed ? s.time : i + 1, { 'text-anchor': 'middle' }); });
    var pts = []; vals.forEach(function (v, i) { if (v != null) pts.push([xs(i), Y(v), i, v]); });
    if (pts.length > 1) {
      svg.appendChild(sv('path', { d: 'M' + pts[0][0] + ',' + Y(0) + ' ' + pts.map(function (p) { return 'L' + p[0] + ',' + p[1]; }).join(' ') + ' L' + pts[pts.length - 1][0] + ',' + Y(0) + ' Z', fill: col, opacity: .14 }));
      svg.appendChild(sv('path', { d: 'M' + pts.map(function (p) { return p[0] + ',' + p[1]; }).join(' L'), fill: 'none', stroke: col, 'stroke-width': 2 }));
    }
    pts.forEach(function (p) {
      var sel = p[2] === S.si, c = sv('circle', { cx: p[0], cy: p[1], r: sel ? 5 : 3, fill: sel ? col : 'var(--panel)', stroke: col, 'stroke-width': 2, style: 'cursor:pointer' });
      var tt = sv('title'); tt.textContent = ss[p[2]].label + ': ' + fmt(p[3], 2) + '%'; c.appendChild(tt); c.onclick = function () { S.si = p[2]; renderAll(); requestDetail(); }; svg.appendChild(c);
      if (sel) svgText(svg, Math.min(Math.max(p[0], 22), W - 22), p[1] - 9, fmt(p[3], 1) + '%', { 'text-anchor': 'middle', 'class': 'lab' });
    });
    f('sparkNote').textContent = (timed ? 'Share of total signal by timepoint.' : 'Share of total signal by sample, in sample order.');
  }
  function renderIso() {
    var d = curDetail(), svg = f('iso'); svg.innerHTML = ''; if (!d) return;
    var W = 300, H = 130, L = 26, R = 8, T = 16, B = 26, n = d.iso.length; if (!n) { f('isoNote').textContent = 'No isotope pattern available.'; return; }
    var bw = (W - L - R) / n, Y = function (a) { return T + (1 - a) * (H - T - B); };
    [0, .5, 1].forEach(function (v) { svg.appendChild(sv('line', { x1: L, x2: W - R, y1: Y(v), y2: Y(v), 'class': 'grid' })); svgText(svg, L - 4, Y(v) + 3, Math.round(v * 100), { 'text-anchor': 'end' }); });
    var apex = d.iso.reduce(function (a, b) { return b[2] > a[2] ? b : a; });
    d.iso.forEach(function (p, i) {
      var x = L + i * bw + bw * .18, w = bw * .64, r = sv('rect', { x: x, y: Y(p[2]), width: w, height: Y(0) - Y(p[2]), fill: i === 0 ? 'var(--warn)' : 'var(--accent)', rx: 1.5 });
      var tt = sv('title'); tt.textContent = 'M+' + p[0] + '  m/z ' + fmt(p[1], 4) + '  ' + fmt(p[2] * 100, 1) + '%'; r.appendChild(tt); svg.appendChild(r);
      svgText(svg, x + w / 2, H - 12, '+' + p[0], { 'text-anchor': 'middle' });
    });
    svgText(svg, L + bw * .5, Y(d.iso[0][2]) - 4, Math.round(d.iso[0][2] * 100) + '%', { 'text-anchor': 'middle', 'class': 'lab' });
    f('isoZ').textContent = 'z = ' + d.z0 + '−';
    f('isoNote').innerHTML = 'The monoisotopic peak (M+0, amber) is <b>' + Math.round(d.iso[0][2] * 100) + '%</b> of the most intense peak (M+' + apex[0] + '). Matching on M+0 alone loses signal for strands this size.';
  }
  function renderComp() {
    if (!hasSamples()) return; var si = S.si; f('compSample').textContent = S.pl.samples[si].label;
    var tot = CLASSES.map(function (c) { return { c: c, v: S.pl.species.filter(function (s) { return s.cls === c.key; }).reduce(function (a, s) { return a + (val(s.pct, si) || 0); }, 0) }; }).filter(function (t) { return t.v > 0; });
    f('stack').innerHTML = tot.map(function (t) { return '<span style="--c:' + t.c.color + ';width:' + t.v + '%" title="' + esc(t.c.label) + ' ' + fmt(t.v, 1) + '%"></span>'; }).join('');
    f('comp').innerHTML = tot.map(function (t) { return '<div><i style="--c:' + t.c.color + '"></i><span>' + esc(t.c.label) + '</span><span>' + fmt(t.v, 1) + '%</span></div>'; }).join('');
  }
  function renderDetail() { renderMirror(); renderCard(); renderIso(); }
  function renderAll() { renderFrame(); renderRibbon(); renderTable(); renderSpark(); renderComp(); renderDetail(); }

  /* ---------- messages ---------- */
  function P(m) { return typeof m === 'string' ? JSON.parse(m) : m; }
  function onLoad(json) {
    var pl = P(json); var prev = S.sel;
    S.pl = pl; S.detail = {};
    var ids = pl.species.map(function (s) { return s.id; });
    if (!prev || ids.indexOf(prev) < 0) { var par = pl.species.filter(function (s) { return s.cls === 'parent'; })[0]; S.sel = par ? par.id : ids[0]; }
    if (S.si >= pl.samples.length) S.si = Math.max(0, Math.floor(pl.samples.length / 2));
    if (S.view !== 'studio') { S.si = Math.floor(pl.samples.length / 2); renderShell(); renderToggles(); }
    renderChips(); renderAll(); requestDetail();
  }
  function onEmpty(json) {
    var o = P(json); S.pl = null; S.sel = null;
    renderEmpty(o.has_library ? { title: 'Studio could not be built from the current results', body: 'Check the R console for the Studio message, or try the demo data.' }
      : { title: 'No oligo defined yet', body: 'Studio shows the sequence map, predicted species, MS/MS and signal for the active oligo. Generate a library first, or explore with the demo data.' });
  }
  function onDetail(json) {
    var d = P(json); S.detail[d.id + '|' + (d.sample || '')] = d; if (S.pl && d.id === S.sel && (d.sample || '') === sampleId()) renderDetail();
  }
  var registered = false;
  function register() {
    if (registered) return;
    S.root = document.getElementById('oms-studio'); if (!S.root) return;
    try { var t = localStorage.getItem('oms-theme'); if (t) S.root.setAttribute('data-theme', t); } catch (e) {}
    Shiny.addCustomMessageHandler('oms-studio-load', onLoad);
    Shiny.addCustomMessageHandler('oms-studio-detail', onDetail);
    Shiny.addCustomMessageHandler('oms-studio-empty', onEmpty);
    registered = true;
  }
  if (window.Shiny && Shiny.addCustomMessageHandler) { document.addEventListener('DOMContentLoaded', register); if (document.readyState !== 'loading') register(); }
  else window.addEventListener('load', register);
})();
