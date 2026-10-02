/* EI-MS Predictor — cliente */
(function () {
  'use strict';
  const $ = (s, r = document) => r.querySelector(s);
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const sub = (f) => esc(f).replace(/(\d+)/g, '<sub>$1</sub>');
  let last = null; let syncingFromJsme = false;

  /* ---------- tema ---------- */
  const root = document.documentElement;
  try { const t = localStorage.getItem('eims-theme'); if (t) root.dataset.theme = t; } catch (_) { /* sin almacenamiento */ }
  $('#themeBtn').addEventListener('click', () => {
    const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
    root.dataset.theme = dark ? 'light' : 'dark';
    try { localStorage.setItem('eims-theme', root.dataset.theme); } catch (_) { /* noop */ }
  });

  /* ---------- sincronización JSME ↔ SMILES ---------- */
  const input = $('#smiles');
  window.addEventListener('jsme-change', () => {
    if (!window.jsmeApplet) return;
    syncingFromJsme = true; input.value = window.jsmeApplet.smiles() || ''; syncingFromJsme = false;
  });
  let t = null;
  function pushToJsme(s) {
    if (!window.jsmeApplet || !s) return;
    try { window.jsmeApplet.readGenericMolecularInput(s); } catch (_) { /* SMILES parcial */ }
  }
  input.addEventListener('input', () => { if (syncingFromJsme) return; clearTimeout(t); t = setTimeout(() => pushToJsme(input.value.trim()), 500); });
  window.addEventListener('jsme-ready', () => { if (input.value) pushToJsme(input.value); });
  $('#examples').addEventListener('change', (e) => {
    const v = e.target.value; if (!v) return;
    input.value = v; pushToJsme(v); e.target.selectedIndex = 0; run();
  });
  $('#form').addEventListener('submit', (e) => { e.preventDefault(); run(); });

  /* ---------- llamada a la API ---------- */
  async function run() {
    let smiles = input.value.trim();
    if (!smiles && window.jsmeApplet) smiles = window.jsmeApplet.smiles();
    if (!smiles) return setStatus('Dibuje una estructura o escriba un SMILES.', 'err');
    const btn = $('#go'); btn.disabled = true; setStatus('Analizando estructura y generando rutas de fragmentación…', 'busy');
    try {
      const r = await fetch('api/predict', { method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles, threshold: Number($('#threshold').value) }) });
      const data = await r.json();
      if (!r.ok) throw new Error(data.error || `Error ${r.status}`);
      last = data; render(data); setStatus('');
    } catch (e) { setStatus(e.message, 'err'); }
    finally { btn.disabled = false; }
  }
  function setStatus(msg, cls) { const s = $('#status'); s.textContent = msg; s.className = 'status ' + (msg ? cls || '' : ''); }

  /* ---------- renderizado ---------- */
  function render(d) {
    $('#results').classList.remove('hidden');
    $('#parentSvg').innerHTML = d.parentSvg || '';
    const rows = [
      ['SMILES canónico', `<code>${esc(d.smiles)}</code>`],
      ['Fórmula molecular', sub(d.formula)],
      ['M⁺• (nominal)', `m/z ${d.nominalMass}`],
      ['Masa monoisotópica', `${d.monoisotopicMass.toFixed(5)} u <span class="muted">(m/z M⁺• ${(d.monoisotopicMass - 0.000549).toFixed(4)})</span>`],
      ['Insaturaciones (RDB)', d.rdb],
      ['Regla del nitrógeno', esc(d.nitrogenRule)],
      ['Pico base predicho', d.basePeak ? `m/z ${d.basePeak.mz} (${sub(d.basePeak.formula)})` : '—'],
      ['Estabilidad de M⁺•', `factor ${d.molecularIon.factor}${d.molecularIon.classes.length ? ' · ' + esc(d.molecularIon.classes.join(', ')) : ''}`],
      ['Motor', `${esc(d.meta.engine)} · RDKit ${esc(d.meta.rdkit)} · ${d.meta.events} rutas evaluadas · ${d.meta.ms} ms`]
    ];
    $('#summary').innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('#disclaimer').textContent = d.meta.disclaimer;
    drawChart(d.spectrum);
    drawTable(d.spectrum);
    drawExplanations(d.explanations);
    $('#refs').innerHTML = d.references.map(r => `<li>${esc(r.text)}</li>`).join('');
    $('#results').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  function drawChart(spec, hlMz) {
    const W = 960, H = 360, m = { l: 52, r: 18, t: 30, b: 44 };
    const iw = W - m.l - m.r, ih = H - m.t - m.b;
    const maxMz = Math.max(...spec.map(p => p.mz));
    const minMz = Math.min(...spec.map(p => p.mz));
    const x0 = Math.max(0, Math.floor((minMz - 5) / 10) * 10), x1 = Math.ceil((maxMz + 8) / 10) * 10;
    const X = v => m.l + (v - x0) / (x1 - x0) * iw; const Y = v => m.t + ih - v / 100 * ih;
    const step = niceStep((x1 - x0) / 10);
    let g = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Espectro de masas predicho"><g class="grid">`;
    for (let y = 0; y <= 100; y += 20) g += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(y)}" y2="${Y(y)}"/>`;
    g += `</g><g class="axis">`;
    for (let y = 0; y <= 100; y += 20) g += `<text x="${m.l - 8}" y="${Y(y) + 4}" text-anchor="end">${y}</text>`;
    g += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(0)}" y2="${Y(0)}"/>`;
    for (let v = Math.ceil(x0 / step) * step; v <= x1; v += step) g += `<line x1="${X(v)}" x2="${X(v)}" y1="${Y(0)}" y2="${Y(0) + 5}"/><text x="${X(v)}" y="${Y(0) + 18}" text-anchor="middle">${v}</text>`;
    g += `<text class="axis-title" x="${m.l + iw / 2}" y="${H - 6}" text-anchor="middle">m/z</text>`;
    g += `<text class="axis-title" transform="translate(14 ${m.t + ih / 2}) rotate(-90)" text-anchor="middle">Abundancia relativa (%)</text></g>`;
    // barras (sticks)
    for (const p of spec) {
      const iso = p.assignments[0] && p.assignments[0].isotope;
      g += `<line class="stick${iso ? ' iso' : ''}${p.mz === hlMz ? ' hl' : ''}" x1="${X(p.mz)}" x2="${X(p.mz)}" y1="${Y(0)}" y2="${Y(Math.max(p.intensity, 0.4))}"/>`;
    }
    // etiquetas selectivas sin colisión
    const placed = [];
    for (const p of spec.slice().sort((a, b) => b.intensity - a.intensity)) {
      if (p.intensity < 6 || placed.length >= 14) break;
      const x = X(p.mz);
      if (placed.some(px => Math.abs(px - x) < 26)) continue;
      placed.push(x); g += `<text class="lbl" x="${x}" y="${Y(p.intensity) - 6}">${p.mz}</text>`;
    }
    // zonas de interacción
    const hw = Math.max(5, iw / (x1 - x0) / 2);
    spec.forEach((p, i) => { g += `<rect class="hit" data-i="${i}" x="${X(p.mz) - hw}" y="${m.t}" width="${hw * 2}" height="${ih}"/>`; });
    g += '</svg>';
    const el = $('#chart'); el.innerHTML = g;
    const tip = $('#tip');
    el.querySelectorAll('.hit').forEach(r => {
      r.addEventListener('mousemove', (e) => {
        const p = spec[+r.dataset.i];
        tip.innerHTML = `<b>m/z ${p.mz}</b> · ${p.intensity.toFixed(1)} %<br>` + p.assignments.slice(0, 3).map(a =>
          `${sub(a.formula)}${a.isotope ? ` <span class="muted">(isótopo M+${a.offset})</span>` : ''} <span class="muted">${Math.round(a.share * 100)}%</span>`).join('<br>');
        tip.hidden = false; tip.style.left = Math.min(e.clientX + 14, innerWidth - 290) + 'px'; tip.style.top = (e.clientY + 14) + 'px';
      });
      r.addEventListener('mouseleave', () => { tip.hidden = true; });
      r.addEventListener('click', () => {
        const p = spec[+r.dataset.i]; const card = document.querySelector(`.ion[data-mz="${p.mz}"]`);
        if (card) { card.open = true; card.scrollIntoView({ behavior: 'smooth', block: 'center' }); }
      });
    });
  }
  function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw))); const n = raw / p; return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * p; }

  function drawTable(spec) {
    $('#peakTable tbody').innerHTML = spec.map(p => {
      const a = p.assignments[0] || {};
      return `<tr><td>${p.mz}</td><td class="num">${p.intensity.toFixed(1)}</td><td>${p.assignments.slice(0, 3).map(x => sub(x.formula) + (x.isotope ? ` (M+${x.offset})` : '')).join(', ')}</td><td>${a.isotope ? 'isotópico' : 'monoisotópico'}</td></tr>`;
    }).join('');
  }
  document.querySelectorAll('.seg button').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('.seg button').forEach(x => x.classList.toggle('on', x === b));
    $('#chartView').hidden = b.dataset.view !== 'chart'; $('#tableView').hidden = b.dataset.view !== 'table';
  }));

  function drawExplanations(list) {
    $('#explanations').innerHTML = list.map((x, i) => {
      const oe = x.ionType.startsWith('OE');
      const steps = x.pathway.map((s, k) => `
        <div class="step">
          <h4>${k + 1}. ${esc(s.ruleName)}</h4>
          <div class="eq">
            ${k === 0 && s.rule !== 'M' ? '<span class="mono">M⁺•</span><span class="op">→</span>' : k > 0 ? `<span class="mono">m/z ${x.pathway[k - 1].ionMz}</span><span class="op">→</span>` : ''}
            ${s.ionSvg ? `<figure class="mol">${s.ionSvg}<figcaption>${sub(s.ionFormula)}${s.oddElectron ? '⁺•' : '⁺'} · m/z ${s.ionMz}</figcaption></figure>` : `<span class="mono"><b>${sub(s.ionFormula)}${s.oddElectron ? '⁺•' : '⁺'}</b> (m/z ${s.ionMz})</span>`}
            ${s.neutralFormula ? `<span class="op">+</span>${s.neutralSvg ? `<figure class="mol small">${s.neutralSvg}<figcaption>${sub(s.neutralFormula)}</figcaption></figure>` : `<span class="mono">${sub(s.neutralFormula)}</span>`}<span class="muted small">(${esc(s.neutralLabel || 'neutro')})</span>` : ''}
          </div>
          <ol>${s.mechanism.map(t => `<li>${esc(t)}</li>`).join('')}</ol>
          ${s.ionSmiles ? `<div class="small muted">SMILES del ion (esqueleto): <code>${esc(s.ionSmiles)}</code>${s.oddElectron ? ' ⁺•' : ''}</div>` : ''}
        </div>`).join('');
      return `<details class="ion" data-mz="${x.mz}" ${i < 3 ? 'open' : ''}>
        <summary><span class="mz">m/z ${x.mz}</span><span class="bar" title="${x.intensity}%"><i style="width:${Math.min(100, x.intensity)}%"></i></span>
          <span class="title">${sub(x.formula)}${oe ? '⁺•' : '⁺'} — ${esc(x.pathway[x.pathway.length - 1].ruleName)}</span>
          <span class="badge ${oe ? 'oe' : 'ee'}">${esc(x.ionType)}</span>
          <span class="muted small">${x.intensity.toFixed(1)} % · m/z exacto ${x.mzExact}${x.lossFromM ? ` · M−${x.lossFromM}` : ''}</span></summary>
        <div class="ion-body">
          <div><div class="mol">${x.parentHighlightSvg || ''}</div><p class="small muted">RDB = ${x.rdb} (${Number.isInteger(x.rdb) ? 'entero → ion impar-electrónico' : 'semientero → ion par-electrónico'})</p></div>
          <div>${steps}${x.alternatives.length ? `<div class="alt">Rutas competitivas que originan la misma composición: ${x.alternatives.map(a => `${esc(a.ruleName)} (${Math.round(a.share * 100)} %)`).join('; ')}</div>` : ''}</div>
        </div></details>`;
    }).join('');
    document.querySelectorAll('.ion > summary').forEach(s => s.addEventListener('click', () => { const d = s.parentElement; if (!d.open && last) drawChart(last.spectrum, +d.dataset.mz); }));
  }

  /* ---------- exportación ---------- */
  function download(name, text, type) {
    const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name;
    document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }
  document.querySelectorAll('[data-export]').forEach(b => b.addEventListener('click', () => {
    if (!last) return; const kind = b.dataset.export; const base = 'eims_' + last.formula;
    if (kind === 'print') return window.print();
    if (kind === 'csv') {
      const csv = ['mz,intensidad_rel,composicion,isotopico'].concat(last.spectrum.map(p => `${p.mz},${p.intensity},${p.assignments.map(a => a.formula).join('|')},${p.assignments[0].isotope ? 1 : 0}`)).join('\n');
      return download(base + '.csv', csv, 'text/csv');
    }
    if (kind === 'json') {
      const clean = JSON.parse(JSON.stringify(last, (k, v) => (/svg$/i.test(k) ? undefined : v)));
      return download(base + '.json', JSON.stringify(clean, null, 2), 'application/json');
    }
    if (kind === 'jdx') {
      // JCAMP-DX 5.01 para espectros de masas (McDonald & Wilks, Appl. Spectrosc. 1988, 42, 151; Lampen et al., Appl. Spectrosc. 1994, 48, 1545)
      const now = new Date();
      const lines = ['##TITLE=Espectro EI predicho: ' + last.smiles, '##JCAMP-DX=5.01', '##DATA TYPE=MASS SPECTRUM', '##DATA CLASS=PEAK TABLE',
        '##ORIGIN=EI-MS rule-based predictor (prediccion in silico)', '##OWNER=public domain', `##LONGDATE=${now.toISOString().slice(0, 10).replace(/-/g, '/')}`,
        `##MOLFORM=${last.formula.replace(/(\d+)/g, ' $1 ').trim()}`, `##MW=${last.nominalMass}`, '##SPECTROMETER/DATA SYSTEM=in silico', '##.SPECTROMETER TYPE=EI', '##.IONIZATION ENERGY=70 eV',
        '##XUNITS=M/Z', '##YUNITS=RELATIVE ABUNDANCE', `##NPOINTS=${last.spectrum.length}`, '##PEAK TABLE=(XY..XY)'];
      for (let i = 0; i < last.spectrum.length; i += 5) lines.push(last.spectrum.slice(i, i + 5).map(p => `${p.mz},${Math.round(p.intensity * 9.99)}`).join(' '));
      lines.push('##END=');
      return download(base + '.jdx', lines.join('\n'), 'chemical/x-jcamp-dx');
    }
  }));
})();
