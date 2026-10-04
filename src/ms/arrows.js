'use strict';
/**
 * Dibujo de mecanismos con flechas curvas ("electron pushing") sobre las estructuras SVG de RDKit.
 *
 * Convenciones (IUPAC Gold Book, «curved arrow notation»; McLafferty & Tureček 1993, cap. 4; Gross 2017, §6.1):
 *   - Flecha de media cabeza (anzuelo, ↷): desplazamiento de UN electrón (escisión homolítica, α-escisión,
 *     transferencia de H radicalaria, McLafferty, retro-Diels–Alder).
 *   - Flecha completa (⇒): desplazamiento de un PAR de electrones (escisión heterolítica / inductiva).
 *   - Sitio de ionización: carga ⊕ y electrón desapareado • sobre el par libre del heteroátomo (ionización n),
 *     sobre el enlace π (ionización π) o sobre el enlace σ (ionización σ).
 *
 * Especificación de puntos (índices del grafo de la especie precursora):
 *   ['a', i]        átomo i
 *   ['m', i, j]     punto medio entre los átomos i y j (enlace existente o enlace que se forma)
 *   ['H', i, k]     hidrógeno explícito sobre el átomo i, dibujado en dirección al átomo k (transferencias de H)
 *   ['mH', i, k, j] punto medio entre ese H (del átomo i, orientado a k) y el átomo j
 *   ['mHi', i, k]   punto medio del enlace i–H
 * Flecha: { t: 'fish' | 'full', from: punto, to: punto }
 * Sitio:  { kind: 'n' | 'pi' | 'sigma', a: i, b?: j, label?: '+•' | '+' }
 */

const COL = { fish: '#c2410c', full: '#1d4ed8', site: '#b91c1c', note: '#6b7280', h: '#374151' };

/** Coordenadas 2D (unidades RDKit) desde un molblock V2000/V3000. */
function molblockCoords(mb) {
  const lines = mb.split('\n');
  const out = [];
  const v3 = lines.findIndex(l => l.includes('M  V30 BEGIN ATOM'));
  if (v3 >= 0) {
    for (let i = v3 + 1; i < lines.length && !lines[i].includes('END ATOM'); i++) {
      const p = lines[i].trim().split(/\s+/); // M V30 idx El x y z ...
      out.push({ x: +p[4], y: +p[5] });
    }
    return out;
  }
  const n = parseInt(lines[3].slice(0, 3), 10);
  for (let i = 0; i < n; i++) {
    const l = lines[4 + i];
    out.push({ x: parseFloat(l.slice(0, 10)), y: parseFloat(l.slice(10, 20)) });
  }
  return out;
}

/** Segmentos de enlace del SVG de RDKit: [{bond, a, b, x1, y1, x2, y2}]. */
function svgBondSegments(svg) {
  const re = /<path class='bond-(\d+) atom-(\d+) atom-(\d+)' d='M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)'/g;
  const segs = []; let m;
  while ((m = re.exec(svg))) segs.push({ bond: +m[1], a: +m[2], b: +m[3], x1: +m[4], y1: +m[5], x2: +m[6], y2: +m[7] });
  return segs;
}

/**
 * Ajuste por mínimos cuadrados de la transformación coordenadas RDKit → píxeles SVG
 * (px = s·x + tx; py = −s·y + ty), usando los extremos de los enlaces entre carbonos sin rótulo.
 */
function fitTransform(coords, segs, labeled) {
  const byBond = new Map();
  for (const s of segs) {
    if (!byBond.has(s.bond)) byBond.set(s.bond, []);
    byBond.get(s.bond).push(s);
  }
  const pairs = [];
  for (const [, ss] of byBond) {
    const { a, b } = ss[0];
    if (!coords[a] || !coords[b]) continue;
    // extremos de la línea principal: el primer M y el último L del primer trazo continuo
    const first = ss[0]; let last = ss[0];
    for (const s of ss.slice(1)) if (Math.abs(s.x1 - last.x2) < 0.6 && Math.abs(s.y1 - last.y2) < 0.6) last = s;
    if (!labeled[a]) pairs.push([coords[a], { x: first.x1, y: first.y1 }]);
    if (!labeled[b]) pairs.push([coords[b], { x: last.x2, y: last.y2 }]);
  }
  if (pairs.length < 2) return null;
  // Resolver s, tx, ty con flip vertical
  const solve = (flip) => {
    const n = pairs.length; let sx = 0, sy = 0, sX = 0, sY = 0, sxx = 0, sxX = 0;
    for (const [c, p] of pairs) {
      const x = c.x, y = flip ? -c.y : c.y;
      sx += x; sy += y; sX += p.x; sY += p.y; sxx += x * x + y * y; sxX += x * p.x + y * p.y;
    }
    const den = sxx - (sx * sx + sy * sy) / n;
    if (Math.abs(den) < 1e-9) return null;
    const s = (sxX - (sx * sX + sy * sY) / n) / den;
    const tx = (sX - s * sx) / n, ty = (sY - s * sy) / n;
    let err = 0;
    for (const [c, p] of pairs) { const y = flip ? -c.y : c.y; err += Math.hypot(s * c.x + tx - p.x, s * y + ty - p.y); }
    return { s, tx, ty, flip, err: err / n };
  };
  const cands = [solve(true), solve(false)].filter(Boolean).sort((u, v) => u.err - v.err);
  const best = cands[0];
  if (!best || best.err > 4 || best.s <= 0) return null;
  return (c) => ({ x: best.s * c.x + best.tx, y: best.s * (best.flip ? -c.y : c.y) + best.ty, scale: best.s });
}

const fmt = (v) => (+v).toFixed(1);

/** Punta de flecha (completa o de media cabeza) en el extremo `p`, con dirección de llegada (dx, dy). */
function arrowHead(p, dx, dy, type, color) {
  const L = Math.hypot(dx, dy) || 1; const ux = dx / L, uy = dy / L;
  const len = 9.0, wid = 4.6;
  const bx = p.x - ux * len, by = p.y - uy * len;
  const nx = -uy, ny = ux;
  if (type === 'full') {
    return `<path d='M ${fmt(p.x)},${fmt(p.y)} L ${fmt(bx + nx * wid)},${fmt(by + ny * wid)} L ${fmt(bx - nx * wid)},${fmt(by - ny * wid)} Z' fill='${color}' stroke='none'/>`;
  }
  // media cabeza (un solo barbo) — anzuelo
  return `<path d='M ${fmt(p.x)},${fmt(p.y)} L ${fmt(bx + nx * wid * 1.2)},${fmt(by + ny * wid * 1.2)} L ${fmt(bx)},${fmt(by)} Z' fill='${color}' stroke='none'/>`;
}

/**
 * Superpone flechas, sitio de ionización (⊕ y •), hidrógenos explícitos y rótulos de átomo sobre el SVG.
 * @param svg     SVG de RDKit de la especie precursora (ya dibujado con sus coordenadas)
 * @param P       función índice-local → {x, y} en píxeles
 * @param ctx     { nbrs: índice-local → [índices-locales], labels: índice-local → texto, centroid }
 */
function overlay(svg, P, ctx, mech) {
  const parts = [];
  const bondLen = ctx.bondLen || 30;
  const hPos = new Map();
  const outward = (i) => {
    const p = P(i); const nb = ctx.nbrs[i] || [];
    let vx = 0, vy = 0;
    for (const j of nb) { const q = P(j); vx += p.x - q.x; vy += p.y - q.y; }
    if (!nb.length) { vx = p.x - ctx.centroid.x; vy = p.y - ctx.centroid.y; }
    const L = Math.hypot(vx, vy);
    if (L < 1e-6) return { x: 0, y: -1 };
    return { x: vx / L, y: vy / L };
  };
  const Hof = (i, k) => {
    const key = i + ':' + k;
    if (hPos.has(key)) return hPos.get(key);
    const p = P(i); let d;
    if (ctx.hFixed && ctx.hFixed.atom === i) { const h = { x: ctx.hFixed.p.x, y: ctx.hFixed.p.y, from: p }; hPos.set(key, h); return h; }
    if (k != null && P(k)) { const q = P(k); const L = Math.hypot(q.x - p.x, q.y - p.y) || 1; d = { x: (q.x - p.x) / L, y: (q.y - p.y) / L }; }
    else d = outward(i);
    // mezcla con la dirección externa para no superponer enlaces
    const o = outward(i); let vx = d.x * 0.65 + o.x * 0.35, vy = d.y * 0.65 + o.y * 0.35; const L = Math.hypot(vx, vy) || 1; vx /= L; vy /= L;
    const h = { x: p.x + vx * bondLen * 0.62, y: p.y + vy * bondLen * 0.62, from: p };
    hPos.set(key, h);
    return h;
  };
  const pt = (s) => {
    if (!s) return null;
    try {
      switch (s[0]) {
        case 'a': return P(s[1]);
        case 'm': { const a = P(s[1]), b = P(s[2]); return a && b ? { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 } : null; }
        case 'H': return Hof(s[1], s[2]);
        case 'mH': { const h = Hof(s[1], s[2]); const b = P(s[3]); return b ? { x: (h.x + b.x) / 2, y: (h.y + b.y) / 2 } : null; }
        case 'mHi': { const h = Hof(s[1], s[2]); return { x: (h.x + h.from.x) / 2, y: (h.y + h.from.y) / 2 }; }
        default: return null;
      }
    } catch (_) { return null; }
  };
  // hidrógenos explícitos que intervienen
  for (const a of mech.arrows || []) for (const s of [a.from, a.to]) if (s && (s[0] === 'H' || s[0] === 'mH' || s[0] === 'mHi')) Hof(s[1], s[2]);
  for (const h of hPos.values()) {
    const ux = h.x - h.from.x, uy = h.y - h.from.y; const L = Math.hypot(ux, uy) || 1;
    parts.push(`<line x1='${fmt(h.from.x + ux / L * 4)}' y1='${fmt(h.from.y + uy / L * 4)}' x2='${fmt(h.x - ux / L * 6)}' y2='${fmt(h.y - uy / L * 6)}' stroke='${COL.h}' stroke-width='1.6'/>`);
    parts.push(`<text x='${fmt(h.x)}' y='${fmt(h.y + 4.5)}' font-size='13' font-family='sans-serif' text-anchor='middle' fill='${COL.h}'>H</text>`);
  }
  // sitio de ionización: ⊕ y electrón desapareado •
  const site = mech.site;
  if (site && P(site.a)) {
    let c; let d;
    if (site.b != null && P(site.b)) {
      const a = P(site.a), b = P(site.b); c = { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 };
      let nx = -(b.y - a.y), ny = b.x - a.x; const L = Math.hypot(nx, ny) || 1; nx /= L; ny /= L;
      if ((c.x + nx - ctx.centroid.x) * nx + (c.y + ny - ctx.centroid.y) * ny < 0) { nx = -nx; ny = -ny; }
      d = { x: nx, y: ny };
    } else { c = P(site.a); d = outward(site.a); }
    const off = site.b != null ? 15 : 19;
    const q = { x: c.x + d.x * off, y: c.y + d.y * off };
    const tx = -d.y, ty = d.x;
    const label = site.label || '+•';
    if (label.includes('+')) {
      const cp = { x: q.x - tx * 6, y: q.y - ty * 6 };
      parts.push(`<circle cx='${fmt(cp.x)}' cy='${fmt(cp.y)}' r='5.6' fill='white' stroke='${COL.site}' stroke-width='1.3'/>`);
      parts.push(`<path d='M ${fmt(cp.x - 3.2)},${fmt(cp.y)} L ${fmt(cp.x + 3.2)},${fmt(cp.y)} M ${fmt(cp.x)},${fmt(cp.y - 3.2)} L ${fmt(cp.x)},${fmt(cp.y + 3.2)}' stroke='${COL.site}' stroke-width='1.4'/>`);
    }
    if (label.includes('•')) {
      const dp = { x: q.x + tx * 6, y: q.y + ty * 6 };
      parts.push(`<circle cx='${fmt(dp.x)}' cy='${fmt(dp.y)}' r='2.7' fill='${COL.site}'/>`);
    }
    if (site.kind === 'sigma') parts.push(`<text x='${fmt(q.x)}' y='${fmt(q.y - 9)}' font-size='9' font-family='sans-serif' text-anchor='middle' fill='${COL.site}'>&#963;</text>`);
  }
  // radicales y cargas adicionales (p. ej., ion distónico)
  for (const r of mech.radicals || []) {
    if (!P(r)) continue; const p = P(r); const d = outward(r);
    parts.push(`<circle cx='${fmt(p.x + d.x * 12)}' cy='${fmt(p.y + d.y * 12)}' r='2.7' fill='${COL.site}'/>`);
  }
  // flechas curvas
  const used = [];
  for (const a of mech.arrows || []) {
    let s = pt(a.from), e = pt(a.to);
    if (!s || !e) continue;
    const dx = e.x - s.x, dy = e.y - s.y; const L = Math.hypot(dx, dy);
    if (L < 2) continue;
    // acortar el destino si es un átomo (no tapar el rótulo)
    if (a.to[0] === 'a' || a.to[0] === 'H') { e = { x: e.x - dx / L * 7, y: e.y - dy / L * 7 }; }
    if (a.from[0] === 'a') { s = { x: s.x + dx / L * 6, y: s.y + dy / L * 6 }; }
    const mx = (s.x + e.x) / 2, my = (s.y + e.y) / 2;
    let nx = -dy / L, ny = dx / L;
    //const h = Math.max(13, Math.min(34, L * 0.5));

 // 1 e⁻ (homolítica): arco más pronunciado; 2 e⁻: curvatura moderada
    const fish = a.t !== 'full';
    const h = fish
      ? Math.max(22, Math.min(44, L * 0.9))   // separación máx. ≈ 11–22 px
      : Math.max(13, Math.min(34, L * 0.5));


    
    // curvar hacia el lado más despejado (lejos de átomos, enlaces y flechas ya dibujadas)
    const clear = (x, y) => {
      let m = Infinity;
      for (const q of ctx.obstacles) m = Math.min(m, Math.hypot(q.x - x, q.y - y));
      for (const q of used) m = Math.min(m, Math.hypot(q.x - x, q.y - y) + 4);
      return m;
    };
    const side = (sg) => {
      let m = Infinity;
      for (const t of [0.3, 0.5, 0.7]) {
        const bx = (1 - t) * (1 - t) * s.x + 2 * t * (1 - t) * (mx + sg * nx * h) + t * t * e.x;
        const by = (1 - t) * (1 - t) * s.y + 2 * t * (1 - t) * (my + sg * ny * h) + t * t * e.y;
        m = Math.min(m, clear(bx, by));
      }
      return m;
    };
    if (a.bend === 1 || a.bend === -1) { if (a.bend === -1) { nx = -nx; ny = -ny; } } else if (side(-1) > side(1) + 0.5) { nx = -nx; ny = -ny; }
    const cx = mx + nx * h, cy = my + ny * h;
    used.push({ x: 0.25 * s.x + 0.5 * cx + 0.25 * e.x, y: 0.25 * s.y + 0.5 * cy + 0.25 * e.y });
    const col = a.t === 'full' ? COL.full : COL.fish;
    parts.push(`<path d='M ${fmt(s.x)},${fmt(s.y)} Q ${fmt(cx)},${fmt(cy)} ${fmt(e.x)},${fmt(e.y)}' fill='none' stroke='${col}' stroke-width='1.4' stroke-linecap='round'/>`);
    parts.push(arrowHead(e, e.x - cx, e.y - cy, a.t, col));
  }
  // rótulos de átomo (índices del dibujo de la molécula original)
  if (ctx.labels) {
    for (const [i, txt] of Object.entries(ctx.labels)) {
      const p = P(+i); if (!p) continue; const d = outward(+i);
      parts.push(`<text x='${fmt(p.x - d.x * 9 + 6)}' y='${fmt(p.y - d.y * 9 + 10)}' font-size='8.5' font-family='sans-serif' fill='${COL.note}'>${txt}</text>`);
    }
  }
  const g = `<g class='mech-overlay'>${parts.join('')}</g>`;
  return svg.replace(/<\/svg>\s*$/, g + '</svg>');
}

/**
 * Genera el SVG del precursor con el mecanismo superpuesto.
 * @param RDKit  módulo RDKit MinimalLib
 * @param sp     Species precursora (átomos del grafo sp.g)
 * @param mech   { site, arrows, radicals, highlightBonds: [[i,j]] }
 */
function mechanismSvg(RDKit, sp, mech, { w = 380, h = 280 } = {}) {
  if (!sp || sp.formulaOverride || !sp.atoms || !sp.atoms.size || !mech) return null;
  const js = sp.toRDKitJSON({ neutralSkeleton: true });
  if (!js) return null;
  const order = [...sp.atoms].sort((a, b) => a - b);
  const loc = new Map(order.map((a, i) => [a, i]));
  let mol = null;
  try {
    const lay = layoutMol(RDKit, js, mech.closeRing ? mech.closeRing.map(i => loc.get(i)) : null, mech.bridgeH !== false);
    if (!lay) return null;
    mol = lay.mol;
    const coords = molblockCoords(mol.get_molblock());
    const molJson = JSON.parse(mol.get_json());
    const m0 = molJson.molecules[0]; const dA = molJson.defaults.atom;
    const labeled = m0.atoms.map((a, i) => (a.z ?? dA.z) !== 6 || (ctxDegree(m0, i) === 0));
    const hb = (mech.highlightBonds || []).map(([i, j]) => [loc.get(i), loc.get(j)]).filter(([i, j]) => i != null && j != null);
    const bondIdx = []; const bondCol = {};
    (m0.bonds || []).forEach((b, k) => {
      if (hb.some(([i, j]) => (b.atoms[0] === i && b.atoms[1] === j) || (b.atoms[0] === j && b.atoms[1] === i))) { bondIdx.push(k); bondCol[k] = [0.98, 0.75, 0.55]; }
    });
    const svg = bondIdx.length
      ? mol.get_svg_with_highlights(JSON.stringify({ width: w, height: h, bonds: bondIdx, atoms: [], highlightBondColors: bondCol, highlightBondWidthMultiplier: 12, padding: 0.12 }))
      : mol.get_svg_with_highlights(JSON.stringify({ width: w, height: h, padding: 0.12 }));
    const T = fitTransform(coords, svgBondSegments(svg), labeled);
    if (!T) return null;
    const px = coords.map(T);
    const P = (gi) => { const k = loc.get(gi); return k == null ? null : px[k]; };
    const nbrs = {};
    for (const gi of order) nbrs[gi] = [];
    for (const b of (m0.bonds || [])) { const [i, j] = b.atoms; nbrs[order[i]].push(order[j]); nbrs[order[j]].push(order[i]); }
    for (const e of sp.extraBonds || []) { if (nbrs[e.a] && nbrs[e.b]) { nbrs[e.a].push(e.b); nbrs[e.b].push(e.a); } }
    const centroid = px.reduce((s, p) => ({ x: s.x + p.x / px.length, y: s.y + p.y / px.length }), { x: 0, y: 0 });
    let bl = 0, nb = 0;
    for (const b of (m0.bonds || [])) { const p = px[b.atoms[0]], q = px[b.atoms[1]]; bl += Math.hypot(p.x - q.x, p.y - q.y); nb++; }
    const labels = {};
    for (const gi of order) { const a = sp.g.atoms[gi]; labels[gi] = String(a.orig ?? gi); }
    const obstacles = [...px];
    for (const b of (m0.bonds || [])) { const p = px[b.atoms[0]], q = px[b.atoms[1]]; for (const t of [0.25, 0.5, 0.75]) obstacles.push({ x: p.x + (q.x - p.x) * t, y: p.y + (q.y - p.y) * t }); }
    const hFixed = lay.hXY && mech.closeRing ? { atom: mech.closeRing[1], p: T(lay.hXY) } : null;
    return overlay(svg, P, { nbrs, labels, centroid, obstacles, bondLen: nb ? bl / nb : 30, hFixed }, mech);
  } catch (_) {
    return null;
  } finally {
    if (mol) try { mol.delete(); } catch (_) { /* noop */ }
  }
}
/**
 * Crea la molécula con coordenadas 2D. Si se indica `close` = [x, d] (aceptor, donador de H), las coordenadas se
 * calculan sobre una copia en la que el H transferido se representa con un pseudoátomo puente x···H···d; así el
 * estado de transición cíclico de seis miembros (McLafferty, eliminaciones 1,4, ciclaciones) aparece plegado como
 * en los libros de texto. Las coordenadas se trasladan después a la molécula original (sin modificar su química).
 * Devuelve { mol, hXY } donde hXY son las coordenadas RDKit del H puente (o null).
 */
function layoutMol(RDKit, js, close, bridgeH = true) {
  const plain = () => {
    const m = RDKit.get_mol(JSON.stringify(js));
    if (!m || !m.is_valid()) { if (m) m.delete(); return null; }
    if (typeof m.set_new_coords === 'function') m.set_new_coords(true);
    return { mol: m, hXY: null };
  };
  if (!close || close.some(i => i == null) || close[0] === close[1]) return plain();
  const js2 = JSON.parse(JSON.stringify(js));
  const M = js2.molecules[0]; const nA = M.atoms.length;
  const compensate = (i) => {
    const at = M.atoms[i]; const h = at.impHs ?? js2.defaults.atom.impHs;
    if (h > 0) at.impHs = h - 1; else at.chg = (at.chg ?? 0) + 1;
  };
  if (bridgeH) {
    M.atoms.push({ z: 0, impHs: 0 });
    M.bonds = [...(M.bonds || []), { bo: 1, atoms: [close[0], nA] }, { bo: 1, atoms: [close[1], nA] }];
  } else {
    M.bonds = [...(M.bonds || []), { bo: 1, atoms: [close[0], close[1]] }];
  }
  compensate(close[0]); compensate(close[1]);
  let ring = null; let base = null;
  try {
    ring = RDKit.get_mol(JSON.stringify(js2));
    if (!ring || !ring.is_valid()) throw new Error('anillo no válido');
    if (typeof ring.set_new_coords === 'function') ring.set_new_coords(true);
    const rc = molblockCoords(ring.get_molblock());
    const b = plain(); if (!b) return null; base = b.mol;
    const mb = base.get_molblock().split('\n');
    if (mb.some(l => l.includes('V3000'))) return { mol: base, hXY: null };
    const n = parseInt(mb[3].slice(0, 3), 10);
    for (let i = 0; i < n; i++) {
      const l = mb[4 + i];
      mb[4 + i] = rc[i].x.toFixed(4).padStart(10, ' ') + rc[i].y.toFixed(4).padStart(10, ' ') + l.slice(20);
    }
    const m = RDKit.get_mol(mb.join('\n'), JSON.stringify({ removeHs: false }));
    if (!m || !m.is_valid()) { if (m) m.delete(); return { mol: base, hXY: null }; }
    base.delete(); base = null;
    return { mol: m, hXY: bridgeH ? rc[nA] : null };
  } catch (_) {
    if (base) return { mol: base, hXY: null };
    return plain();
  } finally {
    if (ring) try { ring.delete(); } catch (_) { /* noop */ }
  }
}
function ctxDegree(m0, i) { return (m0.bonds || []).filter(b => b.atoms.includes(i)).length; }

/** Corchetes y «⁺•» para iones impar-electrónicos dibujados con su esqueleto neutro. */
function bracketIon(svg, label = '+•') {
  if (!svg) return svg;
  const mw = /width='(\d+)px'/.exec(svg); const mh = /height='(\d+)px'/.exec(svg);
  if (!mw || !mh) return svg;
  const W = +mw[1], H = +mh[1];
  const x1 = 6, x2 = W - 16, y1 = 10, y2 = H - 10;
  const s = '#374151';
  const g = `<g class='mech-bracket'><path d='M ${x1 + 6},${y1} L ${x1},${y1} L ${x1},${y2} L ${x1 + 6},${y2} M ${x2 - 6},${y1} L ${x2},${y1} L ${x2},${y2} L ${x2 - 6},${y2}' fill='none' stroke='${s}' stroke-width='1.4'/>` +
    `<circle cx='${x2 + 6}' cy='${y1 + 4}' r='4.6' fill='white' stroke='${COL.site}' stroke-width='1.2'/><path d='M ${x2 + 3.4},${y1 + 4} L ${x2 + 8.6},${y1 + 4} M ${x2 + 6},${y1 + 1.4} L ${x2 + 6},${y1 + 6.6}' stroke='${COL.site}' stroke-width='1.2'/>` +
    (label.includes('•') ? `<circle cx='${x2 + 7}' cy='${y1 + 15}' r='2.5' fill='${COL.site}'/>` : '') + '</g>';
  return svg.replace(/<\/svg>\s*$/, g + '</svg>');
}

module.exports = { mechanismSvg, bracketIon, fitTransform, molblockCoords };
