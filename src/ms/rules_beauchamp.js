'use strict';
/**
 * Reglas de fragmentación EI incorporadas a partir de:
 *   Beauchamp, P. "Basics of Mass Spectroscopy" (capítulo del Spectroscopy Workbook),
 *   Cal Poly Pomona — secciones "Special patterns of fragmentation from organic functional groups"
 *   y hoja resumen "Common fragmentation patterns in mass spectroscopy" (MS_chapter.pdf).
 * Contrastadas con McLafferty & Tureček (1993) y Gross (2017).
 *
 * Patrones que no estaban en la versión 1.1 del motor:
 *   B1. Migración del doble/triple enlace en el radical-catión antes de fragmentarse
 *       (alquenos §2, alquinos §2) → escisiones alílicas/propargílicas y McLafferty de los isómeros.
 *   B2. Alquinos terminales: [M−1]⁺ por pérdida de H• propargílico (alquinos §1).
 *   B3. Alcanos: eliminación de un alcano R–H (y fragmentos alquénicos resultantes) (alcanos §8, p. 19).
 *   B4. Eliminación de moléculas pequeñas: ROH (éteres), RSH (sulfuros), NH₃/RNH₂ (aminas) (hoja resumen, p. 44).
 *   B7. −(H₂O + C₂H₄) / −(H₂S + C₂H₄) en 1-alcanoles/1-alcanotioles.
 * Ya presentes en v1.1 y por tanto NO duplicados aquí: α-escisión anular de cicloalcanoles/aminas (ringAlpha),
 * pérdida de RCOOH/alqueno en ésteres (McLafferty lado alcoxilo), efecto orto.
 *   B6. Fragmentación alquénica consecutiva de los iones alqueno formados por eliminación
 *       (−H₂O, −H₂S, −HX, −ROH, −RSH, −RH…) y por McLafferty de alquenos ("This can lead to alkene fragmentations").
 */
const { Species } = require('./species');
const { MolGraph } = require('../chem/molgraph');
const R = require('./rules');
const { L, B, ev, split, allAtoms, complement, benzylicAllylic, mclafferty } = R;
const { radicalClass } = require('./stability');

/* ------------------------------------------------------------------ */
/* B1. Migración del enlace múltiple C=C / C≡C                         */
/* ------------------------------------------------------------------ */
const MIGR_W = [1, 0.5, 0.35]; // migración extensa: "can migrate through the skeleton to almost any conceivable position"

function shiftCandidates(g, onlyBond = null, avoid = null) {
  const out = [];
  for (const b of g.bonds) {
    if (onlyBond != null && b.idx !== onlyBond) continue;
    if (b.aromatic || b.inRing || (b.order !== 2 && b.order !== 3)) continue;
    if (!g.isC(b.a) || !g.isC(b.b)) continue;
    for (const [p, q] of [[b.a, b.b], [b.b, b.a]]) {
      // q no debe tener otro enlace múltiple (evita alenos/dienos conjugados)
      if (g.atoms[q].nbrs.some(n => n.bond !== b.idx && g.bonds[n.bond].order > 1)) continue;
      if (g.atoms[p].nbrs.some(n => n.bond !== b.idx && (g.bonds[n.bond].order > 1 || g.bonds[n.bond].aromatic))) continue;
      for (const n of g.atoms[q].nbrs) {
        const r = n.atom;
        if (r === p || r === avoid || n.bond === b.idx) continue;
        const bqr = g.bonds[n.bond];
        if (bqr.order !== 1 || bqr.inRing || !g.isC(r) || g.atoms[r].aromatic || !g.isSp3(r)) continue;
        if (g.atoms[r].h < (b.order === 3 ? 2 : 1)) continue;
        out.push({ bond: b.idx, order: b.order, p, q, r, bqr: n.bond });
      }
    }
  }
  return out;
}

function applyShift(g, c) {
  const g2 = g.clone();
  const dh = c.order === 3 ? 2 : 1;
  g2.bonds[c.bond].order = 1;
  g2.bonds[c.bqr].order = c.order;
  g2.atoms[c.p].h += dh;
  g2.atoms[c.r].h -= dh;
  return g2;
}

/** Eventos alílicos/propargílicos y McLafferty que involucran el enlace múltiple recién formado (q=r). */
function eventsOnNewBond(g2, q, r) {
  const nb = new Set([q, r]);
  const touches = (bi) => {
    const bd = g2.bonds[bi];
    return [bd.a, bd.b].some(e => !nb.has(e) && g2.atoms[e].nbrs.some(n => nb.has(n.atom)));
  };
  const out = [];
  for (const e of benzylicAllylic(g2)) if ((e.rule === 'allylic' || e.rule === 'propargylic') && e.cleaved.length && touches(e.cleaved[0])) out.push(e);
  for (const e of mclafferty(g2)) {
    const t = (e.ion.tags.mclafferty || (e.neutral && e.neutral.tags.mclafferty));
    if (t && (t.kind === 'alqueno' || t.kind === 'alquino') && nb.has(t.x) && nb.has(t.y)) out.push(e);
  }
  return out;
}

function migration(g) {
  const out = [];
  const seenG = new Set();
  const sig = (gx) => gx.bonds.map(b => b.order).join('') + '|' + gx.atoms.map(a => a.h).join('');
  const recurse = (gx, depth, c, path) => {
    const g2 = applyShift(gx, c);
    const s = sig(g2); if (seenG.has(s)) return; seenG.add(s);
    const w = MIGR_W[depth];
    const mult = c.order === 3 ? 'triple' : 'doble';
    const sym = c.order === 3 ? '≡' : '=';
    const hopTxt = [...path, `${L(g, c.p)}${sym}${L(g, c.q)} → ${L(g, c.q)}${sym}${L(g, c.r)}`].join('; ');
    for (const e of eventsOnNewBond(g2, c.q, c.r)) {
      out.push(Object.assign(e, {
        rule: 'migr_' + e.rule,
        ruleName: `Migración del ${mult} enlace + ${e.ruleName.charAt(0).toLowerCase() + e.ruleName.slice(1)}`,
        score: e.score * w,
        steps: [
          `Isomerización del radical-catión antes de fragmentarse: el enlace ${mult} migra por desplazamientos de H (1,2/1,3) a lo largo de la cadena (${hopTxt}). ` +
          'Por ello los isómeros de posición con el mismo esqueleto dan espectros muy parecidos (Beauchamp, alquenos/alquinos §2).',
          ...e.steps
        ],
        refs: [...new Set([...(e.refs || []), 'BEAU', 'MT'])]
      }));
    }
    if (depth < 2) for (const c2 of shiftCandidates(g2, c.bqr, c.q)) recurse(g2, depth + 1, c2, [...path, `${L(g, c.p)}${sym}${L(g, c.q)} → ${L(g, c.q)}${sym}${L(g, c.r)}`]);
  };
  for (const c of shiftCandidates(g)) recurse(g, 1, c, []);
  return out;
}

/* ------------------------------------------------------------------ */
/* B2. Pérdida de H• propargílico ([M−1]⁺)                              */
/* ------------------------------------------------------------------ */
function propargylicH(g) {
  const out = [];
  for (const b of g.bonds) {
    if (b.order !== 3 || !g.isC(b.a) || !g.isC(b.b)) continue;
    const terminal = g.atoms[b.a].h > 0 || g.atoms[b.b].h > 0;
    for (const s of [b.a, b.b]) {
      for (const n of g.atoms[s].nbrs) {
        const c = n.atom;
        if (n.bond === b.idx || !g.isC(c) || !g.isSp3(c) || g.atoms[c].h === 0) continue;
        const ion = new Species(g, allAtoms(g), { dH: new Map([[c, -1]]), chg: new Map([[c, 1]]), tags: { propargyl: true } });
        out.push(ev({
          rule: 'propargylH', ruleName: 'Pérdida de H• propargílico ([M−1]⁺, alquinos)', ion, neutral: null, neutralFormula: { H: 1 }, neutralLabel: 'H•',
          cleaved: [], oe: false, score: terminal ? 0.35 : 0.12,
          steps: [
            `Ionización π del triple enlace ${L(g, b.a)}≡${L(g, b.b)}.`,
            `Ruptura homolítica de un enlace C–H propargílico en ${L(g, c)}: se forma el catión propargilo/alenilo [M−H]⁺ estabilizado por resonancia.`,
            'Los alquinos terminales muestran M⁺• débil o ausente y un [M−1]⁺ apreciable.'
          ],
          refs: ['BEAU', 'MT']
        }));
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* B3. Alcanos: eliminación de R–H → radical-catión alqueno             */
/* ------------------------------------------------------------------ */
const ALKENE_SIZE_W = { 2: 0.3, 3: 0.7, 4: 1.0, 5: 0.8, 6: 0.6 };
function alkaneRH(g) {
  const out = [];
  if (g.atoms.some(a => a.el !== 'C') || g.rings.length || g.bonds.some(b => b.order !== 1) || g.n < 5) return out;
  const seen = new Set();
  for (const b of g.bonds) {
    for (const [ca, cg] of [[b.a, b.b], [b.b, b.a]]) {
      const deg = g.atoms[ca].nbrs.length;
      if (deg < 2) continue;
      const parts = split(g, b.idx, ca); if (!parts) continue;
      const [ionAt, neuAt] = parts;
      const cbs = g.atoms[ca].nbrs.filter(n => n.atom !== cg && g.atoms[n.atom].h > 0).sort((x, y) => g.atoms[x.atom].h - g.atoms[y.atom].h);
      if (!cbs.length) continue;
      const cb = cbs[0];
      const key = [...ionAt].sort((x, y) => x - y).join(',');
      if (seen.has(key)) continue; seen.add(key);
      const nIon = ionAt.size;
      const ion = new Species(g, ionAt, { dH: new Map([[cb.atom, -1]]), bo: new Map([[cb.bond, 2]]), oddElectron: true, depictAsNeutral: true, tags: { alkeneIon: true } });
      const neu = new Species(g, neuAt, { dH: new Map([[cg, 1]]) });
      const wDeg = deg >= 4 ? 2.0 : deg === 3 ? 1.8 : 0.12;
      const wN = neuAt.size === 1 ? 0.15 : neuAt.size <= 3 ? 0.5 : 1;
      out.push(ev({
        rule: 'alkaneRH', ruleName: 'Eliminación de un alcano R–H (alcanos ramificados → ion alqueno)', ion, neutral: neu,
        neutralLabel: `alcano ${neu.formulaString}`, cleaved: [b.idx], oe: true,
        score: 0.5 * wDeg * (ALKENE_SIZE_W[nIon] || 0.4) * wN,
        steps: [
          'Ionización σ del esqueleto (M⁺• deslocalizado).',
          `En el punto de ramificación ${L(g, ca)} se rompe ${B(g, ca, cg)} con transferencia de un H desde ${L(g, cb.atom)} al fragmento saliente (análogo a la deshidratación de alcoholes).`,
          `Se elimina el alcano neutro ${neu.formulaString} y queda el radical-catión alqueno ${ion.formulaString}⁺• (masa par, m/z ${ion.nominal}), que a su vez origina fragmentos alquénicos (alílicos, McLafferty).`,
          'Explica iones pares CₙH₂ₙ⁺• (m/z 56, 70, 84…) en alcanos; puede ser incluso el pico base (3,4-dimetilhexano, m/z 56).'
        ],
        refs: ['BEAU', 'MT']
      }));
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* B4. Eliminación de ROH, RSH, NH₃/RNH₂, RCOOH                         */
/* ------------------------------------------------------------------ */
function eliminationHXR(g) {
  const out = []; const seen = new Set();
  for (const X of g.atoms) {
    if (X.aromatic || !['O', 'S', 'N'].includes(X.el) || !g.isSp3(X.idx)) continue;
    let type = null; let base = 0;
    const carbNb = X.nbrs.filter(n => g.carbonylO(n.atom) != null);
    if (X.el === 'O' && X.h === 0 && X.nbrs.length === 2 && X.nbrs.every(n => g.isC(n.atom) && !g.atoms[n.atom].aromatic)) {
      if (!carbNb.length) { type = 'ether'; base = 0.3; } // ésteres: ya cubiertos por el McLafferty del lado alcoxilo (v1.1)
    } else if (X.el === 'S' && X.h === 0 && X.nbrs.length === 2 && X.nbrs.every(n => g.isC(n.atom) && !g.atoms[n.atom].aromatic)) {
      type = 'sulfide'; base = 0.3;
    } else if (X.el === 'N' && !carbNb.length && X.nbrs.every(n => g.isC(n.atom) && !g.atoms[n.atom].aromatic)) {
      type = X.h === 2 ? 'amine1' : 'amine'; base = X.h === 2 ? 0.12 : 0.06;
    }
    if (!type) continue;
    for (const nx of X.nbrs) {
      const ca = nx.atom;
      if (g.carbonylO(ca) != null || !g.isSp3(ca) || g.bonds[nx.bond].inRing) continue;
      const side = g.sideOf(nx.bond, X.idx, allAtoms(g));
      const cb = g.atoms[ca].nbrs.find(n => n.atom !== X.idx && !side.has(n.atom) && g.isC(n.atom) && g.isSp3(n.atom) && g.atoms[n.atom].h > 0);
      if (!cb) continue;
      const ionAt = complement(g, side);
      const key = [...ionAt].sort((a, b) => a - b).join(',');
      if (seen.has(key)) continue; seen.add(key);
      const ion = new Species(g, ionAt, { dH: new Map([[cb.atom, -1]]), bo: new Map([[cb.bond, 2]]), oddElectron: true, depictAsNeutral: true, tags: { alkeneIon: true } });
      const neu = new Species(g, side, { dH: new Map([[X.idx, 1]]) });
      const nC = [...ionAt].filter(i => g.isC(i)).length;
      const lab = { ester: 'ácido carboxílico', ether: 'alcohol', sulfide: 'tiol', amine1: 'NH₃', amine: 'amina' }[type];
      const nm = { ester: 'Eliminación de RCOOH (ésteres de alquilo)', ether: 'Eliminación de ROH (éteres)', sulfide: 'Eliminación de RSH (sulfuros)', amine1: 'Eliminación de NH₃ (aminas primarias)', amine: 'Eliminación de una amina RNH₂/R₂NH' }[type];
      out.push(ev({
        rule: 'elim_' + type, ruleName: nm, ion, neutral: neu, neutralLabel: `${lab} ${neu.formulaString}`, cleaved: [nx.bond], oe: true,
        score: base * (nC >= 3 ? 1 : 0.15),
        steps: [
          `Ionización en el par libre de ${L(g, X.idx)}.`,
          type === 'ester'
            ? `Transferencia de un H β (${L(g, cb.atom)}) del grupo alcoxilo al oxígeno carbonílico (estado de transición cíclico de 6 miembros, tipo McLafferty).`
            : `Transferencia de un H desde ${L(g, cb.atom)} (posición β, 1,2; o más alejada, 1,4) al heteroátomo ${L(g, X.idx)}.`,
          `Ruptura de ${B(g, ca, X.idx)} y expulsión de la molécula neutra ${neu.formulaString} ([M−${g.nominalMass - ion.nominal}]⁺•); queda el radical-catión alqueno ${ion.formulaString}⁺• (masa par si no hay N).`,
          'Igual que la pérdida de H₂O en alcoholes, el ion alqueno resultante genera luego fragmentos alílicos y de McLafferty.'
        ],
        refs: ['BEAU', 'MT']
      }));
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* B6. Fragmentación alquénica de iones alqueno secundarios            */
/* ------------------------------------------------------------------ */
function derivedAlkeneEvents(RDKit, e) {
  const sp = e.ion;
  if (!sp || sp.formulaOverride || !sp.atoms.size) return [];
  const js = sp.toRDKitJSON({ neutralSkeleton: true });
  if (!js) return [];
  let mol = null; let g2;
  try {
    mol = RDKit.get_mol(JSON.stringify(js));
    if (!mol || !mol.is_valid()) return [];
    g2 = new MolGraph(mol.get_json());
  } catch (_) { return []; } finally { if (mol) try { mol.delete(); } catch (_) { /* */ } }
  const order = [...sp.atoms].sort((a, b) => a - b);
  if (g2.n !== order.length) return [];
  g2.atoms.forEach((a, i) => { const o = sp.g.atoms[order[i]]; if (o && o.el === a.el) a.orig = o.orig ?? order[i]; });
  let evs = [];
  try { evs.push(...benzylicAllylic(g2).filter(x => x.rule === 'allylic' || x.rule === 'propargylic')); } catch (_) { /* */ }
  try { evs.push(...mclafferty(g2).filter(x => { const t = x.ion.tags.mclafferty || (x.neutral && x.neutral.tags.mclafferty); return t && t.kind === 'alqueno'; })); } catch (_) { /* */ }
  try { evs.push(...migration(g2)); } catch (_) { /* */ }
  evs = evs.filter(x => { try { return x.ion.nominal > 0 && x.score > 0; } catch (_) { return false; } });
  if (!evs.length) return [];
  const S = evs.reduce((s, x) => s + x.score, 0);
  return evs.map(x => Object.assign(x, {
    k: 0.75 * x.score / (S + 1.2),
    generation: e.generation + 1, parent: e, derived: true, cleaved: [],
    ruleName: 'Fragmentación del ion alqueno: ' + x.ruleName.charAt(0).toLowerCase() + x.ruleName.slice(1),
    steps: [`El radical-catión alqueno ${sp.formulaString}⁺• (m/z ${sp.nominal}) se comporta como un alqueno ionizado (índices de átomos referidos a la molécula original).`, ...x.steps],
    refs: [...new Set([...(x.refs || []), 'BEAU'])]
  }));
}

/* ------------------------------------------------------------------ */
/* B7. Alcoholes y tioles primarios: −(H₂O + C₂H₄) / −(H₂S + C₂H₄)      */
/* ------------------------------------------------------------------ */
function waterEthylene(g) {
  const out = [];
  for (const X of g.atoms) {
    if (!['O', 'S'].includes(X.el) || X.h !== 1 || X.nbrs.length !== 1) continue;
    const ca = X.nbrs[0].atom;
    if (!g.isC(ca) || g.atoms[ca].aromatic || !g.isSp3(ca) || g.atoms[ca].h !== 2) continue;
    // cadena lineal Cα–Cβ–Cγ–Cδ–Cε (se requieren ≥ 5 C para que el ion remanente sea ≥ C₃)
    let chain = [ca]; let ok = true;
    while (chain.length < 5) {
      const last = chain[chain.length - 1];
      const nx = g.atoms[last].nbrs.find(n => !chain.includes(n.atom) && n.atom !== X.idx && g.isC(n.atom) && g.isSp3(n.atom) && !g.bonds[n.bond].inRing);
      if (!nx) { ok = false; break; }
      chain.push(nx.atom);
    }
    if (!ok) continue;
    const f = g.formula; const el = X.el;
    const fx = R.clean({ ...f, C: f.C - 2, H: f.H - 6, [el]: f[el] - 1 });
    out.push(ev({
      rule: 'lossH2O_C2H4', ruleName: el === 'O' ? 'Pérdida consecutiva de H₂O y C₂H₄ (alcoholes primarios)' : 'Pérdida consecutiva de H₂S y C₂H₄ (tioles primarios)',
      ion: new Species(g, new Set(), { formula: fx, oddElectron: true }), neutral: null,
      neutralFormula: el === 'O' ? { C: 2, H: 6, O: 1 } : { C: 2, H: 6, S: 1 }, neutralLabel: `H₂${el} + C₂H₄`, oe: true, cleaved: [],
      score: el === 'O' ? 0.55 : 0.6,
      steps: [
        `Ionización en ${L(g, X.idx)} y eliminación 1,4 de H₂${el} (H del carbono δ, estado de transición de 6 miembros) → ion radical distónico/ciclobutano.`,
        'El ion [M−18]⁺• (o [M−34]⁺•) expulsa etileno (escisión del anillo/escisión alílica tras migración) formando el radical-catión alqueno [M−H₂X−C₂H₄]⁺• de masa par.',
        `Explica el pico (a menudo base) m/z ${(g.nominalMass - (el === 'O' ? 46 : 62))} en ${el === 'O' ? '1-alcanoles' : '1-alcanotioles'} ≥ C₅ (p. ej., 1-hexanol y 1-hexanotiol: m/z 56).`
      ],
      refs: el === 'O' ? ['ALCOH', 'MT', 'BEAU'] : ['MT', 'BEAU']
    }));
  }
  return out;
}

const BEAU_RULES = [migration, propargylicH, alkaneRH, eliminationHXR, waterEthylene];

module.exports = { BEAU_RULES, derivedAlkeneEvents, migration };
