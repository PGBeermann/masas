'use strict';
/**
 * Motor de predicción de espectros EI-MS (70 eV) basado en reglas mecanísticas.
 *
 * Flujo:
 *   SMILES → RDKit (validación, canonización, aromaticidad, anillos) → MolGraph
 *   → reglas primarias sobre M⁺• → fragmentaciones secundarias (≤ 2 generaciones)
 *   → agregación por composición elemental → patrón isotópico → normalización (pico base = 100)
 *   → explicación mecanística de los iones principales.
 *
 * Las intensidades son estimaciones SEMICUANTITATIVAS (orden de magnitud/ranking), no un cálculo
 * cinético RRKM/QET. Sirven para docencia e hipótesis de asignación, no sustituyen la comparación
 * con espectros de referencia (NIST/Wiley).
 */
const { MolGraph } = require('../chem/molgraph');
const { formulaToString, monoMass, nominalMass, ELECTRON_MASS, rdb } = require('../chem/elements');
const { isotopePattern } = require('../chem/isotopes');
const { Species } = require('./species');
const { PRIMARY_RULES, clean } = require('./rules');
const { EXT_RULES } = require('./rules_ext');
const { BEAU_RULES, derivedAlkeneEvents } = require('./rules_beauchamp');
const { secondary } = require('./secondary');
const REFS = require('./references');

const MAX_HEAVY = 60;

function largestFragment(smiles) {
  const parts = smiles.split('.').filter(Boolean);
  return parts.sort((a, b) => b.length - a.length)[0];
}

/** Factor de estabilidad del ion molecular según la clase de compuesto (McLafferty & Tureček 1993, tabla 3.2). */
function molecularIonFactor(g) {
  const arom = g.atoms.some(a => a.aromatic);
  const nRing = g.rings.length;
  let mf = arom ? 1.25 : nRing ? 1.4 : 1.0;
  const notes = [];
  let pen = 1;
  const nCdb = g.bonds.filter(b => b.order === 2 && !b.aromatic && g.isC(b.a) && g.isC(b.b)).length;
  if (!arom && nCdb >= 2) mf *= 1.5;
  for (const a of g.atoms) {
    if (a.el === 'O' && a.h === 1 && a.nbrs.length === 1 && !g.atoms[a.nbrs[0].atom].aromatic && g.carbonylO(a.nbrs[0].atom) == null) {
      const deg = g.atoms[a.nbrs[0].atom].nbrs.filter(n => g.isC(n.atom)).length;
      pen *= [0.12, 0.12, 0.06, 0.015][deg]; notes.push('alcohol');
    }
    if (a.el === 'N' && !a.aromatic && g.isSp3(a.idx) && !a.nbrs.some(n => g.carbonylO(n.atom) != null || g.atoms[n.atom].aromatic)) { pen *= 0.4; notes.push('amina'); }
    if (a.el === 'O' && a.nbrs.length === 2 && a.nbrs.every(n => g.isC(n.atom) && g.carbonylO(n.atom) == null) && !a.aromatic) { pen *= 0.25; notes.push('éter'); }
    if (a.el === 'C' && !a.aromatic && g.carbonylO(a.idx) != null) {
      const fam = a.nbrs.some(n => g.atoms[n.atom].el === 'O' && g.bonds[n.bond].order === 1) ? 0.35 : a.h ? 0.5 : 0.7;
      pen *= fam; notes.push('carbonilo');
    }
    if (['Cl', 'Br', 'I'].includes(a.el) && a.nbrs.length && !g.atoms[a.nbrs[0].atom].aromatic) { pen *= { Cl: 0.3, Br: 0.25, I: 0.8 }[a.el]; notes.push('haluro de alquilo'); }
    if (a.el === 'N' && g.tripleBondPartner(a.idx, 'C') != null && !a.nbrs.some(n => g.atoms[n.atom].aromatic)) { pen *= 0.15; notes.push('nitrilo alifático'); }
    if (a.el === 'C' && a.h === 1 && g.tripleBondPartner(a.idx, 'C') != null) { pen *= 0.1; notes.push('alquino terminal'); }
    if (a.el === 'C' && !a.aromatic && g.isSp3(a.idx) && !g.inRing(a.idx)) {
      const cdeg = a.nbrs.filter(n => g.isC(n.atom)).length;
      if (cdeg === 4) pen *= 0.12; else if (cdeg === 3) pen *= 0.45;
    }
  }
  if (arom) pen = Math.sqrt(pen);
  const nSp3 = g.atoms.filter(a => a.el === 'C' && g.isSp3(a.idx) && !g.inRing(a.idx)).length;
  mf *= Math.exp(-(arom ? 0.06 : 0.25) * Math.max(0, nSp3 - 3));
  return { factor: mf * pen, notes: [...new Set(notes)] };
}

/** Reglas adicionales para núcleos aromáticos sin sustituyentes fragmentables. */
function aromaticCore(g) {
  const out = [];
  const onlyCH = g.atoms.every(a => a.el === 'C' && a.aromatic);
  const arom = g.atoms.filter(a => a.aromatic).length;
  if (!arom) return out;
  const f = g.formula;
  const hasArH = g.atoms.some(a => a.aromatic && a.h > 0);
  if (hasArH) out.push({ rule: 'arH', ruleName: 'Pérdida de H• aromático', ion: new Species(g, new Set(), { formula: clean({ ...f, H: f.H - 1 }) }), neutral: null,
    neutralFormula: { H: 1 }, neutralLabel: 'H•', oe: false, score: onlyCH ? 0.35 : 0.08, cleaved: [], generation: 1,
    steps: ['Ionización π del anillo aromático (M⁺• muy estable).', 'Ruptura C–H del anillo con formación de un catión arilo [M−H]⁺ (proceso de alta energía, intensidad moderada).'], refs: ['MT'] });
  if (onlyCH) out.push({ rule: 'arC2H2', ruleName: 'Expulsión de acetileno desde el anillo', ion: new Species(g, new Set(), { formula: clean({ ...f, C: f.C - 2, H: f.H - 2 }), oddElectron: true }), neutral: null,
    neutralFormula: { C: 2, H: 2 }, neutralLabel: 'C₂H₂', oe: true, score: 0.3, cleaved: [], generation: 1,
    steps: ['El radical-catión aromático se isomeriza (apertura del anillo) y expulsa C₂H₂ ([M−26]⁺•), típico de hidrocarburos aromáticos (benceno 78 → 52).'], refs: ['MT', 'GROSS'] });
  return out;
}

function predict(RDKit, smilesIn, opts = {}) {
  const t0 = Date.now();
  const topN = opts.topExplain || 12;
  if (!smilesIn || typeof smilesIn !== 'string') throw new Error('Debe proporcionar un SMILES.');
  const smiles = largestFragment(smilesIn.trim());
  const mol = RDKit.get_mol(smiles);
  if (!mol || !mol.is_valid()) { if (mol) mol.delete(); throw new Error('SMILES inválido o no sanitizable por RDKit.'); }
  let g;
  try {
    g = new MolGraph(mol.get_json());
    if (g.n > MAX_HEAVY) throw new Error(`La molécula excede el límite de ${MAX_HEAVY} átomos pesados.`);
    if (g.atoms.reduce((s, a) => s + a.chg, 0) !== 0 || g.atoms.some(a => a.rad !== 0)) throw new Error('Introduzca una molécula neutra de capa cerrada (carga neta 0, sin radicales).');
    if (g.n < 2) throw new Error('Se requieren al menos dos átomos pesados.');
    const canonical = mol.get_smiles();
    const desc = JSON.parse(mol.get_descriptors());
    const M = new Species(g, new Set(g.atoms.map(a => a.idx)), { oddElectron: true, depictAsNeutral: true });

    // ---------- reglas primarias ----------
    let events = [];
    for (const rule of [...PRIMARY_RULES, ...EXT_RULES, ...BEAU_RULES]) {
      try { events.push(...rule(g)); } catch (e) { /* regla no aplicable */ }
    }
    events.push(...aromaticCore(g));
    events = events.filter(e => { try { return e.ion.nominal > 0 && e.score > 0; } catch (_) { return false; } });

    // ---------- secundarias (2 generaciones) ----------
    const all = [];
    let frontier = events;
    for (let gen = 1; gen <= 4 && frontier.length; gen++) {
      const next = [];
      for (const e of frontier) {
        all.push(e);
        if (gen >= 4) continue;
        let kids = [];
        try { kids = secondary(e); } catch (_) { kids = []; }
        // iones alqueno OE⁺• formados por eliminación o McLafferty → fragmentación alquénica consecutiva (Beauchamp)
        if (e.ion.tags && e.ion.tags.alkeneIon && !e.derived && gen <= 2) {
          try { kids.push(...derivedAlkeneEvents(RDKit, e)); } catch (_) { /* no aplicable */ }
        }
        let kSum = 0;
        for (const c of kids) { c.score = e.score * c.k; kSum += c.k; next.push(c); }
        if (kSum) e.score *= Math.max(0.35, 1 - 0.5 * kSum);
      }
      frontier = next;
    }

    // ---------- agregación por composición ----------
    const byFormula = new Map();
    for (const e of all) {
      let fs; try { fs = e.ion.formulaString; } catch (_) { continue; }
      if (!byFormula.has(fs)) byFormula.set(fs, { formula: e.ion.formula, formulaString: fs, score: 0, events: [] });
      const agg = byFormula.get(fs);
      agg.score += e.score; agg.events.push(e);
    }
    const mInfo = molecularIonFactor(g);
    const Mfs = M.formulaString;
    const mEvent = {
      rule: 'M', ruleName: 'Ion molecular M⁺•', ion: M, neutral: null, oe: true, generation: 0, cleaved: [],
      steps: [
        'Impacto electrónico (70 eV): M + e⁻ → M⁺• + 2e⁻. Se remueve el electrón de menor energía de ionización: n (N > S > O ≈ I > Br > Cl) > π (aromático, C=C) > σ.',
        `Estabilidad estimada del M⁺•: factor ${mInfo.factor.toFixed(2)}${mInfo.notes.length ? ' (clases: ' + mInfo.notes.join(', ') + ')' : ''}. Los sistemas aromáticos/conjugados dan M⁺• intensos; alcoholes, aminas, éteres y alcanos ramificados, débiles.`
      ],
      refs: ['MT', 'GROSS'], score: 0
    };
    const fragMax = Math.max(0.0001, ...[...byFormula.values()].map(a => a.score));
    mEvent.score = Math.max(fragMax * 0.005, 1.1 * mInfo.factor);
    if (!byFormula.has(Mfs)) byFormula.set(Mfs, { formula: M.formula, formulaString: Mfs, score: 0, events: [] });
    byFormula.get(Mfs).score += mEvent.score; byFormula.get(Mfs).events.unshift(mEvent);

    // ---------- espectro con isótopos ----------
    const peaks = new Map();
    for (const agg of byFormula.values()) {
      const nom = nominalMass(agg.formula);
      agg.nominal = nom; agg.mz = monoMass(agg.formula) - ELECTRON_MASS;
      const pat = isotopePattern(agg.formula);
      const p0 = pat[0].prob;
      for (const p of pat) {
        const mz = nom + p.offset; const inten = agg.score * p.prob / p0;
        if (!peaks.has(mz)) peaks.set(mz, { mz, raw: 0, contrib: [] });
        const pk = peaks.get(mz); pk.raw += inten;
        pk.contrib.push({ formula: agg.formulaString, isotope: p.offset > 0, offset: p.offset, raw: inten });
      }
    }
    const maxRaw = Math.max(...[...peaks.values()].map(p => p.raw));
    const minI = opts.threshold ?? 0.5;
    const spectrum = [...peaks.values()].map(p => ({
      mz: p.mz, intensity: +(100 * p.raw / maxRaw).toFixed(2),
      assignments: p.contrib.sort((a, b) => b.raw - a.raw).map(c => ({ formula: c.formula, isotope: c.isotope, offset: c.offset, share: +(c.raw / p.raw).toFixed(3) }))
    })).filter(p => p.intensity >= minI && p.mz >= 12).sort((a, b) => a.mz - b.mz);
    const base = spectrum.reduce((a, b) => (b.intensity > a.intensity ? b : a), spectrum[0]);

    // ---------- explicación de los iones principales ----------
    const ranked = [...byFormula.values()].map(a => ({ ...a, intensity: 100 * a.score / maxRaw })).sort((a, b) => b.intensity - a.intensity);
    const top = ranked.filter(a => a.intensity >= 1).slice(0, topN);
    if (!top.find(a => a.formulaString === Mfs)) top.push(ranked.find(a => a.formulaString === Mfs));
    const parentSvgFor = (atoms, bonds) => {
      try {
        const set = new Set(atoms); const cut = new Set(bonds || []);
        const keep = g.bonds.filter(b => set.has(b.a) && set.has(b.b) && !cut.has(b.idx)).map(b => b.idx);
        const hb = {}; keep.forEach(i => { hb[i] = [0.36, 0.62, 0.95]; }); cut.forEach(i => { hb[i] = [0.95, 0.45, 0.15]; });
        return mol.get_svg_with_highlights(JSON.stringify({ width: 300, height: 220, addAtomIndices: true, atoms: [...set], bonds: [...keep, ...cut],
          highlightColour: [0.36, 0.62, 0.95], highlightBondColors: hb, highlightBondWidthMultiplier: 10 }));
      } catch (_) { return null; }
    };
    const usedRefs = new Set(['ISO', 'NIST', 'STEV']);
    const explanations = top.map(a => {
      const evs = a.events.slice().sort((x, y) => y.score - x.score);
      const best = evs[0];
      const pathway = [];
      for (let e = best; e; e = e.parent) pathway.unshift(e);
      const steps = pathway.map((e, i) => {
        (e.refs || []).forEach(r => usedRefs.add(r));
        let ionDep = null; let neuDep = null;
        if (!e.ion.tags.noDepict) ionDep = e.ion.depict(RDKit);
        if (e.neutral && !e.neutral.tags.noDepict) neuDep = e.neutral.depict(RDKit, { w: 160, h: 120 });
        let neutralF = null;
        try { neutralF = e.neutral ? e.neutral.formulaString : e.neutralFormula ? formulaToString(e.neutralFormula) : null; } catch (_) { /* */ }
        return {
          generation: e.generation, rule: e.rule, ruleName: e.ruleName,
          ionFormula: e.ion.formulaString, ionMz: e.ion.nominal, oddElectron: !!e.oe,
          neutralFormula: neutralF, neutralLabel: e.neutralLabel || null,
          ionSmiles: ionDep ? ionDep.smiles : null, ionSvg: ionDep ? ionDep.svg : null,
          neutralSmiles: neuDep ? neuDep.smiles : null, neutralSvg: neuDep ? neuDep.svg : null,
          mechanism: e.steps, refs: e.refs
        };
      });
      const hl = best.ion.atoms && best.ion.atoms.size ? best.ion.atoms : (pathway[0].ion.atoms || new Set());
      const r = rdb(a.formula);
      return {
        formula: a.formulaString, mz: a.nominal, mzExact: +a.mz.toFixed(4), intensity: +a.intensity.toFixed(1),
        ionType: Number.isInteger(r) ? 'OE⁺• (impar-electrónico)' : 'EE⁺ (par-electrónico)', rdb: r,
        lossFromM: g.nominalMass - a.nominal,
        parentHighlightSvg: a.formulaString === Mfs ? parentSvgFor([]) : parentSvgFor(hl, pathway[0].cleaved),
        pathway: steps,
        alternatives: Object.values(evs.slice(1).reduce((acc, e) => { (acc[e.ruleName] = acc[e.ruleName] || { ruleName: e.ruleName, share: 0 }).share += e.score / a.score; return acc; }, {})).filter(x => x.ruleName !== best.ruleName).map(x => ({ ruleName: x.ruleName, share: +x.share.toFixed(2) })).slice(0, 4)
      };
    });

    const nN = g.formula.N || 0;
    const result = {
      input: smilesIn, smiles: canonical,
      formula: g.formulaString, nominalMass: g.nominalMass, monoisotopicMass: +g.monoisotopicMass.toFixed(5),
      rdkitExactMW: desc.exactmw, rdb: g.rdb,
      nitrogenRule: `Nº de N = ${nN} → M nominal ${nN % 2 ? 'impar' : 'par'} (${g.nominalMass}).`,
      molecularIon: { mz: g.nominalMass, factor: +mInfo.factor.toFixed(3), classes: mInfo.notes },
      basePeak: base ? { mz: base.mz, formula: base.assignments[0].formula } : null,
      spectrum, explanations,
      parentSvg: parentSvgFor([]),
      references: [...usedRefs].map(k => ({ key: k, text: REFS[k] })).filter(r => r.text),
      meta: { engine: 'EI-MS rule-based predictor v1.2', rdkit: RDKit.version(), events: all.length, ms: Date.now() - t0,
        disclaimer: 'Intensidades semicuantitativas derivadas de reglas mecanísticas (EI 70 eV). Validar con espectros de referencia (NIST/Wiley).' }
    };
    return result;
  } finally {
    mol.delete();
  }
}

module.exports = { predict };
