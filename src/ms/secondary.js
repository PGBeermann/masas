'use strict';
/**
 * Fragmentaciones secundarias (consecutivas) de iones par-electrónicos (EE⁺) y algunos OE⁺•.
 * Regla del electrón par (Karni & Mandelbaum, Org. Mass Spectrom. 15 (1980) 53): los iones EE⁺
 * se fragmentan preferentemente expulsando moléculas neutras (CO, alquenos, C₂H₂, H₂), no radicales.
 */
const { Species } = require('./species');
const { cationClass } = require('./stability');
const { clean, L, B } = require('./rules');
const { formulaToString } = require('../chem/elements');

function child(parent, o) {
  return Object.assign({ steps: [], refs: [], children: [], generation: parent.generation + 1, parent, oe: false, cleaved: [] }, o);
}
const fo = (g, f, tags = {}) => new Species(g, new Set(), { formula: clean(f), tags });

function secondary(e) {
  const out = []; const ion = e.ion; const g = ion.g;
  const f = ion.formula; const fs = ion.formulaString;

  // (1) Acilio → R⁺ + CO
  if (ion.tags.acylium && !ion.formulaOverride) {
    const { c, o } = ion.tags.acylium;
    const r = g.atoms[c].nbrs.find(n => n.atom !== o && ion.atoms.has(n.atom));
    if (r) {
      const atoms = new Set(ion.atoms); atoms.delete(c); atoms.delete(o);
      const dH = new Map(ion.dH); dH.delete(c);
      const bo = new Map(ion.bo); bo.delete(g.bondBetween(c, o));
      const chg = new Map(ion.chg); chg.delete(o);
      let tags = {};
      const cc = cationClass(g, atoms, r.atom);
      if (g.atoms[r.atom].el === 'C') chg.set(r.atom, 1); else tags.noDepict = true;
      const nf = g.formulaOf(atoms, dH);
      if (nf.C && nf.H === 2 * nf.C + 1 && Object.keys(nf).length === 2) tags.alkyl = nf.C;
      if (cc.cls === 'aryl') tags.aryl = true;
      const ni = new Species(g, atoms, { dH, bo, chg, tags });
      const k = cc.cls === 'aryl' ? 0.55 : Math.min(0.6, 0.6 * cc.score);
      out.push(child(e, { rule: 'acyliumCO', ruleName: 'Descarbonilación del ion acilio (−CO)', ion: ni, neutral: null, neutralFormula: { C: 1, O: 1 }, neutralLabel: 'CO', k,
        steps: [`El ion acilio ${fs}⁺ expulsa monóxido de carbono (molécula neutra estable, regla del electrón par) por ruptura heterolítica de ${B(g, c, r.atom)} (⇒).`,
          `Se forma ${cc.label} ${ni.formulaString}⁺ (Δm = 28). ${cc.cls === 'aryl' ? 'Ejemplo clásico: C₆H₅CO⁺ (m/z 105) → C₆H₅⁺ (m/z 77).' : ''}`],
        refs: ['MT', 'GROSS'] }));
    }
  }

  // (2) Reacción del onio: iminio/oxonio con N/O-alquilo → pérdida de alqueno con transferencia de H
  if (ion.tags.onium && !ion.formulaOverride) {
    const { het, c } = ion.tags.onium;
    for (const n of g.atoms[het].nbrs) {
      if (n.atom === c || !ion.atoms.has(n.atom) || !g.isC(n.atom) || !g.isSp3(n.atom) || g.bonds[n.bond].inRing) continue;
      const side = g.sideOf(n.bond, n.atom, ion.atoms);
      if (!side || side.has(het)) continue;
      const beta = g.atoms[n.atom].nbrs.find(m => side.has(m.atom) && g.isC(m.atom) && g.isSp3(m.atom) && (g.atoms[m.atom].h + (ion.dH.get(m.atom) || 0)) > 0);
      if (!beta) continue;
      const atoms = new Set([...ion.atoms].filter(a => !side.has(a)));
      const dH = new Map(ion.dH); dH.set(het, (dH.get(het) || 0) + 1);
      const ni = new Species(g, atoms, { dH, bo: new Map(ion.bo), chg: new Map(ion.chg), tags: { onium: { het, c } } });
      const neu = new Species(g, side, { dH: new Map([[beta.atom, -1]]), bo: new Map([[g.bondBetween(n.atom, beta.atom), 2]]) });
      const k = { N: 0.3, O: 0.8, S: 0.35 }[g.atoms[het].el] || 0.3;
      out.push(child(e, { rule: 'onium', ruleName: 'Reacción del onio (pérdida de alqueno desde el ion iminio/oxonio)', ion: ni, neutral: neu, neutralLabel: `alqueno ${neu.formulaString}`, k,
        steps: [`El ion ${fs}⁺ (EE) con un sustituyente alquilo sobre ${L(g, het)} transfiere un H β (${L(g, beta.atom)}) al heteroátomo (estado de transición de 4 miembros o vía complejo ion–neutro).`,
          `Ruptura heterolítica de ${B(g, het, n.atom)} y eliminación del alqueno neutro ${neu.formulaString}; se forma ${ni.formulaString}⁺.`,
          'Proceso típico de aminas y éteres (p. ej., trietilamina m/z 86 → 58; dietil éter m/z 59 → 31).'],
        refs: ['AMINES', 'MT'] }));
    }
  }

  // (3) Cationes alquilo
  if (ion.tags.alkyl) {
    const n = f.C;
    if (n >= 2 && n <= 5) {
      out.push(child(e, { rule: 'lossH2', ruleName: 'Pérdida de H₂ en carbocationes alquilo', ion: fo(g, { C: n, H: 2 * n - 1 }, { allylic: true }), neutral: null, neutralFormula: { H: 2 }, neutralLabel: 'H₂',
        k: { 2: 0.55, 3: 0.3, 4: 0.25, 5: 0.2 }[n],
        steps: [`El catión ${fs}⁺ elimina H₂ (1,2-eliminación) generando el catión alquenilo/alilo CₙH₂ₙ₋₁⁺ (m/z ${12 * n + 2 * n - 1}).`], refs: ['MT'] }));
    }
    if (n === 4) {
      out.push(child(e, { rule: 'lossCH4', ruleName: 'Pérdida de CH₄ en el catión butilo', ion: fo(g, { C: 3, H: 5 }, { allylic: true }), neutral: null, neutralFormula: { C: 1, H: 4 }, neutralLabel: 'CH₄', k: 0.25,
        steps: ['El catión C₄H₉⁺ (m/z 57), tras isomerizarse al catión terc-butilo/sec-butilo, elimina metano generando el catión alilo C₃H₅⁺ (m/z 41).'], refs: ['MT'] }));
    }
    if (n >= 5) {
      for (const m of [3, 4]) {
        out.push(child(e, { rule: 'alkylAlkene', ruleName: 'Pérdida de alqueno en carbocationes alquilo (reordenamiento por desplazamientos de hidruro)', ion: fo(g, { C: m, H: 2 * m + 1 }, { alkyl: m }), neutral: null,
          neutralFormula: { C: n - m, H: 2 * (n - m) }, neutralLabel: `C${n - m}H${2 * (n - m)}`, k: m === 3 ? 0.22 : 0.2,
          steps: [`El catión ${fs}⁺ se isomeriza (desplazamientos 1,2-H⁻ y alquilo) a estructuras secundarias/terciarias y elimina un alqueno neutro C${n - m}H${2 * (n - m)}, originando C${m}H${2 * m + 1}⁺ (m/z ${14 * m + 1}).`], refs: ['MT', 'GROSS'] }));
      }
    }
  }

  // (4) Iones aromáticos característicos
  const simple = (from, toF, lossF, lossLabel, k, txt) => {
    if (fs !== from) return;
    out.push(child(e, { rule: 'arom', ruleName: `Degradación de ${from}⁺ (−${lossLabel})`, ion: fo(g, toF), neutral: null, neutralFormula: lossF, neutralLabel: lossLabel, k, steps: [txt], refs: ['MT', 'TROP'] }));
  };
  simple('C7H7', { C: 5, H: 5 }, { C: 2, H: 2 }, 'C₂H₂', 0.14, 'El ion tropilio C₇H₇⁺ (m/z 91) pierde acetileno dando el catión ciclopentadienilo C₅H₅⁺ (m/z 65).');
  simple('C5H5', { C: 3, H: 3 }, { C: 2, H: 2 }, 'C₂H₂', 0.45, 'C₅H₅⁺ (m/z 65) pierde acetileno dando el catión ciclopropenilo C₃H₃⁺ (m/z 39), el carbocatión aromático más pequeño.');
  simple('C6H5', { C: 4, H: 3 }, { C: 2, H: 2 }, 'C₂H₂', 0.45, 'El catión fenilo C₆H₅⁺ (m/z 77) pierde acetileno: C₄H₃⁺ (m/z 51), serie aromática 77 → 51.');
  simple('C6H5O', { C: 5, H: 5 }, { C: 1, O: 1 }, 'CO', 0.7, 'El ion fenoxilo C₆H₅O⁺ (m/z 93) expulsa CO dando C₅H₅⁺ (m/z 65).');
  simple('C7H5O', { C: 6, H: 5 }, { C: 1, O: 1 }, 'CO', 0, '');
  return out.filter(x => x.k > 0);
}

module.exports = { secondary };
