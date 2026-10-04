'use strict';
/**
 * Fragmentaciones secundarias (consecutivas) de iones par-electrónicos (EE⁺) y algunos OE⁺•.
 * Regla del electrón par (Karni & Mandelbaum, Org. Mass Spectrom. 15 (1980) 53): los iones EE⁺
 * se fragmentan preferentemente expulsando moléculas neutras (CO, alquenos, C₂H₂, H₂), no radicales.
 */
const { Species } = require('./species');
const { cationClass } = require('./stability');
const { clean, L, B, FH, FL, A, Mid } = require('./rules');
const { formulaToString } = require('../chem/elements');
const { orthoEffect } = require('./rules_ext');

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
      const dist = !!ion.tags.distonic;
      const ni = new Species(g, atoms, { dH, bo, chg, rad: new Map(ion.rad), tags: dist ? { ...tags, noDepict: true, alkeneOE: true } : tags, oddElectron: dist });
      const k = dist ? (atoms.size >= 3 ? 1.4 : 0.8) : cc.cls === 'aryl' ? 0.55 : Math.min(0.6, 0.6 * cc.score);
      out.push(child(e, { mech: { site: { kind: 'n', a: o, label: dist ? '+' : '+' }, arrows: [FL(Mid(c, r.atom), A(c))], highlightBonds: [[c, r.atom]] }, rule: 'acyliumCO', ruleName: dist ? 'Descarbonilación del ion acilio distónico (−CO)' : 'Descarbonilación del ion acilio (−CO)', ion: ni, oe: dist, neutral: null, neutralFormula: { C: 1, O: 1 }, neutralLabel: 'CO', k,
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
      out.push(child(e, { mech: { site: { kind: 'n', a: het, label: '+' }, arrows: [FL(Mid(het, n.atom), ['mH', beta.atom, het, het]), FL(['mHi', beta.atom, het], Mid(n.atom, beta.atom))], highlightBonds: [[het, n.atom]] }, rule: 'onium', ruleName: 'Reacción del onio (pérdida de alqueno desde el ion iminio/oxonio)', ion: ni, neutral: neu, neutralLabel: `alqueno ${neu.formulaString}`, k,
        steps: [`El ion ${fs}⁺ (EE) con un sustituyente alquilo sobre ${L(g, het)} transfiere un H β (${L(g, beta.atom)}) al heteroátomo (estado de transición de 4 miembros o vía complejo ion–neutro).`,
          `Ruptura heterolítica de ${B(g, het, n.atom)} y eliminación del alqueno neutro ${neu.formulaString}; se forma ${ni.formulaString}⁺.`,
          'Proceso típico de aminas y éteres (p. ej., trietilamina m/z 86 → 58; dietil éter m/z 59 → 31).'],
        refs: ['AMINES', 'MT'] }));
    }
  }

  // (2b) Onio con alquilo sobre el carbono (alcoholes): R–CH=OH⁺ → CH₂=OH⁺ + alqueno
  if (ion.tags.onium && !ion.formulaOverride && !ion.oddElectron) {
    const { het, c } = ion.tags.onium;
    const hetEl = g.atoms[het].el;
    for (const n of g.atoms[c].nbrs) {
      if (n.atom === het || !ion.atoms.has(n.atom) || !g.isC(n.atom) || !g.isSp3(n.atom) || g.bonds[n.bond].inRing) continue;
      const side = g.sideOf(n.bond, n.atom, ion.atoms);
      if (!side || side.has(c)) continue;
      const beta = g.atoms[n.atom].nbrs.find(m => side.has(m.atom) && g.isC(m.atom) && g.isSp3(m.atom) && (g.atoms[m.atom].h + (ion.dH.get(m.atom) || 0)) > 0);
      if (!beta) continue;
      const atoms = new Set([...ion.atoms].filter(a => !side.has(a)));
      const dH = new Map(ion.dH); dH.set(c, (dH.get(c) || 0) + 1);
      const ni = new Species(g, atoms, { dH, bo: new Map(ion.bo), chg: new Map(ion.chg), tags: { onium: { het, c } } });
      const neu = new Species(g, side, { dH: new Map([[beta.atom, -1]]), bo: new Map([[g.bondBetween(n.atom, beta.atom), 2]]) });
      const k = { O: 1.0, N: 0.12, S: 0.2 }[hetEl] || 0.1;
      out.push(child(e, { mech: { site: { kind: 'n', a: het, label: '+' }, arrows: [FL(Mid(c, n.atom), ['mH', beta.atom, c, c]), FL(['mHi', beta.atom, c], Mid(n.atom, beta.atom))], highlightBonds: [[c, n.atom]] }, rule: 'oniumC', ruleName: 'Reacción del onio desde el carbono (pérdida de alqueno en R–CH=XH⁺)', ion: ni, neutral: neu, neutralLabel: `alqueno ${neu.formulaString}`, k,
        steps: [`El ion ${fs}⁺ (EE, ${hetEl === 'O' ? 'oxonio' : hetEl === 'N' ? 'iminio' : 'tionio'}) transfiere un H β de la cadena (${L(g, beta.atom)}) al carbono cargado ${L(g, c)} vía un complejo ion–neutro [R⁺ / CH₂=XH].`,
          `Se rompe ${B(g, c, n.atom)} y se elimina el alqueno ${neu.formulaString}; queda ${ni.formulaString}⁺ (p. ej., 2-butanol: m/z 59 → 31, CH₂=OH⁺).`],
        refs: ['ALCOH', 'MT', 'EVEN'] }));
    }
  }

  // (2c) Efecto orto en iones radicales secundarios (p. ej., aspirina 180 → 138 → 120)
  if (ion.oddElectron && !ion.formulaOverride && ion.atoms.size && e.generation <= 2) {
    for (const o of orthoEffect(g, ion.atoms, ion.dH)) {
      out.push(child(e, { ...o, k: Math.min(1.3, 0.4 * o.base) }));
    }
  }

  // (2d) Cadenas de pérdidas neutras tras reordenamientos
  const fchain = (to, lossF, lossLabel, k, txt, tags = {}, oe = false) => out.push(child(e, { rule: 'chain', ruleName: `Pérdida de ${lossLabel} desde ${fs}${ion.oddElectron ? '⁺•' : '⁺'}`,
    ion: fo(g, to, tags), neutral: null, neutralFormula: lossF, neutralLabel: lossLabel, k, oe, steps: [txt], refs: ['MT', 'EVEN'] }));
  if (ion.tags.orthoKetene) fchain({ ...f, C: f.C - 1, O: f.O - 1 }, { C: 1, O: 1 }, 'CO', 0.6, `El ion ceteno ${fs}⁺• formado por efecto orto expulsa CO (contracción del anillo), p. ej., salicilatos m/z 120 → 92.`, { lossCO2: (f.O || 0) >= 2 }, true);
  if (ion.tags.lossCO2 && f.O >= 1) fchain({ ...f, C: f.C - 1, O: f.O - 1 }, { C: 1, O: 1 }, 'CO', 0.3, `Segunda expulsión de CO desde ${fs}⁺• (m/z 92 → 64).`, {}, true);
  if (ion.tags.orthoNitro) fchain({ ...f, C: f.C - 1, O: f.O - 1 }, { C: 1, O: 1 }, 'CO', 0.5, `El ion ${fs}⁺ (o-nitrotolueno −•OH) expulsa CO: m/z 120 → 92.`, { lossHCN: true });
  if (ion.tags.lossHCN && f.N >= 1) fchain({ ...f, C: f.C - 1, N: f.N - 1, H: f.H - 1 }, { C: 1, H: 1, N: 1 }, 'HCN', 0.7, `${fs}⁺ elimina HCN dando el catión C₅H₅⁺ (m/z 65).`);
  if (ion.tags.retroImide && f.O >= 1) fchain({ ...f, C: f.C - 1, O: f.O - 1 }, { C: 1, O: 1 }, 'CO', 1.3, `El ion ${fs}⁺• (tras la retro-reacción) expulsa CO del carbonilo remanente (cafeína 137 → 109).`, { imideHCN: 2 }, true);
  if (ion.tags.imideHCN && f.N >= 1) fchain({ ...f, C: f.C - 1, N: f.N - 1, H: f.H - 1 }, { C: 1, H: 1, N: 1 }, 'HCN', 0.6, `Expulsión de HCN desde el anillo imidazólico de ${fs}⁺• (cafeína 109 → 82 → 55).`, ion.tags.imideHCN > 1 ? { imideHCN: ion.tags.imideHCN - 1 } : {}, true);
  if (ion.tags.cycloalkene) fchain({ ...f, C: f.C - 1, H: f.H - 3 }, { C: 1, H: 3 }, '•CH₃', 0.5, `El radical-catión cicloalqueno ${fs}⁺• formado por deshidratación pierde •CH₃ (ciclohexanol 82 → 67), como el ciclohexeno.`);
  if (ion.oddElectron && f.C >= 3 && f.C <= 5 && f.H === 2 * f.C && Object.keys(f).length === 2) fchain({ C: f.C, H: f.H - 1 }, { H: 1 }, 'H•', 0.4, `El radical-catión ${fs}⁺• pierde H• alílico dando el catión alilo ${f.C === 3 ? 'C₃H₅⁺ (m/z 41)' : ''}.`, { allylic: true });

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

  // (3b) Cationes alilo/propargilo CₙH₂ₙ₋₁⁺ → −H₂ (series 41→39, 55→53, 69→67; Beauchamp, alcanos §7 y alquinos §3)
  if ((ion.tags.allyl || ion.tags.allylic) && f.C >= 3 && f.C <= 5 && f.H === 2 * f.C - 1 && Object.keys(f).length === 2) {
    const n = f.C;
    out.push(child(e, { rule: 'allylH2', ruleName: 'Pérdida de H₂ en cationes alilo (→ C₃H₃⁺ ciclopropenilo y homólogos)', ion: fo(g, { C: n, H: 2 * n - 3 }), neutral: null,
      neutralFormula: { H: 2 }, neutralLabel: 'H₂', k: n === 3 ? 0.3 : 0.15,
      steps: [`El catión ${fs}⁺ (m/z ${14 * n - 1}) elimina H₂ y genera C${n}H${2 * n - 3}⁺ (m/z ${14 * n - 3}); para n = 3 se obtiene el catión ciclopropenilo aromático (m/z 39), presente en casi todos los espectros EI.`],
      refs: ['BEAU', 'MT'] }));
  }

  // (3c) Doble reordenamiento de McLafferty en cetonas con H γ en ambas cadenas (Beauchamp, carbonilos, 4-octanona)
  const mt = ion.tags.mclafferty;
  if (mt && typeof mt === 'object' && mt.kind === 'carbonyl' && !ion.formulaOverride && e.generation === 1) {
    const { x, y, al } = mt;
    for (const n2 of g.atoms[y].nbrs) {
      const al2 = n2.atom;
      if (al2 === x || al2 === al || !ion.atoms.has(al2) || !g.isC(al2) || !g.isSp3(al2)) continue;
      for (const nb of g.atoms[al2].nbrs) {
        const be = nb.atom;
        if (be === y || !ion.atoms.has(be) || !g.isC(be) || !g.isSp3(be) || g.bonds[nb.bond].inRing) continue;
        const ga = g.atoms[be].nbrs.find(m => m.atom !== al2 && ion.atoms.has(m.atom) && g.isC(m.atom) && g.isSp3(m.atom) && g.atoms[m.atom].h > 0);
        if (!ga) continue;
        const side = g.sideOf(nb.bond, be, ion.atoms);
        if (!side || side.has(y)) continue;
        const atoms = new Set([...ion.atoms].filter(a => !side.has(a)));
        const dH = new Map(ion.dH); dH.set(al2, (dH.get(al2) || 0) + 1);
        const ni = new Species(g, atoms, { dH, bo: new Map(ion.bo), oddElectron: true, depictAsNeutral: true, tags: {} });
        const neu = new Species(g, side, { dH: new Map([[ga.atom, -1]]), bo: new Map([[g.bondBetween(be, ga.atom), 2]]) });
        out.push(child(e, { rule: 'mcl2', ruleName: 'Segundo reordenamiento de McLafferty (doble McLafferty)', ion: ni, neutral: neu, neutralLabel: `alqueno ${neu.formulaString}`, k: 0.15, oe: true,
          steps: [`El ion enólico ${fs}⁺• conserva una segunda cadena con H γ (${L(g, ga.atom)}) sobre el carbono ${L(g, y)}.`,
            `Nueva transferencia 1,5 de H (estado de transición de 6 miembros) y escisión β de ${B(g, al2, be)}: se elimina ${neu.formulaString} y se forma ${ni.formulaString}⁺• (m/z ${ni.nominal}, masa par).`,
            'En cetonas dialquílicas con H γ en ambos lados se observan los dos McLafferty simples y el doble (p. ej., 4-octanona: 128 → 100/86 → 58).'],
          refs: ['BEAU', 'MCL59', 'MT'] }));
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
