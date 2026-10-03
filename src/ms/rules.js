'use strict';
/**
 * Reglas de fragmentación unimolecular de iones M+• en EI (70 eV).
 * Cada regla recorre el grafo molecular y genera "eventos" de fragmentación con:
 *   ion (Species), neutro (Species), puntaje heurístico, pasos mecanísticos (texto) y referencias.
 *
 * Convenciones de flechas (McLafferty & Tureček 1993, cap. 4):
 *   ↷ (anzuelo, "fishhook") = desplazamiento de un electrón (escisión homolítica / iniciada por el sitio radical, "α")
 *   ⇒ (flecha completa)      = desplazamiento de un par de electrones (escisión heterolítica / inductiva, "i")
 */
const { Species } = require('./species');
const { cationClass, radicalClass, ONIUM } = require('./stability');

const L = (g, i) => (g.label ? g.label(i) : `${g.atoms[i].el}${i}`); // etiqueta de átomo (índice del dibujo original)
const B = (g, i, j) => `${L(g, i)}–${L(g, j)}`;

function allAtoms(g) { return new Set(g.atoms.map(a => a.idx)); }
function complement(g, set) { const s = new Set(); for (const a of g.atoms) if (!set.has(a.idx)) s.add(a.idx); return s; }

/** Escisión de un enlace acíclico: devuelve [ladoA, ladoB] o null si el enlace está en un anillo. */
function split(g, bondIdx, keepAtom) {
  const all = allAtoms(g);
  const side = g.sideOf(bondIdx, keepAtom, all);
  if (side.size === all.size) return null;
  return [side, complement(g, side)];
}

function ev(o) {
  return Object.assign({ score: 0, steps: [], refs: [], children: [], generation: 1 }, o);
}

/* ------------------------------------------------------------------ */
/* 1. α-Escisión iniciada por el sitio radical en heteroátomos sp3     */
/* ------------------------------------------------------------------ */
function alphaHetero(g) {
  const out = [];
  for (const X of g.atoms) {
    if (!['N', 'O', 'S', 'Cl', 'Br', 'I'].includes(X.el) || X.aromatic) continue;
    if (!g.isSp3(X.idx)) continue;
    const xi = X.idx;
    const isHal = ['Cl', 'Br', 'I'].includes(X.el);
    // atenuación cuando el heteroátomo está conjugado (éster, amida, fenol/anisol)
    const conj = X.nbrs.some(n => g.atoms[n.atom].aromatic || g.carbonylO(n.atom) != null);
    const fConj = conj ? 0.35 : 1.0;
    const fHal = isHal ? 0.15 : 1.0;
    for (const na of X.nbrs) {
      const ca = na.atom;
      if (!g.isC(ca) || g.atoms[ca].aromatic || !g.isSp3(ca)) continue;
      const kind = { N: 'amina', O: X.h ? 'alcohol' : 'éter', S: X.h ? 'tiol' : 'sulfuro', Cl: 'cloruro', Br: 'bromuro', I: 'yoduro' }[X.el];
      // (a) pérdida de un radical R• unido a Cα
      for (const nr of g.atoms[ca].nbrs) {
        if (nr.atom === xi) continue;
        if (g.bonds[nr.bond].inRing) continue;
        const parts = split(g, nr.bond, ca); if (!parts) continue;
        const [ionAt, neuAt] = parts;
        const xc = g.bondBetween(xi, ca);
        const ion = new Species(g, ionAt, { bo: new Map([[xc, 2]]), chg: new Map([[xi, 1]]),
          tags: { onium: { het: xi, c: ca } } });
        const neu = new Species(g, neuAt, { rad: new Map([[nr.atom, 1]]) });
        const rc = radicalClass(g, neuAt, nr.atom);
        // Regla de Stevenson: si el radical saliente es bencílico/alílico (EI baja), la carga tiende a quedarse en él
        const fStev = (X.el === 'O' || X.el === 'S') && /bencilo/.test(rc.label) ? 0.3 : 1;
        const ionS = ONIUM[X.el] * fConj * fHal * fStev;
        out.push(ev({
          rule: 'alpha', ruleName: `α-Escisión (sitio radical en ${X.el}, ${kind})`, ion, neutral: neu,
          neutralLabel: rc.label, cleaved: [nr.bond], oe: false,
          score: ionS * rc.score,
          steps: [
            `Ionización: se remueve un electrón no enlazante del ${X.el} (${L(g, xi)}), el sitio de menor energía de ionización → M⁺• con carga y radical localizados en ${L(g, xi)}.`,
            `El electrón desapareado de ${L(g, xi)} se aparea con un electrón del enlace σ ${B(g, ca, nr.atom)} (↷) formando el enlace π ${L(g, xi)}=${L(g, ca)}.`,
            `El otro electrón del enlace ${B(g, ca, nr.atom)} migra a ${L(g, nr.atom)} (↷): se expulsa ${rc.label} y queda el ${X.el === 'N' ? 'ion iminio' : X.el === 'O' ? 'ion oxonio' : X.el === 'S' ? 'ion tionio' : 'ion halonio'} par-electrónico, estabilizado por el par libre del heteroátomo.`,
            `Regla de Stevenson / pérdida del radical mayor: entre α-escisiones competitivas se favorece expulsar el radical alquilo más grande y estable.`
          ],
          refs: X.el === 'N' ? ['MT', 'AMINES', 'GROSS'] : X.el === 'O' ? ['MT', 'ALCOH', 'GROSS'] : isHal ? ['MT', 'HALIDES'] : ['MT', 'GROSS']
        }));
      }
      // (b) pérdida de H• desde Cα  (M − 1)
      const ringX = g.bonds[na.bond].inRing;
      if (g.atoms[ca].h > 0 && !isHal) {
        const xc = g.bondBetween(xi, ca);
        const ion = new Species(g, allAtoms(g), { dH: new Map([[ca, -1]]), bo: new Map([[xc, 2]]), chg: new Map([[xi, 1]]),
          tags: { onium: { het: xi, c: ca } } });
        out.push(ev({
          rule: 'alphaH', ruleName: `α-Escisión con pérdida de H• (${X.el})`, ion,
          neutral: null, neutralFormula: { H: 1 }, neutralLabel: 'H•', cleaved: [], oe: false,
          score: ONIUM[X.el] * fConj * 0.02 * (ringX ? 8 : 1),
          steps: [
            `Ionización en el par libre de ${L(g, xi)} (M⁺•).`,
            ...(ringX ? [`En heterociclos saturados la α-escisión de enlaces del anillo no cambia la masa (abre el anillo); por eso la pérdida de H• desde ${L(g, ca)} genera un ion iminio/oxonio cíclico [M−1]⁺ intenso (p. ej., piperidina m/z 84).`] : []),
            `Escisión homolítica (↷) del enlace C–H en ${L(g, ca)}: el radical forma el enlace π ${L(g, xi)}=${L(g, ca)} y se expulsa H•, generando [M−H]⁺.`,
            'La pérdida de H• es poco favorable frente a radicales alquilo (H• es el radical menos estable), por lo que [M−1]⁺ suele ser de baja intensidad salvo en aminas y aldehídos.'
          ],
          refs: ['MT', 'GROSS']
        }));
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 2. Grupo carbonilo: α-escisión (ion acilio) y escisión inductiva   */
/* ------------------------------------------------------------------ */
function carbonyl(g) {
  const out = [];
  for (const C of g.atoms) {
    if (C.el !== 'C' || C.aromatic) continue;
    const o = g.carbonylO(C.idx); if (o == null) continue;
    const ci = C.idx; const co = g.bondBetween(ci, o);
    const others = C.nbrs.filter(n => n.atom !== o);
    const fam = classifyCarbonyl(g, ci, o);
    // aldehído: pérdida de H•
    if (C.h > 0) {
      const ion = new Species(g, allAtoms(g), { dH: new Map([[ci, -1]]), bo: new Map([[co, 3]]), chg: new Map([[o, 1]]),
        tags: { acylium: { c: ci, o } } });
      const aromaticNb = others.some(n => g.atoms[n.atom].aromatic);
      out.push(ev({
        rule: 'alphaCO_H', ruleName: 'α-Escisión en aldehído (pérdida de H•)', ion, neutral: null, neutralFormula: { H: 1 },
        neutralLabel: 'H•', cleaved: [], oe: false, score: 2.3 * (aromaticNb ? 0.9 : 0.15),
        steps: [
          `Ionización por remoción de un electrón n del oxígeno carbonílico ${L(g, o)}.`,
          `El radical en ${L(g, o)} forma un segundo enlace π con ${L(g, ci)} (↷) y se rompe homolíticamente el enlace C–H aldehídico: se expulsa H• y se forma el ion acilio R–C≡O⁺ ([M−1]⁺).`,
          aromaticNb ? 'En aldehídos aromáticos el ion [M−1]⁺ es muy intenso por conjugación del acilio con el anillo.' : 'En aldehídos alifáticos [M−1]⁺ es débil; también aparece HCO⁺ (m/z 29).'
        ],
        refs: ['MT', 'SILV']
      }));
    }
    for (const nr of others) {
      const r = nr.atom; if (g.bonds[nr.bond].inRing) continue;
      const parts = split(g, nr.bond, ci); if (!parts) continue;
      const [ionAt, neuAt] = parts;
      // (a) α-escisión → ion acilio
      const ion = new Species(g, ionAt, { bo: new Map([[co, 3]]), chg: new Map([[o, 1]]), tags: { acylium: { c: ci, o } } });
      const neu = new Species(g, neuAt, { rad: new Map([[r, 1]]) });
      const cc = cationClass(g, ionAt, ci); const rc = radicalClass(g, neuAt, r);
      const fHCO = ionAt.size === 2 && g.atoms[ci].h === 1 ? 0.25 : 1; // HCO⁺: EI(HCO•) alta frente a R•
      out.push(ev({
        rule: 'alphaCO', ruleName: `α-Escisión de ${fam} (formación de ion acilio)`, ion, neutral: neu, neutralLabel: rc.label,
        cleaved: [nr.bond], oe: false, score: fHCO * cc.score * rc.score * ([...ionAt].some(i => g.atoms[i].aromatic && g.bondBetween(i, ci) != null) ? 1.5 : 1),
        steps: [
          `Ionización: remoción de un electrón no enlazante (n) del oxígeno carbonílico ${L(g, o)} → M⁺• (radical-catión localizado en O).`,
          `El electrón desapareado del oxígeno forma un nuevo enlace con ${L(g, ci)} (↷), mientras el enlace σ ${B(g, ci, r)} se rompe homolíticamente (↷).`,
          `Se expulsa ${rc.label} (${neu.formulaString}) y se obtiene el ${cc.label}, estabilizado por resonancia R–C⁺=O ↔ R–C≡O⁺.`,
          'Es la fragmentación dominante de cetonas, ésteres, ácidos y amidas; la competencia entre los dos sustituyentes sigue la regla de Stevenson y la estabilidad del radical expulsado.'
        ],
        refs: ['MT', 'GROSS', 'SILV']
      }));
      // (b) escisión inductiva: la carga migra a R (R⁺ + •C(=O)R')
      if (g.isC(r)) {
        const ion2 = new Species(g, neuAt, { chg: new Map([[r, 1]]), tags: alkylTag(g, neuAt) });
        const neu2 = new Species(g, ionAt, { rad: new Map([[ci, 1]]) });
        const cc2 = cationClass(g, neuAt, r);
        out.push(ev({
          rule: 'iCO', ruleName: `Escisión inductiva (i) en ${fam}`, ion: ion2, neutral: neu2, neutralLabel: 'radical acilo',
          cleaved: [nr.bond], oe: false, score: cc2.score * 0.85 * 0.45,
          steps: [
            `Ionización en el oxígeno carbonílico ${L(g, o)}.`,
            `El sitio cargado atrae (efecto inductivo, ⇒) el par de electrones del enlace ${B(g, ci, r)}: ruptura heterolítica.`,
            `La carga positiva queda en ${L(g, r)} como ${cc2.label}; se pierde un radical acilo neutro (${neu2.formulaString}•).`,
            'Compite con la α-escisión; es relevante cuando R⁺ es un catión estabilizado (3°, bencílico, alílico).'
          ],
          refs: ['MT', 'GROSS']
        }));
      }
    }
  }
  return out;
}
function classifyCarbonyl(g, ci, o) {
  const nb = g.atoms[ci].nbrs.filter(n => n.atom !== o).map(n => g.atoms[n.atom]);
  const hasO = nb.some(a => a.el === 'O'); const hasN = nb.some(a => a.el === 'N');
  if (nb.some(a => ['F', 'Cl', 'Br', 'I'].includes(a.el))) return 'haluro de ácido';
  if (hasO) return nb.some(a => a.el === 'O' && a.h > 0) ? 'ácido carboxílico' : 'éster';
  if (hasN) return 'amida';
  if (g.atoms[ci].h > 0) return 'aldehído';
  return 'cetona';
}

/* ------------------------------------------------------------------ */
/* 3. Escisión inductiva C–X (haluros, éteres, nitro)                  */
/* ------------------------------------------------------------------ */
function inductive(g) {
  const out = [];
  for (const b of g.bonds) {
    if (b.order !== 1 || b.inRing) continue;
    for (const [c, x] of [[b.a, b.b], [b.b, b.a]]) {
      if (!g.isC(c)) continue;
      const X = g.atoms[x];
      if (g.carbonylO(c) != null) continue;
      let f = 0; let kind = '';
      if (['Cl', 'Br', 'I', 'F'].includes(X.el)) { f = 2.5; kind = 'haluro'; }
      else if (X.el === 'N' && X.nbrs.filter(n => g.atoms[n.atom].el === 'O').length >= 2) { f = 4.0; kind = 'nitrocompuesto'; }
      else if (X.el === 'O' && !g.atoms[c].aromatic && g.isSp3(c)) { f = 0.3; kind = X.h ? 'alcohol' : 'éter'; }
      else if (X.el === 'S' && !g.atoms[c].aromatic && g.isSp3(c)) { f = 0.45; kind = 'tioéter/tiol'; }
      else if (X.el === 'N' && !g.atoms[c].aromatic && g.isSp3(c)) { f = 0.12; kind = 'amina'; }
      if (!f) continue;
      const parts = split(g, b.idx, c); if (!parts) continue;
      const [ionAt, neuAt] = parts;
      const cc = cationClass(g, ionAt, c); const rc = radicalClass(g, neuAt, x);
      const ion = new Species(g, ionAt, { chg: new Map([[c, 1]]), tags: alkylTag(g, ionAt) });
      const neu = new Species(g, neuAt, { rad: new Map([[x, 1]]) });
      out.push(ev({
        rule: 'inductive', ruleName: `Escisión inductiva (i) del enlace C–${X.el} (${kind})`, ion, neutral: neu, neutralLabel: rc.label,
        cleaved: [b.idx], oe: false, score: f * cc.score * rc.score,
        steps: [
          `Ionización en el heteroátomo ${L(g, x)} (electrón n), sitio de carga del M⁺•.`,
          `El heteroátomo electronegativo atrae el par de electrones del enlace ${B(g, c, x)} (⇒, escisión heterolítica inducida por la carga).`,
          `Se forma ${cc.label} en ${L(g, c)} y se expulsa ${rc.label}.`,
          kind === 'haluro' ? 'La tendencia a perder X• sigue I > Br > Cl >> F, en orden inverso a la fuerza del enlace C–X.' :
            kind === 'nitrocompuesto' ? 'En nitroarenos la pérdida de •NO₂ ([M−46]⁺) es característica.' :
              'La escisión inductiva es más competitiva cuando el carbocatión resultante es secundario, terciario o bencílico.'
        ],
        refs: kind === 'haluro' ? ['MT', 'HALIDES', 'GROSS'] : ['MT', 'GROSS']
      }));
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 4. Iones onio cíclicos de 5 y 3 miembros (Cl, Br, I, S)             */
/*    McLafferty 1962; Beauchamp, MS chapter: halógenos §5, tioles §4   */
/* ------------------------------------------------------------------ */
const CYCLO_ONIUM = {
  5: { Br: 1.9, Cl: 1.4, I: 0.2, S: 0.9 },
  3: { Br: 0.12, Cl: 0.1, I: 0.05, S: 0.55 }
};
function haloniumCyclic(g) {
  const out = []; const seen = new Set();
  for (const X of g.atoms) {
    if (!['Cl', 'Br', 'I', 'S'].includes(X.el) || X.aromatic) continue;
    if (X.el !== 'S' && X.nbrs.length !== 1) continue;
    if (X.el === 'S' && (!g.isSp3(X.idx) || X.nbrs.length > 2 || X.nbrs.some(n => !g.isC(n.atom)))) continue;
    for (const nx of X.nbrs) {
      const c1 = nx.atom;
      if (!g.isC(c1) || !g.isSp3(c1) || g.atoms[c1].aromatic || g.atoms[c1].h < 2 || g.bonds[nx.bond].inRing) continue;
      const paths = [];
      const walk = (path) => {
        if (path.length === 3 || path.length === 5) paths.push(path.slice());
        if (path.length === 5) return;
        const last = path[path.length - 1];
        for (const n of g.atoms[last].nbrs) {
          if (path.includes(n.atom) || n.atom === X.idx) continue;
          if (!g.isC(n.atom) || !g.isSp3(n.atom) || g.bonds[n.bond].inRing) continue;
          path.push(n.atom); walk(path); path.pop();
        }
      };
      walk([c1]);
      for (const p of paths) {
        const ring = p.length === 5 ? 5 : 3;
        const cr = p[p.length - 2]; const cx = p[p.length - 1];  // cr cierra el anillo con X; se rompe cr–cx
        const bi = g.bondBetween(cr, cx);
        const parts = split(g, bi, cr); if (!parts) continue;
        const [ionAt, neuAt] = parts;
        if (!ionAt.has(X.idx)) continue;
        const key = ring + ':' + [...ionAt].sort((a, b) => a - b).join(',');
        if (seen.has(key)) continue; seen.add(key);
        const rc = radicalClass(g, neuAt, cx);
        const ion = new Species(g, ionAt, { extraBonds: [{ a: X.idx, b: cr, order: 1 }], chg: new Map([[X.idx, 1]]), tags: { cyclicOnium: ring } });
        const neu = new Species(g, neuAt, { rad: new Map([[cx, 1]]) });
        const nm = { Br: 'bromonio', Cl: 'cloronio', I: 'yodonio', S: 'sulfonio' }[X.el];
        const prod = ring === 5 ? (X.el === 'S' ? 'tiolanio (tetrahidrotiofenio)' : 'tetrahidrohalogenonio') : (X.el === 'S' ? 'tiiranio' : 'halogenonio de 3 miembros (tipo etilenhalonio)');
        out.push(ev({
          rule: ring === 5 ? 'halonium' : 'onium3', ruleName: `Ciclación con desplazamiento: ion ${nm} cíclico de ${ring} miembros`,
          ion, neutral: neu, neutralLabel: rc.label, cleaved: [bi], oe: false,
          score: CYCLO_ONIUM[ring][X.el] * rc.score,
          steps: [
            `Ionización en el ${X.el === 'S' ? 'azufre' : 'halógeno'} ${L(g, X.idx)} (electrón n).`,
            ring === 5
              ? `El radical-catión de ${L(g, X.idx)} ataca intramolecularmente a ${L(g, cr)} (estado de transición de 5 miembros, desplazamiento rC).`
              : `El par libre de ${L(g, X.idx)} puentea al carbono β ${L(g, cr)} (escisión β asistida, anillo de 3 miembros).`,
            `Se rompe homolíticamente ${B(g, cr, cx)}, se expulsa ${rc.label} y se forma el ion ${prod} ${ion.formulaString}⁺ (m/z ${ion.nominal}).`,
            ring === 5
              ? 'Pico característico de cadenas de ≥ C₅ con Cl, Br, I o S: m/z 91/93 (Cl), 135/137 (Br), 183 (I), 89 + R (S).'
              : 'Ion diagnóstico menor: m/z 63/65 (Cl), 107/109 (Br), 155 (I), 61 + R (S).'
          ],
          refs: X.el === 'S' ? ['MT', 'BEAU'] : ['HALIDES', 'MT', 'BEAU']
        }));
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 5. Escisión bencílica (→ ion tropilio) y 6. alílica               */
/* ------------------------------------------------------------------ */
function benzylicAllylic(g) {
  const out = [];
  for (const cb of g.atoms) {
    if (cb.el !== 'C' || cb.aromatic || !g.isSp3(cb.idx)) continue;
    const arNb = cb.nbrs.find(n => g.atoms[n.atom].aromatic && g.isC(n.atom) && !g.bonds[n.bond].inRing);
    let allylNb = null;
    if (!arNb) {
      for (const n of cb.nbrs) {
        if (!g.isC(n.atom) || g.atoms[n.atom].aromatic) continue;
        const d = g.atoms[n.atom].nbrs.find(m => m.atom !== cb.idx && g.bonds[m.bond].order === 2 && g.isC(m.atom) && !g.bonds[m.bond].aromatic)
          || g.atoms[n.atom].nbrs.find(m => m.atom !== cb.idx && g.bonds[m.bond].order === 3 && g.isC(m.atom));
        if (d) { allylNb = { atom: n.atom, dbl: d.atom, triple: g.bonds[g.bondBetween(n.atom, d.atom)].order === 3 }; break; }
      }
    }
    const anchor = arNb ? arNb.atom : allylNb ? allylNb.atom : null;
    if (anchor == null) continue;
    const isBenz = !!arNb;
    const isProp = !isBenz && allylNb.triple;
    for (const nr of cb.nbrs) {
      if (nr.atom === anchor || g.bonds[nr.bond].inRing) continue;
      const parts = split(g, nr.bond, cb.idx); if (!parts) continue;
      const [ionAt, neuAt] = parts;
      const cc = cationClass(g, ionAt, cb.idx); const rc = radicalClass(g, neuAt, nr.atom);
      const ion = new Species(g, ionAt, { chg: new Map([[cb.idx, 1]]), tags: { benzyl: isBenz, allyl: !isBenz && !isProp, propargyl: isProp } });
      const neu = new Species(g, neuAt, { rad: new Map([[nr.atom, 1]]) });
      const isTrop = isBenz && ion.formulaString === 'C7H7';
      out.push(ev({
        rule: isBenz ? 'benzylic' : isProp ? 'propargylic' : 'allylic',
        ruleName: isBenz ? 'Escisión bencílica (β al anillo) → ion bencilo/tropilio' : isProp ? 'Escisión propargílica (→ C₃H₃⁺ y homólogos)' : 'Escisión alílica',
        ion, neutral: neu, neutralLabel: rc.label, cleaved: [nr.bond], oe: false,
        score: isBenz ? 2.4 * Math.max(cc.score, 2.5) * Math.pow(rc.score, 0.3) : isProp ? Math.max(cc.score, 0.9) * 0.6 * Math.pow(rc.score, 0.3) : Math.max(cc.score, 1.15) * 0.8 * rc.score,
        steps: isBenz ? [
          `Ionización: remoción de un electrón π del anillo aromático (sistema de menor EI) → M⁺• deslocalizado.`,
          `Escisión del enlace ${B(g, cb.idx, nr.atom)}, β respecto al anillo (posición bencílica): el catión bencílico en ${L(g, cb.idx)} queda estabilizado por resonancia con el anillo; se expulsa ${rc.label}.`,
          isTrop ? 'El catión bencilo C₇H₇⁺ (m/z 91) se isomeriza al ion tropilio aromático (6 electrones π, regla de Hückel), responsable de la gran intensidad de m/z 91 en alquilbencenos.'
            : 'El catión bencílico sustituido puede expandirse a un ion tropilio sustituido.'
        ] : isProp ? [
          `Ionización: remoción de un electrón π del triple enlace C≡C (${L(g, anchor)}≡${L(g, allylNb.dbl)}).`,
          `Escisión del enlace propargílico ${B(g, cb.idx, nr.atom)} (β al triple enlace): se forma el catión propargilo (HC≡C–CH₂⁺ ↔ H₂C=C=CH⁺) y se expulsa ${rc.label}.`,
          'La forma alenilo coloca la carga sobre un carbono sp (más electronegativo), por lo que la estabilización es menor que en el alilo; el ion C₃H₃⁺ (m/z 39) puede reordenarse al catión ciclopropenilo aromático (2 e⁻ π). Serie homóloga 39, 53, 67… (Δ = 14, CH₂).'
        ] : [
          `Ionización: remoción de un electrón π del doble enlace C=C (${L(g, anchor)}=${L(g, allylNb.dbl)}).`,
          `Escisión del enlace alílico ${B(g, cb.idx, nr.atom)}: se forma un catión alilo deslocalizado (C=C–C⁺ ↔ ⁺C–C=C) y se expulsa ${rc.label}.`,
          'La migración del doble enlace en el ion radical puede generar varios cationes alílicos isómeros; los alquenos dan series CₙH₂ₙ₋₁⁺ (m/z 41, 55, 69…).'
        ],
        refs: isBenz ? ['MT', 'TROP', 'GROSS'] : ['MT', 'GROSS', 'BEAU']
      }));
    }
    // pérdida de H• bencílico (p. ej., tolueno → m/z 91)
    if (isBenz && cb.h > 0) {
      const ion = new Species(g, allAtoms(g), { dH: new Map([[cb.idx, -1]]), chg: new Map([[cb.idx, 1]]), tags: { benzyl: true } });
      out.push(ev({
        rule: 'benzylicH', ruleName: 'Pérdida de H• bencílico → ion tropilio', ion, neutral: null, neutralFormula: { H: 1 }, neutralLabel: 'H•',
        cleaved: [], oe: false, score: (cb.nbrs.length === 1 ? 2.2 : 0.1) * (ewgOnRing(g, arNb.atom) ? 0.12 : 1),
        steps: [
          'Ionización π del anillo aromático.',
          `Ruptura del enlace C–H bencílico en ${L(g, cb.idx)} con formación del catión bencilo que se reordena a tropilio (C₇H₆R⁺).`,
          'En el tolueno este proceso origina el pico base m/z 91.'
        ],
        refs: ['TROP', 'MT']
      }));
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 7. Escisión σ de enlaces C–C (alcanos y cadenas alquílicas)         */
/* ------------------------------------------------------------------ */
const ALKYL_SIZE_W = { 1: 0.08, 2: 0.45, 3: 1.0, 4: 1.0, 5: 0.65, 6: 0.45, 7: 0.35, 8: 0.3 };
function isSaturatedAlkyl(f) {
  const k = Object.keys(f).filter(x => f[x]);
  return k.every(x => x === 'C' || x === 'H') && f.H === 2 * f.C + 1;
}
function alkylTag(g, atoms) {
  const f = g.formulaOf(atoms);
  return isSaturatedAlkyl(f) ? { alkyl: f.C } : {};
}
function sigmaCC(g) {
  const out = [];
  const hetero = g.atoms.some(a => a.el !== 'C');
  const arom = g.atoms.some(a => a.aromatic);
  const base = hetero || arom ? 0.25 : 1.2;
  for (const b of g.bonds) {
    if (b.order !== 1 || b.inRing || b.aromatic) continue;
    if (!g.isC(b.a) || !g.isC(b.b) || !g.isSp3(b.a) || !g.isSp3(b.b)) continue;
    for (const [c, r] of [[b.a, b.b], [b.b, b.a]]) {
      const parts = split(g, b.idx, c); if (!parts) continue;
      const [ionAt, neuAt] = parts;
      const f = g.formulaOf(ionAt);
      const cc = cationClass(g, ionAt, c); const rc = radicalClass(g, neuAt, r);
      const sat = isSaturatedAlkyl(f);
      const w = sat ? (ALKYL_SIZE_W[f.C] || 0.25) : 0.5;
      const ion = new Species(g, ionAt, { chg: new Map([[c, 1]]), tags: sat ? { alkyl: f.C } : {} });
      const neu = new Species(g, neuAt, { rad: new Map([[r, 1]]) });
      out.push(ev({
        rule: 'sigma', ruleName: 'Escisión σ de enlace C–C', ion, neutral: neu, neutralLabel: rc.label, cleaved: [b.idx], oe: false,
        score: base * cc.score * rc.score * w,
        steps: [
          `Ionización por remoción de un electrón σ (enlaces C–C, en ausencia de electrones n o π de menor EI): la carga y el radical quedan deslocalizados en el esqueleto.`,
          `Ruptura del enlace ${B(g, c, r)} (σ → ${cc.label} + ${rc.label}).`,
          'La escisión se favorece en los puntos de ramificación (formación de carbocationes 2° y 3°); los alcanos lineales generan la serie CₙH₂ₙ₊₁⁺ (m/z 29, 43, 57, 71…) con máximo en C₃–C₄.'
        ],
        refs: ['MT', 'GROSS', 'SILV']
      }));
    }
  }
  return out;
}

/**
 * Energía de ionización vertical aproximada (eV) de un alqueno/aleno neutro, para repartir la carga
 * entre los dos fragmentos de un McLafferty de alqueno según la regla de Stevenson.
 * Valores de referencia (NIST WebBook): eteno 10.51, propeno 9.73, 1-buteno 9.55, 1-penteno 9.49,
 * 1-hexeno 9.44; cada sustituyente alquilo adicional sobre los C sp² la reduce ≈ 0.3–0.4 eV
 * (E-2-buteno 9.10, isobuteno 9.22, 2-metil-2-buteno 8.69). Aleno 9.69.
 */
function alkeneIE(g, atoms, a, b) {
  const nC = [...atoms].filter(i => g.isC(i)).length;
  const base = { 2: 10.51, 3: 9.73, 4: 9.55, 5: 9.49, 6: 9.44 }[nC] ?? 9.40;
  const subs = [a, b].reduce((s, c) => s + g.atoms[c].nbrs.filter(n => atoms.has(n.atom) && n.atom !== a && n.atom !== b && g.isC(n.atom)).length, 0);
  return base - 0.35 * Math.max(0, subs - (nC >= 3 ? 1 : 0));
}
/** Fracción de la carga retenida en el fragmento 1 (Stevenson, reparto tipo Boltzmann con kT efectivo 0,25 eV). */
const stevensonShare = (ie1, ie2) => 1 / (1 + Math.exp((ie1 - ie2) / 0.25));

/* ------------------------------------------------------------------ */
/* 8. Reordenamiento de McLafferty (transferencia γ-H, escisión β)     */
/* ------------------------------------------------------------------ */
function mclafferty(g) {
  const out = []; const seen = new Set();
  const acceptors = [];
  for (const b of g.bonds) {
    const A = g.atoms[b.a], Bt = g.atoms[b.b];
    if (b.aromatic) continue;
    if (b.order === 2 && ((A.el === 'O' && Bt.el === 'C') || (A.el === 'C' && Bt.el === 'O'))) {
      const [x, y] = A.el === 'O' ? [b.a, b.b] : [b.b, b.a];
      acceptors.push({ x, y, kind: 'carbonyl', base: 1.8, bond: b.idx });
    } else if (b.order === 2 && ((A.el === 'N' && Bt.el === 'C') || (A.el === 'C' && Bt.el === 'N'))) {
      const [x, y] = A.el === 'N' ? [b.a, b.b] : [b.b, b.a];
      acceptors.push({ x, y, kind: 'imina', base: 2.0, bond: b.idx });
    } else if (b.order === 3 && ((A.el === 'N' && Bt.el === 'C') || (A.el === 'C' && Bt.el === 'N'))) {
      const [x, y] = A.el === 'N' ? [b.a, b.b] : [b.b, b.a];
      acceptors.push({ x, y, kind: 'nitrilo', base: 1.3, bond: b.idx });
    } else if (b.order === 2 && A.el === 'C' && Bt.el === 'C') {
      acceptors.push({ x: b.a, y: b.b, kind: 'alqueno', base: 2.2, bond: b.idx });
      acceptors.push({ x: b.b, y: b.a, kind: 'alqueno', base: 2.2, bond: b.idx });
    } else if (b.order === 3 && A.el === 'C' && Bt.el === 'C') {
      // alquinos: McLafferty "tipo alqueno" con formación de un radical-catión aleno de masa par (Beauchamp, alquinos §5)
      acceptors.push({ x: b.a, y: b.b, kind: 'alquino', base: 0.75, bond: b.idx });
      acceptors.push({ x: b.b, y: b.a, kind: 'alquino', base: 0.75, bond: b.idx });
    }
  }
  // aceptor aromático: Y = C ipso, X = C orto
  for (const a of g.atoms) {
    if (!a.aromatic || a.el !== 'C') continue;
    const exo = a.nbrs.find(n => !g.atoms[n.atom].aromatic && g.isC(n.atom));
    if (!exo) continue;
    for (const n of a.nbrs) if (g.atoms[n.atom].aromatic && g.atoms[n.atom].h > 0 || (g.atoms[n.atom].aromatic && g.atoms[n.atom].el === 'C'))
      acceptors.push({ x: n.atom, y: a.idx, kind: 'aromático', base: 1.5, bond: n.bond, aromatic: true });
  }
  for (const acc of acceptors) {
    const { x, y } = acc;
    for (const na of g.atoms[y].nbrs) {
      const al = na.atom; if (al === x) continue;
      if (g.bonds[na.bond].order !== 1 || g.atoms[al].aromatic) continue;
      if (!['C', 'O', 'N'].includes(g.atoms[al].el)) continue;
      if (g.atoms[al].el === 'C' && !g.isSp3(al)) continue;
      for (const nb of g.atoms[al].nbrs) {
        const be = nb.atom; if (be === y) continue;
        if (!g.isC(be) || !g.isSp3(be) || g.bonds[nb.bond].inRing) continue;
        for (const ng of g.atoms[be].nbrs) {
          const ga = ng.atom; if (ga === al) continue;
          const gaHet = ['N', 'O'].includes(g.atoms[ga].el) && g.atoms[ga].h > 0 && g.isSp3(ga);
          if (gaHet && !(acc.aromatic || acc.kind === 'alqueno')) continue;
          if (!gaHet && (!g.isC(ga) || g.atoms[ga].h === 0 || g.atoms[ga].aromatic)) continue;
          if (!g.isSp3(ga)) continue;
          const parts = split(g, nb.bond, al); if (!parts) continue;
          const [ionAt, neuAt] = parts;
          if (!ionAt.has(x)) continue;
          const key = acc.kind + ':' + [...ionAt].sort().join(',');
          if (seen.has(key)) continue; seen.add(key);
          const bo = new Map();
          if (!acc.aromatic) {
            bo.set(acc.bond, acc.kind === 'nitrilo' || acc.kind === 'alquino' ? 2 : 1);
            bo.set(na.bond, 2);
          }
          const ion = new Species(g, ionAt, { dH: new Map([[x, 1]]), bo, oddElectron: true, depictAsNeutral: true,
            tags: { mclafferty: { x, y, al, kind: acc.kind }, noDepict: !!acc.aromatic, alkeneIon: acc.kind === 'alqueno' && !gaHet } });
          const bg = g.bondBetween(be, ga);
          const neu = new Species(g, neuAt, { dH: new Map([[ga, -1]]), bo: new Map([[bg, 2]]) });
          const famTxt = acc.kind === 'carbonyl' ? classifyCarbonyl(g, y, x) : acc.kind;
          // Éster con el H γ en la cadena alcoxílica (α = O del éster)
          const alkoxySide = acc.kind === 'carbonyl' && g.atoms[al].el === 'O';
          const acylConj = alkoxySide && g.atoms[y].nbrs.some(n => n.atom !== al && n.atom !== x && (g.atoms[n.atom].aromatic || g.doubleBondPartner(n.atom) != null));
          let score = acc.base;
          if (alkoxySide) score = acylConj ? 1.0 : 0.15;
          if (gaHet) score = 2.2;
          // McLafferty de alqueno/alquino: ambos fragmentos son alquenos (o aleno) de EI parecida → reparto de carga por Stevenson
          // (Beauchamp, alquenos §4: 1-hepteno → 56 [C₄H₈⁺•, 100 %] > 42 [C₃H₆⁺•, 55 %])
          let shareIon = null;
          if ((acc.kind === 'alqueno' || acc.kind === 'alquino') && !gaHet) {
            const ieIon = acc.kind === 'alquino' ? 9.69 - 0.3 * Math.max(0, [...ionAt].filter(i => g.isC(i)).length - 3) : alkeneIE(g, ionAt, y, al);
            const ieNeu = alkeneIE(g, neuAt, be, ga);
            shareIon = stevensonShare(ieIon, ieNeu);
          }
          const total = score;
          if (shareIon != null) score = total * shareIon;
          if (alkoxySide && g.atoms[be].h > 0) {
            // McLafferty + 1: doble transferencia de H → ácido protonado RC(OH)2⁺ (EE)
            const ionP = new Species(g, ionAt, { dH: new Map([[x, 1], [al, 1]]), bo: new Map([[acc.bond, 1]]), chg: new Map([[y, 1]]), tags: { mclPlus1: true } });
            const neuP = new Species(g, neuAt, { dH: new Map([[ga, -1], [be, -1]]), bo: new Map([[bg0(g, be, ga), 2]]), rad: new Map([[be, 1]]) });
            out.push(ev({
              rule: 'mclPlus1', ruleName: 'Reordenamiento de doble H («McLafferty + 1») en ésteres', ion: ionP, neutral: neuP,
              neutralLabel: `radical alquenilo ${neuP.formulaString}•`, cleaved: [nb.bond], oe: false,
              score: acylConj ? 0.45 : 0.4,
              steps: [
                `Ionización en el oxígeno carbonílico ${L(g, x)}.`,
                `Primera transferencia: el H γ de ${L(g, ga)} pasa a ${L(g, x)} por un estado de transición de seis miembros (como en el McLafferty).`,
                `Segunda transferencia: un H de ${L(g, be)} migra al oxígeno alcoxílico ${L(g, al)} (vía complejo ion–neutro) mientras se rompe ${B(g, al, be)}.`,
                `Se expulsa el radical alquenilo ${neuP.formulaString}• y se forma el ácido protonado R–C(OH)₂⁺ (ion par-electrónico, masa impar: [RCOOH + H]⁺; m/z 61 en acetatos, 123 en benzoatos).`,
                'Es característico de ésteres de alcoholes con ≥ 2 carbonos y suele superar al McLafferty simple en ésteres alifáticos.'
              ],
              refs: ['MT', 'GROSS', 'SILV']
            }));
          }
          if (gaHet) ion.tags.benzylHtransfer = true;
          const ionTxt = acc.kind === 'carbonyl' ? 'radical-catión enólico' : acc.kind === 'aromático' ? 'radical-catión metilenciclohexadieno (isotolueno)' : acc.kind === 'alquino' ? 'radical-catión aleno' : 'radical-catión par-másico';
          out.push(ev({
            rule: 'mclafferty', ruleName: gaHet ? `Escisión bencílica/alílica con transferencia de H desde ${g.atoms[ga].el} (McLafferty ${famTxt})` : `Reordenamiento de McLafferty (${famTxt}${alkoxySide ? ', lado alcoxilo' : ''})`,
            ion, neutral: neu, neutralLabel: gaHet ? `${g.atoms[ga].el === 'O' ? 'aldehído/cetona' : 'imina'} ${neu.formulaString}` : `alqueno ${neu.formulaString}`,
            cleaved: [nb.bond], oe: true, score,
            steps: [
              `Ionización en ${acc.kind === 'carbonyl' ? `el par libre del oxígeno ${L(g, x)}` : acc.kind === 'aromático' ? 'el sistema π aromático' : `el sistema π ${L(g, x)}=${L(g, y)}`}.`,
              `Estado de transición cíclico de seis miembros: ${L(g, x)}=${L(g, y)}–${L(g, al)}–${L(g, be)}–${L(g, ga)}–H.`,
              `Transferencia 1,5 del hidrógeno γ (en ${L(g, ga)}) al sitio radical ${L(g, x)} (↷).`,
              `Escisión β del enlace ${B(g, al, be)} (↷): se elimina la molécula neutra ${neu.formulaString} (alqueno) y se forma un ${ionTxt} de masa par (ion impar-electrónico, OE⁺•).`,
              'Los iones OE⁺• de masa par (sin N) son diagnósticos de reordenamientos; requiere un H en posición γ accesible.',
              ...(shareIon != null ? [`Ambos productos son alquenos: la carga se reparte según su energía de ionización (regla de Stevenson); fracción estimada en ${ion.formulaString}⁺•: ${(100 * shareIon).toFixed(0)} %.`] : [])
            ],
            refs: ['MCL59', 'MT', 'GROSS']
          }));
          // carga retenida en el alqueno (minoritaria)
          const ion2 = neu.clone({ oddElectron: true, depictAsNeutral: true, tags: { mclaffertyAlkene: true, alkeneIon: !gaHet } });
          const neu2 = ion.clone();
          out.push(ev({
            rule: 'mclaffertyAlk', ruleName: 'McLafferty con retención de carga en el alqueno', ion: ion2, neutral: neu2,
            neutralLabel: `${acc.kind === 'carbonyl' ? 'enol' : 'neutro'} ${neu2.formulaString}`, cleaved: [nb.bond], oe: true,
            score: shareIon != null ? total * (1 - shareIon) : alkoxySide ? (neuAt.size >= 3 ? 0.45 * Math.min(neuAt.size, 6) / 2 : 0.01) : score * 0.08 * Math.min(neuAt.size, 6) / 2,
            steps: [
              'Mismo estado de transición de seis miembros que el reordenamiento de McLafferty.',
              ...(shareIon != null ? [`Fragmentos alqueno/alqueno: la carga se reparte por la regla de Stevenson; fracción estimada en ${ion2.formulaString}⁺•: ${(100 * (1 - shareIon)).toFixed(0)} % (p. ej., 1-hepteno: C₄H₈⁺• m/z 56 = pico base, C₃H₆⁺• m/z 42 ≈ 55 %).`] : []),
              `La carga queda en el fragmento alqueno ${ion2.formulaString}⁺• cuando su energía de ionización es comparable o menor que la del enol (regla de Stevenson).`
            ],
            refs: shareIon != null ? ['STEV', 'MT', 'BEAU'] : ['STEV', 'MT']
          }));
        }
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 9. Retro-Diels–Alder (ciclohexenos)                                 */
/* ------------------------------------------------------------------ */
function retroDielsAlder(g) {
  const out = [];
  g.rings.forEach(ring => {
    if (ring.length !== 6) return;
    const n = ring.length;
    for (let k = 0; k < n; k++) {
      const a1 = ring[k], a2 = ring[(k + 1) % n];
      const b12 = g.bond(a1, a2);
      const isDouble = b12 && ((b12.order === 2 && !b12.aromatic) || b12.aromatic);
      if (!isDouble) continue;
      const a3 = ring[(k + 2) % n], a4 = ring[(k + 3) % n], a5 = ring[(k + 4) % n], a6 = ring[(k + 5) % n];
      const restAll = [a3, a4, a5, a6];
      if (restAll.some(i => g.atoms[i].aromatic)) continue;
      // el resto del anillo debe ser saturado (salvo el enlace a1=a2)
      const ringBonds = [[a2, a3], [a3, a4], [a4, a5], [a5, a6], [a6, a1]].map(([i, j]) => g.bond(i, j));
      if (ringBonds.some(bd => !bd || bd.order !== 1)) continue;
      if (b12.aromatic && !(g.atoms[a1].aromatic && g.atoms[a2].aromatic)) continue;
      const b34 = g.bondBetween(a3, a4), b56 = g.bondBetween(a5, a6);
      const comps = g.componentsWithout([b34, b56]);
      if (comps.length !== 2) continue;
      const diene = comps.find(c => c.has(a1)); const ene = comps.find(c => c.has(a4));
      if (!diene || !ene || diene === ene) continue;
      const boD = new Map(); if (!b12.aromatic) boD.set(b12.idx, 1);
      boD.set(g.bondBetween(a2, a3), 2); boD.set(g.bondBetween(a6, a1), 2);
      const dieneSp = new Species(g, diene, { bo: b12.aromatic ? new Map() : boD, oddElectron: true, depictAsNeutral: true, tags: { rda: 'diene', noDepict: b12.aromatic } });
      const eneSp = new Species(g, ene, { bo: new Map([[g.bondBetween(a4, a5), 2]]), oddElectron: true, depictAsNeutral: true, tags: { rda: 'ene' } });
      const dieneShare = ene.size <= 2 && [...ene].every(i => g.atoms[i].el === 'C') ? 0.95 : Math.min(0.85, 0.6 + 0.05 * (diene.size - ene.size));
      const steps = (ionName, neuName) => [
        `Ionización: remoción de un electrón π del doble enlace ${L(g, a1)}=${L(g, a2)} del anillo de ciclohexeno.`,
        `Reacción retro-Diels–Alder (cicloreversión [4+2] formal, por pasos en el ion radical): se rompen los enlaces alílicos ${B(g, a3, a4)} y ${B(g, a5, a6)}.`,
        `Productos: ${ionName} (ion radical, OE⁺•) + ${neuName} neutro. La carga se reparte según la regla de Stevenson (normalmente en el fragmento dieno).`
      ];
      out.push(ev({ rule: 'rda', ruleName: 'Retro-Diels–Alder (carga en el dieno)', ion: dieneSp, neutral: eneSp, neutralLabel: `dienófilo ${eneSp.formulaString}`,
        cleaved: [b34, b56], oe: true, score: 1.6 * dieneShare, steps: steps(`dieno ${dieneSp.formulaString}`, `dienófilo ${eneSp.formulaString}`), refs: ['RDA', 'MT', 'GROSS'] }));
      out.push(ev({ rule: 'rda', ruleName: 'Retro-Diels–Alder (carga en el dienófilo)', ion: eneSp.clone(), neutral: dieneSp.clone(), neutralLabel: `dieno ${dieneSp.formulaString}`,
        cleaved: [b34, b56], oe: true, score: 1.6 * (1 - dieneShare), steps: steps(`dienófilo ${eneSp.formulaString}`, `dieno ${dieneSp.formulaString}`), refs: ['RDA', 'MT'] }));
    }
  });
  return out;
}

/* ------------------------------------------------------------------ */
/* 10. Eliminaciones de moléculas neutras pequeñas desde M+•           */
/* ------------------------------------------------------------------ */
function neutralLosses(g) {
  const out = [];
  const M = allAtoms(g);
  const onceRules = new Set(['lossHCN', 'phenolCO', 'phenolHCO', 'lossCH2O', 'lossNO']); const used = new Set();
  const mk = (o) => { if (onceRules.has(o.rule)) { if (used.has(o.rule)) return; used.add(o.rule); } out.push(ev(Object.assign({ oe: true, cleaved: [] }, o))); };
  const minus = (atoms) => { const s = new Set(M); atoms.forEach(a => s.delete(a)); return s; };

  for (const X of g.atoms) {
    // H2O de alcoholes / H2S de tioles (eliminación 1,4 / 1,2)
    if ((X.el === 'O' || X.el === 'S') && X.h === 1 && X.nbrs.length === 1) {
      const ca = X.nbrs[0].atom;
      if (!g.atoms[ca].aromatic && g.isSp3(ca)) {
        const cb = g.atoms[ca].nbrs.find(n => n.atom !== X.idx && g.isC(n.atom) && g.atoms[n.atom].h > 0 && g.isSp3(n.atom));
        if (cb) {
          const nC = g.atoms.filter(a => a.el === 'C').length;
          const ion = new Species(g, minus([X.idx]), { dH: new Map([[cb.atom, -1]]), bo: new Map([[g.bondBetween(ca, cb.atom), 2]]), oddElectron: true, depictAsNeutral: true, tags: { cycloalkene: g.inRing(ca), alkeneIon: !g.inRing(ca) } });
          mk({ rule: 'lossH2O', ruleName: X.el === 'O' ? 'Eliminación de H₂O (alcohol)' : 'Eliminación de H₂S (tiol)', ion, neutral: null,
            neutralFormula: X.el === 'O' ? { H: 2, O: 1 } : { H: 2, S: 1 }, neutralLabel: X.el === 'O' ? 'H₂O' : 'H₂S',
            score: (X.el === 'O' ? 0.7 : 0.4) * (nC >= 4 ? 1.0 : 0.5) * (g.inRing(ca) ? 2.5 : [1, 1, 0.3, 0.12][Math.min(3, g.atoms[ca].nbrs.filter(n => g.isC(n.atom)).length)]),
            steps: [
              `Ionización en el par libre de ${L(g, X.idx)}.`,
              'Transferencia de H desde un carbono δ (1,4) o β (1,2) al heteroátomo a través de un estado de transición cíclico (preferentemente de 6 miembros).',
              `Eliminación de ${X.el === 'O' ? 'H₂O ([M−18]⁺•)' : 'H₂S ([M−34]⁺•)'} y formación del radical-catión alqueno/cicloalcano.`,
              'Con frecuencia el M⁺• de alcoholes es débil o ausente y [M−18]⁺• es el ion de mayor masa observable.'
            ], refs: ['ALCOH', 'MT'] });
        }
      }
    }
    // HX de haluros de alquilo
    if (['Cl', 'Br', 'F'].includes(X.el) && X.nbrs.length === 1) {
      const ca = X.nbrs[0].atom;
      if (!g.atoms[ca].aromatic && g.isSp3(ca)) {
        const cb = g.atoms[ca].nbrs.find(n => n.atom !== X.idx && g.isC(n.atom) && g.atoms[n.atom].h > 0 && g.isSp3(n.atom));
        if (cb) {
          const ion = new Species(g, minus([X.idx]), { dH: new Map([[cb.atom, -1]]), bo: new Map([[g.bondBetween(ca, cb.atom), 2]]), oddElectron: true, depictAsNeutral: true, tags: { alkeneIon: true } });
          mk({ rule: 'lossHX', ruleName: `Eliminación de H${X.el}`, ion, neutral: null, neutralFormula: { H: 1, [X.el]: 1 }, neutralLabel: `H${X.el}`,
            score: { Cl: 0.45, Br: 0.12, F: 0.7 }[X.el],
            steps: [`Ionización en el halógeno ${L(g, X.idx)}.`, `Abstracción de un H vecinal y eliminación de H${X.el} neutro, dando el radical-catión del alqueno ([M−H${X.el}]⁺•).`],
            refs: ['HALIDES', 'MT'] });
        }
      }
    }
    // Fenoles: −CO y −HCO•
    if (X.el === 'O' && X.h === 1 && X.nbrs.length === 1 && g.atoms[X.nbrs[0].atom].aromatic) {
      const f = g.formula;
      const fCO = { ...f, C: f.C - 1, O: f.O - 1 };
      mk({ rule: 'phenolCO', ruleName: 'Pérdida de CO en fenoles', ion: new Species(g, M, { formula: clean(fCO), oddElectron: true }), neutral: null,
        neutralFormula: { C: 1, O: 1 }, neutralLabel: 'CO', score: 0.55 * (g.atoms.filter(a => a.el !== 'C' && a.el !== 'H').length > 1 ? 0.3 : 1),
        steps: ['Ionización π del anillo fenólico.', 'Tautomerización del ion radical a la forma ciclohexadienona (ceto).', 'Contracción del anillo con expulsión de CO ([M−28]⁺•) → radical-catión ciclopentadieno (C₅H₆⁺• en fenol, m/z 66).'],
        refs: ['MT', 'SILV'] });
      const fHCO = { ...fCO, H: fCO.H - 1 };
      mk({ rule: 'phenolHCO', ruleName: 'Pérdida de HCO• en fenoles', ion: new Species(g, M, { formula: clean(fHCO) }), neutral: null, oe: false,
        neutralFormula: { C: 1, H: 1, O: 1 }, neutralLabel: 'HCO•', score: 0.45,
        steps: ['Tras la pérdida de CO, el ion ciclopentadieno expulsa H• dando el catión ciclopentadienilo (C₅H₅⁺, m/z 65 en fenol); globalmente [M−29]⁺.'],
        refs: ['MT'] });
    }
    // Anilinas y piridinas: −HCN
    if (X.el === 'N' && ((X.aromatic && X.nbrs.length === 2 && X.h === 0) || (!X.aromatic && X.h === 2 && X.nbrs.length === 1 && g.atoms[X.nbrs[0].atom].aromatic))) {
      const f = g.formula; const fx = clean({ ...f, C: f.C - 1, N: f.N - 1, H: f.H - 1 });
      if (fx.C > 0 && fx.H >= 0) mk({ rule: 'lossHCN', ruleName: 'Pérdida de HCN (aminas aromáticas / azinas)', ion: new Species(g, M, { formula: fx, oddElectron: true }), neutral: null,
        neutralFormula: { H: 1, C: 1, N: 1 }, neutralLabel: 'HCN', score: (X.aromatic ? X.rings.some(ri => g.rings[ri].length === 6) && g.rings.length === 1 : true) ? 0.9 : 0.25,
        steps: ['Ionización π/n del sistema aromático nitrogenado.', 'Apertura/contracción del anillo con expulsión de HCN neutro ([M−27]⁺•), característica de anilinas, piridinas e indoles.'],
        refs: ['MT', 'SILV'] });
    }
    // Nitroarenos: −NO (ion fenoxilo)
    if (X.el === 'N' && X.nbrs.filter(n => g.atoms[n.atom].el === 'O').length >= 2) {
      const ar = X.nbrs.find(n => g.atoms[n.atom].aromatic);
      if (ar) {
        const f = g.formula; const fx = clean({ ...f, N: f.N - 1, O: f.O - 1 });
        mk({ rule: 'lossNO', ruleName: 'Reordenamiento nitro–nitrito y pérdida de NO•', ion: new Species(g, M, { formula: fx }), neutral: null, oe: false,
          neutralFormula: { N: 1, O: 1 }, neutralLabel: 'NO•', score: 0.7,
          steps: ['Ionización en el grupo nitro.', 'Isomerización Ar–NO₂⁺• → Ar–O–N=O⁺• (nitro–nitrito).', 'Expulsión de NO• ([M−30]⁺) con formación del ion arilóxido ArO⁺.'],
          refs: ['MT', 'GROSS'] });
      }
    }
    // Pérdida de ceteno en acetanilidas / acetatos de arilo (y enoles acetilados)
    if ((X.el === 'N' || X.el === 'O') && X.nbrs.some(n => g.atoms[n.atom].aromatic)) {
      for (const n of X.nbrs) {
        const c = n.atom; const o = g.carbonylO(c); if (o == null) continue;
        const me = g.atoms[c].nbrs.find(m => m.atom !== o && m.atom !== X.idx && g.isC(m.atom) && g.atoms[m.atom].h > 0 && g.isSp3(m.atom));
        if (!me) continue;
        const acyl = g.sideOf(n.bond, c, M);
        const ion = new Species(g, minus([...acyl]), { dH: new Map([[X.idx, 1]]), oddElectron: true, depictAsNeutral: true });
        const neu = new Species(g, acyl, { dH: new Map([[me.atom, -1]]), bo: new Map([[g.bondBetween(c, me.atom), 2]]) });
        mk({ rule: 'ketene', ruleName: 'Eliminación de ceteno (CH₂=C=O)', ion, neutral: neu, neutralLabel: `ceteno ${neu.formulaString}`, score: 3.5,
          steps: [`Ionización en ${L(g, X.idx)} o en el anillo aromático.`, `Transferencia de H desde ${L(g, me.atom)} hacia ${L(g, X.idx)} a través de un estado de transición de 4 miembros (o vía el anillo, 6 miembros).`,
            `Ruptura de ${B(g, X.idx, c)} y expulsión de ceteno neutro; se forma el radical-catión de la anilina/fenol ([M−42]⁺•).`],
          refs: ['MT', 'SILV'] });
      }
    }
    // Aril metil éteres: −CH2O
    if (X.el === 'O' && X.nbrs.length === 2) {
      const ar = X.nbrs.find(n => g.atoms[n.atom].aromatic); const me = X.nbrs.find(n => g.isC(n.atom) && g.atoms[n.atom].h === 3);
      if (ar && me) {
        const f = g.formula; const fx = clean({ ...f, C: f.C - 1, O: f.O - 1, H: f.H - 2 });
        mk({ rule: 'lossCH2O', ruleName: 'Pérdida de formaldehído en anisoles', ion: new Species(g, M, { formula: fx, oddElectron: true }), neutral: null,
          neutralFormula: { C: 1, H: 2, O: 1 }, neutralLabel: 'CH₂O', score: 0.7,
          steps: ['Ionización en el anillo o en el O del éter.', 'Transferencia de H del metilo al anillo (estado de transición de 4 miembros) y eliminación de CH₂O ([M−30]⁺•).'],
          refs: ['MT'] });
        const ion = new Species(g, minus([me.atom]), { formula: clean({ ...f, C: f.C - 1, H: f.H - 3 }), tags: { aryloxy: true } });
        mk({ rule: 'arylOMe', ruleName: 'Pérdida de •CH₃ en aril metil éteres', ion, neutral: null, oe: false, neutralFormula: { C: 1, H: 3 }, neutralLabel: '•CH₃', score: 0.8,
          steps: ['Ionización en el O del éter (conjugado con el anillo).', 'Ruptura homolítica O–CH₃ con formación del ion fenoxilo/oxociclohexadienilo ArO⁺ ([M−15]⁺).'],
          refs: ['MT'] });
      }
    }
  }
  return out;
}
function bg0(g, a, b) { return g.bondBetween(a, b); }
/** ¿Hay un sustituyente atractor (NO₂, C=O, C≡N) en el anillo aromático que contiene `ar`? */
function ewgOnRing(g, ar) {
  const ring = g.rings.find(r => r.includes(ar)); if (!ring) return false;
  return ring.some(i => g.atoms[i].nbrs.some(n => {
    const a = g.atoms[n.atom]; if (a.aromatic) return false;
    if (a.el === 'N' && a.nbrs.filter(m => g.atoms[m.atom].el === 'O').length >= 2) return true;
    if (a.el === 'C' && (g.carbonylO(n.atom) != null || g.tripleBondPartner(n.atom, 'N') != null)) return true;
    return false;
  }));
}
function clean(f) { const o = {}; for (const [k, v] of Object.entries(f)) if (v > 0) o[k] = v; return o; }

/* ------------------------------------------------------------------ */
/* 11. Cicloalcanos/cicloalquenos: apertura del anillo                  */
/* ------------------------------------------------------------------ */
function carbocycleOpening(g) {
  const out = [];
  const done = new Set();
  for (const ring of g.rings) {
    if (ring.length < 5 || ring.length > 7) continue;
    if (ring.some(i => g.atoms[i].aromatic || g.atoms[i].el !== 'C')) continue;
    const fused = ring.some(i => g.atoms[i].rings.length > 1);
    if (fused) continue;
    // anillos con heteroátomo exocíclico (cetonas, alcoholes, aminas cíclicas) → regla ringAlpha
    if (ring.some(i => g.atoms[i].nbrs.some(n => !g.isC(n.atom)))) continue;
    const hasDb = ring.some((i, k) => { const b = g.bond(i, ring[(k + 1) % ring.length]); return b && b.order === 2; });
    const f = g.formula; const key = hasDb ? 'ene' : 'ane'; if (done.has(key)) continue; done.add(key);
    const ringTxt = `${ring.length} miembros`;
    out.push({ rule: 'ringCH3', ruleName: 'Apertura del anillo y pérdida de •CH₃ (cicloalcano/cicloalqueno)', ion: new Species(g, new Set(), { formula: clean({ ...f, C: f.C - 1, H: f.H - 3 }) }),
      neutral: null, neutralFormula: { C: 1, H: 3 }, neutralLabel: '•CH₃', oe: false, cleaved: [], generation: 1, score: hasDb ? 2.2 : 0.5,
      steps: [`Ionización ${hasDb ? 'π del doble enlace' : 'σ del anillo'} (anillo de ${ringTxt}).`,
        `${hasDb ? 'Escisión alílica' : 'Escisión σ'} de un enlace del anillo: se forma un ion distónico (radical y carga separados) sin cambio de masa.`,
        'Transferencia de hidrógeno (1,4/1,5) que genera un radical metilo terminal y escisión que expulsa •CH₃ ([M−15]⁺); en el ciclohexeno origina el pico base m/z 67.'],
      refs: ['MT', 'GROSS'] });
    out.push({ rule: 'ringC2H4', ruleName: 'Apertura del anillo y pérdida de etileno', ion: new Species(g, new Set(), { formula: clean({ ...f, C: f.C - 2, H: f.H - 4 }), oddElectron: true }),
      neutral: null, neutralFormula: { C: 2, H: 4 }, neutralLabel: 'C₂H₄', oe: true, cleaved: [], generation: 1, score: hasDb ? 0.2 : 2.6,
      steps: ['Ionización σ del anillo y apertura a un ion radical distónico.', 'Segunda escisión con expulsión de etileno neutro ([M−28]⁺•); en el ciclohexano es el pico base m/z 56.'],
      refs: ['MT', 'GROSS'] });
  }
  return out;
}

const PRIMARY_RULES = [carbocycleOpening, alphaHetero, carbonyl, inductive, haloniumCyclic, benzylicAllylic, sigmaCC, mclafferty, retroDielsAlder, neutralLosses];

module.exports = { PRIMARY_RULES, isSaturatedAlkyl, clean, L, B, ev, split, allAtoms, complement, classifyCarbonyl, ewgOnRing, benzylicAllylic, mclafferty };
