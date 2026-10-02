'use strict';
/**
 * Reglas ampliadas (v1.1):
 *   12. α-Escisión en anillos con ion distónico (cetonas, alcoholes, aminas, éteres y sulfuros cíclicos)
 *   13. Retro-reacción en imidas/ureas cíclicas (pérdida de R–N=C=O; purinas, uracilos, barbitúricos)
 *   14. Efecto orto (eliminación de ROH/H₂O y de •OH entre sustituyentes vecinos)
 *   15. Pérdida de ROH en ésteres (formación de ion ceteno)
 *
 * Fundamentos: McLafferty & Tureček (1993) §8.3–8.6 (anillos, ion distónico), §8.9 (efecto orto);
 * Gross (2017) §6.5–6.11; Schwarz, H. Top. Curr. Chem. 1978, 73, 231 (efectos orto);
 * Yates, B. F.; Bouma, W. J.; Radom, L. Tetrahedron 1986, 42, 6225 (iones distónicos).
 */
const { Species } = require('./species');
const { radicalClass } = require('./stability');
const { ev, L, B, clean } = require('./rules');

const allOf = (g) => new Set(g.atoms.map(a => a.idx));

/** Recorrido del anillo empezando en `start` y alejándose de `avoid` (vecino en el anillo). */
function ringPath(ring, start, avoid) {
  const n = ring.length; const i = ring.indexOf(start); const j = ring.indexOf(avoid);
  const dir = j === (i + 1) % n ? -1 : 1;
  const P = [];
  for (let k = 0; k < n; k++) P.push(ring[((i + dir * k) % n + n) % n]);
  return P; // P[n-1] === avoid
}
function twoSides(g, bonds, a, b, subset) {
  const comps = g.componentsWithout(bonds, subset || allOf(g));
  const ca = comps.find(c => c.has(a)); const cb = comps.find(c => c.has(b));
  if (!ca || !cb || ca === cb) return null;
  if (comps.length !== 2) return null;
  return [ca, cb];
}

/* ------------------------------------------------------------------ */
/* 12. α-Escisión en anillos (ion distónico)                           */
/* ------------------------------------------------------------------ */
function ringAlpha(g) {
  const out = []; const seen = new Set();
  const push = (key, e) => { if (seen.has(key)) return; seen.add(key); out.push(ev(e)); };
  for (const R of g.rings) {
    const n = R.length;
    if (n < 5 || n > 7 || R.some(i => g.atoms[i].aromatic)) continue;
    for (const cc of R) {
      if (!g.isC(cc)) continue;
      const sites = [];
      const o = g.carbonylO(cc);
      if (o != null && !R.includes(o)) sites.push({ type: 'a', het: o, hetBond: g.bondBetween(cc, o), order: 3 });
      if (g.isSp3(cc)) {
        for (const nb of g.atoms[cc].nbrs) {
          const X = g.atoms[nb.atom];
          if (!['O', 'N', 'S'].includes(X.el) || X.aromatic || !g.isSp3(X.idx) || g.bonds[nb.bond].order !== 1) continue;
          if (X.nbrs.some(m => m.atom !== cc && (g.atoms[m.atom].aromatic || g.carbonylO(m.atom) != null))) continue;
          sites.push({ type: R.includes(X.idx) ? 'c' : 'b', het: X.idx, hetBond: nb.bond, order: 2 });
        }
      }
      for (const site of sites) {
        const cbs = g.atoms[cc].nbrs.map(m => m.atom).filter(a => R.includes(a) && a !== site.het && g.isC(a));
        for (const cb of cbs) {
          const P = ringPath(R, cc, cb);
          if (site.type === 'c' && P[1] !== site.het) continue;
          const bCb = g.bondBetween(cc, cb);
          if (g.bonds[bCb].order !== 1) continue;
          const charge = { bo: [[site.hetBond, site.order]], chg: [[site.het, 1]] };
          const typeTxt = { a: 'cetona cíclica', b: site.het != null && g.atoms[site.het].el === 'N' ? 'amina en anillo carbocíclico' : 'alcohol/éter en anillo carbocíclico', c: g.atoms[site.het].el === 'N' ? 'amina cíclica' : g.atoms[site.het].el === 'O' ? 'éter cíclico' : 'sulfuro cíclico' }[site.type];
          const ionName = site.type === 'a' ? 'ion acilio' : g.atoms[site.het].el === 'N' ? 'ion iminio' : g.atoms[site.het].el === 'O' ? 'ion oxonio' : 'ion tionio';
          const step1 = [
            site.type === 'a' ? `Ionización en el par libre del oxígeno carbonílico ${L(g, site.het)}.` : `Ionización en el par libre de ${L(g, site.het)}.`,
            `α-Escisión del enlace del anillo ${B(g, cc, cb)} (↷): el anillo se abre sin cambio de masa y se forma un ion distónico — la carga queda como ${ionName} en ${L(g, cc)} y el radical en ${L(g, cb)}.`
          ];
          // ---- Ruta B: transferencia de H al radical + β-escisión → pérdida de radical alquilo (EE)
          const d = site.type === 'c' ? 2 : 1;
          const donor = P[d];
          if (d + 2 <= n - 1 && g.isC(donor) && g.atoms[donor].h > 0 && (n - 1 - d) >= 3) {
            const bBeta = g.bondBetween(P[d + 1], P[d + 2]);
            const bDb = g.bondBetween(P[d], P[d + 1]);
            const sides = g.bonds[bBeta].order === 1 && g.bonds[bDb].order === 1 ? twoSides(g, [bCb, bBeta], cc, cb) : null;
            if (sides) {
              const [ionAt, neuAt] = sides;
              const ion = new Species(g, ionAt, { dH: new Map([[donor, -1]]), bo: new Map([...charge.bo, [bDb, 2]]), chg: new Map(charge.chg),
                tags: site.type === 'a' ? { acylium: { c: cc, o: site.het } } : { onium: { het: site.het, c: cc } } });
              const neu = new Species(g, neuAt, { dH: new Map([[cb, 1]]), rad: new Map([[P[d + 2], 1]]) });
              const rc = radicalClass(g, neuAt, P[d + 2]);
              const ts = n - d + 1; // miembros del estado de transición (átomos de la cadena donador…radical + H)
              const base = { a: 2.4, b: 1.9, c: 0.3 }[site.type] * (ts >= 6 ? 1 : 0.8);
              push('B' + ion.key(), {
                rule: 'ringAlphaB', ruleName: `α-Escisión en anillo + transferencia de H + β-escisión (${typeTxt})`, ion, neutral: neu, neutralLabel: rc.label,
                cleaved: [bCb, bBeta], oe: false, score: base * rc.score,
                steps: [...step1,
                  `Transferencia 1,${ts - 1} de un H desde ${L(g, donor)} al radical primario en ${L(g, cb)} (estado de transición de ${ts} miembros): el radical migra a ${L(g, donor)}, posición estabilizada por el grupo cargado vecino.`,
                  `β-Escisión (↷) del enlace ${B(g, P[d + 1], P[d + 2])}: se forma el doble enlace ${L(g, P[d])}=${L(g, P[d + 1])} conjugado con el sitio de carga y se expulsa ${rc.label} (${neu.formulaString}•).`,
                  site.type === 'a' ? 'En la ciclohexanona y la ciclopentanona origina el pico base m/z 55 (CH₂=CH–C≡O⁺).' :
                    site.type === 'b' ? 'En el ciclohexanol origina el pico base m/z 57 (CH₂=CH–CH=OH⁺); en la ciclohexilamina, m/z 56.' :
                      'En la piperidina origina m/z 56 (CH₂=CH–CH=NH₂⁺).'],
                refs: ['MT', 'GROSS', 'DISTONIC']
              });
            }
          }
          // ---- Ruta A: β-escisión desde el radical inicial → pérdida de alqueno (ion distónico OE)
          const minIdx = site.type === 'c' ? 2 : 1;
          if (n - 3 >= minIdx) {
            const bA = g.bondBetween(P[n - 2], P[n - 3]);
            const bDb = g.bondBetween(P[n - 1], P[n - 2]);
            const sides = g.bonds[bA].order === 1 && g.bonds[bDb].order === 1 ? twoSides(g, [bCb, bA], cc, cb) : null;
            if (sides) {
              const [ionAt, neuAt] = sides;
              const ion = new Species(g, ionAt, { bo: new Map(charge.bo), chg: new Map(charge.chg), rad: new Map([[P[n - 3], 1]]), oddElectron: true,
                tags: site.type === 'a' ? { acylium: { c: cc, o: site.het }, distonic: true } : { distonic: true } });
              const neu = new Species(g, neuAt, { bo: new Map([[bDb, 2]]) });
              const base = { a: 1.3, b: 0.25, c: n === 5 && g.atoms[site.het].el === 'N' ? 1.3 : 0.22 }[site.type];
              push('A' + ion.key(), {
                rule: 'ringAlphaA', ruleName: `α-Escisión en anillo + β-escisión con pérdida de alqueno (${typeTxt})`, ion, neutral: neu, neutralLabel: `alqueno ${neu.formulaString}`,
                cleaved: [bCb, bA], oe: true, score: base,
                steps: [...step1,
                  `β-Escisión (↷) del enlace ${B(g, P[n - 2], P[n - 3])} iniciada por el radical en ${L(g, cb)}: se expulsa el alqueno ${neu.formulaString} (${L(g, cb)}=${L(g, P[n - 2])}).`,
                  `Queda un ion radical distónico ${ion.formulaString}⁺• (OE, masa ${ion.nominal % 2 ? 'impar' : 'par'}), con la carga en ${L(g, cc)}/${L(g, site.het)} y el radical en ${L(g, P[n - 3])}.`,
                  site.type === 'a' ? 'En cetonas cíclicas este ion pierde CO a continuación (ciclohexanona 98 → 70 → 42).' : 'En la pirrolidina origina el pico base m/z 43.'],
                refs: ['MT', 'GROSS', 'DISTONIC']
              });
            }
          }
          // ---- Ruta C (heteroátomo en el anillo): expulsión de R₂C=X neutro, carga en el resto hidrocarbonado
          if (site.type === 'c') {
            const X = site.het; const bX = g.bondBetween(X, P[2]);
            const sides = twoSides(g, [bCb, bX], cc, cb);
            if (sides) {
              const [neuAt, ionAt] = sides;
              const ion = new Species(g, ionAt, { oddElectron: true, tags: { noDepict: true, alkeneOE: true } });
              const neu = new Species(g, neuAt, { bo: new Map([[site.hetBond, 2]]) });
              const base = { O: 1.3, N: 0.12, S: 0.6 }[g.atoms[X].el];
              push('C' + ion.key(), {
                rule: 'ringAlphaC', ruleName: `Apertura del anillo con expulsión de ${neu.formulaString} (${typeTxt})`, ion, neutral: neu, neutralLabel: `neutro ${neu.formulaString}`,
                cleaved: [bCb, bX], oe: true, score: base,
                steps: [...step1,
                  `Escisión del segundo enlace ${B(g, X, P[2])} (⇒, inductiva): se expulsa ${neu.formulaString} neutro (${g.atoms[X].el === 'O' ? 'aldehído' : g.atoms[X].el === 'N' ? 'imina' : 'tioaldehído'}) y la carga queda en el fragmento hidrocarbonado ${ion.formulaString}⁺• (regla de Stevenson: EI(CH₂O) = 10,9 eV > EI(C₃H₆)).`,
                  'En el tetrahidrofurano origina el pico base m/z 42.'],
                refs: ['MT', 'STEV']
              });
            }
          }
        }
      }
    }
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* 13. Retro-reacción de imidas/ureas cíclicas: pérdida de R–N=C=O      */
/* ------------------------------------------------------------------ */
function isocyanateLoss(g) {
  const best = new Map();
  for (const R of g.rings) {
    if (R.length !== 6) continue;
    for (let k = 0; k < 6; k++) {
      for (const dir of [1, -1]) {
        const nI = R[k]; const cI = R[(k + dir + 6) % 6];
        if (g.atoms[nI].el !== 'N' || !g.isC(cI)) continue;
        const o = g.carbonylO(cI); if (o == null || R.includes(o)) continue;
        const a = R[(k - dir + 6) % 6]; const b = R[(k + 2 * dir + 6) % 6];
        const bNa = g.bondBetween(nI, a); const bCb = g.bondBetween(cI, b);
        const sides = twoSides(g, [bNa, bCb], nI, a);
        if (!sides) continue;
        const [neuAt, ionAt] = sides;
        if (neuAt.size > 6) continue;
        const ion = new Species(g, ionAt, { oddElectron: true, tags: { noDepict: true, retroImide: true } });
        const neu = new Species(g, neuAt, { bo: new Map([[g.bondBetween(nI, cI), 2]]) });
        const fs = ion.formulaString;
        const e = ev({
          rule: 'isocyanate', ruleName: 'Retro-reacción en imida/urea cíclica: pérdida de isocianato R–N=C=O', ion, neutral: neu,
          neutralLabel: `isocianato ${neu.formulaString}`, cleaved: [bNa, bCb], oe: true, score: 1.3,
          steps: [
            `Ionización en el sistema n/π del anillo de pirimidinadiona (o urea/imida cíclica).`,
            `Ruptura concertada (formalmente retro-[2+2+2]/retro-Diels–Alder) de los enlaces ${B(g, nI, a)} y ${B(g, cI, b)}.`,
            `Se expulsa el isocianato neutro ${neu.formulaString} (unidad ${L(g, nI)}–${L(g, cI)}=O) y queda el ion radical ${fs}⁺•.`,
            'Es la fragmentación característica de xantinas (cafeína 194 → 137), uracilos y barbitúricos; le siguen pérdidas consecutivas de CO y HCN.'
          ],
          refs: ['MT', 'GROSS']
        });
        if (!best.has(fs) || best.get(fs).score < e.score) best.set(fs, e);
      }
    }
  }
  return [...best.values()];
}

/* ------------------------------------------------------------------ */
/* 14. Efecto orto                                                     */
/* ------------------------------------------------------------------ */
/**
 * Busca pares de sustituyentes en posición orto (donador de H / grupo aceptor) dentro del
 * subconjunto `atoms` (con ajustes dH). Devuelve eventos con su puntaje base.
 */
function orthoEffect(g, atoms = allOf(g), dH = new Map()) {
  const out = []; const seen = new Set();
  const H = (i) => g.atoms[i].h + (dH.get(i) || 0);
  for (const c1 of atoms) {
    const A1 = g.atoms[c1]; if (!A1.aromatic || A1.el !== 'C') continue;
    for (const dn of A1.nbrs) {
      const D = dn.atom; const Da = g.atoms[D];
      if (!atoms.has(D) || Da.aromatic || g.bonds[dn.bond].order !== 1 || H(D) <= 0) continue;
      const dKind = ['O', 'N', 'S'].includes(Da.el) ? 'het' : (Da.el === 'C' && g.isSp3(D)) ? 'C' : null;
      if (!dKind) continue;
      for (const rn of A1.nbrs) {
        const c2 = rn.atom; if (!atoms.has(c2) || !g.atoms[c2].aromatic || c2 === D) continue;
        if (!g.rings.some(r => r.includes(c1) && r.includes(c2))) continue;
        for (const an of g.atoms[c2].nbrs) {
          const Ac = an.atom; const Aa = g.atoms[Ac];
          if (!atoms.has(Ac) || Aa.aromatic) continue;
          // (i) éster / ácido orto: eliminación de ROH (H₂O en ácidos)
          const oc = g.carbonylO(Ac);
          if (Aa.el === 'C' && oc != null && atoms.has(oc)) {
            const oe = Aa.nbrs.find(m => m.atom !== oc && m.atom !== c2 && g.atoms[m.atom].el === 'O' && atoms.has(m.atom) && g.bonds[m.bond].order === 1);
            if (!oe) continue;
            const side = g.sideOf(oe.bond, oe.atom, atoms);
            if (!side || side.has(D) || side.size > 8) continue;
            const ionAt = new Set([...atoms].filter(a => !side.has(a)));
            const ndH = new Map(dH); ndH.set(D, (ndH.get(D) || 0) - 1);
            const ion = new Species(g, ionAt, { dH: ndH, oddElectron: true, tags: { noDepict: true, orthoKetene: true } });
            const nH = new Map(); for (const a of side) if (dH.has(a)) nH.set(a, dH.get(a)); nH.set(oe.atom, (nH.get(oe.atom) || 0) + 1);
            const neu = new Species(g, side, { dH: nH });
            const key = 'i' + ion.formulaString; if (seen.has(key)) continue; seen.add(key);
            const lost = neu.formulaString === 'H2O' ? 'H₂O' : neu.formulaString === 'CH4O' ? 'CH₃OH' : neu.formulaString;
            out.push({ rule: 'orthoROH', ruleName: `Efecto orto: eliminación de ${lost}`, ion, neutral: neu, neutralLabel: lost, oe: true, cleaved: [oe.bond],
              base: dKind === 'het' ? (side.size === 1 ? 2.6 : 4.5) : 1.6,
              steps: [
                `Ionización en el anillo aromático o en el grupo ${dKind === 'het' ? Da.el + '–H' : 'alquilo'} orto.`,
                `El H de ${L(g, D)} (sustituyente en ${L(g, c1)}) se transfiere al oxígeno ${L(g, oe.atom)} del grupo carboxílico vecino en ${L(g, c2)} mediante un estado de transición de seis miembros, posible sólo por la proximidad orto.`,
                `Ruptura de ${B(g, Ac, oe.atom)} y eliminación de ${lost} neutro: se forma un ion radical tipo ceteno/o-quinona metiduro ${ion.formulaString}⁺• (OE).`,
                'Los isómeros meta y para no pueden dar esta eliminación; el ion [M−ROH]⁺• permite asignar la sustitución orto (p. ej., salicilato de metilo 152 → 120, pico base).'
              ], refs: ['ORTHO', 'MT'] });
          }
          // (ii) nitro orto a CH / OH / NH: pérdida de •OH
          if (Aa.el === 'N') {
            const os = Aa.nbrs.filter(m => g.atoms[m.atom].el === 'O' && atoms.has(m.atom));
            if (os.length < 2) continue;
            const ionAt = new Set(atoms); ionAt.delete(os[0].atom);
            const ndH = new Map(dH); ndH.set(D, (ndH.get(D) || 0) - 1);
            const ion = new Species(g, ionAt, { dH: ndH, tags: { noDepict: true, orthoNitro: true } });
            const key = 'ii' + ion.formulaString; if (seen.has(key)) continue; seen.add(key);
            out.push({ rule: 'orthoOH', ruleName: 'Efecto orto en nitroarenos: pérdida de •OH', ion, neutral: null, neutralFormula: { O: 1, H: 1 }, neutralLabel: '•OH', oe: false, cleaved: [],
              base: dKind === 'C' ? 3.2 : 0.6,
              steps: [
                `Ionización en el grupo nitro (${L(g, Ac)}).`,
                `Un oxígeno del NO₂ abstrae un H de ${L(g, D)} (sustituyente orto en ${L(g, c1)}) a través de un estado de transición de seis miembros.`,
                `Expulsión de •OH ([M−17]⁺) y formación de un ion par-electrónico cíclico (tipo antranilo/benzisoxazolio) ${ion.formulaString}⁺.`,
                'Diagnóstico de o-nitrotoluenos (m/z 120, pico base); ausente en los isómeros meta y para.'
              ], refs: ['ORTHO', 'MT'] });
          }
        }
      }
    }
  }
  return out;
}
function orthoPrimary(g) {
  return orthoEffect(g).map(e => ev({ ...e, score: e.base }));
}

/* ------------------------------------------------------------------ */
/* 15. Pérdida de ROH en ésteres (ion ceteno)                           */
/* ------------------------------------------------------------------ */
function esterROH(g) {
  const out = [];
  const all = allOf(g);
  for (const C of g.atoms) {
    if (C.el !== 'C' || C.aromatic) continue;
    const o = g.carbonylO(C.idx); if (o == null) continue;
    const oe = C.nbrs.find(n => n.atom !== o && g.atoms[n.atom].el === 'O' && g.bonds[n.bond].order === 1 && g.atoms[n.atom].nbrs.length === 2);
    const ca = C.nbrs.find(n => n.atom !== o && g.isC(n.atom) && g.isSp3(n.atom) && g.atoms[n.atom].h > 0 && !g.bonds[n.bond].inRing);
    if (!oe || !ca || g.bonds[oe.bond].inRing) continue;
    const side = g.sideOf(oe.bond, oe.atom, all); if (!side || side.has(C.idx)) continue;
    const ionAt = new Set([...all].filter(a => !side.has(a)));
    const ion = new Species(g, ionAt, { dH: new Map([[ca.atom, -1]]), bo: new Map([[ca.bond, 2]]), oddElectron: true, depictAsNeutral: true });
    const neu = new Species(g, side, { dH: new Map([[oe.atom, 1]]) });
    out.push(ev({ rule: 'esterROH', ruleName: `Pérdida de ${neu.formulaString} (alcohol) en ésteres`, ion, neutral: neu, neutralLabel: `alcohol ${neu.formulaString}`,
      oe: true, cleaved: [oe.bond], score: 0.12,
      steps: [`Ionización en el oxígeno carbonílico ${L(g, o)}.`,
        `Transferencia de un H α (${L(g, ca.atom)}) al oxígeno alcoxílico ${L(g, oe.atom)} (estado de transición de 4 miembros o vía enol).`,
        `Eliminación del alcohol ${neu.formulaString} y formación del radical-catión ceteno R–CH=C=O⁺• ([M−ROH]⁺•, normalmente débil salvo en ésteres orto-sustituidos).`],
      refs: ['MT'] }));
  }
  return out;
}

const EXT_RULES = [ringAlpha, isocyanateLoss, orthoPrimary, esterROH];
module.exports = { EXT_RULES, orthoEffect };
