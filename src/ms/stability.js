'use strict';
/**
 * Escalas heurísticas de estabilidad de cationes y radicales usadas para estimar
 * la abundancia relativa de los iones (modelo semicuantitativo).
 *
 * Fundamentos:
 *  - Regla de Stevenson–Audier: en una escisión simple, la carga queda preferentemente en el
 *    fragmento de menor energía de ionización (Stevenson, D. P. Discuss. Faraday Soc. 10 (1951) 35;
 *    Audier, H. E. Org. Mass Spectrom. 2 (1969) 283).
 *  - Estabilidad de carbocationes (afinidades por hidruro en fase gas): metilo << 1° < 2° < 3° ≈ alilo
 *    < bencilo/tropilio, y estabilización por heteroátomos n-donadores N > S > O > halógeno
 *    (iminio > tionio > oxonio). McLafferty & Tureček, "Interpretation of Mass Spectra", 4.ª ed.,
 *    University Science Books (1993), cap. 4; Gross, J. H. "Mass Spectrometry: A Textbook", 3.ª ed., Springer (2017), cap. 6.
 *  - Regla de la "pérdida del radical alquilo mayor": en α-escisiones competitivas se pierde
 *    preferentemente el radical más grande (McLafferty & Tureček 1993, §4.6).
 * Los valores numéricos son factores de ponderación calibrados con espectros EI de referencia
 * (NIST/EPA/NIH Mass Spectral Library) y NO son energías termoquímicas.
 */

const ONIUM = { N: 3.2, S: 2.3, O: 1.7, Cl: 0.5, Br: 0.5, I: 0.4, F: 0.1, P: 1.5, Se: 2.0 };

function carbonNbrs(g, atoms, i) { return g.atoms[i].nbrs.filter(n => atoms.has(n.atom) && g.atoms[n.atom].el === 'C'); }

/** Clasifica el centro catiónico en el átomo c del fragmento `atoms`. */
function cationClass(g, atoms, c) {
  const a = g.atoms[c];
  if (a.el !== 'C') {
    // catión centrado en heteroátomo (p. ej., X+ tras pérdida de R•): poco estable salvo excepciones
    return { cls: 'hetero', score: 0.25, label: `catión centrado en ${a.el}` };
  }
  if (a.aromatic) return { cls: 'aryl', score: 0.4, label: 'catión arilo (fenilo)' };
  const o = g.carbonylO(c);
  if (o != null && atoms.has(o)) {
    const het = a.nbrs.find(n => atoms.has(n.atom) && n.atom !== o && ['O', 'N'].includes(g.atoms[n.atom].el));
    if (het) return { cls: 'acylium-het', score: 0.9, label: g.atoms[het.atom].el === 'O' ? 'ion alcoxicarbonilo (+O≡C–OR)' : 'ion carbamoilo (+O≡C–NR₂)' };
    return { cls: 'acylium', score: 2.3, label: 'ion acilio (R–C≡O⁺)' };
  }
  let best = null;
  for (const n of a.nbrs) {
    if (!atoms.has(n.atom)) continue;
    const nb = g.atoms[n.atom]; const bnd = g.bonds[n.bond];
    if (bnd.order !== 1) continue;
    if (ONIUM[nb.el] && !nb.aromatic) {
      const cand = { cls: 'onium', het: n.atom, score: ONIUM[nb.el],
        label: nb.el === 'N' ? 'ion iminio (C=N⁺)' : nb.el === 'O' ? 'ion oxonio (C=O⁺)' : nb.el === 'S' ? 'ion tionio (C=S⁺)' : `ion ${nb.el}-onio` };
      if (!best || cand.score > best.score) best = cand;
    }
    if (ONIUM[nb.el] && nb.aromatic && nb.el !== 'C') {
      const cand = { cls: 'onium', het: n.atom, score: ONIUM[nb.el] * 0.6, label: 'catión estabilizado por heteroátomo aromático' };
      if (!best || cand.score > best.score) best = cand;
    }
    if (nb.el === 'C' && nb.aromatic) {
      const cand = { cls: 'benzyl', score: 2.5, label: 'catión bencilo ⇄ ion tropilio' };
      if (!best || cand.score > best.score) best = cand;
    }
    if (nb.el === 'C' && !nb.aromatic) {
      const dbl = nb.nbrs.find(m => m.atom !== c && atoms.has(m.atom) && g.bonds[m.bond].order === 2 && g.atoms[m.atom].el === 'C');
      if (dbl) {
        const cand = { cls: 'allyl', score: 1.15, label: 'catión alilo (estabilizado por resonancia)' };
        if (!best || cand.score > best.score) best = cand;
      }
    }
  }
  if (best) return best;
  const deg = carbonNbrs(g, atoms, c).length;
  return [
    { cls: 'methyl', score: 0.04, label: 'catión metilo' },
    { cls: 'primary', score: 0.4, label: 'carbocatión primario' },
    { cls: 'secondary', score: 0.95, label: 'carbocatión secundario' },
    { cls: 'tertiary', score: 1.45, label: 'carbocatión terciario' }
  ][Math.min(deg, 3)];
}

/** Estabilidad del radical neutro expulsado con el electrón desapareado en r. */
function radicalClass(g, atoms, r) {
  if (r == null) return { score: 0.12, label: 'H•' };
  const a = g.atoms[r];
  const heavy = atoms.size;
  const sizeBonus = 1 + 0.07 * Math.min(heavy - 1, 8);
  let s; let label;
  switch (a.el) {
    case 'I': s = 1.7; label = 'I•'; break;
    case 'Br': s = 1.15; label = 'Br•'; break;
    case 'Cl': s = 0.6; label = 'Cl•'; break;
    case 'F': s = 0.04; label = 'F•'; break;
    case 'O': s = heavy === 1 ? 0.18 : 0.8; label = heavy === 1 ? 'HO•' : 'radical alcoxilo/aciloxilo'; break;
    case 'N': {
      const nO = a.nbrs.filter(n => atoms.has(n.atom) && g.atoms[n.atom].el === 'O').length;
      if (nO >= 2) { s = 1.0; label = 'NO₂•'; } else { s = heavy === 1 ? 0.2 : 0.35; label = 'radical amino'; }
      break;
    }
    case 'S': s = 0.6; label = 'radical tiilo'; break;
    case 'C': {
      if (a.aromatic) { s = 0.08; label = 'radical arilo'; break; }
      const o = g.carbonylO(r);
      if (o != null && atoms.has(o)) {
        const oh = a.nbrs.find(n => n.atom !== o && atoms.has(n.atom) && g.atoms[n.atom].el === 'O');
        s = oh ? 0.6 : 0.85; label = oh ? 'radical carboxilo/alcoxicarbonilo' : 'radical acilo'; break;
      }
      if (g.tripleBondPartner(r, 'N') != null) { s = 0.3; label = '•CN'; break; }
      const cnb = a.nbrs.filter(n => atoms.has(n.atom));
      if (cnb.some(n => g.atoms[n.atom].aromatic)) { s = 1.5; label = 'radical bencilo'; break; }
      if (cnb.some(n => ['N', 'O', 'S'].includes(g.atoms[n.atom].el))) { s = 1.15; label = 'radical α-heterosustituido'; break; }
      const deg = cnb.filter(n => g.atoms[n.atom].el === 'C').length;
      s = [0.12, 1.0, 1.15, 1.3][Math.min(deg, 3)];
      label = ['•CH₃', 'radical alquilo primario', 'radical alquilo secundario', 'radical alquilo terciario'][Math.min(deg, 3)];
      break;
    }
    default: s = 0.4; label = `radical centrado en ${a.el}`;
  }
  return { score: s * sizeBonus, label };
}

module.exports = { cationClass, radicalClass, ONIUM };
