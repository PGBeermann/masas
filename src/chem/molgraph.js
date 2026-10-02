'use strict';
/**
 * Grafo molecular construido a partir de RDKit.js (RDKit MinimalLib, WASM).
 * RDKit se encarga de: parseo/validación SMILES, sanitización, aromaticidad,
 * percepción de anillos (SSSR) y canonización. Este módulo expone un grafo
 * de átomos pesados con H implícitos, sobre el cual operan las reglas de fragmentación.
 */
const { Z_TO_SYMBOL, formulaToString, monoMass, nominalMass, rdb } = require('./elements');

class MolGraph {
  constructor(json, rdkitMol) {
    const doc = typeof json === 'string' ? JSON.parse(json) : json;
    const dA = doc.defaults.atom; const dB = doc.defaults.bond;
    const mol = doc.molecules[0];
    const ext = (mol.extensions || []).find(e => e.name === 'rdkitRepresentation') || {};
    const aromA = new Set(ext.aromaticAtoms || []);
    const aromB = new Set(ext.aromaticBonds || []);
    this.rings = (ext.atomRings || []).map(r => r.slice());
    this.atoms = mol.atoms.map((a, i) => {
      const z = a.z ?? dA.z;
      const sym = Z_TO_SYMBOL[z];
      if (!sym) throw new Error(`Elemento con Z=${z} no soportado por el predictor.`);
      return { idx: i, z, el: sym, h: a.impHs ?? dA.impHs, chg: a.chg ?? dA.chg, rad: a.nRad ?? dA.nRad,
        aromatic: aromA.has(i), nbrs: [], rings: [] };
    });
    this.bonds = (mol.bonds || []).map((b, i) => ({ idx: i, a: b.atoms[0], b: b.atoms[1], order: b.bo ?? dB.bo,
      aromatic: aromB.has(i), inRing: false }));
    for (const b of this.bonds) {
      this.atoms[b.a].nbrs.push({ atom: b.b, bond: b.idx });
      this.atoms[b.b].nbrs.push({ atom: b.a, bond: b.idx });
    }
    this.rings.forEach((r, ri) => {
      r.forEach(ai => this.atoms[ai].rings.push(ri));
      for (let k = 0; k < r.length; k++) {
        const bi = this.bondBetween(r[k], r[(k + 1) % r.length]);
        if (bi != null) this.bonds[bi].inRing = true;
      }
    });
    this.rdkitMol = rdkitMol || null;
  }

  get n() { return this.atoms.length; }
  bondBetween(i, j) { const x = this.atoms[i].nbrs.find(n => n.atom === j); return x ? x.bond : null; }
  bond(i, j) { const b = this.bondBetween(i, j); return b == null ? null : this.bonds[b]; }
  other(bond, i) { return bond.a === i ? bond.b : bond.a; }
  isC(i) { return this.atoms[i].el === 'C'; }
  inRing(i) { return this.atoms[i].rings.length > 0; }
  /** sp3: sin enlaces múltiples ni aromaticidad. */
  isSp3(i) {
    const a = this.atoms[i];
    return !a.aromatic && a.nbrs.every(n => this.bonds[n.bond].order === 1);
  }
  heavyDegree(i, within) { return this.atoms[i].nbrs.filter(n => !within || within.has(n.atom)).length; }
  doubleBondPartner(i, el) {
    for (const n of this.atoms[i].nbrs) {
      const b = this.bonds[n.bond];
      if (b.order === 2 && !b.aromatic && (!el || this.atoms[n.atom].el === el)) return n.atom;
    }
    return null;
  }
  tripleBondPartner(i, el) {
    for (const n of this.atoms[i].nbrs) {
      const b = this.bonds[n.bond];
      if (b.order === 3 && (!el || this.atoms[n.atom].el === el)) return n.atom;
    }
    return null;
  }
  /** ¿El carbono i es carbonílico (C=O)? Devuelve el índice del O o null. */
  carbonylO(i) { return this.isC(i) ? this.doubleBondPartner(i, 'O') : null; }

  /** Componentes conexos tras eliminar los enlaces indicados (por índice). */
  componentsWithout(removedBonds, subset) {
    const removed = new Set(removedBonds);
    const allowed = subset || new Set(this.atoms.map(a => a.idx));
    const seen = new Set(); const comps = [];
    for (const start of allowed) {
      if (seen.has(start)) continue;
      const comp = new Set([start]); const stack = [start]; seen.add(start);
      while (stack.length) {
        const cur = stack.pop();
        for (const nb of this.atoms[cur].nbrs) {
          if (removed.has(nb.bond) || !allowed.has(nb.atom) || seen.has(nb.atom)) continue;
          seen.add(nb.atom); comp.add(nb.atom); stack.push(nb.atom);
        }
      }
      comps.push(comp);
    }
    return comps;
  }
  /** Lado de un enlace que contiene al átomo `side` tras romperlo (dentro de subset). */
  sideOf(bondIdx, side, subset, extraRemoved = []) {
    const comps = this.componentsWithout([bondIdx, ...extraRemoved], subset);
    return comps.find(c => c.has(side));
  }

  /** Fórmula de un subconjunto de átomos con ajustes de H por átomo. */
  formulaOf(atomSet, dH = new Map()) {
    const f = {};
    let h = 0;
    for (const i of atomSet) {
      const a = this.atoms[i];
      f[a.el] = (f[a.el] || 0) + 1;
      h += a.h + (dH.get(i) || 0);
    }
    if (h > 0) f.H = (f.H || 0) + h;
    if (h < 0) throw new Error('Conteo de H negativo en fragmento');
    return f;
  }

  get formula() { return this.formulaOf(new Set(this.atoms.map(a => a.idx))); }
  get formulaString() { return formulaToString(this.formula); }
  get monoisotopicMass() { return monoMass(this.formula); }
  get nominalMass() { return nominalMass(this.formula); }
  get rdb() { return rdb(this.formula); }
}

module.exports = { MolGraph };
