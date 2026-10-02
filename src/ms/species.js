'use strict';
/**
 * Especie química (ion o neutro) definida como subconjunto de átomos del grafo
 * de la molécula original + modificaciones locales (H transferidos, órdenes de enlace,
 * cargas, radicales). Permite calcular fórmula/masa y generar SMILES/SVG con RDKit.
 */
const { formulaToString, monoMass, nominalMass, ELECTRON_MASS, rdb } = require('../chem/elements');

class Species {
  constructor(g, atoms, o = {}) {
    this.g = g;
    this.atoms = atoms instanceof Set ? atoms : new Set(atoms || []);
    this.dH = o.dH || new Map();          // átomo -> ΔH
    this.bo = o.bo || new Map();          // índice de enlace -> nuevo orden
    this.chg = o.chg || new Map();        // átomo -> carga formal
    this.rad = o.rad || new Map();        // átomo -> nº electrones desapareados
    this.extraBonds = o.extraBonds || []; // [{a,b,order}] enlaces nuevos (p. ej. ion halonio cíclico)
    this.formulaOverride = o.formula || null; // especies descritas sólo por fórmula (reordenamientos profundos)
    this.tags = o.tags || {};
    this.label = o.label || '';
    this.oddElectron = !!o.oddElectron;
    this.depictAsNeutral = !!o.depictAsNeutral; // iones-radical deslocalizados: se dibuja el esqueleto neutro y se rotula "+•"
  }
  get formula() { return this.formulaOverride || this.g.formulaOf(this.atoms, this.dH); }
  get formulaString() { return formulaToString(this.formula); }
  get nominal() { return nominalMass(this.formula); }
  get mono() { return monoMass(this.formula); }
  /** m/z exacto del catión monocargado. */
  get mzExact() { return this.mono - ELECTRON_MASS; }
  get rdb() { return rdb(this.formula); }
  key() { return this.formulaString + '|' + [...this.atoms].sort((a, b) => a - b).join(',') + '|' + [...this.dH].sort().join(';'); }

  clone(o = {}) {
    return new Species(this.g, new Set(this.atoms), {
      dH: new Map(this.dH), bo: new Map(this.bo), chg: new Map(this.chg), rad: new Map(this.rad),
      extraBonds: this.extraBonds.slice(), formula: this.formulaOverride, tags: { ...this.tags },
      label: this.label, oddElectron: this.oddElectron, depictAsNeutral: this.depictAsNeutral, ...o
    });
  }

  /** Construye el JSON commonchem/rdkitjson del subgrafo. */
  toRDKitJSON({ neutralSkeleton = false } = {}) {
    if (this.formulaOverride || !this.atoms.size) return null;
    const g = this.g; const map = new Map(); const atoms = [];
    for (const i of [...this.atoms].sort((a, b) => a - b)) {
      map.set(i, atoms.length);
      const a = g.atoms[i];
      const rec = { z: a.z, impHs: Math.max(0, a.h + (this.dH.get(i) || 0)) };
      if (neutralSkeleton && a.chg) rec.chg = a.chg;
      if (!neutralSkeleton) {
        const c = (this.chg.get(i) || 0) + (a.chg || 0); if (c) rec.chg = c;
        const r = this.rad.get(i); if (r) rec.nRad = r;
      }
      atoms.push(rec);
    }
    const bonds = [];
    for (const b of g.bonds) {
      if (!this.atoms.has(b.a) || !this.atoms.has(b.b)) continue;
      const order = this.bo.has(b.idx) ? this.bo.get(b.idx) : b.order;
      if (order === 0) continue;
      bonds.push({ bo: order, atoms: [map.get(b.a), map.get(b.b)] });
    }
    for (const e of this.extraBonds) bonds.push({ bo: e.order || 1, atoms: [map.get(e.a), map.get(e.b)] });
    return {
      commonchem: { version: 10 },
      defaults: { atom: { z: 6, impHs: 0, chg: 0, nRad: 0, isotope: 0, stereo: 'unspecified' }, bond: { bo: 1, stereo: 'unspecified' } },
      molecules: [{ atoms, bonds }]
    };
  }

  /** SMILES y SVG mediante RDKit (null si no es representable). */
  depict(RDKit, { svg = true, w = 220, h = 160 } = {}) {
    const tryBuild = (neutral) => {
      const js = this.toRDKitJSON({ neutralSkeleton: neutral });
      if (!js) return null;
      let mol = null;
      try {
        mol = RDKit.get_mol(JSON.stringify(js));
        if (!mol || !mol.is_valid()) { if (mol) mol.delete(); return null; }
        const out = { smiles: mol.get_smiles() };
        if (svg) out.svg = mol.get_svg(w, h);
        mol.delete();
        return out;
      } catch (e) { if (mol) try { mol.delete(); } catch (_) { /* noop */ } return null; }
    };
    return tryBuild(this.depictAsNeutral) || tryBuild(true);
  }
}

module.exports = { Species };
