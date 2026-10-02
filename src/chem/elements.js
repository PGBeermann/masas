'use strict';
/**
 * Datos isotópicos (masas monoisotópicas y abundancias naturales).
 * Fuentes:
 *  - Wang, M. et al. "The AME 2020 atomic mass evaluation (II)". Chinese Physics C 45 (2021) 030003.
 *  - Meija, J. et al. "Isotopic compositions of the elements 2013 (IUPAC Technical Report)".
 *    Pure Appl. Chem. 88 (2016) 293–306 (abundancias representativas).
 */
const ELECTRON_MASS = 0.000548579909; // u (CODATA 2018)

const ISOTOPES = {
  H:  [{ m: 1, mass: 1.00782503223, ab: 0.999885 }, { m: 2, mass: 2.01410177812, ab: 0.000115 }],
  B:  [{ m: 10, mass: 10.01293695, ab: 0.199 }, { m: 11, mass: 11.00930536, ab: 0.801 }],
  C:  [{ m: 12, mass: 12.0, ab: 0.9893 }, { m: 13, mass: 13.00335483507, ab: 0.0107 }],
  N:  [{ m: 14, mass: 14.00307400443, ab: 0.99636 }, { m: 15, mass: 15.00010889888, ab: 0.00364 }],
  O:  [{ m: 16, mass: 15.99491461957, ab: 0.99757 }, { m: 17, mass: 16.9991317565, ab: 0.00038 }, { m: 18, mass: 17.99915961286, ab: 0.00205 }],
  F:  [{ m: 19, mass: 18.99840316273, ab: 1.0 }],
  Si: [{ m: 28, mass: 27.97692653465, ab: 0.92223 }, { m: 29, mass: 28.9764946649, ab: 0.04685 }, { m: 30, mass: 29.973770136, ab: 0.03092 }],
  P:  [{ m: 31, mass: 30.97376199842, ab: 1.0 }],
  S:  [{ m: 32, mass: 31.9720711744, ab: 0.9499 }, { m: 33, mass: 32.9714589098, ab: 0.0075 }, { m: 34, mass: 33.967867004, ab: 0.0425 }, { m: 36, mass: 35.96708071, ab: 0.0001 }],
  Cl: [{ m: 35, mass: 34.968852682, ab: 0.7576 }, { m: 37, mass: 36.965902602, ab: 0.2424 }],
  Br: [{ m: 79, mass: 78.9183376, ab: 0.5069 }, { m: 81, mass: 80.9162897, ab: 0.4931 }],
  I:  [{ m: 127, mass: 126.9044719, ab: 1.0 }],
  Se: [{ m: 74, mass: 73.9224759, ab: 0.0089 }, { m: 76, mass: 75.9192137, ab: 0.0937 }, { m: 77, mass: 76.9199142, ab: 0.0763 }, { m: 78, mass: 77.9173091, ab: 0.2377 }, { m: 80, mass: 79.9165218, ab: 0.4961 }, { m: 82, mass: 81.9166995, ab: 0.0873 }]
};

const Z_TO_SYMBOL = { 1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 34: 'Se', 35: 'Br', 53: 'I' };
const SYMBOL_TO_Z = Object.fromEntries(Object.entries(Z_TO_SYMBOL).map(([z, s]) => [s, Number(z)]));

/** Isótopo más abundante (define la masa monoisotópica). */
function mono(sym) {
  const iso = ISOTOPES[sym];
  if (!iso) throw new Error(`Elemento no soportado: ${sym}`);
  return iso.reduce((a, b) => (b.ab > a.ab ? b : a));
}

/** Orden de Hill para fórmulas. */
function hillOrder(formula) {
  const keys = Object.keys(formula).filter(k => formula[k] > 0);
  const hasC = keys.includes('C');
  return keys.sort((a, b) => {
    if (hasC) {
      if (a === 'C') return -1; if (b === 'C') return 1;
      if (a === 'H') return -1; if (b === 'H') return 1;
    }
    return a.localeCompare(b);
  });
}

function formulaToString(formula) {
  return hillOrder(formula).map(k => k + (formula[k] > 1 ? formula[k] : '')).join('');
}

/** Masa monoisotópica de una fórmula neutra (u). */
function monoMass(formula) {
  let m = 0;
  for (const [el, n] of Object.entries(formula)) if (n) m += mono(el).mass * n;
  return m;
}

/** Masa nominal (enteros del isótopo más abundante). */
function nominalMass(formula) {
  let m = 0;
  for (const [el, n] of Object.entries(formula)) if (n) m += mono(el).m * n;
  return m;
}

/** Grados de insaturación (anillos + dobles enlaces), RDB = C − H/2 − X/2 + N/2 + 1. */
function rdb(formula) {
  const c = (formula.C || 0) + (formula.Si || 0);
  const h = (formula.H || 0) + (formula.F || 0) + (formula.Cl || 0) + (formula.Br || 0) + (formula.I || 0);
  const n = (formula.N || 0) + (formula.P || 0) + (formula.B || 0);
  return c - h / 2 + n / 2 + 1;
}

module.exports = { ELECTRON_MASS, ISOTOPES, Z_TO_SYMBOL, SYMBOL_TO_Z, mono, formulaToString, monoMass, nominalMass, rdb, hillOrder };
