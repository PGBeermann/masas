'use strict';
/**
 * Patrón isotópico por convolución polinomial agrupada en masa nominal
 * (método de Yergey, J. A. Int. J. Mass Spectrom. Ion Phys. 52 (1983) 337–349).
 */
const { ISOTOPES, mono } = require('./elements');


/**
 * Devuelve [{offset, rel}] con offset = masa nominal − masa nominal monoisotópica, rel en % del pico máximo.
 */
function isotopePattern(formula, maxOffset = 8) {
  // distribución como arreglo de probabilidades por desplazamiento nominal
  let dist = [1];
  for (const [el, n] of Object.entries(formula)) {
    if (!n) continue;
    const iso = ISOTOPES[el]; const m0 = mono(el).m;
    const single = [];
    for (const i of iso) { const d = i.m - m0; single[d - Math.min(...iso.map(x => x.m - m0))] = (single[d - Math.min(...iso.map(x => x.m - m0))] || 0) + i.ab; }
    const minD = Math.min(...iso.map(x => x.m - m0));
    // potencia por cuadrados sucesivos
    let base = { arr: single, shift: minD };
    let res = { arr: [1], shift: 0 };
    let k = n;
    while (k > 0) {
      if (k & 1) res = mul(res, base, maxOffset);
      base = mul(base, base, maxOffset);
      k >>= 1;
    }
    dist = mul({ arr: dist, shift: 0 }, res, maxOffset);
    dist = shiftTo(dist);
  }
  const max = Math.max(...dist);
  return dist.map((p, i) => ({ offset: i, rel: (p / max) * 100, prob: p })).filter(x => x.rel >= 0.05);
}
function mul(a, b, maxOffset) {
  const arr = [];
  for (let i = 0; i < a.arr.length; i++) { if (!a.arr[i]) continue;
    for (let j = 0; j < b.arr.length; j++) { if (!b.arr[j]) continue;
      const k = i + j; if (k + a.shift + b.shift > maxOffset + 2) continue; arr[k] = (arr[k] || 0) + a.arr[i] * b.arr[j]; } }
  return { arr, shift: a.shift + b.shift };
}
function shiftTo(d) {
  // los elementos con isótopo más ligero que el monoisotópico (B, Se) generan shift negativo; lo recortamos al origen
  const out = [];
  for (let i = 0; i < d.arr.length; i++) { const off = i + d.shift; if (off >= 0) out[off] = (out[off] || 0) + (d.arr[i] || 0); }
  for (let i = 0; i < out.length; i++) out[i] = out[i] || 0;
  return out;
}

module.exports = { isotopePattern };
