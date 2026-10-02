'use strict';
const init = require('@rdkit/rdkit'); const { predict } = require('../src/ms/engine'); const refs = require('./reference');
init().then(R => {
  let tot = 0, hit = 0, basehit = 0;
  for (const c of refs) {
    let r; try { r = predict(R, c.smiles); } catch (e) { console.log(c.name, 'ERROR', e.message, e.stack.split('\n')[1]); continue; }
    const top = r.spectrum.slice().sort((a, b) => b.intensity - a.intensity).slice(0, 8).map(p => `${p.mz}:${p.intensity.toFixed(0)}`).join(' ');
    const refTop = Object.entries(c.ref).sort((a, b) => b[1] - a[1]).map(([m, i]) => `${m}:${i}`).join(' ');
    const get = mz => (r.spectrum.find(p => p.mz == mz) || { intensity: 0 }).intensity;
    const refKeys = Object.keys(c.ref);
    const found = refKeys.filter(m => get(m) >= 2).length; tot += refKeys.length; hit += found;
    const refBase = refKeys.sort((a, b) => c.ref[b] - c.ref[a])[0]; if (r.basePeak && r.basePeak.mz == refBase) basehit++;
    console.log(`${c.name.padEnd(20)} ${r.formula.padEnd(10)} pred: ${top}\n${''.padEnd(31)} ref : ${refTop}   [${found}/${refKeys.length}]`);
  }
  console.log(`\nPicos de referencia predichos (≥2%): ${hit}/${tot}  |  Pico base correcto: ${basehit}/${refs.length}`);
});
