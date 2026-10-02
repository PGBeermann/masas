'use strict';
const init = require('@rdkit/rdkit'); const { predict } = require('../src/ms/engine');
const set = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'CC(=O)Oc1ccccc1C(=O)O', 'CCN(CC)CC(=O)Nc1c(C)cccc1C',
  'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'C[Si](C)(C)OC', 'CCOP(=O)(OCC)OCC', 'CSC', 'c1ccc2ccccc2c1', 'CC1=CCC(CC1)C(=C)C',
  'O=C1CCCCC1', 'N#Cc1ccccc1', 'CCCCCCCCCCCCCCCC', 'ClC(Cl)(Cl)Cl', 'FC(F)(F)c1ccccc1', 'CC(=O)OCC', 'OCC(O)CO', 'c1ccncc1', 'C1CCC2CCCCC2C1', 'CCOC(=O)C=Cc1ccccc1',
  'Ic1ccccc1', 'CC(C)(C)c1ccccc1', 'C=CC=C', 'CC#N', 'C[N+](C)(C)C.[Cl-]', 'O=C(O)c1ccccc1O', 'CCCCCCC=O', 'NCCc1ccc(O)c(O)c1', 'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O', 'C1CCCCC1'];
init().then(R => {
  let ok = 0;
  for (const s of set) {
    try { const t = Date.now(); const r = predict(R, s); ok++;
      const top = r.spectrum.slice().sort((a, b) => b.intensity - a.intensity).slice(0, 6).map(p => p.mz + ':' + p.intensity.toFixed(0)).join(' ');
      console.log(`${s.padEnd(50).slice(0,50)} ${r.formula.padEnd(10)} ${String(Date.now() - t).padStart(5)} ms  ${top}  expl=${r.explanations.length}`);
    } catch (e) { console.log(s, 'ERROR:', e.message); }
  }
  console.log(`OK ${ok}/${set.length}`);
});
