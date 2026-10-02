'use strict';
// Compuestos de referencia con picos EI principales (aprox., NIST WebBook SRD 69)
module.exports = [
  { name: '2-hexanona', smiles: 'CCCCC(C)=O', ref: { 43: 100, 58: 60, 100: 8, 71: 8, 85: 5 } },
  { name: '3-pentanona', smiles: 'CCC(=O)CC', ref: { 57: 100, 29: 80, 86: 30 } },
  { name: 'tolueno', smiles: 'Cc1ccccc1', ref: { 91: 100, 92: 75, 65: 10, 39: 8 } },
  { name: 'etilbenceno', smiles: 'CCc1ccccc1', ref: { 91: 100, 106: 30, 77: 8, 51: 8 } },
  { name: 'butirofenona', smiles: 'CCCC(=O)c1ccccc1', ref: { 105: 100, 120: 45, 77: 40, 148: 15, 51: 15 } },
  { name: 'trietilamina', smiles: 'CCN(CC)CC', ref: { 86: 100, 58: 30, 101: 18, 30: 20 } },
  { name: '1-bromobutano', smiles: 'CCCCBr', ref: { 57: 100, 41: 55, 29: 40, 136: 8, 138: 8 } },
  { name: '1-bromohexano', smiles: 'CCCCCCBr', ref: { 135: 100, 137: 95, 43: 80, 41: 50, 55: 40, 85: 20 } },
  { name: 'ciclohexeno', smiles: 'C1CC=CCC1', ref: { 67: 100, 82: 40, 54: 60, 41: 30 } },
  { name: 'n-octano', smiles: 'CCCCCCCC', ref: { 43: 100, 57: 35, 85: 25, 71: 20, 41: 30, 29: 20, 114: 6 } },
  { name: '2-butanol', smiles: 'CCC(C)O', ref: { 45: 100, 59: 20, 31: 20, 74: 1 } },
  { name: 'benzoato de metilo', smiles: 'COC(=O)c1ccccc1', ref: { 105: 100, 77: 70, 136: 30, 51: 30 } },
  { name: 'benzaldehído', smiles: 'O=Cc1ccccc1', ref: { 105: 100, 106: 95, 77: 90, 51: 45 } },
  { name: 'dietil éter', smiles: 'CCOCC', ref: { 31: 100, 59: 60, 74: 60, 45: 30, 29: 50 } },
  { name: 'ácido butanoico', smiles: 'CCCC(=O)O', ref: { 60: 100, 73: 25, 42: 15, 45: 15, 88: 3 } },
  { name: 'acetanilida', smiles: 'CC(=O)Nc1ccccc1', ref: { 93: 100, 135: 50, 43: 30, 66: 15, 65: 15 } },
  { name: 'nitrobenceno', smiles: 'O=[N+]([O-])c1ccccc1', ref: { 77: 100, 123: 60, 51: 55, 93: 15, 65: 15 } },
  { name: 'fenol', smiles: 'Oc1ccccc1', ref: { 94: 100, 66: 30, 65: 25, 39: 15 } },
  { name: 'benceno', smiles: 'c1ccccc1', ref: { 78: 100, 77: 20, 52: 20, 51: 20 } },
  { name: '2-metilpropan-2-ol', smiles: 'CC(C)(C)O', ref: { 59: 100, 31: 30, 41: 15, 43: 12 } }
];
