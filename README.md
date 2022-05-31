# MXFP (Macromolecule eXtended FingerPrint)

Open-source version of the MXFP fingerprint based on Python and the RDKit


```
#import required packages
from mxfp import mxfp
from rdkit import Chem

#convert SMILES to rdchem.Mol object
polymyxin_b2_smiles = 'C[C@H]([C@H]1C(=O)NCC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCN)CCN)CC(C)C)CC2=CC=CC=C2)CCN)NC(=O)[C@H](CCN)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCN)NC(=O)CCCCC(C)C)O'
polymyxin_b2_mol = Chem.MolFromSmiles(polymyxin_b2_smiles)

#Initialize MXFP calculator class
MXFP = mxfp.MXFPCalculator()

#Calculate MXFP
polymyxin_b2_mxfp = MXFP.calc_mxfp(polymyxin_b2_mol)

print(polymyxin_b2_mxfp)
```
