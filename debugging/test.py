import cProfile
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from mxfp import MXFPCalculator

MXFP = MXFPCalculator()

df = pd.read_csv('/home/joule/Desktop/Projects/MXFP_python/debugging/ChEMBL_short24.1.mxfp', header=None, sep='\s+|;', engine='python')
df = df[df.columns[[0, 1]]]
df.columns = ['SMILES', 'Name']
df = df.head(100000)

mollist = []
for i in range(len(df.SMILES)):
    smiles = df['SMILES'].iloc[i]
    mol = Chem.MolFromSmiles(smiles)
    mollist.append(mol)

def calc_mxfp(mollist):
    mxfp_list = MXFP.get_manyxfp(mollist)
    return mxfp_list

cProfile.run('calc_mxfp(mollist)')