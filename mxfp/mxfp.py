import numpy as np
from math import exp
import random

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops

random.seed(42)

DISTANCE_BINS = [0, 1, 2, 3, 4, 5, 6, 7.1, 8.4, 9.9, 11.6, 13.7, 16.2,
             19.1, 22.6, 26.6, 31.4, 37.1, 43.7, 51.6, 60.9, 71.8, 84.8,
             100.0, 118, 139.3, 164.4, 193.9, 228.9, 270, 318.7]

MXFP_SMARTS = {
'HAT' : ['[*]'],
'HYD' : ['[C]', '[a]', '[S]'],
'ARO' : ['[a]'],
'HBA' : ['[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'],
'HBD' : ['[!$([#6,H0,-,-2,-3])]'],
'POS' : ['[*+]'],
'NEG' : ['[*-]'],
}

LABELS = MXFP_SMARTS.keys()

GAUSSROWS = []
for distance in range(1000):
    if distance == 0:
        gaussrow = [exp(-0.5*(((distance_bin)/(0.09))**2)) for distance_bin in DISTANCE_BINS]
    else:
        gaussrow = [exp(-0.5*(((distance_bin-distance)/(distance*0.09))**2)) for distance_bin in DISTANCE_BINS]
    GAUSSROWS.append(gaussrow)


class MXFPCalculator():

    def __init__(self, distance_bins=DISTANCE_BINS, mxfp_smarts=MXFP_SMARTS, labels=LABELS, gaussrows=GAUSSROWS, dimensionality='2D'):
        self.distance_bins = DISTANCE_BINS 
        self.mxfp_smarts = MXFP_SMARTS
        self.labels = LABELS
        self.gaussrows = GAUSSROWS
        self.dimensionality = dimensionality


    def get_propmat(self, mol):

        """
        Calculates which of the pharmacophore features are attributed
        to the atoms contained in the query molecule.

        Returns a list of the length of number of atoms in the query
        molecule. The i'th element of the list contains a list of pharmacophore
        labels (str) assigned to that atom.
        """
    
        propmat = np.zeros((len(self.labels), mol.GetNumAtoms()))

        for i, label in enumerate(self.labels):
            category = self.mxfp_smarts[label]
            for smarts in category:
                pattern = Chem.MolFromSmarts(smarts)
                matched = False
                for match in mol.GetSubstructMatches(pattern, uniquify=True):
                    for j in match:
                        propmat[i, j] = 1
                    matched = True
                if matched: break
            
        return propmat


    def get_aplist(self, mol, property):
        """
        Determines all atom-pairs belonging to one of the
        selected pharmacophore category.

        Returns a list of n topological distances (int) for n atom-pairs.
        """

        propmat = self.get_propmat(mol)

        if self.dimensionality == '3D':
            distmat = AllChem.Get3DDistanceMatrix(mol)
        else:
            distmat = rdmolops.GetDistanceMatrix(mol)

        for i in range(0, mol.GetNumAtoms()):
            for j in range(i, mol.GetNumAtoms()):
                if propmat[property, i] and propmat[property, j] == 1:
                    yield int(distmat[i, j])


    def get_gausslist(self, mol, property):
        """
        Calculates the gaussian value for each atom pair
        as described in the publication.

        Returns a matrix of n atom pairs and k gaussian values (float).
        """

        gausslist = np.array([self.gaussrows[distance] for distance in self.get_aplist(mol, property)])
        gausssum = np.sum(gausslist.T, axis=0)

        return gausslist, gausssum


    def get_natoms(self, mol, property):
        """
        Calculates the number of atoms belonging to one of 
        the pharmacophore categories. 

        Returns the number of atoms (int)
        """

        propmat = self.get_propmat(mol)
        
        num_atoms = np.count_nonzero(propmat[property] == 1)

        return num_atoms


    def calc_property_mxfp(self, mol, property):

        gausslist, gausssum = self.get_gausslist(mol, property)
        natoms = self.get_natoms(mol, property)

        property_mxfp = np.zeros(len(self.distance_bins))

        for i in range(len(self.distance_bins)):
            if natoms != 0:
                property_mxfp[i] = np.floor((100/(natoms**1.5))*np.sum(gausslist[:, i]/gausssum))
        
        return property_mxfp


    def calc_mxfp(self, mol):
        """
        Calculates the MXFP for one molecule. 

        Returns the MXFP vector of 217 values (int)
        """

        mxfp = np.concatenate([self.calc_property_mxfp(mol, i) for i in range(len(self.labels))], axis=0)
        mxfp = mxfp.astype(int)
    
        return mxfp


    def get_manyxfp(self, mollist):
        """
        Calculates the MXFP for all molecules in a list. 

        Returns a list of MXFP vectors of 217 values (int)
        """

        manyxfp = [self.calc_mxfp(mol) for mol in mollist]

        return manyxfp