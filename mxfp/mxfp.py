import numpy as np
from math import exp

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, rdchem


DISTANCE_BINS = [0, 1, 2, 3, 4, 5, 6, 7.1, 8.4, 9.9, 11.6, 13.7, 16.2,
             19.1, 22.6, 26.6, 31.4, 37.1, 43.7, 51.6, 60.9, 71.8, 84.8,
             100.0, 118, 139.3, 164.4, 193.9, 228.9, 270, 318.7]

MXFP_SMARTS = {
'HAC' : ['[*]'],
'HYD' : ['[C]', '[a]', '[S]'],
'ARO' : ['[a]'],
'HBA' : ['[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]'],
'HBD' : ['[!$([#6,H0,-,-2,-3])]'],
'POS' : ['[*+]'],
'NEG' : ['[*-]'],
}

class MXFPCalculator:

    def __init__(self, dimensionality:str='2D', max_dist:int=1000, categories:list=['HAC', 'HYD', 'ARO', 'HBA', 'HBD', 'POS', 'NEG']) -> None:
        """
        MXFP calculator class

        Parameters 
        ----------
        dimensionality : '2D'
            Default coordinates are 2-dimensional. If you are working with 3-dimensional coordinates then set parameter to '3D'. 
        
        max_dist: 1000
            Maximum atom pair distance found in the molecule. Lower values speed up calculations.

        """
        self.distance_bins = DISTANCE_BINS 
        self.mxfp_smarts = MXFP_SMARTS
        self.dimensionality = dimensionality
        self.max_dist = max_dist
        self.labels = categories
    
        gaussrows = []
        for distance in range(self.max_dist):
            if distance == 0:
                gaussrow = [exp(-0.5*(((distance_bin)/(0.09))**2)) for distance_bin in DISTANCE_BINS]
            else:
                gaussrow = [exp(-0.5*(((distance_bin-distance)/(distance*0.09))**2)) for distance_bin in DISTANCE_BINS]
            gaussrows.append(gaussrow)
        
        self.gaussrows = gaussrows

    def get_property_matrix(self, mol:rdchem.Mol) -> np.ndarray:
        """ 

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.

        Returns
        -------
        property matrix : np.ndarray 
            (7 pharmacophore categories) x (number of heavy atoms) matrix.\n 
            Each heavy atom is represented by a column of 7 binary values,\n
            which state whether the selected atom possesses the property (1) or does not possess the property (0).

        """

        propmat = np.zeros((len(self.labels), mol.GetNumAtoms()))

        try: 
            for i, label in enumerate(self.labels):
                category = self.mxfp_smarts[label]
                for smarts in category:
                    pattern = Chem.MolFromSmarts(smarts)
                    matched = False
                    for match in mol.GetSubstructMatches(pattern, uniquify=True):
                        for j in match:
                            propmat[i, j] = 1
                        matched = True
                    if matched: 
                        break
            return propmat
        except:
            print("Oopsie! Something went wrong while trying to generate the property matrix.\n"
                "Have you checked if the input file format and/or the input molecule is valid?")


    def yield_ap(self, mol:rdchem.Mol, property:int):
        """

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.
        
        property : int
            Integer value (0-6) designating the pharmacophore property for which to find all atom pairs. \n
            0 = Heavy Atom Count (HAC)\n
            1 = Hydrophobic Atoms (HYD)\n
            2 = Aromatic Atoms (ARO)\n
            3 = H-bond Acceptors (HBA)\n
            4 = H-bond Donors (HBD)\n
            5 = Positively Charged Atoms (POS)\n
            6 = Negatively Charged Atoms (NEG)\n

        Returns
        -------
         atom pairs : Generator
            Returns a generator object containing all pairs of atoms displaying the selected pharmacophore property.

        """

        propmat = self.get_property_matrix(mol)

        if self.dimensionality == '3D':
            distmat = AllChem.Get3DDistanceMatrix(mol)
        else:
            distmat = rdmolops.GetDistanceMatrix(mol)

        try:
            for i in range(0, mol.GetNumAtoms()):
                for j in range(i, mol.GetNumAtoms()):
                    if propmat[property, i] and propmat[property, j] == 1:
                        yield int(distmat[i, j])
        except:
            print(f"Uh-Oh! Couldn't calculate the atom pairs for property {property}")


    def get_gausslist(self, mol:rdchem.Mol, property:int) -> tuple:
        """

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.

        property : int
            Integer value (0-6) designating the pharmacophore property for which to find all atom pairs. \n
            0 = Heavy Atom Count (HAC)\n
            1 = Hydrophobic Atoms (HYD)\n
            2 = Aromatic Atoms (ARO)\n
            3 = H-bond Acceptors (HBA)\n
            4 = H-bond Donors (HBD)\n
            5 = Positively Charged Atoms (POS)\n
            6 = Negatively Charged Atoms (NEG)\n

        Returns
        -------
        gausslist : np.ndarray 
            (n atom pairs) x (31 Gaussian values for each atom pair) matrix.\n
            Gaussian values are calculated beforehand.

        gaussum : int
            List containing the sum of the 31 Gaussian values for each atom pair. 

        """

        try: 
            gausslist = np.array([self.gaussrows[distance] for distance in self.yield_ap(mol, property)])
            gausssum = np.sum(gausslist.T, axis=0)
            return gausslist, gausssum
        except:
            print(f"Yikes! Couldn't generate the list and the sum of the Gaussian values for property {property}")


    def count_atoms(self, mol:rdchem.Mol, property:int) -> int:
        """

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.

        property : int
            Integer value (0-6) designating the pharmacophore property for which to find all atom pairs. \n
            0 = Heavy Atom Count (HAC)\n
            1 = Hydrophobic Atoms (HYD)\n
            2 = Aromatic Atoms (ARO)\n
            3 = H-bond Acceptors (HBA)\n
            4 = H-bond Donors (HBD)\n
            5 = Positively Charged Atoms (POS)\n
            6 = Negatively Charged Atoms (NEG)\n

        Returns
        -------
        num_atoms : int 
            Count of all atoms that possess the selected pharmacophore property.

        """

        propmat = self.get_property_matrix(mol)

        try:
            num_atoms = np.count_nonzero(propmat[property] == 1)
            return num_atoms
        except:
            print(f"Ooof! Couldn't count the atoms with property {property}")


    def calc_property_mxfp(self, mol:rdchem.Mol, property:int) -> np.ndarray:
        """

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.

        property : int
            Integer value (0-6) designating the pharmacophore property for which to find all atom pairs. \n
            0 = Heavy Atom Count (HAC)\n
            1 = Hydrophobic Atoms (HYD)\n
            2 = Aromatic Atoms (ARO)\n
            3 = H-bond Acceptors (HBA)\n
            4 = H-bond Donors (HBD)\n
            5 = Positively Charged Atoms (POS)\n
            6 = Negatively Charged Atoms (NEG)\n

        Returns
        -------
        property_mxfp : np.ndarray
            MXFP values for the selected property (31 values)

        """
        
        gausslist, gausssum = self.get_gausslist(mol, property)
        natoms = self.count_atoms(mol, property)

        property_mxfp = np.zeros(len(self.distance_bins))

        try:
            for i in range(len(self.distance_bins)):
                if natoms != 0:
                    property_mxfp[i] = np.floor((100/(natoms**1.5))*np.sum(gausslist[:, i]/gausssum))
            return property_mxfp
        except:
            print(f"Yuck! Something didn't work as expected. Couldn't generate MXFP values for property {property}")
        

    def mxfp_from_mol(self, mol:rdchem.Mol) -> np.ndarray:
        """
        Calculates the MXFP for a selected molecule. 

        Parameters 
        ----------
        mol : rdkit.Chem.rdchem.Mol
            RDKit mol object.

        Returns
        -------
        mxfp : np.ndarray
            MXFP for the selected molecule (217 values).
        """

        try:
            mxfp = np.concatenate([self.calc_property_mxfp(mol, i) for i in range(len(self.labels))], axis=0)
            mxfp = mxfp.astype(int)
            return mxfp
        except:
            print("Ouch! Couldn't calculate MXFP for this molecule.")


    def mxfp_from_smiles(self, smiles:str) -> np.ndarray:
            """
            Calculates the MXFP for a selected molecule. 

            Parameters 
            ----------
            smiles : str
                A SMILES string of .

            Returns
            -------
            mxfp : np.ndarray
                MXFP for the selected molecule (217 values).
            """

            try:
                mol = Chem.MolFromSmiles(smiles)
                mxfp = np.concatenate([self.calc_property_mxfp(mol, i) for i in range(len(self.labels))], axis=0)
                mxfp = mxfp.astype(int)
                return mxfp
            except:
                print("Ouch! Couldn't calculate MXFP for this molecule. Have you checked if the SMILES string is valid?")