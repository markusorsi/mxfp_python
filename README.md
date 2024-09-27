# **MXFP** (**M**acromolecule e**X**tended **F**inger**P**rint)

<img src="https://img.shields.io/pypi/v/mxfp?color=success&label=Version&style=flat-square"/> <img src="https://img.shields.io/badge/Python-3.11-blue?style=flat-square"/> <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square"/>

## Theory

MXFP (Macromolecule eXtended FingerPrint) is a 217-dimensional fuzzy fingerprint that encodes atom pairs of seven pharmacophore groups, making it suitable for the comparison of large molecules and scaffold hopping. Each atom is assigned to one or more of the following pharmacophore categories:

1. Heavy Atoms (HAC)
2. Hydrophobic Atoms (HYD)
3. Aromatic Atoms (ARO)
4. H-bond Acceptors (HBA)
5. H-bond Donors (HBD)
6. Positively Charged Atoms (POS)
7. Negatively Charged Atoms (NEG)

For each pharmacophore category, all possible atom pairs are determined and converted to a Gaussian with an 18% width centered at the atom pair topological (2D) or Euclidean (3D) distance. This Gaussian is sampled at 31 distance bins \((d_{i})\) spanning from \(d_{0} = 0\) to \(d_{30} = 317.8\) bonds at exponentially increasing intervals. The Gaussian value \(g_{jk}(d_{i})\) for an atom pair with distance \(d_{jk}\) is calculated as follows:

$$ g_{jk}(d_{i}) = e^{- \frac{1}{2} \left( \frac{d_{i} - d_{jk}}{d_{jk} \cdot 0.09} \right)^2 } $$

Each of the obtained 31 Gaussian values is normalized to the sum of all 31 values, \(s_{jk}\), to ensure that every atom pair contributes equally to the fingerprint.

$$ s_{jk}  = \sum_{i = 0}^{30} g_{jk}(d_{i}) $$

The sum of normalized Gaussian contributions from all atom pairs of one pharmacophore category at distance \(d_{i}\) is normalized by the total number of atoms in that category \(N_{c}\) raised to the power of 1.5 to reduce the sensitivity of the fingerprint to molecule size. This value is multiplied by 100 and rounded to yield the final fingerprint bit value \(v_{ci}\).

$$ v_{ci} = \frac{100}{N_{c}^{1.5}} \sum_{j = 1} \sum_{k = 1} \frac{g_{jk}(d_{i})}{s_{jk}} $$

The 31 fingerprint bit values from each of the 7 atom categories are concatenated, yielding the 217-dimensional fingerprint vector.

## Getting Started

### Prerequisites

You will need the following prerequisites: 

- [Python 3.6](https://www.python.org)
- [NumPy](https://numpy.org)
- [RDKit](https://www.rdkit.org)

## Installing MXFP

There are several ways to get started using MXFP:

#### **Installing via GitHub**

To obtain a local copy of the project, clone the GitHub repository:

```bash
git clone https://github.com/reymond-group/mxfp_python.git
```

#### **Installing via Conda**

To create a ready-to-use Conda environment, download the `mxfp.yml` file from the repository and run the following commands:

```bash
conda env create -f mxfp.yml
conda activate mxfp
```

#### **Installing via pip**

To install MXFP in an existing Conda environment, activate the environment and install MXFP via pip:

```bash
conda activate my_environment
pip install mxfp
```

## Using MXFP

In your Python script or Jupyter notebook:

1. Import the required libraries (RDKit, MXFP).
2. Convert SMILES to an `rdchem.Mol` object with RDKit (optional).
3. Initialize the `MXFPCalculator` class.
4. Calculate the MXFP of your molecule either from the `rdchem.Mol` object or directly from SMILES.

```python
# Import the required libraries (RDKit, MXFP)
from rdkit import Chem
from mxfp import MXFPCalculator

# Convert SMILES to rdchem.Mol object with RDKit (optional)
polymyxin_b2_smiles = 'C[C@H]([C@H]1C(=O)NCC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCN)CCN)CC(C)C)CC2=CC=CC=C2)CCN)NC(=O)[C@H](CCN)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCN)NC(=O)CCCCC(C)C)O'
polymyxin_b2_mol = Chem.MolFromSmiles(polymyxin_b2_smiles)

# Initialize the MXFPCalculator class
calculator = MXFPCalculator()

# Calculate MXFP of your molecule from rdchem.Mol object or directly from SMILES
polymyxin_b2_mxfp = calculator.mxfp_from_mol(polymyxin_b2_mol)  # from rdchem.Mol object
polymyxin_b2_mxfp = calculator.mxfp_from_smiles(polymyxin_b2_smiles)  # from SMILES
```

If you are working with 3D coordinates and wish to use Euclidean atom-pair distances instead of topological distances, initialize the `MXFPCalculator` class with the `dimensionality='3D'` parameter. Note that you will not be able to calculate MXFP from a SMILES string if you use the '3D' option, so you need to provide an `rdchem.Mol` object.

```python
# Initialize the MXFPCalculator class for 3D calculations
calculator_3d = MXFPCalculator(dimensionality='3D')
```

## License

This project is licensed under the [MIT License](https://choosealicense.com/licenses/mit/).

## Contact

<img src="https://img.shields.io/twitter/follow/reymondgroup?style=social"/> 
<img src="https://img.shields.io/twitter/follow/markusorsi?style=social"/>