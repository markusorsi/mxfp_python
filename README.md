# **MXFP** (**M**acromolecule e**X**tended **F**inger**P**rint)

<img src="https://img.shields.io/pypi/v/mxfp?color=success&label=Version&style=flat-square"/>
<img src="https://img.shields.io/badge/Python-3.9-blue?style=flat-square"/>
<img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square"/>

<br>

MXFP is a 217D fingerprint counting atom-pairs using a fuzzy approach to assign atom-pairs to distance bins.

For the MXFP presented here, we compute exact topological distances between atoms, which is suitable for any molecule. Furthermore, we use seven atom categories by additionally computing aromatic (AR), H-bond donor (HBD), and H-bond acceptor atoms (HBA), which are important to differentiate molecules such as polycyclic aromatic hydrocarbons, oligosaccharides, and oligonucleotides. As for 3DP and 2DP, we do not consider cross-category atom pairs in MXFP. Each atom pair is converted to a Gaussian of 18 % width centered at the atom pair topological distance, which is the shortest path between the two atoms counted in bonds. This Gaussian is then sampled at 31 distances d<sub>i</sub> spanning from d<sub>0</sub> = 0 to d<sub>30</sub> = 317.8 bonds at exponentially increasing intervals. 


The sampled Gaussian values are normalized and added to the MXFP distance bins for the corresponding atom-pair category, and distance bins of each category are normalized to size (see Methods and Equation 1 for details). Sampling atom-pair Gaussians at exponentially increasing distances allows to describe molecules up to a very large size using only a limited number of dimensions in the atom pair fingerprint. The approach furthermore partly erases differences between atom pairs separated by a similar number of bonds at large distances, which favors the perception of global molecular shape over structural detail.

<br>

## Getting started

### Prerequisites

You will need following prerequisites: 

* [Python 3.9](https://www.python.org)
* [NumPy](https://numpy.org)
* [RDKit](https://www.rdkit.org)

<br>

## Installing MXFP

There are several ways in which you can get started using MXFP:

<br>

### **Installing via GitHub**

To have local copy of the project clone the GitHub repository as follows:

```console
git clone https://github.com/reymond-group/mxfp_python.git
```

<br>

### **Installing via Conda**

To create a ready to use Conda environment download the mxfp.yml file from the repository and run following commands:

```console
conda env create -f mxfp.yml
```

```console
conda activate mxfp
```

<br>

### **Installing via pip**

To install mxfp on an existing Conda environment activate the environment and install mxfp via pip using following commands:

```console
conda activate myenv
```

```console
pip install mxfp
```
<br>

## Using MXFP

In your python script:

1. Import the required libraries (RDKit, MXFP)
2. Convert SMILES to rdchem.Mol object with RDKit
3. Initialize the MXFPCalculator class
4. Calculate MXFP of your molecule

```python
#Import the required libraries (RDKit, MXFP)
from rdkit import Chem
from mxfp import mxfp

#Convert SMILES to rdchem.Mol object with RDKit
polymyxin_b2_smiles = 'C[C@H]([C@H]1C(=O)NCC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N1)CCN)CCN)CC(C)C)CC2=CC=CC=C2)CCN)NC(=O)[C@H](CCN)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCN)NC(=O)CCCCC(C)C)O'
polymyxin_b2_mol = Chem.MolFromSmiles(polymyxin_b2_smiles)

#Initialize the MXFPCalculator class
MXFP = mxfp.MXFPCalculator()

#Calculate MXFP of your molecule
polymyxin_b2_mxfp = MXFP.calc_mxfp(polymyxin_b2_mol)
print(polymyxin_b2_mxfp)
```

## License
[MIT](https://choosealicense.com/licenses/mit/)

## Contact

<img src="https://img.shields.io/twitter/follow/reymondgroup?style=social"/> 
<img src="https://img.shields.io/twitter/follow/markusorsi?style=social"/>