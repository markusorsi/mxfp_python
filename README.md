# MXFP (Macromolecule eXtended FingerPrint)

Open source python version of MXFP based on the RDKit package. 

## Getting started

### Prerequisites

You will need following prerequisites: 

* [Python 3.9](https://www.python.org)
* [NumPy](https://numpy.org)
* [RDKit](https://www.rdkit.org)

## Installing MXFP

There are several ways in which you can get started using MXFP:

* ### Installing via GitHub

To have local copy of the project clone the GitHub repository as follows:

```console
git clone https://github.com/reymond-group/mxfp_python.git
cd mxfp_python
```

* ### Installing via Conda 

To create a ready to use Conda environment download the mxfp.yml file from the repository and run following commands:

```console
conda env create -f mxfp.yml
conda activate mxfp
```

* ### Installing via pip

To install mxfp on an existing Conda environment activate the environment and install mxfp via pip using following commands:

```console
conda activate myenv
pip install mxfp
```

## Using MXFP

In your python script:

* Import the required libraries (RDKit, MXFP)
* Convert SMILES to rdchem.Mol object with RDKit
* Initialize the MXFPCalculator class
* Calculate MXFP of your molecule

<br>

```python
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

## License
[MIT](https://choosealicense.com/licenses/mit/)
