# SASA 
This program aims to calculate the Solvent-Accessible Surface Area (SASA) of a protein.

## Setup your environment 

Clone the repository: 

```bash
git clone https://github.com/souadyoujilabadi/SASA_project.git
```

Move to the new directory:

```bash
cd SASA_project
```

Install [anaconda](https://www.anaconda.com/products/distribution).

Install [mamba](https://github.com/mamba-org/mamba):

```bash
conda install mamba -n base -c conda-forge
```

Create the `SASA` conda environment: 

```
mamba env create --file environment.yml
```

Load the `SASA` conda environment:

```
conda activate env_SASA
```

To deactivate an active environment, use:

```
conda deactivate 
```

# Example of use: 

```bash
python3 bin/SASA.py -pdb data/pdb/1bzv.pdb -n 92
```

Results using this example and program: 

3749.25 ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)

9.92 PERCENTAGE OF ACCESSIBILITY