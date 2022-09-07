# Solvent-accessible surface area (SASA) of a protein 

## Setup your environment

Clone the repository: 

```bash
git clone https://github.com/souadyoujilabadi/projects.git
```

Create the `SASA` conda environment: 

```
conda env create -n env_SASA --file environment.yml 
```

Load the `SASA` conda environment:

```
conda activate env_SASA
```

To deactivate an active environment, use:

```
conda deactivate 
```

#conda remove -n env_SASA --all
#conda env list



#conda env export -n env_SASA --no-builds | grep -v "^prefix:" > environment.yml  



Example of use: 
python3 SASA.py -pdb "1bzv.pdb" -n 1000

Results using this example and program: 
Au mois avec cette example : moyenne

Use already existing libraries 
