VScreenML_v2
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/gandrianov/vscreenml_v2/workflows/CI/badge.svg)](https://github.com/gandrianov/vscreenml_v2/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/gandrianov/VScreenML_v2/branch/master/graph/badge.svg)](https://codecov.io/gh/gandrianov/VScreenML_v2/branch/master)

Modification of original VScreenML module (https://www.pnas.org/content/117/31/18477)

Example of script

```python
import pandas as pd
import sys, json, glob

from vscreenml_v2 import pyrosetta_wrapper
from vscreenml_v2 import rdkit_wrapper
from vscreenml_v2 import rfscore_wrapper
from vscreenml_v2 import binana_wrapper
from oddt.toolkits.extras.rdkit import MolToPDBQTBlock


if __name__ == "__main__":

    pdb = sys.argv[1] #"minimizations/mini_complex_1a28_STR_A_active.pdb"
    params = sys.argv[2] #"LigandAlignmentsPDBParams/mini_complex_1a28_STR_A_active.params"

    features = {"name":pdb.split("/")[-1].split(".")[0]}

     # init pyrosetta
    pyrosetta_wrapper.InitPyRosetta(params + " -mute all")
    bound_pdb = pyrosetta_wrapper.LoadPDB(pdb)
    ligand_id = pyrosetta_wrapper.GetLigandResidueIndex(bound_pdb)
    unbound = pyrosetta_wrapper.MoveLigandFromActiveSite(bound_pdb)
    ligand = pyrosetta_wrapper.ExtractLigandFromPose(bound_pdb, ligand_id)
    protein = pyrosetta_wrapper.ExtractProteinFromPose(bound_pdb, ligand_id)

    sfxn = pyrosetta_wrapper.GetScoringFunction()

    features.update(pyrosetta_wrapper.GetLigandSasa(bound_pdb, protein, ligand,
                                                    ligand_id, sfxn))

    features.update(pyrosetta_wrapper.GetEnergyTerms(bound_pdb, unbound, sfxn))
    features.update(pyrosetta_wrapper.GetMinimizeAppFeatures(bound_pdb, unbound))

    ligand_pdb = pyrosetta_wrapper.GetPDBStringFromPose(ligand)
    protein_pdb = pyrosetta_wrapper.GetPDBStringFromPose(protein)

    ligand_pdb = rdkit_wrapper.ReadLigandPDB(ligand_pdb, params=params)
    protein_pdb = rdkit_wrapper.ReadProteinPDB(protein_pdb)

    ligand_pdbqt = MolToPDBQTBlock(ligand_pdb,
                                   flexible=True,
                                   addHs=False,
                                   computeCharges=True)

    protein_pdbqt = MolToPDBQTBlock(protein_pdb,
                                    flexible=False,
                                    addHs=False,
                                    computeCharges=True)

    features.update(binana_wrapper.GetBinanaFeatures(ligand_pdbqt, protein_pdbqt))
    features.update(rfscore_wrapper.GetRFscore(ligand_pdb, protein_pdb))

    features.update(rdkit_wrapper.CalcRDKitFeatures(ligand_pdb))

    features = pd.DataFrame(features)
    features.to_csv(sys.argv[3], index=False)

    #with open(sys.argv[3], "w") as fwr:
    #    json = json.dumps(features)
    #    fwr.write(json)
    #    fwr.close()

```

### Copyright

Copyright (c) 2021, Grigorii Andrianov


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
