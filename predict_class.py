import argparse

import pickle 
import pandas as pd
import xgboost as xgb

from vscreenml_v2 import pyrosetta_wrapper
from vscreenml_v2 import rdkit_wrapper
from vscreenml_v2 import rfscore_wrapper
from vscreenml_v2 import binana_wrapper
from oddt.toolkits.extras.rdkit import MolToPDBQTBlock

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", required=True)
    parser.add_argument("-params", required=True)
    parser.add_argument("-xgbmodel", required=True)
    parser.add_argument("-output", required=True)

    return parser.parse_args()


def calc_features(pdb, params):

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

    return pd.DataFrame(features)


def predict_class(features, model):

    with open(model, "rb") as f:
        model = pickle.load(f)

        features["predicted_class"] = model.predict(features.drop("name", axis=1))
        features["predicted_proba"]  = model.predict_proba(features.drop(["name", "predicted_class"], axis=1))[:,1]

        return features


if __name__ == "__main__":

    args = args()
    features = calc_features(args.s, args.params)
    features = predict_class(features, model)

    features.to_csv(args.output, index=False)
