#!/usr/bin/env python
import argparse
import pandas as pd
import xgboost as xgb

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-features", required=True)
    parser.add_argument("-output", required=True)
    parser.add_argument("-model", default="DUDE", choices=["DUDE", "DUDEZ", "DeepCoy"])
    
    return parser.parse_args()

def prepare_features(features):

    required_features = ["TotalExposedSasa", "TotalBSA", "InterfaceHydrophobicSasa", "InterfacePolarSasa", "InteractionScore",
                         "FaAtrInteraction", "FaRepInteraction", "FaSolInteraction", "FaElecInteraction", "HBondBbScInteraction",
                         "HBondScInteraction", "GenBonded", "HBInterface", "TotalPackStat", "InterfaceUnsat", "SideFlexAlpha",
                         "SideFlexBeta", "SideFlexOther", "BackFlexAlpha", "BackFlexBeta", "BackFlexOther", "PiPi", "TStacking", 
                         "CationPi", "SaltBridge", "TotalElec", "TotalHBond", "TotalHphobics", "6.6", "6.7", "6.8", "6.9", "6.15",
                         "6.16", "6.17", "6.35", "6.53", "7.6", "7.7", "7.8", "7.9", "7.15", "7.16", "7.17", "7.35", "7.53", "8.6",
                         "8.7", "8.8", "8.9", "8.15", "8.16", "8.17", "8.35", "8.53", "16.6", "16.7", "16.8", "16.9", "16.15", 
                         "16.16", "16.17", "16.35", "16.53", "FCsp3", "NumHAcceptors", "NumHDonors", "MolLogP", "TPSA", "VdwSA"]

    for f in required_features:
        if f not in features.columns:
            raise(f"Feature {f} is not presented in the file")

    return features[required_features].values

if __name__ == "__main__":

    args = args()

    if args.model == "DUDE":
        model_fname = f"{__file__.replace('predict_class.py','')}../data/dude_features/vscreenml_dude.json"
    elif args.model == "DUDEZ":
        model_fname = f"{__file__.replace('predict_class.py','')}../data/dudez_features/vscreenml_dudez.json"
    elif args.model == "DeepCoy":
        model_fname = f"{__file__.replace('predict_class.py','')}../data/deepcoy_features/vscreenml_deepcoy.json"

    data = pd.read_csv(args.features)
    features = prepare_features(data)

    model = xgb.XGBClassifier(use_label_encoder=False)
    model.load_model(model_fname)

    data["Predicted_Class"] = model.predict(features)
    data["VScreenML_Score"] = model.predict_proba(features)[:,1]

    data.to_csv(args.output, index=False)
