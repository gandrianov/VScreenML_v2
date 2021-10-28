from .binana.interactions import get_all_interactions
from .binana.load_ligand_receptor import from_texts

def GetBinanaFeatures(ligand, protein):

    ligand, protein = from_texts(ligand, protein)
    all_interacts = get_all_interactions(ligand, protein)
    
    for key in all_interacts["cat_pi"]["counts"].keys():
        all_interacts["pi_pi"]["counts"][key] = all_interacts["cat_pi"]["counts"][key]

    data = {}

    d =  all_interacts["active_site_flexibility"]["counts"]
    data["SIDE_flex_ALPHA"] = d.get("SIDECHAIN_ALPHA", 0.0)
    data["SIDE_flex_BETA"] = d.get("SIDECHAIN_BETA", 0.0)
    data["SIDE_flex_OTHER"] = d.get("SIDECHAIN_OTHER", 0.0)
    data["BACK_flex_ALPHA"] = d.get("BACKBONE_ALPHA", 0.0)
    data["BACK_flex_BETA"] = d.get("BACKBONE_BETA", 0.0)
    data["BACK_flex_OTHER"] = d.get("BACKBONE_OTHER", 0.0)

    d = all_interacts["electrostatic_energies"]["counts"]
    data["TotalElec"] = sum(d.values())

    d = all_interacts["hydrogen_bonds"]["counts"]
    data["TotalHbond"] = sum(d.values())

    d = all_interacts["hydrophobics"]["counts"]
    data["Hphobe_contact"] = sum(d.values())

    d = all_interacts["hydrophobics"]["counts"]
    data["Hphobe_contact"] = sum(d.values())

    data["Pi_Pi"] = all_interacts["pi_pi"]["counts"].get("pi_stacking", 0.0)
    data["T-stack"] = all_interacts["pi_pi"]["counts"].get("T_stacking", 0.0)

    d = all_interacts["cat_pi"]["counts"]
    data["Cation-pi"] = sum(d.values())

    d = all_interacts["salt_bridges"]["counts"]
    data["Salt_Bridge"] = sum(d.values())

    return data