from oddt.toolkits.rdk import Molecule
from oddt.scoring.descriptors import close_contacts_descriptor

LigandChemElements = [6, 7, 8, 9, 15, 16, 17, 35, 53]
ProteinChemElements = [6, 7, 8, 16]

def GetRFscore(ligand, protein, cutoff=12):

	ligand = Molecule(ligand)
	protein = Molecule(protein)
	
	descriptors = close_contacts_descriptor(protein, cutoff=cutoff, 
											protein_types=ProteinChemElements, 
											ligand_types=LigandChemElements, 
											aligned_pairs=False)

	features = descriptors.build(ligand).tolist()[0]
	return dict(zip(descriptors.titles, features))