from oddt.toolkits.rdk import Molecule
from oddt.scoring.descriptors import close_contacts_descriptor

LIGAND_ATOMS = [6, 7, 8, 9, 15, 16, 17, 35, 53]
PROTEIN_ATOMS = [6, 7, 8, 16]

class RFScoreCalculator:
    def __init__(self, ligand_atoms=LIGAND_ATOMS, protein_atoms=PROTEIN_ATOMS, cutoff=12):
        self.cutoff = cutoff
        self.ligand_atoms = ligand_atoms
        self.protein_atoms = protein_atoms


    def GetScore(self, ligand_mol, protein_mol):
        ligand_mol = Molecule(ligand_mol)
        protein_mol = Molecule(protein_mol)

        descriptors = close_contacts_descriptor(protein_mol,
                                                cutoff=self.cutoff, 
                                                protein_types=self.protein_atoms, 
                                                ligand_types=self.ligand_atoms, 
                                                aligned_pairs=False)

        features = descriptors.build(ligand_mol).tolist()[0]
        return dict(zip(descriptors.titles, features))