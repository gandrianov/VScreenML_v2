from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3, CalcTPSA
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
from rdkit.Chem.rdFreeSASA import CalcSASA, classifyAtoms
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdchem import AtomMonomerType

class RDKitCalculator:
    def __init__(self):
        pass

    def CalcFCsp3(self, mol):
        return CalcFractionCSP3(mol)

    def CalcNumHAcceptors(self, mol):
        return NumHAcceptors(mol)

    def CalcNumHDonors(self, mol):
        return NumHDonors(mol)

    def CalcMolLogP(self, mol):
        return MolLogP(mol)

    def CalcTPSA(self, mol):
        return CalcTPSA(mol)

    def CalcVdwSA(self, mol):
        
        mol_copy = Chem.Mol(mol)

        for a in mol_copy.GetAtoms():
            info = a.GetMonomerInfo()
            if info:
                info.SetMonomerType(AtomMonomerType.UNKNOWN)

        radii = classifyAtoms(mol_copy)
        succ = AllChem.EmbedMolecule(mol_copy)

        if succ != -1:
            return CalcSASA(mol_copy, radii, confIdx=-1)
        else:
            return None

