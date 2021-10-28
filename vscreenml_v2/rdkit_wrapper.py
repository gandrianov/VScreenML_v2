from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3, CalcTPSA
from rdkit.Chem.Lipinski import NumHAcceptors, NumHDonors
from rdkit.Chem.rdFreeSASA import CalcSASA, classifyAtoms
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdchem import AtomMonomerType

def reset_aromaticity(mol, idx):

    atom = mol.GetAtomWithIdx(idx)

    atom.SetIsAromatic(False)

    for bond in atom.GetBonds():
        bond.SetIsAromatic(False)
        if bond.GetBondType() is Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    return mol

def correct_charge(mol, atom_id):
    DEFAULTS = {"O":2, "N":3}

    atom = mol.GetAtomWithIdx(atom_id)

    symbol = atom.GetSymbol()
    if symbol in DEFAULTS:
        bonds = atom.GetBonds()
        num_bonds = len(bonds)
        valences = [bond.GetBondTypeAsDouble() for bond in bonds]
        charge = sum(valences) - DEFAULTS[symbol]
        atom.SetFormalCharge(int(charge))

    return mol


def solve_problem(mol):

    problems = Chem.rdmolops.DetectChemistryProblems(mol)

    for problem in problems:
        try:
            problem_atoms = [problem.GetAtomIdx()]
        except:
            problem_atoms = problem.GetAtomIndices()

        problem_type = problem.GetType()

        if problem_type == "AtomValenceException":
            for idx in problem_atoms:
                mol = correct_charge(mol, idx)
        if problem_type == "AtomKekulizeException":
            # print("AtomKekulizeException", problem_atoms)
            for idx in problem_atoms:
                mol = reset_aromaticity(mol, idx)

        if problem_type == "KekulizeException":
            for idx, _ in enumerate(mol.GetAtoms()):
                mol = correct_charge(mol, idx)

            #print(write_smi(mol, kekuleSmiles=False))

            # print(problem_type, problem_atoms)
            # print(write_smi(mol))
            # print()

    if len(problems) != 0:
        mol = solve_problem(mol)

    Chem.SanitizeMol(mol)
    Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)

    return mol

def parse_bonds(filename):

    bond_types = { "1":Chem.BondType.SINGLE,
                   "2":Chem.BondType.DOUBLE,
                   "3":Chem.BondType.TRIPLE,
                   "4":Chem.BondType.AROMATIC }

    bonds = []

    with open(filename, "r") as f:
        f = f.read().split("\n")

    for row in f:
        if row.startswith("BOND_TYPE"):
            row = row.replace(" RING", "")
            row = [elem for elem in row.split(" ")[1:] if elem != ""]
            row = [row[0], row[1], row[-1].replace("#ORGBND", "")]
            row[-1] = bond_types[row[-1]]
            bonds.append(row)

    return bonds


def ReadLigandPDB(mol, params):

    mol = mol[0:mol.find("CONECT")]

    mol = Chem.MolFromPDBBlock(mol,
                              removeHs=False,
                              proximityBonding=False,
                              sanitize=False)

    atom_names = {}

    for i, atom in enumerate(mol.GetAtoms()):
        atom_name = atom.GetPDBResidueInfo().GetName()
        atom_names[atom_name.replace(" ","")] = i

    bonds = parse_bonds(params)

    mol_editable = Chem.RWMol(mol)

    for bond in bonds:
        if bond[0] in atom_names and bond[1] in atom_names:
            begin_atom_name = bond[0]
            end_atom_name = bond[1]
            bond_type = bond[2]

            begin_idx = atom_names[begin_atom_name]
            end_idx = atom_names[end_atom_name]

            mol_editable.AddBond(begin_idx, end_idx, order=bond_type)
            mol = mol_editable.GetMol()
    
    mol = solve_problem(mol)

    for i, _ in enumerate(mol.GetAtoms()):
        mol = correct_charge(mol, i)

    Chem.SanitizeMol(mol)

    for i, atom in enumerate(mol.GetAtoms()):
        atom.UpdatePropertyCache()

    return mol

def CalcRDKitFeatures(mol):

    for a in mol.GetAtoms():
        info = a.GetMonomerInfo()
        if info:
            info.SetMonomerType(AtomMonomerType.UNKNOWN)

    data = {}

    Chem.SanitizeMol(mol)
    
    data["fsp3"] = CalcFractionCSP3(mol)
    data["acceptorcount"] = NumHAcceptors(mol)
    data["donorcount"] = NumHDonors(mol)
    data["logp"] = MolLogP(mol)
    data["polarsurfacearea"] = CalcTPSA(mol)

    AllChem.EmbedMolecule(mol)
    radii = classifyAtoms(mol)
    data["vdwsa"] = CalcSASA(mol, radii)
    
    return data



def ReadProteinPDB(pdb_string):

    # templates for assigning bonds for input protein structure

    class Aminoacid:
        def __init__(self, name):
            self.name = name
            self.atoms = {}

    bond_types = {1:Chem.BondType.SINGLE,
                  2:Chem.BondType.DOUBLE,
                  3:Chem.BondType.TRIPLE,
                  4:Chem.BondType.AROMATIC}

    aas_structure = {
        "GLY":{"N" : [("CA", 1), ("H",1)],
               "CA": [("C", 1), ("1HA", 1), ("2HA", 1)],
               "C" : [("O",2)]},
        "ALA":{"N":  [("CA", 1),("H",1)],
               "CA": [("CB",1),("C",1),("HA",1)],
               "C":  [("O",2)],
               "CB": [("1HB",1), ("2HB",1), ("3HB",1)]},
        "ARG":{"N":  [("CA", 1),("H",1)],
               "CA": [("CB",1), ("C",1),("HA",1)],
               "C":  [("O",2)],"CB":[("CG",1),("1HB",1), ("2HB",1)],
               "CG": [("CD",1),("1HG",1),("2HG",1)],
               "CD": [("NE",1),("1HD",1),("2HD",1)],
               "NE": [("CZ",1), ("HE",1)],
               "CZ": [("NH1",2), ("NH2",1)],
               "NH1":[("1HH1",1), ("2HH1",1)],
               "NH2":[("2HH2",1), ("1HH2",1)]},
        "ASN":{"N":  [("CA", 1),("H",1)],
               "CA": [("CB",1), ("C",1), ("HA",1)], "C":[("O",2)], 
               "CB": [("CG",1), ("1HB",1), ("2HB",1)], 
               "CG": [("OD1",2), ("ND2",1)], "ND2":[("1HD2",1), ("2HD2",1)]},
        "ASP":{"N":  [("CA", 1), ("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],
               "C":  [("O",2)],
               "CB": [("CG",1),("1HB",1), ("2HB",1)], "CG": [("OD1",2), ("OD2",1)]},
        "GLU":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1),("C",1), ("HA",1)],
               "C":  [("O",2)],
               "CB": [("CG",1),("1HB",1),("2HB",1)],
               "CG": [("CD",1),("1HG",1),("2HG",1)],
               "CD": [("OE1",2),("OE2",1)]},
        "GLN":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],
               "C":  [("O",2)],
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)],
               "CG": [("CD",1), ("1HG",1), ("2HG",1)] ,
               "CD": [("OE1",2), ("NE2",1)], 
               "NE2":[("1HE2",1), ("2HE2",1)]},
        "HIS":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("ND1",1), ("CD2",2)], 
               "ND1":[("CE1",2)], 
               "CE1":[("NE2",1), ("HE1",1)], 
               "NE2":[("CD2",1), ("HE2",1)], "CD2":[("HD2",1)]},
        "LEU":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],
               "C":  [("O",2)],
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)],
               "CG": [("CD1",1), ("CD2",1), ("HG",1)],
               "CD2":[("1HD2",1), ("2HD2",1), ("3HD2",1)],
               "CD1":[("1HD1",1), ("2HD1",1), ("3HD1",1)]},
        "LYS":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],
               "C":  [("O",2)],
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)],
               "CG": [("CD",1), ("1HG",1), ("2HG",1)],
               "CD": [("CE",1), ("1HD",1), ("2HD",1)],
               "CE": [("NZ",1), ("1HE",1), ("2HE",1)], 
               "NZ": [("1HZ",1), ("2HZ",1), ("3HZ",1)]},        
        "MET":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("SD",1), ("1HG",1), ("2HG",1)],
               "SD": [("CE",1)], 
               "CE": [("1HE",1), ("2HE",1), ("3HE",1)]},
        "PHE":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("CD1",2),("CD2",1)],
               "CD1":[("CE1",1), ("HD1",1)],
               "CE1":[("CZ",2), ("HE1",1)],
               "CZ": [("CE2",1), ("HZ",1)], 
               "CE2":[("CD2",1), ("HE2",1)], 
               "CD2":[("HD2",1)]},
        "PRO":{"N":  [("CA", 1),("CD",1)],
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("CD",1), ("1HG",1), ("2HG",1)], 
               "CD": [("1HD",1), ("2HD",1)]},
        "TRP":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("CD1",2), ("CD2",1)], 
               "CD1":[("NE1",1), ("HD1",1)], 
               "NE1":[("CE2",1), ("HE1",1)], 
               "CE2":[("CD2",2), ("CZ2",1)], 
               "CZ2":[("CH2",2), ("HZ2",1)],
               "CH2":[("CZ3",1), ("HH2",1)],
               "CZ3":[("CE3",2), ("HZ3",1)], 
               "CE3":[("CD2",1), ("HE3",1)]},        
        "TYR":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG",1),  ("1HB",1), ("2HB",1)], 
               "CG": [("CD1",2), ("CD2",1)], 
               "CD1":[("CE1",1), ("HD1",1)], 
               "CD2":[("HD2",1)], 
               "CE1":[("CZ",2), ("HE1",1)], 
               "CZ": [("CE2",1), ("OH",1)], 
               "CE2":[("CD2",2), ("HE2",1)], 
               "OH": [("HH",1)]},
        "CYS":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("SG",1),  ("1HB",1), ("2HB",1)], 
               "SG": [("HG",1)]},
        "SER":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("OG",1),  ("1HB",1), ("2HB",1)], 
               "OG": [("HG",1)]},
        "THR":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("OG1",1), ("CG2",1), ("HB",1)],  
               "OG1":[("HG1",1)], 
               "CG2":[("1HG2",1), ("2HG2",1), ("3HG2",1)]},
        "ILE":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG1",1), ("CG2",1), ("HB",1)],  
               "CG1":[("CD1",1),  ("1HG1",1), ("2HG1",1)], 
               "CD1":[("1HD1",1), ("2HD1",1), ("3HD1",1)], 
               "CG2":[("1HG2",1), ("2HG2",1), ("3HG2",1)]},
        "VAL":{"N":  [("CA", 1),("H",1)], 
               "CA": [("CB",1), ("C",1), ("HA",1)],   
               "C":  [("O",2)], 
               "CB": [("CG1",1), ("CG2",1), ("HB",1)],  
               "CG1":[("1HG1",1), ("2HG1",1), ("3HG1",1)], 
               "CG2":[("1HG2",1), ("2HG2",1), ("3HG2",1)]}
    }


    # drop any bonds annotations
    pdb_string = pdb_string[0:pdb_string.find("CONECT")]

    # read PDB to RDKit.Mol
    mol = Chem.MolFromPDBBlock(pdb_string,
                               removeHs=False,
                               proximityBonding=False,
                               sanitize=False)


    # parse residue names and number in RDKit.Mol
    atoms = {}

    for i, atom in enumerate(mol.GetAtoms()):

        atom_name = atom.GetPDBResidueInfo().GetName().strip()
        residue_num = int(atom.GetPDBResidueInfo().GetResidueNumber())
        residue_name = atom.GetPDBResidueInfo().GetResidueName()

        if residue_num not in atoms:
            atoms[residue_num] = Aminoacid(residue_name)

        atoms[residue_num].atoms[atom_name] = i


    # create writeable RDKit.Mol
    rwmol = Chem.RWMol(mol)
    
    # assign bonds based on amino acid templates
    for i, (residue_num, aminoacid) in enumerate(atoms.items()):
        residue_name = aminoacid.name
        for atom_name, mol_index in aminoacid.atoms.items():
            # connect backbone N of aa with backbone C of aa - 1 
            if atom_name == "N":
                # need for skipping first amino acid in the sequence
                if residue_num - 1 in atoms:
                    rwmol.AddBond(mol_index, 
                                  atoms[residue_num - 1].atoms["C"],
                                  Chem.BondType.SINGLE)

            # check if amino acid in templates
            if atom_name in aas_structure[residue_name]:
                # iterate thru all atoms in template
                for a, bond_order in aas_structure[residue_name][atom_name]:
                    # skip extra terminal H (usually Rosetta have NH3+ 
                    # in the beginning)
                    
                    if i == 0 and a == "H":
                        continue
                    

                    if a not in atoms[residue_num].atoms:
                        continue
                    
                    rwmol.AddBond(mol_index, atoms[residue_num].atoms[a],
                                  bond_types[bond_order])


    # connect terminal N and Hs
    for i in range(3):
        rwmol.AddBond(atoms[min(atoms)].atoms["N"],
                      atoms[min(atoms)].atoms[f"{i+1}H"],
                      Chem.BondType.SINGLE)

    # connect terminal C and extra oxigen
    rwmol.AddBond(atoms[max(atoms)].atoms["C"], atoms[max(atoms)].atoms["OXT"], Chem.BondType.SINGLE)

    # get read-only RDKit.Mol
    mol = rwmol.GetMol()
    mol = solve_problem(mol)
    mol.UpdatePropertyCache()
   
    return mol

