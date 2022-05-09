from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import AtomMonomerType

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

BOND_TYPES = {1:Chem.BondType.SINGLE,
              2:Chem.BondType.DOUBLE,
              3:Chem.BondType.TRIPLE,
              4:Chem.BondType.AROMATIC}

AMINOACID_CONNECTIVITY = {"GLY":{"N" : [("CA", 1), ("H",1)],
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
                                 "ND1":[("CE1",2), ("HD1",1)], 
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
                                 "CG2":[("1HG2",1), ("2HG2",1), ("3HG2",1)]}}


def ResetAtomAromaticity(mol, atom_id):

    """Reset aromaticity flags on target atom and on bonds attached to the atom"""

    non_aromatic_mol = Chem.Mol(mol)

    atom = non_aromatic_mol.GetAtomWithIdx(atom_id)

    atom.SetIsAromatic(False)

    for bond in atom.GetBonds():
        bond.SetIsAromatic(False)
        if bond.GetBondType() is Chem.rdchem.BondType.AROMATIC:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    return non_aromatic_mol


def AddCorrectCharge(mol, atom_id):

    """Replace formal charge based on valence"""

    DEFAULT_CHARGE = {"O":2, "N":3}

    recharged_mol = Chem.Mol(mol)

    atom = recharged_mol.GetAtomWithIdx(atom_id)
    symbol = atom.GetSymbol()

    if symbol in DEFAULT_CHARGE:
        bonds = atom.GetBonds()
        valences = [bond.GetBondTypeAsDouble() for bond in bonds]
        charge = sum(valences) - DEFAULT_CHARGE[symbol]
        atom.SetFormalCharge(int(charge))

    return recharged_mol


def SolveMolProblem(mol):

    """Detects possible chemistry problems and tries recursively solve them"""

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        problem = e.cause
        problem_type = problem.GetType()

        if problem_type == "AtomValenceException":
            mol = AddCorrectCharge(mol, problem.GetAtomIdx())

        elif problem_type == "AtomKekulizeException":
            mol = ResetAtomAromaticity(mol, problem.GetAtomIdx())

    for atom in mol.GetAtoms():
        try:
            atom.UpdatePropertyCache(strict=True)
        except Exception as e:
            problem = e.cause
            problem_type = problem.GetType()

            if problem_type == "AtomValenceException":
                mol = AddCorrectCharge(mol, problem.GetAtomIdx())

            elif problem_type == "AtomKekulizeException":
                mol = ResetAtomAromaticity(mol, problem.GetAtomIdx())

    problems = Chem.rdmolops.DetectChemistryProblems(mol)

    for problem in problems:
        try:
            # problem related only to one atom
            problem_atoms = [problem.GetAtomIdx()]
        except:
            # problem related to multiple atoms
            problem_atoms = problem.GetAtomIndices()

        problem_type = problem.GetType()

        # atom = mol.GetAtomWithIdx(problem_atoms[0])
        # atom_name = atom.GetPDBResidueInfo().GetName().strip()
        # residue_num = str(atom.GetPDBResidueInfo().GetResidueNumber())
        # residue_name = atom.GetPDBResidueInfo().GetResidueName().strip()
        # chain = atom.GetPDBResidueInfo().GetChainId().strip()

        # print(atom_name, residue_num, residue_name, chain)

        if problem_type == "AtomValenceException":
            for idx in problem_atoms:
                mol = AddCorrectCharge(mol, idx)

        if problem_type == "KekulizeException":
            for idx, _ in enumerate(mol.GetAtoms()):
                mol = AddCorrectCharge(mol, idx)

        if problem_type == "AtomKekulizeException":
            for idx in problem_atoms:
                mol = ResetAtomAromaticity(mol, idx)

    # run again to check if problems disappear or not
    if len(problems) != 0:
        mol = SolveMolProblem(mol)

    # clean mol and assign correct aromaticity
    mol.UpdatePropertyCache(strict=True)
    Chem.SanitizeMol(mol)

    # Chem.rdmolops.SetAromaticity(mol, model=Chem.rdmolops.AromaticityModel.AROMATICITY_MDL)
    # Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
    
    return mol


def ReadBondsOrderFromParams(params):

    """Restore bonds order from params file"""

    BOND_TYPES = {"1":Chem.BondType.SINGLE,
                  "2":Chem.BondType.DOUBLE,
                  "3":Chem.BondType.TRIPLE,
                  "4":Chem.BondType.AROMATIC}

    bonds = []

    for row in params.split("\n"):
        if row.startswith("BOND_TYPE"):
            row = row.replace(" RING", "")
            row = [elem for elem in row.split(" ")[1:] if elem != ""]
            bond = {"A1":row[0], "A2":row[1], "ORDER":BOND_TYPES[row[-1].replace("#ORGBND", "")]}
            bonds.append(bond)

    return bonds


def GetPDBAtomNames(mol):

    """Extracts PDB atoms names"""

    names = {}

    for i, atom in enumerate(mol.GetAtoms()):
        name = atom.GetPDBResidueInfo().GetName()
        names[name.strip()] = i

    return names


def ReadPDBParamsFile(pdb, params):

    """Read PDB and params file to build RDKit.Chem.Mol object"""
    mol = Chem.MolFromPDBBlock(pdb, removeHs=False, proximityBonding=False, sanitize=False)

    atom_names = GetPDBAtomNames(mol)
    bonds = ReadBondsOrderFromParams(params)

    mol_editable = Chem.RWMol(mol)

    # delete old bonds
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        mol_editable.RemoveBond(a1, a2)

    # add new bonds from params file

    for b in bonds:
        a1, a2, bond_order = b["A1"], b["A2"], b["ORDER"]

        if a1 not in atom_names or a2 not in atom_names:
            continue

        a1_idx = atom_names[a1]
        a2_idx = atom_names[a2]
        mol_editable.AddBond(a1_idx, a2_idx, order=bond_order)
    
    # clean if it has any chemical problems
    mol = mol_editable.GetMol()
    mol = SolveMolProblem(mol)

    return mol


def ReadPDBFile(pdb):

    mol = Chem.MolFromPDBBlock(pdb, removeHs=False, proximityBonding=False, sanitize=False)

    # parse residue names and number in RDKit.Mol
    structure = {}

    for i, atom in enumerate(mol.GetAtoms()):

        atom_name = atom.GetPDBResidueInfo().GetName().strip()
        residue_num = str(atom.GetPDBResidueInfo().GetResidueNumber())
        residue_name = atom.GetPDBResidueInfo().GetResidueName().strip()
        chain = atom.GetPDBResidueInfo().GetChainId().strip()

        if chain not in structure:
            structure[chain] = {}

        if residue_name + residue_num not in structure[chain]:
            structure[chain][residue_name + residue_num] = {}

        if atom_name not in structure[chain][residue_name + residue_num]:
            structure[chain][residue_name + residue_num][atom_name] = i

    # create writeable RDKit.Mol
    mol_editable = Chem.RWMol(mol)

    # delete old bonds
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()

        mol_editable.RemoveBond(a1, a2)
    
    for chain, residues in structure.items():

        previous_residue = None
        
        for residue, atoms in residues.items():

            terminal_residue = "1H" in atoms

            if residue[:3] not in AMINOACID_CONNECTIVITY:
                print(f"{residue} is not presented in canonical amino acids set")
                continue

            for atom, i in atoms.items():

                # Add bonds for capped hydrogens for N-terminal amino acid
                if atom in ["1H", "2H", "3H"]:
                    mol_editable.AddBond(i, structure[chain][residue]["N"], BOND_TYPES[1])

                # connect terminal C and extra oxigen
                if atom == "OXT":
                    mol_editable.AddBond(i, structure[chain][residue]["C"], BOND_TYPES[1])
                    
                # connect backbone N of aa with backbone C of aa - 1 
                if atom == "N":
                    if previous_residue is not None:

                        new_residue_number = int(residue[3:])
                        old_residue_number = int(previous_residue[3:])

                        # addd C - N bonding only if numbering is sequential
                        if new_residue_number - old_residue_number == 1:
                            mol_editable.AddBond(i, structure[chain][previous_residue]["C"], BOND_TYPES[1])

                    previous_residue = residue

                if atom in AMINOACID_CONNECTIVITY[residue[:3]]:
                    for a2_name, bond in AMINOACID_CONNECTIVITY[residue[:3]][atom]:
                        # skip default N-hydrogen bonding if amino acid is terminal
                        if a2_name == "H" and terminal_residue:
                            continue

                        if a2_name in structure[chain][residue]:
                            mol_editable.AddBond(i, structure[chain][residue][a2_name], BOND_TYPES[bond])


    # clean mol
    mol = mol_editable.GetMol()
    mol = SolveMolProblem(mol)
    mol.UpdatePropertyCache()
    Chem.SanitizeMol(mol)
   
    return mol

