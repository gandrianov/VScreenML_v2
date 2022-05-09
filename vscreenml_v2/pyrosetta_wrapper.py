import pyrosetta

# import utilities
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.core.pose import renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover

def MoveLigandFromActiveSite(pose, residue_id=None):

    """Moves target residue on 1M angstroms"""
    
    unbound_pose = pose.clone()

    if residue_id is None:
        # if residue is not specified, then pick last jump (chain)
        # in most cases last jump is target ligand (depends how PDB of complex was prepared)
        jump_id = unbound_pose.num_jump()
    else:
        # define residue_id belongs to what jumps
        n_residues = pose.size()
        jump = unbound_pose.fold_tree().get_residue_edge(residue_id)
        jump_id = jump.label() #e.is_jump()            
        
    trans_mover = RigidBodyTransMover(unbound_pose, jump_id)
    trans_mover.step_size(1000000.0)
    trans_mover.apply(unbound_pose)

    return unbound_pose


def ExtractResidueFromPose(pose, residue_id):

    """Extracts target residue from pose. Useful for extracting ligand from pose"""

    ligand = pose.clone()
    n_residues = pose.size()
    
    # if target residue is first
    if residue_id == 1:
        delete_region(ligand, 2, n_residues)

    # if target residue is last
    elif residue_id == n_residues:
        delete_region(ligand, 1, n_residues - 1)

    # if target residue is in the middle
    else:
        delete_region(ligand, residue_id + 1, n_residues)
        delete_region(ligand, 1, residue_id - 1)

    return ligand


def DeleteResidueFromPose(pose, residue_id):

    """Deletes target residue from pose """

    cleaned_pose = pose.clone()
    delete_region(cleaned_pose, residue_id, residue_id)

    return cleaned_pose


def GetPDBStringFromPose(pose):

    """Prepares PDB string"""
    
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    
    return buffer.str()


def GetResidueIndex(pose, name="LG1"):

    """Searches target residue by name"""

    for i in range(1, pose.size()+1):
        if pose.residue(i).name() == name:
            return i


def LoadPDBFile(pdbfile):

    """Reads PDB file and converts it pose object"""

    pose = pyrosetta.Pose()

    with open(pdbfile, "r") as pdb:
        pdb = pdb.read()
        pdb = ValidatePDBConsistency(pdb)
        pose_from_pdbstring(pose, pdb)

    return pose


def GetScoringFunction():

    """Wrapper to get scoring function"""

    sfxn = get_score_function()
    return sfxn


def InitPyRosetta(args):
    pyrosetta.init(args)


def RenumberResiduesInPose(pose):

    renumbered_pose = pose.clone()

    for i in range(renumbered_pose.size()):
        renumbered_pose.pdb_info().number(i+1, i+1)
        renumbered_pose.pdb_info().icode(i+1, " ")

    return renumbered_pose


def ValidatePDBConsistency(pdb_rows):

    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    pdb_rows = pdb_rows.split("\n")

    chains = [r[21] for r in pdb_rows if r.startswith("ATOM")]
    chains = set(chains)    

    het_chain = list(set(letters).difference(chains))[0]

    for i, row in enumerate(pdb_rows):
        if row.startswith("HETNAM"):
            row = row[:15] + het_chain + row[16:]
        elif row.startswith("HETATM"):
            row = row[:21] + het_chain + row[22:]

        pdb_rows[i] = row

    pdb_rows = "\n".join(pdb_rows)

    return pdb_rows
    