#!/usr/bin/env python
import argparse
import pyrosetta
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver
from pyrosetta.rosetta.core.scoring.constraints import add_fa_constraints_from_cmdline

# Basic usage of minimizer:
# python pyrosetta_minimizer.py -ligand <pdb> -params <params> -protein <pdb> -output <pdb>
# Usage of minimizer with constraints:
# python pyrosetta_minimizer.py -ligand <pdb> -params <params> -protein <pdb> -output <pdb> -constraints <cst>
# Usage of minimizer with pymol observer:
# python pyrosetta_minimizer.py -ligand <pdb> -params <params> -protein <pdb> -output <pdb> -pymol-observer true

def args():

    """Wrapper for input arguments"""
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-ligand", type=str, required=True,
                        help="PDB file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-params", type=str, required=True,
                        help="Params file of ligand generated by generic_potential/mol2genparams.py")
    parser.add_argument("-protein", type=str, required=True,
                        help="PDB file of input protein structure")
    parser.add_argument("-output", type=str, required=True,
                        help="File name for minimized protein-ligand complex")
    parser.add_argument("-init-args", type=str, default="-in:auto_setup_metals -packing:no_optH",
                        help="Argument for initializing PyRosetta state")
    parser.add_argument("-constraints", type=str, required=False,
                        help="File for adding AtomPair constraints")
    parser.add_argument("-minimizer_fn", type=str, default="lbfgs_armijo_nonmonotone",
                        help="Name of minimization function")
    parser.add_argument("-minimizer_tolerance", type=float, default=0.000001,
                        help="Tolerance of minimization")
    parser.add_argument("-pymol-observer", type=bool, default=False,
                        help="Send minimization frame to PyMOL")

    args = parser.parse_args()

    return args


def MinimizeComplex(pose, sfxn, **kwargs):

    """Minimizes merged protein-ligand complex"""
    
    # all manipulations were done on cloned pose
    minimized_pose = pose.clone()

    if args.constraints is not None:
        # add constransts to pose
        add_fa_constraints_from_cmdline(minimized_pose, sfxn)
        # add weight to score function
        sfxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 1.0)
        print(sfxn)
    
    # create move map 
    mmap = MoveMap()
    # allow flexibility on side chains
    mmap.set_chi(True)
    # allow flexibility on back bones
    mmap.set_bb(True)
    # allow flexibility on chains
    mmap.set_jump(True)
    
    # minimizer
    minimizer = MinMover(movemap_in=mmap,
                         scorefxn_in=sfxn,
                         min_type_in=kwargs["minimizer_fn"],
                         tolerance_in=kwargs["minimizer_tolerance"],
                         use_nb_list_in=True,
                         deriv_check_in=False,
                         deriv_check_verbose_in=False)

    minimizer.max_iter(10000)


    if args.pymol_observer == True:
        pyobs = AddPyMOLObserver(minimized_pose, keep_history=True, update_interval=5)
        pyobs.attach(minimized_pose)

    # minimize pose
    minimizer.apply(minimized_pose)
    
    return minimized_pose


def MoveLigandFromActiveSite(pose):

    """Creates unbound pose by moving ligand from protein on 1M angstrom"""
    
    # all manipulations were done on cloned pose
    unbound_pose = pose.clone()

    # create mover focused on last chain (usually ligand of interest)
    trans_mover = RigidBodyTransMover(unbound_pose, unbound_pose.num_jump())
    trans_mover.step_size(1000000.0)

    # move ligand
    trans_mover.apply(unbound_pose)
    
    return unbound_pose

def GetPDBStringFromPose(pose):

    """Creates PDB string from pose"""
    
    # buffer for writing
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    
    return buffer.str()


def MergePDBsToOneComplex(protein_file, ligand_file):

    """Merges protein and ligand pdbs into one text"""
    
    protein_pdb = open(args.protein, "r").read()
    ligand_pdb = open(args.ligand, "r").read()

    bound_pdb = protein_pdb + ligand_pdb

    return bound_pdb



if __name__ == "__main__":

    args = args()

    args.init_args += f" -beta -extra_res_fa {args.params}"

    if args.constraints is not None:
        args.init_args += f" -constraints:cst_fa_file {args.constraints}"

    pyrosetta.init(args.init_args)

    # scoring function
    sfxn = get_score_function()
    # create pose for minimization
    bound_pose = pyrosetta.rosetta.core.pose.Pose()
    # create merged pdb file of ligand and protein
    bound_pdb = MergePDBsToOneComplex(args.ligand, args.protein)
    # put merged pdb to pose
    pose_from_pdbstring(bound_pose, pdbcontents=bound_pdb)
    
    # minimize bound pose
    bound_pose = MinimizeComplex(bound_pose, sfxn, **vars(args))
    # move ligand from protein on 1M angstroms to some direction
    unbound_pose = MoveLigandFromActiveSite(bound_pose)
    # calculate interaction energy
    interaction_energy = sfxn(bound_pose) - sfxn(unbound_pose)

    # convert pose to pdb string
    bound_pdb = GetPDBStringFromPose(bound_pose)
    # add interaction score to end of pdb string
    bound_pdb += f"\nInteraction energy: {interaction_energy}\n" 

    # write to file
    with open(args.output, "w") as fwr:
        fwr.write(bound_pdb)
        fwr.close()
