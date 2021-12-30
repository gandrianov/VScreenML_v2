# TODO:
# - add params file through residue topology

# import pyrosetta
import pyrosetta

# import utilities
from pyrosetta.rosetta.utility import vector1_double
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.core.pose import renumber_pdbinfo_based_on_conf_chains
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover

# energy terms
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol, gen_bonded
from pyrosetta.rosetta.core.scoring import fa_elec, hbond_bb_sc, hbond_sc

# calculators
from pyrosetta.rosetta.core.pose.metrics import CalculatorFactory 
from pyrosetta.rosetta.core.scoring import calc_per_res_hydrophobic_sasa # metrics 
from pyrosetta.rosetta.core.pose.metrics.simple_calculators import SasaCalculatorLegacy 
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import NumberHBondsCalculator
from pyrosetta.rosetta.protocols.pose_metric_calculators import PackstatCalculator
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import BuriedUnsatisfiedPolarsCalculator


def MinimizeComplex(pose, sfxn, minimizer_type="lbfgs_armijo_nonmonotone",
                    minimizer_tolerance=0.00001, use_nblist=True):
    
    mpose = pose.clone()
    
    options = MinimizerOptions(minimizer_type, minimizer_tolerance)
    options.nblist_auto_update(use_nblist)
    
    mmap = MoveMap()
    mmap.set_chi(True)
    mmap.set_bb(True)
    mmap.set_jump(True)
    
    minimizer = MinMover()
    minimizer.run(mpose, mmap, sfxn, options)
    
    return mpose


def MoveLigandFromActiveSite(pose):
    
    ppose = pose.clone()
    trans_mover = RigidBodyTransMover(ppose, ppose.num_jump())
    trans_mover.step_size(1000000.0)
    trans_mover.apply(ppose)
    
    return ppose

def ExtractLigandFromPose(pose, res_id):

    ppose = pose.clone()
    ppose_size = ppose.size()
    
    if res_id == 1:
        delete_region(ppose, 2, ppose_size)
    elif res_id == ppose_size:
        delete_region(ppose, 1, res_id - 1)
    else:
        delete_region(ppose, res_id + 1, ppose_size)
        delete_region(ppose, 1, res_id - 1)

    return ppose


def ExtractProteinFromPose(pose, res_id):

    ppose = pose.clone()
    delete_region(ppose, res_id, res_id)

    return ppose


def GetPDBStringFromPose(pose):
    
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    
    return buffer.str()


def GetLigandResidueIndex(pose):

    for i in range(1, pose.size()+1):
        if not pose.residue(i).is_protein():
            return i

    return None

def InitPyRosetta(params):
    pyrosetta.init(f"-beta -extra_res_fa {params}")


def LoadPDB(pdbfile):

    def PreprocessPDB(pdb_string):
        # available letters
        letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" + "0123456789"

        # collect chains with protein from PDB file
        chains = [r[21] for r in pdb_string.split("\n") if r.startswith("ATOM")]
        chains = set(chains)    

        print(letters, chains, set(letters).difference(chains))

        het_chain = list(set(letters).difference(chains))[0]

        new_pdb_string = []

        for r in pdb_string.split("\n"):
            if r.startswith("HETNAM"):
                r = r[:15] + het_chain + r[16:]
            elif r.startswith("HETATM"):
                r = r[:21] + het_chain + r[22:]
            new_pdb_string.append(r)

        return "\n".join(new_pdb_string)


    pdb_string = open(pdbfile, "r").read()
    # preprocess chains ids in pdb
    pdb_string = PreprocessPDB(pdb_string) 
    # create Pose object for PDB
    pose = pyrosetta.Pose()
    # add PDB to the pose
    pose_from_pdbstring(pose, pdb_string)

    return pose


def GetScoringFunction():
    sfxn = get_score_function()
    return sfxn

def GetLigandSasa(bound, protein, ligand, ligand_id, sfxn):

    data = {}
    calc_factory = CalculatorFactory.Instance()
    
    sasa_calc_name = "sasa"
    if not calc_factory.check_calculator_exists(sasa_calc_name):
        sasa_calc = SasaCalculatorLegacy()
        calc_factory.register_calculator(sasa_calc_name, sasa_calc)

    bound_sasa = float(bound.print_metric(sasa_calc_name, "total_sasa"))
    protein_sasa = float(protein.print_metric(sasa_calc_name, "total_sasa"))
    ligand_sasa = float(ligand.print_metric(sasa_calc_name, "total_sasa"))    
    
    data["total_exposed_sasa"] = 1 - ( ( (protein_sasa + ligand_sasa) - bound_sasa) * 0.5 / ligand_sasa)
    
    probe_radius = 1.4
    complex_rsd_sasa = vector1_double()
    separated_rsd_sasa = vector1_double()
    complex_rsd_hsasa = vector1_double()
    separated_rsd_hsasa = vector1_double()
    
    calc_per_res_hydrophobic_sasa(bound, complex_rsd_sasa, complex_rsd_hsasa, probe_radius, False)
    calc_per_res_hydrophobic_sasa(protein, separated_rsd_sasa, separated_rsd_hsasa, probe_radius, False)

    complex_sasa = 0.0
    complex_polar_sasa = 0.0
    complex_hydrophobic_sasa = 0.0
    separated_sasa = 0.0
    separated_polar_sasa = 0.0
    separated_hydrophobic_sasa = 0.0
    
    for j in range(1, bound.size()+1):
        if j == ligand_id:
            continue
        complex_sasa += complex_rsd_sasa[j]
        complex_hydrophobic_sasa += complex_rsd_hsasa[j]
        separated_sasa += separated_rsd_sasa[j]
        separated_hydrophobic_sasa += separated_rsd_hsasa[j]
        
    complex_polar_sasa = complex_sasa - complex_hydrophobic_sasa
    separated_polar_sasa = separated_sasa - separated_hydrophobic_sasa

    interface_sasa = abs(complex_sasa - separated_sasa)

    if interface_sasa != 0.0:
        data["interface_hydrophobic_sasa"] = abs(complex_hydrophobic_sasa - separated_hydrophobic_sasa) / interface_sasa
        data["interface_polar_sasa"] = abs(complex_polar_sasa - separated_polar_sasa) / interface_sasa
    else:
        data["interface_hydrophobic_sasa"] = 0.0
        data["interface_polar_sasa"] = 0.0

    return data

def GetEnergyTerms(bound, unbond, sfxn):
    data = {}
    
    bou_emap = bound.energies().total_energies()
    ubo_emap = unbond.energies().total_energies()
    
    data["interaction_score"] = sfxn(bound) - sfxn(unbond)
    data["fa_atr"] = bou_emap[fa_atr] - ubo_emap[fa_atr]
    data["fa_rep"] = bou_emap[fa_rep] - ubo_emap[fa_rep]
    data["fa_sol"] = bou_emap[fa_sol] - ubo_emap[fa_sol]
    data["fa_elec"] = bou_emap[fa_elec] - ubo_emap[fa_elec]
    data["hbond_bb_sc"] = bou_emap[hbond_bb_sc] - ubo_emap[hbond_bb_sc]
    data["hbond_sc"] = bou_emap[hbond_sc] - ubo_emap[hbond_sc]
    data["gen_bonded"] = bou_emap[gen_bonded]
    
    return data
    
def GetMinimizeAppFeatures(bound, unbound):
    data = {}
    
    calc_factory = CalculatorFactory.Instance()

    sasa_calc_name = "sasa"
    if not calc_factory.check_calculator_exists(sasa_calc_name):
        sasa_calc = SasaCalculatorLegacy()
        calc_factory.register_calculator(sasa_calc_name, sasa_calc)

    bound_sasa = float(bound.print_metric(sasa_calc_name, "total_sasa"))
    unbound_sasa = float(unbound.print_metric(sasa_calc_name, "total_sasa"))
    data["total_bsa"] = unbound_sasa - bound_sasa

    hbond_calc_name = "hbond"
    if not calc_factory.check_calculator_exists(hbond_calc_name):
        hb_calc = NumberHBondsCalculator()
        calc_factory.register_calculator(hbond_calc_name, hb_calc)

    bound_hb = float(bound.print_metric(hbond_calc_name, "all_Hbonds"))
    unbound_hb = float(unbound.print_metric(hbond_calc_name, "all_Hbonds"))
    data["interface_hb"] = bound_hb - unbound_hb

    packstat_calc_name = "packstat"
    if not calc_factory.check_calculator_exists(packstat_calc_name):
        packstat_calc = PackstatCalculator()
        calc_factory.register_calculator(packstat_calc_name, packstat_calc)

    # packstat has problem with PDB's insertion codes
    # so we need to renumber all residues
    renumber_pdbinfo_based_on_conf_chains(bound)
    renumber_pdbinfo_based_on_conf_chains(unbound)

    bound_packstat = float(bound.print_metric(packstat_calc_name, "total_packstat"))
    unbound_packstat = float(unbound.print_metric(packstat_calc_name, "total_packstat"))
    data["total_packstats"] = bound_packstat - unbound_packstat

    burunsat_calc_name = "burunsat"
    if not calc_factory.check_calculator_exists(burunsat_calc_name):
        burunsat_calc = BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name)
        calc_factory.register_calculator(burunsat_calc_name, burunsat_calc)

    bound_unsat = float(bound.print_metric(burunsat_calc_name, "all_bur_unsat_polars"))
    unbound_unsat = float(unbound.print_metric(burunsat_calc_name, "all_bur_unsat_polars"))
    data["interface_unsat"] = bound_unsat - unbound_unsat
    
    return data




