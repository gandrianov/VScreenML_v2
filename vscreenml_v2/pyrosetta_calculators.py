from pyrosetta_wrapper import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.utility import vector1_double
from pyrosetta.rosetta.core.pose.metrics import CalculatorFactory 
from pyrosetta.rosetta.core.scoring import calc_per_res_hydrophobic_sasa 
from pyrosetta.rosetta.core.pose.metrics.simple_calculators import SasaCalculatorLegacy 
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import NumberHBondsCalculator
from pyrosetta.rosetta.protocols.pose_metric_calculators import PackstatCalculator
from pyrosetta.rosetta.protocols.simple_pose_metric_calculators import BuriedUnsatisfiedPolarsCalculator

class SasaCalculator:

    def __init__(self, sasa_calc_name="sasa"):

        self.sasa_calc = SasaCalculatorLegacy()
        self.sasa_calc_name = sasa_calc_name

        self.calc_factory = CalculatorFactory.Instance()

        if not self.calc_factory.check_calculator_exists(self.sasa_calc_name):
            self.calc_factory.register_calculator(self.sasa_calc_name, self.sasa_calc)

    def GetTotalExposedSasa(self, bound_pose, unbound_pose, ligand_pose):

        bound_sasa = bound_pose.print_metric(self.sasa_calc_name, "total_sasa")
        bound_sasa = float(bound_sasa)
        
        unbound_sasa = unbound_pose.print_metric(self.sasa_calc_name, "total_sasa")
        unbound_sasa = float(unbound_sasa)

        ligand_sasa = ligand_pose.print_metric(self.sasa_calc_name, "total_sasa")
        ligand_sasa = float(ligand_sasa)

        return 1 - ( ( (unbound_sasa + ligand_sasa) - bound_sasa) * 0.5 / ligand_sasa)

    def GetTotalBSA(self, bound_pose, unbound_pose):

        bound_sasa = bound_pose.print_metric(self.sasa_calc_name, "total_sasa")
        bound_sasa = float(bound_sasa)
        unbound_sasa = unbound_pose.print_metric(self.sasa_calc_name, "total_sasa")
        unbound_sasa = float(unbound_sasa)

        return unbound_sasa - bound_sasa
    
    def GetInterfaceSasa(self, bound_pose, unbound_pose, ligand_id, probe_radius=1.4):

        complex_rsd_sasa = vector1_double()
        separated_rsd_sasa = vector1_double()

        complex_rsd_hsasa = vector1_double()
        separated_rsd_hsasa = vector1_double()

        calc_per_res_hydrophobic_sasa(bound_pose, complex_rsd_sasa, complex_rsd_hsasa, probe_radius, False)
        calc_per_res_hydrophobic_sasa(unbound_pose, separated_rsd_sasa, separated_rsd_hsasa, probe_radius, False)

        complex_sasa = sum(complex_rsd_sasa) - complex_rsd_sasa[ligand_id]
        complex_hydrophobic_sasa = sum(complex_rsd_hsasa) - complex_rsd_hsasa[ligand_id]

        separated_sasa = sum(separated_rsd_sasa) - separated_rsd_sasa[ligand_id]
        separated_hydrophobic_sasa = sum(separated_rsd_hsasa) - separated_rsd_hsasa[ligand_id]

        complex_polar_sasa = complex_sasa - complex_hydrophobic_sasa
        separated_polar_sasa = separated_sasa - separated_hydrophobic_sasa

        interface_sasa = abs(complex_sasa - separated_sasa)

        if interface_sasa != 0.0:
            interface_hydrophobic_sasa = abs(complex_hydrophobic_sasa - separated_hydrophobic_sasa) / interface_sasa
            interface_polar_sasa = abs(complex_polar_sasa - separated_polar_sasa) / interface_sasa
        else:
            interface_hydrophobic_sasa = 0.0
            interface_polar_sasa = 0.0

        return interface_hydrophobic_sasa, interface_polar_sasa


class EnergyCalculator:
    def __init__(self):
        self.sfxn = GetScoringFunction()

    def CalcEnergies(self, bound_pose, unbound_pose):
        self.bound_total_energy = self.sfxn(bound_pose)
        self.unbound_total_energy = self.sfxn(unbound_pose)
        self.bound_energies = bound_pose.energies().total_energies()
        self.unbound_energies = unbound_pose.energies().total_energies()

    def GetInteractionScore(self):
        return self.bound_total_energy - self.unbound_total_energy

    def GetFaAtrInteraction(self):
        return self.bound_energies[fa_atr] - self.unbound_energies[fa_atr]

    def GetFaRepInteraction(self):
        return self.bound_energies[fa_rep] - self.unbound_energies[fa_rep]

    def GetFaSolInteraction(self):
        return self.bound_energies[fa_sol] - self.unbound_energies[fa_sol]

    def GetFaElecInteraction(self):
        return self.bound_energies[fa_elec] - self.unbound_energies[fa_elec]

    def GetHBondBbScInteraction(self):
        return self.bound_energies[hbond_bb_sc] - self.unbound_energies[hbond_bb_sc]

    def GetHBondScInteraction(self):
        return self.bound_energies[hbond_sc] - self.unbound_energies[hbond_sc]

    def GetGenBonded(self):
        return self.bound_energies[gen_bonded] 


class PoseMetricCalculator:
    def __init__(self, sasa_calc_name="sasa", hbond_calc_name="hbond", packstat_calc_name="packstat", burunsat_calc_name="burunsat", packstat_oversample=100):

        self.sasa_calc = SasaCalculatorLegacy()
        self.sasa_calc_name = sasa_calc_name

        self.hb_calc = NumberHBondsCalculator()
        self.hbond_calc_name = hbond_calc_name

        self.packstat_calc = PackstatCalculator(oversample=packstat_oversample, remove_nonprotein_res=False)
        self.packstat_calc_name = packstat_calc_name

        self.burunsat_calc = BuriedUnsatisfiedPolarsCalculator(self.sasa_calc_name, self.hbond_calc_name)
        self.burunsat_calc_name = burunsat_calc_name

        self.calc_factory = CalculatorFactory.Instance()

        if not self.calc_factory.check_calculator_exists(self.sasa_calc_name):
            self.calc_factory.register_calculator(self.sasa_calc_name, self.sasa_calc)

        if not self.calc_factory.check_calculator_exists(self.hbond_calc_name):
            self.calc_factory.register_calculator(self.hbond_calc_name, self.hb_calc)

        if not self.calc_factory.check_calculator_exists(self.packstat_calc_name):  
            self.calc_factory.register_calculator(self.packstat_calc_name, self.packstat_calc)
            
        if not self.calc_factory.check_calculator_exists(self.burunsat_calc_name):
            self.calc_factory.register_calculator(self.burunsat_calc_name, self.burunsat_calc)


    def GetHBInterface(self, bound_pose, unbound_pose):
        bound_hb = bound_pose.print_metric(self.hbond_calc_name, "all_Hbonds")
        bound_hb = float(bound_hb)
        unbound_hb = unbound_pose.print_metric(self.hbond_calc_name, "all_Hbonds")
        unbound_hb = float(unbound_hb)
        
        return bound_hb - unbound_hb

    def GetTotalPackStat(self, bound_pose, unbound_pose):

        bound_packstat = bound_pose.print_metric(self.packstat_calc_name, "total_packstat")
        unbound_packstat = unbound_pose.print_metric(self.packstat_calc_name, "total_packstat")
        
        bound_packstat = float(bound_packstat)
        unbound_packstat = float(unbound_packstat)

        return bound_packstat - unbound_packstat

    def GetInterfaceUnsat(self, bound_pose, unbound_pose):

    	# calculates number of polar atoms that became buried during binding with ligand
    	# detailed explanation - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7971855/pdf/pcbi.1008061.pdf

        bound_unsat = bound_pose.print_metric(self.burunsat_calc_name, "all_bur_unsat_polars")
        bound_unsat = float(bound_unsat)

        unbound_unsat = unbound_pose.print_metric(self.burunsat_calc_name, "all_bur_unsat_polars")
        unbound_unsat = float(unbound_unsat)

        return bound_unsat - unbound_unsat


