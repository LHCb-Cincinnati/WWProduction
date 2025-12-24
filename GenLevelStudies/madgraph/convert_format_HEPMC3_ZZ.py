""" Program to convert HEPMC3 files to a ROOT Tree
"""

# Imports
import sys
import pdb
import pyhepmc
import ROOT
import numpy as np
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at
import PDFReweighter

def get_nonradiative_decay(particle):
    if (particle.children[0].pid == particle.pid):
        return(get_nonradiative_decay(particle.children[0]))
    else:
        return(particle)

def check_parents(particle, bad_particle):
    if particle.parents:
        if (particle == bad_particle):
            return(False)
        else:
            return(check_parents(particle.parents[0], bad_particle))
    else:
        return(True)

def fill_kinematic_array(particle, array):
    mom = particle.momentum
    array[:] = (mom.px, mom.py, mom.pz, mom.pt(), 
                np.sqrt(mom.px *mom.px + mom.py * mom.py + mom.pz * mom.pz), 
                mom.eta(), mom.e, mom.phi(), mom.m(), particle.pid, 
                particle.status)
    return(array)

# NNPDF31 Reweighter
pdfrwgt_NNPDF31 = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="NNPDF31_lo_as_0118"
)
# CT18NLO Reweighter
pdfrwgt_CT18LO = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="CT18LO"
)
# MSHT20LO Reweighter
pdfrwgt_MSHT20LO = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="MSHT20lo_as130"
)

# Files
ifile_name = "/data/home/ganowak/MG5_aMC_v2_9_25/ZZ_LO/Events/run_02_decayed_1/tag_1_pythia8_events.hepmc" 
ofile_name = "Test1.root"
# Setting up Path
ww_path = at.find_WW_path()

# Inputs
Z_id = 23
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
jet_array = np.array([-5,-3,-1])

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:status/F'
wgt_str = "AUX_MUR05_MUF05/F:AUX_MUR05_MUF10/F:AUX_MUR05_MUF20/F:AUX_MUR10_MUF05/F:AUX_MUR10_MUF10/F:AUX_MUR10_MUF20/F:AUX_MUR20_MUF05/F:AUX_MUR20_MUF10/F:AUX_MUR20_MUF20/F"
initial_conditions_str = 'x1/F:id1/F:x2/F:id2/F'
pdfrwgt_str = 'nnpdf31lo/F:ct18lo/F:msht20lo/F'
evt_array = np.array([0], dtype=np.float32)
wgt_array = np.array([0]*9, dtype=np.float32)
Z1_array = np.array([0]*11, dtype=np.float32)
Z1_lminus_array = np.array([0]*11, dtype=np.float32)
Z1_lplus_array = np.array([0]*11, dtype=np.float32)
Z2_array = np.array([0]*11, dtype=np.float32)
Z2_lminus_array = np.array([0]*11, dtype=np.float32)
Z2_lplus_array = np.array([0]*11, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
pdfrwgt_array = np.array([0]*3, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/madgraph/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch(
    'InitialConditions', initial_conditions_array, initial_conditions_str
)
tree.Branch('pdfReweight', pdfrwgt_array, pdfrwgt_str)
tree.Branch('Weights', wgt_array, wgt_str)
tree.Branch('Z1', Z1_array, var_str)
tree.Branch('Z1_lminus', Z1_lminus_array, var_str)
tree.Branch('Z1_lplus', Z1_lplus_array, var_str)
tree.Branch('Z2', Z2_array, var_str)
tree.Branch('Z2_lminus', Z2_lminus_array, var_str)
tree.Branch('Z2_lplus', Z2_lplus_array, var_str)

# Variables
i = 0
counter1 = 0
counter2 = 0

with pyhepmc.open(ifile_name) as f:
    for event in f:
        evt_array[:] = i
        # Find relevant particles
        initial_parton_list = [particle for particle in event.particles if particle.status==21]
        Z_list = [particle for particle in event.particles if particle.pid ==23]
        Z1 = Z_list[0]
        Z_list.reverse()
        for Z in Z_list:
            if check_parents(Z, Z1):
                Z2 = Z
        Z1 = get_nonradiative_decay(Z1)
        Z2 = get_nonradiative_decay(Z2)

        # Find leptonic daughters
        for daughter in Z1.children:
            if daughter.pid in lepton_pid_array:
                Z1_lminus = daughter
            elif daughter.pid in -1*lepton_pid_array:
                Z1_lplus = daughter
        for daughter in Z2.children:
            if daughter.pid in lepton_pid_array:
                Z2_lminus = daughter
            elif daughter.pid in -1*lepton_pid_array:
                Z2_lplus = daughter

        # Get Initial Conditions Info
        if len(initial_parton_list) != 2:
            print("More than 2 initial partons")
            print(initial_parton_list)
        parton1 = initial_parton_list[0]
        parton2 = initial_parton_list[1]
        x1 = parton1.momentum.p3mod() / 6500
        id1 = parton1.pid
        x2 = parton2.momentum.p3mod() / 6500
        id2 = parton2.pid
        nnpdf_reweight = pdfrwgt_NNPDF31.calculate_reweighting_factor(
            x1, x2, id1, id2, sqrt_s=13e3
        )
        ct18lo_reweight = pdfrwgt_CT18LO.calculate_reweighting_factor(
            x1, x2, id1, id2, sqrt_s=13e3
        )
        msht20lo_reweight = pdfrwgt_MSHT20LO.calculate_reweighting_factor(
            x1, x2, id1, id2, sqrt_s=13e3
        )

        # Fill arrays
        pdfrwgt_array[:] = (nnpdf_reweight, ct18lo_reweight, msht20lo_reweight)
        initial_conditions_array[:] = (x1, id1, x2, id2)

        Z1_array = fill_kinematic_array(Z, Z1_array)
        Z1_lminus_array = fill_kinematic_array(Z1_lminus, Z1_lminus_array)
        Z1_lplus_array = fill_kinematic_array(Z1_lplus, Z1_lplus_array)
        Z2_array = fill_kinematic_array(Z, Z2_array)
        Z2_lminus_array = fill_kinematic_array(Z2_lminus, Z2_lminus_array)
        Z2_lplus_array = fill_kinematic_array(Z2_lplus, Z2_lplus_array)
        wgt_array[:] = event.weights[1:10]
        tree.Fill()
        i = i+1

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
tree.Print()
file.Write()
