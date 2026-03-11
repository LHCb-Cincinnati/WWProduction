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

def fill_kinematic_array(particle, array):
    mom = particle.momentum
    array[:] = (mom.px, mom.py, mom.pz, mom.pt(), 
                np.sqrt(mom.px *mom.px + mom.py * mom.py + mom.pz * mom.pz), 
                mom.eta(), mom.e, mom.phi(), mom.m(), particle.pid, 
                particle.status)
    return(array)

def fill_scalewgt_array(event, weight_name_list):
    array = np.array([0]*len(weight_name_list), dtype=np.float32)
    for index, weight_name in enumerate(weight_name_list):
        weight_index = event.weight_names.index(weight_name)
        array[index] = event.weights[weight_index] 
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
ifile_name = "/data/home/ganowak/MG5_aMC_v2_9_25/tW_LO/Events/run_05_decayed_1/tag_1_pythia8_events.hepmc" 
# ofile_name = "tW_MG5_LO_CMTS09_mu10.root"
ofile_name = "Test.root"
# Setting up Path
ww_path = at.find_WW_path()

# Inputs
W_pid = -24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
weight_name_list = ([
    "MUF=0.5_MUR=0.5_PDF=10770_MERGING=0.000",
    "MUF=0.5_MUR=1.0_PDF=10770_MERGING=0.000",
    "MUF=0.5_MUR=2.0_PDF=10770_MERGING=0.000",
    "MUF=1.0_MUR=0.5_PDF=10770_MERGING=0.000",
    "Weight",
    "MUF=1.0_MUR=2.0_PDF=10770_MERGING=0.000",
    "MUF=2.0_MUR=0.5_PDF=10770_MERGING=0.000",
    "MUF=2.0_MUR=1.0_PDF=10770_MERGING=0.000",
    "MUF=2.0_MUR=2.0_PDF=10770_MERGING=0.000",
])

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:status/F'
wgt_str = "AUX_MUR05_MUF05/F:AUX_MUR05_MUF10/F:AUX_MUR05_MUF20/F:AUX_MUR10_MUF05/F:AUX_MUR10_MUF10/F:AUX_MUR10_MUF20/F:AUX_MUR20_MUF05/F:AUX_MUR20_MUF10/F:AUX_MUR20_MUF20/F"
initial_conditions_str = 'x1/F:id1/F:x2/F:id2/F'
pdfrwgt_str = 'nnpdf31lo/F:ct18lo/F:msht20lo/F'
evt_array = np.array([0], dtype=np.float32)
wgt_array = np.array([0]*9, dtype=np.float32)
target_particle_array = np.array([0]*11, dtype=np.float32)
target_lepton_array = np.array([0]*11, dtype=np.float32)
target_neutrino_array = np.array([0]*11, dtype=np.float32)
target_antiparticle_array = np.array([0]*11, dtype=np.float32)
target_antilepton_array = np.array([0]*11, dtype=np.float32)
target_antineutrino_array = np.array([0]*11, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
pdfrwgt_array = np.array([0]*3, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/madgraph/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('Weights', wgt_array, wgt_str)
tree.Branch(
    'InitialConditions', initial_conditions_array, initial_conditions_str
)
tree.Branch('pdfReweight', pdfrwgt_array, pdfrwgt_str)
tree.Branch('TargetParticle', target_particle_array, var_str)
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetNeutrino', target_neutrino_array, var_str)
tree.Branch('TargetAntiParticle', target_antiparticle_array, var_str)
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
tree.Branch('TargetAntiNeutrino', target_antineutrino_array, var_str)

# Variables
i = 0
counter1 = 0
counter2 = 0

with pyhepmc.open(ifile_name) as f:
    for event in f:
        evt_array[:] = i
        # Find relevant particles
        initial_parton_list = []
        for particle in event.particles:
            if particle.pid==W_pid:
                wminus = particle
            elif particle.pid==-1*W_pid:
                wplus = particle
            elif particle.status == 21:
                initial_parton_list.append(particle)
        # Find top b-daughter and w-daughter
        if (len(wplus.children) == 0) or (len(wminus.children) == 0):
            continue
        wplus = get_nonradiative_decay(wplus)
        wminus = get_nonradiative_decay(wminus)
        for daughter in wminus.children:
            if daughter.pid in lepton_pid_array:
                lepton_plus = daughter
            elif daughter.pid in -1*neutrino_pid_array:
                neutrino_plus = daughter
        for daughter in wplus.children:
            if daughter.pid in -1*lepton_pid_array:
                lepton_minus = daughter
            elif daughter.pid in neutrino_pid_array:
                neutrino_minus = daughter

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
        target_particle_array = fill_kinematic_array(wplus, target_particle_array)
        target_lepton_array = fill_kinematic_array(lepton_minus, target_lepton_array)
        target_neutrino_array = fill_kinematic_array(neutrino_minus, target_neutrino_array)
        target_antiparticle_array = fill_kinematic_array(wminus, target_antiparticle_array)
        target_antilepton_array = fill_kinematic_array(lepton_plus, target_antilepton_array)
        target_antineutrino_array = fill_kinematic_array(neutrino_plus, target_antineutrino_array)
        wgt_array[:] = fill_scalewgt_array(event, weight_name_list)
        tree.Fill()
        i = i+1

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
tree.Print()
file.Write()
