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

# NNPDF31LO Reweighter
num_NNPDF31LO_pdfsets = 101
pdfrwgt_list_NNPDF31LO = [0] * num_NNPDF31LO_pdfsets
nnpdf31lo_bra_str = ""
for i in range(num_NNPDF31LO_pdfsets):
    pdfrwgt_list_NNPDF31LO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="NNPDF31_lo_as_0118",
        pdf_set_to_num = i
    )
    nnpdf31lo_bra_str += f"PDFMember{i}Weight/F:"
nnpdf31lo_bra_str = nnpdf31lo_bra_str[:-1]

# CT18LO Reweighter
num_CT18LO_pdfsets = 1
pdfrwgt_list_CT18LO = [0] * num_CT18LO_pdfsets
ct18lo_bra_str = ""
for i in range(num_CT18LO_pdfsets):
    pdfrwgt_list_CT18LO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="CT18LO",
        pdf_set_to_num = i
    )
    ct18lo_bra_str += f"PDFMember{i}Weight/F:"
ct18lo_bra_str = ct18lo_bra_str[:-1]

# MSHT20LO Reweighter
num_MSHT20LO_pdfsets = 61
pdfrwgt_list_MSHT20LO = [0] * num_MSHT20LO_pdfsets
msht20lo_bra_str = ""
for i in range(num_MSHT20LO_pdfsets):
    pdfrwgt_list_MSHT20LO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="MSHT20lo_as130",
        pdf_set_to_num = i
    )
    msht20lo_bra_str += f"PDFMember{i}Weight/F:"
msht20lo_bra_str = msht20lo_bra_str[:-1]

# Files
ifile_name = "/data/home/ganowak/MG5_aMC_v2_9_25/WZ_LO/Events/run_16_decayed_1/tag_1_pythia8_events.hepmc" 
ofile_name = "WZ_MG5_LO_CMTS09_mu10.root"
# ofile_name = "Test.root"
# Setting up Path
ww_path = at.find_WW_path()

# Inputs
Z_id = 23
W_id = 24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
jet_array = np.array([-5,-3,-1])
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
Z_array = np.array([0]*11, dtype=np.float32)
W_array = np.array([0]*11, dtype=np.float32)
Z_lminus_array = np.array([0]*11, dtype=np.float32)
Z_lplus_array = np.array([0]*11, dtype=np.float32)
W_lepton_array = np.array([0]*11, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
pdfrwgt_array = np.array([0]*3, dtype=np.float32)
nnpdf31_wgt_array = np.array([0]*num_NNPDF31LO_pdfsets, dtype=np.float32) 
ct18lo_wgt_array = np.array([0]*num_CT18LO_pdfsets, dtype=np.float32) 
msht20lo_wgt_array = np.array([0]*num_MSHT20LO_pdfsets, dtype=np.float32) 

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/madgraph/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch(
    'InitialConditions', initial_conditions_array, initial_conditions_str
)
tree.Branch('pdfReweight', pdfrwgt_array, pdfrwgt_str)
tree.Branch('NNPDF31LO_Members', nnpdf31_wgt_array, nnpdf31lo_bra_str)
tree.Branch('CT18LO_Members', ct18lo_wgt_array, ct18lo_bra_str)
tree.Branch('MSHT20LO_Members', msht20lo_wgt_array, msht20lo_bra_str)
tree.Branch('Weights', wgt_array, wgt_str)
tree.Branch('Z', Z_array, var_str)
tree.Branch('W', W_array, var_str)
tree.Branch('Z_lminus', Z_lminus_array, var_str)
tree.Branch('Z_lplus', Z_lplus_array, var_str)
tree.Branch('W_lepton', W_lepton_array, var_str)

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
            if particle.pid==Z_id:
                Z = particle
            elif abs(particle.pid)==W_id:
                W = particle
            elif particle.status == 21:
                initial_parton_list.append(particle)
        Z = get_nonradiative_decay(Z)
        W = get_nonradiative_decay(W)
        # Find leptonic daughters
        for daughter in Z.children:
            if daughter.pid in lepton_pid_array:
                Z_lminus = daughter
            elif daughter.pid in -1*lepton_pid_array:
                Z_lplus = daughter
        for daughter in W.children:
            if daughter.pid in lepton_pid_array:
                W_lepton = daughter
            elif daughter.pid in -1*lepton_pid_array:
                W_lepton = daughter

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
        # Fill large PDF member arrays
        for i in range(num_NNPDF31LO_pdfsets):
            nnpdf31_wgt_array[i] = pdfrwgt_list_NNPDF31LO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )
        for i in range(num_CT18LO_pdfsets):
            ct18lo_wgt_array[i] = pdfrwgt_list_CT18LO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )
        for i in range(num_MSHT20LO_pdfsets):
            msht20lo_wgt_array[i] = pdfrwgt_list_MSHT20LO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )

        # Fill arrays
        pdfrwgt_array[:] = (
            nnpdf31_wgt_array[0], ct18lo_wgt_array[0], msht20lo_wgt_array[0]
        )
        initial_conditions_array[:] = (x1, id1, x2, id2)
        Z_array = fill_kinematic_array(Z, Z_array)
        W_array = fill_kinematic_array(W, W_array)
        Z_lminus_array = fill_kinematic_array(Z_lminus, Z_lminus_array)
        Z_lplus_array = fill_kinematic_array(Z_lplus, Z_lplus_array)
        W_lepton_array = fill_kinematic_array(W_lepton, W_lepton_array)
        wgt_array[:] = fill_scalewgt_array(event, weight_name_list)
        tree.Fill()
        i = i+1

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
tree.Print()
file.Write()
