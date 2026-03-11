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

# NNPDF31NLO Reweighter
num_NNPDF31NLO_pdfsets = 101
pdfrwgt_list_NNPDF31NLO = [0] * num_NNPDF31NLO_pdfsets
nnpdf31nlo_bra_str = ""
for i in range(num_NNPDF31NLO_pdfsets):
    pdfrwgt_list_NNPDF31NLO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="NNPDF31_nlo_as_0118",
        pdf_set_to_num = i
    )
    nnpdf31nlo_bra_str += f"PDFMember{i}Weight/F:"
nnpdf31nlo_bra_str = nnpdf31nlo_bra_str[:-1]

# CT18NLO Reweighter
num_CT18NLO_pdfsets = 59
pdfrwgt_list_CT18NLO = [0] * num_CT18NLO_pdfsets
ct18nlo_bra_str = ""
for i in range(num_CT18NLO_pdfsets):
    pdfrwgt_list_CT18NLO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="CT18NLO",
        pdf_set_to_num = i
    )
    ct18nlo_bra_str += f"PDFMember{i}Weight/F:"
ct18nlo_bra_str = ct18nlo_bra_str[:-1]

# MSHT20NLO Reweighter
num_MSHT20NLO_pdfsets = 65
pdfrwgt_list_MSHT20NLO = [0] * num_MSHT20NLO_pdfsets
msht20nlo_bra_str = ""
for i in range(num_MSHT20NLO_pdfsets):
    pdfrwgt_list_MSHT20NLO[i] = PDFReweighter.PDFReweighter(
        pdf_set_from="CT09MCS",
        pdf_set_to="MSHT20nlo_as118",
        pdf_set_to_num = i
    )
    msht20nlo_bra_str += f"PDFMember{i}Weight/F:"
msht20nlo_bra_str = msht20nlo_bra_str[:-1]

# Files
ifile_name = "/data/home/ganowak/MG5_aMC_v2_9_25/Z_tautau_NLO/Events/run_13/events_PYTHIA8_0.hepmc" 
ofile_name = "DFDY_NLO_InvMassGT500.root"
# Setting up Path
ww_path = at.find_WW_path()

# Inputs
tau_pid = 15
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:status/F'
wgt_str = "AUX_MUR05_MUF05/F:AUX_MUR05_MUF10/F:AUX_MUR05_MUF20/F:AUX_MUR10_MUF05/F:AUX_MUR10_MUF10/F:AUX_MUR10_MUF20/F:AUX_MUR20_MUF05/F:AUX_MUR20_MUF10/F:AUX_MUR20_MUF20/F"
initial_conditions_str = 'x1/F:id1/F:x2/F:id2/F'
pdfrwgt_str = 'nnpdf31lo/F:ct18lo/F:msht20lo/F:nnpdf31nlo/F:ct18nlo/F:msht20nlo/F'
evt_array = np.array([0], dtype=np.float32)
wgt_array = np.array([0]*9, dtype=np.float32)
target_particle_array = np.array([0]*11, dtype=np.float32)
target_lepton_array = np.array([0]*11, dtype=np.float32)
target_neutrino_array = np.array([0]*11, dtype=np.float32)
target_antiparticle_array = np.array([0]*11, dtype=np.float32)
target_antilepton_array = np.array([0]*11, dtype=np.float32)
target_antineutrino_array = np.array([0]*11, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
pdfrwgt_array = np.array([0]*6, dtype=np.float32)
nnpdf31lo_wgt_array = np.array([0]*num_NNPDF31LO_pdfsets, dtype=np.float32) 
ct18lo_wgt_array = np.array([0]*num_CT18LO_pdfsets, dtype=np.float32) 
msht20lo_wgt_array = np.array([0]*num_MSHT20LO_pdfsets, dtype=np.float32) 
nnpdf31nlo_wgt_array = np.array([0]*num_NNPDF31NLO_pdfsets, dtype=np.float32) 
ct18nlo_wgt_array = np.array([0]*num_CT18NLO_pdfsets, dtype=np.float32) 
msht20nlo_wgt_array = np.array([0]*num_MSHT20NLO_pdfsets, dtype=np.float32) 

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
tree.Branch('NNPDF31LO_Members', nnpdf31lo_wgt_array, nnpdf31lo_bra_str)
tree.Branch('CT18LO_Members', ct18lo_wgt_array, ct18lo_bra_str)
tree.Branch('MSHT20LO_Members', msht20lo_wgt_array, msht20lo_bra_str)
tree.Branch('NNPDF31NLO_Members', nnpdf31nlo_wgt_array, nnpdf31nlo_bra_str)
tree.Branch('CT18NLO_Members', ct18nlo_wgt_array, ct18nlo_bra_str)
tree.Branch('MSHT20NLO_Members', msht20nlo_wgt_array, msht20nlo_bra_str)
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
        # Find top and anti-top in event.particles
        #for particle in event.particles:
        #    if particle.pid==tau_pid:
        #        tau = particle
        #    elif particle.pid==-1*tau_pid:
        #        antitau = particle
        tau_list = [particle for particle in event.particles if particle.pid == tau_pid]
        antitau_list = [particle for particle in event.particles if particle.pid == -1*tau_pid]
        # Find relevant particles
        initial_parton_list = []
        for particle in event.particles:
            if particle.status == 21:
                initial_parton_list.append(particle)
        antitau = get_nonradiative_decay(antitau_list[0])
        tau = get_nonradiative_decay(tau_list[0])
        for daughter in tau.children:
            if daughter.pid in lepton_pid_array:
                lepton_minus = daughter
            elif daughter.pid in -1*neutrino_pid_array:
                neutrino_minus = daughter
        for daughter in antitau.children:
            if daughter.pid in -1*lepton_pid_array:
                lepton_plus = daughter
            elif daughter.pid in neutrino_pid_array:
                neutrino_plus = daughter

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
            nnpdf31lo_wgt_array[i] = pdfrwgt_list_NNPDF31LO[i].calculate_reweighting_factor(
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
        for i in range(num_NNPDF31NLO_pdfsets):
            nnpdf31nlo_wgt_array[i] = pdfrwgt_list_NNPDF31NLO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )
        for i in range(num_CT18NLO_pdfsets):
            ct18nlo_wgt_array[i] = pdfrwgt_list_CT18NLO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )
        for i in range(num_MSHT20NLO_pdfsets):
            msht20nlo_wgt_array[i] = pdfrwgt_list_MSHT20NLO[i].calculate_reweighting_factor(
                x1, x2, id1, id2, sqrt_s=13e3
            )
        # Fill arrays
        pdfrwgt_array[:] = (
            nnpdf31lo_wgt_array[0], ct18lo_wgt_array[0], msht20lo_wgt_array[0],
            nnpdf31nlo_wgt_array[0], ct18nlo_wgt_array[0], msht20nlo_wgt_array[0]
        )
        initial_conditions_array[:] = (x1, id1, x2, id2)
        target_particle_array = fill_kinematic_array(tau_list[0], target_particle_array)
        target_lepton_array = fill_kinematic_array(lepton_minus, target_lepton_array)
        target_neutrino_array = fill_kinematic_array(neutrino_minus, target_neutrino_array)
        target_antiparticle_array = fill_kinematic_array(antitau_list[0], target_antiparticle_array)
        target_antilepton_array = fill_kinematic_array(lepton_plus, target_antilepton_array)
        target_antineutrino_array = fill_kinematic_array(neutrino_plus, target_antineutrino_array)
        wgt_array[:] = event.weights[1:10]
        tree.Fill()
        i = i+1

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
tree.Print()
file.Write()
