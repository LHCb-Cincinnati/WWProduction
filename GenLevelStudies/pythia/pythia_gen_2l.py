# Imports
# STL Packages
import sys
import pdb
import os
# Scikit Packages
import numpy as np
# HEP Packages
import pythia8
import ROOT
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Functions
def deltaR(vec1, vec2):
    delta_r = np.sqrt(
        (vec1.phi() - vec2.phi())**2
        + (vec1.eta() - vec2.eta())**2
    )
    return(delta_r)

def fill_UE_array(array, pythia, lepton_index):
    lepton_vec = pythia.event[lepton_index].p()
    spT_sum = 0
    num_good_cones = 0
    isGoodCone_array = [0]*7
    spT_array = [0]*7
    for i in range(7):
        cone_spT = 0
        is_good_cone = True
        lepton_vec.rot(0, (2*np.pi) / 6)
        for particle in pythia.event:
            if ((deltaR(particle.p(), lepton_vec) < 0.50) and (particle.pT() < 20)):
                cone_spT += particle.pT()
            elif (particle.pT() > 20):
                is_good_cone = False
                isGoodCone_array[i] = False
                spT_array[i] = 0
                break
        isGoodCone_array[i] = True
        spT_array[i] = cone_spT
        num_good_cones +=1
        spT_sum += cone_spT
    array[:] = (
        spT_sum,
        num_good_cones,
        spT_array[0],
        isGoodCone_array[0],
        spT_array[1],
        isGoodCone_array[1],
        spT_array[2],
        isGoodCone_array[2],
        spT_array[3],
        isGoodCone_array[3],
        spT_array[4],
        isGoodCone_array[4],
        spT_array[5],
        isGoodCone_array[5]
    )
    return(array)

# Setting up Path
ww_path = at.find_WW_path()

# Inputs
card_file_name = "ttbarProduction.cmnd"
ofile_name = "Test.root"
target_pid = -24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
# jet_array = np.array([-5,-3,-1])
jet_array = np.array([])

# Read in Card File 
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/defaults.cmnd")
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)
if "PYTHIA_NUMBER_EVENTS" in os.environ:
    pythia.readString("Main:numberOfEvents = " + os.environ["PYTHIA_NUMBER_EVENTS"])

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:charge/F:status/F'
ue_str = "spTAve/F:NumGoodCones/F:Cone1spT/F:Cone1isGoodCone/F:Cone2spT/F:Cone2isGoodCone/F:Cone3spT/F:Cone3isGoodCone/F:Cone4spT/F:Cone4isGoodCone/F:Cone5spT/F:Cone5isGoodCone/F:Cone6spT/F:Cone6isGoodCone/F:"
initial_conditions_str = 'x1/F:id1/F:x2/F:id2/F'
evt_array = np.array([0], dtype=np.float32)
target_particle_array = np.array([0]*12, dtype=np.float32)
born_target_lepton_array = np.array([0]*12, dtype=np.float32)
bare_target_lepton_array = np.array([0]*12, dtype=np.float32)
target_lepton_UE_array = np.array([0]*14, dtype=np.float32)
target_neutrino_array = np.array([0]*12, dtype=np.float32)
target_antiparticle_array = np.array([0]*12, dtype=np.float32)
born_target_antilepton_array = np.array([0]*12, dtype=np.float32)
bare_target_antilepton_array = np.array([0]*12, dtype=np.float32)
target_antilepton_UE_array = np.array([0]*14, dtype=np.float32)
target_antineutrino_array = np.array([0]*12, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
if jet_array.any():
    target_jet_array = np.array([0]*12, dtype=np.float32)
    target_antijet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch(
    'InitialConditions', initial_conditions_array, initial_conditions_str
)
tree.Branch('TargetParticle', target_particle_array, var_str)
tree.Branch('BornTargetLepton', born_target_lepton_array, var_str)
tree.Branch('BareTargetLepton', bare_target_lepton_array, var_str)
tree.Branch('TargetLeptonUE', target_lepton_UE_array, ue_str)
tree.Branch('TargetNeutrino', target_neutrino_array, var_str)
tree.Branch('TargetAntiParticle', target_antiparticle_array, var_str)
tree.Branch('BornTargetAntiLepton', born_target_antilepton_array, var_str)
tree.Branch('BareTargetAntiLepton', bare_target_antilepton_array, var_str)
tree.Branch('TargetAntiLeptonUE', target_antilepton_UE_array, ue_str)
tree.Branch('TargetAntiNeutrino', target_antineutrino_array, var_str)
if jet_array.any():
    tree.Branch('TargetJet', target_jet_array, var_str)
    tree.Branch('TargetAntiJet', target_antijet_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent
    itarget_particle = None
    itarget_antiparticle = None
    ilepton_plus = None
    ilepton_minus = None
    ibare_lepton_plus = None
    ibare_lepton_minus = None
    ineutrino_plus = None
    ineutrino_minus = None
    ijet = None
    iantijet = None
    for index, particle in enumerate(pythia.event):
        if particle.id()==target_pid:
            itarget_particle = index
        elif particle.id()==-1*target_pid:
            itarget_antiparticle = index
    if itarget_particle is None or itarget_antiparticle is None:
        print(f"WARNING: could not find both target particles in event {iEvent}; event skipped.")
        continue
    for idaughter in at.get_target_children(itarget_particle, pythia.event, child_index_list=[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in lepton_pid_array:
            ilepton_plus = idaughter
        elif daughter.id() in -1*neutrino_pid_array:
            ineutrino_plus = idaughter
        elif daughter.id() in jet_array:
            ijet = idaughter
    for idaughter in at.get_target_children(itarget_antiparticle, pythia.event, child_index_list=[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in -1*lepton_pid_array:
            ilepton_minus = idaughter
        elif daughter.id() in neutrino_pid_array:
            ineutrino_minus = idaughter
        elif daughter.id() in -1*jet_array:
            iantijet = idaughter

    required_indices = [
        ilepton_plus, ilepton_minus, ineutrino_plus, ineutrino_minus
    ]
    if any(index is None for index in required_indices):
        print(f"WARNING: could not find all target leptons/neutrinos in event {iEvent}; event skipped.")
        continue

    ibare_lepton_plus = at.get_nonradiative_decay(pythia.event, ilepton_plus)
    ibare_lepton_minus = at.get_nonradiative_decay(pythia.event, ilepton_minus)
    if ibare_lepton_plus is None or ibare_lepton_minus is None:
        print(f"WARNING: ambiguous bare lepton history in event {iEvent}; event skipped.")
        continue

    # Get Initial Conditions Info
    x1pdf = pythia.infoPython().x1pdf()
    x2pdf = pythia.infoPython().x2pdf()
    id1pdf = pythia.infoPython().id1pdf()
    id2pdf = pythia.infoPython().id2pdf()
    initial_conditions_array[:] = (x1pdf, id1pdf, x2pdf, id2pdf)
    target_particle_array = at.fill_array(target_particle_array, pythia.event,
                                        itarget_particle)
    born_target_lepton_array = at.fill_array(born_target_lepton_array,
                                        pythia.event, ilepton_plus)
    bare_target_lepton_array = at.fill_array(bare_target_lepton_array,
                                        pythia.event, ibare_lepton_plus)
    target_lepton_UE_array = fill_UE_array(target_lepton_UE_array, pythia, ibare_lepton_plus)
    target_neutrino_array = at.fill_array(target_neutrino_array, pythia.event,
                                        ineutrino_plus)
    target_antiparticle_array = at.fill_array(target_antiparticle_array, pythia.event,
                                        itarget_antiparticle)
    born_target_antilepton_array = at.fill_array(born_target_antilepton_array,
                                        pythia.event, ilepton_minus)
    bare_target_antilepton_array = at.fill_array(bare_target_antilepton_array,
                                        pythia.event, ibare_lepton_minus)
    target_antilepton_UE_array = fill_UE_array(target_antilepton_UE_array, pythia, ibare_lepton_minus)
    target_antineutrino_array = at.fill_array(target_antineutrino_array, pythia.event,
                                        ineutrino_minus)
    if jet_array.any():
        target_jet_array = at.fill_array(target_jet_array, pythia.event,
                                        ijet)
        target_antijet_array = at.fill_array(target_antijet_array, pythia.event,
                                        iantijet)
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
