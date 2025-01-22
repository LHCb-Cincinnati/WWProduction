# Imports
# STL Packages
import sys
import pdb
# Scikit Packages
import numpy as np
# HEP Packages
import pythia8
import ROOT
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

def deltaR(vec1, vec2):
    delta_r = np.sqrt(
        (vec1.phi() - vec2.phi())**2
        + (vec1.eta() - vec2.eta())**2
    )
    return(delta_r)

def fill_UE_array(array, event, lepton_index):
    lepton_vec = event[lepton_index].p()
    spT_sum = 0
    num_good_cones = 0
    isGoodCone_array = [0]*6
    maxpT_array = [0]*6
    spT_array = [0]*6
    for i in range(6):
        cone_spT = 0
        cone_max_pT = 0
        isGoodCone_array[i] = True
        lepton_vec.rot(0, (2*np.pi) / 6)
        for particle in event:
            # Find nearby particles
            if ((deltaR(particle.p(), lepton_vec) < 0.50) and (particle.isFinal())):
                # Find max pT object
                if particle.pT() > cone_max_pT:
                    cone_max_pT = particle.pT()
                cone_spT += particle.pT()
                # Find if particles spoil the cone.
                if (particle.pT() > 20):
                    isGoodCone_array[i] = False
        spT_array[i] = cone_spT
        maxpT_array[i] = cone_max_pT
    good_cone_indices_list = np.where(isGoodCone_array)[0]
    for i in range(5):
        if i in good_cone_indices_list:
            spT_sum += spT_array[i]
            num_good_cones += 1
    if num_good_cones == 0:
        spTAve = 0
    else:
        spTAve = spT_sum / num_good_cones
    array[:] = (
        spTAve,
        num_good_cones,
        max(maxpT_array),
        spT_array[0],
        isGoodCone_array[0],
        maxpT_array[0],
        spT_array[1],
        isGoodCone_array[1],
        maxpT_array[1],
        spT_array[2],
        isGoodCone_array[2],
        maxpT_array[2],
        spT_array[3],
        isGoodCone_array[3],
        maxpT_array[3],
        spT_array[4],
        isGoodCone_array[4],
        maxpT_array[4],
        spT_array[5],
        isGoodCone_array[5],
        maxpT_array[5],
    )
    return(array)

def fill_u_array(array, event, lepton_index):
    iw = event[lepton_index].mother1()
    u = (
        (
            (event[lepton_index].px() * event[iw].px())
            + (event[lepton_index].py() * event[iw].py())
        )
        / event[lepton_index].pT()
    )
    array[:] = (u)
    return(array)

def fill_iso_array(array, event, lepton_index):
    cone_spT = 0
    lepton_vec = event[lepton_index].p()
    for particle in event:
        # Find nearby particles
        if ((deltaR(particle.p(), lepton_vec) < 0.50) and (deltaR(particle.p(), lepton_vec) > 0.05) and (particle.isFinal())):
            cone_spT += particle.pT()
    array[:] = (cone_spT)
    return(array)

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/data/home/ganowak/WWProduction/pythia8311/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "DrellYanProduction.cmnd"
ofile_name = "DrellYanProductionUE.root"
target_pid = -24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
# jet_array = np.array([-5,-3,-1])
jet_array = np.array([])

# Read in Card File 
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:charge/F:status/F'
ue_str = "spTAve/F:NumGoodCones/F:MaxpT/F:Cone1spT/F:Cone1isGoodCone/F:Cone1MaxPt/F:Cone2spT/F:Cone2isGoodCone/F:Cone2MaxPt/F:Cone3spT/F:Cone3isGoodCone/F:Cone3MaxPt/F:Cone4spT/F:Cone4isGoodCone/F:Cone4MaxPt/F:Cone5spT/F:Cone5isGoodCone/F:Cone5MaxPt/F:Cone6spT/F:Cone6isGoodCone/F:Cone6MaxPt/F:Cone7spT/F:Cone7isGoodCone/F:Cone7MaxPt/F"
evt_array = np.array([0], dtype=np.float32)
target_lepton_array = np.array([0]*12, dtype=np.float32)
target_lepton_UE_array = np.array([0]*21, dtype=np.float32)
target_lepton_u = np.array([0], dtype=np.float32)
target_lepton_isolation = np.array([0], dtype=np.float32)
target_antilepton_array = np.array([0]*12, dtype=np.float32)
target_antilepton_UE_array = np.array([0]*21, dtype=np.float32)
target_antilepton_u = np.array([0], dtype=np.float32)
target_antilepton_isolation = np.array([0], dtype=np.float32)
if jet_array.any():
    target_jet_array = np.array([0]*12, dtype=np.float32)
    target_jet_UE_array = np.array([0]*21, dtype=np.float32)
    target_antijet_array = np.array([0]*12, dtype=np.float32)
    target_antijet_UE_array = np.array([0]*21, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetLeptonUE', target_lepton_UE_array, ue_str)
tree.Branch("TargetLeptonu", target_lepton_u, "u/F")
tree.Branch("TargetLeptonIso", target_lepton_isolation, "spT/F")
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
tree.Branch('TargetAntiLeptonUE', target_antilepton_UE_array, ue_str)
tree.Branch("TargetAntiLeptonu", target_antilepton_u, "u/F")
tree.Branch("TargetAntiLeptonIso", target_antilepton_isolation, "spT/F")
if jet_array.any():
    tree.Branch('TargetJet', target_jet_array, var_str)
    tree.Branch('TargetJetUE', target_jet_UE_array, ue_str)
    tree.Branch('TargetAntiJet', target_antijet_array, var_str)
    tree.Branch('TargetAntiJetUE', target_antijet_UE_array, ue_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    ilepton_minus = 0
    ilepton_plus = 0
    evt_array[:] = iEvent
    # Loop through particles to find leptons and jets 
    for index, particle in enumerate(pythia.event):
        # if particle.statusAbs() == 23:
        if particle.idAbs() == 15:
            part_children = at.get_target_children(index, pythia.event, child_index_list=[])
            for ichild in part_children:
                if pythia.event[ichild].id() in lepton_pid_array:
                    ilepton_minus = ichild
                elif pythia.event[ichild].id() in -1*lepton_pid_array:
                    ilepton_plus = ichild
                elif particle.id() in jet_array:
                    ijet = index
                elif particle.id() in -1*jet_array:
                    iantijet = index

    # Fill Arrays
    target_lepton_array = at.fill_array(
        target_lepton_array, pythia.event, ilepton_minus
    )
    target_lepton_UE_array = fill_UE_array(
        target_lepton_UE_array, pythia.event, ilepton_minus
    )
    target_lepton_u = fill_u_array(target_lepton_u, pythia.event, ilepton_minus)
    target_lepton_isolation = fill_iso_array(target_lepton_isolation, pythia.event, ilepton_minus)
    target_antilepton_array = at.fill_array(
        target_antilepton_array, pythia.event, ilepton_plus
    )
    target_antilepton_UE_array = fill_UE_array(
        target_antilepton_UE_array, pythia.event, ilepton_plus
    )
    target_antilepton_u = fill_u_array(target_antilepton_u, pythia.event, ilepton_plus)
    target_antilepton_isolation = fill_iso_array(target_antilepton_isolation, pythia.event, ilepton_plus)
    if jet_array.any():
        target_jet_array = at.fill_array(target_jet_array, pythia.event,
                                        ijet)
        target_jet_UE_array = fill_UE_array(target_jet_UE_array, pythia.event,
                                        ijet)
        target_antijet_array = at.fill_array(target_antijet_array, pythia.event,
                                        iantijet)
        target_antijet_UE_array = fill_UE_array(target_antijet_UE_array, pythia.event,
                                        iantijet)
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
