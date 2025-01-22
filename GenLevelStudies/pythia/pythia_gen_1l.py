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

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/data/home/ganowak/WWProduction/pythia8311/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "ZJets.cmnd"
ofile_name = "Test.root"
target_pid = 23

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
target_lepton_array = np.array([0]*12, dtype=np.float32)
target_antilepton_array = np.array([0]*12, dtype=np.float32)
target_jet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
tree.Branch('TargetJet', target_jet_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    counter1+=1
    evt_array[:] = iEvent
    target_indices_list = []
    for index, particle in enumerate(pythia.event):
        if particle.statusAbs() == 23: # Look for outgoing particles of hard process
                    target_indices_list.append(index)
    
    for iparticle in target_indices_list:
        if pythia.event[iparticle].id() == 13 or pythia.event[iparticle].id() == 11:
            ilepton_minus = iparticle
            lepton_minus = pythia.event[iparticle]
        elif pythia.event[iparticle].id() == -13 or pythia.event[iparticle].id() == -11:
            ilepton_plus = iparticle
            lepton_plus = pythia.event[iparticle]
        else:
            itarget_jet = iparticle
            jet = pythia.event[iparticle]

    good_event = (
        (lepton_plus.pT() > 20)
        & (lepton_minus.pT() > 20)
        & ((jet.pT() > 20))
        & (abs(lepton_minus.eta()) < 2.5)
        & (abs(lepton_plus.eta()) < 2.5)
        & (abs(jet.eta()) < 5)
    )
    if good_event:
        counter2+=1
    target_lepton_array = at.fill_array(target_lepton_array, pythia.event,
                                        ilepton_plus)
    target_antilepton_array = at.fill_array(target_antilepton_array, pythia.event,
                                        ilepton_minus)
    target_jet_array = at.fill_array(target_jet_array, pythia.event,
                                    itarget_jet)

    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
