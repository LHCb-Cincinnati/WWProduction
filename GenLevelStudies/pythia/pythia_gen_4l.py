# Imports
# STL Packages
import sys
import pdb
# Sikit Packages
import numpy as np
# HEP Packages
import pythia8
import ROOT
# Personal Packages
sys.path.append(".")
import AnalysisTools as at

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/Applications/pythia8310/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "WZProduction.cmnd"
ofile_name = "Test.root"

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(ww_path+ "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
l1_array = np.array([0]*12, dtype=np.float32)
l2_array = np.array([0]*12, dtype=np.float32)
l3_array = np.array([0]*12, dtype=np.float32)
l4_array = np.array([0]*12, dtype=np.float32)
lepton_array_list = [
    l1_array,
    l2_array,
    l3_array,
    l4_array
] 
lepton_list = [0]*4

# Set up ROOT
file = ROOT.TFile.Open(ww_path+ "/GenLevelStudies/pythia/" +ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('l1', l1_array, var_str)
tree.Branch('l2', l2_array, var_str)
tree.Branch('l3', l3_array, var_str)
tree.Branch('l4', l4_array, var_str)

# Variables
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    counter1+=1
    evt_array[:] = iEvent
    target_indices_list = []
    # Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.statusAbs() == 23:
            target_indices_list.append(index)
    # Find Target Particles Loop
    for index in range(len(target_indices_list)):
        lepton_array_list[index] =at.fill_array(
            lepton_array_list[index],
            pythia.event, 
            target_indices_list[index]
        )
        lepton_list[index] = target_indices_list[index]
    # Find good events
    good_event = (
        (lplus.pT() > 20)
        & (lminus.pT() > 20)
        & (lepton_three.pT() > 20)
        & (abs(lminus.eta()) < 2.5)
        & (abs(lplus.eta()) < 2.5)
        & (abs(lepton_three.eta()) < 2.5)
    )
    if good_event:
        counter2+=1
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()