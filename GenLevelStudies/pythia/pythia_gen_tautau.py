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
card_file_name = "MuMuProduction.cmnd"
ofile_name = "Test.root"

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
tau_array = np.array([0]*12, dtype=np.float32)
taubar_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('lepton', tau_array, var_str)
tree.Branch('leptonbar', taubar_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    counter1+=1
    for index, particle in enumerate(pythia.process):
        if particle.id()==13:
            itau = index
        elif particle.id()==-13:
            itaubar = index
    tau_array = at.fill_array(
        tau_array,
        pythia.process,
        itau
    )
    taubar_array = at.fill_array(
        taubar_array,
        pythia.process,
        itaubar
    )
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
file.Write()
pythia.stat()
tree.Print()
