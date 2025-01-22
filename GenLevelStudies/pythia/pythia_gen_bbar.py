# Imports
# STL Packages
import sys
import pdb
import math
# Scikit Packages
import numpy as np
# HEP Packages
import pythia8
import ROOT
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Functions
def isBHadron(part):
    if math.floor(part.idAbs()/100) == 5:
        return(True)
    else:
        return(False)

def get_daughter(part):
    for index in range(part.daughter2() - part.daughter1()):
        daughter = pythia.event[index + part.daughter1()]
        print(daughter.id())
        if isBHadron(daughter):
            get_daughter(daughter)
        elif (daughter.idAbs() == 11) | (daughter.idAbs() == 13):
            return(index + part.daughter1())
        else:
            return(False)

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/data/home/ganowak/WWProduction/pythia8311/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "bbarStrongProduction.cmnd"
ofile_name = "bbarStrongProduction.root"

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
targetlepton_array = np.array([0]*12, dtype=np.float32)
targetantilepton_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('TargetLepton', targetlepton_array, var_str)
tree.Branch('TargetAntiLepton', targetantilepton_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    lminus_index = 0
    lminus_pT = 0
    lplus_index = 0
    lplus_pT = 0
    if not pythia.next():
        continue
    
    for index, part in enumerate(pythia.event):
        part_pt = part.pT()
        if (part.id() == 11) | (part.id() == 13):
            if part_pt > lminus_pT:
                lminus_index = index
                lminus_pT = part_pt
        if (part_pt == -11) | (part.id() == -13):
            if part_pt > lplus_pT:
                lplus_index = index
                lplus_pT = part_pt
    
    targetlepton_array = at.fill_array(
        targetlepton_array,
        pythia.event,
        lminus_index
    )
    targetantilepton_array = at.fill_array(
        targetantilepton_array,
        pythia.event,
        lplus_index
    )
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
file.Write()
pythia.stat()
tree.Print()
