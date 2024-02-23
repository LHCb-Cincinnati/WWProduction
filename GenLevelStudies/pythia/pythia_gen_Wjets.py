# Imports
import sys
import pdb
import yaml

import numpy as np

import pythia8
import ROOT
import AnalysisTools as at

# Functions
def fill_array(array, event, index):
    particle = event[index]
    px = particle.px()
    py = particle.py()
    pz = particle.pz()
    p = np.sqrt(px*px + py*py + pz*pz)
    array[:] = (px, py, pz,
                particle.pT(), p, particle.eta(), particle.e(),
                particle.phi(), particle.m0(),particle.id(),
                particle.charge(), particle.status())
    return(array)

def GetTargetChildren(target_index, event, child_index_list):
    target = event[target_index]
    for idaughter in target.daughterList():
        child_index_list.append(idaughter)
        if abs(event[idaughter].id()) == 24:
            child_index_list = GetTargetChildren(idaughter, event,
                                                child_index_list)
    return(child_index_list)

def check_mother(event, particle, pid):
    imother = particle.mother1()
    mother = event[imother]
    if mother.id() == pid:
        return(True)
    elif mother.id() == particle.id():
        return(check_mother(event, mother, pid))
    else:
        return(False)

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/Applications/pythia8310/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "WjetsProduction.cmnd"
ofile_name = "WjetsProduction.root"

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
target_jet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetJet', target_jet_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent
    target_indices_list = []
    # Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.statusAbs() == 23: # Look for outgoing particles of hard process
            target_indices_list.append(index)
    # Find Target Particles Loop
    for itarget in target_indices_list:
        target_particle = pythia.event[itarget]
        itarget_mother = target_particle.mother1()
        if (pythia.event[itarget_mother].idAbs() == 24):
            if target_particle.idAbs() in [11, 13]:
                target_lepton_array = fill_array(target_lepton_array, pythia.event, itarget)
        else:
            target_jet_array = fill_array(target_jet_array, pythia.event, itarget)
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
