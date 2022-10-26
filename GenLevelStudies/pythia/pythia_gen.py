# Imports
import sys
import pdb
import yaml

import numpy as np

import pythia8
import ROOT

# Functions
def fill_array(array, event, index):
    particle = event[index]
    px = particle.px()
    py = particle.py()
    pz = particle.pz()
    p = np.sqrt(px*px + py*py + pz*pz)
    array[:] = (particle.pT(), p, particle.eta(), particle.e(),
                particle.phi(), particle.m0(),particle.id(),
                particle.charge(), particle.status())
    return(array)

#  Set up makefile configuration.
cfg = open("Makefile.inc")
lib = "/Applications/pythia8307/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Parse Inputs
with open('config.yaml', 'r') as file:
    config_dict = yaml.safe_load(file)
card_file_name = config_dict['card_file']
ofile_name = config_dict['output_file']
additional_particles_dict = config_dict['additional_particles']

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
wplus_array = np.array([0]*9, dtype=np.float32)
lplus_array = np.array([0]*9, dtype=np.float32)
nplus_array = np.array([0]*9, dtype=np.float32)
wminus_array = np.array([0]*9, dtype=np.float32)
lminus_array = np.array([0]*9, dtype=np.float32)
nminus_array = np.array([0]*9, dtype=np.float32)
additional_particles_array = [np.array([0]*9, dtype=np.float32)
                            *len(additional_particles_dict)]


# Set up ROOT
file = ROOT.TFile.Open(ofile_name,"RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('WeakBosonPlus', wplus_array, var_str)
tree.Branch('LeptonPlus', lplus_array, var_str)
tree.Branch('NeutrinoPlus', nplus_array, var_str)
tree.Branch('WeakBosonMinus', wminus_array, var_str)
tree.Branch('LeptonMinus', lminus_array, var_str)
tree.Branch('NeutrinoMinus', nminus_array, var_str)
for array, key in zip(additional_particles_array, additional_particles_dict.keys()):
    tree.Branch(additional_particles_dict[key], nminus_array, var_str)

# Testing

for iEvent in range(10):
    if not pythia.next():
        continue
    evt_array[:] = iEvent
    for index, particle in enumerate(pythia.event):
        if particle.id()==5:
            print('Bottom')
        if particle.id()==6:
            print('Top')
        if particle.id()==24:
            iW_plus = index
            ilepton_plus = particle.daughter1()
            ineutrino_plus = particle.daughter2()
        elif particle.id()==-24:
            iW_minus = index
            ilepton_minus = particle.daughter1()
            ineutrino_minus = particle.daughter2()
    wplus_array = fill_array(wplus_array, pythia.event, iW_plus)
    lplus_array = fill_array(lplus_array, pythia.event, ilepton_plus)
    nplus_array = fill_array(nplus_array, pythia.event, ineutrino_plus)
    wminus_array = fill_array(wminus_array, pythia.event, iW_minus)
    lminus_array = fill_array(lminus_array, pythia.event, ilepton_minus)
    nminus_array = fill_array(nminus_array, pythia.event, ineutrino_minus)
    tree.Fill()

tree.Print()
file.Write()
