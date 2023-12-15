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

#  Set up makefile configuration.
cfg = open("Makefile.inc")
lib = "/Applications/pythia8307/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "WjetsProduction.cmnd"
ofile_name = "WjetsProduction.root"
target_pid = -6
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
jet_array = np.array([-5,-3,-1])

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F'
evt_array = np.array([0], dtype=np.float32)
target_particle_array = np.array([0]*12, dtype=np.float32)
target_lepton_array = np.array([0]*12, dtype=np.float32)
target_neutrino_array = np.array([0]*12, dtype=np.float32)
target_antiparticle_array = np.array([0]*12, dtype=np.float32)
target_antilepton_array = np.array([0]*12, dtype=np.float32)
target_antineutrino_array = np.array([0]*12, dtype=np.float32)
if jet_array.any():
    target_jet_array = np.array([0]*12, dtype=np.float32)
    target_antijet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ofile_name,"RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('TargetParticle', target_particle_array, var_str)
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetNeutrino', target_neutrino_array, var_str)
tree.Branch('TargetAntiParticle', target_antiparticle_array, var_str)
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
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
    for index, particle in enumerate(pythia.event):
        if particle.id()==target_pid:
            itarget_particle = index
        elif particle.id()==-1*target_pid:
            itarget_antiparticle = index
    
    for idaughter in GetTargetChildren(itarget_particle, pythia.event,[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in lepton_pid_array:
            ilepton_plus = idaughter
        elif daughter.id() in -1*neutrino_pid_array:
            ineutrino_plus = idaughter
        elif daughter.id() in jet_array:
            ijet = idaughter
    for idaughter in GetTargetChildren(itarget_antiparticle, pythia.event,[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in -1*lepton_pid_array:
            ilepton_minus = idaughter
        elif daughter.id() in neutrino_pid_array:
            ineutrino_minus = idaughter
        elif daughter.id() in -1*jet_array:
            iantijet = idaughter

    target_particle_array = fill_array(target_particle_array, pythia.event,
                                        itarget_particle)
    target_lepton_array = fill_array(target_lepton_array, pythia.event,
                                        ilepton_plus)
    target_neutrino_array = fill_array(target_neutrino_array, pythia.event,
                                        ineutrino_plus)
    target_antiparticle_array = fill_array(target_antiparticle_array, pythia.event,
                                        itarget_antiparticle)
    target_antilepton_array = fill_array(target_antilepton_array, pythia.event,
                                        ilepton_minus)
    target_antineutrino_array = fill_array(target_antineutrino_array, pythia.event,
                                        ineutrino_minus)
    # if jet_array.any():
    #     target_jet_array = fill_array(target_jet_array, pythia.event,
    #                                     ijet)
    #     target_antijet_array = fill_array(target_antijet_array, pythia.event,
    #                                     iantijet)
    # if check_mother(pythia.event, pythia.event[iW_plus], 6):
    #     counter1+=1
    # if check_mother(pythia.event, pythia.event[iW_minus], -6):
    #     counter2+=1
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
