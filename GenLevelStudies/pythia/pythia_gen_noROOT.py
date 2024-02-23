# Imports
import sys
import pdb

import numpy as np

import pythia8
import awkward as ak
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
card_file_name = "ttbarProduction.cmnd"
ofile_name = "Test.root"
target_pid = -6
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
jet_array = np.array([-5,-3,-1])

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

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
target_jet_array = np.array([0]*12, dtype=np.float32)
target_antijet_array = np.array([0]*12, dtype=np.float32)
array_builder = ak.ArrayBuilder()

# Testing
counter1 = 0
counter2 = 0

# Event Loop

for iEvent in range(nEvent):
    with array_builder.record():    
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
        if jet_array.any():
            target_jet_array = fill_array(target_jet_array, pythia.event,
                                            ijet)
            target_antijet_array = fill_array(target_antijet_array, pythia.event,
                                            iantijet)
        array_builder.field("Particle").append( target_particle_array)
        array_builder.field("Lepton").append(target_lepton_array)
        array_builder.field("Neutrino").append(target_neutrino_array)
        array_builder.field("Jet").append(target_jet_array)
        array_builder.field("AntiParticle").append(target_antiparticle_array)
        array_builder.field("AntiLepton").append(target_antilepton_array)
        array_builder.field("AntiNeutrino").append(target_antineutrino_array)
        array_builder.field("Antijet").append(target_antijet_array)
        # if check_mother(pythia.event, pythia.event[iW_plus], 6):
        #     counter1+=1
        # if check_mother(pythia.event, pythia.event[iW_minus], -6):
        #     counter2+=1
data_array = array_builder.snapshot()
print(data_array.Lepton)
print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
