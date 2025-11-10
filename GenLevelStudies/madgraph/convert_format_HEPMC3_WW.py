""" Program to convert HEPMC3 files to a ROOT Tree
"""

# Imports
import sys
import pdb
import pyhepmc
import ROOT
import numpy as np
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

def get_nonradiative_decay(particle):
    if (particle.children[0].pid == particle.pid):
        return(get_nonradiative_decay(particle.children[0]))
    else:
        return(particle)

def fill_kinematic_array(particle, array):
    mom = particle.momentum
    array[:] = (mom.px, mom.py, mom.pz, mom.pt(), 
                np.sqrt(mom.px *mom.px + mom.py * mom.py + mom.pz * mom.pz), 
                mom.eta(), mom.e, mom.phi(), mom.m(), particle.pid, 
                particle.status)
    return(array)

# Files
ifile_name = "/data/home/ganowak/MG5_aMC_v2_9_25/WW_NLO_GaussStatement_1JetMatching/Events/run_01/events_PYTHIA8_0.hepmc" 
ofile_name = "Test1.root"
# Setting up Path
ww_path = at.find_WW_path()

# Inputs
W_pid = -24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:status/F'
wgt_str = "AUX_MUR05_MUF05/F:AUX_MUR05_MUF10/F:AUX_MUR05_MUF20/F:AUX_MUR10_MUF05/F:AUX_MUR10_MUF10/F:AUX_MUR10_MUF20/F:AUX_MUR20_MUF05/F:AUX_MUR20_MUF10/F:AUX_MUR20_MUF20/F"
evt_array = np.array([0], dtype=np.float32)
wgt_array = np.array([0]*9, dtype=np.float32)
target_particle_array = np.array([0]*11, dtype=np.float32)
target_lepton_array = np.array([0]*11, dtype=np.float32)
target_neutrino_array = np.array([0]*11, dtype=np.float32)
target_antiparticle_array = np.array([0]*11, dtype=np.float32)
target_antilepton_array = np.array([0]*11, dtype=np.float32)
target_antineutrino_array = np.array([0]*11, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/madgraph/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch('Weights', wgt_array, wgt_str)
tree.Branch('TargetParticle', target_particle_array, var_str)
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetNeutrino', target_neutrino_array, var_str)
tree.Branch('TargetAntiParticle', target_antiparticle_array, var_str)
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
tree.Branch('TargetAntiNeutrino', target_antineutrino_array, var_str)

# Variables
i = 0
counter1 = 0
counter2 = 0

with pyhepmc.open(ifile_name) as f:
    for event in f:
        evt_array[:] = i
        # Find top and anti-top in event.particles
        for particle in event.particles:
            if particle.pid==W_pid:
                wminus = particle
            elif particle.pid==-1*W_pid:
                wplus = particle
        # Find top b-daughter and w-daughter
        if (len(wplus.children) == 0) or (len(wminus.children) == 0):
            continue
        wplus = get_nonradiative_decay(wplus)
        wminus = get_nonradiative_decay(wminus)
        for daughter in wminus.children:
            if daughter.pid in lepton_pid_array:
                lepton_plus = daughter
            elif daughter.pid in -1*neutrino_pid_array:
                neutrino_plus = daughter
        for daughter in wplus.children:
            if daughter.pid in -1*lepton_pid_array:
                lepton_minus = daughter
            elif daughter.pid in neutrino_pid_array:
                neutrino_minus = daughter

        target_particle_array = fill_kinematic_array(wplus, target_particle_array)
        target_lepton_array = fill_kinematic_array(lepton_minus, target_lepton_array)
        target_neutrino_array = fill_kinematic_array(neutrino_minus, target_neutrino_array)
        target_antiparticle_array = fill_kinematic_array(wminus, target_antiparticle_array)
        target_antilepton_array = fill_kinematic_array(lepton_plus, target_antilepton_array)
        target_antineutrino_array = fill_kinematic_array(neutrino_plus, target_antineutrino_array)
        wgt_array[:] = event.weights[1:10]
        tree.Fill()
        i = i+1

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
tree.Print()
file.Write()
