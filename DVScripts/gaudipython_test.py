# Imports
import sys
import pdb
import inspect
from collections import OrderedDict

import numpy as np
import ROOT

import GaudiPython
from GaudiConf import IOHelper
from Configurables import DaVinci
from Configurables import ApplicationMgr 

# Configure DaVinci
DaVinci().InputType = 'DST'
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2016'
DaVinci().Simulation = True
# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
# DaVinci().CondDBtag = 'sim-20170721-2-vc-md100'
# DaVinci().DDDBtag = 'dddb-20170721-3'
DaVinci().CondDBtag = 'sim-20161124-2-vc-md100'
DaVinci().DDDBtag = 'dddb-20150724'

# Use the local input data
WW_data = '00057387_00000003_3.AllStreams.dst'
tautau_data = '00144995_00000138_7.AllStreams.dst'
IOHelper().inputFiles([
    '~/WWProduction/Data/TestFiles/' + WW_data
], clear=True)

# Containers for ROOT Tree
tEvt_num = np.array([0], dtype=np.float32)
# Truth Containers
tFew_lepton_flag = np.array([0], dtype=np.float32)
tMany_lepton_flag = np.array([0], dtype=np.float32)
tLeading_lepton_array = np.array(5*[0], dtype=np.float32)
tTrailing_lepton_array = np.array(5*[0], dtype=np.float32)
# Reco Containers
rElectron_array = np.array(5*[0], dtype=np.float32)
rMuon_array = np.array(5*[0], dtype=np.float32)

# Create ROOT Tree
ofile = ROOT.TFile('~/WWProduction/Data/DVTuples/ofile.root', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')
tree.Branch('tEvt_num', tEvt_num, 'tEvt_num/F')
# Truth Branches
tree.Branch('tFew_lepton_flag', tFew_lepton_flag, 'tFew_lepton_flag/F')
tree.Branch('tMany_lepton_flag',tMany_lepton_flag, 'tMany_lepton_flag/F')
tree.Branch('tLeading_lepton_array', tLeading_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('tTrailing_lepton_array', tTrailing_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
# Reco Branches
tree.Branch('rElectron_array', rElectron_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('rMuon_array', rMuon_array, 'px/F:py/F:pz/F:e/F:pid/F')

# Set up DaVinci Objects
import DVoption_Sequences as DVSequences
from DVoption_tools import DeltaRMatching

dilepton_seq = DVSequences.GetDiLeptonSequence()
# mc_particle_flow = DVSequences.GetMCParticleFlow()
# mc_jet_builder = DVSequences.GetMCJetBuilder(mc_particle_flow)
# hlt_particle_flow = DVSequences.GetHLTParticleFlow(dilepton_seq)
# hlt_jet_builder = DVSequences.GetHLTJetBuilder(hlt_particle_flow)
DaVinci().appendToMainSequence([dilepton_seq])
# DaVinci().appendToMainSequence([dilepton_seq, mc_particle_flow,
#                                 mc_jet_builder, hlt_particle_flow,
#                                 hlt_jet_builder])

# Configure GaudiPython
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

genTool = gaudi.toolsvc().create(
    'DaVinciSmartAssociator',
    interface = 'IParticle2MCWeightedAssociator')

# Run
evtnum = 0
gaudi.run(1)
evtmax = 100

# Testing
counter = 0

while (evtnum<evtmax):
# while (bool(tes['/Event'])):
    # Putting this gaud.run(1) here skips the first event in preference of
    # simplicity.
    gaudi.run(1) 
    evtnum += 1

    # Initialize Arrays
    tEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
    tFew_lepton_flag[:] = np.array(0, dtype=np.float32)
    tMany_lepton_flag[:] = np.array(0, dtype=np.float32)

    # Truth Stuff
    truth_particles = tes['MC/Particles']
    truth_highptleptons_dict = {mc_particle.pt():mc_particle
                    for mc_particle in truth_particles
                    if (((abs(mc_particle.particleID().pid()) == 13) 
                        or (abs(mc_particle.particleID().pid()) == 11))
                        and (mc_particle.pt() > 10000))}
    if (len(truth_highptleptons_dict)>2):
        tMany_lepton_flag[:] = np.array(1, dtype=np.float32)
    elif (len(truth_highptleptons_dict)<2):
        # pdb.set_trace()
        tFew_lepton_flag[:] = np.array(1, dtype=np.float32)
        tree.GetBranch("tEvt_num").Fill()
        tree.GetBranch("tFew_lepton_flag").Fill()
        tree.GetBranch("tMany_lepton_flag").Fill()
        continue
    high_pt_vals = sorted(truth_highptleptons_dict.keys(), reverse=True)[:2]
    leading_lepton = truth_highptleptons_dict[high_pt_vals[0]]
    trailing_lepton = truth_highptleptons_dict[high_pt_vals[1]]
    tLeading_lepton_array[:] = np.array((
                    leading_lepton.momentum().X(),
                    leading_lepton.momentum().Y(),
                    leading_lepton.momentum().Z(),
                    leading_lepton.momentum().E(),
                    leading_lepton.particleID().pid()), 
                    dtype=np.float32)
    tTrailing_lepton_array[:] = np.array((
                    trailing_lepton.momentum().X(),
                    trailing_lepton.momentum().Y(),
                    trailing_lepton.momentum().Z(),
                    trailing_lepton.momentum().E(),
                    trailing_lepton.particleID().pid()), 
                    dtype=np.float32)

    # Fill MC Branches
    tree.GetBranch("tEvt_num").Fill()
    tree.GetBranch("tFew_lepton_flag").Fill()
    tree.GetBranch("tMany_lepton_flag").Fill()
    tree.GetBranch("tLeading_lepton_array").Fill()
    tree.GetBranch("tTrailing_lepton_array").Fill()

    # Reco Stuff
    # If collection is not filled, skip the loop
    if (not bool(tes['Phys/StdAllLooseMuons/Particles']) or not bool(tes['Phys/StdAllLooseElectrons/Particles'])):
        continue

    candidates =  tes[dilepton_seq.outputLocation()]
    # mc_pfs = tes[mc_particle_flow.Output]
    # mc_jets = tes[mc_jet_builder.Output]
    # hlt_pfs = tes[hlt_particle_flow.Output]
    # hlt_jets = tes[hlt_jet_builder.Output]
    # reg_jets = tes['Phys/StdHltJets/Particles']

    for index in range(len(candidates)):
        daughter_id_list = [abs(daughter.particleID().pid()) for daughter in candidates[index].daughters()]
        muon_index = daughter_id_list.index(13)
        electron_index = daughter_id_list.index(11)
        # test_ele = genTool.relatedMCPs(candidates[index].daughters()[electron_index])
        # test_muon = genTool.relatedMCPs(candidates[index].daughters()[muon_index])
        # pdb.set_trace()

        rMuon_array[:] = np.array((candidates[index].daughters()[muon_index].momentum().X(),
                            candidates[index].daughters()[muon_index].momentum().Y(),
                            candidates[index].daughters()[muon_index].momentum().Z(),
                            candidates[index].daughters()[muon_index].momentum().E(),
                            candidates[index].daughters()[muon_index].particleID().pid()),
                            dtype=np.float32)
        rElectron_array[:] = np.array((candidates[index].daughters()[electron_index].momentum().X(),
                            candidates[index].daughters()[electron_index].momentum().Y(),
                            candidates[index].daughters()[electron_index].momentum().Z(),
                            candidates[index].daughters()[electron_index].momentum().E(),
                            candidates[index].daughters()[electron_index].particleID().pid()),
                            dtype=np.float32)
        tree.GetBranch("rElectron_array").Fill()
        tree.GetBranch("rMuon_array").Fill()


tree.SetEntries(-1)
ofile.Write()
ofile.Close()

