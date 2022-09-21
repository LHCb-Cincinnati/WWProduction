# Imports
import sys
import pdb

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
DaVinci().CondDBtag = 'sim-20161124-2-vc-md100'
DaVinci().DDDBtag = 'dddb-20150724'

# Use the local input data
IOHelper().inputFiles([
    '~/WWProduction/Data/TestFiles/00057387_00000003_3.AllStreams.dst'
], clear=True)

# Containers for ROOT Tree
evt_num = np.array([0], dtype=np.float32)
electron_array = np.array(4*[0], dtype=np.float32)
muon_array = np.array(4*[0], dtype=np.float32)

# Create ROOT Tree
ofile = ROOT.TFile('~/WWProduction/Data/DVTuples/ofile.root', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')
tree.Branch('evt_num', evt_num, 'evt_num/F')
tree.Branch('electron_array', electron_array, 'electron_array[' + str(4) + ']/F')
tree.Branch('muon_array', muon_array, 'muon_array[' + str(4) + ']/F')

# Set up DaVinci Objects
import DVoption_Sequences as DVSequences
dilepton_seq = DVSequences.GetDiLeptonSequence()
mc_particle_flow = DVSequences.GetMCParticleFlow()
mc_jet_builder = DVSequences.GetMCJetBuilder(mc_particle_flow)
hlt_particle_flow = DVSequences.GetHLTParticleFlow(dilepton_seq)
hlt_jet_builder = DVSequences.GetHLTJetBuilder(hlt_particle_flow)
DaVinci().appendToMainSequence([dilepton_seq, mc_particle_flow,
                                mc_jet_builder, hlt_particle_flow,
                                hlt_jet_builder])

# Configure GaudiPython
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

# Run
evtnum = 0
gaudi.run(1)
evtmax = 100
while evtnum<evtmax:
# while (bool(tes['/Event'])):
    # Putting this gaud.run(1) here skips the first event in preference of
    # simplicity.
    gaudi.run(1) 
    evtnum += 1

    # If collection is not filled, skip the loop
    if (not bool(tes['Phys/StdAllLooseMuons/Particles']) or not bool(tes['Phys/StdAllLooseElectrons/Particles'])):
        continue

    candidates =  tes[dilepton_seq.outputLocation()]
    mc_pfs = tes[mc_particle_flow.Output]
    mc_jets = tes[mc_jet_builder.Output]
    hlt_pfs = tes[hlt_particle_flow.Output]
    hlt_jets = tes[hlt_jet_builder.Output]
    for index in range(len(candidates)):
        pdb.set_trace()
        evt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
        daughter_id_list = [abs(daughter.particleID().pid()) for daughter in candidates[index].daughters()]
        muon_index = daughter_id_list.index(13)
        electron_index = daughter_id_list.index(11)
        muon_array[:] = np.array((candidates[index].daughters()[muon_index].momentum().X(),
                            candidates[index].daughters()[muon_index].momentum().Y(),
                            candidates[index].daughters()[muon_index].momentum().Z(),
                            candidates[index].daughters()[muon_index].momentum().E()), dtype=np.float32)
        electron_array[:] = np.array((candidates[index].daughters()[electron_index].momentum().X(),
                            candidates[index].daughters()[electron_index].momentum().Y(),
                            candidates[index].daughters()[electron_index].momentum().Z(),
                            candidates[index].daughters()[electron_index].momentum().E()), dtype=np.float32)
    
        tree.Fill()

ofile.Write()
ofile.Close()

