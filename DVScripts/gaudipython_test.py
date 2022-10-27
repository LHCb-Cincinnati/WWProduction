# Imports
# Standard Library Packages
import sys
import yaml
import pdb
import inspect
# Standard Packages
import numpy as np
# Specific Packages
import ROOT
import GaudiPython
from GaudiConf import IOHelper
from Configurables import DaVinci
from Configurables import ApplicationMgr 
# Personal Packages
import DVoption_Sequences as DVSequences
from DVoption_tools import DeltaRMatching

# Parse Inputs
with open('DVScripts/config.yaml', 'r') as file:
    config_dict = yaml.safe_load(file)

# Configure DaVinci
DaVinci().InputType = 'DST'
DaVinci().PrintFreq = config_dict['PrintFrequency']
DaVinci().DataType = config_dict['DataType']
DaVinci().Simulation = True
# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().CondDBtag = config_dict['CondDBtag']
DaVinci().DDDBtag = config_dict['DDDBtag']

# Use the local input data
IOHelper().inputFiles([
    '~/WWProduction/Data/TestFiles/' + config_dict['DataFile']
], clear=True)

# Containers for ROOT Tree
# Truth Containers
tEvt_num = np.array([0], dtype=np.float32)
tDecay_process_array = np.array(3*[0], dtype=np.float32)
tLeading_lepton_array = np.array(5*[0], dtype=np.float32)
tTrailing_lepton_array = np.array(5*[0], dtype=np.float32)
# Reco Containers
rEvt_num = np.array([0], dtype=np.float32)
rLeading_lepton_array = np.array(5*[0], dtype=np.float32)
rTrailing_lepton_array = np.array(5*[0], dtype=np.float32)
rMuonid_array = np.array(2*[0], dtype=np.float32)
rElectronid_array = np.array(2*[0], dtype=np.float32)

# Create ROOT Tree
ofile_name = config_dict['OutputFile']
ofile = ROOT.TFile(f'~/WWProduction/Data/DVTuples/{ofile_name}', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')

# Truth Branches
tree.Branch('tEvt_num', tEvt_num, 'tEvt_num/F')
tree.Branch('tLeading_lepton_array', tLeading_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('tTrailing_lepton_array', tTrailing_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('tDecay_process_array', tDecay_process_array, 'mumu/F:mue/F:ee/F')
# Reco Branches
tree.Branch('rEvt_num', rEvt_num, 'rEvt_num/F')
tree.Branch('rLeading_lepton_array', rLeading_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('rTrailing_lepton_array', rTrailing_lepton_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('rMuonid_array', rMuonid_array, 'leading/F:trailing/F')
tree.Branch('rElectronid_array', rElectronid_array, 'leading/F:trailing/F')

# Set up DaVinci Objects
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
evtmax = config_dict['EvtMax']
counter = 0

while (evtnum<evtmax):
# while bool(tes['/Event']):
    evtnum += 1

    # Initialize Truth Arrays
    tEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
    tDecay_process_array[:] = np.array(3*[0], dtype=np.float32)

    # Truth Stuff
    truth_particles = tes['MC/Particles']
    W_list = [particle for particle in truth_particles
             if ((abs(particle.particleID().pid())==6)
             and particle.originVertex().type()==1)]
    pdb.set_trace()
    if len(W_list)!=2:
        pdb.set_trace()
    W1_decay_dict = {particle_ref.pt():particle_ref.target()
                for particle_ref in W_list[0].endVertices()[0].products()
                    if (abs(particle_ref.particleID().pid())==11
                        or abs(particle_ref.particleID().pid())==13)
                        and (W_list[0].particleID().pid()
                            / particle_ref.particleID().pid() < 0)}
    W2_decay_dict = {particle_ref.pt():particle_ref.target()
            for particle_ref in W_list[1].endVertices()[0].products()
                if (abs(particle_ref.particleID().pid())==11
                    or abs(particle_ref.particleID().pid())==13)
                    and (W_list[1].particleID().pid()
                        / particle_ref.particleID().pid() < 0)}
    if (not W1_decay_dict) or (not W2_decay_dict):
        gaudi.run(1)
        continue
    high_pT_W1_decay_product = W1_decay_dict[max(W1_decay_dict.keys())]
    high_pT_W2_decay_product = W2_decay_dict[max(W2_decay_dict.keys())]
    Ws_decay_dict = {high_pT_W1_decay_product.pt(): high_pT_W1_decay_product,
                    high_pT_W2_decay_product.pt(): high_pT_W2_decay_product}
    truth_leading_lepton = Ws_decay_dict[max(Ws_decay_dict.keys())]
    truth_trailing_lepton = Ws_decay_dict[min(Ws_decay_dict.keys())]

    if ((abs(truth_leading_lepton.particleID().pid()) == 13) and (abs(truth_trailing_lepton.particleID().pid()) == 13)):
        tDecay_process_array[:] = np.array([1,0,0], dtype=np.float32)
    elif ((abs(truth_leading_lepton.particleID().pid()) == 11) and (abs(truth_trailing_lepton.particleID().pid()) == 13)
        or (abs(truth_leading_lepton.particleID().pid()) == 13) and (abs(truth_trailing_lepton.particleID().pid()) == 11)):
        tDecay_process_array[:] = np.array([0,1,0], dtype=np.float32)
    elif ((abs(truth_leading_lepton.particleID().pid()) == 11) and (abs(truth_trailing_lepton.particleID().pid()) == 11)):
        tDecay_process_array[:] = np.array([0,0,1], dtype=np.float32)
    
    tLeading_lepton_array[:] = np.array((
                    truth_leading_lepton.momentum().X(),
                    truth_leading_lepton.momentum().Y(),
                    truth_leading_lepton.momentum().Z(),
                    truth_leading_lepton.momentum().E(),
                    truth_leading_lepton.particleID().pid()), 
                    dtype=np.float32)
    tTrailing_lepton_array[:] = np.array((
                    truth_trailing_lepton.momentum().X(),
                    truth_trailing_lepton.momentum().Y(),
                    truth_trailing_lepton.momentum().Z(),
                    truth_trailing_lepton.momentum().E(),
                    truth_trailing_lepton.particleID().pid()), 
                    dtype=np.float32)

    
    # Fill MC Branches
    tree.GetBranch("tEvt_num").Fill()
    tree.GetBranch("tLeading_lepton_array").Fill()
    tree.GetBranch("tTrailing_lepton_array").Fill()
    tree.GetBranch("tDecay_process_array").Fill()

    # Reco Stuff
    # If collection is not filled, skip the loop
    if (not bool(tes['Phys/StdAllLooseMuons/Particles']) or not bool(tes['Phys/StdAllLooseElectrons/Particles'])):
        gaudi.run(1) 
        continue

    candidates =  tes[dilepton_seq.outputLocation()]
    # mc_pfs = tes[mc_particle_flow.Output]
    # mc_jets = tes[mc_jet_builder.Output]
    # hlt_pfs = tes[hlt_particle_flow.Output]
    # hlt_jets = tes[hlt_jet_builder.Output]
    # reg_jets = tes['Phys/StdHltJets/Particles']

    for index in range(len(candidates)):

        # Initialize Reco Arrays
        rEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
        rMuonid_array[:] = np.array(2*[0], dtype=np.float32)
        rElectronid_array[:] = np.array(2*[0], dtype=np.float32)

        reco_leptons_dict = {particle.pt():particle for particle
                            in candidates[index].daughters()}
        reco_leading_lepton = reco_leptons_dict[max(reco_leptons_dict.keys())]
        reco_trailing_lepton = reco_leptons_dict[min(reco_leptons_dict.keys())]
        # test_ele = genTool.relatedMCPs(candidates[index].daughters()[electron_index])
        # test_muon = genTool.relatedMCPs(candidates[index].daughters()[muon_index])

        rLeading_lepton_array[:] = np.array((
                            reco_leading_lepton.momentum().X(),
                            reco_leading_lepton.momentum().Y(),
                            reco_leading_lepton.momentum().Z(),
                            reco_leading_lepton.momentum().E(),
                            reco_leading_lepton.particleID().pid()),
                            dtype=np.float32)
        rTrailing_lepton_array[:] = np.array((
                            reco_trailing_lepton.momentum().X(),
                            reco_trailing_lepton.momentum().Y(),
                            reco_trailing_lepton.momentum().Z(),
                            reco_trailing_lepton.momentum().E(),
                            reco_trailing_lepton.particleID().pid()),
                            dtype=np.float32)
        if (not bool(reco_leading_lepton.proto().muonPID()) or not bool(reco_trailing_lepton.proto().muonPID())):
            rMuonid_array[:] = np.array((0,0), dtype=np.float32)
            counter +=1
        else:
            rMuonid_array[:] = np.array((reco_leading_lepton.proto().muonPID().IsMuon(),
                                        reco_trailing_lepton.proto().muonPID().IsMuon()),
                                        dtype=np.float32)

        tree.GetBranch("rEvt_num").Fill()                                
        tree.GetBranch("rLeading_lepton_array").Fill()
        tree.GetBranch("rTrailing_lepton_array").Fill()
        tree.GetBranch("rMuonid_array").Fill()
    gaudi.run(1) 


print(f'Count: {counter}')
tree.SetEntries(-1)
ofile.Write()
ofile.Close()