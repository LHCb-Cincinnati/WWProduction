# Imports
# Standard Library Packages
import sys
import yaml
import pdb
import inspect
# Standard Packages
import numpy as np
# HEP Packages
import ROOT
# LHCb Packages
import GaudiPython
from GaudiConf import IOHelper
from Configurables import DaVinci
from Configurables import ApplicationMgr 
# Personal Packages
import DVoption_Sequences as DVSequences
import DVoption_tools as tools

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
rLeading_lepton_id_array = np.array(9*[0], dtype=np.float32)
rTrailing_lepton_id_array = np.array(9*[0], dtype=np.float32)
rLeading_lepton_deltaRmatch = np.array([0], dtype=np.float32)
rTrailing_lepton_deltaRmatch = np.array([0], dtype=np.float32)

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
tree.Branch('rLeading_lepton_id_array', rLeading_lepton_id_array, 'ip/F:isMuon/F:isMuonLoose/F:isMuonTight/F:EcalE/F:HcalE/F:PRSE/F:TrackP/F:ConepT/F')
tree.Branch('rTrailing_lepton_id_array', rTrailing_lepton_id_array, 'ip/F:isMuon/F:isMuonLoose/F:isMuonTight/F:EcalE/F:HcalE/F:PRSE/F:TrackP/F:ConepT/F')
tree.Branch('rLeading_lepton_deltaRmatch', rLeading_lepton_deltaRmatch, 'rLeading_lepton_deltaRmatch/F')
tree.Branch('rTrailing_lepton_deltaRmatch', rTrailing_lepton_deltaRmatch, 'rTrailing_lepton_deltaRmatch/F')


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
SourceParticles = config_dict['SourceParticles']
counter = 0

# while (evtnum<evtmax):
while bool(tes['/Event']) and (evtnum<evtmax):
    evtnum += 1

    # Initialize Truth Arrays
    tEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
    tDecay_process_array[:] = np.array(3*[0], dtype=np.float32)

    # Truth Stuff
    truth_particles = tes['MC/Particles']
    decay_dict = {}
    for key in SourceParticles:
        source_pid = config_dict['SourceParticles'][key]['SourcePID']
        antiparticle_source = config_dict['SourceParticles'][key]['AntiParticleSource']
        num_min_products = config_dict['SourceParticles'][key]['NumMinProducts']
        target_pid_array = np.array(config_dict['SourceParticles'][key]['TargetPIDs'])
        particle = tools.FindDecayProduct(truth_particles, source_pid, num_min_products, target_pid_array)
        if antiparticle_source:
            antiparticle = tools.FindDecayProduct(truth_particles, -1*source_pid, num_min_products, -1*target_pid_array)    
        decay_dict[particle.pt()] = particle
        decay_dict[antiparticle.pt()] = antiparticle
    truth_leading_lepton = decay_dict[max(decay_dict.keys())]
    truth_trailing_lepton = decay_dict[min(decay_dict.keys())]

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
    if (not bool(tes['Phys/StdAllLooseMuons/Particles'])
        or not bool(tes['Phys/StdAllLooseElectrons/Particles'])
        or not bool(tes[dilepton_seq.outputLocation()])):
        gaudi.run(1) 
        continue

    candidates =  tes[dilepton_seq.outputLocation()]
    # mc_pfs = tes[mc_particle_flow.Output]
    # mc_jets = tes[mc_jet_builder.Output]
    # hlt_pfs = tes[hlt_particle_flow.Output]
    # hlt_jets = tes[hlt_jet_builder.Output]
    # reg_jets = tes['Phys/StdHltJets/Particles']

    candidates_dilepton_pT_sum = [sum([candidate[1].momentum().pt(),
                                    candidate[2].momentum().pt()])
                                    for candidate in candidates]
    best_candidate_index = candidates_dilepton_pT_sum.index(
                            max(candidates_dilepton_pT_sum))

    # pdb.set_trace()

    # Initialize Reco Arrays
    rEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
    rLeading_lepton_id_array[:] = np.array(9*[0], dtype=np.float32)
    rTrailing_lepton_id_array[:] = np.array(9*[0], dtype=np.float32)

    reco_leptons_dict = {particle.pt():particle for particle
                        in candidates[best_candidate_index].daughters()}
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
    rLeading_lepton_deltaRmatch[:] = np.array([tools.DeltaR(reco_leading_lepton, truth_leading_lepton)], dtype=np.float32)
    rTrailing_lepton_deltaRmatch[:] = np.array([tools.DeltaR(reco_trailing_lepton, truth_trailing_lepton)], dtype=np.float32)
    # This is a catch-all for the one or two times the muonPID object
    # is undefined.
    if (not bool(reco_leading_lepton.proto().muonPID()) or not bool(reco_trailing_lepton.proto().muonPID())):
        counter +=1
    else:
        # Just a reminder the lepton id array is currently filled as
        # impact parameter
        # isMuon
        # isMuonLoose
        # isMuonTight
        # energy deposited in ECal
        # energy deposited in HCal
        # Energy deposited in PRS
        # Track Momentum
        # pT in a cone of angle ...
        rLeading_lepton_id_array[:] = np.array((0,
                                    reco_leading_lepton.proto().muonPID().IsMuon(),
                                    reco_leading_lepton.proto().muonPID().IsMuonLoose(),
                                    reco_leading_lepton.proto().muonPID().IsMuonTight(),
                                    reco_leading_lepton.proto().extraInfo()[332],
                                    reco_leading_lepton.proto().extraInfo()[333],
                                    reco_leading_lepton.proto().extraInfo()[349],
                                    0,
                                    0),
                                    dtype=np.float32)
        rTrailing_lepton_id_array[:] = np.array((0,
                                    reco_trailing_lepton.proto().muonPID().IsMuon(),
                                    reco_trailing_lepton.proto().muonPID().IsMuonLoose(),
                                    reco_trailing_lepton.proto().muonPID().IsMuonTight(),
                                    reco_trailing_lepton.proto().extraInfo()[332],
                                    reco_trailing_lepton.proto().extraInfo()[333],
                                    reco_trailing_lepton.proto().extraInfo()[349],
                                    0,
                                    0),
                                    dtype=np.float32)

    tree.GetBranch("rEvt_num").Fill()                                
    tree.GetBranch("rLeading_lepton_array").Fill()
    tree.GetBranch("rTrailing_lepton_array").Fill()
    tree.GetBranch("rLeading_lepton_id_array").Fill()
    tree.GetBranch("rTrailing_lepton_id_array").Fill()
    tree.GetBranch("rLeading_lepton_deltaRmatch").Fill()
    tree.GetBranch("rTrailing_lepton_deltaRmatch").Fill()

    # Continue on with the event loop
    gaudi.run(1) 


print(f'Count: {counter}')
tree.SetEntries(-1)
ofile.Write()
ofile.Close()