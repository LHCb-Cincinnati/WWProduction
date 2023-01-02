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
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2016'
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
tMuon_array = np.array(5*[0], dtype=np.float32)
# Reco Containers
rEvt_num = np.array([0], dtype=np.float32)
rMuon_array = np.array(5*[0], dtype=np.float32)
rMuon_id_array = np.array(10*[0], dtype=np.float32)
rMuon_deltaRmatch = np.array([0], dtype=np.float32)

# Create ROOT Tree
ofile_name = 'MuonEfficiencyTuplettbar.root'
ofile = ROOT.TFile(f'~/WWProduction/Data/DVTuples/{ofile_name}', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')

# Truth Branches
tree.Branch('tEvt_num', tEvt_num, 'tEvt_num/F')
tree.Branch('tMuon_array', tMuon_array, 'px/F:py/F:pz/F:e/F:pid/F')
# Reco Branches
tree.Branch('rEvt_num', rEvt_num, 'rEvt_num/F')
tree.Branch('rMuon_array', rMuon_array, 'px/F:py/F:pz/F:e/F:pid/F')
tree.Branch('rMuon_id_array', rMuon_id_array, 'ip/F:ipchi2/F:isMuon/F:isMuonLoose/F:isMuonTight/F:EcalE/F:HcalE/F:PRSE/F:TrackP/F:ConepT/F')
tree.Branch('rMuon_deltaRmatch', rMuon_deltaRmatch, 'rLeading_lepton_deltaRmatch/F')

# Configure GaudiPython
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()
pvrTool = gaudi.toolsvc().create(
    "GenericParticle2PVRelator<_p2PVWithIPChi2, "
    "OfflineDistanceCalculatorName>/P2PVWithIPChi2",
    interface = "IRelatedPVFinder")
dstTool = gaudi.toolsvc().create(
    "LoKi::TrgDistanceCalculator",
    interface = "IDistanceCalculator")

# Run
evtnum = 0
gaudi.run(1)
evtmax = 1000000
counter = 0

# while (evtnum<evtmax):
while bool(tes['/Event']) and (evtnum<evtmax):
    evtnum += 1

    # Initialize Truth Arrays
    tEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)

    # Truth Stuff
    truth_particles = tes['MC/Particles']
    decay_dict = {}
    truth_muon_dict = {particle.pt():particle for particle
                        in truth_particles
                        if abs(particle.particleID().pid())==13}
    truth_leading_muon = truth_muon_dict[max(truth_muon_dict.keys())]
    tMuon_array[:] = np.array((
                    truth_leading_muon.momentum().X(),
                    truth_leading_muon.momentum().Y(),
                    truth_leading_muon.momentum().Z(),
                    truth_leading_muon.momentum().E(),
                    truth_leading_muon.particleID().pid()), 
                    dtype=np.float32)
    
    # Fill MC Branches
    tree.GetBranch("tEvt_num").Fill()
    tree.GetBranch("tMuon_array").Fill()

    # Reco Stuff
    # If collection is not filled, skip the loop
    if (not bool(tes['Phys/StdAllLooseMuons/Particles'])):
        gaudi.run(1) 
        continue

    # Initialize Reco Arrays
    rEvt_num[:] = np.array(tes['DAQ/ODIN'].eventNumber(), dtype=np.float32)
    rMuon_id_array[:] = np.array(10*[0], dtype=np.float32)
    
    candidates =  tes['Phys/StdAllLooseMuons/Particles']
    primary_verts = tes['Rec/Vertex/Primary']
    reco_leptons_dict = {particle.pt():particle for particle in candidates}
    reco_leading_lepton = reco_leptons_dict[max(reco_leptons_dict.keys())]

    rMuon_array[:] = np.array((
                        reco_leading_lepton.momentum().X(),
                        reco_leading_lepton.momentum().Y(),
                        reco_leading_lepton.momentum().Z(),
                        reco_leading_lepton.momentum().E(),
                        reco_leading_lepton.particleID().pid()),
                        dtype=np.float32)
    rMuon_deltaRmatch[:] = np.array([tools.DeltaR(reco_leading_lepton, truth_leading_muon)], dtype=np.float32)
    # This is a catch-all for the one or two times the muonPID object
    # is undefined.
    if (not bool(reco_leading_lepton.proto().muonPID())):
        counter +=1
    else:
        # Find IP
        Muon_ip = np.array([0], dtype=np.double)
        Muon_ipchi2 = np.array([0], dtype=np.double)
        leading_lepton_primary_vert = pvrTool.relatedPV(reco_leading_lepton,
                                                        primary_verts)
        dstTool.distance(reco_leading_lepton, leading_lepton_primary_vert,
                         Muon_ip, Muon_ipchi2)

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
        rMuon_id_array[:] = np.array((Muon_ip,
                                    Muon_ipchi2,
                                    reco_leading_lepton.proto().muonPID().IsMuon(),
                                    reco_leading_lepton.proto().muonPID().IsMuonLoose(),
                                    reco_leading_lepton.proto().muonPID().IsMuonTight(),
                                    reco_leading_lepton.proto().extraInfo()[332],
                                    reco_leading_lepton.proto().extraInfo()[333],
                                    reco_leading_lepton.proto().extraInfo()[331],
                                    reco_leading_lepton.proto().extraInfo()[504],
                                    0),
                                    dtype=np.float32)

    tree.GetBranch("rEvt_num").Fill()                                
    tree.GetBranch("rMuon_array").Fill()
    tree.GetBranch("rMuon_id_array").Fill()
    tree.GetBranch("rMuon_deltaRmatch").Fill()

    # Continue on with the event loop
    gaudi.run(1) 


print(f'Count: {counter}')
tree.SetEntries(-1)
ofile.Write()
ofile.Close()