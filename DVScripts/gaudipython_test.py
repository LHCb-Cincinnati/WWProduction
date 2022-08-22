# Imports
import sys
import pdb
from array import array

import numpy as np
import ROOT

import GaudiPython
from GaudiConf import IOHelper
from Configurables import DaVinci
from Configurables import CombineParticles
from Configurables import ApplicationMgr 
from PhysConf.Selections import Selection
from PhysConf.Selections import SelectionSequence
from PhysConf.Selections import CombineSelection
from StandardParticles import StdAllLooseMuons as Muons


# Configure DaVinci
DaVinci().InputType = 'DST'
DaVinci().TupleFile = '~/WWProduction/Data/DVTuples/Test.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = '2016'
DaVinci().Simulation = True
# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().EvtMax = -1
DaVinci().CondDBtag = 'sim-20161124-2-vc-md100'
DaVinci().DDDBtag = 'dddb-20150724'



# Use the local input data
IOHelper().inputFiles([
    '~/WWProduction/Data/TestFiles/00057387_00000003_3.AllStreams.dst'
], clear=True)

# Containers for ROOT Tree
#evt_num = np.array([0])
evt_num = array('f', [ 1.5 ])
l_plus_array = np.array(4*[0], dtype=np.float32)
l_minus_array = np.array(4*[0], dtype=np.float32)

# Create ROOT Tree
ofile = ROOT.TFile('~/WWProduction/Data/DVTuples/ofile.root', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')
tree.Branch('evt_num', evt_num, 'evt_num/F')
#tree.Branch('EventNumber', evt_num, 'EventNumber/D')
tree.Branch('l_plus_array', l_plus_array, 'l_plus_array[' + str(4) + ']/F')
tree.Branch('l_minus_array', l_minus_array, 'l_minus_array[' + str(4) + ']/F')


# dilepton Object Cuts
# Cuts on the individual muons, dilepton object, and parent particle.
# All the cuts here are fake to accept almost any lepton pair.
dilepton_decay_products = {
    'mu+': '(PT > 5000*MeV) & (ETA < 5) & (ETA > 2)',
    'mu-': '(PT > 5000*MeV) & (ETA < 5) & (ETA > 2)'
}
dilepton_comb = '(AM > 100*MeV)'
#dilepton_comb = "15*GeV<AMAXCHILD(MAXTREE('mu+'==ABSID,PT)"
dilepton_mother = '(MM > 100*MeV)'


# dilepton Object Definition
# Decay: Intermediate -> mu+ mu-
# Where Intermediate is a fake parent particle created only to get information
# on the leptons in the event.
dilepton_sel = CombineSelection(
    'dilepton_sel',
    [Muons],
    DecayDescriptor='Intermediate -> mu- mu+',
    DaughtersCuts=dilepton_decay_products,
    CombinationCut=dilepton_comb,
    MotherCut=dilepton_mother
)

dilepton_seq = SelectionSequence('dilepton_Seq', TopSelection=dilepton_sel)
DaVinci().appendToMainSequence([dilepton_seq.sequence()])

# Configure GaudiPython
gaudi = GaudiPython.AppMgr()
tes   = gaudi.evtsvc()

# Run
evtmax = 100
#evtmax = float('inf')
evtnum = 0
bevents = 0
gaudi.run(1)
#pdb.set_trace()
while evtnum < evtmax:
#while bool(tes['/Event']):
    # if not bool(tes['/Event']):
    #     print(f'EVENT HERE {evt_num}')
    #     pdb.set_trace()
    #     bevents+=1
    #     break
    # print(type(evt_num))
    # print(evt_num)
    if not bool(tes['/Event']):
        print(f'EVENT HERE {evt_num}')
        #pdb.set_trace()
        bevents+=1
    #print(bool(tes['/Event']))

    evtnum += 1
    if not bool(tes['Phys/StdAllLooseMuons/Particles']):
        gaudi.run(1)
        continue
    candidates =  tes[dilepton_seq.outputLocation()]
    for index in range(len(candidates)):
        print(candidates)
        evt_num[0] = tes['DAQ/ODIN'].eventNumber()
        daughter_id_list = [daughter.particleID().pid() for daughter in candidates[index].daughters()]
        l_minus_index = daughter_id_list.index(13)
        l_plus_index = daughter_id_list.index(-13)
        l_minus_array[:] = np.array((candidates[index].daughters()[l_minus_index].momentum().X(),
                            candidates[index].daughters()[l_minus_index].momentum().Y(),
                            candidates[index].daughters()[l_minus_index].momentum().Z(),
                            candidates[index].daughters()[l_minus_index].momentum().E()), dtype=np.float32)
        l_plus_array[:] = np.array((candidates[index].daughters()[l_plus_index].momentum().X(),
                            candidates[index].daughters()[l_plus_index].momentum().Y(),
                            candidates[index].daughters()[l_plus_index].momentum().Z(),
                            candidates[index].daughters()[l_plus_index].momentum().E()), dtype=np.float32)
        print(l_minus_array)                           
        print(l_plus_array)
        print(evt_num)
        tree.Fill()
    gaudi.run(1)

tree.Print()
ofile.Write()
ofile.Close()

