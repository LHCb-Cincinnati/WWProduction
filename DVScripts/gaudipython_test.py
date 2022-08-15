# Imports
import sys
import pdb

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
evt_num = array('d', [0])
l_plus_array = np.array(4*[0])
l_minus_array = np.array(4*[0])

# Create ROOT Tree
ofile = ROOT.TFile('~/WWProduction/Data/DVTuples/ofile.root', 'RECREATE')
tree = ROOT.TTree('Tree', 'Tree')
tree.Branch('EventNumber', evt_num, 'EventNumber/D')
tree.Branch('l_plus_array', l_plus_array, 'l_plus_array[' + str(4) + ']/D')
tree.Branch('l_minus_array', l_minus_array, 'l_minus_array[' + str(4) + ']/D')


# dilepton Object Cuts
# Cuts on the individual muons, dilepton object, and parent particle.
# All the cuts here are fake to accept almost any lepton pair.
dilepton_decay_products = {
    'mu+': 'PT > 5000*MeV',
    'mu-': 'PT > 5000*MeV'
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
evtmax = 1000
evtnum = 0
while evtnum < evtmax:
    gaudi.run(1)
    evtnum += 1
    prts =  tes[dilepton_seq.outputLocation()]
    if prts:
        evt_num[0] = tes['DAQ/ODIN'].eventNumber()
        l_minus_array[:] = (prts[0].daughters()[0].momentum().X(),
                            prts[0].daughters()[0].momentum().Y(),
                            prts[0].daughters()[0].momentum().Z(),
                            prts[0].daughters()[0].momentum().E())
        l_plus_array[:] = (prts[0].daughters()[1].momentum().X(),
                           prts[0].daughters()[1].momentum().Y(),
                           prts[0].daughters()[1].momentum().Z(),
                           prts[0].daughters()[1].momentum().E())
        tree.Fill()

ofile.Write()
ofile.Close()

