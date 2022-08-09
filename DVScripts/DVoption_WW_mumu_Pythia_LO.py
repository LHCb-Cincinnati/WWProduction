# Imports
from GaudiConf import IOHelper

from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import MCDecayTreeTuple
from Configurables import DaVinci
from Configurables import CheckPV
from Configurables import TupleToolKinematic
from Configurables import TupleToolPid
from Configurables import TupleToolPrimaries
from Configurables import CombineParticles

from PhysConf.Selections import Selection
from PhysConf.Selections import SelectionSequence
from PhysConf.Selections import CombineSelection

from StandardParticles import StdAllLooseMuons as Muons


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

# Decay Tree Configuration
dtt = DecayTreeTuple('Tuple')
dtt.Inputs = dilepton_seq.outputLocations()
dtt.Decay = 'Intermediate -> ^mu+ ^mu-'
dtt.addBranches({
    'Intermediate': 'Intermediate -> mu+ mu-',
    'lplus': 'Intermediate -> ^mu+ mu-',
    'lminus': 'Intermediate -> mu+ ^mu-',
})

# Configure DaVinci
DaVinci().UserAlgorithms += [dilepton_seq.sequence(), dtt]
DaVinci().InputType = 'DST'
DaVinci().TupleFile = '~/WWProduction/Data/DVTuples/WW_mumu_Pythia_LO_Tuple.root'
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
