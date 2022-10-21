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
DaVinci().CondDBtag = 'sim-20170721-2-vc-md100'
DaVinci().DDDBtag = 'dddb-20170721-3'
# DaVinci().CondDBtag = 'sim-20161124-2-vc-md100'
# DaVinci().DDDBtag = 'dddb-20150724'

# Use the local input data
WW_data = '00057387_00000003_3.AllStreams.dst'
tautau_data = '00144995_00000138_7.AllStreams.dst'
IOHelper().inputFiles([
    '~/WWProduction/Data/TestFiles/' + tautau_data
], clear=True)


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
ele_counter = 0
muon_counter = 0
hadron_counter = 0

# while (evtnum<evtmax):
while bool(tes['/Event']):
    evtnum += 1

    # Truth Stuff
    # print('New Event:')
    truth_particles = tes['MC/Particles']
    tau_list = [particle for particle in truth_particles if particle.particleID().pid()==15]
    decay_list = [particle.particleID().pid() for particle in tau_list[0].endVertices()[0].products()]
    # print(decay_list)
    if -12 in decay_list:
        ele_counter+=1
    elif 13 in decay_list:
        muon_counter+=1
    else:
        hadron_counter+=1
    gaudi.run(1) 

print(f'Electron Count: {ele_counter}')
print(f'Muon Count: {muon_counter}')
print(f'Hadron Count: {hadron_counter}')