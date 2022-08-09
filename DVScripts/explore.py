import sys

import GaudiPython as GP
from GaudiConf import IOHelper
from Configurables import DaVinci

dv = DaVinci()
dv.DataType = '2016'
dv.Simulation = True

# Retrieve file path to open as the last command line argument
inputFiles = [sys.argv[-1]]
IOHelper('ROOT').inputFiles(inputFiles)

appMgr = GP.AppMgr()
evt = appMgr.evtsvc()

appMgr.run(1)

def advance(decision):
    """Advance until stripping decision is true, returns
    number of events by which we advanced"""
    n = 0
    while evt['/Event']:
        reports = evt['/Event/Strip/Phys/DecReports']
        if reports.hasDecisionName('Stripping{}Decision'.format(decision)):
            break

        appMgr.run(1)
        n += 1
    return n

def get_candidates(stripping_line):
    candidates = evt['/Event/AllStreams/Phys/{}/Particles'.format(stripping_line)]
    return(candidates)

def quick_explore(decision):
    n=0
    print_count = 0
    while (evt['/Event'] and print_count<100):
        reports = evt['/Event/Strip/Phys/DecReports']
        if reports.hasDecisionName('Stripping{}Decision'.format(decision)):
            candidates = get_candidates(decision)

            # Print the candidates for the event
            print("###########################\n")
            print(f"Event: {n}.  Number of Candidates: {len(candidates)}\n")
            for candidate in candidates:
                print(candidate)
            print_count+=1
            print("###########################\n")
        appMgr.run(1)
        n+=1
    return n
