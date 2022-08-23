# Imports
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

sys.path.append('/Users/nowakg/OneDrive - University of Cincinnati/Research/WeakBosonDecay/Analysis') # This is Awesome!
import analysis_tools as at

# Argument Parser
parser = argparse.ArgumentParser(description='''Process files and settings for analysis.
                                Required arguments are: input_files.
                                Optional arguments are: cross_section, testing.
                                Ex: python3 analysis.py test.root -c 53.5 -p''')
parser.add_argument('input_files',type=open, nargs='+',
                    help='The file or files to be anayzed. Input files should be root files.')
parser.add_argument('-c', '--cross_section', type=float, default=False, nargs='*', 
                    help='''The cross section used to normalize histograms in fb.
                    This option should either not be specified or have exactly
                    as many arguments as input files.''')
parser.add_argument('-t', '--testing', default=True, action='store_true',
                    help='''Flag to indicate if a new output folder should be created for 
                    this analysis.  -t means that all plots will go in folder labelled Test.''')
parser.add_argument('-p', '--production', dest='testing', action='store_false',
                    help='''Flag to indicate if a new output folder should be created for this 
                    analysis.  -p means that all plots will be put into a new output folder with
                    the same name as the input file.''')
args = parser.parse_args()
print(args)

# Make sure the cross-sections have as many aguments as files.
if (not args.cross_section):
    args.cross_section = [False] * len(args.input_files)
elif (len(args.cross_section) != len(args.input_files)):
    raise RuntimeError('''The number of given cross sections is different from
                       the number of given input files.  Please either don't
                       specify a cross section arugment or specify as many
                       cross sections as there are input files.''')

# Parse inputs
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['AnalysisTree'].arrays()

# # Create Vectors

lminus_vec = vector.zip({
    'px': tree['LeptonMinus'].px,
    'py': tree['LeptonMinus'].py,
    'pz': tree['LeptonMinus'].pz,
    'e': tree['LeptonMinus'].energy,
})

lplus_vec = vector.zip({
    'px': tree['LeptonPlus'].px,
    'py': tree['LeptonPlus'].py,
    'pz': tree['LeptonPlus'].pz,
    'e': tree['LeptonPlus'].energy,
})
dilepton_vec = lminus_vec + lplus_vec

# Cuts
lepton_eta_cut = ((lminus_vec.eta > 2) & (lminus_vec.eta<5)
                    & (lplus_vec.eta > 2) & (lplus_vec.eta<5))
high_pT_lepton_cut = ((lminus_vec.pt>15) | (lplus_vec.pt>15))
low_pT_lepton_cut = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
muon_pid_cut = ((tree['LeptonMinus'].id==13) & (tree['LeptonPlus'].id==-13))
invariant_mass_cut = (dilepton_vec.m>10)
lepton_cuts = lepton_eta_cut&high_pT_lepton_cut&low_pT_lepton_cut&muon_pid_cut
drellyan_cuts = lepton_eta_cut&low_pT_lepton_cut&invariant_mass_cut

# Calculate Quantitites
greater_pt_array = (lminus_vec.pt > lplus_vec.pt)
lesser_pt_array = -1*(greater_pt_array-1)
leading_lepton_pT_array = greater_pt_array * lminus_vec.pt + lesser_pt_array * lplus_vec.pt
trailing_lepton_pT_array = lesser_pt_array * lminus_vec.pt + greater_pt_array * lplus_vec.pt
delta_phi_array = np.abs(lminus_vec.deltaphi(lplus_vec))
delta_eta_array = np.abs(lminus_vec.deltaeta(lplus_vec))
delta_r_array = np.abs(lminus_vec.deltaR(lplus_vec))

# # Create output directory if it does not yet exist
# # and change current directory to output directory
#file_path = at.create_folder_path('WW_GenLevelPythia' + '.root', args.testing)

file_path = 'DY_GenLevelPythia'
#os.mkdir(file_path)
os.chdir(file_path)

# scale_factor
#cross_section = 3038 # in fb (WW)
cross_section = 35.53E6 # in fb (DY)
nevents = 100000
scale_factor = cross_section / nevents


# Plots
at.create_hist(dilepton_vec.m[drellyan_cuts], 'DiLepton Mass', bins=50, range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(dilepton_vec.pt[drellyan_cuts], 'Lepton Pair pT', yscale='log', bins=50, range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(leading_lepton_pT_array[drellyan_cuts], 'Leading Lepton pT', yscale='log', bins=50, range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(trailing_lepton_pT_array[drellyan_cuts], 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(lminus_vec.pt[drellyan_cuts], 'lminus pT', bins=50, yscale='log', range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(lplus_vec.pt[drellyan_cuts], 'lplus pT', yscale='log', bins=50, range=(0,150),
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(delta_phi_array[drellyan_cuts], 'Delta Phi', bins=50,
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(delta_r_array[drellyan_cuts], 'Delta R', bins=50,
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))
at.create_hist(delta_eta_array[drellyan_cuts], 'Delta Eta', bins=50,
            weights=[scale_factor]*len(dilepton_vec[drellyan_cuts]))