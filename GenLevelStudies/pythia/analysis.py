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

# Functions
# def get_px(pT_array, eta, phi):
#     from numpy import cos
#     return(pT_array*cos(phi))

# def get_py(pT_array, eta, phi):
#     from numpy import sin
#     return(pT_array*sin(phi))

# def get_pz(pT_array, eta, phi):
#     from numpy import sinh
#     return(pT_array*sinh(eta))

# Argument Parser
parser = argparse.ArgumentParser(description='Process files and settings for analysis.')
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
# lminus_array = tree['LeptonMinus']
# lplus_array = tree['LeptonPlus']
# dilepton_xmom = (get_px(lminus_array.pT, lminus_array.eta, lminus_array.phi)
#                 + get_px(lplus_array.pT, lminus_array.eta, lplus_array.phi))
# dilepton_ymom = (get_py(lminus_array.pT, lminus_array.eta, lminus_array.phi)
#                 + get_py(lplus_array.pT, lminus_array.eta, lplus_array.phi))
# dilepton_zmom = (get_pz(lminus_array.pT, lminus_array.eta, lminus_array.phi)
#                 + get_pz(lplus_array.pT, lminus_array.eta, lplus_array.phi))
# dilepton_pT_array = np.sqrt(dilepton_xmom**2 + dilepton_ymom**2)
# dilepton_momentum2 = dilepton_xmom**2 + dilepton_ymom**2 + dilepton_zmom**2
# dilepton_energy2 = (lminus_array.energy + lminus_array.energy)**2
# dilepton_mass_array = np.sqrt(dilepton_energy2 - dilepton_momentum2)

# Cuts
lepton_eta_cut = ((lminus_vec.eta > 2) & (lminus_vec.eta<5)
                    & (lplus_vec.eta > 2) & (lplus_vec.eta<5))
high_pT_lepton_cut = ((lminus_vec.pt>15) | (lplus_vec.pt>15))
low_pT_lepton_cut = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
muon_pid_cut = ((tree['LeptonMinus'].id==13) & (tree['LeptonPlus'].id==-13))
lepton_cuts = lepton_eta_cut&high_pT_lepton_cut&low_pT_lepton_cut&muon_pid_cut

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
file_path = at.create_folder_path('WW_GenLevelPythia' + '.root', args.testing)
os.chdir(file_path)

# # Plots
# at.create_hist(dilepton_vec.m, 'DiLepton Mass', bins=50, range=(0,150000),
#             weights=at.calculate_weights(cross_section, dilepton_vec))
# at.create_hist(dilepton_vec.pt, 'Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(cross_section, dilepton_vec))
# at.create_hist(leading_lepton_pT_array, 'Leading Lepton pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(cross_section, leading_lepton_pT_array))
# at.create_hist(trailing_lepton_pT_array, 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(cross_section, trailing_lepton_pT_array))
# at.create_hist(lminus_vec.pt, 'lminus pT', bins=50, yscale='log', range=(0,150000),
#             weights=at.calculate_weights(cross_section, lminus_vec))
# at.create_hist(lplus_vec.pt, 'lplus pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(cross_section, lplus_vec))
# at.create_hist(delta_phi_array, 'Delta Phi', bins=50,
#             weights=at.calculate_weights(cross_section, delta_phi_array))
# at.create_hist(delta_r_array, 'Delta R', bins=50,
#             weights=at.calculate_weights(cross_section, delta_r_array))
# at.create_hist(delta_eta_array, 'Delta Eta', bins=50,
#             weights=at.calculate_weights(cross_section, delta_eta_array))