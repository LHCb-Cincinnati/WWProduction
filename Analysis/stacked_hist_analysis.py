# Imports
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

import analysis_tools as at


# Parse inputs
args = at.parse_user_input(sys.argv)
if len(args.input_files) < 2:
    raise RuntimeError('Not enough input files provided.')

file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section # Cross section list in fb

# Open the file
ifile_list = [uproot.open(file) for file in file_name_list]
tree_list = [file['Tuple/DecayTree'].arrays() for file in ifile_list]

# Create Indices
first_lepton_key = 'muon'
second_lepton_key = 'electron'

# Create Vectors
for index in range(len(tree_list)):
    lminus_vec_list[index] = vector.zip({
        'px': tree_list[index][first_lepton_key + '_PX'],
        'py': tree_list[index][first_lepton_key + '_PY'],
        'pz': tree_list[index][first_lepton_key + '_PZ'],
        'e': tree_list[index][first_lepton_key + '_PE'],
    })

    lplus_vec_list[index] = vector.zip({
        'px': tree_list[index][second_lepton_key + '_PX'],
        'py': tree_list[index][second_lepton_key + '_PY'],
        'pz': tree_list[index][second_lepton_key + '_PZ'],
        'e': tree_list[index][second_lepton_key + '_PE'],
    })

    dilepton_vec[index] = lminus_vec[index] + lplus_vec[index]


    # Calculate Quantitites
    greater_pt_array[index] = (lminus_vec[index].pt > lplus_vec[index].pt)
    lesser_pt_array[index] = -1*(greater_pt_array[index]-1)
    leading_lepton_pT_array[index] = greater_pt_array[index] * lminus_vec[index].pt + lesser_pt_array[index] * lplus_vec[index].pt
    trailing_lepton_pT_array[index] = lesser_pt_array[index] * lminus_vec[index].pt + greater_pt_array[index] * lplus_vec[index].pt
    delta_phi_array[index] = np.abs(lminus_vec[index].deltaphi(lplus_vec[index]))
    delta_eta_array[index] = np.abs(lminus_vec[index].deltaeta(lplus_vec[index]))
    delta_r_array[index] = np.abs(lminus_vec[index].deltaR(lplus_vec[index]))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path('StackedHist.root', args.testing)
os.chdir(file_path)

# Plots
# at.create_hist(dilepton_vec.m, 'DiLepton Mass', bins=50, range=(0,150000),
#             weights=at.calculate_weights(dilepton_vec, cross_section))
# at.create_hist(dilepton_vec.pt, 'Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(dilepton_vec, cross_section))
# at.create_hist(leading_lepton_pT_array, 'Leading Lepton pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(leading_lepton_pT_array, cross_section))
# at.create_hist(trailing_lepton_pT_array, 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(trailing_lepton_pT_array, cross_section))
# at.create_hist(lminus_vec.pt, 'Muon pT', bins=50, yscale='log', range=(0,150000),
#             weights=at.calculate_weights(lminus_vec, cross_section))
# at.create_hist(lplus_vec.pt, 'Electron pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(lplus_vec, cross_section))
# at.create_hist(delta_phi_array, 'Delta Phi', bins=50,
#             weights=at.calculate_weights(delta_phi_array, cross_section))
# at.create_hist(delta_r_array, 'Delta R', bins=50,
#             weights=at.calculate_weights(delta_r_array, cross_section))
# at.create_hist(delta_eta_array, 'Delta Eta', bins=50,
#             weights=at.calculate_weights(delta_eta_array, cross_section))