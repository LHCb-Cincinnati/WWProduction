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
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tuple/DecayTree'].arrays()

# Create Indices
first_lepton_key = 'muon'
second_lepton_key = 'electron'

# Create Vectors
lminus_vec = vector.zip({
    'px': tree[first_lepton_key + '_PX'],
    'py': tree[first_lepton_key + '_PY'],
    'pz': tree[first_lepton_key + '_PZ'],
    'e': tree[first_lepton_key + '_PE'],
})

lplus_vec = vector.zip({
    'px': tree[second_lepton_key + '_PX'],
    'py': tree[second_lepton_key + '_PY'],
    'pz': tree[second_lepton_key + '_PZ'],
    'e': tree[second_lepton_key + '_PE'],
})
dilepton_vec = lminus_vec + lplus_vec


# Calculate Quantitites
greater_pt_array = (lminus_vec.pt > lplus_vec.pt)
lesser_pt_array = -1*(greater_pt_array-1)
leading_lepton_pT_array = greater_pt_array * lminus_vec.pt + lesser_pt_array * lplus_vec.pt
trailing_lepton_pT_array = lesser_pt_array * lminus_vec.pt + greater_pt_array * lplus_vec.pt
delta_phi_array = np.abs(lminus_vec.deltaphi(lplus_vec))
delta_eta_array = np.abs(lminus_vec.deltaeta(lplus_vec))
delta_r_array = np.abs(lminus_vec.deltaR(lplus_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
at.create_hist(dilepton_vec.m, 'DiLepton Mass', bins=50, range=(0,150000),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(dilepton_vec.pt, 'Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(leading_lepton_pT_array, 'Leading Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(leading_lepton_pT_array, cross_section))
at.create_hist(trailing_lepton_pT_array, 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(trailing_lepton_pT_array, cross_section))
#at.create_hist(lminus_vec.pt, 'lminus pT', bins=50, yscale='log', range=(0,150000),
#            weights=at.calculate_weights(lminus_vec, cross_section))
#at.create_hist(lplus_vec.pt, 'lplus pT', yscale='log', bins=50, range=(0,150000),
#            weights=at.calculate_weights(lplus_vec, cross_section))
at.create_hist(delta_phi_array, 'Delta Phi', bins=50,
            weights=at.calculate_weights(delta_phi_array, cross_section))
at.create_hist(delta_r_array, 'Delta R', bins=50,
            weights=at.calculate_weights(delta_r_array, cross_section))
at.create_hist(delta_eta_array, 'Delta Eta', bins=50,
            weights=at.calculate_weights(delta_eta_array, cross_section))