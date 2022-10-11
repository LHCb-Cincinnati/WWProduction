# Imports
import sys
import os
import argparse 
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

import analysis_tools as at


# Parse inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tree']

# Create Indices
truth_leading_lepton_key = 'tLeading_lepton_array'
truth_trailing_lepton_key = 'tTrailing_lepton_array'

# Create Vectors
tLeading_lepton_vec = vector.zip({
    'px': tree[truth_leading_lepton_key].array()['px'],
    'py': tree[truth_leading_lepton_key].array()['py'],
    'pz': tree[truth_leading_lepton_key].array()['pz'],
    'e': tree[truth_leading_lepton_key].array()['e']
})

tTrailing_lepton_vec = vector.zip({
    'px': tree[truth_trailing_lepton_key].array()['px'],
    'py': tree[truth_trailing_lepton_key].array()['py'],
    'pz': tree[truth_trailing_lepton_key].array()['pz'],
    'e': tree[truth_trailing_lepton_key].array()['e']
})
tDilepton_vec = tLeading_lepton_vec + tTrailing_lepton_vec


# Calculate Quantitites
# greater_pt_array = (lminus_vec.pt > lplus_vec.pt)
# lesser_pt_array = -1*(greater_pt_array-1)
# leading_lepton_pT_array = greater_pt_array * lminus_vec.pt + lesser_pt_array * lplus_vec.pt
# trailing_lepton_pT_array = lesser_pt_array * lminus_vec.pt + greater_pt_array * lplus_vec.pt
tDelta_phi_array = np.abs(tLeading_lepton_vec.deltaphi(tTrailing_lepton_vec))
tDelta_eta_array = np.abs(tLeading_lepton_vec.deltaeta(tTrailing_lepton_vec))
tDelta_r_array = np.abs(tLeading_lepton_vec.deltaR(tTrailing_lepton_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
at.create_hist(tDilepton_vec.m, 'Truth DiLepton Mass', bins=50, range=(0,170000),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tDilepton_vec.pt, 'Truth Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tLeading_lepton_vec.pt, 'Truth Leading Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tLeading_lepton_vec.pt, cross_section))
at.create_hist(tTrailing_lepton_vec.pt, 'Truth Trailing Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tTrailing_lepton_vec.pt, cross_section))
at.create_hist(tLeading_lepton_vec.eta, 'Truth Leading Lepton Eta', bins=50,
            weights=at.calculate_weights(tLeading_lepton_vec, cross_section))
at.create_hist(tTrailing_lepton_vec.eta, 'Truth Trailing Lepton Eta', bins=50,
            weights=at.calculate_weights(tTrailing_lepton_vec, cross_section))
at.create_hist(tDelta_phi_array, 'Truth Delta Phi', bins=50,
            weights=at.calculate_weights(tDelta_phi_array, cross_section))
at.create_hist(tDelta_r_array, 'Truth Delta R', bins=50,
            weights=at.calculate_weights(tDelta_r_array, cross_section))
at.create_hist(tDelta_eta_array, 'Truth Delta Eta', bins=50,
            weights=at.calculate_weights(tDelta_eta_array, cross_section))