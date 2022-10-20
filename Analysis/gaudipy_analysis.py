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
reco_leading_lepton_key = 'rLeading_lepton_array'
reco_trailing_lepton_key = 'rTrailing_lepton_array'

# Masks
process_mask = (tree['tDecay_process_array'].array()['mue'] == 1)
good_reco_event = [evt for evt in tEvt_num if evt in rEvt_num and process_mask[list(tEvt_num).index(evt)]]

full_mask = process_mask

# Create Vectors
tLeading_lepton_vec = vector.awk(tree[truth_leading_lepton_key].array())
tTrailing_lepton_vec = vector.awk(tree[truth_trailing_lepton_key].array())
tDilepton_vec = tLeading_lepton_vec + tTrailing_lepton_vec
rLeading_lepton_vec = vector.awk(tree[reco_leading_lepton_key].array())
rTrailing_lepton_vec = vector.awk(tree[reco_trailing_lepton_key].array())
rDilepton_vec = rLeading_lepton_vec + rTrailing_lepton_vec
pdb.set_trace()

# Mask Vectors and Extract Lepton Info
tLeading_lepton_vec = tLeading_lepton_vec[full_mask]
tTrailing_lepton_vec = tTrailing_lepton_vec[full_mask]
tDilepton_vec = tDilepton_vec[full_mask]
tMuon_vec = ak.where((abs(tLeading_lepton_vec.pid)==13), tLeading_lepton_vec, tTrailing_lepton_vec)
tElectron_vec = ak.where((abs(tLeading_lepton_vec.pid)==11), tLeading_lepton_vec, tTrailing_lepton_vec)

# Calculate Quantitites
tDelta_phi_array = np.abs(tMuon_vec.deltaphi(tElectron_vec))
tDelta_eta_array = np.abs(tMuon_vec.deltaeta(tElectron_vec))
tDelta_r_array = np.abs(tMuon_vec.deltaR(tElectron_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
at.create_hist(tDilepton_vec.m, 'Truth DiLepton Mass', bins=50, range=(0,200000),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tDilepton_vec.pt, 'Truth Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tMuon_vec.pt, 'Truth Muon pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tMuon_vec.pt, cross_section))
at.create_hist(tElectron_vec.pt, 'Truth Electron pT', yscale='log', bins=50, range=(0,150000),
            weights=at.calculate_weights(tElectron_vec.pt, cross_section))
at.create_hist(tMuon_vec.eta, 'Truth Muon Eta', bins=50,
            weights=at.calculate_weights(tMuon_vec, cross_section))
at.create_hist(tElectron_vec.eta, 'Truth Electron Eta', bins=50,
            weights=at.calculate_weights(tElectron_vec, cross_section))
at.create_hist(tDelta_phi_array, 'Truth Delta Phi', bins=50,
            weights=at.calculate_weights(tDelta_phi_array, cross_section))
at.create_hist(tDelta_r_array, 'Truth Delta R', bins=50,
            weights=at.calculate_weights(tDelta_r_array, cross_section))
at.create_hist(tDelta_eta_array, 'Truth Delta Eta', bins=50,
            weights=at.calculate_weights(tDelta_eta_array, cross_section))