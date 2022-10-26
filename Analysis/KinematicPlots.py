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

# Create Keys
truth_leading_lepton_key = 'tLeading_lepton_array'
truth_trailing_lepton_key = 'tTrailing_lepton_array'
reco_leading_lepton_key = 'rLeading_lepton_array'
reco_trailing_lepton_key = 'rTrailing_lepton_array'

# Create Vectors
tLeading_lepton_vec = vector.awk(tree[truth_leading_lepton_key].array())
tTrailing_lepton_vec = vector.awk(tree[truth_trailing_lepton_key].array())
tDilepton_vec = tLeading_lepton_vec + tTrailing_lepton_vec
rLeading_lepton_vec = vector.awk(tree[reco_leading_lepton_key].array())
rTrailing_lepton_vec = vector.awk(tree[reco_trailing_lepton_key].array())
rDilepton_vec = rLeading_lepton_vec + rTrailing_lepton_vec

# Masks
tEvt_num = tree['tEvt_num'].array().to_list()
rEvt_num = tree['rEvt_num'].array().to_list()
truth_lepton_acceptance_mask = ((tLeading_lepton_vec.eta>1.594)
                                & (tTrailing_lepton_vec.eta<1.594))
process_mask = (tree['tDecay_process_array'].array()['mue'] == 1)
good_reco_event_mask = [process_mask[tEvt_num.index(evt)] for evt in rEvt_num]
full_truth_mask = process_mask & truth_lepton_acceptance_mask
high_pT_lepton_mask = ((tTrailing_lepton_vec.pt>15) | (tTrailing_lepton_vec.pt>15))
one_lepton_acc_mask = ((tLeading_lepton_vec.eta>1.594) 
                        | (tTrailing_lepton_vec.eta>1.594))


# Mask Vectors and Extract Lepton Info
# Truth Vectors
tLeading_lepton_vec = tLeading_lepton_vec
tTrailing_lepton_vec = tTrailing_lepton_vec
tDilepton_vec = tDilepton_vec
# tMuon_vec = ak.where((abs(tLeading_lepton_vec.pid)==13), tLeading_lepton_vec, tTrailing_lepton_vec)
# tElectron_vec = ak.where((abs(tLeading_lepton_vec.pid)==11), tLeading_lepton_vec, tTrailing_lepton_vec)
# Reco Vectors
# rLeading_lepton_vec = rLeading_lepton_vec[good_reco_event_mask]
# rTrailing_lepton_vec = rTrailing_lepton_vec[good_reco_event_mask]
# rDilepton_vec = rDilepton_vec[good_reco_event_mask]
# rMuon_vec = ak.where((abs(rLeading_lepton_vec.pid)==13), rLeading_lepton_vec, rTrailing_lepton_vec)
# rElectron_vec = ak.where((abs(rLeading_lepton_vec.pid)==11), rLeading_lepton_vec, rTrailing_lepton_vec)

# Calculate Quantitites
# # Truth Quantities
tDelta_phi_array = np.abs(tLeading_lepton_vec.deltaphi(tTrailing_lepton_vec))
tDelta_eta_array = np.abs(tLeading_lepton_vec.deltaeta(tTrailing_lepton_vec))
tDelta_r_array = np.abs(tLeading_lepton_vec.deltaR(tTrailing_lepton_vec))
# # Reco Quantities
# rDelta_phi_array = np.abs(rMuon_vec.deltaphi(rElectron_vec))
# rDelta_eta_array = np.abs(rMuon_vec.deltaeta(rElectron_vec))
# rDelta_r_array = np.abs(rMuon_vec.deltaR(rElectron_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Testing

# # Plots
# Testing
at.create_hist(tLeading_lepton_vec.pt / 1000, 'Truth Leading Lepton Pt', bins=100, range=(0,200),
                density=True)
at.create_hist(tTrailing_lepton_vec.pt / 1000, 'Truth Trailing Lepton Pt', bins=100, range=(0,200),
                density=True)
at.create_hist(tLeading_lepton_vec.eta, 'Truth Leading Lepton Eta', bins=50, range=(-7,7),
                density=True)
at.create_hist(tTrailing_lepton_vec.eta, 'Truth Trailing Lepton Eta', bins=50, range=(-7,7),
                density=True)
at.create_hist(tDilepton_vec.m / 1000, 'Truth DiLepton Mass', bins=100, range=(0,500),
                density=True)
at.create_hist(tDilepton_vec.pt / 1000, 'Truth DiLepton pT', bins=100, range=(0,200),
                density=True)
at.create_hist(tDelta_phi_array, 'Truth Delta Phi', bins=50,
                density=True)
at.create_hist(tDelta_r_array, 'Truth Delta R', bins=50, range=(0, 7),
                density=True)
at.create_hist(tDelta_eta_array, 'Truth Delta Eta', bins=50, range=(0, 3),
                density=True)

# # Truth Plots
# at.create_hist(tDilepton_vec.m, 'Truth DiLepton Mass', bins=50, range=(0,200000),
#                weights=at.calculate_weights(tDilepton_vec, cross_section))
# at.create_hist(tDilepton_vec.m, 'Truth DiLepton Mass', bins=50, range=(0,200000),
#             weights=at.calculate_weights(tDilepton_vec, cross_section))
# at.create_hist(tDilepton_vec.pt, 'Truth Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(tDilepton_vec, cross_section))
# at.create_hist(tMuon_vec.pt, 'Truth Muon pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(tMuon_vec.pt, cross_section))
# at.create_hist(tElectron_vec.pt, 'Truth Electron pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(tElectron_vec.pt, cross_section))
# at.create_hist(tMuon_vec.eta, 'Truth Muon Eta', bins=50, range=(2, 5),
#             weights=at.calculate_weights(tMuon_vec, cross_section))
# at.create_hist(tElectron_vec.eta, 'Truth Electron Eta', bins=50, range=(2, 5),
#             weights=at.calculate_weights(tElectron_vec, cross_section))
# at.create_hist(tDelta_phi_array, 'Truth Delta Phi', bins=50,
#             weights=at.calculate_weights(tDelta_phi_array, cross_section))
# at.create_hist(tDelta_r_array, 'Truth Delta R', bins=50, range=(0, 4),
#             weights=at.calculate_weights(tDelta_r_array, cross_section))
# at.create_hist(tDelta_eta_array, 'Truth Delta Eta', bins=50, range=(0, 3),
#             weights=at.calculate_weights(tDelta_eta_array, cross_section))
# # Reco Plots
# at.create_hist(rDilepton_vec.m, 'Reco DiLepton Mass', bins=50, range=(0,200000),
#             weights=at.calculate_weights(rDilepton_vec, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rDilepton_vec.pt, 'Reco Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(rDilepton_vec, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rMuon_vec.pt, 'Reco Muon pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(rMuon_vec.pt, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rElectron_vec.pt, 'Reco Electron pT', yscale='log', bins=50, range=(0,150000),
#             weights=at.calculate_weights(rElectron_vec.pt, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rMuon_vec.eta, 'Reco Muon Eta', bins=50, range=(2, 5),
#             weights=at.calculate_weights(rMuon_vec, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rElectron_vec.eta, 'Reco Electron Eta', bins=50, range=(2, 5),
#             weights=at.calculate_weights(rElectron_vec, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rDelta_phi_array, 'Reco Delta Phi', bins=50,
#             weights=at.calculate_weights(rDelta_phi_array, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rDelta_r_array, 'Reco Delta R', bins=50, range=(0, 4),
#             weights=at.calculate_weights(rDelta_r_array, cross_section, nsim=len(tDelta_phi_array)))
# at.create_hist(rDelta_eta_array, 'Reco Delta Eta', bins=50, range=(0, 3),
#             weights=at.calculate_weights(rDelta_eta_array, cross_section, nsim=len(tDelta_phi_array)))