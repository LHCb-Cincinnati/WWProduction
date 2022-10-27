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
from hepunits.units import keV, MeV, GeV

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
truth_lepton_acceptance_mask = ((tLeading_lepton_vec.eta>1.596)
                                & (tTrailing_lepton_vec.eta>1.596))
process_mask = (tree['tDecay_process_array'].array()['mue'] == 1)
good_reco_event_mask = [process_mask[tEvt_num.index(evt)] for evt in rEvt_num]
high_pT_lepton_mask = ((tTrailing_lepton_vec.pt>15) | (tTrailing_lepton_vec.pt>15))
one_lepton_acc_mask = ((tLeading_lepton_vec.eta>1.596) 
                        | (tTrailing_lepton_vec.eta>1.596))
full_truth_mask = process_mask & truth_lepton_acceptance_mask & high_pT_lepton_mask


# Mask Vectors and Extract Lepton Info
# Truth Vectors
tLeading_lepton_vec = tLeading_lepton_vec[full_truth_mask]
tTrailing_lepton_vec = tTrailing_lepton_vec[full_truth_mask]
tDilepton_vec = tDilepton_vec[full_truth_mask]
tMuon_vec = ak.where((abs(tLeading_lepton_vec.pid)==13), tLeading_lepton_vec, tTrailing_lepton_vec)
tElectron_vec = ak.where((abs(tLeading_lepton_vec.pid)==11), tLeading_lepton_vec, tTrailing_lepton_vec)
# Reco Vectors
# rLeading_lepton_vec = rLeading_lepton_vec[good_reco_event_mask]
# rTrailing_lepton_vec = rTrailing_lepton_vec[good_reco_event_mask]
# rDilepton_vec = rDilepton_vec[good_reco_event_mask]
# rMuon_vec = ak.where((abs(rLeading_lepton_vec.pid)==13), rLeading_lepton_vec, rTrailing_lepton_vec)
# rElectron_vec = ak.where((abs(rLeading_lepton_vec.pid)==11), rLeading_lepton_vec, rTrailing_lepton_vec)

# Calculate Quantitites
# # Truth Quantities
tDelta_phi_array = np.abs(tMuon_vec.deltaphi(tElectron_vec))
tDelta_eta_array = np.abs(tMuon_vec.deltaeta(tElectron_vec))
tDelta_r_array = np.abs(tMuon_vec.deltaR(tElectron_vec))
# # Reco Quantities
# rDelta_phi_array = np.abs(rMuon_vec.deltaphi(rElectron_vec))
# rDelta_eta_array = np.abs(rMuon_vec.deltaeta(rElectron_vec))
# rDelta_r_array = np.abs(rMuon_vec.deltaR(rElectron_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Testing

# Plots
# Truth Plots
at.create_hist(tDilepton_vec.m / GeV, 'Truth DiLepton Mass', bins=50, range=(0,500),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tDilepton_vec.pt / GeV, 'Truth DiLepton pT Log', yscale='log', bins=50, range=(0,500),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tDilepton_vec.pt / GeV, 'Truth DiLepton pT Linear', bins=50, range=(0,500),
            weights=at.calculate_weights(tDilepton_vec, cross_section))
at.create_hist(tMuon_vec.pt / GeV, 'Truth Muon pT Log', yscale='log', bins=50, range=(0,150),
            weights=at.calculate_weights(tMuon_vec.pt, cross_section))
at.create_hist(tElectron_vec.pt / GeV, 'Truth Electron pT Log', yscale='log', bins=50, range=(0,150),
            weights=at.calculate_weights(tElectron_vec.pt, cross_section))
at.create_hist(tMuon_vec.pt / GeV, 'Truth Muon pT Linear', bins=50, range=(0,150),
            weights=at.calculate_weights(tMuon_vec.pt, cross_section))
at.create_hist(tElectron_vec.pt / GeV, 'Truth Electron pT Linear', bins=50, range=(0,150),
            weights=at.calculate_weights(tElectron_vec.pt, cross_section))
at.create_hist(tMuon_vec.eta, 'Truth Muon Eta', bins=50, range=(1.5, 7),
            weights=at.calculate_weights(tMuon_vec, cross_section))
at.create_hist(tElectron_vec.eta, 'Truth Electron Eta', bins=50, range=(1.5, 7),
            weights=at.calculate_weights(tElectron_vec, cross_section))
at.create_hist(tDelta_phi_array, 'Truth Delta Phi', bins=50,
            weights=at.calculate_weights(tDelta_phi_array, cross_section))
at.create_hist(tDelta_r_array, 'Truth Delta R', bins=50, range=(0, 5),
            weights=at.calculate_weights(tDelta_r_array, cross_section))
at.create_hist(tDelta_eta_array, 'Truth Delta Eta', bins=50, range=(0, 3),
            weights=at.calculate_weights(tDelta_eta_array, cross_section))
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