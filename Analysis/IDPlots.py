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
reco_leading_lepton_key = 'rLeading_lepton_id_array'
reco_trailing_lepton_key = 'rTrailing_lepton_id_array'

# Create Vectors
rLeading_lepton_id_array = tree[reco_leading_lepton_key].array()
rTrailing_lepton_id_array = tree[reco_trailing_lepton_key].array()

# Masks


# Mask Vectors and Extract Lepton Info
# Truth Vectors
# tLeading_lepton_vec = tLeading_lepton_vec[full_truth_mask]
# tTrailing_lepton_vec = tTrailing_lepton_vec[full_truth_mask]
# tDilepton_vec = tDilepton_vec[full_truth_mask]
# tMuon_vec = ak.where((abs(tLeading_lepton_vec.pid)==13), tLeading_lepton_vec, tTrailing_lepton_vec)
# tElectron_vec = ak.where((abs(tLeading_lepton_vec.pid)==11), tLeading_lepton_vec, tTrailing_lepton_vec)
# # Reco Vectors
# # rLeading_lepton_vec = rLeading_lepton_vec[good_reco_event_mask]
# # rTrailing_lepton_vec = rTrailing_lepton_vec[good_reco_event_mask]
# # rDilepton_vec = rDilepton_vec[good_reco_event_mask]
# # rMuon_vec = ak.where((abs(rLeading_lepton_vec.pid)==13), rLeading_lepton_vec, rTrailing_lepton_vec)
# # rElectron_vec = ak.where((abs(rLeading_lepton_vec.pid)==11), rLeading_lepton_vec, rTrailing_lepton_vec)

# Calculate Quantitites

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Testing

# Plots
# Reco Plots
at.create_hist(rLeading_lepton_id_array.ip, 'Leading Lepton IP', bins=10, range=(0,0.05))
