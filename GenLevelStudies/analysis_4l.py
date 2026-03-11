# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector
sys.path.append(".")
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
cross_section = args.cross_section # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tree'].arrays()

# Create Vectors
Z1_lminus_vec = vector.zip({
    'px': tree['Z1_lminus'].px,
    'py': tree['Z1_lminus'].py,
    'pz': tree['Z1_lminus'].pz,
    'e': tree['Z1_lminus'].e,
    'pid': tree['Z1_lminus'].pid
})

Z1_lplus_vec = vector.zip({
    'px': tree['Z1_lplus'].px,
    'py': tree['Z1_lplus'].py,
    'pz': tree['Z1_lplus'].pz,
    'e': tree['Z1_lplus'].e,
    'pid': tree['Z1_lplus'].pid
})
Z2_lminus_vec = vector.zip({
    'px': tree['Z2_lminus'].px,
    'py': tree['Z2_lminus'].py,
    'pz': tree['Z2_lminus'].pz,
    'e': tree['Z2_lminus'].e,
    'pid': tree['Z2_lminus'].pid
})

Z2_lplus_vec = vector.zip({
    'px': tree['Z2_lplus'].px,
    'py': tree['Z2_lplus'].py,
    'pz': tree['Z2_lplus'].pz,
    'e': tree['Z2_lplus'].e,
    'pid': tree['Z2_lplus'].pid
})



# Masks
# Z1_lminus and Z2_lplus pairs
Z1lmZ2lp_mue_decay_mask = (
    (abs(Z1_lminus_vec.pid)!=abs(Z2_lplus_vec.pid))
    & (Z1_lminus_vec.pid*Z2_lplus_vec.pid < 0)
)
Z1lmZ2lp_tight_acc_mask = (
    (Z1_lminus_vec.eta>2.2)
    & (Z1_lminus_vec.eta<4.4)
    & (Z2_lplus_vec.eta>2.2)
    & (Z2_lplus_vec.eta<4.4)
)
Z1lmZ2lp_high_pT_mask = (
    (Z1_lminus_vec.pt>20)
    & (Z2_lplus_vec.pt>20)
)
Z1lmZ2lp_masks = (
    Z1lmZ2lp_mue_decay_mask
    & Z1lmZ2lp_tight_acc_mask
    & Z1lmZ2lp_high_pT_mask
)
# Z2_lminus and Z1_lplus pairs
Z2lmZ1lp_mue_decay_mask = (
    (abs(Z2_lminus_vec.pid)!=abs(Z1_lplus_vec.pid))
    & (Z2_lminus_vec.pid*Z1_lplus_vec.pid < 0)
)
Z2lmZ1lp_tight_acc_mask = (
    (Z2_lminus_vec.eta>2.2)
    & (Z2_lminus_vec.eta<4.4)
    & (Z1_lplus_vec.eta>2.2)
    & (Z1_lplus_vec.eta<4.4)
)
Z2lmZ1lp_high_pT_mask = (
    (Z2_lminus_vec.pt>20)
    & (Z1_lplus_vec.pt>20)
)
Z2lmZ1lp_masks = (
    Z2lmZ1lp_mue_decay_mask
    & Z2lmZ1lp_tight_acc_mask
    & Z2lmZ1lp_high_pT_mask
)
# Apply Masks
print(f"Total number of events: {len(Z1_lminus_vec)}")
print(f"Mu-E Decay Mode: {sum((Z1lmZ2lp_mue_decay_mask) | (Z2lmZ1lp_mue_decay_mask))}")
print(f"Tight Eta Acc: {sum((Z1lmZ2lp_tight_acc_mask) | (Z2lmZ1lp_tight_acc_mask))}")
print(f"High pT: {sum((Z1lmZ2lp_high_pT_mask) | (Z2lmZ1lp_high_pT_mask))}")
print(f"LHCb Acceptance: {sum((Z1lmZ2lp_masks)|(Z2lmZ1lp_masks))}")
