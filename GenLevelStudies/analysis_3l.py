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
lminus_vec = vector.zip({
    'px': tree['Z_lminus'].px,
    'py': tree['Z_lminus'].py,
    'pz': tree['Z_lminus'].pz,
    'e': tree['Z_lminus'].e,
    'pid': tree['Z_lminus'].pid
})

lplus_vec = vector.zip({
    'px': tree['Z_lplus'].px,
    'py': tree['Z_lplus'].py,
    'pz': tree['Z_lplus'].pz,
    'e': tree['Z_lplus'].e,
    'pid': tree['Z_lplus'].pid
})

thirdl_vec = vector.zip({
    'px': tree['W_lepton'].px,
    'py': tree['W_lepton'].py,
    'pz': tree['W_lepton'].pz,
    'e': tree['W_lepton'].e,
    'pid': tree['W_lepton'].pid
})

dilepton_vec = lminus_vec + lplus_vec
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

# Masks
# general masks
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt>15)
                         | (thirdl_vec.eta>1.596) & (thirdl_vec.pt>15))
# lminus and wlepton pair masks
lmwl_mue_decay_mask = (
    (abs(lminus_vec.pid)!=abs(thirdl_vec.pid))
    & (lminus_vec.pid*thirdl_vec.pid < 0)
)
lmwl_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lmwl_high_pT_mask = (
    (lminus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lmwl_masks = (
    lmwl_mue_decay_mask
    & lmwl_tight_acc_mask
    & lmwl_high_pT_mask
)
# lplus and wlepton pair masks
lpwl_mue_decay_mask = (
    (abs(lplus_vec.pid)!=abs(thirdl_vec.pid))
    & (lplus_vec.pid*thirdl_vec.pid < 0)
)
lpwl_tight_acc_mask = (
    (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lpwl_high_pT_mask = (
    (lplus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lpwl_masks = (
    lpwl_mue_decay_mask
    & lpwl_tight_acc_mask
    & lpwl_high_pT_mask
)
# Apply Masks
print(f"Total number of events: {len(lminus_vec)}")
print(f"Gauss Cuts: {sum((one_lepton_gauss_mask))}")
print(f"Mu-E Decay Mode: {sum((lpwl_mue_decay_mask) | (lmwl_mue_decay_mask))}")
print(f"Tight Eta Acc: {sum((lpwl_tight_acc_mask) | (lmwl_tight_acc_mask))}")
print(f"High pT: {sum((lpwl_high_pT_mask) | (lmwl_high_pT_mask))}")
print(f"LHCb Acceptance: {sum((lmwl_masks)|(lpwl_masks))}")
