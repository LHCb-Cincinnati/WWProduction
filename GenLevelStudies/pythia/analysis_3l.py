# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

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
    'px': tree['lminus'].px,
    'py': tree['lminus'].py,
    'pz': tree['lminus'].pz,
    'e': tree['lminus'].energy,
    'pid': tree['lminus'].id
})

lplus_vec = vector.zip({
    'px': tree['lplus'].px,
    'py': tree['lplus'].py,
    'pz': tree['lplus'].pz,
    'e': tree['lplus'].energy,
    'pid': tree['lplus'].id
})

thirdl_vec = vector.zip({
    'px': tree['LeptonThree'].px,
    'py': tree['LeptonThree'].py,
    'pz': tree['LeptonThree'].pz,
    'e': tree['LeptonThree'].energy,
    'pid': tree['LeptonThree'].id
})

dilepton_vec = lminus_vec + lplus_vec
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

# Masks
all_lepton_tight_acc_mask = ((lminus_vec.eta>2.2)
                                & (lminus_vec.eta<4.4)
                                & (lplus_vec.eta>2.2)
                                & (lplus_vec.eta<4.4)
                                & (thirdl_vec.eta>2.2)
                                & (thirdl_vec.eta<4.4)) 
high_pT_lepton_mask = ((lminus_vec.pt>20) & (lplus_vec.pt>20) & (thirdl_vec.pt>20))
low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5) & (thirdl_vec.pt>5))
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
one_lepton_tight_acc_mask = ((lminus_vec.eta>2) | (lplus_vec.eta>2))
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt>15)
                         | (thirdl_vec.eta>1.596) & (thirdl_vec.pt>15))
invariant_mass_mask = (dilepton_vec.m>10)
# lepton_mask = both_lepton_tight_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask

# Apply Masks
# dilepton_vec = dilepton_vec[one_lepton_gauss_mask&low_pT_lepton_mask]

print(f"Total number of events: {len(lminus_vec)}")
print(f"Gauss Cuts: {sum((one_lepton_gauss_mask))}")
print(f"DaVinci Cut: {sum((low_pT_lepton_mask))}")
print(f"LHCb Acceptance: {sum((all_lepton_tight_acc_mask)&(high_pT_lepton_mask))}")
print(f"Total Cuts: {sum((one_lepton_gauss_mask&low_pT_lepton_mask))}")