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
lepton_vec = vector.zip({
    'px': tree['TargetLepton'].px,
    'py': tree['TargetLepton'].py,
    'pz': tree['TargetLepton'].pz,
    'e': tree['TargetLepton'].energy,
    'pid': tree['TargetLepton'].id
})

jet_vec = vector.zip({
    'px': tree['TargetJet'].px,
    'py': tree['TargetJet'].py,
    'pz': tree['TargetJet'].pz,
    'e': tree['TargetJet'].energy,
    'pid': tree['TargetJet'].id
})

# Masks
both_product_loose_acc_mask = ((lepton_vec.eta>2)
                                & (lepton_vec.eta<5)
                                & (jet_vec.eta>2)
                                & (jet_vec.eta<5)) 
both_product_tight_acc_mask = ((lepton_vec.eta>2.2)
                                & (lepton_vec.eta<4.4)
                                & (jet_vec.eta>2.2)
                                & (jet_vec.eta<4.4)) 
high_pT_product_mask = ((lepton_vec.pt>20) & (jet_vec.pt>20))
low_pT_product_mask = ((lepton_vec.pt>5) & (jet_vec.pt>5))
one_product_gauss_mask = ((lepton_vec.eta>1.596) & (lepton_vec.pt>17))
product_mask = both_product_tight_acc_mask&high_pT_product_mask&low_pT_product_mask

# Print
print(f"Total number of events: {len(lepton_vec)}")
print(f"Gauss Cuts: {sum((one_product_gauss_mask))}")
print(f"DaVinci Cut: {sum((low_pT_product_mask))}")
print(f"Total Cuts: {sum((one_product_gauss_mask&low_pT_product_mask))}")
print(f"Loose Acceptance Cuts: {sum(both_product_loose_acc_mask)}")
print(f"Tight Acceptance Cuts: {sum(product_mask)}")