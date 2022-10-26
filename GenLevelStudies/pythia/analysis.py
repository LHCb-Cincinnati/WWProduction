# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector
sys.path.append('../../Analysis') # This is Awesome!
import analysis_tools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['AnalysisTree'].arrays()

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['LeptonMinus'].px,
    'py': tree['LeptonMinus'].py,
    'pz': tree['LeptonMinus'].pz,
    'e': tree['LeptonMinus'].energy,
    'pid': tree['LeptonMinus'].id
})

lplus_vec = vector.zip({
    'px': tree['LeptonPlus'].px,
    'py': tree['LeptonPlus'].py,
    'pz': tree['LeptonPlus'].pz,
    'e': tree['LeptonPlus'].energy,
    'pid': tree['LeptonPlus'].id
})
dilepton_vec = lminus_vec + lplus_vec
muon_vec = ak.where((abs(lminus_vec.pid)==13), lminus_vec, lplus_vec)
electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

# Masks
both_lepton_acc_mask = ((lminus_vec.eta>1.596) & (lplus_vec.eta>1.596))
high_pT_lepton_mask = ((lminus_vec.pt>15) | (lplus_vec.pt>15))
low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
mumu_decay_mask = ((tree['LeptonMinus'].id==13) & (tree['LeptonPlus'].id==-13))
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
invariant_mass_mask = (dilepton_vec.m>10)
lepton_mask = both_lepton_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mumu_decay_mask
drellyan_mask = both_lepton_acc_mask&low_pT_lepton_mask&invariant_mass_mask

# Apply Masks
leading_lepton_vec = leading_lepton_vec[one_lepton_acc_mask&high_pT_lepton_mask]
trailing_lepton_vec = trailing_lepton_vec[one_lepton_acc_mask&high_pT_lepton_mask]

# Calculate Quantitites
delta_phi_array = np.abs(leading_lepton_vec.deltaphi(trailing_lepton_vec))
delta_eta_array = np.abs(leading_lepton_vec.deltaeta(trailing_lepton_vec))
delta_r_array = np.abs(leading_lepton_vec.deltaR(trailing_lepton_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

pdb.set_trace()
# Plots
at.create_hist(leading_lepton_vec.pt, 'StandAlone Leading Lepton Pt', bins=100, range=(0,200),
                density=True)
at.create_hist(trailing_lepton_vec.pt, 'StandAlone Trailing Lepton Pt', bins=100, range=(0,200),
                density=True)
at.create_hist(leading_lepton_vec.eta, 'StandAlone Leading Lepton Eta', bins=50, range=(-7,7),
                density=True)
at.create_hist(trailing_lepton_vec.eta, 'StandAlone Trailing Lepton Eta', bins=50, range=(-7,7),
                density=True)
at.create_hist(dilepton_vec.m, 'StandAlone DiLepton Mass', bins=100, range=(0,500),
                density=True)
at.create_hist(dilepton_vec.pt, 'StandAlone DiLepton pT', bins=100, range=(0,200),
                density=True)
at.create_hist(delta_phi_array, 'StandAlone Delta Phi', bins=50,
                density=True)
at.create_hist(delta_r_array, 'StandAlone Delta R', bins=50, range=(0, 4),
                density=True)
at.create_hist(delta_eta_array, 'StandAlone Delta Eta', bins=50, range=(0, 3),
                density=True)