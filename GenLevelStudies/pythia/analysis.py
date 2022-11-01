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
tree = ifile['Tree'].arrays()

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['TargetLepton'].px,
    'py': tree['TargetLepton'].py,
    'pz': tree['TargetLepton'].pz,
    'e': tree['TargetLepton'].energy,
    'pid': tree['TargetLepton'].id
})

lplus_vec = vector.zip({
    'px': tree['TargetAntiLepton'].px,
    'py': tree['TargetAntiLepton'].py,
    'pz': tree['TargetAntiLepton'].pz,
    'e': tree['TargetAntiLepton'].energy,
    'pid': tree['TargetAntiLepton'].id
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
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
invariant_mass_mask = (dilepton_vec.m>10)
lepton_mask = both_lepton_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask
drellyan_mask = both_lepton_acc_mask&low_pT_lepton_mask&invariant_mass_mask

# Apply Masks
muon_vec = muon_vec[lepton_mask]
electron_vec = electron_vec[lepton_mask]
dilepton_vec = dilepton_vec[lepton_mask]

# Calculate Quantitites
delta_phi_array = np.abs(muon_vec.deltaphi(electron_vec))
delta_eta_array = np.abs(muon_vec.deltaeta(electron_vec))
delta_r_array = np.abs(muon_vec.deltaR(electron_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
at.create_hist(dilepton_vec.m, 'Standalone DiLepton Mass', bins=50, range=(0,500),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Log', yscale='log', bins=50, range=(0,500),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Linear', bins=50, range=(0,500),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(muon_vec.pt, 'Standalone Muon pT Log', yscale='log', bins=50, range=(0,150),
            weights=at.calculate_weights(muon_vec.pt, cross_section))
at.create_hist(electron_vec.pt, 'Standalone Electron pT Log', yscale='log', bins=50, range=(0,150),
            weights=at.calculate_weights(electron_vec.pt, cross_section))
at.create_hist(muon_vec.pt, 'Standalone Muon pT Linear', bins=50, range=(0,150),
            weights=at.calculate_weights(muon_vec.pt, cross_section))
at.create_hist(electron_vec.pt, 'Standalone Electron pT Linear', bins=50, range=(0,150),
            weights=at.calculate_weights(electron_vec.pt, cross_section))
at.create_hist(muon_vec.eta, 'Standalone Muon Eta', bins=50, range=(1.5, 5),
            weights=at.calculate_weights(muon_vec, cross_section))
at.create_hist(electron_vec.eta, 'Standalone Electron Eta', bins=50, range=(1.5, 5),
            weights=at.calculate_weights(electron_vec, cross_section))
at.create_hist(delta_phi_array, 'Standalone Delta Phi', bins=50,
            weights=at.calculate_weights(delta_phi_array, cross_section))
at.create_hist(delta_r_array, 'Standalone Delta R', bins=50, range=(0, 5),
            weights=at.calculate_weights(delta_r_array, cross_section))
at.create_hist(delta_eta_array, 'Standalone Delta Eta', bins=50, range=(0, 3),
            weights=at.calculate_weights(delta_eta_array, cross_section))