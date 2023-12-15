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
both_lepton_loose_acc_mask = ((lminus_vec.eta>2)
                                & (lminus_vec.eta<5)
                                & (lplus_vec.eta>2)
                                & (lplus_vec.eta<5)) 
both_lepton_tight_acc_mask = ((lminus_vec.eta>2.2)
                                & (lminus_vec.eta<4.4)
                                & (lplus_vec.eta>2.2)
                                & (lplus_vec.eta<4.4)) 
high_pT_lepton_mask = ((lminus_vec.pt>20) & (lplus_vec.pt>20))
low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
one_lepton_tight_acc_mask = ((lminus_vec.eta>2) | (lplus_vec.eta>2))
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt>15))
invariant_mass_mask = (dilepton_vec.m>10)
lepton_mask = both_lepton_tight_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask

# # Apply Masks
# muon_vec = muon_vec[lepton_mask&mue_decay_mask]
# electron_vec = electron_vec[lepton_mask&mue_decay_mask]
# dilepton_vec = dilepton_vec[lepton_mask&mue_decay_mask]

# # Testing Masks
# # lminus_vec = lminus_vec[one_lepton_loose_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask]
# # lplus_vec = lplus_vec[one_lepton_loose_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask]
# lplus_tight_acc_mask = ((lplus_vec.eta>2) & (lplus_vec.eta<5))
# lminus_tight_acc_mask = ((lminus_vec.eta>2) & (lminus_vec.eta<5))
# #jet_tight_acc_mask = ((jet_vec.eta>2) & (jet_vec.eta<5))
# # antijet_tight_acc_mask = ((antijet_vec.eta>2) & (antijet_vec.eta<5))
# one_lepton_tight_acc_mask = (((lminus_vec.eta>2)
#                         & (lminus_vec.eta<5))
#                         | ((lplus_vec.eta>2)
#                         & (lplus_vec.eta<5)))
# mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
#                     | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
# both_lepton_loose_acc_mask = ((lminus_vec.eta>1.596)
#                                 & (lplus_vec.eta>1.596))
# one_jet_one_lepton_mask = ((lminus_tight_acc_mask&antijet_tight_acc_mask) | (lplus_tight_acc_mask&jet_tight_acc_mask))
# two_jet_two_lepton_mask = ((lminus_tight_acc_mask&antijet_tight_acc_mask) & (lplus_tight_acc_mask&jet_tight_acc_mask))

print(f"Total number of events: {len(lminus_vec)}")
print(f"Gauss Cuts: {sum((one_lepton_gauss_mask))}")
print(f"DaVinci Cut: {sum((low_pT_lepton_mask))}")
print(f"Total Cuts: {sum((one_lepton_gauss_mask&low_pT_lepton_mask))}")
print(f"Loose Acceptance Cuts: {sum(both_lepton_loose_acc_mask&mue_decay_mask)}")
print(f"Tight Acceptance Cuts: {sum(lepton_mask)}")
pdb.set_trace()

# Calculate Quantitites
delta_phi_array = np.abs(muon_vec.deltaphi(electron_vec))
delta_eta_array = np.abs(muon_vec.deltaeta(electron_vec))
delta_r_array = np.abs(muon_vec.deltaR(electron_vec))
mt2_array = 2 * muon_vec.pt * electron_vec.pt * (1 + np.cos(muon_vec.deltaphi(electron_vec)))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
# 2D eta plots
# fig, axs = plt.subplots()
# hist, x_edges, y_edges = np.histogram2d(leading_lepton_vec.eta.to_numpy(),
#                                         jet_vec.eta.to_numpy(),
#                                         bins=(14,14),
#                                         range = ((-7,7),(-7,7)))
# X, Y = np.meshgrid(x_edges, y_edges)
# color_mesh = axs.pcolormesh(X, Y, hist)
# fig.colorbar(color_mesh, ax=axs)
# plt.title('ttbar Leading Lepton Eta vs Jet Eta')
# plt.xlabel('Leading Lepton Eta')
# plt.ylabel('Jet Eta')
# fig.savefig('ttbar_LeadingEtaVsJetEta' + '.png')

at.create_hist(mt2_array, 'Standalone MT2', bins=50, range=(0,5000),
            weights=at.calculate_weights(dilepton_vec, cross_section))

at.create_hist(dilepton_vec.m, 'Standalone DiLepton Mass', bins=50, range=(0,500),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Log', yscale='log', bins=50, range=(0,150),
            weights=at.calculate_weights(dilepton_vec, cross_section))
at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Linear', bins=50, range=(0,150),
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