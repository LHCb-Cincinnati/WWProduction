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
jets = args.jets

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
both_lepton_loose_acc_mask = ((lminus_vec.eta>1.596)
                                & (lplus_vec.eta>1.596))
both_lepton_tight_acc_mask = ((lminus_vec.eta>2)
                                & (lminus_vec.eta<5)
                                & (lplus_vec.eta>2)
                                & (lplus_vec.eta<5)) 
high_pT_lepton_mask = ((lminus_vec.pt>15) | (lplus_vec.pt>15))
low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_tight_acc_mask = (((lminus_vec.eta>2)
                            & (lminus_vec.eta<5))
                            | ((lplus_vec.eta>2)
                            & (lplus_vec.eta<5)))
one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
invariant_mass_mask = (dilepton_vec.m>10)
lepton_mask = both_lepton_tight_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask

# Apply Masks
print("One Lepton in Acc", sum(one_lepton_tight_acc_mask))
print("One high pT lepton", sum(high_pT_lepton_mask))
print("Both Leptons in Acc", sum(both_lepton_tight_acc_mask))
print("Muon Electron Decay Mode", sum(mue_decay_mask))
print("LHCb Online Selection", sum(one_lepton_tight_acc_mask&high_pT_lepton_mask))
print("LHCb Offline selection", sum(both_lepton_tight_acc_mask&high_pT_lepton_mask&mue_decay_mask))
print("Muon Electron Decay Mode and both leptons in Acc", sum(mue_decay_mask&both_lepton_tight_acc_mask))
print("LHCb Offline selection", sum(both_lepton_tight_acc_mask&high_pT_lepton_mask&mue_decay_mask&low_pT_lepton_mask))
print("LHCb Offline selection", sum(both_lepton_tight_acc_mask&high_pT_lepton_mask&mue_decay_mask&low_pT_lepton_mask&(muon_vec.pt>15)))
pdb.set_trace()
muon_vec = muon_vec[lepton_mask]
electron_vec = electron_vec[lepton_mask]
dilepton_vec = dilepton_vec[lepton_mask]

# Jets
if jets:
    jet_vec = vector.zip({
        'px': tree['TargetJet'].px,
        'py': tree['TargetJet'].py,
        'pz': tree['TargetJet'].pz,
        'e': tree['TargetJet'].energy,
        'pid': tree['TargetJet'].id
    })

    antijet_vec = vector.zip({
        'px': tree['TargetAntiJet'].px,
        'py': tree['TargetAntiJet'].py,
        'pz': tree['TargetAntiJet'].pz,
        'e': tree['TargetAntiJet'].energy,
        'pid': tree['TargetAntiJet'].id
    })

# Calculate Quantitites
delta_phi_array = np.abs(muon_vec.deltaphi(electron_vec))
delta_eta_array = np.abs(muon_vec.deltaeta(electron_vec))
delta_r_array = np.abs(muon_vec.deltaR(electron_vec))

pdb.set_trace()
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

# at.create_hist(dilepton_vec.m, 'Standalone DiLepton Mass', bins=50, range=(0,500),
#             weights=at.calculate_weights(dilepton_vec, cross_section))
# at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Log', yscale='log', bins=50, range=(0,500),
#             weights=at.calculate_weights(dilepton_vec, cross_section))
# at.create_hist(dilepton_vec.pt, 'Standalone DiLepton pT Linear', bins=50, range=(0,500),
#             weights=at.calculate_weights(dilepton_vec, cross_section))
# at.create_hist(muon_vec.pt, 'Standalone Muon pT Log', yscale='log', bins=50, range=(0,150),
#             weights=at.calculate_weights(muon_vec.pt, cross_section))
# at.create_hist(electron_vec.pt, 'Standalone Electron pT Log', yscale='log', bins=50, range=(0,150),
#             weights=at.calculate_weights(electron_vec.pt, cross_section))
# at.create_hist(muon_vec.pt, 'Standalone Muon pT Linear', bins=50, range=(0,150),
#             weights=at.calculate_weights(muon_vec.pt, cross_section))
# at.create_hist(electron_vec.pt, 'Standalone Electron pT Linear', bins=50, range=(0,150),
#             weights=at.calculate_weights(electron_vec.pt, cross_section))
# at.create_hist(muon_vec.eta, 'Standalone Muon Eta', bins=50, range=(1.5, 5),
#             weights=at.calculate_weights(muon_vec, cross_section))
# at.create_hist(electron_vec.eta, 'Standalone Electron Eta', bins=50, range=(1.5, 5),
#             weights=at.calculate_weights(electron_vec, cross_section))
# at.create_hist(delta_phi_array, 'Standalone Delta Phi', bins=50,
#             weights=at.calculate_weights(delta_phi_array, cross_section))
# at.create_hist(delta_r_array, 'Standalone Delta R', bins=50, range=(0, 5),
#             weights=at.calculate_weights(delta_r_array, cross_section))
# at.create_hist(delta_eta_array, 'Standalone Delta Eta', bins=50, range=(0, 3),
#             weights=at.calculate_weights(delta_eta_array, cross_section))