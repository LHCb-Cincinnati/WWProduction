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
reco_leading_lepton_id_key = 'rTrailing_lepton_id_array'
reco_trailing_lepton_id_key = 'rTrailing_lepton_id_array'

# Create Vectors/Arrays
# Truth Vectors
tLeading_lepton_vec = vector.awk(tree[truth_leading_lepton_key].array())
tTrailing_lepton_vec = vector.awk(tree[truth_trailing_lepton_key].array())
tDilepton_vec = tLeading_lepton_vec + tTrailing_lepton_vec
tMuon_vec = ak.where((abs(tLeading_lepton_vec.pid)==13),
                    tLeading_lepton_vec, tTrailing_lepton_vec)
tElectron_vec = ak.where((abs(tLeading_lepton_vec.pid)==11),
                        tLeading_lepton_vec, tTrailing_lepton_vec)
# Reconstructed Vectors
rLeading_lepton_vec = vector.awk(tree[reco_leading_lepton_key].array())
rTrailing_lepton_vec = vector.awk(tree[reco_trailing_lepton_key].array())
rDilepton_vec = rLeading_lepton_vec + rTrailing_lepton_vec
rLeading_lepton_id_array = tree[reco_leading_lepton_id_key].array()
rTrailing_lepton_id_array = tree[reco_trailing_lepton_id_key].array()
rMuon_vec = ak.where((abs(rLeading_lepton_vec.pid)==13),
                    rLeading_lepton_vec, rTrailing_lepton_vec)
rElectron_vec = ak.where((abs(rLeading_lepton_vec.pid)==11),
                        rLeading_lepton_vec, rTrailing_lepton_vec)
rMuon_id_array = ak.where((abs(rLeading_lepton_vec.pid)==13),
                    rLeading_lepton_id_array, rTrailing_lepton_id_array)
rElectron_id_array = ak.where((abs(rLeading_lepton_vec.pid)==11),
                    rLeading_lepton_id_array, rTrailing_lepton_id_array)

# Masks
tEvt_num = tree['tEvt_num'].array()
rEvt_num = tree['rEvt_num'].array()
truth_lepton_acceptance_mask = ((tLeading_lepton_vec.eta>2)
                                & (tLeading_lepton_vec.eta<5)
                                & (tTrailing_lepton_vec.eta>2)
                                & (tTrailing_lepton_vec.eta<5))
process_mask = (tree['tDecay_process_array'].array()['mue'] == 1)
good_reco_event_mask = [process_mask[np.where(tEvt_num==evt)[0][0]] for evt in rEvt_num]
full_truth_mask = process_mask & truth_lepton_acceptance_mask

# Mask Vectors
# Truth Vectors
tEvt_num = tEvt_num[full_truth_mask]
tLeading_lepton_vec = tLeading_lepton_vec[full_truth_mask]
tTrailing_lepton_vec = tTrailing_lepton_vec[full_truth_mask]
tDilepton_vec = tDilepton_vec[full_truth_mask]
tMuon_vec = tMuon_vec[full_truth_mask]
tElectron_vec = tElectron_vec[full_truth_mask]
# Reco Vectors
rEvt_num = rEvt_num[good_reco_event_mask]
rLeading_lepton_vec = rLeading_lepton_vec[good_reco_event_mask]
rTrailing_lepton_vec = rTrailing_lepton_vec[good_reco_event_mask]
rLeading_lepton_id_array = rLeading_lepton_id_array[good_reco_event_mask]
rTrailing_lepton_id_array = rTrailing_lepton_id_array[good_reco_event_mask]
rDilepton_vec = rDilepton_vec[good_reco_event_mask]
rMuon_vec = rMuon_vec[good_reco_event_mask]
rElectron_vec = rElectron_vec[good_reco_event_mask]
rMuon_id_array = rMuon_id_array[good_reco_event_mask]
rElectron_id_array = rElectron_id_array[good_reco_event_mask]

for i in range(len(rEvt_num)):
    if i==3:
        continue
    index = np.where(tEvt_num==rEvt_num[i])[0][0]
    at.DeltaRMatchBool(rMuon_vec[i], tMuon_vec[index])
pdb.set_trace()

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path('EfficiencyPlots.root', args.testing)
os.chdir(file_path)

# Testing
# Truth Hist
truth_hist, tx_edges, ty_edges = np.histogram2d(
                            tMuon_vec.pt.to_numpy(),
                            tMuon_vec.eta.to_numpy(),
                            bins=(10, 3),
                            range=((0, 150000), (2, 5)),
                            weights=at.calculate_weights(tDilepton_vec, 
                                                        cross_section))
# Loose Reco Hist
loose_reco_hist, rx_edges, ry_edges = np.histogram2d(
                            rMuon_vec.pt.to_numpy(),
                            rMuon_vec.eta.to_numpy(),
                            bins=(10, 3),
                            range=((0, 150000), (2, 5)),
                            weights=at.calculate_weights(rDilepton_vec, 
                                                        cross_section,
                                                        nsim=len(tDilepton_vec)))
# Tight Reco Hist
tight_reco_hist, rx_edges, ry_edges = np.histogram2d(
                            rTightmuon_vec.pt.to_numpy(),
                            rTightmuon_vec.eta.to_numpy(),
                            bins=(10, 3),
                            range=((0, 150000), (2, 5)),
                            weights=at.calculate_weights(rTightmuon_vec, 
                                                        cross_section,
                                                        nsim=len(tDilepton_vec)))


# Dividing Histograms
loose_eff_hist = loose_reco_hist/truth_hist
loose_eff_hist = np.nan_to_num(loose_eff_hist)
tight_eff_hist = loose_reco_hist/truth_hist
tight_eff_hist = np.nan_to_num(tight_eff_hist)


# Plots
# Truth Muon Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, truth_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Muon Histogram')
plt.xlabel('Muon pT (MeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMuonHistogram' + '.png')
# Loose Reco Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, loose_reco_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Reco Muon Histogram')
plt.xlabel('Muon pT (MeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('RecoMuonHistogram' + '.png')
# Loose Muon ID Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, loose_eff_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Loose Muon ID Efficiency Plot')
plt.xlabel('Muon pT (MeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('LooseMuonIDEfficiencyPlot' + '.png')
# Tight Muon ID Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, tight_eff_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Tight Muon ID Efficiency Plot')
plt.xlabel('Muon pT (MeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TightMuonIDEfficiencyPlot' + '.png')