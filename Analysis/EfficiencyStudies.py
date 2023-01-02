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
truth_leading_lepton_key = 'tMuon_array'
reco_leading_lepton_key = 'rMuon_array'
reco_leading_lepton_id_key = 'rMuon_id_array'
rLeading_lepton_deltaRmatch_key = 'rMuon_deltaRmatch'

# Create Vectors/Arrays
# Truth Vectors
tEvt_num = tree['tEvt_num'].array()
tLeading_lepton_vec = vector.awk(tree[truth_leading_lepton_key].array())
# Reconstructed Vectors
rEvt_num = tree['rEvt_num'].array()
rLeading_lepton_vec = vector.awk(tree[reco_leading_lepton_key].array())
rLeading_lepton_id_array = tree[reco_leading_lepton_id_key].array()
rLeading_lepton_deltaRmatch_array = tree[rLeading_lepton_deltaRmatch_key].array()

# Masks
truth_lepton_acceptance_mask = (tLeading_lepton_vec.eta>1.596)
rgood_muon_match_mask = (rLeading_lepton_deltaRmatch_array<0.1)
# rgood_muon_match_evts = [evt for index,evt in enumerate(rEvt_num) if rgood_muon_match_mask[index]]
# tevt_mask = [np.size(np.where(tEvt_num==evt))>0 for evt in rgood_muon_match_evts]

# Mask Vectors
# Truth Vectors
tEvt_num_full = tEvt_num
tEvt_num = tEvt_num
tLeading_lepton_vec = tLeading_lepton_vec
# Reco Vectors
rEvt_num = rEvt_num[rgood_muon_match_mask]
rLeading_lepton_vec = rLeading_lepton_vec[rgood_muon_match_mask]
rLeading_lepton_id_array = rLeading_lepton_id_array[rgood_muon_match_mask]
rLeading_lepton_deltaRmatch_array = rLeading_lepton_deltaRmatch_array[rgood_muon_match_mask]

# Calculate Quantities
# Debug for loop
# for evt in rEvt_num:
#     tindex = np.where(tEvt_num==evt)[0][0]
#     truth_match = tLeading_lepton_vec[tindex]

truth_matched_vec = vector.awk([tLeading_lepton_vec[np.where(tEvt_num==evt)[0][0]] for evt in rEvt_num])
# pdb.set_trace()

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path('EfficiencyPlots.root', args.testing)
os.chdir(file_path)

# Testing
# Truth Hist
truth_hist, tx_edges, ty_edges = np.histogram2d(
                            tLeading_lepton_vec.rho.to_numpy() / 1000,
                            tLeading_lepton_vec.eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))
# Loose Reco Hist
loose_reco_hist, rx_edges, ry_edges = np.histogram2d(
                            truth_matched_vec.rho.to_numpy() / 1000,
                            truth_matched_vec.eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))
# Tight Reco Hist
tight_reco_hist, rx_edges, ry_edges = np.histogram2d(
                            truth_matched_vec[(rLeading_lepton_id_array.isMuon==1)].rho.to_numpy() / 1000,
                            truth_matched_vec[(rLeading_lepton_id_array.isMuon==1)].eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))


# Dividing Histograms
loose_eff_hist = loose_reco_hist/truth_hist
loose_eff_hist = np.nan_to_num(loose_eff_hist)
tight_eff_hist = tight_reco_hist/truth_hist
tight_eff_hist = np.nan_to_num(tight_eff_hist)


# Plots
# Truth Muon Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, truth_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Muon Histogram')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMuonHistogram' + '.png')
# Loose Reco Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, loose_reco_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Matched Muon Histogram')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMatchedMuonHistogram' + '.png')
# Loose Muon ID Plot Zoomed in
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, loose_eff_hist.transpose(), shading='flat', vmin=0.80, vmax=1)
fig.colorbar(color_mesh, ax=axs)
plt.title('Loose Muon ID Efficiency Plot Zoomed In')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('LooseMuonIDEfficiencyPlotZoom' + '.png')
# Tight Muon ID Plot Zoomed in
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, tight_eff_hist.transpose(), shading='flat', vmin=0.80, vmax=1)
fig.colorbar(color_mesh, ax=axs)
plt.title('Tight Muon ID Efficiency Plot Zoomed In')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TightMuonIDEfficiencyPlotZoom' + '.png')
# Loose Muon ID Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, loose_eff_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Loose Muon ID Efficiency Plot')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('LooseMuonIDEfficiencyPlot' + '.png')
# Tight Muon ID Plot
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, tight_eff_hist.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Tight Muon ID Efficiency Plot')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TightMuonIDEfficiencyPlot' + '.png')