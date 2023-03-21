# Imports
import sys
import os
import argparse 
import pdb

import numpy as np
import matplotlib as mpl
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

rwmass_eta_cut = ((rLeading_lepton_vec[rgood_muon_match_mask].eta>2.2)&(rLeading_lepton_vec[rgood_muon_match_mask].eta<4.4))
rwmass_pt_cut = ((rLeading_lepton_vec[rgood_muon_match_mask].pt>28000)&(rLeading_lepton_vec[rgood_muon_match_mask].eta<52000))
rwmass_mask = rwmass_eta_cut&rwmass_pt_cut
#rwmass_mask = [True]*len(rwmass_eta_cut)
twmass_eta_cut = ((tLeading_lepton_vec.eta>2.2)&(tLeading_lepton_vec.eta<4.4))
twmass_pt_cut = ((tLeading_lepton_vec.pt>28000))
twmass_mask = twmass_eta_cut&twmass_pt_cut
#twmass_mask = [True]*len(tLeading_lepton_vec.pt)


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

# Truth Hist no cuts
truth_hist_no_cuts, tx_edges, ty_edges = np.histogram2d(
                            tLeading_lepton_vec.rho.to_numpy() / 1000,
                            tLeading_lepton_vec.eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))
# Truth Hist with Cuts
truth_hist_wmass_cuts, tx_edges, ty_edges = np.histogram2d(
                            tLeading_lepton_vec[twmass_mask].rho.to_numpy() / 1000,
                            tLeading_lepton_vec[twmass_mask].eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))

# Reco Hist no Cuts
reco_hist_no_cuts, rx_edges, ry_edges = np.histogram2d(
                            truth_matched_vec.rho.to_numpy() / 1000,
                            truth_matched_vec.eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))
# Reco Hist with Cuts
reco_hist_wmass_cuts, rx_edges, ry_edges = np.histogram2d(
                            truth_matched_vec[rwmass_mask].rho.to_numpy() / 1000,
                            truth_matched_vec[rwmass_mask].eta.to_numpy(),
                            bins=(8, 5),
                            range=((0, 160), (1, 6)))


# Dividing Histograms
acc_hist = reco_hist_no_cuts/truth_hist_no_cuts
acc_hist = np.nan_to_num(acc_hist)
eff_hist = reco_hist_wmass_cuts/truth_hist_wmass_cuts
eff_hist = np.nan_to_num(eff_hist)


# Plots
#plt.rcParams.update({'font.size': 16})
# Acceptance Plots
at.create_imposed_hist([tLeading_lepton_vec.eta, truth_matched_vec.eta], "Muon Eta", bins=50, range=(0, 7), label=['MC', 'Reco'], color=['red', 'green'])
at.create_imposed_hist([tLeading_lepton_vec.pt / 1000, truth_matched_vec.rho / 1000], "Muon PT", bins=50, range=(0, 150), label=['MC', 'Reco'], color=['red', 'green'])
at.create_imposed_hist([tLeading_lepton_vec.eta[twmass_mask], truth_matched_vec.eta[rwmass_mask]], "Muon Eta (W Mass Cuts)", bins=50, range=(0, 7), label=['MC', 'Reco'], color=['red', 'green'])
at.create_imposed_hist([tLeading_lepton_vec.pt[twmass_mask] / 1000, truth_matched_vec.rho[rwmass_mask] / 1000], "Muon PT (W Mass Cuts)", bins=50, range=(0, 150), label=['MC', 'Reco'], color=['red', 'green'])

# Truth Muon Plot no Cuts
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, truth_hist_no_cuts.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Muon Histogram (No Cuts)')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMuonHistogramNoCuts' + '.png')
# Truth Muon Plot with Cuts
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, truth_hist_wmass_cuts.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Muon Histogram (Wmass Cuts)')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMuonHistogramWmassCuts' + '.png')
# Reco Plot no Cuts
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, reco_hist_wmass_cuts.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Matched Muon Histogram (No Cuts)')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMatchedMuonHistogramNoCuts' + '.png')
# Reco Plot with Cuts
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
color_mesh = axs.pcolormesh(X, Y, reco_hist_wmass_cuts.transpose(), shading='flat')
fig.colorbar(color_mesh, ax=axs)
plt.title('Truth Matched Muon Histogram (Wmass Cuts)')
plt.xlabel('Muon pT (GeV)')
plt.ylabel('Muon pseudorapidity ')
fig.savefig('TruthMatchedMuonHistogramWmassCuts' + '.png')

# Exam Plots
fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
divnorm = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=0.80, vmax=1)
colors_below90 = plt.cm.coolwarm(np.linspace(0, 0.4, 256))
colors_above90 = plt.cm.coolwarm(np.linspace(0.6, 1, 256))
all_colors = np.vstack((colors_below90, colors_above90))
efficiency_map = mpl.colors.LinearSegmentedColormap.from_list(
    'efficiency_map', all_colors)
color_mesh = axs.pcolormesh(X, Y, acc_hist.transpose(), shading='flat', norm=divnorm, cmap=efficiency_map)
fig.colorbar(color_mesh, ax=axs, ticks=[0,0.2,0.4,0.6,0.8,0.85, 0.9, 0.95, 1])
plt.title('Muon ID Acceptance Plot')
plt.xlabel('Muon $p_T$ (GeV)')
plt.ylabel('Muon $\\eta$')
fig.savefig('MuonIDAcceptancePlot' + '.png', dpi=600)

fig, axs = plt.subplots()
X, Y = np.meshgrid(tx_edges, ty_edges)
divnorm = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=0.80, vmax=1)
colors_below90 = plt.cm.coolwarm(np.linspace(0, 0.4, 256))
colors_above90 = plt.cm.coolwarm(np.linspace(0.6, 1, 256))
all_colors = np.vstack((colors_below90, colors_above90))
efficiency_map = mpl.colors.LinearSegmentedColormap.from_list(
    'efficiency_map', all_colors)
color_mesh = axs.pcolormesh(X, Y, eff_hist.transpose(), shading='flat', norm=divnorm, cmap=efficiency_map)
fig.colorbar(color_mesh, ax=axs, ticks=[0,0.2,0.4,0.6,0.8,0.85, 0.9, 0.95, 1])
plt.title('Muon ID Efficiency Plot')
plt.xlabel('Muon $p_T$ (GeV)')
plt.ylabel('Muon $\\eta$')
fig.savefig('MuonIDEfficiencyPlot' + '.png', dpi=600)