# Imports
import sys
import os
import pdb
import pickle
# Scikit Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages
import boost_histogram as bh
import awkward as ak
import uproot
import vector
from hepunits.units import MeV, GeV
# Personal Packages 
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
ofile_name = args.output
cross_section = args.cross_section # Cross section in fb
luminosity = args.luminosity # Luminosity in fb^-1

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tree'].arrays()

# Define histograms
muon_eta_deltar_hist = bh.Histogram(
    bh.axis.Regular(12, 2, 5),
    bh.axis.Regular(50, 0, 5),
    storage=bh.storage.Weight()
)
# leadinglepton_id_eta_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 2.2, 4.4), storage=bh.storage.Weight())

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['TargetLepton'].px * GeV,
    'py': tree['TargetLepton'].py * GeV,
    'pz': tree['TargetLepton'].pz * GeV,
    'e': tree['TargetLepton'].e * GeV,
    'pid': tree['TargetLepton'].pid
})

lplus_vec = vector.zip({
    'px': tree['TargetAntiLepton'].px * GeV,
    'py': tree['TargetAntiLepton'].py * GeV,
    'pz': tree['TargetAntiLepton'].pz * GeV,
    'e': tree['TargetAntiLepton'].e * GeV,
    'pid': tree['TargetAntiLepton'].pid
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
both_lepton_tight_acc_mask = ((np.abs(lminus_vec.eta)<6) & (np.abs(lplus_vec.eta)<6)) 
                                # & (lminus_vec.eta<4.4)
                                # & (lplus_vec.eta>2.2)
high_pT_lepton_mask = ((lminus_vec.pt / GeV >20) & (lplus_vec.pt / GeV >20))
low_pT_lepton_mask = ((lminus_vec.pt / GeV >5) & (lplus_vec.pt / GeV >5))
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
one_lepton_tight_acc_mask = ((lminus_vec.eta>2) | (lplus_vec.eta>2))
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt / GeV >15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt / GeV >15))
invariant_mass_mask = (dilepton_vec.m / GeV >100)
lepton_mask = (
    both_lepton_tight_acc_mask
    & high_pT_lepton_mask
    & low_pT_lepton_mask
    & mue_decay_mask
)

# Apply Masks
muon_vec = muon_vec[lepton_mask]
electron_vec = electron_vec[lepton_mask]
leading_lepton_vec = leading_lepton_vec[lepton_mask]
trailing_lepton_vec = trailing_lepton_vec[lepton_mask]
dilepton_vec = dilepton_vec[lepton_mask]

# Print Info
print(f"Total number of events: {len(lminus_vec)}")
print(f"Gauss Cuts: {sum((one_lepton_gauss_mask))}")
print(f"Mu-E Decay Mode: {sum((mue_decay_mask))}")
print(f"Gauss Cuts and Mu-E Decay Mode: {sum((one_lepton_gauss_mask&mue_decay_mask))}")
print(f"High pT Cuts: {sum((high_pT_lepton_mask))}")
print(f"Tight Acceptance Eta Cuts: {sum(both_lepton_tight_acc_mask)}")
print(f"High Mass Cut: {sum(invariant_mass_mask)}")
print(f"High Mass Cut and Tight Acceptance: {sum(invariant_mass_mask&lepton_mask)}")
print(f"Total Cuts: {sum(lepton_mask)}")
if args.debug:
    exit()

# Calculate Quantitites
delta_phi_array = np.abs(muon_vec.deltaphi(electron_vec))
delta_eta_array = np.abs(muon_vec.deltaeta(electron_vec))
delta_r_array = np.abs(muon_vec.deltaR(electron_vec))

# Fill histograms
# Identified Rough Binned
muon_eta_deltar_hist.fill(
    muon_vec.eta,
    delta_r_array
)

# Scale histograms
if not cross_section:
    scale_factor = 1
else:
    scale_factor = cross_section

# Scale histograms
muon_eta_deltar_hist *= scale_factor

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
# Eta Eta Yield 2D Hist
fig, axs = plt.subplots()
colormesh = plt.pcolormesh(
    muon_eta_deltar_hist.axes[0].edges,
    muon_eta_deltar_hist.axes[1].edges,
    muon_eta_deltar_hist.view().value.T
)
colorbar = plt.colorbar(colormesh)
colorbar.ax.set_ylabel("Yield", rotation=-90)
axs.set_xlabel("$\eta_{\mu}$")
axs.set_ylabel("$\Delta r_{e \mu}$")
save_str = ''.join("MuonEtaDeltarHist".split())
plt.savefig(save_str + '.png')
plt.close()

# # Save histograms
# os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
# pickle_dict = {"LeadingLeptonpTRough": [leadinglepton_id_pT_rghbin_hist]}
# with open(ofile_name + ".pkl", "wb") as f:
#     pickle.dump(pickle_dict, f)
