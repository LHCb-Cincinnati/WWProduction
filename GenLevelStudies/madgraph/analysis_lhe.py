# Imports
# STL Packages
import pdb
import os
import sys
import pickle
# Scikit Packages
import matplotlib.pyplot as plt
import numpy as np
# HEP Packages
import awkward as ak
import vector
import pylhe
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
ifile_name =  args.input_files[0].name
cross_section = args.cross_section # Cross section in fb

# pylhe SETUP
# # Ensures that pylhe reads in data as an akward array.
# pylhe.register_awkward() 
# Read in madgraph data.  
arr = pylhe.to_awkward(pylhe.read_lhe_with_attributes(ifile_name))

# Get indicies of interesting particles.
muon_event_mask = np.zeros(len(arr['particles']['id']), dtype=bool)
electron_event_mask = np.zeros(len(arr['particles']['id']), dtype=bool)
for index, array in enumerate(arr['particles']['id']):
    if 13 in abs(array):
        muon_event_mask[index] = True
    if 11 in abs(array):
        electron_event_mask[index] = True
id_array = arr['particles']['id'][((muon_event_mask)&(electron_event_mask))]
muon_index_array = [np.where(abs(array)==13)[0][0] for array in id_array]
electron_index_array = [np.where(abs(array)==11)[0][0] for array in id_array]
arr = arr[((muon_event_mask)&(electron_event_mask))]
muon_vec = arr['particles']['vector'][range(len(arr)), muon_index_array]
muon_vec["pid"] = arr['particles']['id'][range(len(arr)), muon_index_array]
electron_vec = arr['particles']['vector'][range(len(arr)), electron_index_array]
electron_vec["pid"] = arr['particles']['id'][range(len(arr)), electron_index_array]

# Cuts
muon_eta_mask = (muon_vec.eta>2.2) & (muon_vec.eta<4.4)
muon_pt_mask = (muon_vec.pt>20)
electron_eta_mask = (electron_vec.eta>2.2) & (electron_vec.eta<4.4)
electron_pt_mask = (electron_vec.pt>20)
cut_mask = (
    muon_eta_mask
    & muon_pt_mask
    & electron_eta_mask
    & electron_pt_mask
)

# Redefining arrays with cuts
muon_vec = muon_vec[cut_mask]
electron_vec = electron_vec[cut_mask]
# Calculate Additional Arrays
dilepton_vec = muon_vec + electron_vec
leading_lepton_vec = np.where(
    (muon_vec.pt > electron_vec.pt),
    muon_vec,
    electron_vec
)
trailing_lepton_vec = np.where(
    (muon_vec.pt < electron_vec.pt),
    muon_vec,
    electron_vec
)
deltaphi_array = np.abs(muon_vec.deltaphi(electron_vec))
deltaeta_array = np.abs(muon_vec.deltaeta(electron_vec))
deltar_array = np.abs(muon_vec.deltaR(electron_vec))

# Cuts
print(f"Number of Events in Tight Acceptance: {len(muon_vec)}")
# Plots
at.create_folder_path(args.output, args.testing)
os.chdir(at.find_WW_path() + "/GenLevelStudies/Figures/" + args.output)
at.create_stair(muon_vec.pt, "Muon pT", normalize=True, bins=50, range=(20, 150))
at.create_stair(muon_vec.eta, "Muon eta", normalize=True, bins=50, range=(2.2, 4.4))
at.create_stair(electron_vec.pt, "Electron pT", normalize=True, bins=50, range=(20, 150))
at.create_stair(electron_vec.eta, "Electron eta", normalize=True, bins=50, range=(2.2, 4.4))
at.create_stair(leading_lepton_vec.pt, "Leading Lepton pT", normalize=True, bins=50, range=(20, 150))
at.create_stair(leading_lepton_vec.eta, "Leading Lepton eta", normalize=True, bins=50, range=(2.2, 4.4))
at.create_stair(trailing_lepton_vec.pt, "Trailing Lepton pT", normalize=True, bins=50, range=(20, 150))
at.create_stair(trailing_lepton_vec.eta, "Trailing Lepton eta", normalize=True, bins=50, range=(2.2, 4.4))
at.create_stair(dilepton_vec.pt, "DiLepton pT", normalize=True, bins=50, range=(20, 150))
at.create_stair(dilepton_vec.m, "DiLepton mass", normalize=True, bins=50, range=(20, 300))
at.create_stair(deltaphi_array, "Delta Phi", normalize=True, bins=50, range=(0, 3.14))
at.create_stair(deltaeta_array, "Delta Eta", normalize=True, bins=50, range=(0, 3))
at.create_stair(deltar_array, "Delta R", normalize=True, bins=50, range=(0, 5))

# Histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "LeadingLeptonpT": [leading_lepton_vec.pt],
    "LeadingLeptonEta": [leading_lepton_vec.eta],
    "TrailingLeptonpT": [trailing_lepton_vec.pt],
    "TrailingLeptonEta": [trailing_lepton_vec.eta],
    "ElectronpT": [electron_vec.pt],
    "ElectronEta": [electron_vec.eta],
    "MuonpT": [muon_vec.pt],
    "MuonEta": [muon_vec.eta],
    "DiLeptonpT": [dilepton_vec.pt],
    "DiLeptonMass": [dilepton_vec.m],
    "DeltaPhi": [deltaphi_array],
    "DeltaEta": [deltaeta_array],
    "DeltaR": [deltar_array]
}
with open(args.output + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)