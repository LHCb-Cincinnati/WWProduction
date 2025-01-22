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
# Identified Rough Binning
leadinglepton_id_pT_rghbin_hist = bh.Histogram(bh.axis.Regular(13, 20, 150), storage=bh.storage.Weight())
leadinglepton_id_eta_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 2.2, 4.4), storage=bh.storage.Weight())
trailinglepton_id_pT_rghbin_hist = bh.Histogram(bh.axis.Regular(13, 20, 150), storage=bh.storage.Weight())
trailinglepton_id_eta_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 2.2, 4.4), storage=bh.storage.Weight())
electron_id_pT_rghbin_hist = bh.Histogram(bh.axis.Regular(13, 20, 150), storage=bh.storage.Weight())
electron_id_eta_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 2.2, 4.4), storage=bh.storage.Weight())
muon_id_pT_rghbin_hist = bh.Histogram(bh.axis.Regular(13, 20, 150), storage=bh.storage.Weight())
muon_id_eta_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 2.2, 4.4), storage=bh.storage.Weight())
dilepton_id_pT_rghbin_hist = bh.Histogram(bh.axis.Regular(13, 20, 150), storage=bh.storage.Weight())
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
deltaphi_id_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 0, 3.14), storage=bh.storage.Weight())
deltar_id_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 0, 5), storage=bh.storage.Weight())
deltaeta_id_rghbin_hist = bh.Histogram(bh.axis.Regular(6, 0, 3), storage=bh.storage.Weight())
# Identified Fine Binning
leadinglepton_id_pT_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
leadinglepton_id_eta_finebin_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
trailinglepton_id_pT_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
trailinglepton_id_eta_finebin_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
electron_id_pT_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
electron_id_eta_finebin_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
muon_id_pT_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
muon_id_eta_finebin_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
dilepton_id_pT_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
deltaphi_id_finebin_hist = bh.Histogram(bh.axis.Regular(50, 0, 3.14), storage=bh.storage.Weight())
deltar_id_finebin_hist = bh.Histogram(bh.axis.Regular(50, 0, 5), storage=bh.storage.Weight())
deltaeta_id_finebin_hist = bh.Histogram(bh.axis.Regular(50, 0, 3), storage=bh.storage.Weight())

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
both_lepton_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
) 
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
leadinglepton_id_pT_rghbin_hist.fill(leading_lepton_vec.pt / GeV)
leadinglepton_id_eta_rghbin_hist.fill(leading_lepton_vec.eta)
trailinglepton_id_pT_rghbin_hist.fill(trailing_lepton_vec.pt / GeV)
trailinglepton_id_eta_rghbin_hist.fill(trailing_lepton_vec.eta)
electron_id_pT_rghbin_hist.fill(electron_vec.pt / GeV)
electron_id_eta_rghbin_hist.fill(electron_vec.eta)
muon_id_pT_rghbin_hist.fill(muon_vec.pt / GeV)
muon_id_eta_rghbin_hist.fill(muon_vec.eta)
dilepton_id_pT_rghbin_hist.fill(dilepton_vec.pt / GeV)
dilepton_id_mass_rghbin_hist.fill(dilepton_vec.m / GeV)
deltaphi_id_rghbin_hist.fill(delta_phi_array)
deltar_id_rghbin_hist.fill(delta_r_array)
deltaeta_id_rghbin_hist.fill(delta_eta_array)
# Identified Fine Binned
leadinglepton_id_pT_finebin_hist.fill(leading_lepton_vec.pt / GeV)
leadinglepton_id_eta_finebin_hist.fill(leading_lepton_vec.eta)
trailinglepton_id_pT_finebin_hist.fill(trailing_lepton_vec.pt / GeV)
trailinglepton_id_eta_finebin_hist.fill(trailing_lepton_vec.eta)
electron_id_pT_finebin_hist.fill(electron_vec.pt / GeV)
electron_id_eta_finebin_hist.fill(electron_vec.eta)
muon_id_pT_finebin_hist.fill(muon_vec.pt / GeV)
muon_id_eta_finebin_hist.fill(muon_vec.eta)
dilepton_id_pT_finebin_hist.fill(dilepton_vec.pt / GeV)
dilepton_id_mass_finebin_hist.fill(dilepton_vec.m / GeV)
deltaphi_id_finebin_hist.fill(delta_phi_array)
deltar_id_finebin_hist.fill(delta_r_array)
deltaeta_id_finebin_hist.fill(delta_eta_array)

# Scale histograms
if not cross_section:
    scale_factor = 1
else:
    scale_factor = cross_section
# Identified Rough Binned
leadinglepton_id_pT_rghbin_hist *= scale_factor
leadinglepton_id_eta_rghbin_hist *= scale_factor
trailinglepton_id_pT_rghbin_hist *= scale_factor
trailinglepton_id_eta_rghbin_hist *= scale_factor
electron_id_pT_rghbin_hist *= scale_factor
electron_id_eta_rghbin_hist *= scale_factor
muon_id_pT_rghbin_hist *= scale_factor
muon_id_eta_rghbin_hist *= scale_factor
dilepton_id_pT_rghbin_hist *= scale_factor
dilepton_id_mass_rghbin_hist *= scale_factor
deltaphi_id_rghbin_hist *= scale_factor
deltar_id_rghbin_hist *= scale_factor
deltaeta_id_rghbin_hist *= scale_factor
# Identified Fine Binned
leadinglepton_id_pT_finebin_hist *= scale_factor
leadinglepton_id_eta_finebin_hist *= scale_factor
trailinglepton_id_pT_finebin_hist *= scale_factor
trailinglepton_id_eta_finebin_hist *= scale_factor
electron_id_pT_finebin_hist *= scale_factor
electron_id_eta_finebin_hist *= scale_factor
muon_id_pT_finebin_hist *= scale_factor
muon_id_eta_finebin_hist *= scale_factor
dilepton_id_pT_finebin_hist *= scale_factor
dilepton_id_mass_finebin_hist *= scale_factor
deltaphi_id_finebin_hist *= scale_factor
deltar_id_finebin_hist *= scale_factor
deltaeta_id_finebin_hist *= scale_factor

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
# Identified Rough Binned
at.create_stair(leadinglepton_id_pT_rghbin_hist, "Leading Lepton pT Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(leadinglepton_id_pT_rghbin_hist, "Leading Lepton pT MC Log Rough Binning",
                yscale='log', luminosity=luminosity)
at.create_stair(leadinglepton_id_eta_rghbin_hist, "Leading Lepton Eta Rough Binning",
                luminosity=luminosity)
at.create_stair(trailinglepton_id_pT_rghbin_hist, "Trailing Lepton pT Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(trailinglepton_id_pT_rghbin_hist, "Trailing Lepton pT Log Rough Binning",
                yscale='log', luminosity=luminosity)
at.create_stair(trailinglepton_id_eta_rghbin_hist, "Trailing Lepton Eta Rough Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_pT_rghbin_hist, "DiLepton pT Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_pT_rghbin_hist, "DiLepton pT Log Rough Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Log Rough Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(deltaphi_id_rghbin_hist, "Delta Phi Rough Binning",
                luminosity=luminosity)
at.create_stair(deltar_id_rghbin_hist, "Delta R Rough Binning",
                luminosity=luminosity)
at.create_stair(deltaeta_id_rghbin_hist, "Delta Eta Rough Binning",
                luminosity=luminosity)
# Identified Fine Binned
at.create_stair(leadinglepton_id_pT_finebin_hist, "Leading Lepton pT Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(leadinglepton_id_pT_finebin_hist, "Leading Lepton pT MC Log Fine Binning",
                yscale='log', luminosity=luminosity)
at.create_stair(leadinglepton_id_eta_finebin_hist, "Leading Lepton Eta Fine Binning",
                luminosity=luminosity)
at.create_stair(trailinglepton_id_pT_finebin_hist, "Trailing Lepton pT Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(trailinglepton_id_pT_finebin_hist, "Trailing Lepton pT Log Fine Binning",
                yscale='log', luminosity=luminosity)
at.create_stair(trailinglepton_id_eta_finebin_hist, "Trailing Lepton Eta Fine Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_pT_finebin_hist, "DiLepton pT Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_pT_finebin_hist, "DiLepton pT Log Fine Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Log Fine Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(deltaphi_id_finebin_hist, "Delta Phi Fine Binning",
                luminosity=luminosity)
at.create_stair(deltar_id_finebin_hist, "Delta R Fine Binning",
                luminosity=luminosity)
at.create_stair(deltaeta_id_finebin_hist, "Delta Eta Fine Binning",
                luminosity=luminosity)

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {"LeadingLeptonpTRough": [leadinglepton_id_pT_rghbin_hist],
                "LeadingLeptonEtaRough": [leadinglepton_id_eta_rghbin_hist],
                "TrailingLeptonpTRough": [trailinglepton_id_pT_rghbin_hist],
                "TrailingLeptonEtaRough": [trailinglepton_id_eta_rghbin_hist],
                "Electron pT": [electron_id_pT_rghbin_hist],
                "Electron Eta": [electron_id_eta_rghbin_hist],
                "MuonpTRough": [muon_id_pT_rghbin_hist],
                "MuonEtaRough": [muon_id_eta_rghbin_hist],
                "DiLeptonpTRough": [dilepton_id_pT_rghbin_hist],
                "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
                "DeltaPhiRough": [deltaphi_id_rghbin_hist],
                "DeltaRRough": [deltar_id_rghbin_hist],
                "DeltaEtaRough": [deltaeta_id_rghbin_hist],
                "LeadingLeptonpTFine": [leadinglepton_id_pT_finebin_hist],
                "LeadingLeptonEtaFine": [leadinglepton_id_eta_finebin_hist],
                "TrailingLeptonpTFine": [trailinglepton_id_pT_finebin_hist],
                "TrailingLeptonEtaFine": [trailinglepton_id_eta_finebin_hist],
                "ElectronpTFine": [electron_id_pT_finebin_hist],
                "ElectronEtaFine": [electron_id_eta_finebin_hist],
                "MuonpTFine": [muon_id_pT_finebin_hist],
                "MuonEtaFine": [muon_id_eta_finebin_hist],
                "DiLeptonpTFine": [dilepton_id_pT_finebin_hist],
                "DileptonMassFine": [dilepton_id_mass_finebin_hist],
                "DeltaPhiFine": [deltaphi_id_finebin_hist],
                "DeltaRFine": [deltar_id_finebin_hist],
                "DeltaEtaFine": [deltaeta_id_finebin_hist],}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
