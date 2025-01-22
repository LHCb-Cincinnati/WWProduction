
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
muon_pT_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
muon_UE_spT_hist = bh.Histogram(bh.axis.Regular(50, 0, 20), storage=bh.storage.Weight())
muon_num_cones_used_hist = bh.Histogram(bh.axis.Variable(np.linspace(-0.5, 6.5, 8)), storage=bh.storage.Weight())
muon_UE_Cone1spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_UE_Cone2spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_UE_Cone3spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_UE_Cone4spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_UE_Cone5spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_UE_Cone6spT_hist = bh.Histogram(bh.axis.Regular(30, 0, 15), storage=bh.storage.Weight())
muon_isolation_hist = bh.Histogram(bh.axis.Regular(20, 0, 20), storage=bh.storage.Weight())
muon_netisolation_hist = bh.Histogram(bh.axis.Regular(20, 0, 20), storage=bh.storage.Weight())
muon_u_hist = bh.Histogram(bh.axis.Regular(50, -100, 200), storage=bh.storage.Weight())
muon_u_lowpT_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_u_highpT_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_u_loweta_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_u_higheta_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_u_lowspT_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_u_highspT_hist = bh.Histogram(bh.axis.Regular(30, -100, 200), storage=bh.storage.Weight())
muon_eta_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_denom_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_num_hist = bh.Histogram(bh.axis.Regular(50, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_lowu_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_lowu_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_midu_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_midu_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_highu_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_highu_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_lowUEspT_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_lowUEspT_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_midUEspT_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_midUEspT_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_highUEspT_denom_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_eta_highUEspT_num_hist = bh.Histogram(bh.axis.Regular(8, 2.2, 4.4), storage=bh.storage.Weight())
muon_pT_denom_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
muon_pT_num_hist = bh.Histogram(bh.axis.Regular(50, 20, 150), storage=bh.storage.Weight())
muon_UEspT_denom_hist = bh.Histogram(bh.axis.Regular(50, 0, 20), storage=bh.storage.Weight())
muon_UEspT_num_hist = bh.Histogram(bh.axis.Regular(50, 0, 20), storage=bh.storage.Weight())
muon_u_denom_hist = bh.Histogram(bh.axis.Regular(50, -100, 200), storage=bh.storage.Weight())
muon_u_num_hist = bh.Histogram(bh.axis.Regular(50, -100, 200), storage=bh.storage.Weight())
muon_spT_u_hist = bh.Histogram(
    bh.axis.Regular(20, 0, 20),
    bh.axis.Regular(20, -100, 200),
    storage=bh.storage.Weight()
)
muon_eta_spT_hist = bh.Histogram(
    bh.axis.Regular(20, 2.2, 4.4),
    bh.axis.Regular(20, 0, 20),
    storage=bh.storage.Weight()
)
muon_eta_spT_lowpT_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, 0, 20),
    storage=bh.storage.Weight()
)
muon_eta_spT_highpT_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, 0, 20),
    storage=bh.storage.Weight()
)
muon_eta_spT_denom_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, 0, 20),
    storage=bh.storage.Weight()
)
muon_eta_spT_num_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, 0, 20),
    storage=bh.storage.Weight()
)
muon_eta_u_hist = bh.Histogram(
    bh.axis.Regular(20, 2.2, 4.4),
    bh.axis.Regular(20, -100, 200),
    storage=bh.storage.Weight()
)
muon_eta_u_lowpT_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, -100, 200),
    storage=bh.storage.Weight()
)
muon_eta_u_highpT_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, -100, 200),
    storage=bh.storage.Weight()
)
muon_eta_u_denom_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, -100, 200),
    storage=bh.storage.Weight()
)
muon_eta_u_num_hist = bh.Histogram(
    bh.axis.Regular(10, 2.2, 4.4),
    bh.axis.Regular(10, -100, 200),
    storage=bh.storage.Weight()
)

# Extra Hists
lepton_deltar_jet_hist = bh.Histogram(bh.axis.Regular(100, 0, 10), storage=bh.storage.Weight())
lepton_deltar_antijet_hist = bh.Histogram(bh.axis.Regular(100, 0, 10), storage=bh.storage.Weight())
antilepton_deltar_jet_hist = bh.Histogram(bh.axis.Regular(100, 0, 10), storage=bh.storage.Weight())
antilepton_deltar_antijet_hist = bh.Histogram(bh.axis.Regular(100, 0, 10), storage=bh.storage.Weight())

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
muon_UE_arr = ak.where((abs(lminus_vec.pid)==13), tree["TargetLeptonUE"], tree["TargetAntiLeptonUE"])
muon_u_arr = ak.where((abs(lminus_vec.pid)==13), tree["TargetLeptonu"], tree["TargetAntiLeptonu"])
muon_iso_arr = ak.where((abs(lminus_vec.pid)==13), tree["TargetLeptonIso"], tree["TargetAntiLeptonIso"])
electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)
# Masks
# Muon Masks
muon_acc = (
    (muon_vec.eta > 2.2)
    & (muon_vec.eta < 4.4)
    & (muon_vec.pt > 20)
    # & (electron_vec.eta < 4.4)
    # & (electron_vec.eta > 2.2)
    # & (electron_vec.pt > 20)
)
# Apply Masks
muon_tree = tree[muon_acc]
muon_vec = muon_vec[muon_acc]
muon_UE_arr = muon_UE_arr[muon_acc]
muon_iso_arr = muon_iso_arr[muon_acc]
muon_u_arr = muon_u_arr[muon_acc]

# Calculate Quantitites
muon_netiso_arr = muon_iso_arr - muon_UE_arr.spTAve

# New Masks 
muon_iso_cut = (muon_netiso_arr < 2.5)
muon_lowu_cut = (
    (muon_u_arr < -40)
    & (muon_u_arr > -45) 
)
muon_midu_cut = (
    (muon_u_arr > 25)
    & (muon_u_arr < 30) 
)
muon_highu_cut = (
    (muon_u_arr > 95)
    & (muon_u_arr < 100) 
)
muon_lowUEspT_cut = (
    (muon_UE_arr.spTAve > 0.5)
    & (muon_UE_arr.spTAve < 0.75) 
)
muon_midUEspT_cut = (
    (muon_UE_arr.spTAve > 3)
    & (muon_UE_arr.spTAve < 3.5) 
)
muon_highUEspT_cut = (
    (muon_UE_arr.spTAve > 7)
    & (muon_UE_arr.spTAve < 7.5) 
)

# Fill histograms
muon_pT_hist.fill(muon_vec.pt / GeV) 
muon_UE_spT_hist.fill(muon_UE_arr.spTAve)
muon_num_cones_used_hist.fill(muon_UE_arr.NumGoodCones)
muon_UE_Cone1spT_hist.fill(muon_UE_arr.Cone1spT)
muon_UE_Cone2spT_hist.fill(muon_UE_arr.Cone2spT)
muon_UE_Cone3spT_hist.fill(muon_UE_arr.Cone3spT)
muon_UE_Cone4spT_hist.fill(muon_UE_arr.Cone4spT)
muon_UE_Cone5spT_hist.fill(muon_UE_arr.Cone5spT)
muon_UE_Cone6spT_hist.fill(muon_UE_arr.Cone6spT)
muon_netisolation_hist.fill(muon_netiso_arr)
muon_isolation_hist.fill(muon_iso_arr)
muon_u_hist.fill(muon_u_arr)
muon_u_lowpT_hist.fill(muon_u_arr[(muon_vec.pt / GeV < 50)])
muon_u_highpT_hist.fill(muon_u_arr[(muon_vec.pt / GeV > 50)])
muon_u_loweta_hist.fill(muon_u_arr[(muon_vec.eta < 3)])
muon_u_higheta_hist.fill(muon_u_arr[(muon_vec.eta > 3.5)])
muon_u_lowspT_hist.fill(muon_u_arr[(muon_UE_arr.spTAve < 5)])
muon_u_highspT_hist.fill(muon_u_arr[(muon_UE_arr.spTAve > 5)])
muon_spT_u_hist.fill(muon_UE_arr.spTAve, muon_u_arr)
muon_eta_spT_lowpT_hist.fill(muon_vec.eta[(muon_vec.pt / GeV > 50)], muon_UE_arr.spTAve[(muon_vec.pt / GeV > 50)])
muon_eta_spT_highpT_hist.fill(muon_vec.eta[(muon_vec.pt / GeV < 50)], muon_UE_arr.spTAve[(muon_vec.pt / GeV < 50)])
muon_eta_spT_hist.fill(muon_vec.eta, muon_UE_arr.spTAve)
muon_eta_spT_denom_hist.fill(muon_vec.eta, muon_UE_arr.spTAve)
muon_eta_spT_num_hist.fill(muon_vec[muon_iso_cut].eta, muon_UE_arr[muon_iso_cut].spTAve)
muon_eta_u_hist.fill(muon_vec.eta, muon_u_arr)
muon_eta_u_lowpT_hist.fill(muon_vec.eta[(muon_vec.pt / GeV < 50)], muon_u_arr[(muon_vec.pt / GeV < 50)])
muon_eta_u_highpT_hist.fill(muon_vec.eta[(muon_vec.pt / GeV > 50)], muon_u_arr[(muon_vec.pt / GeV > 50)])
muon_eta_u_denom_hist.fill(muon_vec.eta, muon_u_arr)
muon_eta_u_num_hist.fill(muon_vec[muon_iso_cut].eta, muon_u_arr[muon_iso_cut])
muon_eta_hist.fill(muon_vec.eta)
muon_eta_denom_hist.fill(muon_vec.eta)
muon_eta_num_hist.fill(muon_vec[muon_iso_cut].eta)
muon_eta_lowu_denom_hist.fill(muon_vec[muon_lowu_cut].eta)
muon_eta_lowu_num_hist.fill(muon_vec[muon_iso_cut&muon_lowu_cut].eta)
muon_eta_midu_denom_hist.fill(muon_vec[muon_midu_cut].eta)
muon_eta_midu_num_hist.fill(muon_vec[muon_iso_cut&muon_midu_cut].eta)
muon_eta_highu_denom_hist.fill(muon_vec[muon_highu_cut].eta)
muon_eta_highu_num_hist.fill(muon_vec[muon_iso_cut&muon_highu_cut].eta)
muon_eta_lowUEspT_denom_hist.fill(muon_vec[muon_lowUEspT_cut].eta)
muon_eta_lowUEspT_num_hist.fill(muon_vec[muon_iso_cut&muon_lowUEspT_cut].eta)
muon_eta_midUEspT_denom_hist.fill(muon_vec[muon_midUEspT_cut].eta)
muon_eta_midUEspT_num_hist.fill(muon_vec[muon_iso_cut&muon_midUEspT_cut].eta)
muon_eta_highUEspT_denom_hist.fill(muon_vec[muon_highUEspT_cut].eta)
muon_eta_highUEspT_num_hist.fill(muon_vec[muon_iso_cut&muon_highUEspT_cut].eta)
muon_pT_denom_hist.fill(muon_vec.pt / GeV)
muon_pT_num_hist.fill(muon_vec[muon_iso_cut].pt / GeV)
muon_UEspT_denom_hist.fill(muon_UE_arr.spTAve)
muon_UEspT_num_hist.fill(muon_UE_arr[muon_iso_cut].spTAve)
muon_u_denom_hist.fill(muon_u_arr)
muon_u_num_hist.fill(muon_u_arr[muon_iso_cut])

# Divide histograms
muon_eta_eff_hist = at.divide_bh_histograms(
    muon_eta_num_hist,
    muon_eta_denom_hist,
    binomial_error = True
)
muon_eta_lowu_eff_hist = at.divide_bh_histograms(
    muon_eta_lowu_num_hist,
    muon_eta_lowu_denom_hist,
    binomial_error = True
)
muon_eta_midu_eff_hist = at.divide_bh_histograms(
    muon_eta_midu_num_hist,
    muon_eta_midu_denom_hist,
    binomial_error = True
)
muon_eta_highu_eff_hist = at.divide_bh_histograms(
    muon_eta_highu_num_hist,
    muon_eta_highu_denom_hist,
    binomial_error = True
)
muon_eta_lowUEspT_eff_hist = at.divide_bh_histograms(
    muon_eta_lowUEspT_num_hist,
    muon_eta_lowUEspT_denom_hist,
    binomial_error = True
)
muon_eta_midUEspT_eff_hist = at.divide_bh_histograms(
    muon_eta_midUEspT_num_hist,
    muon_eta_midUEspT_denom_hist,
    binomial_error = True
)
muon_eta_highUEspT_eff_hist = at.divide_bh_histograms(
    muon_eta_highUEspT_num_hist,
    muon_eta_highUEspT_denom_hist,
    binomial_error = True
)
muon_UEspT_eff_hist = at.divide_bh_histograms(
    muon_UEspT_num_hist,
    muon_UEspT_denom_hist,
    binomial_error = True
)
muon_pT_eff_hist = at.divide_bh_histograms(
    muon_pT_num_hist,
    muon_pT_denom_hist,
    binomial_error = True
)
muon_u_eff_hist = at.divide_bh_histograms(
    muon_u_num_hist,
    muon_u_denom_hist,
    binomial_error = True
)
muon_eta_u_eff_hist = at.divide_bh_histograms(
    muon_eta_u_num_hist,
    muon_eta_u_denom_hist,
    binomial_error = True
)
muon_eta_spT_eff_hist = at.divide_bh_histograms(
    muon_eta_spT_num_hist,
    muon_eta_spT_denom_hist,
    binomial_error = True
)

# # Extra Hists
# lepton_deltar_jet_hist.fill(
#     np.sqrt(
#         (lminus_vec.eta - tree["TargetJet"].eta) **2
#         + (lminus_vec.eta - tree["TargetJet"].phi) **2
#     )
# )
# lepton_deltar_antijet_hist.fill(
#     np.sqrt(
#         (lminus_vec.eta - tree["TargetAntiJet"].eta) **2
#         + (lminus_vec.eta - tree["TargetAntiJet"].phi) **2
#     )
# )
# antilepton_deltar_jet_hist.fill(
#     np.sqrt(
#         (lplus_vec.eta - tree["TargetJet"].eta) **2
#         + (lplus_vec.eta - tree["TargetJet"].phi) **2
#     )
# )
# antilepton_deltar_antijet_hist.fill(
#     np.sqrt(
#         (lplus_vec.eta - tree["TargetAntiJet"].eta) **2
#         + (lplus_vec.eta - tree["TargetAntiJet"].phi) **2
#     )
# )

# Print Info
if args.debug:
    exit()

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
# Identified Rough Binned
at.create_stair(
    muon_pT_hist, 
    "Muon pT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_eta_hist, 
    "Muon Eta",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_spT_hist, 
    "Muon UE spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_num_cones_used_hist, 
    "Muon Num Good Cones",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone1spT_hist, 
    "Muon UE Cone 1 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone2spT_hist, 
    "Muon UE Cone 2 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone3spT_hist, 
    "Muon UE Cone 3 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone4spT_hist, 
    "Muon UE Cone 4 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone5spT_hist, 
    "Muon UE Cone 5 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_UE_Cone6spT_hist, 
    "Muon UE Cone 6 spT",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_isolation_hist, 
    "Muon Isolation",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_netisolation_hist, 
    "Muon Net Isolation",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_u_hist, 
    "Muon u",
    normalize=True,
    luminosity=luminosity
)
at.create_stair(
    muon_u_lowpT_hist,
    "Muon u Low pT",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_u_highpT_hist,
    "Muon u High pT",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_u_loweta_hist,
    "Muon u Low Eta",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_u_higheta_hist,
    "Muon u High Eta",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_u_lowspT_hist,
    "Muon u Low spT",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_u_highspT_hist,
    "Muon u High spT",
    normalize = True,
    luminosity = luminosity
)
at.create_stair(
    muon_eta_eff_hist, 
    "Muon Isolation Efficiency (Eta)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_lowu_eff_hist, 
    "Muon Isolation Efficiency (Eta Lowu)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_midu_eff_hist, 
    "Muon Isolation Efficiency (Eta Midu)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_highu_eff_hist, 
    "Muon Isolation Efficiency (Eta Highu)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_lowUEspT_eff_hist, 
    "Muon Isolation Efficiency (Eta LowUEspT)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_midUEspT_eff_hist, 
    "Muon Isolation Efficiency (Eta MidUEspT)",
    luminosity=luminosity
)
at.create_stair(
    muon_eta_highUEspT_eff_hist, 
    "Muon Isolation Efficiency (Eta HighUEspT)",
    luminosity=luminosity
)
at.create_stair(
    muon_UEspT_eff_hist, 
    "Muon Isolation Efficiency (UEspT)",
    luminosity=luminosity
)
at.create_stair(
    muon_pT_eff_hist, 
    "Muon Isolation Efficiency (pT)",
    luminosity=luminosity
)
at.create_stair(
    muon_u_eff_hist, 
    "Muon Isolation Efficiency (u)",
    luminosity=luminosity
)
at.plot_2dhist(
    muon_spT_u_hist,
    "Muon spT u",
    normalize = True
)
at.plot_2dhist(
    muon_eta_spT_hist,
    "Muon Eta spT",
    normalize = True
)
at.plot_2dhist(
    muon_eta_spT_lowpT_hist,
    "Muon Eta spT lowPT",
    normalize = True
)
at.plot_2dhist(
    muon_eta_spT_highpT_hist,
    "Muon Eta spT high pT",
    normalize = True
)
at.plot_2dhist(
    muon_eta_u_hist,
    "Muon Eta u",
    normalize = True
)
at.plot_2dhist(
    muon_eta_u_lowpT_hist,
    "Muon Eta u LowpT",
    normalize = True
)
at.plot_2dhist(
    muon_eta_u_highpT_hist,
    "Muon Eta u highpT",
    normalize = True
)
at.plot_2dhist(
    muon_eta_u_eff_hist,
    "Muon Eta u Isolation Eff",
    vmin = 0,
    vmax = 1
)
at.plot_2dhist(
    muon_eta_spT_eff_hist,
    "Muon Eta spT Isolation Eff",
    vmin = 0,
    vmax = 1
)

# # Extra Hists
# at.create_stair(
#     lepton_deltar_antijet_hist,
#     "Lepton AntiJet DeltaR",
#     normalize = True
# )
# at.create_stair(
#     lepton_deltar_jet_hist,
#     "Lepton Jet DeltaR",
#     normalize = True
# )
# at.create_stair(
#     antilepton_deltar_antijet_hist,
#     "AntiLepton AntiJet DeltaR",
#     normalize = True
# )
# at.create_stair(
#     antilepton_deltar_jet_hist,
#     "AntiLepton Jet DeltaR",
#     normalize = True
# )

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "MuonpT": [muon_pT_hist],
    "MuonEta": [muon_eta_hist],
    "MuonUEspT": [muon_UE_spT_hist],
    "MuonNumGoodCones": [muon_num_cones_used_hist],
    "MuonUECone1spT": [muon_UE_Cone1spT_hist],
    "MuonUECone2spT": [muon_UE_Cone2spT_hist],
    "MuonUECone3spT": [muon_UE_Cone3spT_hist],
    "MuonUECone4spT": [muon_UE_Cone4spT_hist],
    "MuonUECone5spT": [muon_UE_Cone5spT_hist],
    "MuonUECone6spT": [muon_UE_Cone6spT_hist],
    "MuonIso": [muon_isolation_hist],
    "MuonnetIso": [muon_netisolation_hist],
    "Muonu": [muon_u_hist],
    "MuonuLowpT": [muon_u_lowpT_hist],
    "MuonuHighpT": [muon_u_highpT_hist],
    "MuonuLowEta": [muon_u_loweta_hist],
    "MuonuHighEta": [muon_u_higheta_hist],
    "MuonuLowspT": [muon_u_lowspT_hist],
    "MuonuHighspT": [muon_u_highspT_hist],
    "MuonEtaIsoEff": [muon_eta_eff_hist],
    "MuonEtaLowuIsoEff": [muon_eta_lowu_eff_hist],
    "MuonEtaMiduIsoEff": [muon_eta_midu_eff_hist],
    "MuonEtaHighuIsoEff": [muon_eta_highu_eff_hist],
    "MuonEtaLowUEspTIsoEff": [muon_eta_lowUEspT_eff_hist],
    "MuonEtaMidUEspTIsoEff": [muon_eta_midUEspT_eff_hist],
    "MuonEtaHighUEspTIsoEff": [muon_eta_highUEspT_eff_hist],
    "MuonpTIsoEff": [muon_pT_eff_hist],
    "MuonUEspTIsoEff": [muon_UEspT_eff_hist],
    "MuonuEff": [muon_u_eff_hist]
}

with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)