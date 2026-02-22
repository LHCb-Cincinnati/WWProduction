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
file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section
ofile_name = args.output
luminosity = args.luminosity # Luminosity in fb^-1

# Get Tree
with uproot.open(file_name_list[0]) as f:
    # Get the Tree
    tree = f['Tree'].arrays()

# Scale Factor
if args.cross_section == False:
    scale_factor = 1
else:
    scale_factor = args.cross_section

# DFDY invariant mass cuts
invmass_cut = 0

# Binning Scheme
ttbar_bins_list = list(np.linspace(0, 100, 4)) + [150, 200, 2000]
bins_list = ttbar_bins_list 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_pdfreweight_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_pdfreweight_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
weight_name_list = tree["pdfReweight"].fields
weighthist_dict = {}
for weight_name in weight_name_list:
    weighthist_dict[weight_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )

# Create Vectors
Z1_lminus_vec = vector.zip({
    'px': tree['Z1_lminus'].px,
    'py': tree['Z1_lminus'].py,
    'pz': tree['Z1_lminus'].pz,
    'e': tree['Z1_lminus'].e,
    'pid': tree['Z1_lminus'].pid
})

Z1_lplus_vec = vector.zip({
    'px': tree['Z1_lplus'].px,
    'py': tree['Z1_lplus'].py,
    'pz': tree['Z1_lplus'].pz,
    'e': tree['Z1_lplus'].e,
    'pid': tree['Z1_lplus'].pid
})
Z2_lminus_vec = vector.zip({
    'px': tree['Z2_lminus'].px,
    'py': tree['Z2_lminus'].py,
    'pz': tree['Z2_lminus'].pz,
    'e': tree['Z2_lminus'].e,
    'pid': tree['Z2_lminus'].pid
})

Z2_lplus_vec = vector.zip({
    'px': tree['Z2_lplus'].px,
    'py': tree['Z2_lplus'].py,
    'pz': tree['Z2_lplus'].pz,
    'e': tree['Z2_lplus'].e,
    'pid': tree['Z2_lplus'].pid
})

# Masks
# Z1_lminus and Z2_lplus pairs
Z1lmZ2lp_mue_decay_mask = (
    (abs(Z1_lminus_vec.pid)!=abs(Z2_lplus_vec.pid))
    & (Z1_lminus_vec.pid*Z2_lplus_vec.pid < 0)
)
Z1lmZ2lp_tight_acc_mask = (
    (Z1_lminus_vec.eta>2.2)
    & (Z1_lminus_vec.eta<4.4)
    & (Z2_lplus_vec.eta>2.2)
    & (Z2_lplus_vec.eta<4.4)
)
Z1lmZ2lp_high_pT_mask = (
    (Z1_lminus_vec.pt>20)
    & (Z2_lplus_vec.pt>20)
)
Z1lmZ2lp_masks = (
    Z1lmZ2lp_mue_decay_mask
    & Z1lmZ2lp_tight_acc_mask
    & Z1lmZ2lp_high_pT_mask
)
# Z2_lminus and Z1_lplus pairs
Z2lmZ1lp_mue_decay_mask = (
    (abs(Z2_lminus_vec.pid)!=abs(Z1_lplus_vec.pid))
    & (Z2_lminus_vec.pid*Z1_lplus_vec.pid < 0)
)
Z2lmZ1lp_tight_acc_mask = (
    (Z2_lminus_vec.eta>2.2)
    & (Z2_lminus_vec.eta<4.4)
    & (Z1_lplus_vec.eta>2.2)
    & (Z1_lplus_vec.eta<4.4)
)
Z2lmZ1lp_high_pT_mask = (
    (Z2_lminus_vec.pt>20)
    & (Z1_lplus_vec.pt>20)
)
Z2lmZ1lp_masks = (
    Z2lmZ1lp_mue_decay_mask
    & Z2lmZ1lp_tight_acc_mask
    & Z2lmZ1lp_high_pT_mask
)
# Apply Masks
z2lmz1lp_dilepton_vec = Z2_lminus_vec[Z2lmZ1lp_masks] + Z1_lplus_vec[Z2lmZ1lp_masks]
z2lpz1lm_dilepton_vec = Z2_lplus_vec[Z1lmZ2lp_masks] + Z1_lminus_vec[Z1lmZ2lp_masks]

# Print Info
print(f"Total Cuts: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_rghbin_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_profilehist.fill(
    z2lpz1lm_dilepton_vec.m, weight=scale_factor, sample=z2lpz1lm_dilepton_vec.m
)
dilepton_id_mass_pdfreweight_profilehist.fill(
    z2lmz1lp_dilepton_vec.m, weight=scale_factor, sample=z2lmz1lp_dilepton_vec.m
)
for weight_name in tree["pdfReweight"].fields:
    weighthist_dict[weight_name].fill(
        z2lmz1lp_dilepton_vec.m, 
        weight=tree["pdfReweight"][weight_name][Z2lmZ1lp_masks]*scale_factor
    )
    weighthist_dict[weight_name].fill(
        z2lpz1lm_dilepton_vec.m, 
        weight=tree["pdfReweight"][weight_name][Z1lmZ2lp_masks]*scale_factor
    )


# Print Statements:
print(f"Unweighted Events: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks)}")
print(f"Weighted Events: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks) * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks) * scale_factor}")


if args.debug:
    pdb.set_trace()
    exit()

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Log Rough Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Log Fine Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Linear K-Factor Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Log K-Factor Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_profilehist, "DiLepton Mass Linear K-Factor Binning Profile Hist",
                luminosity=luminosity)

# Special pdf plots
for weight_name in tree["pdfReweight"].fields:
    at.create_stacked_stair(
        [weighthist_dict[weight_name], dilepton_id_mass_pdfreweight_hist],
        f"DiLeptonMasspdfReweightBinning_{weight_name}",
        [weight_name, "CM09MTS"]
    )

fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    dilepton_id_mass_pdfreweight_hist.view().value, 
    edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
    label="CM09MTS"
)
for weight_name in tree["pdfReweight"].fields:
    axs.stairs(
        weighthist_dict[weight_name].view().value, 
        edges=weighthist_dict[weight_name].axes[0].edges,
        label=weight_name
    )
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    # root_file["DileptonKFactorFine"] = upper_env_hist

# # Save histograms
# os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
# pickle_dict = {
#     "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
#     "DileptonMassFine": [dilepton_id_mass_finebin_hist],
#     "DileptonKFactorFine": [dilepton_id_mass_pdfreweight_hist],
#     # "DileptonKFactorFine": [upper_env_hist],
#     "DileptonKFactorProfile": [dilepton_id_mass_pdfreweight_profilehist],
# }
# with open(ofile_name + ".pkl", "wb") as f:
#     pickle.dump(pickle_dict, f)
