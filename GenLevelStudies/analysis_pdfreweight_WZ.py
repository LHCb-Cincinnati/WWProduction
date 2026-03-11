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

# Variables
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
nnpdf31lo_hist_dict = {}
for member_name in tree["NNPDF31LO_Members"].fields:
    nnpdf31lo_hist_dict[member_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
msht20lo_hist_dict = {}
for member_name in tree["MSHT20LO_Members"].fields:
    msht20lo_hist_dict[member_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )


# Create Vectors
lminus_vec = vector.zip({
    'px': tree['Z_lminus'].px,
    'py': tree['Z_lminus'].py,
    'pz': tree['Z_lminus'].pz,
    'e': tree['Z_lminus'].e,
    'pid': tree['Z_lminus'].pid
})

lplus_vec = vector.zip({
    'px': tree['Z_lplus'].px,
    'py': tree['Z_lplus'].py,
    'pz': tree['Z_lplus'].pz,
    'e': tree['Z_lplus'].e,
    'pid': tree['Z_lplus'].pid
})

thirdl_vec = vector.zip({
    'px': tree['W_lepton'].px,
    'py': tree['W_lepton'].py,
    'pz': tree['W_lepton'].pz,
    'e': tree['W_lepton'].e,
    'pid': tree['W_lepton'].pid
})

# Masks
# general masks
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt>15)
                         | (thirdl_vec.eta>1.596) & (thirdl_vec.pt>15))
# lminus and wlepton pair masks
lmwl_mue_decay_mask = (
    (abs(lminus_vec.pid)!=abs(thirdl_vec.pid))
    & (lminus_vec.pid*thirdl_vec.pid < 0)
)
lmwl_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lmwl_high_pT_mask = (
    (lminus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lmwl_masks = (
    lmwl_mue_decay_mask
    # & lmwl_tight_acc_mask
    & lmwl_high_pT_mask
)
# lplus and wlepton pair masks
lpwl_mue_decay_mask = (
    (abs(lplus_vec.pid)!=abs(thirdl_vec.pid))
    & (lplus_vec.pid*thirdl_vec.pid < 0)
)
lpwl_tight_acc_mask = (
    (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lpwl_high_pT_mask = (
    (lplus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lpwl_masks = (
    lpwl_mue_decay_mask
    # & lpwl_tight_acc_mask
    & lpwl_high_pT_mask
)
# Apply Masks
lpwl_dilepton_vec = lplus_vec[lpwl_masks] + thirdl_vec[lpwl_masks]
lmwl_dilepton_vec = lplus_vec[lmwl_masks] + thirdl_vec[lmwl_masks]

# Print Info
print(f"Total Cuts: {sum(lpwl_masks | lmwl_masks)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_rghbin_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_pdfreweight_profilehist.fill(
    lmwl_dilepton_vec.m, weight=scale_factor, sample=lmwl_dilepton_vec.m
)
dilepton_id_mass_pdfreweight_profilehist.fill(
    lpwl_dilepton_vec.m, weight=scale_factor, sample=lpwl_dilepton_vec.m
)
for weight_name in tree["pdfReweight"].fields:
    weighthist_dict[weight_name].fill(
        lpwl_dilepton_vec.m, 
        weight=tree["pdfReweight"][weight_name][lpwl_masks]*scale_factor
    )
    weighthist_dict[weight_name].fill(
        lmwl_dilepton_vec.m, 
        weight=tree["pdfReweight"][weight_name][lmwl_masks]*scale_factor
    )
for member_name in tree["NNPDF31LO_Members"].fields:
    nnpdf31lo_hist_dict[member_name].fill(
        lpwl_dilepton_vec.m, 
        weight=tree["NNPDF31LO_Members"][member_name][lpwl_masks]
    )
    nnpdf31lo_hist_dict[member_name].fill(
        lmwl_dilepton_vec.m, 
        weight=tree["NNPDF31LO_Members"][member_name][lmwl_masks]*scale_factor
    )
for member_name in tree["MSHT20LO_Members"].fields:
    msht20lo_hist_dict[member_name].fill(
        lpwl_dilepton_vec.m, 
        weight=tree["MSHT20LO_Members"][member_name][lpwl_masks]*scale_factor
    )
    msht20lo_hist_dict[member_name].fill(
        lmwl_dilepton_vec.m, 
        weight=tree["MSHT20LO_Members"][member_name][lmwl_masks]*scale_factor
    )

# Print Statements:
print(f"Unweighted Events: {sum(lpwl_masks | lmwl_masks)}")
print(f"Weighted Events: {sum(lpwl_masks | lmwl_masks) * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * sum(lpwl_masks | lmwl_masks) * scale_factor}")


if args.debug:
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
# NNPDF31LO RMS Plots
at.create_stacked_stair(
    [*at.calc_pdf_rms(nnpdf31lo_hist_dict), weighthist_dict["nnpdf31lo"]],
    f"DiLeptonMass_PDFReweight_NNPD31LO",
    ["Upper Envelope", "Lower Envelope", "Central Value"],
    yscale="log"
)
# MSHT20LO RMS Plots
at.create_stacked_stair(
    [*at.calc_pdf_rms(msht20lo_hist_dict), weighthist_dict["msht20lo"]],
    f"DiLeptonMass_PDFReweight_MSHT20LO",
    ["Upper Envelope", "Lower Envelope", "Central Value"],
    yscale="log"
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
axs.set_ylabel("$ \\frac{d \\sigma}{M_{e \\mu}}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    for weight_name in tree["pdfReweight"].fields:
        root_file[f"pdfReweight_{weight_name}"] =  weighthist_dict[f"{weight_name}"]

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonpdfReweightFine": [dilepton_id_mass_pdfreweight_hist],
    "DileptonpdfReweightProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
# Store reweight histograms
# for weight_name in tree["pdfReweight"].fields:
#     f"pdfReweight_{weight_name}": [weighthist_dict[f"{weight_name}_Reweight"]]

with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
