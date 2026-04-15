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
    & lmwl_tight_acc_mask
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
    & lpwl_tight_acc_mask
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
        weight=tree["NNPDF31LO_Members"][member_name][lpwl_masks]*scale_factor
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

# RMS calculation for individual pdf families
nnpdf31lo_hist_dict = at.calc_pdf_rms(nnpdf31lo_hist_dict)
msht20lo_hist_dict = at.calc_pdf_rms(msht20lo_hist_dict)
# Special treatment for k-factor binning histogram w/ pdf uncertainties
# Create mean NLO PDF histogram from individual families
nnpdf31lo_view = weighthist_dict["nnpdf31lo"].view()
ct18lo_view = weighthist_dict["ct18lo"].view()
msht20lo_view = weighthist_dict["msht20lo"].view()
lo_mean_hist = weighthist_dict["nnpdf31lo"].copy()
lo_mean_view = lo_mean_hist.view()

# Per-bin mean of the three NLO histograms
values_stack = np.stack(
    [
        nnpdf31lo_view.value,
        ct18lo_view.value,
        msht20lo_view.value,
    ],
    axis=0,
)
mean_values = np.mean(values_stack, axis=0)
lo_mean_view.value = mean_values
sq_diff = (values_stack - mean_values) ** 2
rms_squared = np.mean(sq_diff, axis=0)
lo_mean_view.variance = np.sqrt(rms_squared)
weighthist_dict["lo_mean"] = lo_mean_hist
# Create RMS NLO PDF histogram
lo_rmslow_hist = weighthist_dict["lo_mean"].copy()
lo_rmslow_hist.view().value = lo_rmslow_hist.view().value - lo_rmslow_hist.view().variance
lo_rmshigh_hist = weighthist_dict["lo_mean"].copy()
lo_rmshigh_hist.view().value = lo_rmshigh_hist.view().value + lo_rmshigh_hist.view().variance

# Reweight Histograms
pdf_rwgt_hist_central = at.divide_bh_histograms(
    weighthist_dict["lo_mean"], 
    dilepton_id_mass_pdfreweight_hist, 
    error_type="pass_through"
)
pdf_rwgt_hist_lowerRMS = weighthist_dict["lo_mean"].copy()
pdf_rwgt_hist_lowerRMS.view().value = (
    pdf_rwgt_hist_lowerRMS.view().value
    - pdf_rwgt_hist_lowerRMS.view().variance
)
pdf_rwgt_hist_lowerRMS = at.divide_bh_histograms(
    pdf_rwgt_hist_lowerRMS, 
    dilepton_id_mass_pdfreweight_hist, 
    error_type="pass_through"
)
pdf_rwgt_hist_upperRMS = weighthist_dict["lo_mean"].copy()
pdf_rwgt_hist_upperRMS.view().value = (
    pdf_rwgt_hist_upperRMS.view().value
    + pdf_rwgt_hist_upperRMS.view().variance
)
pdf_rwgt_hist_upperRMS = at.divide_bh_histograms(
    pdf_rwgt_hist_upperRMS, 
    dilepton_id_mass_pdfreweight_hist, 
    error_type="pass_through"
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
at.create_stair(weighthist_dict["lo_mean"], "DiLepton Mass Linear K-Factor Binning Mean PDF",
                luminosity=luminosity)
at.create_stair(weighthist_dict["lo_mean"], "DiLepton Mass Linear K-Factor Binning Mean PDF",
                yscale="log", luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_central, "PDF Reweight Central Value",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_central, "PDF Reweight Central Value", yscale="log",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_lowerRMS, "PDF Reweight Lower RMS",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_lowerRMS, "PDF Reweight Lower RMS", yscale="log",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_upperRMS, "PDF Reweight Upper RMS",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_upperRMS, "PDF Reweight Upper RMS", yscale="log",
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
# PDF Reweight Plots
# Mean RMS Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    pdf_rwgt_hist_central.view().value, 
    edges=pdf_rwgt_hist_central.axes[0].edges,
    label="Central Value",
    color="black",
    zorder=3
)
axs.stairs(
    pdf_rwgt_hist_lowerRMS.view().value, 
    edges=pdf_rwgt_hist_lowerRMS.axes[0].edges,
    label="Lower Envelope",
    color="black"
)
axs.stairs(
    pdf_rwgt_hist_upperRMS.view().value, 
    edges=pdf_rwgt_hist_upperRMS.axes[0].edges,
    label="Upper Envelope",
    color="black"
)
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels, loc="lower center")
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Scales.png')
plt.close()
# NNPDF31LO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
nnpdf31lo_lowerRMS_hist, nnpdf31lo_central_hist, nnpdf31lo_upperRMS_hist = at.calc_rms_ratio_hists(nnpdf31lo_hist_dict["PDFMember0Weight"]) 
axs.stairs(
    nnpdf31lo_central_hist.view().value, 
    edges=nnpdf31lo_central_hist.axes[0].edges,
    label="Central Value",
    color="black",
    zorder=3
)
axs.stairs(
    nnpdf31lo_lowerRMS_hist.view().value, 
    edges=nnpdf31lo_lowerRMS_hist.axes[0].edges,
    label="Lower Envelope",
    color="black"
)
axs.stairs(
    nnpdf31lo_upperRMS_hist.view().value, 
    edges=nnpdf31lo_upperRMS_hist.axes[0].edges,
    label="Upper Envelope",
    color="black"
)
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels, loc="lower center")
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_NNPDF31LO.png')
plt.close()
# MSHT20LO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
msht20lo_lowerRMS_hist, msht20lo_central_hist, msht20lo_upperRMS_hist = at.calc_rms_ratio_hists(msht20lo_hist_dict["PDFMember0Weight"]) 
axs.stairs(
    msht20lo_central_hist.view().value, 
    edges=msht20lo_central_hist.axes[0].edges,
    label="Central Value",
    color="black",
    zorder=3
)
axs.stairs(
    msht20lo_lowerRMS_hist.view().value, 
    edges=msht20lo_lowerRMS_hist.axes[0].edges,
    label="Lower Envelope",
    color="black"
)
axs.stairs(
    msht20lo_upperRMS_hist.view().value, 
    edges=msht20lo_upperRMS_hist.axes[0].edges,
    label="Upper Envelope",
    color="black"
)
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels, loc="lower center")
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_MSHT20LO.png')
plt.close()
# NNPDF31LO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    weighthist_dict["lo_mean"].view().value, 
    edges=weighthist_dict["lo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black",
    zorder=3
)
axs.bar(
    x=weighthist_dict["lo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["lo_mean"].view().variance, 
    bottom=weighthist_dict["lo_mean"].view().value - weighthist_dict["lo_mean"].view().variance, 
    width=weighthist_dict["lo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean RMS Envelope",
    zorder=2
)
axs.stairs(
    nnpdf31lo_hist_dict["PDFMember0Weight"].view().value, 
    edges=nnpdf31lo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
axs.bar(
    x=nnpdf31lo_hist_dict["PDFMember0Weight"].axes[0].centers, 
    height= 2*nnpdf31lo_hist_dict["PDFMember0Weight"].view().variance, 
    bottom=(
        nnpdf31lo_hist_dict["PDFMember0Weight"].view().value 
        - nnpdf31lo_hist_dict["PDFMember0Weight"].view().variance
    ), 
    width=nnpdf31lo_hist_dict["PDFMember0Weight"].axes[0].widths, 
    linewidth=0, 
    color="green", 
    alpha=0.25, 
    label="NNPDF31LO RMS Envelope"
)
axs.stairs(
    weighthist_dict["msht20lo"].view().value, 
    edges=weighthist_dict["msht20lo"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.stairs(
    weighthist_dict["ct18lo"].view().value, 
    edges=weighthist_dict["ct18lo"].axes[0].edges,
    label="CT18LO",
    color="blue"
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 3.5))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_NNPDF31LO.png')
plt.close()
# MSHT20LO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    weighthist_dict["lo_mean"].view().value, 
    edges=weighthist_dict["lo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black",
    zorder=3
)
axs.bar(
    x=weighthist_dict["lo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["lo_mean"].view().variance, 
    bottom=weighthist_dict["lo_mean"].view().value - weighthist_dict["lo_mean"].view().variance, 
    width=weighthist_dict["lo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean RMS Envelope",
    zorder=2
)
axs.stairs(
    msht20lo_hist_dict["PDFMember0Weight"].view().value, 
    edges=msht20lo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.bar(
    x=msht20lo_hist_dict["PDFMember0Weight"].axes[0].centers, 
    height= 2*msht20lo_hist_dict["PDFMember0Weight"].view().variance, 
    bottom=(
        msht20lo_hist_dict["PDFMember0Weight"].view().value 
        - msht20lo_hist_dict["PDFMember0Weight"].view().variance
    ), 
    width=msht20lo_hist_dict["PDFMember0Weight"].axes[0].widths, 
    linewidth=0, 
    color="red", 
    alpha=0.25, 
    label="MSHT20LO RMS Envelope"
)
axs.stairs(
    weighthist_dict["nnpdf31lo"].view().value, 
    edges=weighthist_dict["nnpdf31lo"].axes[0].edges,
    label="NNPDF31LO",
    color="green"
)
axs.stairs(
    weighthist_dict["ct18lo"].view().value, 
    edges=weighthist_dict["ct18lo"].axes[0].edges,
    label="CT18LO",
    color="blue"
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 3.5))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_MSHT20LO.png')
plt.close()
# General PDF Variation Plot 
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    dilepton_id_mass_pdfreweight_hist.view().value, 
    edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
    label="CT09MCS"
)
for weight_name in tree["pdfReweight"].fields:
    axs.stairs(
        weighthist_dict[weight_name].view().value, 
        edges=weighthist_dict[weight_name].axes[0].edges,
        label=weight_name.upper()
    )
axs.stairs(
    weighthist_dict["lo_mean"].view().value, 
    edges=weighthist_dict["lo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black"
)
axs.bar(
    x=weighthist_dict["lo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["lo_mean"].view().variance, 
    bottom=weighthist_dict["lo_mean"].view().value - weighthist_dict["lo_mean"].view().variance, 
    width=weighthist_dict["lo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean w/ RMS Envelope"
)
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    root_file["DiLeptonMass_LOMean"] = weighthist_dict["lo_mean"]
    for weight_name in tree["pdfReweight"].fields:
        root_file[f"DiLeptonMass_{weight_name}"] = weighthist_dict[f"{weight_name}"]
    root_file["PDFReweight_Mean_Central"] = pdf_rwgt_hist_central
    root_file["PDFReweight_Mean_LowerRMS"] = pdf_rwgt_hist_lowerRMS
    root_file["PDFReweight_Mean_UpperRMS"] = pdf_rwgt_hist_upperRMS
    root_file["PDFReweight_MSHT20LO_Central"] = msht20lo_central_hist
    root_file["PDFReweight_MSHT20LO_LowerRMS"] = msht20lo_lowerRMS_hist
    root_file["PDFReweight_MSHT20LO_UpperRMS"] = msht20lo_upperRMS_hist
    root_file["PDFReweight_NNPDF31LO_Central"] = nnpdf31lo_central_hist
    root_file["PDFReweight_NNPDF31LO_LowerRMS"] = nnpdf31lo_lowerRMS_hist
    root_file["PDFReweight_NNPDF31LO_UpperRMS"] = nnpdf31lo_upperRMS_hist

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonpdfReweightFine": [dilepton_id_mass_pdfreweight_hist],
    "DiLeptonMass_LOMean": [weighthist_dict["lo_mean"]],
    "PDFReweight_Mean_Central": [pdf_rwgt_hist_central],
    "PDFReweight_Mean_LowerRMS": [pdf_rwgt_hist_lowerRMS],
    "PDFReweight_Mean_UpperRMS": [pdf_rwgt_hist_upperRMS],
    "PDFReweight_MSHT20LO_Central": [msht20lo_central_hist],
    "PDFReweight_MSHT20LO_LowerRMS": [msht20lo_lowerRMS_hist],
    "PDFReweight_MSHT20LO_UpperRMS": [msht20lo_upperRMS_hist],
    "PDFReweight_NNPDF31LO_Central": [nnpdf31lo_central_hist],
    "PDFReweight_NNPDF31LO_LowerRMS": [nnpdf31lo_lowerRMS_hist],
    "PDFReweight_NNPDF31LO_UpperRMS": [nnpdf31lo_upperRMS_hist],
    "DileptonpdfReweightProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
# Store reweight histograms
# for weight_name in tree["pdfReweight"].fields:
#     f"pdfReweight_{weight_name}": [weighthist_dict[f"{weight_name}_Reweight"]]

with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
