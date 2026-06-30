
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
import mplhep as hep
# Personal Packages
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
weights_bool = True
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
pdf4lhc_hist_dict = {}
for member_name in tree["PDF4LHC21_Members"].fields:
    pdf4lhc_hist_dict[member_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['TargetAntiLepton'].px * GeV,
    'py': tree['TargetAntiLepton'].py * GeV,
    'pz': tree['TargetAntiLepton'].pz * GeV,
    'e': tree['TargetAntiLepton'].e * GeV,
    'pid': tree['TargetAntiLepton'].pid
})
lplus_vec = vector.zip({
    'px': tree['TargetLepton'].px * GeV,
    'py': tree['TargetLepton'].py * GeV,
    'pz': tree['TargetLepton'].pz * GeV,
    'e': tree['TargetLepton'].e * GeV,
    'pid': tree['TargetLepton'].pid
})
dilepton_vec = lminus_vec + lplus_vec
muon_vec = ak.where((abs(lminus_vec.pid)==13), lminus_vec, lplus_vec)
electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

# Masks
both_lepton_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
)
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
deltar_mask = (np.abs(lminus_vec.deltaR(lplus_vec)) > 0.1)
high_pT_lepton_mask = ((electron_vec.pt / GeV >30) & (muon_vec.pt / GeV >25))
lepton_mask = (
    both_lepton_tight_acc_mask
    & high_pT_lepton_mask
    & mue_decay_mask
    & deltar_mask
)

# Apply Masks
# Note here that muon/electron vec and dilepton vec have different lengths.
muon_vec = muon_vec[mue_decay_mask & high_pT_lepton_mask & deltar_mask]
electron_vec = electron_vec[mue_decay_mask & high_pT_lepton_mask & deltar_mask]
dilepton_vec = dilepton_vec[lepton_mask]

# Print Info
print(f"Total Cuts: {sum(lepton_mask)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_pdfreweight_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_pdfreweight_profilehist.fill(
    dilepton_vec.m / GeV, weight=scale_factor, sample=dilepton_vec.m / GeV
)
for weight_name in tree["pdfReweight"].fields:
    weighthist_dict[weight_name].fill(
        dilepton_vec.m / GeV,
        weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
    )
for member_name in tree["NNPDF31LO_Members"].fields:
    nnpdf31lo_hist_dict[member_name].fill(
        dilepton_vec.m / GeV,
        weight=tree["NNPDF31LO_Members"][member_name][lepton_mask]*scale_factor
    )
for member_name in tree["MSHT20LO_Members"].fields:
    msht20lo_hist_dict[member_name].fill(
        dilepton_vec.m / GeV,
        weight=tree["MSHT20LO_Members"][member_name][lepton_mask]*scale_factor
    )
for member_name in tree["PDF4LHC21_Members"].fields:
    pdf4lhc_hist_dict[member_name].fill(
        dilepton_vec.m / GeV,
        weight=tree["PDF4LHC21_Members"][member_name][lepton_mask]*scale_factor
    )

# Store statistical variances before rewriting
variance_stack = np.stack(
    [
        weighthist_dict["nnpdf31lo"].view().variance,
        weighthist_dict["ct18lo"].view().variance,
        weighthist_dict["msht20lo"].view().variance,
    ],
    axis=0,
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
# Create NNPDF31LO RMS LO PDF histogram
nnpdf31lo_rmslow_hist = nnpdf31lo_hist_dict["PDFMember0Weight"].copy()
nnpdf31lo_rmslow_hist.view().value = nnpdf31lo_rmslow_hist.view().value - nnpdf31lo_rmslow_hist.view().variance
nnpdf31lo_rmslow_hist.view().variance = weighthist_dict["nnpdf31lo"].view().variance
nnpdf31lo_rmshigh_hist = nnpdf31lo_hist_dict["PDFMember0Weight"].copy()
nnpdf31lo_rmshigh_hist.view().value = nnpdf31lo_rmshigh_hist.view().value + nnpdf31lo_rmshigh_hist.view().variance
nnpdf31lo_rmshigh_hist.view().variance = weighthist_dict["nnpdf31lo"].view().variance
nnpdf31lo_rmslow_hist.view().value[-1] = 0.0
# Create MSHT20LO RMS LO PDF histogram
msht20lo_rmslow_hist = msht20lo_hist_dict["PDFMember0Weight"].copy()
msht20lo_rmslow_hist.view().value = msht20lo_rmslow_hist.view().value - msht20lo_rmslow_hist.view().variance
msht20lo_rmslow_hist.view().variance = weighthist_dict["msht20lo"].view().variance
msht20lo_rmshigh_hist = msht20lo_hist_dict["PDFMember0Weight"].copy()
msht20lo_rmshigh_hist.view().value = msht20lo_rmshigh_hist.view().value + msht20lo_rmshigh_hist.view().variance
msht20lo_rmshigh_hist.view().variance = weighthist_dict["msht20lo"].view().variance
msht20lo_rmslow_hist.view().value[-1] = 0.0
pdf4lhc_central_hist, pdf4lhc_lowerenv_hist, pdf4lhc_upperenv_hist = at.calc_pdf_mc_envelope(
    pdf4lhc_hist_dict
)

# Reweight Histograms
lo_mean_hist_wstaterrs = weighthist_dict["lo_mean"].copy()
lo_mean_hist_wstaterrs.view().variance = np.sum(variance_stack) / 9
pdf_rwgt_hist_central = at.divide_bh_histograms(
    lo_mean_hist_wstaterrs,
    dilepton_id_mass_pdfreweight_hist
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
pdf4lhc_rwgt_hist_central = at.divide_bh_histograms(
    pdf4lhc_central_hist,
    dilepton_id_mass_pdfreweight_hist
)

# Print Statements:
print(f"Unweighted Events: {len(dilepton_vec)}")
print(f"Weighted Events: {len(dilepton_vec) * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * len(dilepton_vec) * scale_factor}")

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
at.create_stair(weighthist_dict["lo_mean"], "DiLepton Mass Linear K-Factor Binning Mean PDF",
                luminosity=luminosity)
at.create_stair(weighthist_dict["lo_mean"], "DiLepton Mass Linear K-Factor Binning Mean PDF",
                yscale="log", luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_central, "PDF Reweight Central Value",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_central, "PDF Reweight Central Value Log", yscale="log",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_lowerRMS, "PDF Reweight Lower RMS",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_lowerRMS, "PDF Reweight Lower RMS Log", yscale="log",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_upperRMS, "PDF Reweight Upper RMS",
                luminosity=luminosity)
at.create_stair(pdf_rwgt_hist_upperRMS, "PDF Reweight Upper RMS Log", yscale="log",
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
# Mean PDF Reweight Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    pdf_rwgt_hist_central.view().value,
    edges=pdf_rwgt_hist_central.axes[0].edges,
#     label="Central Value",
    color="black",
    zorder=3
)
axs.errorbar(
    pdf_rwgt_hist_central.axes[0].centers,
    pdf_rwgt_hist_central.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = np.sqrt(pdf_rwgt_hist_central.view().variance),
    label="RMS Deviation"
)
# axs.stairs(
#     pdf_rwgt_hist_lowerRMS.view().value,
#     edges=pdf_rwgt_hist_lowerRMS.axes[0].edges,
#     label="Lower Envelope",
#     color="black"
# )
# axs.stairs(
#     pdf_rwgt_hist_upperRMS.view().value,
#     edges=pdf_rwgt_hist_upperRMS.axes[0].edges,
#     label="Upper Envelope",
#     color="black"
# )
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Reweight Factor")
# hist_handles, hist_labels = axs.get_legend_handles_labels()
# axs.legend(hist_handles, hist_labels, loc="lower center")
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Reweight.png')
plt.close()
# PDF4LHC Reweight Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    pdf4lhc_rwgt_hist_central.view().value,
    edges=pdf4lhc_rwgt_hist_central.axes[0].edges,
    color="black",
    zorder=3
)
axs.errorbar(
    pdf4lhc_rwgt_hist_central.axes[0].centers,
    pdf4lhc_rwgt_hist_central.view().value,
    ecolor="black",
    linestyle="",
    yerr=np.sqrt(pdf4lhc_rwgt_hist_central.view().variance),
    label="Statistical Uncertainty"
)
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF4LHC Reweight Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('PDF4LHCReweightHistogram.png')
plt.close()
# Mean PDF Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
meanpdflo_lowerratio_hist, meanpdflo_centralratio_hist, meanpdflo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["lo_mean"], lo_rmslow_hist, lo_rmshigh_hist
)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=meanpdflo_centralratio_hist.axes[0].centers,
    height=(meanpdflo_upperratio_hist.view().value - meanpdflo_lowerratio_hist.view().value),
    bottom=meanpdflo_lowerratio_hist.view().value,
    width=meanpdflo_centralratio_hist.axes[0].widths,
    linewidth=0,
    color="grey",
    alpha=0.25,
    label="RMS Deviation"
)
# axs.stairs(
#     meanpdflo_central_hist.view().value,
#     edges=meanpdflo_central_hist.axes[0].edges,
#     label="Central Value",
#     color="black",
#     zorder=3
# )
# axs.stairs(
#     meanpdflo_rmslow_hist.view().value,
#     edges=meanpdflo_rmslow_hist.axes[0].edges,
#     label="Lower Envelope",
#     color="black"
# )
# axs.stairs(
#     meanpdflo_rmshigh_hist.view().value,
#     edges=meanpdflo_rmshigh_hist.axes[0].edges,
#     label="Upper Envelope",
#     color="black"
# )
axs.set_xlim((0, 300))
axs.set_ylim((0, 4))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_MeanPDF.png')
plt.close()
# NNPDF31LO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
nnpdf31lo_lowerratio_hist, nnpdf31lo_centralratio_hist, nnpdf31lo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["nnpdf31lo"], nnpdf31lo_rmslow_hist, nnpdf31lo_rmshigh_hist
)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=nnpdf31lo_centralratio_hist.axes[0].centers,
    height=(nnpdf31lo_upperratio_hist.view().value - nnpdf31lo_lowerratio_hist.view().value),
    bottom=nnpdf31lo_lowerratio_hist.view().value,
    width=nnpdf31lo_centralratio_hist.axes[0].widths,
    linewidth=0,
    color="green",
    alpha=0.25,
    label="RMS Deviation"
)
# axs.stairs(
#     nnpdf31lo_central_hist.view().value,
#     edges=nnpdf31lo_central_hist.axes[0].edges,
#     label="Central Value",
#     color="black",
#     zorder=3
# )
# axs.stairs(
#     nnpdf31lo_rmslow_hist.view().value,
#     edges=nnpdf31lo_rmslow_hist.axes[0].edges,
#     label="Lower Envelope",
#     color="black"
# )
# axs.stairs(
#     nnpdf31lo_rmshigh_hist.view().value,
#     edges=nnpdf31lo_rmshigh_hist.axes[0].edges,
#     label="Upper Envelope",
#     color="black"
# )
axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_NNPDF31LO.png')
plt.close()
# MSHT20LO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
msht20lo_lowerratio_hist, msht20lo_centralratio_hist, msht20lo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["msht20lo"], msht20lo_rmslow_hist, msht20lo_rmshigh_hist
)
pdf4lhc_lowerratio_hist, pdf4lhc_centralratio_hist, pdf4lhc_upperratio_hist = at.calc_rms_ratio_hists(
    pdf4lhc_central_hist, pdf4lhc_lowerenv_hist, pdf4lhc_upperenv_hist
)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=msht20lo_centralratio_hist.axes[0].centers,
    height=(msht20lo_upperratio_hist.view().value - msht20lo_lowerratio_hist.view().value),
    bottom=msht20lo_lowerratio_hist.view().value,
    width=msht20lo_centralratio_hist.axes[0].widths,
    linewidth=0,
    color="red",
    alpha=0.25,
    label="RMS Deviation"
)
# axs.stairs(
#     msht20lo_central_hist.view().value,
#     edges=msht20lo_central_hist.axes[0].edges,
#     label="Central Value",
#     color="black",
#     zorder=3
# )
# axs.stairs(
#     msht20lo_lowerRMS_hist.view().value,
#     edges=msht20lo_lowerRMS_hist.axes[0].edges,
#     label="Lower Envelope",
#     color="black"
# )
# axs.stairs(
#     msht20lo_upperRMS_hist.view().value,
#     edges=msht20lo_upperRMS_hist.axes[0].edges,
#     label="Upper Envelope",
#     color="black"
# )
# axs.bar(
#     x=msht20lo_upperRMS_hist.axes[0].centers,
#     height=(msht20lo_upperRMS_hist.view().value - msht20lo_lowerRMS_hist.view().value),
#     bottom=msht20lo_upperRMS_hist.view().value,
#     width=msht20lo_upperRMS_hist.axes[0].widths,
#     # align='edge',
#     linewidth=0,
#     color="red",
#     alpha=0.25,
#     # zorder=-1,
#     label="NNPDF31NLO Envelope"
# )
axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_MSHT20LO.png')
plt.close()
# PDF4LHC Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=pdf4lhc_centralratio_hist.axes[0].centers,
    height=(pdf4lhc_upperratio_hist.view().value - pdf4lhc_lowerratio_hist.view().value),
    bottom=pdf4lhc_lowerratio_hist.view().value,
    width=pdf4lhc_centralratio_hist.axes[0].widths,
    linewidth=0,
    color="purple",
    alpha=0.25,
    label="PDF4LHC 68% CL Envelope"
)
axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('DiLeptonMassPDFVariations_Ratio_PDF4LHC.png')
plt.close()
# NNPDF31LO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    dilepton_id_mass_pdfreweight_hist.view().value,
    edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
    label="CT09MCS",
    color="grey"
)
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
axs.stairs(
    pdf4lhc_central_hist.view().value,
    edges=pdf4lhc_central_hist.axes[0].edges,
    label="PDF4LHC21",
    color="purple"
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 5.0))
axs.set_title("")
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
    dilepton_id_mass_pdfreweight_hist.view().value,
    edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
    label="CT09MCS",
    color="grey"
)
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
axs.stairs(
    pdf4lhc_central_hist.view().value,
    edges=pdf4lhc_central_hist.axes[0].edges,
    label="PDF4LHC21",
    color="purple"
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 5.0))
axs.set_title("")
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
    label="CT09MCS",
    color="grey"
)
axs.stairs(
    msht20lo_hist_dict["PDFMember0Weight"].view().value,
    edges=msht20lo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.stairs(
    weighthist_dict["ct18lo"].view().value,
    edges=weighthist_dict["ct18lo"].axes[0].edges,
    label="CT18LO",
    color="blue"
)
axs.stairs(
    weighthist_dict["nnpdf31lo"].view().value,
    edges=weighthist_dict["nnpdf31lo"].axes[0].edges,
    label="NNPDF31LO",
    color="green"
)
axs.stairs(
    pdf4lhc_central_hist.view().value,
    edges=pdf4lhc_central_hist.axes[0].edges,
    label="PDF4LHC21",
    color="purple"
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
    label="Mean RMS Envelope"
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 5.0))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_MeanPDF.png')
plt.close()
# PDF4LHC Replica Variation Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    dilepton_id_mass_pdfreweight_hist.view().value,
    edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
    label="CT09MCS",
    color="grey"
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
axs.stairs(
    weighthist_dict["msht20lo"].view().value,
    edges=weighthist_dict["msht20lo"].axes[0].edges,
    label="MSHT20LO",
    color="red"
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
    label="Mean RMS Envelope"
)
axs.stairs(
    pdf4lhc_central_hist.view().value,
    edges=pdf4lhc_central_hist.axes[0].edges,
    label="PDF4LHC21",
    color="purple",
    zorder=3
)
axs.bar(
    x=pdf4lhc_central_hist.axes[0].centers,
    height=(pdf4lhc_upperenv_hist.view().value - pdf4lhc_lowerenv_hist.view().value),
    bottom=pdf4lhc_lowerenv_hist.view().value,
    width=pdf4lhc_central_hist.axes[0].widths,
    linewidth=0,
    color="purple",
    alpha=0.25,
    label="PDF4LHC21 68% CL Envelope",
    zorder=2
)
axs.set_xlim((0, 300))
axs.set_ylim((0, 5.0))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('DiLeptonMassPDFVariations_PDF4LHCReplicas.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    root_file["DiLeptonMass_LOMean"] = weighthist_dict["lo_mean"]
    for weight_name in tree["pdfReweight"].fields:
        root_file[f"DiLeptonMass_{weight_name}"] = weighthist_dict[f"{weight_name}"]
    root_file["PDFReweight_Mean_Central"] = pdf_rwgt_hist_central
    root_file["PDFReweight_PDF4LHC_Central"] = pdf4lhc_rwgt_hist_central
    root_file["PDFRatio_Mean_LowerRMS"] = meanpdflo_lowerratio_hist
    root_file["PDFRatio_Mean_UpperRMS"] = meanpdflo_upperratio_hist
    root_file["PDFRatio_MSHT20LO_Central"] = msht20lo_centralratio_hist
    root_file["PDFRatio_MSHT20LO_LowerRMS"] = msht20lo_lowerratio_hist
    root_file["PDFRatio_MSHT20LO_UpperRMS"] = msht20lo_upperratio_hist
    root_file["PDFRatio_NNPDF31LO_Central"] = nnpdf31lo_centralratio_hist
    root_file["PDFRatio_NNPDF31LO_LowerRMS"] = nnpdf31lo_lowerratio_hist
    root_file["PDFRatio_NNPDF31LO_UpperRMS"] = nnpdf31lo_upperratio_hist
    root_file["DiLeptonMass_PDF4LHC_Central"] = pdf4lhc_central_hist
    root_file["DiLeptonMass_PDF4LHC_LowerEnv"] = pdf4lhc_lowerenv_hist
    root_file["DiLeptonMass_PDF4LHC_UpperEnv"] = pdf4lhc_upperenv_hist
    root_file["PDFRatio_PDF4LHC_Central"] = pdf4lhc_centralratio_hist
    root_file["PDFRatio_PDF4LHC_LowerEnv"] = pdf4lhc_lowerratio_hist
    root_file["PDFRatio_PDF4LHC_UpperEnv"] = pdf4lhc_upperratio_hist

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonpdfReweightFine": [dilepton_id_mass_pdfreweight_hist],
    "DiLeptonMass_LOMean": [weighthist_dict["lo_mean"]],
    "PDFReweight_Mean_Central": [pdf_rwgt_hist_central],
    "PDFReweight_PDF4LHC_Central": [pdf4lhc_rwgt_hist_central],
    "PDFRatio_Mean_LowerRMS": [meanpdflo_lowerratio_hist],
    "PDFRatio_Mean_UpperRMS": [meanpdflo_upperratio_hist],
    "PDFRatio_MSHT20LO_Central": [msht20lo_centralratio_hist],
    "PDFRatio_MSHT20LO_LowerRMS": [msht20lo_lowerratio_hist],
    "PDFRatio_MSHT20LO_UpperRMS": [msht20lo_upperratio_hist],
    "PDFRatio_NNPDF31LO_Central": [nnpdf31lo_centralratio_hist],
    "PDFRatio_NNPDF31LO_LowerRMS": [nnpdf31lo_lowerratio_hist],
    "PDFRatio_NNPDF31LO_UpperRMS": [nnpdf31lo_upperratio_hist],
    "DiLeptonMass_PDF4LHC_Central": [pdf4lhc_central_hist],
    "DiLeptonMass_PDF4LHC_LowerEnv": [pdf4lhc_lowerenv_hist],
    "DiLeptonMass_PDF4LHC_UpperEnv": [pdf4lhc_upperenv_hist],
    "PDFRatio_PDF4LHC_Central": [pdf4lhc_centralratio_hist],
    "PDFRatio_PDF4LHC_LowerEnv": [pdf4lhc_lowerratio_hist],
    "PDFRatio_PDF4LHC_UpperEnv": [pdf4lhc_upperratio_hist],
    "DileptonpdfReweightProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
# Store reweight histograms
# for weight_name in tree["pdfReweight"].fields:
#     pickle_dict[f"pdfReweight_{weight_name}"] = [weighthist_dict[f"{weight_name}_Reweight"]]

with open(ofile_name + ".pkl", "wb") as f:
  pickle.dump(pickle_dict, f)
