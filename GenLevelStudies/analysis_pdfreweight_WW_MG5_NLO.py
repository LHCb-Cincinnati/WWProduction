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

# Functions
def plot_pdf_rms(hist_dict, pdf_str):
    """ hist should be in the order upper, lower, central
    """

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist = hist_dict["PDFMember0Weight"]
    axs.stairs(
        hist.view().value, 
        edges=hist.axes[0].edges,
        label = pdf_str.upper() + " Central Value"
    )
    axs.bar(
        hist.axes[0].centers,
        hist.view().value,
        width = hist.axes[0].widths,
        yerr = hist.view().variance,
        fill = False,
        label = pdf_str.upper() + " RMS Envelope"
    )
    axs.set_xlim((0, 300))
    # axs.set_ylim((0, 35))
    axs.set_title("PDF Variations")
    axs.set_xlabel("$M_{e \\mu} (GeV)$")
    axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}}$")
    fig.savefig(f'DiLeptonMass_PDFVariations_{pdf_str}.png')
    plt.close()


# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section
ofile_name = args.output
luminosity = args.luminosity # Luminosity in fb^-1

# Get Tree
tree_iterator = uproot.iterate(
    [file_name + ":Tree"for file_name in file_name_list],
    step_size=45000
)

# Scale Factor
if args.cross_section == False:
    scale_factor = 1
else:
    scale_factor = args.cross_section

# Variables
unweighted_event_counter = 0
unweighted_event_counter = 0
num_nnpdf31lo_members = 101
num_msht20lo_members = 61
num_nnpdf31nlo_members = 101
num_ct18nlo_members = 59
num_msht20nlo_members = 65
num_pdf4lhc_members = 101

# Binning Scheme
bins_list = list(np.linspace(0, 100, 4)) + [150, 200, 2000] 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_pdfreweight_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_pdfreweight_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
weight_name_list = ["nnpdf31lo", "ct18lo", "msht20lo", "nnpdf31nlo", "ct18nlo", "msht20nlo", "pdf4lhc21"]
weighthist_dict = {}
for weight_name in weight_name_list:
    weighthist_dict[weight_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
nnpdf31lo_hist_dict = {}
for i in range(num_nnpdf31lo_members):
    nnpdf31lo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
msht20lo_hist_dict = {}
for i in range(num_msht20lo_members):
    msht20lo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
nnpdf31nlo_hist_dict = {}
for i in range(num_nnpdf31nlo_members):
    nnpdf31nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
ct18nlo_hist_dict = {}
for i in range(num_ct18nlo_members):
    ct18nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
msht20nlo_hist_dict = {}
for i in range(num_msht20nlo_members):
    msht20nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
pdf4lhc_hist_dict = {}
for i in range(num_pdf4lhc_members):
    pdf4lhc_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )

# Loop over files in small steps
for tree in tree_iterator:
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
    for member_name in tree["NNPDF31NLO_Members"].fields:
        nnpdf31nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["NNPDF31NLO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["CT18NLO_Members"].fields:
        ct18nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["CT18NLO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["MSHT20NLO_Members"].fields:
        msht20nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["MSHT20NLO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["PDF4LHC21_Members"].fields:
        pdf4lhc_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["PDF4LHC21_Members"][member_name][lepton_mask]*scale_factor
        )

    unweighted_event_counter += len(dilepton_vec)
    if args.debug:
        break

# Store statistical variances before rewriting
variance_stack = np.stack(
    [
        weighthist_dict["nnpdf31lo"].view(flow=True).variance,
        weighthist_dict["ct18lo"].view(flow=True).variance,
        weighthist_dict["msht20lo"].view(flow=True).variance,
    ],
    axis=0,
)
# RMS calculation for individual pdf families
nnpdf31nlo_hist_dict = at.calc_pdf_rms(nnpdf31nlo_hist_dict)
msht20nlo_hist_dict = at.calc_pdf_rms(msht20nlo_hist_dict)
ct18nlo_hist_dict = at.calc_pdf_rms(ct18nlo_hist_dict)
# Create mean NLO PDF histogram from individual families
nnpdf31nlo_view = weighthist_dict["nnpdf31nlo"].view(flow=True)
ct18nlo_view = weighthist_dict["ct18nlo"].view(flow=True)
msht20nlo_view = weighthist_dict["msht20nlo"].view(flow=True)
nlo_mean_hist = weighthist_dict["nnpdf31nlo"].copy()
nlo_mean_view = nlo_mean_hist.view(flow=True)

# Per-bin mean of the three NLO histograms
values_stack = np.stack(
    [
        nnpdf31nlo_view.value,
        ct18nlo_view.value,
        msht20nlo_view.value,
    ],
    axis=0,
)
mean_values = np.mean(values_stack, axis=0)
nlo_mean_view.value = mean_values
nlo_mean_view.variance = np.sqrt(np.mean(np.abs(values_stack-mean_values)**2, axis=0))
weighthist_dict["nlo_mean"] = nlo_mean_hist
# Create RMS NLO PDF histogram
nlo_rmslow_hist = weighthist_dict["nlo_mean"].copy()
nlo_rmslow_hist.view(flow=True).value = nlo_rmslow_hist.view(flow=True).value - nlo_rmslow_hist.view(flow=True).variance
nlo_rmshigh_hist = weighthist_dict["nlo_mean"].copy()
nlo_rmshigh_hist.view(flow=True).value = nlo_rmshigh_hist.view(flow=True).value + nlo_rmshigh_hist.view(flow=True).variance
# Create NNPDF31LO RMS NLO PDF histogram
nnpdf31nlo_rmslow_hist = nnpdf31nlo_hist_dict["PDFMember0Weight"].copy()
nnpdf31nlo_rmslow_hist.view(flow=True).value = nnpdf31nlo_rmslow_hist.view(flow=True).value - nnpdf31nlo_rmslow_hist.view(flow=True).variance
nnpdf31nlo_rmslow_hist.view(flow=True).variance = weighthist_dict["nnpdf31nlo"].view(flow=True).variance
nnpdf31nlo_rmshigh_hist = nnpdf31nlo_hist_dict["PDFMember0Weight"].copy()
nnpdf31nlo_rmshigh_hist.view(flow=True).value = nnpdf31nlo_rmshigh_hist.view(flow=True).value + nnpdf31nlo_rmshigh_hist.view(flow=True).variance
nnpdf31nlo_rmshigh_hist.view(flow=True).variance = weighthist_dict["nnpdf31nlo"].view(flow=True).variance
nnpdf31nlo_rmslow_hist.view(flow=True).value[-1] = 0.0
# Create MSHT20LO RMS NLO PDF histogram
msht20nlo_rmslow_hist = msht20nlo_hist_dict["PDFMember0Weight"].copy()
msht20nlo_rmslow_hist.view(flow=True).value = msht20nlo_rmslow_hist.view(flow=True).value - msht20nlo_rmslow_hist.view(flow=True).variance
msht20nlo_rmslow_hist.view(flow=True).variance = weighthist_dict["msht20nlo"].view(flow=True).variance
msht20nlo_rmshigh_hist = msht20nlo_hist_dict["PDFMember0Weight"].copy()
msht20nlo_rmshigh_hist.view(flow=True).value = msht20nlo_rmshigh_hist.view(flow=True).value + msht20nlo_rmshigh_hist.view(flow=True).variance
msht20nlo_rmshigh_hist.view(flow=True).variance = weighthist_dict["msht20nlo"].view(flow=True).variance
msht20nlo_rmslow_hist.view(flow=True).value[-1] = 0.0
# Create MSHT20LO RMS NLO PDF histogram
ct18nlo_rmslow_hist = ct18nlo_hist_dict["PDFMember0Weight"].copy()
ct18nlo_rmslow_hist.view(flow=True).value = ct18nlo_rmslow_hist.view(flow=True).value - ct18nlo_rmslow_hist.view(flow=True).variance
ct18nlo_rmslow_hist.view(flow=True).variance = weighthist_dict["ct18nlo"].view(flow=True).variance
ct18nlo_rmshigh_hist = ct18nlo_hist_dict["PDFMember0Weight"].copy()
ct18nlo_rmshigh_hist.view(flow=True).value = ct18nlo_rmshigh_hist.view(flow=True).value + ct18nlo_rmshigh_hist.view(flow=True).variance
ct18nlo_rmshigh_hist.view(flow=True).variance = weighthist_dict["ct18nlo"].view(flow=True).variance
ct18nlo_rmslow_hist.view(flow=True).value[-1] = 0.0
pdf4lhc_central_hist, pdf4lhc_lowerenv_hist, pdf4lhc_upperenv_hist = at.calc_pdf_mc_envelope(
    pdf4lhc_hist_dict
)

# Reweight Histograms
nlo_mean_hist_wstaterrs = weighthist_dict["nlo_mean"].copy()
nlo_mean_hist_wstaterrs.view(flow=True).variance = np.sum(variance_stack) / 9
pdf_rwgt_hist_central = at.divide_bh_histograms(
    nlo_mean_hist_wstaterrs, 
    dilepton_id_mass_pdfreweight_hist
)
pdf_rwgt_hist_lowerRMS = weighthist_dict["nlo_mean"].copy()
pdf_rwgt_hist_lowerRMS.view(flow=True).value = (
    pdf_rwgt_hist_lowerRMS.view(flow=True).value
    - pdf_rwgt_hist_lowerRMS.view(flow=True).variance
)
pdf_rwgt_hist_lowerRMS = at.divide_bh_histograms(
    pdf_rwgt_hist_lowerRMS, 
    dilepton_id_mass_pdfreweight_hist, 
    error_type="pass_through"
)
pdf_rwgt_hist_upperRMS = weighthist_dict["nlo_mean"].copy()
pdf_rwgt_hist_upperRMS.view(flow=True).value = (
    pdf_rwgt_hist_upperRMS.view(flow=True).value
    + pdf_rwgt_hist_upperRMS.view(flow=True).variance
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

# Divide Hists
pdfrwgt_hist = at.divide_bh_histograms(weighthist_dict["nlo_mean"], dilepton_id_mass_pdfreweight_hist)
pdfrwgt_hist.view(flow=True).variance = (
    weighthist_dict["nlo_mean"].view(flow=True).variance 
    / dilepton_id_mass_pdfreweight_hist.view(flow=True).value
)

# Print Statements:
print(f"Unweighted Events: {unweighted_event_counter}")
print(f"Weighted Events: {unweighted_event_counter * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name].view(flow=True).value.sum() / dilepton_id_mass_pdfreweight_hist.view(flow=True).value.sum() * unweighted_event_counter * scale_factor}")
print(f"Mean PDF: {weighthist_dict['nlo_mean'].view(flow=True).value.sum() / dilepton_id_mass_pdfreweight_hist.view(flow=True).value.sum() * unweighted_event_counter * scale_factor}")

if args.debug:
    pass
    # exit()

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
at.create_stair(pdfrwgt_hist, "DiLepton Mass Linear PDF RWT Hist",
                luminosity=luminosity)

# Special pdf plots
for weight_name in tree["pdfReweight"].fields:
    at.create_stacked_stair(
        [weighthist_dict[weight_name], dilepton_id_mass_pdfreweight_hist],
        f"DiLeptonMasspdfReweightBinning_{weight_name}",
        [weight_name, "CT09MCS"]
    )
# NNPDF31LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(nnpdf31lo_hist_dict),
    "NNPD31LO"
)
# MSHT20LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(msht20lo_hist_dict), "MSHT20LO"
)
# NNPDF31NLO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(nnpdf31nlo_hist_dict),
    "NNPD31NLO"
)
# CT18NLO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(ct18nlo_hist_dict),
    "CT18NLO"
)
# MSHT20LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(msht20nlo_hist_dict),
    "MSHT20NLO"
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
    dilepton_id_mass_pdfreweight_profilehist.view().value,
    pdf_rwgt_hist_central.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = pdf_rwgt_hist_central.view().variance,
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
meanpdfnlo_lowerratio_hist, meanpdfnlo_centralratio_hist, meanpdfnlo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["nlo_mean"], nlo_rmslow_hist, nlo_rmshigh_hist
) 
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=meanpdfnlo_centralratio_hist.axes[0].centers, 
    height=(meanpdfnlo_upperratio_hist.view().value - meanpdfnlo_lowerratio_hist.view().value), 
    bottom=meanpdfnlo_lowerratio_hist.view().value, 
    width=meanpdfnlo_centralratio_hist.axes[0].widths, 
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
# axs.set_ylim((0, 4.0))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_MeanPDF.png')
plt.close()
# NNPDF31NLO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
nnpdf31nlo_lowerratio_hist, nnpdf31nlo_centralratio_hist, nnpdf31nlo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["nnpdf31nlo"], nnpdf31nlo_rmslow_hist, nnpdf31nlo_rmshigh_hist
) 
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=nnpdf31nlo_centralratio_hist.axes[0].centers, 
    height=(nnpdf31nlo_upperratio_hist.view().value - nnpdf31nlo_lowerratio_hist.view().value), 
    bottom=nnpdf31nlo_lowerratio_hist.view().value, 
    width=nnpdf31nlo_centralratio_hist.axes[0].widths, 
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
# axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_NNPDF31LO.png')
plt.close()
# MSHT20NLO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
msht20nlo_lowerratio_hist, msht20nlo_centralratio_hist, msht20nlo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["msht20nlo"], msht20nlo_rmslow_hist, msht20nlo_rmshigh_hist
) 
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=msht20nlo_centralratio_hist.axes[0].centers, 
    height=(msht20nlo_upperratio_hist.view().value - msht20nlo_lowerratio_hist.view().value), 
    bottom=msht20nlo_lowerratio_hist.view().value, 
    width=msht20nlo_centralratio_hist.axes[0].widths, 
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
# axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_ct18LO.png')
plt.close()
# CT18NLO Ratio Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
ct18nlo_lowerratio_hist, ct18nlo_centralratio_hist, ct18nlo_upperratio_hist = at.calc_rms_ratio_hists(
    weighthist_dict["ct18nlo"], ct18nlo_rmslow_hist, ct18nlo_rmshigh_hist
) 
pdf4lhc_lowerratio_hist, pdf4lhc_centralratio_hist, pdf4lhc_upperratio_hist = at.calc_rms_ratio_hists(
    pdf4lhc_central_hist, pdf4lhc_lowerenv_hist, pdf4lhc_upperenv_hist
)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=ct18nlo_centralratio_hist.axes[0].centers, 
    height=(ct18nlo_upperratio_hist.view().value - ct18nlo_lowerratio_hist.view().value), 
    bottom=ct18nlo_lowerratio_hist.view().value, 
    width=ct18nlo_centralratio_hist.axes[0].widths, 
    linewidth=0, 
    color="red", 
    alpha=0.25, 
    label="RMS Deviation"
)
# axs.stairs(
#     ct18lo_central_hist.view().value, 
#     edges=ct18lo_central_hist.axes[0].edges,
#     label="Central Value",
#     color="black",
#     zorder=3
# )
# axs.stairs(
#     ct18lo_lowerRMS_hist.view().value, 
#     edges=ct18lo_lowerRMS_hist.axes[0].edges,
#     label="Lower Envelope",
#     color="black"
# )
# axs.stairs(
#     ct18lo_upperRMS_hist.view().value, 
#     edges=ct18lo_upperRMS_hist.axes[0].edges,
#     label="Upper Envelope",
#     color="black"
# )
# axs.bar(
#     x=ct18lo_upperRMS_hist.axes[0].centers, 
#     height=(ct18lo_upperRMS_hist.view().value - ct18lo_lowerRMS_hist.view().value), 
#     bottom=ct18lo_upperRMS_hist.view().value, 
#     width=ct18lo_upperRMS_hist.axes[0].widths, 
#     # align='edge', 
#     linewidth=0, 
#     color="red", 
#     alpha=0.25, 
#     # zorder=-1, 
#     label="NNPDF31NLO Envelope"
# )
# axs.set_ylim((0, 4.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_Ratio_CT18NLO.png')
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
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("PDF Deviation to the Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('DiLeptonMassPDFVariations_Ratio_PDF4LHC.png')
plt.close()
# Many Family PDF Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
# Plot NNPDF31NLO
axs.stairs(
    weighthist_dict["nnpdf31nlo"].view().value, 
    edges=weighthist_dict["nnpdf31nlo"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
# Plot CT18NLO
axs.stairs(
    weighthist_dict["ct18nlo"].view().value, 
    edges=weighthist_dict["ct18nlo"].axes[0].edges,
    label="CT18NLO",
    color="blue"
)
# Plot MSHT20NLO
axs.stairs(
    weighthist_dict["msht20nlo"].view().value, 
    edges=weighthist_dict["msht20nlo"].axes[0].edges,
    label="MSHT20NLO",
    color="red"
)
# Plot NLO Mean
axs.stairs(
    weighthist_dict["nlo_mean"].view().value, 
    edges=weighthist_dict["nlo_mean"].axes[0].edges,
    color="black"
)
axs.stairs(
    pdf4lhc_central_hist.view().value,
    edges=pdf4lhc_central_hist.axes[0].edges,
    label="PDF4LHC21",
    color="purple"
)

axs.bar(
    x=weighthist_dict["nlo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["nlo_mean"].view().variance, 
    bottom=weighthist_dict["nlo_mean"].view().value - weighthist_dict["nlo_mean"].view().variance, 
    width=weighthist_dict["nlo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    zorder=-1, 
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
# NNPDF31NLO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    weighthist_dict["nlo_mean"].view().value, 
    edges=weighthist_dict["nlo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black",
    zorder=3
)
axs.bar(
    x=weighthist_dict["nlo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["nlo_mean"].view().variance, 
    bottom=weighthist_dict["nlo_mean"].view().value - weighthist_dict["nlo_mean"].view().variance, 
    width=weighthist_dict["nlo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean RMS Envelope",
    zorder=2
)
axs.stairs(
    nnpdf31nlo_hist_dict["PDFMember0Weight"].view().value, 
    edges=nnpdf31nlo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
axs.bar(
    x=nnpdf31nlo_hist_dict["PDFMember0Weight"].axes[0].centers, 
    height= 2*nnpdf31nlo_hist_dict["PDFMember0Weight"].view().variance, 
    bottom=(
        nnpdf31nlo_hist_dict["PDFMember0Weight"].view().value 
        - nnpdf31nlo_hist_dict["PDFMember0Weight"].view().variance
    ), 
    width=nnpdf31nlo_hist_dict["PDFMember0Weight"].axes[0].widths, 
    linewidth=0, 
    color="green", 
    alpha=0.25, 
    label="NNPDF31NLO RMS Envelope"
)
axs.stairs(
    weighthist_dict["msht20nlo"].view().value, 
    edges=weighthist_dict["msht20nlo"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.stairs(
    weighthist_dict["ct18nlo"].view().value, 
    edges=weighthist_dict["ct18nlo"].axes[0].edges,
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
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_NNPDF31NLO.png')
plt.close()
# CT18NLO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    weighthist_dict["nlo_mean"].view().value, 
    edges=weighthist_dict["nlo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black",
    zorder=3
)
axs.bar(
    x=weighthist_dict["nlo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["nlo_mean"].view().variance, 
    bottom=weighthist_dict["nlo_mean"].view().value - weighthist_dict["nlo_mean"].view().variance, 
    width=weighthist_dict["nlo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean RMS Envelope",
    zorder=2
)
axs.stairs(
    nnpdf31nlo_hist_dict["PDFMember0Weight"].view().value, 
    edges=nnpdf31nlo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
axs.bar(
    x=ct18nlo_hist_dict["PDFMember0Weight"].axes[0].centers, 
    height= 2*ct18nlo_hist_dict["PDFMember0Weight"].view().variance, 
    bottom=(
        ct18nlo_hist_dict["PDFMember0Weight"].view().value 
        - ct18nlo_hist_dict["PDFMember0Weight"].view().variance
    ), 
    width=ct18nlo_hist_dict["PDFMember0Weight"].axes[0].widths, 
    linewidth=0, 
    color="blue", 
    alpha=0.25, 
    label="CT18NLO RMS Envelope"
)
axs.stairs(
    weighthist_dict["msht20nlo"].view().value, 
    edges=weighthist_dict["msht20nlo"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.stairs(
    weighthist_dict["ct18nlo"].view().value, 
    edges=weighthist_dict["ct18nlo"].axes[0].edges,
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
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_CT18NLO.png')
plt.close()
# MSHT20NLO RMS Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    weighthist_dict["nlo_mean"].view().value, 
    edges=weighthist_dict["nlo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black",
    zorder=3
)
axs.bar(
    x=weighthist_dict["nlo_mean"].axes[0].centers, 
    height= 2*weighthist_dict["nlo_mean"].view().variance, 
    bottom=weighthist_dict["nlo_mean"].view().value - weighthist_dict["nlo_mean"].view().variance, 
    width=weighthist_dict["nlo_mean"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    label="Mean RMS Envelope",
    zorder=2
)
axs.stairs(
    msht20nlo_hist_dict["PDFMember0Weight"].view().value, 
    edges=nnpdf31nlo_hist_dict["PDFMember0Weight"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
axs.bar(
    x=msht20nlo_hist_dict["PDFMember0Weight"].axes[0].centers, 
    height= 2*msht20nlo_hist_dict["PDFMember0Weight"].view().variance, 
    bottom=(
        msht20nlo_hist_dict["PDFMember0Weight"].view().value 
        - msht20nlo_hist_dict["PDFMember0Weight"].view().variance
    ), 
    width=msht20nlo_hist_dict["PDFMember0Weight"].axes[0].widths, 
    linewidth=0, 
    color="red", 
    alpha=0.25, 
    label="MSHT20NLo RMS Envelope"
)
axs.stairs(
    weighthist_dict["msht20nlo"].view().value, 
    edges=weighthist_dict["msht20nlo"].axes[0].edges,
    label="MSHT20LO",
    color="red"
)
axs.stairs(
    weighthist_dict["ct18nlo"].view().value, 
    edges=weighthist_dict["ct18nlo"].axes[0].edges,
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
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations_MSHT20NLO.png')
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
    weighthist_dict["nnpdf31nlo"].view().value,
    edges=weighthist_dict["nnpdf31nlo"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
axs.stairs(
    weighthist_dict["ct18nlo"].view().value,
    edges=weighthist_dict["ct18nlo"].axes[0].edges,
    label="CT18NLO",
    color="blue"
)
axs.stairs(
    weighthist_dict["msht20nlo"].view().value,
    edges=weighthist_dict["msht20nlo"].axes[0].edges,
    label="MSHT20NLO",
    color="red"
)
axs.stairs(
    weighthist_dict["nlo_mean"].view().value,
    edges=weighthist_dict["nlo_mean"].axes[0].edges,
    label="Mean Central Value",
    color="black"
)
axs.bar(
    x=weighthist_dict["nlo_mean"].axes[0].centers,
    height= 2*weighthist_dict["nlo_mean"].view().variance,
    bottom=weighthist_dict["nlo_mean"].view().value - weighthist_dict["nlo_mean"].view().variance,
    width=weighthist_dict["nlo_mean"].axes[0].widths,
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
    root_file["DiLeptonMass_LOMean"] = weighthist_dict["nlo_mean"]
    for weight_name in tree["pdfReweight"].fields:
        root_file[f"DiLeptonMass_{weight_name}"] = weighthist_dict[f"{weight_name}"]
    root_file["PDFReweight_Mean_Central"] = pdf_rwgt_hist_central
    root_file["PDFReweight_PDF4LHC_Central"] = pdf4lhc_rwgt_hist_central
    root_file["PDFRatio_Mean_LowerRMS"] = meanpdfnlo_lowerratio_hist
    root_file["PDFRatio_Mean_UpperRMS"] = meanpdfnlo_upperratio_hist
    root_file["PDFRatio_MSHT20NLO_Central"] = msht20nlo_centralratio_hist
    root_file["PDFRatio_MSHT20NLO_LowerRMS"] = msht20nlo_lowerratio_hist
    root_file["PDFRatio_MSHT20NLO_UpperRMS"] = msht20nlo_upperratio_hist
    root_file["PDFRatio_NNPDF31NLO_Central"] = nnpdf31nlo_centralratio_hist
    root_file["PDFRatio_NNPDF31NLO_LowerRMS"] = nnpdf31nlo_lowerratio_hist
    root_file["PDFRatio_NNPDF31NLO_UpperRMS"] = nnpdf31nlo_upperratio_hist
    root_file["PDFRatio_CT18NLO_Central"] = ct18nlo_centralratio_hist
    root_file["PDFRatio_CT18NLO_LowerRMS"] = ct18nlo_lowerratio_hist
    root_file["PDFRatio_CT18NLO_UpperRMS"] = ct18nlo_upperratio_hist
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
    "DiLeptonMass_LOMean": [weighthist_dict["nlo_mean"]],
    "PDFReweight_Mean_Central": [pdf_rwgt_hist_central],
    "PDFReweight_PDF4LHC_Central": [pdf4lhc_rwgt_hist_central],
    "PDFRatio_Mean_LowerRMS": [meanpdfnlo_lowerratio_hist],
    "PDFRatio_Mean_UpperRMS": [meanpdfnlo_upperratio_hist],
    "PDFRatio_MSHT20NLO_Central": [msht20nlo_centralratio_hist],
    "PDFRatio_MSHT20NLO_LowerRMS": [msht20nlo_lowerratio_hist],
    "PDFRatio_MSHT20NLO_UpperRMS": [msht20nlo_upperratio_hist],
    "PDFRatio_NNPDF31NLO_Central": [nnpdf31nlo_centralratio_hist],
    "PDFRatio_NNPDF31NLO_LowerRMS": [nnpdf31nlo_lowerratio_hist],
    "PDFRatio_NNPDF31NLO_UpperRMS": [nnpdf31nlo_upperratio_hist],
    "PDFRatio_CT18NLO_Central": [ct18nlo_centralratio_hist],
    "PDFRatio_CT18NLO_LowerRMS": [ct18nlo_lowerratio_hist],
    "PDFRatio_CT18NLO_UpperRMS": [ct18nlo_upperratio_hist],
    "DiLeptonMass_PDF4LHC_Central": [pdf4lhc_central_hist],
    "DiLeptonMass_PDF4LHC_LowerEnv": [pdf4lhc_lowerenv_hist],
    "DiLeptonMass_PDF4LHC_UpperEnv": [pdf4lhc_upperenv_hist],
    "PDFRatio_PDF4LHC_Central": [pdf4lhc_centralratio_hist],
    "PDFRatio_PDF4LHC_LowerEnv": [pdf4lhc_lowerratio_hist],
    "PDFRatio_PDF4LHC_UpperEnv": [pdf4lhc_upperratio_hist],
    "DileptonpdfReweightProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
