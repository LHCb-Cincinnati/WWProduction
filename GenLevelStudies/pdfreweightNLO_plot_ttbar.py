""" Program to compute the pdf variations to the K-factor reweighting for ttbar process.

Program to compute the K-factor reweighting for ttbar and DFDY process.  
The input for this program should be the Lower RMS Histogram 
histogram, central Mean PDF, and the Upper RMS histogram pickle 
files in that order.
"""

# Imports
# Standard Library Imports
import os
import pdb
import pickle
import sys
import logging
import array
from dataclasses import dataclass
# Scikit Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
# HEP Imports
import ROOT
from hepunits.units import MeV, GeV
import awkward as ak
import uproot
import boost_histogram as bh
import iminuit
from iminuit import minimize  # has same interface as scipy.optimize.minimize
from iminuit import Minuit, describe
from iminuit.cost import LeastSquares
# Personal Imports
import AnalysisTools as at

# Exponental Decay function
def exp_decay(x, a, b, c):
    return(a*np.exp(c*x) + b)

# Parse inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
cross_section = args.cross_section # Cross section in fb
luminosity = args.luminosity # Luminosity in fb^-1
file_list = [file.name for file in args.input_files]
ofile_name = args.output
logging.info(f"Arguments: {args}")

# Define Variables

# Retrieve Histos
# Here data are the histograms and index is the index of the file in the
# list of files. So data_list looks like:
# data_list[1] = dictionary first file's histograms
# The dictionary looks like
# dict[Histogram Tag] = Histogram
with open(file_list[0], 'rb') as f:
    hist_dict = pickle.load(f)
# Grab Histograms from File.
profile_hist = hist_dict["DiLeptonMassProfile"][0]
meanpdf_lowerRMS_hist = hist_dict["DileptonKFactorFine_MeanPDF_LowerRMS"][0]
meanpdf_central_hist = hist_dict["DileptonKFactorFine_MeanPDF_Central"][0]
meanpdf_upperRMS_hist = hist_dict["DileptonKFactorFine_MeanPDF_UpperRMS"][0]
pdf4lhc21_lowerRMS_hist = hist_dict["DileptonKFactorFine_PDF4LHC21_LowerRMS"][0]
pdf4lhc21_central_hist = hist_dict["DileptonKFactorFine_PDF4LHC21_Central"][0]
pdf4lhc21_upperRMS_hist = hist_dict["DileptonKFactorFine_PDF4LHC21_UpperRMS"][0]
nnpdf31nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_LowerRMS"][0]
nnpdf31nlo_central_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_Central"][0]
nnpdf31nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_UpperRMS"][0]
msht20nlo_central_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_Central"][0]
msht20nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_LowerRMS"][0]
msht20nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_UpperRMS"][0]
ct18nlo_central_hist = hist_dict["DileptonKFactorFine_CT18NLO_Central"][0]
ct18nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_CT18NLO_LowerRMS"][0]
ct18nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_CT18NLO_UpperRMS"][0]

# Create ratio histograms
# Large final bin ratio to 2 due to large errors.
meanpdf_lowerratio_hist, meanpdf_centralratio_hist, meanpdf_upperratio_hist = at.calc_rms_ratio_hists(
    meanpdf_central_hist, meanpdf_lowerRMS_hist, meanpdf_upperRMS_hist
) 
nnpdf31nlo_lowerratio_hist, nnpdf31nlo_centralratio_hist, nnpdf31nlo_upperratio_hist = at.calc_rms_ratio_hists(
    nnpdf31nlo_central_hist, nnpdf31nlo_lowerRMS_hist, nnpdf31nlo_upperRMS_hist
) 
msht20nlo_lowerratio_hist, msht20nlo_centralratio_hist, msht20nlo_upperratio_hist = at.calc_rms_ratio_hists(
    msht20nlo_central_hist, msht20nlo_lowerRMS_hist, msht20nlo_upperRMS_hist
)
ct18nlo_lowerratio_hist, ct18nlo_centralratio_hist, ct18nlo_upperratio_hist = at.calc_rms_ratio_hists(
    ct18nlo_central_hist, ct18nlo_lowerRMS_hist, ct18nlo_upperRMS_hist
) 
pdf4lhc21_lowerratio_hist, pdf4lhc21_centralratio_hist, pdf4lhc21_upperratio_hist = at.calc_rms_ratio_hists(
    pdf4lhc21_central_hist, pdf4lhc21_lowerRMS_hist, pdf4lhc21_upperRMS_hist
) 

# Fit K-Factor
fit_func = exp_decay

# # Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# K-Factor with scale variations of PDF families
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value"
)
axs.bar(
    x=meanpdf_central_hist.axes[0].centers, 
    height=2*(meanpdf_upperRMS_hist.view().value - meanpdf_central_hist.view().value), 
    bottom=meanpdf_lowerRMS_hist.view().value, 
    width=meanpdf_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    # zorder=-1, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.stairs(
    pdf4lhc21_central_hist.view().value, 
    edges=pdf4lhc21_central_hist.axes[0].edges,
    color="purple",
    label="PDF4LHC21"
)
# Setting axes and legend
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 3.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("K-Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFFamilies.png')
plt.close()
# Specific PDF comparison plots
# NNPDF31NLO
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=meanpdf_central_hist.axes[0].centers, 
    height=2*(meanpdf_upperRMS_hist.view().value - meanpdf_central_hist.view().value), 
    bottom=meanpdf_lowerRMS_hist.view().value, 
    width=meanpdf_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.bar(
    x=nnpdf31nlo_central_hist.axes[0].centers, 
    height=(nnpdf31nlo_upperRMS_hist.view().value - nnpdf31nlo_lowerRMS_hist.view().value), 
    bottom=nnpdf31nlo_lowerRMS_hist.view().value, 
    width=nnpdf31nlo_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="green", 
    alpha=0.25, 
    # zorder=-1, 
    label="NNPDF31NLO RMS Envelope"
)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.stairs(
    pdf4lhc21_central_hist.view().value, 
    edges=pdf4lhc21_central_hist.axes[0].edges,
    color="purple",
    label="PDF4LHC21"
)
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 3.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("K-Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_NNPDF31NLO.png')
plt.close()
# MSHT20NLO
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=meanpdf_central_hist.axes[0].centers, 
    height=2*(meanpdf_upperRMS_hist.view().value - meanpdf_central_hist.view().value), 
    bottom=meanpdf_lowerRMS_hist.view().value, 
    width=meanpdf_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.bar(
    x=msht20nlo_central_hist.axes[0].centers, 
    height=(msht20nlo_upperRMS_hist.view().value - msht20nlo_lowerRMS_hist.view().value), 
    bottom=msht20nlo_lowerRMS_hist.view().value, 
    width=msht20nlo_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="red", 
    alpha=0.25, 
    # zorder=-1, 
    label="MSHT20NLO RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.stairs(
    pdf4lhc21_central_hist.view().value, 
    edges=pdf4lhc21_central_hist.axes[0].edges,
    color="purple",
    label="PDF4LHC21"
)
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 3.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("K-Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_MSHT20NLO.png')
plt.close()
# CT18NLO
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=meanpdf_central_hist.axes[0].centers, 
    height=2*(meanpdf_upperRMS_hist.view().value - meanpdf_central_hist.view().value), 
    bottom=meanpdf_lowerRMS_hist.view().value, 
    width=meanpdf_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.bar(
    x=ct18nlo_central_hist.axes[0].centers, 
    height=(ct18nlo_upperRMS_hist.view().value - ct18nlo_lowerRMS_hist.view().value), 
    bottom=ct18nlo_lowerRMS_hist.view().value, 
    width=ct18nlo_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="blue", 
    alpha=0.25, 
    # zorder=-1, 
    label="CT18NLO RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    pdf4lhc21_central_hist.view().value, 
    edges=pdf4lhc21_central_hist.axes[0].edges,
    color="purple",
    label="PDF4LHC21"
)
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 3.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("K-Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_CT18NLO.png')
plt.close()
# PDF4LHC21
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=meanpdf_central_hist.axes[0].centers, 
    height=2*(meanpdf_upperRMS_hist.view().value - meanpdf_central_hist.view().value), 
    bottom=meanpdf_lowerRMS_hist.view().value, 
    width=meanpdf_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    pdf4lhc21_central_hist.view().value, 
    edges=pdf4lhc21_central_hist.axes[0].edges,
    color="purple",
    label="PDF4LHC21"
)
axs.bar(
    x=pdf4lhc21_central_hist.axes[0].centers, 
    height=(pdf4lhc21_upperRMS_hist.view().value - pdf4lhc21_lowerRMS_hist.view().value), 
    bottom=pdf4lhc21_lowerRMS_hist.view().value, 
    width=pdf4lhc21_central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="purple", 
    alpha=0.25, 
    # zorder=-1, 
    label="PDF4LHC21 68% CL Envelope"
)
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 3.0))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("K-Factor")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_PDF4LHC21.png')
plt.close()
# Plot reweight histogram with errors
# Mean PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
# axs.errorbar(
#     profile_hist.view().value,
#     meanpdf_centralratio_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     # yerr = (meanpdf_upperratio_hist.view().value - meanpdf_lowerratio_hist.view().value) / 2,
#     label="RMS Deviation"
# )
axs.bar(
    x=meanpdf_centralratio_hist.axes[0].centers, 
    height=(meanpdf_upperratio_hist.view().value - meanpdf_lowerratio_hist.view().value), 
    bottom=meanpdf_lowerratio_hist.view().value, 
    width=meanpdf_centralratio_hist.axes[0].widths, 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    label="RMS Deviation"
)
# axs.scatter(
#     profile_hist.view().value,
#     meanpdf_upperratio_hist.view().value,
#     color = "blue",
# #     label="Upper RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     meanpdf_upperratio_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(meanpdf_centralratio_hist.view().variance),
#     label="Upper RMS"
# )
# axs.scatter(
#     profile_hist.view().value,
#     meanpdf_lowerratio_hist.view().value,
#     color = "red",
# #     label="Lower RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     meanpdf_lowerratio_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(meanpdf_centralratio_hist.view().variance),
#     label="Lower RMS"
# )
# Setting axes and legend
# ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0.5, 1.5))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$\\frac{\\sigma(K)}{K}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {meanpdf_minuit['Central'].fval:.1f} / {meanpdf_minuit['Central'].ndof:.0f} = {meanpdf_minuit['Central'].fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_MeanPDF.png')
plt.close()

# NNPDF31NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
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
# axs.errorbar(
#     profile_hist.view().value,
#     nnpdf31nlo_centralratio_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     # yerr = (nnpdf31nlo_upperratio_hist.view().value - nnpdf31nlo_lowerratio_hist.view().value) / 2,
#     label="RMS Deviation"
# )
# axs.scatter(
#     profile_hist.view().value,
#     nnpdf31nlo_upperratio_hist.view().value,
#     color = "blue",
# #     label="Upper RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     nnpdf31nlo_upperratio_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(nnpdf31nlo_upperRMS_hist.view().variance),
#     label="Upper RMS"
# )
# axs.scatter(
#     profile_hist.view().value,
#     nnpdf31nlo_lowerratio_hist.view().value,
#     color = "red",
# #     label="Lower RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     nnpdf31nlo_lowerratio_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(nnpdf31nlo_lowerRMS_hist.view().variance),
#     label="Lower RMS"
# )
# Setting axes and legend
# ymax = max(nnpdf31nlo_upperRMS_hist.view().value)
# axs.set_ylim((0, 3.0))
axs.set_ylim((0.5, 1.5))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$\\frac{\\sigma(K)}{K}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {nnpdf31nlo_minuit['Central'].fval:.1f} / {nnpdf31nlo_minuit['Central'].ndof:.0f} = {nnpdf31nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_NNPDF31NLO.png')
plt.close()

# MSHT20NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
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
# axs.errorbar(
#     profile_hist.view().value,
#     msht20nlo_centralratio_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     # yerr = (msht20nlo_upperratio_hist.view().value - msht20nlo_lowerratio_hist.view().value) / 2,
#     label="RMS Deviation"
# )
# axs.scatter(
#     profile_hist.view().value,
#     msht20nlo_upperratio_hist.view().value,
#     color = "blue",
# #     label="Upper RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     msht20nlo_upperratio_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(msht20nlo_upperRMS_hist.view().variance),
#     label="Upper RMS"
# )
# axs.scatter(
#     profile_hist.view().value,
#     msht20nlo_lowerratio_hist.view().value,
#     color = "red",
# #     label="Lower RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     msht20nlo_lowerratio_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(msht20nlo_lowerRMS_hist.view().variance),
#     label="Lower RMS"
# )
# Setting axes and legend
# ymax = max(msht20nlo_upperRMS_hist.view().value)
axs.set_ylim((0.5, 1.5))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$\\frac{\\sigma(K)}{K}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {msht20nlo_minuit['Central'].fval:.1f} / {msht20nlo_minuit['Central'].ndof:.0f} = {msht20nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_MSHT20NLO.png')
plt.close()

# CT18NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=ct18nlo_centralratio_hist.axes[0].centers, 
    height=(ct18nlo_upperratio_hist.view().value - ct18nlo_lowerratio_hist.view().value), 
    bottom=ct18nlo_lowerratio_hist.view().value, 
    width=ct18nlo_centralratio_hist.axes[0].widths, 
    linewidth=0, 
    color="blue", 
    alpha=0.25, 
    label="RMS Deviation"
)
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_centralratio_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     # yerr = (ct18nlo_upperratio_hist.view().value - ct18nlo_lowerratio_hist.view().value) / 2,
#     label="RMS Deviation"
# )
# axs.scatter(
#     profile_hist.view().value,
#     ct18nlo_upperratio_hist.view().value,
#     color = "blue",
# #     label="Upper RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_upperratio_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(ct18nlo_upperRMS_hist.view().variance),
#     label="Upper RMS"
# )
# axs.scatter(
#     profile_hist.view().value,
#     ct18nlo_lowerratio_hist.view().value,
#     color = "red",
# #     label="Lower RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_lowerratio_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(ct18nlo_lowerRMS_hist.view().variance),
#     label="Lower RMS"
# )
# Setting axes and legend
# ymax = max(ct18nlo_upperRMS_hist.view().value)
axs.set_ylim((0.5, 1.5))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$\\frac{\\sigma(K)}{K}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {ct18nlo_minuit['Central'].fval:.1f} / {ct18nlo_minuit['Central'].ndof:.0f} = {ct18nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_CT18NLO.png')
plt.close()

# PDF4LHC21 PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
axs.bar(
    x=pdf4lhc21_centralratio_hist.axes[0].centers, 
    height=(pdf4lhc21_upperratio_hist.view().value - pdf4lhc21_lowerratio_hist.view().value), 
    bottom=pdf4lhc21_lowerratio_hist.view().value, 
    width=pdf4lhc21_centralratio_hist.axes[0].widths, 
    linewidth=0, 
    color="purple", 
    alpha=0.25, 
    label="68% CL"
)
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_centralratio_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     # yerr = (ct18nlo_upperratio_hist.view().value - ct18nlo_lowerratio_hist.view().value) / 2,
#     label="RMS Deviation"
# )
# axs.scatter(
#     profile_hist.view().value,
#     ct18nlo_upperratio_hist.view().value,
#     color = "blue",
# #     label="Upper RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_upperratio_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(ct18nlo_upperRMS_hist.view().variance),
#     label="Upper RMS"
# )
# axs.scatter(
#     profile_hist.view().value,
#     ct18nlo_lowerratio_hist.view().value,
#     color = "red",
# #     label="Lower RMS"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_lowerratio_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(ct18nlo_lowerRMS_hist.view().variance),
#     label="Lower RMS"
# )
# Setting axes and legend
# ymax = max(ct18nlo_upperRMS_hist.view().value)
axs.set_ylim((0.5, 1.5))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$\\frac{\\sigma(K)}{K}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
# bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
# text_x = bbox.x0 + 0.1
# text_y = bbox.y1 + 0.05  # small offset below the legend
# axs.text(
#     text_x, text_y,
#     f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {ct18nlo_minuit['Central'].fval:.1f} / {ct18nlo_minuit['Central'].ndof:.0f} = {ct18nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
#     transform=axs.transAxes,
#     fontsize=axs.legend().get_texts()[0].get_fontsize(),
#     verticalalignment='top',
#     # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
# )
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_PDF4LHC21.png')
plt.close()

# Create ROOT histogram
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["KFactor_rwgt_MeanPDF_Central"] = meanpdf_central_hist 
    root_file["KFactor_rwgt_MeanPDF_lowerRMS"] = meanpdf_lowerRMS_hist
    root_file["KFactor_rwgt_MeanPDF_upperRMS"] = meanpdf_upperRMS_hist
    root_file["KFactor_rwgt_NNPDF31NLO_Central"] = nnpdf31nlo_central_hist 
    root_file["KFactor_rwgt_MSHT20NLO_Central"] = msht20nlo_central_hist 
    root_file["KFactor_rwgt_CT18NLO_Central"] = ct18nlo_central_hist 
    root_file["KFactor_rwgt_PDF4LHC21_Central"] = pdf4lhc21_central_hist 
    root_file["PDFRatio_rwgt_PDF4LHC21_LowerRMS"] = pdf4lhc21_lowerratio_hist
    root_file["PDFRatio_rwgt_PDF4LHC21_UpperRMS"] = pdf4lhc21_upperratio_hist
    root_file["PDFRatio_rwgt_NNPDF31NLO_LowerRMS"] = nnpdf31nlo_lowerratio_hist
    root_file["PDFRatio_rwgt_NNPDF31NLO_UpperRMS"] = nnpdf31nlo_upperratio_hist
    root_file["PDFRatio_rwgt_MSHT20NLO_LowerRMS"] = msht20nlo_lowerratio_hist
    root_file["PDFRatio_rwgt_MSHT20NLO_UpperRMS"] = msht20nlo_upperratio_hist
    root_file["PDFRatio_rwgt_CT18NLO_LowerRMS"] = ct18nlo_lowerratio_hist
    root_file["PDFRatio_rwgt_CT18NLO_UpperRMS"] = ct18nlo_upperratio_hist
    root_file["PDFRatio_rwgt_MeanPDF_lowerRMS"] = meanpdf_lowerratio_hist
    root_file["PDFRatio_rwgt_MeanPDF_upperRMS"] = meanpdf_upperratio_hist

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "KFactor_rwgt_MeanPDF_Central": [meanpdf_central_hist], 
    "KFactor_rwgt_MeanPDF_lowerRMS": [meanpdf_lowerRMS_hist],
    "KFactor_rwgt_MeanPDF_upperRMS": [meanpdf_upperRMS_hist],
    "KFactor_rwgt_NNPDF31NLO_Central": [nnpdf31nlo_central_hist],
    "KFactor_rwgt_MSHT20NLO_Central": [msht20nlo_central_hist],
    "KFactor_rwgt_CT18NLO_Central": [ct18nlo_central_hist],
    "KFactor_rwgt_PDF4LHC21_Central": [pdf4lhc21_central_hist],
    "PDFRatio_rwgt_PDF4LHC21_LowerRMS": [pdf4lhc21_lowerratio_hist],
    "PDFRatio_rwgt_PDF4LHC21_UpperRMS": [pdf4lhc21_upperratio_hist],
    "PDFRatio_rwgt_NNPDF31NLO_LowerRMS": [nnpdf31nlo_lowerratio_hist],
    "PDFRatio_rwgt_NNPDF31NLO_UpperRMS": [nnpdf31nlo_upperratio_hist],
    "PDFRatio_rwgt_MSHT20NLO_LowerRMS": [msht20nlo_lowerratio_hist],
    "PDFRatio_rwgt_MSHT20NLO_UpperRMS": [msht20nlo_upperratio_hist],
    "PDFRatio_rwgt_CT18NLO_LowerRMS": [ct18nlo_lowerratio_hist],
    "PDFRatio_rwgt_CT18NLO_UpperRMS": [ct18nlo_upperratio_hist],
    "PDFRatio_rwgt_MeanPDF_lowerRMS": [meanpdf_lowerratio_hist],
    "PDFRatio_rwgt_MeanPDF_upperRMS": [meanpdf_upperratio_hist]
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)


