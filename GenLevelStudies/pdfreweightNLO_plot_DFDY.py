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
from collections import namedtuple
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
lowerRMS_hist = hist_dict["DileptonKFactorFine_LowerRMS"][0]
central_hist = hist_dict["DileptonKFactorFine_Central"][0]
upperRMS_hist = hist_dict["DileptonKFactorFine_UpperRMS"][0]
nnpdf31nlo_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO"][0]
msht20nlo_hist = hist_dict["DileptonKFactorFine_MSHT20NLO"][0]
ct18nlo_hist = hist_dict["DileptonKFactorFine_CT18NLO"][0]

# # Fit K-Factor
# fit_func = exp_decay
# # Central Value Fit
# least_squares = LeastSquares(
#     profile_hist.view().value, 
#     central_hist.view().value, 
#     central_hist.view().variance, 
#     fit_func
# )
# m = Minuit(least_squares, a = 0.5, b = 0.88, c = -0.02)  # Starting values for ttbar exponential fit.
# # m.limits = [(0, 20), (0.75, 1), (None, None)] # Bounds for ttbar exponential fit.
# m.migrad()  # finds minimum of least_squares function
# m.hesse()  # accurately computes uncertainties
# print("Central Value Fit Info:")
# print(f"chi^2/n_dof = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}")
# for p, v, e in zip(m.parameters, m.values, m.errors):
#     print(f"{p} = {v:.3f} +- {e:.3f}")
# 
# # Lower Envelope Fit
# lowerRMS_least_squares = LeastSquares(
#     profile_hist.view().value, 
#     lowerRMS_hist.view().value, 
#     lowerRMS_hist.view().variance, 
#     fit_func
# )
# lowerenv_m = Minuit(lowerRMS_least_squares, a = 0.8, b = 1.0, c = -0.01)  # Starting values for ttbar exponential fit.
# # lowerenv_m.limits = [(0.6, 2), (0.4, 1.1), (None, None)] # Bounds for ttbar exponential fit.
# lowerenv_m.limits = [(None, None), (None, None), (None, None)] # Bounds for ttbar exponential fit.
# lowerenv_m.migrad()  # finds minimum of least_squares function
# lowerenv_m.hesse()  # accurately computes uncertainties
# print("Lower Envelope Fit Info:")
# print(f"chi^2/n_dof = {lowerenv_m.fval:.1f} / {lowerenv_m.ndof:.0f} = {lowerenv_m.fmin.reduced_chi2:.1f}")
# for p, v, e in zip(lowerenv_m.parameters, lowerenv_m.values, lowerenv_m.errors):
#     print(f"{p} = {v:.3f} +- {e:.3f}")
# 
# # Upper Envelope Fit
# upperenv_least_squares = LeastSquares(
#     profile_hist.view().value, 
#     upperRMS_hist.view().value, 
#     upperRMS_hist.view().variance, 
#     fit_func
# )
# # upperenv_m = Minuit(upperenv_least_squares, a = 1.5)  # Starting values for DFDY flat fit.
# upperenv_m = Minuit(upperenv_least_squares, a = 0.8, b = 1.0, c = -0.01)  # Starting values for ttbar exponential fit.
# upperenv_m.limits = [(0.6, 2), (0.9, 1.1), (None, None)] # Bounds for ttbar exponential fit.
# upperenv_m.migrad()  # finds minimum of least_squares function
# upperenv_m.hesse()  # accurately computes uncertainties
# print("Upper Envelope Fit Info:")
# print(f"chi^2/n_dof = {upperenv_m.fval:.1f} / {upperenv_m.ndof:.0f} = {upperenv_m.fmin.reduced_chi2:.1f}")
# for p, v, e in zip(upperenv_m.parameters, upperenv_m.values, upperenv_m.errors):
#     print(f"{p} = {v:.3f} +- {e:.3f}")

# # Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# K-Factor with scale variations of PDF families
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=central_hist.axes[0].centers, 
    height=2*(upperRMS_hist.view().value - central_hist.view().value), 
    bottom=lowerRMS_hist.view().value, 
    width=central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_hist.view().value, 
    edges=nnpdf31nlo_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    msht20nlo_hist.view().value, 
    edges=msht20nlo_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    ct18nlo_hist.view().value, 
    edges=ct18nlo_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
# Setting axes and legend
ymax = max(upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
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
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=central_hist.axes[0].centers, 
    height=2*(upperRMS_hist.view().value - central_hist.view().value), 
    bottom=lowerRMS_hist.view().value, 
    width=central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    nnpdf31nlo_hist.view().value, 
    edges=nnpdf31nlo_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.bar(
    x=nnpdf31nlo_hist.axes[0].centers, 
    height=2*(nnpdf31nlo_hist.view().variance), 
    bottom=(nnpdf31nlo_hist.view().value - nnpdf31nlo_hist.view().variance), 
    width=nnpdf31nlo_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="green", 
    alpha=0.25, 
    # zorder=-1, 
    label="NNPDF31NLO Envelope"
)
axs.stairs(
    msht20nlo_hist.view().value, 
    edges=msht20nlo_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.stairs(
    ct18nlo_hist.view().value, 
    edges=ct18nlo_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
ymax = max(upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_NNPDF31NLO.png')
plt.close()
# MSHT20NLO
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=central_hist.axes[0].centers, 
    height=2*(upperRMS_hist.view().value - central_hist.view().value), 
    bottom=lowerRMS_hist.view().value, 
    width=central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    msht20nlo_hist.view().value, 
    edges=msht20nlo_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
axs.bar(
    x=msht20nlo_hist.axes[0].centers, 
    height=2*(msht20nlo_hist.view().variance), 
    bottom=(msht20nlo_hist.view().value - msht20nlo_hist.view().variance), 
    width=msht20nlo_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="red", 
    alpha=0.25, 
    # zorder=-1, 
    label="MSHT20NLO Envelope"
)
axs.stairs(
    nnpdf31nlo_hist.view().value, 
    edges=nnpdf31nlo_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    ct18nlo_hist.view().value, 
    edges=ct18nlo_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
ymax = max(upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_MSHT20NLO.png')
plt.close()
# CT18NLO
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black",
    label="Mean PDF Central Value",
    zorder=3
)
axs.bar(
    x=central_hist.axes[0].centers, 
    height=2*(upperRMS_hist.view().value - central_hist.view().value), 
    bottom=lowerRMS_hist.view().value, 
    width=central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    zorder=2, 
    label="Mean PDF RMS Envelope"
)
axs.stairs(
    ct18nlo_hist.view().value, 
    edges=ct18nlo_hist.axes[0].edges,
    color="blue",
    label="CT18NLO"
)
axs.bar(
    x=ct18nlo_hist.axes[0].centers, 
    height=2*(ct18nlo_hist.view().variance), 
    bottom=(ct18nlo_hist.view().value - ct18nlo_hist.view().variance), 
    width=ct18nlo_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="blue", 
    alpha=0.25, 
    # zorder=-1, 
    label="CT18NLO Envelope"
)
axs.stairs(
    nnpdf31nlo_hist.view().value, 
    edges=nnpdf31nlo_hist.axes[0].edges,
    color="green",
    label="NNPDF31NLO"
)
axs.stairs(
    msht20nlo_hist.view().value, 
    edges=msht20nlo_hist.axes[0].edges,
    color="red",
    label="MSHT20NLO"
)
ymax = max(upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_CT18NLO.png')
plt.close()

# Plot reweight histogram with errors
# fig, axs = plt.subplots()
# plt.subplots_adjust(top=0.85)
# x_fine = np.linspace(central_hist.axes[0].edges[0], central_hist.axes[0].edges[-1], 300)
# axs.stairs(
#     central_hist.view().value, 
#     edges=central_hist.axes[0].edges,
#     color="black"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     central_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     yerr = central_hist.view().variance,
#     label="RMS Envelope"
# )
# axs.plot(
#     x_fine, 
#     fit_func(x_fine, *m.values),
#     color="black", label="Central Value Fit"
# )
# axs.fill_between(
#     x_fine, 
#     fit_func(x_fine, *lowerenv_m.values), 
#     fit_func(x_fine, *upperenv_m.values), 
#     color='blue', alpha=0.3, label='Envelope'
# )
# # axs.errorbar(
# #     profile_hist.view().value,
# #     lowerRMS_hist.view().value,
# #     ecolor = "blue",
# #     linestyle = "",
# #     yerr = lowerRMS_hist.view().variance,
# #     label="Lower Envelope Statistical Error"
# # )
# # axs.errorbar(
# #     profile_hist.view().value,
# #     upperRMS_hist.view().value,
# #     ecolor = "red",
# #     linestyle = "",
# #     yerr = upperRMS_hist.view().variance,
# #     label="Upper Envelope Statistical Error"
# # )
# # Setting axes and legend
# ymax = max(upperRMS_hist.view().value)
# axs.set_ylim((0, 1.15*ymax))
# axs.set_xlim((0, 300))
# axs.set_title("")
# axs.set_xlabel("$M_{e \\mu} (GeV)$")
# axs.set_ylabel("Reweight Factor (NLO/LO)")
# hist_handles, hist_labels = axs.get_legend_handles_labels()
# axs.legend(hist_handles, hist_labels)
# # Setting fit box above legend
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
# # Slightly fancy to remove whitespace
# fig.savefig('KFactorwPDFVariations.png')
# plt.close()

# Create ROOT histogram
# Save Histos in ROOT
# os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
# with uproot.recreate(ofile_name + ".root") as root_file:
#     for name, parameterization in zip(["Central", "UpperRMS", "LowerRMS"],[m, upperenv_m, lowerenv_m]):
#         n_params = len(parameterization.params)
#         param_bins = array.array(
#             'd', [float(x) for x in np.arange(-0.5, n_params, 1)]
#         )
#         param_hist = ROOT.TH2D(
#             f"K-Factor Parameterization Hist {name} PDF",
#             f"K-Factor Parameterization Hist {name} PDF",
#             1,
#             array.array('d', [0,1]),
#             n_params,
#             param_bins
#         )
#         for index in range(n_params):
#             bin_index = param_hist.GetBin(1, index + 1)
#             param_hist.SetBinContent(bin_index, parameterization.values[index])
#             param_hist.SetBinError(bin_index, parameterization.errors[index])
#         root_file[f"rwgt_{name}PDF_hist"] = param_hist

