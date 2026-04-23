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

@dataclass
class FitConfig:
    """Configuration for a single histogram fit."""
    initial_params: dict
    limits: list

    DEFAULT_PARAMS = {"a": 0.5, "b": 0.88, "c": -0.02}
    DEFAULT_LIMITS = [(None, None), (None, None), (None, None)]

    @classmethod
    def default(cls) -> "FitConfig":
        return cls(
            initial_params={"a": 0.5, "b": 0.88, "c": -0.02},
            limits=[(None, None), (None, None), (None, None)],
        )

def fit_histograms(
    profile_hist,
    lower_rms_hist,
    central_hist,
    upper_rms_hist,
    fit_func,
    lower_config = FitConfig(
        initial_params={"a": 0.8, "b": 1, "c": -0.01},
        limits = [(0.6, 2), (0.4, 1.1), (None, None)]
    ),
    central_config = FitConfig.default(),
    upper_config = FitConfig(
        initial_params={"a": 0.8, "b": 1, "c": -0.01},
        limits = [(0.6, 2), (0.9, 1.1), (None, None)]
    )
):
    """
    Fit three boost histograms (lower RMS, central, upper RMS) using iminuit.

    Args:
        profile_hist:    Boost histogram providing the x-axis bin centers.
        lower_rms_hist:  Boost histogram for the lower RMS envelope.
        central_hist:    Boost histogram for the central values.
        upper_rms_hist:  Boost histogram for the upper RMS envelope.
        fit_func:        Callable with signature f(x, a, b, c) to fit against.
        initial_params:  Starting values for the Minuit minimiser.
        limits:          Parameter limits passed to Minuit.

    Returns:
        Tuple of (lower_m, central_m, upper_m) Minuit objects, all post-fit.
    """
    hists = [
        ("Lower RMS", lower_rms_hist, lower_config),
        ("Central",   central_hist,   central_config),
        ("Upper RMS", upper_rms_hist, upper_config),
    ]

    fitted = {}
    for label, hist, config in hists:
        ls = LeastSquares(
            profile_hist.view().value,
            hist.view().value,
            np.sqrt(hist.view().variance),
            fit_func,
        )
        m = Minuit(ls, **config.initial_params)
        m.limits = config.limits
        m.migrad()
        m.hesse()

        print(f"{label} Fit Info:")
        print(f"  chi^2/n_dof = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}")
        for p, v, e in zip(m.parameters, m.values, m.errors):
            print(f"  {p} = {v:.3f} +- {e:.3f}")

        fitted[label] = m

    return fitted

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
nnpdf31nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_LowerRMS"][0]
nnpdf31nlo_central_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_Central"][0]
nnpdf31nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_NNPDF31NLO_UpperRMS"][0]
msht20nlo_central_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_Central"][0]
msht20nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_LowerRMS"][0]
msht20nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_MSHT20NLO_UpperRMS"][0]
ct18nlo_central_hist = hist_dict["DileptonKFactorFine_CT18NLO_Central"][0]
ct18nlo_lowerRMS_hist = hist_dict["DileptonKFactorFine_CT18NLO_LowerRMS"][0]
ct18nlo_upperRMS_hist = hist_dict["DileptonKFactorFine_CT18NLO_UpperRMS"][0]

# Fit K-Factor
fit_func = exp_decay

# Mean PDF Minuit
meanpdf_minuit = fit_histograms(
    profile_hist, meanpdf_lowerRMS_hist, meanpdf_central_hist, meanpdf_upperRMS_hist, fit_func
)
nnpdf31nlo_minuit = fit_histograms(
    profile_hist, nnpdf31nlo_lowerRMS_hist, nnpdf31nlo_central_hist, nnpdf31nlo_upperRMS_hist, fit_func
)
msht20nlo_minuit = fit_histograms(
    profile_hist, msht20nlo_lowerRMS_hist, msht20nlo_central_hist, msht20nlo_upperRMS_hist, fit_func
)
ct18nlo_minuit = fit_histograms(
    profile_hist, ct18nlo_lowerRMS_hist, ct18nlo_central_hist, ct18nlo_upperRMS_hist, fit_func
)

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
# Setting axes and legend
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
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
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
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
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
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
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
fig.savefig('KFactorwPDFFamilies_CT18NLO.png')
plt.close()

# Plot reweight histogram with errors
# Mean PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
x_fine = np.linspace(meanpdf_central_hist.axes[0].edges[0], meanpdf_central_hist.axes[0].edges[-1], 300)
axs.stairs(
    meanpdf_central_hist.view().value, 
    edges=meanpdf_central_hist.axes[0].edges,
    color="black"
)
axs.errorbar(
    profile_hist.view().value,
    meanpdf_central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = (meanpdf_upperRMS_hist.view().value - meanpdf_lowerRMS_hist.view().value) / 2,
    label="RMS Error"
)
axs.plot(
    x_fine, 
    fit_func(x_fine, *meanpdf_minuit["Central"].values),
    color="black", label="Central Value Fit"
)
axs.fill_between(
    x_fine, 
    fit_func(x_fine, *meanpdf_minuit["Lower RMS"].values), 
    fit_func(x_fine, *meanpdf_minuit["Upper RMS"].values), 
    color='blue', alpha=0.3, label='Envelope'
)
# axs.errorbar(
#     profile_hist.view().value,
#     meanpdf_lowerRMS_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = meanpdf_lowerRMS_hist.view().variance,
#     label="Lower Envelope Statistical Error"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     meanpdf_upperRMS_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = meanpdf_upperRMS_hist.view().variance,
#     label="Upper Envelope Statistical Error"
# )
# Setting axes and legend
ymax = max(meanpdf_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
text_x = bbox.x0 + 0.1
text_y = bbox.y1 + 0.05  # small offset below the legend
axs.text(
    text_x, text_y,
    f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {meanpdf_minuit['Central'].fval:.1f} / {meanpdf_minuit['Central'].ndof:.0f} = {meanpdf_minuit['Central'].fmin.reduced_chi2:.1f}",
    transform=axs.transAxes,
    fontsize=axs.legend().get_texts()[0].get_fontsize(),
    verticalalignment='top',
    # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_MeanPDF.png')
plt.close()

# NNPDF31NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
x_fine = np.linspace(nnpdf31nlo_central_hist.axes[0].edges[0], nnpdf31nlo_central_hist.axes[0].edges[-1], 300)
axs.stairs(
    nnpdf31nlo_central_hist.view().value, 
    edges=nnpdf31nlo_central_hist.axes[0].edges,
    color="black"
)
axs.errorbar(
    profile_hist.view().value,
    nnpdf31nlo_central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = (nnpdf31nlo_upperRMS_hist.view().value - nnpdf31nlo_lowerRMS_hist.view().value) / 2,
    label="RMS Error"
)
axs.plot(
    x_fine, 
    fit_func(x_fine, *nnpdf31nlo_minuit["Central"].values),
    color="black", label="Central Value Fit"
)
axs.fill_between(
    x_fine, 
    fit_func(x_fine, *nnpdf31nlo_minuit["Lower RMS"].values), 
    fit_func(x_fine, *nnpdf31nlo_minuit["Upper RMS"].values), 
    color='blue', alpha=0.3, label='Envelope'
)
# axs.errorbar(
#     profile_hist.view().value,
#     nnpdf31nlo_lowerRMS_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = nnpdf31nlo_lowerRMS_hist.view().variance,
#     label="Lower Envelope Statistical Error"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     nnpdf31nlo_upperRMS_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = nnpdf31nlo_upperRMS_hist.view().variance,
#     label="Upper Envelope Statistical Error"
# )
# Setting axes and legend
ymax = max(nnpdf31nlo_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
text_x = bbox.x0 + 0.1
text_y = bbox.y1 + 0.05  # small offset below the legend
axs.text(
    text_x, text_y,
    f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {nnpdf31nlo_minuit['Central'].fval:.1f} / {nnpdf31nlo_minuit['Central'].ndof:.0f} = {nnpdf31nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
    transform=axs.transAxes,
    fontsize=axs.legend().get_texts()[0].get_fontsize(),
    verticalalignment='top',
    # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_NNPDF31NLO.png')
plt.close()

# MSHT20NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
x_fine = np.linspace(msht20nlo_central_hist.axes[0].edges[0], msht20nlo_central_hist.axes[0].edges[-1], 300)
axs.stairs(
    msht20nlo_central_hist.view().value, 
    edges=msht20nlo_central_hist.axes[0].edges,
    color="black"
)
axs.errorbar(
    profile_hist.view().value,
    msht20nlo_central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = (msht20nlo_upperRMS_hist.view().value - msht20nlo_lowerRMS_hist.view().value) / 2,
    label="RMS Error"
)
axs.plot(
    x_fine, 
    fit_func(x_fine, *msht20nlo_minuit["Central"].values),
    color="black", label="Central Value Fit"
)
axs.fill_between(
    x_fine, 
    fit_func(x_fine, *msht20nlo_minuit["Lower RMS"].values), 
    fit_func(x_fine, *msht20nlo_minuit["Upper RMS"].values), 
    color='blue', alpha=0.3, label='Envelope'
)
# axs.errorbar(
#     profile_hist.view().value,
#     msht20nlo_lowerRMS_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = msht20nlo_lowerRMS_hist.view().variance,
#     label="Lower Envelope Statistical Error"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     msht20nlo_upperRMS_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = msht20nlo_upperRMS_hist.view().variance,
#     label="Upper Envelope Statistical Error"
# )
# Setting axes and legend
ymax = max(msht20nlo_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
text_x = bbox.x0 + 0.1
text_y = bbox.y1 + 0.05  # small offset below the legend
axs.text(
    text_x, text_y,
    f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {msht20nlo_minuit['Central'].fval:.1f} / {msht20nlo_minuit['Central'].ndof:.0f} = {msht20nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
    transform=axs.transAxes,
    fontsize=axs.legend().get_texts()[0].get_fontsize(),
    verticalalignment='top',
    # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_MSHT20NLO.png')
plt.close()

# CT18NLO PDF
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
x_fine = np.linspace(ct18nlo_central_hist.axes[0].edges[0], ct18nlo_central_hist.axes[0].edges[-1], 300)
axs.stairs(
    ct18nlo_central_hist.view().value, 
    edges=ct18nlo_central_hist.axes[0].edges,
    color="black"
)
axs.errorbar(
    profile_hist.view().value,
    ct18nlo_central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = (ct18nlo_upperRMS_hist.view().value - ct18nlo_lowerRMS_hist.view().value) / 2,
    label="RMS Error"
)
axs.plot(
    x_fine, 
    fit_func(x_fine, *ct18nlo_minuit["Central"].values),
    color="black", label="Central Value Fit"
)
axs.fill_between(
    x_fine, 
    fit_func(x_fine, *ct18nlo_minuit["Lower RMS"].values), 
    fit_func(x_fine, *ct18nlo_minuit["Upper RMS"].values), 
    color='blue', alpha=0.3, label='Envelope'
)
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_lowerRMS_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = ct18nlo_lowerRMS_hist.view().variance,
#     label="Lower Envelope Statistical Error"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     ct18nlo_upperRMS_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = ct18nlo_upperRMS_hist.view().variance,
#     label="Upper Envelope Statistical Error"
# )
# Setting axes and legend
ymax = max(ct18nlo_upperRMS_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Setting fit box above legend
bbox = axs.legend().get_window_extent().transformed(axs.transAxes.inverted())
text_x = bbox.x0 + 0.1
text_y = bbox.y1 + 0.05  # small offset below the legend
axs.text(
    text_x, text_y,
    f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {ct18nlo_minuit['Central'].fval:.1f} / {ct18nlo_minuit['Central'].ndof:.0f} = {ct18nlo_minuit['Central'].fmin.reduced_chi2:.1f}",
    transform=axs.transAxes,
    fontsize=axs.legend().get_texts()[0].get_fontsize(),
    verticalalignment='top',
    # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwPDFVariations_CT18NLO.png')
plt.close()

# Create ROOT histogram
# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    for family_name, minuit_dict in zip(["MeanPDF", "NNPDF31NLO", "MSHT20NLO", "CT18NLO"], [meanpdf_minuit, nnpdf31nlo_minuit, msht20nlo_minuit, ct18nlo_minuit]):
        for fit_name in ["Central", "Upper RMS", "Lower RMS"]:
            parameterization = minuit_dict[fit_name]
            n_params = len(parameterization.params)
            param_bins = array.array(
                'd', [float(x) for x in np.arange(-0.5, n_params, 1)]
            )
            param_hist = ROOT.TH2D(
                f"K-Factor Parameterization Hist {family_name} {fit_name}",
                f"K-Factor Parameterization Hist {family_name} {fit_name}",
                1,
                array.array('d', [0,1]),
                n_params,
                param_bins
            )
            for index in range(n_params):
                bin_index = param_hist.GetBin(1, index + 1)
                param_hist.SetBinContent(bin_index, parameterization.values[index])
                param_hist.SetBinError(bin_index, parameterization.errors[index])
            root_file[f"rwgt_{family_name}{fit_name}_param_hist"] = param_hist

