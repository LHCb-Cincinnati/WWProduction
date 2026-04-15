""" Program to compute the K-factor reweighting for ttbar and DFDY process.

Program to compute the K-factor reweighting for ttbar and DFDY process.  
The input for this program should be the lower envelope histogram pickle file,
central value histogram pickle file, and the upper envelope histogram pickle 
file in that order.  You may need to change the starting values for each of the
central value, lower/upper env histograms to get them to converge correctly.
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

# Functions
# Flat function
def flat(x, a):
    return(0*x + a)

# Linear function
def linear(x, a, b):
    return(a*x + b)

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
    data = pickle.load(f)

# Grab Histograms from File.
central_hist = data["DileptonKFactorFine_mu10"][0]
lower_hist = data["DileptonKFactorFine_lowerenv"][0]
upper_hist = data["DileptonKFactorFine_upperenv"][0]
profile_hist = data["DiLeptonMassProfile"][0]

# Fit K-Factor
fit_func = exp_decay
# Central Value Fit
least_squares = LeastSquares(
    profile_hist.view().value, 
    central_hist.view().value, 
    np.sqrt(central_hist.view().variance), 
    fit_func
)
# m = Minuit(least_squares, a = 1)  # Starting values for DFDY flat fit.
m = Minuit(least_squares, a = 0.5, b = 0.88, c = -0.02)  # Starting values for ttbar exponential fit.
m.limits = [(0, 20), (0.6, 0.9), (None, None)] # Bounds for ttbar exponential fit.
m.migrad()  # finds minimum of least_squares function
m.hesse()  # accurately computes uncertainties
print("Central Value Fit Info:")
print(f"chi^2/n_dof = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(m.parameters, m.values, m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")

# Lower Envelope Fit
lowerenv_least_squares = LeastSquares(
    profile_hist.view().value, 
    lower_hist.view().value, 
    np.sqrt(lower_hist.view().variance), 
    fit_func
)
# lowerenv_m = Minuit(lowerenv_least_squares, a = 0.8)  # Starting values for DFDY flat fit.
lowerenv_m = Minuit(lowerenv_least_squares, a = 0.5, b = 0.75, c = -0.1)  # Starting values for ttbar exponential fit.
lowerenv_m.limits = [(0.1, 5), (0.5, 0.8), (-0.2, 0)] # Bounds for ttbar exponential fit.
# lowerenv_m = Minuit(least_squares, a = 0.5, b = 0.88, c = -0.02)  # Starting values for ttbar exponential fit.
# lowerenv_m.limits = [(0, 20), (0.6, 1), (None, None)] # Bounds for ttbar exponential fit.
lowerenv_m.migrad()  # finds minimum of least_squares function
lowerenv_m.hesse()  # accurately computes uncertainties
print("Lower Envelope Fit Info:")
print(f"chi^2/n_dof = {lowerenv_m.fval:.1f} / {lowerenv_m.ndof:.0f} = {lowerenv_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(lowerenv_m.parameters, lowerenv_m.values, lowerenv_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")

# Upper Envelope Fit
upperenv_least_squares = LeastSquares(
    profile_hist.view().value, 
    upper_hist.view().value, 
    np.sqrt(upper_hist.view().variance), 
    fit_func
)
# upperenv_m = Minuit(upperenv_least_squares, a = 1.5)  # Starting values for DFDY flat fit.
upperenv_m = Minuit(upperenv_least_squares, a = 0.8, b = 1.0, c = -0.01)  # Starting values for ttbar exponential fit.
upperenv_m.limits = [(0.6, 2), (0.9, 1.1), (None, None)] # Bounds for ttbar exponential fit.
# upperenv_m = Minuit(least_squares, a = 0.5, b = 0.88, c = -0.02)  # Starting values for ttbar exponential fit.
# upperenv_m.limits = [(0, 20), (0.8, 1), (None, None)] # Bounds for ttbar exponential fit.
upperenv_m.migrad()  # finds minimum of least_squares function
upperenv_m.hesse()  # accurately computes uncertainties
print("Upper Envelope Fit Info:")
print(f"chi^2/n_dof = {upperenv_m.fval:.1f} / {upperenv_m.ndof:.0f} = {upperenv_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(upperenv_m.parameters, upperenv_m.values, upperenv_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")

# # Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# K-Factor with scale variations
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
x_fine = np.linspace(central_hist.axes[0].edges[0], central_hist.axes[0].edges[-1], 300)
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black"
)
axs.errorbar(
    profile_hist.view().value, central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = np.sqrt(central_hist.view().variance),
    label="Statistical Error"
)
axs.plot(
    x_fine, 
    fit_func(x_fine, *m.values),
    color="black", label="Central Value Fit"
)
axs.fill_between(
    x_fine, 
    fit_func(x_fine, *lowerenv_m.values), 
    fit_func(x_fine, *upperenv_m.values), 
    color='blue', alpha=0.3, label='Envelope'
)
# axs.errorbar(
#     profile_hist.view().value,
#     lower_hist.view().value,
#     ecolor = "blue",
#     linestyle = "",
#     yerr = np.sqrt(lower_hist.view().variance),
#     label="Lower Envelope Statistical Error"
# )
# axs.errorbar(
#     profile_hist.view().value,
#     upper_hist.view().value,
#     ecolor = "red",
#     linestyle = "",
#     yerr = np.sqrt(upper_hist.view().variance),
#     label="Upper Envelope Statistical Error"
# )
# Setting axes and legend
ymax = max(upper_hist.view().value)
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
    f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}",
    transform=axs.transAxes,
    fontsize=axs.legend().get_texts()[0].get_fontsize(),
    verticalalignment='top',
    # bbox=dict(boxstyle="round,pad=0.3", fc='white', ec='black')
)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwScaleVariations.png')
plt.close()

# Plot reweight histogram with errors
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.bar(
    central_hist.axes[0].centers,
    central_hist.view().value,
    width = central_hist.axes[0].widths,
    # yerr = np.sqrt(kfactor_hist_list[1].view().variance),
    fill = False,
    label = "Reweight Factor"
    )
axs.errorbar(
    profile_hist.view().value,
    central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    # width = kfactor_hist_list[1].axes[0].widths,
    yerr = np.sqrt(central_hist.view().variance),
    label="Statistical Error"
)
# axs.plot(
#     kfactor_hist_list[1].axes[0].edges, 
#     fit_func(kfactor_hist_list[1].axes[0].edges, *m.values),
#     label = "Fit"
# )
ymax = np.max(
    central_hist.view().value 
    + np.sqrt(central_hist.view().variance)
)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 300))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
# hist_handles, hist_labels = axs.get_legend_handles_labels()
# axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwErrors.png')
plt.close()

# Create ROOT histogram
# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    n_params = len(upperenv_m.params)
    param_bins = array.array(
        'd', [float(x) for x in np.arange(-0.5, n_params, 1)]
    )
    # Central Parameter Hist
    central_param_hist = ROOT.TH2D(
        "Central K-Factor Parameterization Hist",
        "Central K-Factor Parameterization Hist",
        1,
        array.array('d', [0,1]),
        n_params,
        param_bins
    )
    for index in range(n_params):
        bin_index = central_param_hist.GetBin(1, index + 1)
        central_param_hist.SetBinContent(bin_index, m.values[index])
        central_param_hist.SetBinError(bin_index, m.errors[index])
    root_file["central_rwgt_hist"] = central_param_hist
    # Lower Parameter Hist
    lower_param_hist = ROOT.TH2D(
        "Lower K-Factor Parameterization Hist",
        "Lower K-Factor Parameterization Hist",
        1,
        array.array('d', [0,1]),
        n_params,
        param_bins
    )
    for index in range(n_params):
        bin_index = lower_param_hist.GetBin(1, index + 1)
        lower_param_hist.SetBinContent(bin_index, lowerenv_m.values[index])
        lower_param_hist.SetBinError(bin_index, lowerenv_m.errors[index])
    root_file["lower_rwgt_hist"] = lower_param_hist
    # Upper Parameter Hist
    upper_param_hist = ROOT.TH2D(
        "Upper K-Factor Parameterization Hist",
        "Upper K-Factor Parameterization Hist",
        1,
        array.array('d', [0,1]),
        n_params,
        param_bins
    )
    for index in range(n_params):
        bin_index = upper_param_hist.GetBin(1, index + 1)
        upper_param_hist.SetBinContent(bin_index, upperenv_m.values[index])
        upper_param_hist.SetBinError(bin_index, upperenv_m.errors[index])
    root_file["upper_rwgt_hist"] = upper_param_hist

