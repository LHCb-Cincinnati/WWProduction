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
from scipy.interpolate import make_interp_spline, BSpline, UnivariateSpline
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

def create_heavyside_hist(hist, x0, y1, y2):
    heavyside_hist = hist.copy()
    for index, center in enumerate(heavyside_hist.axes[0].centers):
        if center < x0:
            heavyside_hist.view().value[index] = y1
            heavyside_hist.view().variance[index] = 1.0
        else:
            heavyside_hist.view().value[index] = y2
            heavyside_hist.view().variance[index] = 1.0
    return(heavyside_hist)

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
central_hist = hist_dict["DileptonKFactorFine_mu10"][0]
lowerenv_hist = hist_dict["DileptonKFactorFine_lowerenv"][0]
upperenv_hist = hist_dict["DileptonKFactorFine_upperenv"][0]
profile_hist = hist_dict["DiLeptonMassProfile"][0]

# Fit Functions
fit_func = flat
# Central Fit
central_least_squares = LeastSquares(
    profile_hist.view().value[1:-1], 
    central_hist.view().value[1:-1], 
    np.sqrt(central_hist.view().variance[1:-1]), 
    fit_func
)
central_m = Minuit(central_least_squares, a = 1.2)  # Starting values for DFDY flat fit.
central_m.migrad()  # finds minimum of least_squares function
central_m.hesse()  # accurately computes uncertainties
print("Central Value Fit Info:")
print(f"chi^2/n_dof = {central_m.fval:.1f} / {central_m.ndof:.0f} = {central_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(central_m.parameters, central_m.values, central_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")
central_heavyside_hist = create_heavyside_hist(
    central_hist, 
    33.0, 
    central_hist.view().value[0], 
    central_m.values[0]
)
# Lower Fit
lower_least_squares = LeastSquares(
    profile_hist.view().value[1:-1], 
    lowerenv_hist.view().value[1:-1], 
    np.sqrt(lowerenv_hist.view().variance[1:-1]), 
    fit_func
)
lower_m = Minuit(lower_least_squares, a = 1.2)  # Starting values for DFDY flat fit.
lower_m.migrad()  # finds minimum of least_squares function
lower_m.hesse()  # accurately computes uncertainties
print("Lower Value Fit Info:")
print(f"chi^2/n_dof = {lower_m.fval:.1f} / {lower_m.ndof:.0f} = {lower_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(lower_m.parameters, lower_m.values, lower_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")
lower_heavyside_hist = create_heavyside_hist(
    lowerenv_hist, 
    33.0, 
    lowerenv_hist.view().value[0], 
    lower_m.values[0]
)
# Upper Fit
upper_least_squares = LeastSquares(
    profile_hist.view().value[1:-1], 
    upperenv_hist.view().value[1:-1], 
    np.sqrt(upperenv_hist.view().variance[1:-1]), 
    fit_func
)
upper_m = Minuit(upper_least_squares, a = 1.2)  # Starting values for DFDY flat fit.
upper_m.migrad()  # finds minimum of least_squares function
upper_m.hesse()  # accurately computes uncertainties
print("Upper Value Fit Info:")
print(f"chi^2/n_dof = {upper_m.fval:.1f} / {upper_m.ndof:.0f} = {upper_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(upper_m.parameters, upper_m.values, upper_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")
upper_heavyside_hist = create_heavyside_hist(
    upperenv_hist, 
    33.0, 
    upperenv_hist.view().value[0], 
    upper_m.values[0]
)

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# K-Factor with scale variations
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    label="Central Value",
    color="black"
)
axs.errorbar(
   profile_hist.view().value,
   central_hist.view().value,
   ecolor = "black",
   linestyle = "",
   yerr = np.sqrt(central_hist.view().variance),
   label="Statistical Error"
)
axs.bar(
    x=central_hist.axes[0].centers, 
    height=2*(upper_heavyside_hist.view().value - lower_heavyside_hist.view().value), 
    bottom=lower_heavyside_hist.view().value, 
    width=central_hist.axes[0].widths, 
    # align='edge', 
    linewidth=0, 
    color="blue", 
    alpha=0.25, 
    # zorder=-1, 
    label="Scale Envelope"
)
# axs.stairs(
#     lowerenv_hist.view().value, 
#     edges=lowerenv_hist.axes[0].edges,
#     color="blue",
#     label="Lower Envelope"
# )
# axs.stairs(
#     upperenv_hist.view().value, 
#     edges=upperenv_hist.axes[0].edges,
#     color="red",
#     label="Upper Envelope"
# )
# # Setting axes and legend
ymax = max(upperenv_hist.view().value)
axs.set_ylim((0, 1.15*ymax))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
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
    fill = False,
    color = "black",
    label = "Reweight Factor"
    )
axs.errorbar(
    profile_hist.view().value,
    central_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = np.sqrt(central_hist.view().variance),
    label="Statistical Error"
)
axs.bar(
    central_heavyside_hist.axes[0].centers,
    central_heavyside_hist.view().value,
    width = central_heavyside_hist.axes[0].widths,
    fill = False,
    color = "orange",
    label = "Fit"
    )
ymax = np.max(
    central_hist.view().value 
    + np.sqrt(central_hist.view().variance)
)
# axs.set_ylim((0, 1.15*ymax))
axs.set_ylim((0, 1.6))
axs.set_xlim((0, 400))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Reweight Factor (NLO/LO)")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('KFactorwErrors.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["central_scale_hist"] = central_heavyside_hist
    root_file["lowerRMS_scale_hist"] = lower_heavyside_hist
    root_file["upperRMS_scale_hist"] = upper_heavyside_hist

# Save histograms
pickle_dict = {
    "central_scale_hist": [central_heavyside_hist],
    "lowerRMS_scale_hist": [lower_heavyside_hist],
    "upperRMS_scale_hist": [upper_heavyside_hist]
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)

