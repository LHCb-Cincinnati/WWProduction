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
def hist_to_function(hist):
    """
    Convert a 1D boost_histogram.Histogram into a callable function
    that returns the histogram value at a given x.
    """
    # Extract bin edges and histogram values
    edges = hist.axes[0].edges
    values = hist.view()

    # Handle potential shape mismatches (e.g. flow bins)
    if len(values) != len(edges) - 1:
        # Attempt to trim or pad appropriately
        values = values[:len(edges) - 1]

    def hist_func(x):
        """
        Given an x-value (or array), return the corresponding histogram value.
        Values outside the histogram range return 0.
        """
        x = np.atleast_1d(x)
        # np.digitize returns the index of the bin (1-based)
        idx = np.digitize(x, edges) - 1
        # Clamp indices to valid bins
        valid = (idx >= 0) & (idx < len(values))
        result = np.zeros_like(x, dtype=float)
        result[valid] = values[idx[valid]].value
        # Return scalar if input was scalar
        return result[0] if np.isscalar(x) else result

    return hist_func
def histogram_points(hist):
    """
    Given a 1D boost-histogram histogram, return two arrays:
      x_vals: sorted array of bin edges and centers
      y_vals: histogram values at those x positions

    The bin value is duplicated at both edges and the center
    of its bin to allow for plotting or interpolation.
    """
    edges = hist.axes[0].edges
    centers = hist.axes[0].centers
    values = hist.view()

    # Build combined arrays for plotting or analysis
    x_vals = []
    y_vals = []

    for i in range(len(values)):
        # For each bin, record left edge, center, and right edge
        center = centers[i]
        right = edges[i + 1]
        val = values[i].value
        if i == 0:
            left = edges[i]
            x_vals.extend([left, center, right])
            y_vals.extend([val, val, val])
        else:
            x_vals.extend([center, right])
            y_vals.extend([val, val])

    # Convert to numpy arrays and sort by x for convenience
    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    order = np.argsort(x_vals)

    return x_vals[order], y_vals[order]

# Flat function
def flat(x, a):
    return(0*x + a)

# Linear function
def linear(x, a, b):
    return(a*x + b)

# Exponental Decay function
def exp_decay(x, a, b, c):
    return(a*np.exp(c*x) + b)

# Start Logger
logging.basicConfig(
    filename='StackedHist.log',
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p'
)
plt.set_loglevel (level = 'warning')
# get the the logger with the name 'PIL'
pil_logger = logging.getLogger('PIL')  
# override the logger logging level to INFO
pil_logger.setLevel(logging.INFO)

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



# # Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# K-Factor with scale variations
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
# lower_func = hist_to_function(kfactor_hist_list[0])
# upper_func = hist_to_function(kfactor_hist_list[2])
# x_sparse = np.linspace(kfactor_hist_list[0].axes[0].edges[0], kfactor_hist_list[0].axes[0].edges[-1], 300)
# x_fine = np.linspace(kfactor_hist_list[0].axes[0].edges[0], kfactor_hist_list[0].axes[0].edges[-1], 200)
# lower_spl = make_interp_spline(x_sparse[:-1], lower_func(x_sparse)[:-1], k=3) 
# lower_env = lower_spl(x_fine)    
# upper_spl = make_interp_spline(x_sparse[:-1], upper_func(x_sparse)[:-1], k=3) 
# upper_env = upper_spl(x_fine)    
axs.stairs(
    central_hist.view().value, 
    edges=central_hist.axes[0].edges,
    color="black"
)
# axs.errorbar(
#     profile_hist.view().value,
#     central_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     yerr = np.sqrt(central_hist.view().variance),
#     label="Statistical Error"
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
axs.stairs(
    lowerenv_hist.view().value, 
    edges=lowerenv_hist.axes[0].edges,
    color="blue",
    label="Lower Envelope"
)
# axs.errorbar(
#    profile_hist.view().value,
#    lowerenv_hist.view().value,
#    ecolor = "blue",
#    linestyle = "",
#    yerr = np.sqrt(lowerenv_hist.view().variance),
#    label="Lower Envelope Statistical Error"
# )
axs.stairs(
    upperenv_hist.view().value, 
    edges=upperenv_hist.axes[0].edges,
    color="red",
    label="Upper Envelope"
)
# axs.errorbar(
#    profile_hist.view().value,
#    upperenv_hist.view().value,
#    ecolor = "red",
#    linestyle = "",
#    yerr = np.sqrt(upperenv_hist.view().variance),
#    label="Upper Envelope Statistical Error"
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
fig.savefig('KFactorwScaleVariations.png')
plt.close()

# Plot reweight histogram with errors
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
spl = UnivariateSpline(
    profile_hist.view().value, 
    central_hist.view().value,
    k=3,
    s=0
)
x_new = np.linspace(0, 400, 1000)
y_smooth = spl(x_new)
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
plt.plot(x_new, y_smooth, '-', label='Smoothed Spline Fit')
# axs.plot(
#     kfactor_hist_list[1].axes[0].edges, 
#     fit_func(kfactor_hist_list[1].axes[0].edges, *m.values),
#     label = "Fit"
# )
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

# Create ROOT histogram
# Save Histos in ROOT
# os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
# with uproot.recreate(ofile_name + ".root") as root_file:
#     n_params = len(upperenv_m.params)
#     param_bins = array.array(
#         'd', [float(x) for x in np.arange(-0.5, n_params, 1)]
#     )
#     param_hist = ROOT.TH2D(
#         "K-Factor Parameterization Hist",
#         "K-Factor Parameterization Hist",
#         1,
#         array.array('d', [0,1]),
#         n_params,
#         param_bins
#     )
#     for index in range(n_params):
#         bin_index = param_hist.GetBin(1, index + 1)
#         param_hist.SetBinContent(bin_index, upperenv_m.values[index])
#         param_hist.SetBinError(bin_index, upperenv_m.errors[index])
#     root_file["rwgt_hist"] = param_hist
