# Imports
# Standard Library Imports
import os
import pdb
import pickle
import sys
import logging
from collections import namedtuple
# Scikit Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
# HEP Imports
from hepunits.units import MeV, GeV
import awkward as ak
import uproot
import boost_histogram as bh
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
data_list = [0] * len(file_list)
label_list = ["Lower Envelope", "Central Value", "Upper Envelope"]
logging.info(f"Label list provided: {label_list}")

# Retrieve Histos
# Here data are the histograms and index is the index of the file in the
# list of files. So data_list looks like:
# data_list[1] = dictionary first file's histograms
# The dictionary looks like
# dict[Histogram Tag] = Histogram
for index, file in enumerate(file_list):
    with open(file, 'rb') as f:
        data = pickle.load(f)
        logging.debug(f"Keys taken from {index} file: {data.keys()}")
        data_list[index] = data

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Loop through the histogram tags available from the first file
for hist_tag in data_list[0].keys():
    # Add histogram to hist_draw_list from first file.
    logging.debug(f"Drawing {hist_tag}")
    hist_draw_list = data_list[0][hist_tag]
    logging.debug(f"Added {data_list[0][hist_tag]} to hist_draw_list")
    # Loop through other files and add other histograms with the
    # same tag to the hist_draw_list
    for file in data_list[1:]:
        # Check if tag is in histogram
        if hist_tag not in file.keys():
            raise RuntimeWarning(f"Tag: {hist_tag} not found in file: {file}.")
        hist_draw_list+=file[hist_tag]
        logging.debug(f"Added {file[hist_tag]} to hist_draw_list")
    # K-Factor with scale variations
    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    lower_func = hist_to_function(hist_draw_list[0])
    upper_func = hist_to_function(hist_draw_list[2])
    x_sparse = np.linspace(hist_draw_list[0].axes[0].edges[0], hist_draw_list[0].axes[0].edges[-1], 300)
    x_fine = np.linspace(hist_draw_list[0].axes[0].edges[0], hist_draw_list[0].axes[0].edges[-1], 200)
    # x_arr, lower_env = histogram_points(hist_draw_list[0])
    # x_arr, upper_env = histogram_points(hist_draw_list[2])
    lower_spl = make_interp_spline(x_sparse[:-1], lower_func(x_sparse)[:-1], k=3) 
    lower_env = lower_spl(x_fine)    
    upper_spl = make_interp_spline(x_sparse[:-1], upper_func(x_sparse)[:-1], k=3) 
    upper_env = upper_spl(x_fine)    
    #axs.plot(hist_draw_list[0].axes[0].centers, hist_draw_list[0].view().value,
    #            label=label_list[1])
    axs.stairs(hist_draw_list[1].view().value, edges=hist_draw_list[1].axes[0].edges,
                label=label_list[1])
    # axs.plot(hist_draw_list[2].axes[0].centers, hist_draw_list[2].view().value,
    #             label=label_list[2])
    axs.fill_between(x_fine, lower_env, upper_env, color='blue', alpha=0.3, label='Envelope')
    ymax = max([max(hist.view().value) for hist in hist_draw_list])
    axs.set_ylim((0, 1.15*ymax))
    axs.set_title("")
    axs.set_xlabel("$M_{e \\mu} (GeV)$")
    axs.set_ylabel("Reweight Factor (NLO/LO)")
    hist_handles, hist_labels = axs.get_legend_handles_labels()
    axs.legend(hist_handles, hist_labels, loc = "lower center")
    # Slightly fancy to remove whitespace
    fig.savefig('KFactorwScaleVariations.png')
    plt.close()

    # Plot reweight histogram with errors
    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    axs.bar(
        hist_draw_list[1].axes[0].centers,
        hist_draw_list[1].view().value,
        width = hist_draw_list[1].axes[0].widths,
        yerr = np.sqrt(hist_draw_list[1].view().variance),
        fill = False,
        )
    ymax = max(hist_draw_list[1].view().value)
    axs.set_ylim((0, 1.15*ymax))
    axs.set_title("")
    axs.set_xlabel("$M_{e \\mu} (GeV)$")
    axs.set_ylabel("Reweight Factor (NLO/LO)")
    # Slightly fancy to remove whitespace
    fig.savefig('KFactorwErrors.png')
    plt.close()
    logging.debug(f"Finished Drawing {hist_tag}")
