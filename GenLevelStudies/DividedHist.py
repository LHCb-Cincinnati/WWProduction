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
import matplotlib
import matplotlib.pyplot as plt
# HEP Imports
import ROOT
from hepunits.units import MeV, GeV
import awkward as ak
import uproot
import vector
import boost_histogram as bh
# Personal Imports
import AnalysisTools as at

# # Start Logger
# logging.basicConfig(
#     filename='StackedHist.log',
#     level=logging.DEBUG,
#     format='%(asctime)s %(levelname)s %(message)s',
#     datefmt='%m/%d/%Y %I:%M:%S %p'
# )
# plt.set_loglevel (level = 'warning')
# # get the the logger with the name 'PIL'
# pil_logger = logging.getLogger('PIL')  
# # override the logger logging level to INFO
# pil_logger.setLevel(logging.INFO)

# Parse inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
cross_section = args.cross_section # Cross section in fb
luminosity = args.luminosity # Luminosity in fb^-1
file_list = [file.name for file in args.input_files]
ofile_name = args.output
# logging.info(f"Arguments: {args}")

# Create ROOT File
folder_path = at.create_folder_path(ofile_name, args.testing)
root_ofile = ROOT.TFile(at.find_WW_path() + "/GenLevelStudies/Histograms/" + ofile_name + ".root", "RECREATE")

# Define Variables
data_list = [0] * len(file_list)
hist_title = "Weights WW Pythia LO to Madgraph NLO"

# Retrieve Histos
# Here data are the histograms and index is the index of the file in the
# list of files. So data_list looks like:
# data_list[1] = dictionary first file's histograms
# The dictionary looks like
# dict[Histogram Tag] = Histogram
for index, file in enumerate(file_list):
    with open(file, 'rb') as f:
        data = pickle.load(f)
        # logging.debug(f"Keys taken from {index} file: {data.keys()}")
        data_list[index] = data

# Save Figures
os.chdir(folder_path)
for hist_tag in data_list[0].keys():
    # Add histogram to hist_draw_list from first file.
    # logging.debug(f"Drawing {hist_tag}")
    hist_draw_list = data_list[0][hist_tag]
    # logging.debug(f"Added {data_list[0][hist_tag]} to hist_draw_list")
    # Loop through other files and add other histograms with the
    # same tag to the hist_draw_list
    for file in data_list[1:]:
        # Check if tag is in histogram
        if hist_tag not in file.keys():
            raise RuntimeWarning(f"Tag: {hist_tag} not found in file: {file}.")
        hist_draw_list+=file[hist_tag]
        # logging.debug(f"Added {file[hist_tag]} to hist_draw_list")
    print(hist_tag)
    div_hist = at.divide_bh_histograms(hist_draw_list[0], hist_draw_list[1],binomial_error=False)        
    # Draw the histogram.  Only the first histograms corresponding to the
    # number of labels given are drawn.
    if len(div_hist.axes) == 1:
        at.create_stair(
            div_hist,
            hist_title,
        )
        at.create_stair(
            div_hist,
            hist_title + " Log",
            yscale="log"
        )
        # Create ROOT hist
        root_hist = ROOT.TH1F(
            "rwgt_hist",
            "rwgt_hist",
            len(div_hist.axes[0].edges) - 1,
            div_hist.axes[0].edges
        )
        for index in range(len(div_hist.axes[0].edges) - 1):
            root_hist.SetBinContent(index + 1, div_hist.view().value[index])
            root_hist.SetBinError(index + 1, div_hist.view().variance[index])
        logging.debug(f"Finished Drawing {hist_tag}")
        root_hist.Write()
    else:
        fig, axs = plt.subplots()
        # Messing with Color Maps
        norm=matplotlib.colors.TwoSlopeNorm(vmin=0.5, vcenter=1, vmax=2)
        cmap = matplotlib.cm.get_cmap("bwr").copy()
        # cmap.Normalize(vim = 0.5, vmax = 2, clip = False)
        colormesh = plt.pcolormesh(
            div_hist.axes[0].edges,
            div_hist.axes[1].edges,
            div_hist.view().value,
            cmap=cmap,
            norm=norm
        )
        axs.set_xlabel("$\eta_e$")
        axs.set_ylabel("$\eta_{\mu}$")
        colorbar = plt.colorbar(colormesh)
        colorbar.ax.set_ylabel("Ratio")
        cmap.set_under("grey")
        # cmap.set_over("grey")
        # cmap.set_bad("white")
        plt.savefig("Test" + ".png")

    # Save histograms
    os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
    pickle_dict = {
        "DiLeptonMass": [div_hist],
    }
    with open(ofile_name + ".pkl", "wb") as f:
        pickle.dump(pickle_dict, f)
    
    # Close ROOT File
    root_ofile.Close()
