# Imports
# Standard Library Imports
import os
import pdb
import pickle
import sys
import logging
# Scikit Imports
import numpy as np
import matplotlib.pyplot as plt
# HEP Imports
from hepunits.units import MeV, GeV
import awkward as ak
import uproot
import vector
# Personal Imports
import AnalysisTools as at

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
label_list = ["No Ws", "MadSpin"]
logging.info(f"Label list provided: {label_list}")
bin_dict= {
    "LeadingLeptonpT": 50,
    "LeadingLeptonEta": 50,
    "TrailingLeptonpT": 50,
    "TrailingLeptonEta":50,
    "ElectronpT": 50,
    "ElectronEta": 50,
    "MuonpT": 50,
    "MuonEta": 50,
    "DiLeptonpT": 50,
    "DiLeptonMass": 50,
    "DeltaPhi": 50,
    "DeltaEta": 50,
    "DeltaR": 50
}
range_dict= {
    "LeadingLeptonpT": (20, 150),
    "LeadingLeptonEta": (2.2, 4.4),
    "TrailingLeptonpT": (20, 150),
    "TrailingLeptonEta": (2.2, 4.4),
    "ElectronpT": (20, 150),
    "ElectronEta": (2.2, 4.4),
    "MuonpT": (20, 150),
    "MuonEta": (2.2, 4.4),
    "DiLeptonpT": (20, 150),
    "DiLeptonMass": (20, 300),
    "DeltaPhi": (0, 3.14),
    "DeltaEta": (0, 3),
    "DeltaR": (0, 5)
}

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
    # Draw the histogram.  Only the first histograms corresponding to the
    # number of labels given are drawn.
    at.create_stacked_stair(
        hist_draw_list[:len(label_list)],
        hist_tag,
        label_list,
        normalize=True,
        bins=bin_dict[hist_tag],
        range=range_dict[hist_tag] 
    )
    # at.create_stacked_stair(hist_draw_list[:len(label_list)], hist_tag + "_Log", label_list,
    #                         yscale="log", luminosity=luminosity, normalize=True)
    logging.debug(f"Finished Drawing {hist_tag}")
