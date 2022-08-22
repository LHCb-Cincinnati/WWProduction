# Imports
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector


# Functions
def calculate_weights(cross_section, array):
    scale_factor = bool(cross_section)*(cross_section / len(array) - 1) + 1
    weights_array = [scale_factor] * len(array)
    return(weights_array)

def create_folder_path(file_name, test_mode_flag):
    if test_mode_flag:
        folder_name = 'Test'
    else:
        folder_name = file_name[:-5]
        folder_name = folder_name.split('/')[-1]
    path = os.environ['HOME'] + '/WWProduction/Data/Figures/' + folder_name
    if not os.path.exists(path):
        os.mkdir(path)
    return(path)

def calculate_hist_stats(hist, bins):
    ''' Calculate count, mean, and variance for a numpy histogram.

    Calculate count, mean, and variance of the given histogram using the
    histogram, not the data in the histogram.

    Args:
    hist (np.array): A numpy histogram that represents some dataset.
    bins (np.array): The bin edge placements for hist.

    Returns:
    hist_count (float): The number of samples in hist.
    hist_mean (float): The average value of hist.
    hist_variance (float): The variance of the samples in hist.

    '''

    hist_count = np.sum(hist)
    hist_mids = 0.5*(bins[1:] + bins[:-1])
    hist_mean = np.average(hist_mids, weights=(hist/hist_count))
    hist_var = np.sqrt(np.average((hist_mids - hist_mean)**2,
                       weights=(hist/hist_count)))
    return(hist_count, hist_mean, hist_var)

def create_hist(array, title, yscale='linear', **kwargs):
    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist, bins, patches = axs.hist(array, **kwargs)
    hist_count, hist_mean, hist_var = calculate_hist_stats(hist, bins)
    plt.yscale(yscale)
    plt.title(title)
    fig_string = (f"Statistics:\n"
                  f"Count: {hist_count:.2f}\n"
                  f"Mean:  {hist_mean:.2f}\n"
                  f"Sigma: {hist_var:.2f}")

    axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))

    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')

# Argument Parser
parser = argparse.ArgumentParser(description='Process files and settings for analysis.')
parser.add_argument('input_files',type=open, nargs='+',
                    help='The file or files to be anayzed. Input files should be root files.')
parser.add_argument('-c', '--cross_section', type=float, default=False, nargs='*', 
                    help='''The cross section used to normalize histograms in fb.
                    This option should either not be specified or have exactly
                    as many arguments as input files.''')
parser.add_argument('-t', '--testing', default=True, action='store_true',
                    help='''Flag to indicate if a new output folder should be created for 
                    this analysis.  -t means that all plots will go in folder labelled Test.''')
parser.add_argument('-p', '--production', dest='testing', action='store_false',
                    help='''Flag to indicate if a new output folder should be created for this 
                    analysis.  -p means that all plots will be put into a new output folder with
                    the same name as the input file.''')
args = parser.parse_args()
print(args)

# Make sure the cross-sections have as many aguments as files.
if (not args.cross_section):
    args.cross_section = [False] * len(args.input_files)
elif (len(args.cross_section) != len(args.input_files)):
    raise RuntimeError('''The number of given cross sections is different from
                       the number of given input files.  Please either don't
                       specify a cross section arugment or specify as many
                       cross sections as there are input files.''')

# Parse inputs
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tuple/DecayTree'].arrays()

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['lminus_PX'],
    'py': tree['lminus_PY'],
    'pz': tree['lminus_PZ'],
    'e': tree['lminus_PE'],
})

lplus_vec = vector.zip({
    'px': tree['lplus_PX'],
    'py': tree['lplus_PY'],
    'pz': tree['lplus_PZ'],
    'e': tree['lplus_PE'],
})
dilepton_vec = lminus_vec + lplus_vec


# Calculate Quantitites
greater_pt_array = (lminus_vec.pt > lplus_vec.pt)
lesser_pt_array = -1*(greater_pt_array-1)
leading_lepton_pT_array = greater_pt_array * lminus_vec.pt + lesser_pt_array * lplus_vec.pt
trailing_lepton_pT_array = lesser_pt_array * lminus_vec.pt + greater_pt_array * lplus_vec.pt
delta_phi_array = np.abs(lminus_vec.deltaphi(lplus_vec))
delta_eta_array = np.abs(lminus_vec.deltaeta(lplus_vec))
delta_r_array = np.abs(lminus_vec.deltaR(lplus_vec))

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = create_folder_path(file_name, args.testing)
os.chdir(file_path)

# Plots
create_hist(dilepton_vec.m, 'DiLepton Mass', bins=50, range=(0,150000),
            weights=calculate_weights(cross_section, dilepton_vec))
create_hist(dilepton_vec.pt, 'Lepton Pair pT', yscale='log', bins=50, range=(0,150000),
            weights=calculate_weights(cross_section, dilepton_vec))
create_hist(leading_lepton_pT_array, 'Leading Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=calculate_weights(cross_section, leading_lepton_pT_array))
create_hist(trailing_lepton_pT_array, 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150000),
            weights=calculate_weights(cross_section, trailing_lepton_pT_array))
create_hist(lminus_vec.pt, 'lminus pT', bins=50, yscale='log', range=(0,150000),
            weights=calculate_weights(cross_section, lminus_vec))
create_hist(lplus_vec.pt, 'lplus pT', yscale='log', bins=50, range=(0,150000),
            weights=calculate_weights(cross_section, lplus_vec))
create_hist(delta_phi_array, 'Delta Phi', bins=50,
            weights=calculate_weights(cross_section, delta_phi_array))
create_hist(delta_r_array, 'Delta R', bins=50,
            weights=calculate_weights(cross_section, delta_r_array))
create_hist(delta_eta_array, 'Delta Eta', bins=50,
            weights=calculate_weights(cross_section, delta_eta_array))