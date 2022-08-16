# Imports
import sys
import os

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector


# Functions
def create_folder_path(file_name):
    folder_name = file_name[:-5]
    folder_name = folder_name.split('/')[-1]
    path = os.environ['HOME'] + '/WWProduction/Data/Figures/' + folder_name
    if not os.path.exists(path):
        os.mkdir(path)
    return(path) 

def create_hist(array, title, yscale='linear', **kwargs):
    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    axs.hist(array, **kwargs)
    plt.yscale(yscale)
    plt.title(title)
    fig_string = (f"Statistics:\n"
                  f"Count: {np.size(array)}\n"
                  f"Mean:  {np.mean(array):.2f}\n"
                  f"Sigma: {np.std(array):.2f}")

    axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))

    # Slightly fancy to remove whitespace
    save_str = ''.join(title.split())
    plt.savefig(save_str + '.png')


# Open the file
file_name =  sys.argv[1]
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
file_path = create_folder_path(file_name)
os.chdir(file_path)

# Plots
create_hist(dilepton_vec.m, 'DiLepton Mass', bins=50, range=(0,150000))
create_hist(dilepton_vec.pt, 'Lepton Pair pT', yscale='log', bins=50, range=(0,150000))
create_hist(leading_lepton_pT_array, 'Leading Lepton pT', yscale='log', bins=50, range=(0,150000))
create_hist(trailing_lepton_pT_array, 'Trailing Lepton pT', yscale='log', bins=50, range=(0,150000))
create_hist(lminus_vec.pt, 'lminus pT', bins=50, yscale='log', range=(0,150000))
create_hist(lplus_vec.pt, 'lplus pT', yscale='log', bins=50, range=(0,150000))
create_hist(delta_phi_array, 'Delta Phi', bins=50)
create_hist(delta_r_array, 'Delta R', bins=50)
create_hist(delta_eta_array, 'Delta Eta', bins=50)
