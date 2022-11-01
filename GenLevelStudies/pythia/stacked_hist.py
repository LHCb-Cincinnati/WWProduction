# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
from sklearn.decomposition import dict_learning
import uproot
import vector
sys.path.append('../../Analysis') # This is Awesome!
import analysis_tools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name =  args.input_files[0].name
cross_section = args.cross_section[0] # Cross section in fb

# Functions
def GetAnalysisArray(tree):
    lminus_vec = vector.zip({
        'px': tree['TargetLepton'].px,
        'py': tree['TargetLepton'].py,
        'pz': tree['TargetLepton'].pz,
        'e': tree['TargetLepton'].energy,
        'pid': tree['TargetLepton'].id
    })

    lplus_vec = vector.zip({
        'px': tree['TargetAntiLepton'].px,
        'py': tree['TargetAntiLepton'].py,
        'pz': tree['TargetAntiLepton'].pz,
        'e': tree['TargetAntiLepton'].energy,
        'pid': tree['TargetAntiLepton'].id
    })

    dilepton_vec = lminus_vec + lplus_vec
    muon_vec = ak.where((abs(lminus_vec.pid)==13), lminus_vec, lplus_vec)
    electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
    leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
    trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

    # Masks
    both_lepton_acc_mask = ((lminus_vec.eta>1.596) & (lplus_vec.eta>1.596))
    high_pT_lepton_mask = ((lminus_vec.pt>15) | (lplus_vec.pt>15))
    low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
    mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                        | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
    one_lepton_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
    invariant_mass_mask = (dilepton_vec.m>10)
    lepton_mask = both_lepton_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask
    drellyan_mask = both_lepton_acc_mask&low_pT_lepton_mask&invariant_mass_mask

    # Apply Masks
    dilepton_vec = dilepton_vec[lepton_mask]
    muon_vec = muon_vec[lepton_mask]
    electron_vec = electron_vec[lepton_mask]
    leading_lepton_vec = leading_lepton_vec[lepton_mask]
    trailing_lepton_vec = trailing_lepton_vec[lepton_mask]

    # Calculate Quantitites
    delta_phi_array = np.abs(muon_vec.deltaphi(electron_vec))
    delta_eta_array = np.abs(muon_vec.deltaeta(electron_vec))
    delta_r_array = np.abs(muon_vec.deltaR(electron_vec))

    array = ak.zip({'delta_phi': delta_phi_array,
                    'delta_eta': delta_eta_array,
                    'delta_r': delta_r_array,
                    'dilepton_vec': dilepton_vec,
                    'muon_vec': muon_vec,
                    'electron_vec': electron_vec,
                    'leading_lepton_vec': leading_lepton_vec,
                    'trailing_lepton_vec': trailing_lepton_vec})
    return(array)

# Open the file
data_list = [0]*len(args.input_files)
weight_list = [0]*len(args.input_files)
for index, file_name in enumerate(args.input_files):
    ifile = uproot.open(file_name.name)
    tree = ifile['Tree'].arrays()
    data_list[index] = GetAnalysisArray(tree)
    weight_list[index] = at.calculate_weights(data_list[index].delta_phi,
                                              args.cross_section[index])

# pdb.set_trace()
# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path("Stacked_Hist_WW_ttbar.root", False)
os.chdir(file_path)
label_list = ["WW", "ttbar"]

# Plots
at.create_stacked_hist([data_list[index].delta_eta for index in range(len(data_list))],
                        "Standalone Delta Eta",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 3),
                        label=label_list)
at.create_stacked_hist([data_list[index].delta_phi for index in range(len(data_list))],
                        "Standalone Delta Phi",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50,
                        label=label_list)
at.create_stacked_hist([data_list[index].delta_r for index in range(len(data_list))],
                        "Standalone Delta R",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 5),
                        label=label_list)
at.create_stacked_hist([data_list[index].dilepton_vec.m for index in range(len(data_list))],
                        "Standalone Dilepton Mass",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 500),
                        label=label_list)
at.create_stacked_hist([data_list[index].dilepton_vec.pt for index in range(len(data_list))],
                        "Standalone Dilepton pT Linear",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 500),
                        label=label_list)
at.create_stacked_hist([data_list[index].dilepton_vec.pt for index in range(len(data_list))],
                        "Standalone Dilepton pT Log",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 500), yscale='log',
                        label=label_list)
at.create_stacked_hist([data_list[index].muon_vec.pt for index in range(len(data_list))],
                        "Standalone Muon pT Linear",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 150),
                        label=label_list)
at.create_stacked_hist([data_list[index].muon_vec.pt for index in range(len(data_list))],
                        "Standalone Muon pT Log",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 150), yscale='log',
                        label=label_list)
at.create_stacked_hist([data_list[index].electron_vec.pt for index in range(len(data_list))],
                        "Standalone Electron pT Linear",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 150),
                        label=label_list)
at.create_stacked_hist([data_list[index].electron_vec.pt for index in range(len(data_list))],
                        "Standalone Electron pT Log",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(0, 150), yscale='log',
                        label=label_list)
at.create_stacked_hist([data_list[index].electron_vec.eta for index in range(len(data_list))],
                        "Standalone Electron Eta",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(1.5, 5),
                        label=label_list)
at.create_stacked_hist([data_list[index].muon_vec.eta for index in range(len(data_list))],
                        "Standalone Muon Eta",
                        weights=[weight_list[index] for index in range(len(data_list))],
                        bins=50, range=(1.5, 5),
                        label=label_list)