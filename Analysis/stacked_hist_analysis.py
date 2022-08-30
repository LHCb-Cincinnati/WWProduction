# Imports
import sys
import os
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

import analysis_tools as at


# Parse inputs
args = at.parse_user_input(sys.argv)
if len(args.input_files) < 2:
    raise RuntimeError('Not enough input files provided.')

file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section # Cross section list in fb

# Create Indices
first_lepton_key = 'muon'
second_lepton_key = 'electron'

# Create Storage DFs
df_dict = {}
column_list = ('DiLepton Mass', 'Lepton Pair pT', 'Leading Lepton pT',
                'Trailing Lepton pT', 'Muon pT', 'Electron pT', 'Delta Phi',
                'Delta R', 'Delta Eta')
for file_name in file_name_list:
    df_dict[file_name] = pd.DataFrame(columns=column_list)


# Open each file and store relevant quantities in a DF
for file_name in file_name_list:
    # Open the file
    ifile = uproot.open(file_name)
    tree = ifile['Tuple/DecayTree'].arrays()

    # Create Vectors
    lminus_vec = vector.zip({
        'px': tree[first_lepton_key + '_PX'],
        'py': tree[first_lepton_key + '_PY'],
        'pz': tree[first_lepton_key + '_PZ'],
        'e': tree[first_lepton_key + '_PE'],
    })

    lplus_vec = vector.zip({
        'px': tree[second_lepton_key + '_PX'],
        'py': tree[second_lepton_key + '_PY'],
        'pz': tree[second_lepton_key + '_PZ'],
        'e': tree[second_lepton_key + '_PE'],
    })

    dilepton_vec = lminus_vec + lplus_vec

    # Calculate Quantitites
    greater_pt_array = (lminus_vec.pt > lplus_vec.pt)
    lesser_pt_array = -1*(greater_pt_array-1)
    leading_lepton_pT_array = (greater_pt_array * lminus_vec.pt
                               + lesser_pt_array * lplus_vec.pt)
    trailing_lepton_pT_array = (lesser_pt_array * lminus_vec.pt
                                + greater_pt_array * lplus_vec.pt)
    dilepton_mass_array = dilepton_vec.m
    dilepton_pT_array = dilepton_vec.pt
    muon_pT_array = lminus_vec.pt
    electron_pT_array = lplus_vec.pt
    delta_phi_array = np.abs(lminus_vec.deltaphi(lplus_vec))
    delta_eta_array = np.abs(lminus_vec.deltaeta(lplus_vec))
    delta_r_array = np.abs(lminus_vec.deltaR(lplus_vec))

#     # Fill Pandas DFs
    df_dict[file_name]['DiLepton Mass'] = dilepton_mass_array.to_numpy()
    df_dict[file_name]['Lepton Pair pT'] = dilepton_pT_array.to_numpy()
    df_dict[file_name]['Leading Lepton pT'] = leading_lepton_pT_array.to_numpy()
    df_dict[file_name]['Trailing Lepton pT'] = trailing_lepton_pT_array.to_numpy()
    df_dict[file_name]['Muon pT'] = muon_pT_array.to_numpy()
    df_dict[file_name]['Electron pT'] = electron_pT_array.to_numpy()
    df_dict[file_name]['Delta Phi'] = delta_phi_array.to_numpy()
    df_dict[file_name]['Delta R'] = delta_r_array.to_numpy()
    df_dict[file_name]['Delta Eta'] = delta_eta_array.to_numpy()

# Create output directory if it does not yet exist
# and change current directory to output directory
file_path = at.create_folder_path('StackedHist.root', args.testing)
os.chdir(file_path)

# Plots
at.create_stacked_hist((df_dict[file_name_list[0]]['DiLepton Mass'], 
                        df_dict[file_name_list[1]]['DiLepton Mass']),
                        'DiLepton Mass', bins=50, range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Lepton Pair pT'], 
                        df_dict[file_name_list[1]]['Lepton Pair pT']),
                        'Lepton Pair pT', bins=50, yscale='log',
                        range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Leading Lepton pT'], 
                        df_dict[file_name_list[1]]['Leading Lepton pT']),
                        'Leading Lepton pT', bins=50, yscale='log',
                        range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Trailing Lepton pT'], 
                        df_dict[file_name_list[1]]['Trailing Lepton pT']),
                        'Trailing Lepton pT', bins=50, yscale='log',
                        range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Muon pT'], 
                        df_dict[file_name_list[1]]['Muon pT']),
                        'Muon pT', bins=50, yscale='log',
                        range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Electron pT'], 
                        df_dict[file_name_list[1]]['Electron pT']),
                        'Electron pT', bins=50, yscale='log',
                        range=(0,150000),
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Delta Eta'], 
                        df_dict[file_name_list[1]]['Delta Eta']),
                        'Delta Eta', bins=50,
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Delta Phi'], 
                        df_dict[file_name_list[1]]['Delta Phi']),
                        'Delta Phi', bins=50,
                        label=(file_name_list[0], file_name_list[1]))
at.create_stacked_hist((df_dict[file_name_list[0]]['Delta R'], 
                        df_dict[file_name_list[1]]['Delta R']),
                        'Delta R', bins=50,
                        label=(file_name_list[0], file_name_list[1]))
