# Imports
# STL Packages
import sys
import os
import pdb
# Scikit Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages
import uproot
import awkward as ak
import vector
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
ifile_name =  args.input_files[0].name
cross_section = args.cross_section # Cross section in fb

# Open Files
ifile = uproot.open(ifile_name)
tree = ifile["LHEF"]

# Collect Initial Arrays
""" Indices 5 and 7 are those of the neutrino.  Indicies
    4 and 6 are those of the leptons.
"""
lepton_index = 6
antilepton_index = 3
neutrino_index = 4
antineutrino_index = 7
cut_string = f"(Particle.Eta[:,{lepton_index}] > 0) & (Particle.Eta[:,{antilepton_index}] > 0)"

# Creating Arrays
lepton_array = vector.array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, lepton_index],
                          'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, lepton_index],
                          'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, lepton_index],
                          'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, lepton_index]})

antilepton_array = vector.array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, antilepton_index],
                          'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, antilepton_index],
                          'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, antilepton_index],
                          'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, antilepton_index]})

neutrino_array = vector.array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, neutrino_index],
                          'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, neutrino_index],
                          'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, neutrino_index],
                          'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, neutrino_index]})

antineutrino_array = vector.array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, antineutrino_index],
                          'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, antineutrino_index],
                          'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, antineutrino_index],
                          'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, antineutrino_index]})

# Cuts
lepton_cut_mask = (lepton_array.eta>2) & (lepton_array.eta<5)
antilepton_cut_mask = (antilepton_array.eta>2) & (antilepton_array.eta<5)
cut_mask = lepton_cut_mask & antilepton_cut_mask

# Redefining arrays with cuts
lepton_array = lepton_array[cut_mask]
antilepton_array = antilepton_array[cut_mask]
neutrino_array = neutrino_array[cut_mask]
antineutrino_array = antineutrino_array[cut_mask]

# Calculate Additional Arrays
lepton_pair_array = lepton_array + antilepton_array
missing_array = neutrino_array + antineutrino_array

# Calculate Leading and Trailing Lepton arrays
# Gabe trying to be creative for no apparent reason.
greater_pt_array = (lepton_array.pt > antilepton_array.pt)
lesser_pt_array = -1*(greater_pt_array-1)
leading_lepton_pT_array = (greater_pt_array * lepton_array.pt
                           + lesser_pt_array * antilepton_array.pt)
trailing_lepton_pT_array = (lesser_pt_array * lepton_array.pt
                            + greater_pt_array * antilepton_array.pt)

# Calculate delta phi
closer_deltaphi_array = (np.abs(missing_array.deltaphi(lepton_array))
                         < np.abs(missing_array.deltaphi(antilepton_array)))
farther_deltaphi_array = -1*(closer_deltaphi_array-1)
deltaphi_array = np.abs(missing_array.deltaphi(lepton_array) * closer_deltaphi_array
                  + missing_array.deltaphi(antilepton_array) * farther_deltaphi_array)
missing_pt_proj_array = missing_array.pt * np.sin(deltaphi_array)

# Set Up Plotting Director
path = at.create_folder_path(args.output, args.testing)
os.chdir(path)

# Make Plots
at.create_stair(leading_lepton_pT_array, "Test", luminosity=args.luminosity)

