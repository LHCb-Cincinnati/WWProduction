# Imports
from black import E
from pandas import cut

import uproot
import numpy as np
import awkward as ak
import vector
import matplotlib.pyplot as plt
import pylhe

# Functions
def hist_string(array):
    fig_string = (f"Statistics:\n"
                  f"Count: {np.size(array)}\n"
                  f"Mean:  {np.mean(array):.2f}\n"
                  f"Sigma: {np.std(array):.2f}")
    return(fig_string)

# Collect Initial Arrays
""" Indices 5 and 7 are those of the neutrino.  Indicies
    4 and 6 are those of the leptons.
"""

lepton_index = 6
antilepton_index = 3
neutrino_index = 4
antineutrino_index = 7

cut_string = f"(Particle.Eta[:,{lepton_index}] > 0) & (Particle.Eta[:,{antilepton_index}] > 0)"

# # ROOT SETUP
# # Open Files
# ifile = uproot.open("WW_DoubleBosonDecay.root")
# tree = ifile["LHEF"]
# ofile = uproot.recreate("test.root")

# # Creating Arrays
# lepton_array = vector.Array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, lepton_index],
#                           'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, lepton_index],
#                           'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, lepton_index],
#                           'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, lepton_index]})

# antilepton_array = vector.Array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, antilepton_index],
#                           'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, antilepton_index],
#                           'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, antilepton_index],
#                           'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, antilepton_index]})

# neutrino_array = vector.Array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, neutrino_index],
#                           'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, neutrino_index],
#                           'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, neutrino_index],
#                           'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, neutrino_index]})

# antineutrino_array = vector.Array({'px': tree.arrays("Particle.Px", cut_string, how=tuple)[0][:, antineutrino_index],
#                           'py': tree.arrays("Particle.Py", cut_string, how=tuple)[0][:, antineutrino_index],
#                           'pz': tree.arrays("Particle.Pz", cut_string, how=tuple)[0][:, antineutrino_index],
#                           'E': tree.arrays("Particle.E", cut_string, how=tuple)[0][:, antineutrino_index]})


# LHE SETUP
ifile_str = "WW_DoubleBosonDecay_NLO.lhe"
pylhe.register_awkward()
arr = pylhe.to_awkward(pylhe.read_lhe_with_attributes(ifile_str))
lepton_array = arr['particles']['vector'][:,lepton_index]
antilepton_array = arr['particles']['vector'][:,antilepton_index]
neutrino_array = arr['particles']['vector'][:,neutrino_index]
antineutrino_array = arr['particles']['vector'][:,antineutrino_index]

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

leading_lepton_pT_array = greater_pt_array * lepton_array.pt + lesser_pt_array * antilepton_array.pt
trailing_lepton_pT_array = lesser_pt_array * lepton_array.pt + greater_pt_array * antilepton_array.pt

# Calculate delta phi
closer_deltaphi_array = (np.abs(missing_array.deltaphi(lepton_array)) < np.abs(missing_array.deltaphi(antilepton_array)))
farther_deltaphi_array = -1*(closer_deltaphi_array-1)

deltaphi_array = np.abs(missing_array.deltaphi(lepton_array) * closer_deltaphi_array
                  + missing_array.deltaphi(antilepton_array) * farther_deltaphi_array)

missing_pt_proj_array = missing_array.pt * np.sin(deltaphi_array)

# # Making some plots
# fig, axs = plt.subplots(1, 2)
# #First Plot
# axs[0].hist(leading_lepton_pT_array, bins=100, range=(0,150))
# axs[0].text(0.5, 0.75, s=hist_string(leading_lepton_pT_array), transform=axs[0].transAxes,
#              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))
# axs[0].set_title("Leading Lepton pT")
# axs[0].set_xlabel("GeV/c")
# axs[0].set_ylabel("Events/Gev")
# # Second Plot
# axs[1].hist(trailing_lepton_pT_array, bins=100, range=(0,150))
# axs[1].text(0.5, 0.75, hist_string(trailing_lepton_pT_array), transform=axs[1].transAxes,
#              bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0))
# axs[1].set_title("Trailing Lepton pT")
# axs[1].set_xlabel("GeV/c")
# axs[1].set_ylabel("Events/Gev")
# plt.show()


# Create the Output Tree
# ofile["tree"] = {
#     "LeadingLeptonpT": leading_lepton_pT_array,
#     "TrailingLeptonpT": trailing_lepton_pT_array,
#     "LeptonPairET": lepton_pair_array.et,
#     "LeptonPairPt": lepton_pair_array.pt,
#     "LeptonPairInvariantMass": lepton_pair_array.mass,
#     "MissingPT": missing_array.pt,
#     "MissingET": missing_array.et,
#     "DeltaPhi": deltaphi_array,
#     "MissingPTProj": missing_pt_proj_array,
# }

# Close the files
# ifile.close()
# ofile.close()