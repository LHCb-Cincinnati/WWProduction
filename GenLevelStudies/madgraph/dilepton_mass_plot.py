# Imports
# STL Packages
import sys
import os
import pdb
import pickle
# Scikit Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages
from hepunits.units import MeV, GeV
import uproot
import awkward as ak
import vector
import boost_histogram as bh
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
cross_section = args.cross_section # Cross section in fb
luminosity = args.luminosity # Luminosity in fb^-1
file_list = [file.name for file in args.input_files]
ofile_name = args.output

# Histograms
dilepton_mass_hist = bh.Histogram(bh.axis.Regular(60, 0, 300), storage=bh.storage.Weight())

# Iterator over large files
for tree in uproot.iterate(file_list[0] + ":Tree"):
    # Create Vectors
    lminus_vec = vector.zip({
        'px': tree['TargetLepton'].px,
        'py': tree['TargetLepton'].py,
        'pz': tree['TargetLepton'].pz,
        'e': tree['TargetLepton'].e,
        'pid': tree['TargetLepton'].pid
    })

    lplus_vec = vector.zip({
        'px': tree['TargetAntiLepton'].px,
        'py': tree['TargetAntiLepton'].py,
        'pz': tree['TargetAntiLepton'].pz,
        'e': tree['TargetAntiLepton'].e,
        'pid': tree['TargetAntiLepton'].pid
    })
    dilepton_vec = lminus_vec + lplus_vec
    muon_vec = ak.where((abs(lminus_vec.pid)==13), lminus_vec, lplus_vec)
    electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
    leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
    trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

    # Masks
    both_lepton_loose_acc_mask = ((lminus_vec.eta>2)
                                    & (lminus_vec.eta<5)
                                    & (lplus_vec.eta>2)
                                    & (lplus_vec.eta<5)) 
    both_lepton_tight_acc_mask = ((lminus_vec.eta>2.2)
                                    & (lminus_vec.eta<4.4)
                                    & (lplus_vec.eta>2.2)
                                    & (lplus_vec.eta<4.4)) 
    high_pT_lepton_mask = ((lminus_vec.pt>20) & (lplus_vec.pt>20))
    low_pT_lepton_mask = ((lminus_vec.pt>5) & (lplus_vec.pt>5))
    mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                        | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
    one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
    one_lepton_tight_acc_mask = ((lminus_vec.eta>2) | (lplus_vec.eta>2))
    one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                            | (lplus_vec.eta>1.596) & (lplus_vec.pt>15))
    invariant_mass_mask = (dilepton_vec.m>10)
    lepton_mask = both_lepton_tight_acc_mask&high_pT_lepton_mask&low_pT_lepton_mask&mue_decay_mask

    # Apply Masks
    muon_vec = muon_vec[lepton_mask]
    electron_vec = electron_vec[lepton_mask]
    dilepton_vec = dilepton_vec[lepton_mask]

    # Fill Histograms
    dilepton_mass_hist.fill(dilepton_vec.m)

    print(f"Total number of events: {len(lminus_vec)}")
    print(f"Gauss Cuts: {sum((one_lepton_gauss_mask))}")
    print(f"DaVinci Cut: {sum((low_pT_lepton_mask))}")
    print(f"Loose Acceptance Cuts: {sum(both_lepton_loose_acc_mask&mue_decay_mask)}")
    print(f"Tight Acceptance Cuts: {sum(lepton_mask)}")

# Scale histograms
if args.cross_section:
    dilepton_mass_hist*=args.cross_section

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# lepton pT Hist
at.create_stair(
    dilepton_mass_hist,
    "Dilepton Mass",
    luminosity=luminosity
)

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMass": [dilepton_mass_hist],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)