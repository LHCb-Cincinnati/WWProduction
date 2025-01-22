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
import mplhep as hep
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
dilepton_mass_hist = bh.Histogram(
    bh.axis.Variable(
        list(np.linspace(0, 100, 33)) + list(np.linspace(105, 200, 19)) + [220, 240, 260, 280, 300]),
    storage=bh.storage.Weight()
    )
eta_hist = bh.Histogram(
    bh.axis.Regular(12,-6, 6),
    bh.axis.Regular(12,-6, 6),
    storage=bh.storage.Weight()
)

# Counters
num_events = 0
num_tight_acceptance_cuts = 0

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
    lepton_mask = (
        high_pT_lepton_mask
        & both_lepton_tight_acc_mask
        & low_pT_lepton_mask
        & mue_decay_mask
    )

    # Apply Masks
    muon_vec = muon_vec[lepton_mask]
    electron_vec = electron_vec[lepton_mask]
    dilepton_vec = dilepton_vec[lepton_mask]

    # Fill Histograms
    dilepton_mass_hist.fill(dilepton_vec.m)
    eta_hist.fill(electron_vec.eta, muon_vec.eta)

    # Update Counters
    num_events += len(lminus_vec)
    num_tight_acceptance_cuts += sum(lepton_mask)

print(f"Total number of events: {num_events}")
print(f"Tight Acceptance Cuts: {num_tight_acceptance_cuts}")

# Scale histograms
if args.cross_section:
    dilepton_mass_hist*=args.cross_section
    eta_hist*=(args.cross_section * args.luminosity)

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# dilepton Mass Hist
plt.style.use(hep.style.LHCb2)  # or ATLAS/LHCb2
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
plt.stairs(
    dilepton_mass_hist.view().value,
    edges=dilepton_mass_hist.axes[0].edges
)
plt.xlabel("$M_{\mu e}$ (GeV)")
plt.ylabel("$ \\frac{d \sigma}{d M_{\mu e}}$ (GeV/fb)")
hist_count, hist_mean, hist_var = at.calculate_hist_stats(dilepton_mass_hist)
fig_string = (f"Statistics:\n"
                f"Count: {hist_count:.2f}\n"
                f"Mean:  {hist_mean:.2f}\n"
                f"Sigma: {hist_var:.2f}")
axs.text(0.8, 1.02, fig_string, transform=axs.transAxes,
            bbox=dict(facecolor='none', edgecolor='0.7', pad=3.0), fontsize=24)
# if luminosity:
#     axs.text(0, 1.01, "$\mathcal{L} = " + f"{luminosity}" + "fb^{-1}$",
#             transform=axs.transAxes, fontsize=12)
# Slightly fancy to remove whitespace
save_str = ''.join("DifferentialCrossSectioninDiLeptonMass".split())
plt.savefig(save_str + '.png')
plt.close()

# Eta Eta Yield 2D Hist
fig, axs = plt.subplots()
colormesh = plt.pcolormesh(
    eta_hist.axes[0].edges,
    eta_hist.axes[1].edges,
    eta_hist.view().value
)
colorbar = plt.colorbar(colormesh)
colorbar.ax.set_ylabel("Yield", rotation=-90)
axs.set_xlabel("$\eta_e$")
axs.set_ylabel("$\eta_{\mu}$")
save_str = ''.join("EtaEtaYieldHist".split())
plt.savefig(save_str + '.png')
plt.close()

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DifferentialCrossSectionMemu": [dilepton_mass_hist],
    "EtaEtaYield": [eta_hist],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
