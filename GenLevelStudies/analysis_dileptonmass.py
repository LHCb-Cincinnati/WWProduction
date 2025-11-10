# Imports
import sys
import os
import pdb
import pickle
# Scikit Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages
import boost_histogram as bh
import awkward as ak
import uproot
import vector
from hepunits.units import MeV, GeV
# Personal Packages 
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])
weights_bool = True
args = parser.args
file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section
ofile_name = args.output
luminosity = args.luminosity # Luminosity in fb^-1

# Get Tree
with uproot.open(file_name_list[0]) as f:
    # Get the Tree
    tree = f['Tree'].arrays()

# Scale Factor
if args.cross_section == False:
    scale_factor = 1
else:
    scale_factor = args.cross_section

# DFDY invariant mass cuts
invmass_cut = 0

# Binning Scheme
ttbar_bins_list = list(np.linspace(0, 100, 4)) + [150, 200, 2000]
DFDY_bins_list = list(np.linspace(0, 100, 8)) + [200, 2000]
bins_list = ttbar_bins_list

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_kfactorbin_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
eta_hist = bh.Histogram(
    bh.axis.Regular(12,-6, 6),
    bh.axis.Regular(12,-6, 6),
    storage=bh.storage.Weight()
)
if weights_bool:
    weight_name_list = ["AUX_MUR05_MUF05", "AUX_MUR05_MUF10", "AUX_MUR05_MUF20", "AUX_MUR10_MUF05", "AUX_MUR10_MUF10", "AUX_MUR10_MUF20", "AUX_MUR20_MUF05", "AUX_MUR20_MUF10", "AUX_MUR20_MUF20"]
    weighthist_dict = {}
    lower_env_hist = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
    upper_env_hist = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
    for weight_name in weight_name_list:
        weighthist_dict[weight_name] = bh.Histogram(
            bh.axis.Variable(bins_list), storage=bh.storage.Weight()
        )

# Create Vectors
lminus_vec = vector.zip({
    'px': tree['TargetLepton'].px * GeV,
    'py': tree['TargetLepton'].py * GeV,
    'pz': tree['TargetLepton'].pz * GeV,
    'e': tree['TargetLepton'].e * GeV,
    'pid': tree['TargetLepton'].pid
})
lplus_vec = vector.zip({
    'px': tree['TargetAntiLepton'].px * GeV,
    'py': tree['TargetAntiLepton'].py * GeV,
    'pz': tree['TargetAntiLepton'].pz * GeV,
    'e': tree['TargetAntiLepton'].e * GeV,
    'pid': tree['TargetAntiLepton'].pid
})
tau_vec = vector.zip({
    'px': tree['TargetParticle'].px * GeV,
    'py': tree['TargetParticle'].py * GeV,
    'pz': tree['TargetParticle'].pz * GeV,
    'e': tree['TargetParticle'].e * GeV,
    'pid': tree['TargetParticle'].pid
})
antitau_vec = vector.zip({
    'px': tree['TargetAntiParticle'].px * GeV,
    'py': tree['TargetAntiParticle'].py * GeV,
    'pz': tree['TargetAntiParticle'].pz * GeV,
    'e': tree['TargetAntiParticle'].e * GeV,
    'pid': tree['TargetAntiParticle'].pid
})
dilepton_vec = lminus_vec + lplus_vec
ditau_vec = tau_vec + antitau_vec
muon_vec = ak.where((abs(lminus_vec.pid)==13), lminus_vec, lplus_vec)
electron_vec = ak.where((abs(lminus_vec.pid)==11), lminus_vec, lplus_vec)
leading_lepton_vec = ak.where((lplus_vec.pt>lminus_vec.pt), lplus_vec, lminus_vec)
trailing_lepton_vec = ak.where((lplus_vec.pt<lminus_vec.pt), lplus_vec, lminus_vec)

# Masks
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt / GeV >15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt / GeV >15))
both_lepton_loose_acc_mask = ((lminus_vec.eta>2)
                                & (lminus_vec.eta<5)
                                & (lplus_vec.eta>2)
                                & (lplus_vec.eta<5)) 
both_lepton_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
) 
gauss_both_lepton_tight_acc_mask = (
    (lminus_vec[one_lepton_gauss_mask].eta>2.2)
    & (lminus_vec[one_lepton_gauss_mask].eta<4.4)
    & (lplus_vec[one_lepton_gauss_mask].eta>2.2)
    & (lplus_vec[one_lepton_gauss_mask].eta<4.4)
) 
mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                    | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
gauss_mue_decay_mask = (((lminus_vec[one_lepton_gauss_mask].pid==13) & (lplus_vec[one_lepton_gauss_mask].pid==-11))
                    | ((lminus_vec[one_lepton_gauss_mask].pid==11) & (lplus_vec[one_lepton_gauss_mask].pid==-13)))
one_lepton_loose_acc_mask = ((lminus_vec.eta>1.596) | (lplus_vec.eta>1.596))
one_lepton_tight_acc_mask = ((lminus_vec.eta>2) | (lplus_vec.eta>2))
invariant_mass_mask = (dilepton_vec.m / GeV >100)
deltar_mask = (np.abs(lminus_vec.deltaR(lplus_vec)) > 0.1)
gauss_deltar_mask = (np.abs(lminus_vec[one_lepton_gauss_mask].deltaR(lplus_vec[one_lepton_gauss_mask])) > 0.1)
gauss_high_pT_lepton_mask = ((lminus_vec[one_lepton_gauss_mask].pt / GeV >20) & (lplus_vec[one_lepton_gauss_mask].pt / GeV >20))
gauss_total_cuts = (
    gauss_high_pT_lepton_mask
    & gauss_both_lepton_tight_acc_mask
    & gauss_mue_decay_mask
    & gauss_deltar_mask
)
high_pT_lepton_mask = ((lminus_vec.pt / GeV >20) & (lplus_vec.pt / GeV >20))
low_pT_lepton_mask = ((lminus_vec.pt / GeV >5) & (lplus_vec.pt / GeV >5))
tautau_invmass_cut = (ditau_vec.m  / GeV < 500)
lepton_mask = (
    both_lepton_tight_acc_mask
    & high_pT_lepton_mask
    & mue_decay_mask
    & deltar_mask
    # & tautau_invmass_cut
)

# Apply Masks
# Note here that muon/electron vec and dilepton vec have different lengths.
muon_vec = muon_vec[mue_decay_mask & high_pT_lepton_mask & deltar_mask]
electron_vec = electron_vec[mue_decay_mask & high_pT_lepton_mask & deltar_mask]
dilepton_vec = dilepton_vec[lepton_mask]

# Print Info
print(f"Total Cuts: {sum(lepton_mask)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_kfactorbin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
eta_hist.fill(
    muon_vec.eta, electron_vec.eta, weight=scale_factor
)
if weights_bool:
    for weight_name in tree["Weights"].fields:
        weighthist_dict[weight_name].fill(
            dilepton_vec.m / GeV, weight=tree["Weights"][weight_name][lepton_mask]
        )

# Weight Hists
if weights_bool:
    for index in range(len(dilepton_id_mass_kfactorbin_hist.view().value)):
        min_val = weighthist_dict["AUX_MUR10_MUF10"].view()[index].value
        min_var = weighthist_dict["AUX_MUR10_MUF10"].view()[index].variance 
        max_val = weighthist_dict["AUX_MUR10_MUF10"].view()[index].value 
        max_var = weighthist_dict["AUX_MUR10_MUF10"].view()[index].variance
        for weight_name in tree["Weights"].fields:
            hist_view = weighthist_dict[weight_name].view()[index]
            if hist_view.value < min_val:
                min_val = hist_view.value
                min_var = hist_view.variance 
            if hist_view.value > max_val:
                max_val = hist_view.value
                max_var = hist_view.variance 
        lower_env_hist.view()[index] = [min_val, min_var]
        upper_env_hist.view()[index] = [max_val, max_var]

# Print Statements:
print(f"Unweighted Events: {len(dilepton_vec)}")
print(f"Weighted Events: {len(dilepton_vec) * scale_factor}")
if weights_bool:
    print(f"Lower Bound: {sum(lower_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * len(dilepton_vec) * scale_factor}")
    print(f"Upper Bound: {sum(upper_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * len(dilepton_vec) * scale_factor}")

# Divide Hists
    lower_env_hist = at.divide_bh_histograms(lower_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
    upper_env_hist = at.divide_bh_histograms(upper_env_hist, weighthist_dict["AUX_MUR10_MUF10"])

# Multiply Hists
    lower_env_hist = at.multiply_bh_histograms(lower_env_hist, dilepton_id_mass_kfactorbin_hist)
    upper_env_hist = at.multiply_bh_histograms(upper_env_hist, dilepton_id_mass_kfactorbin_hist)

if args.debug:
    pdb.set_trace()
    exit()

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Linear Rough Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_rghbin_hist, "DiLepton Mass Log Rough Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Linear Fine Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_finebin_hist, "DiLepton Mass Log Fine Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_kfactorbin_hist, "DiLepton Mass Linear K-Factor Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_kfactorbin_hist, "DiLepton Mass Log K-Factor Binning", yscale='log',
                luminosity=luminosity)
# Envelope Calculation
if weights_bool:
    at.create_stacked_stair(
        [lower_env_hist, dilepton_id_mass_kfactorbin_hist, upper_env_hist],
        "DiLeptonMassScaleVariations",
        ["Lower Envelope", "Central Value", "Upper Envelope"]
    )
# Eta Eta Yield 2D Hist
fig, axs = plt.subplots()
colormesh = plt.pcolormesh(
    eta_hist.axes[0].edges,
    eta_hist.axes[1].edges,
    eta_hist.view().value
)
colorbar = plt.colorbar(colormesh)
colorbar.ax.set_ylabel("Yield", rotation=-90)
axs.set_ylabel("$\eta_e$")
axs.set_xlabel("$\eta_{\mu}$")
save_str = ''.join("EtaEtaYieldHist".split())
plt.savefig(folder_path + "/" + save_str + '.png')
plt.close()

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    # "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    # "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    # "DileptonKFactorFine": [dilepton_id_mass_kfactorbin_hist],
    "DileptonKFactorFine": [upper_env_hist],
    "EtaEtaYield": [eta_hist]
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
