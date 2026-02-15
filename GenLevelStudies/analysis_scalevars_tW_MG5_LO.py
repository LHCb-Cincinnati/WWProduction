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
bins_list = ttbar_bins_list 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_kfactorbin_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_kfactorbin_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
weight_name_list = ["AUX_MUR05_MUF05", "AUX_MUR05_MUF10", "AUX_MUR05_MUF20", "AUX_MUR10_MUF05", "AUX_MUR10_MUF10", "AUX_MUR10_MUF20", "AUX_MUR20_MUF05", "AUX_MUR20_MUF10", "AUX_MUR20_MUF20"]
weighthist_dict = {}
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
dilepton_id_mass_kfactorbin_profilehist.fill(
    dilepton_vec.m / GeV, weight=scale_factor, sample=dilepton_vec.m / GeV
)
for weight_name in tree["Weights"].fields:
    weighthist_dict[weight_name].fill(
        dilepton_vec.m / GeV, weight=tree["Weights"][weight_name][lepton_mask]
    )

# Print Statements:
print(f"Unweighted Events: {len(dilepton_vec)}")
print(f"Weighted Events: {len(dilepton_vec) * scale_factor}")
for weight_name in tree["pdfReweight"].fields:
print(f"Lower Bound: {sum(lower_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * len(dilepton_vec) * scale_factor}")
print(f"Upper Bound: {sum(upper_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * len(dilepton_vec) * scale_factor}")

# Divide Hists
lower_env_hist = at.divide_bh_histograms(lower_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
upper_env_hist = at.divide_bh_histograms(upper_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
central_hist = at.divide_bh_histograms(weighthist_dict["AUX_MUR10_MUF10"], weighthist_dict["AUX_MUR10_MUF10"])


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
at.create_stair(dilepton_id_mass_kfactorbin_profilehist, "DiLepton Mass Linear K-Factor Binning Profile Hist",
                luminosity=luminosity)
# Scale Variation Plot
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
axs.bar(
    central_hist.axes[0].centers,
    central_hist.view().value,
    width = central_hist.axes[0].widths,
    fill = False,
    label = "Central Value"
)
axs.bar(
    lower_env_hist.axes[0].centers,
    lower_env_hist.view().value,
    width = lower_env_hist.axes[0].widths,
    fill = False,
)
axs.errorbar(
    lower_env_hist.axes[0].centers,
    lower_env_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = np.sqrt(lower_env_hist.view().variance),
    label="Lower Envelope"
)
axs.bar(
    upper_env_hist.axes[0].centers,
    upper_env_hist.view().value,
    width = upper_env_hist.axes[0].widths,
    fill = False,
    label = "Upper Env"
)
axs.errorbar(
    upper_env_hist.axes[0].centers,
    upper_env_hist.view().value,
    ecolor = "black",
    linestyle = "",
    yerr = np.sqrt(upper_env_hist.view().variance),
    label="Upper Envelope"
)
axs.set_xlim((0, 200))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Scale Variation Ratio Compared to Nominal")
# hist_handles, hist_labels = axs.get_legend_handles_labels()
# axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassScaleVariations.png')
plt.close()

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_kfactorbin_hist
    # root_file["DileptonKFactorFine"] = upper_env_hist

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonKFactorFine": [dilepton_id_mass_kfactorbin_hist],
    # "DileptonKFactorFine": [upper_env_hist],
    "DileptonKFactorProfile": [dilepton_id_mass_kfactorbin_profilehist],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
