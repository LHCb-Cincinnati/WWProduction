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
import iminuit
from iminuit import minimize  # has same interface as scipy.optimize.minimize
from iminuit import Minuit, describe
from iminuit.cost import LeastSquares
# Personal Packages 
import AnalysisTools as at

# Flat function
def flat(x, a):
    return(0*x + a)

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
Z1_lminus_vec = vector.zip({
    'px': tree['Z1_lminus'].px,
    'py': tree['Z1_lminus'].py,
    'pz': tree['Z1_lminus'].pz,
    'e': tree['Z1_lminus'].e,
    'pid': tree['Z1_lminus'].pid
})

Z1_lplus_vec = vector.zip({
    'px': tree['Z1_lplus'].px,
    'py': tree['Z1_lplus'].py,
    'pz': tree['Z1_lplus'].pz,
    'e': tree['Z1_lplus'].e,
    'pid': tree['Z1_lplus'].pid
})
Z2_lminus_vec = vector.zip({
    'px': tree['Z2_lminus'].px,
    'py': tree['Z2_lminus'].py,
    'pz': tree['Z2_lminus'].pz,
    'e': tree['Z2_lminus'].e,
    'pid': tree['Z2_lminus'].pid
})

Z2_lplus_vec = vector.zip({
    'px': tree['Z2_lplus'].px,
    'py': tree['Z2_lplus'].py,
    'pz': tree['Z2_lplus'].pz,
    'e': tree['Z2_lplus'].e,
    'pid': tree['Z2_lplus'].pid
})

# Masks
# Z1_lminus and Z2_lplus pairs
Z1lmZ2lp_mue_decay_mask = (
    (abs(Z1_lminus_vec.pid)!=abs(Z2_lplus_vec.pid))
    & (Z1_lminus_vec.pid*Z2_lplus_vec.pid < 0)
)
Z1lmZ2lp_electron_vec = ak.where((abs(Z1_lminus_vec.pid)==11), Z1_lminus_vec, Z2_lplus_vec)
Z1lmZ2lp_muon_vec = ak.where((abs(Z1_lminus_vec.pid)==13), Z1_lminus_vec, Z2_lplus_vec)
Z1lmZ2lp_tight_acc_mask = (
    (Z1_lminus_vec.eta>2.2)
    & (Z1_lminus_vec.eta<4.4)
    & (Z2_lplus_vec.eta>2.2)
    & (Z2_lplus_vec.eta<4.4)
)
Z1lmZ2lp_high_pT_mask = (
    (Z1lmZ2lp_electron_vec.pt>30)
    & (Z1lmZ2lp_muon_vec.pt>25)
)
Z1lmZ2lp_masks = (
    Z1lmZ2lp_mue_decay_mask
    & Z1lmZ2lp_tight_acc_mask
    & Z1lmZ2lp_high_pT_mask
)
# Z2_lminus and Z1_lplus pairs
Z2lmZ1lp_mue_decay_mask = (
    (abs(Z2_lminus_vec.pid)!=abs(Z1_lplus_vec.pid))
    & (Z2_lminus_vec.pid*Z1_lplus_vec.pid < 0)
)
Z2lmZ1lp_electron_vec = ak.where((abs(Z2_lminus_vec.pid)==11), Z2_lminus_vec, Z1_lplus_vec)
Z2lmZ1lp_muon_vec = ak.where((abs(Z2_lminus_vec.pid)==13), Z2_lminus_vec, Z1_lplus_vec)
Z2lmZ1lp_tight_acc_mask = (
    (Z2_lminus_vec.eta>2.2)
    & (Z2_lminus_vec.eta<4.4)
    & (Z1_lplus_vec.eta>2.2)
    & (Z1_lplus_vec.eta<4.4)
)
Z2lmZ1lp_high_pT_mask = (
    (Z2lmZ1lp_electron_vec.pt>30)
    & (Z2lmZ1lp_muon_vec.pt>25)
)
Z2lmZ1lp_masks = (
    Z2lmZ1lp_mue_decay_mask
    & Z2lmZ1lp_tight_acc_mask
    & Z2lmZ1lp_high_pT_mask
)
# Apply Masks
z2lmz1lp_dilepton_vec = Z2_lminus_vec[Z2lmZ1lp_masks] + Z1_lplus_vec[Z2lmZ1lp_masks]
z2lpz1lm_dilepton_vec = Z2_lplus_vec[Z1lmZ2lp_masks] + Z1_lminus_vec[Z1lmZ2lp_masks]

# Print Info
print(f"Total Cuts: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_rghbin_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_hist.fill(z2lmz1lp_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_hist.fill(z2lpz1lm_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_profilehist.fill(
    z2lpz1lm_dilepton_vec.m, weight=scale_factor, sample=z2lpz1lm_dilepton_vec.m
)
dilepton_id_mass_kfactorbin_profilehist.fill(
    z2lmz1lp_dilepton_vec.m, weight=scale_factor, sample=z2lmz1lp_dilepton_vec.m
)
if weights_bool:
    for weight_name in tree["Weights"].fields:
        weighthist_dict[weight_name].fill(
            z2lmz1lp_dilepton_vec.m, weight=tree["Weights"][weight_name][Z2lmZ1lp_masks]
        )
        weighthist_dict[weight_name].fill(
            z2lpz1lm_dilepton_vec.m, weight=tree["Weights"][weight_name][Z1lmZ2lp_masks]
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
print(f"Unweighted Events: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks)}")
print(f"Weighted Events: {sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks) * scale_factor}")
if weights_bool:
    print(f"Lower Bound: {sum(lower_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks) * scale_factor}")
    print(f"Upper Bound: {sum(upper_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * sum(Z2lmZ1lp_masks | Z1lmZ2lp_masks) * scale_factor}")

# Divide Hists
    lower_env_hist = at.divide_bh_histograms(lower_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
    upper_env_hist = at.divide_bh_histograms(upper_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
    central_hist = at.divide_bh_histograms(weighthist_dict["AUX_MUR10_MUF10"], weighthist_dict["AUX_MUR10_MUF10"])

# Multiply Hists
#     lower_env_hist = at.multiply_bh_histograms(lower_env_hist, dilepton_id_mass_kfactorbin_hist)
#     upper_env_hist = at.multiply_bh_histograms(upper_env_hist, dilepton_id_mass_kfactorbin_hist)

# Fit Scale Variations
fit_func = flat
# Lower Envelope Fit
lowerenv_least_squares = LeastSquares(
    dilepton_id_mass_kfactorbin_profilehist.view().value, 
    lower_env_hist.view().value, 
    np.sqrt(lower_env_hist.view().variance), 
    fit_func
)
lowerenv_m = Minuit(lowerenv_least_squares, a=0.85)  
lowerenv_m.migrad()  # finds minimum of least_squares function
lowerenv_m.hesse()  # accurately computes uncertainties
print("Lower Envelope Fit Info:")
print(f"chi^2/n_dof = {lowerenv_m.fval:.1f} / {lowerenv_m.ndof:.0f} = {lowerenv_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(lowerenv_m.parameters, lowerenv_m.values, lowerenv_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")
# Upper Envelope Fit
upperenv_least_squares = LeastSquares(
    dilepton_id_mass_kfactorbin_profilehist.view().value, 
    upper_env_hist.view().value, 
    np.sqrt(upper_env_hist.view().variance), 
    fit_func
)
upperenv_m = Minuit(upperenv_least_squares, a=1.15)  
upperenv_m.migrad()  # finds minimum of least_squares function
upperenv_m.hesse()  # accurately computes uncertainties
print("upper Envelope Fit Info:")
print(f"chi^2/n_dof = {upperenv_m.fval:.1f} / {upperenv_m.ndof:.0f} = {upperenv_m.fmin.reduced_chi2:.1f}")
for p, v, e in zip(upperenv_m.parameters, upperenv_m.values, upperenv_m.errors):
    print(f"{p} = {v:.3f} +- {e:.3f}")

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
plt.axhline(y=1, color='black', linestyle='-', label="Central Value")
# axs.bar(
#     central_hist.axes[0].centers,
#     central_hist.view().value,
#     width = central_hist.axes[0].widths,
#     fill = False,
#     label = "Central Value"
# )
# axs.bar(
#     lower_env_hist.axes[0].centers,
#     lower_env_hist.view().value,
#     width = lower_env_hist.axes[0].widths,
#     fill = False,
# )
# axs.errorbar(
#     lower_env_hist.axes[0].centers,
#     lower_env_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     yerr = np.sqrt(lower_env_hist.view().variance),
#     label="Lower Envelope"
# )
# axs.bar(
#     upper_env_hist.axes[0].centers,
#     upper_env_hist.view().value,
#     width = upper_env_hist.axes[0].widths,
#     fill = False,
#     label = "Upper Env"
# )
# axs.errorbar(
#     upper_env_hist.axes[0].centers,
#     upper_env_hist.view().value,
#     ecolor = "black",
#     linestyle = "",
#     yerr = np.sqrt(upper_env_hist.view().variance),
#     label="Upper Envelope"
# )
axs.bar(
    upper_env_hist.axes[0].centers,
    height=(upper_env_hist.view().value - lower_env_hist.view().value), 
    bottom=lower_env_hist.view().value, 
    width = upper_env_hist.axes[0].widths,
    linewidth=0, 
    color="grey", 
    alpha=0.25, 
    # zorder=-1, 
    label="Scale Variation Envelope"
)
axs.set_xlim((0, 300))
axs.set_ylim((0.8, 1.2))
axs.set_title("")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("Scale Variation Ratio Compared to Nominal")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
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
