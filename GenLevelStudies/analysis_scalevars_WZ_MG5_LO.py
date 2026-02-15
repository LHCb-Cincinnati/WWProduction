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
lminus_vec = vector.zip({
    'px': tree['Z_lminus'].px,
    'py': tree['Z_lminus'].py,
    'pz': tree['Z_lminus'].pz,
    'e': tree['Z_lminus'].e,
    'pid': tree['Z_lminus'].pid
})

lplus_vec = vector.zip({
    'px': tree['Z_lplus'].px,
    'py': tree['Z_lplus'].py,
    'pz': tree['Z_lplus'].pz,
    'e': tree['Z_lplus'].e,
    'pid': tree['Z_lplus'].pid
})

thirdl_vec = vector.zip({
    'px': tree['W_lepton'].px,
    'py': tree['W_lepton'].py,
    'pz': tree['W_lepton'].pz,
    'e': tree['W_lepton'].e,
    'pid': tree['W_lepton'].pid
})

# Masks
# general masks
one_lepton_gauss_mask = ((lminus_vec.eta>1.596) & (lminus_vec.pt>15)
                         | (lplus_vec.eta>1.596) & (lplus_vec.pt>15)
                         | (thirdl_vec.eta>1.596) & (thirdl_vec.pt>15))
# lminus and wlepton pair masks
lmwl_mue_decay_mask = (
    (abs(lminus_vec.pid)!=abs(thirdl_vec.pid))
    & (lminus_vec.pid*thirdl_vec.pid < 0)
)
lmwl_tight_acc_mask = (
    (lminus_vec.eta>2.2)
    & (lminus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lmwl_high_pT_mask = (
    (lminus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lmwl_masks = (
    lmwl_mue_decay_mask
    & lmwl_tight_acc_mask
    & lmwl_high_pT_mask
)
# lplus and wlepton pair masks
lpwl_mue_decay_mask = (
    (abs(lplus_vec.pid)!=abs(thirdl_vec.pid))
    & (lplus_vec.pid*thirdl_vec.pid < 0)
)
lpwl_tight_acc_mask = (
    (lplus_vec.eta>2.2)
    & (lplus_vec.eta<4.4)
    & (thirdl_vec.eta>2.2)
    & (thirdl_vec.eta<4.4)
)
lpwl_high_pT_mask = (
    (lplus_vec.pt>20)
    & (thirdl_vec.pt>20)
)
lpwl_masks = (
    lpwl_mue_decay_mask
    & lpwl_tight_acc_mask
    & lpwl_high_pT_mask
)
# Apply Masks
lpwl_dilepton_vec = lplus_vec[lpwl_masks] + thirdl_vec[lpwl_masks]
lmwl_dilepton_vec = lplus_vec[lmwl_masks] + thirdl_vec[lmwl_masks]

# Print Info
print(f"Total Cuts: {sum(lpwl_masks | lmwl_masks)}")

# Fill histograms
dilepton_id_mass_rghbin_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_rghbin_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_finebin_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_hist.fill(lpwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_hist.fill(lmwl_dilepton_vec.m, weight=scale_factor)
dilepton_id_mass_kfactorbin_profilehist.fill(
    lmwl_dilepton_vec.m, weight=scale_factor, sample=lmwl_dilepton_vec.m
)
dilepton_id_mass_kfactorbin_profilehist.fill(
    lpwl_dilepton_vec.m, weight=scale_factor, sample=lpwl_dilepton_vec.m
)
if weights_bool:
    for weight_name in tree["Weights"].fields:
        weighthist_dict[weight_name].fill(
            lpwl_dilepton_vec.m, weight=tree["Weights"][weight_name][lpwl_masks]
        )
        weighthist_dict[weight_name].fill(
            lmwl_dilepton_vec.m, weight=tree["Weights"][weight_name][lmwl_masks]
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
print(f"Unweighted Events: {sum(lpwl_masks | lmwl_masks)}")
print(f"Weighted Events: {sum(lpwl_masks | lmwl_masks) * scale_factor}")
if weights_bool:
    print(f"Lower Bound: {sum(lower_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * sum(lpwl_masks | lmwl_masks) * scale_factor}")
    print(f"Upper Bound: {sum(upper_env_hist.view().value) / sum(weighthist_dict['AUX_MUR10_MUF10'].view().value) * sum(lpwl_masks | lmwl_masks) * scale_factor}")

# Divide Hists
    lower_env_hist = at.divide_bh_histograms(lower_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
    upper_env_hist = at.divide_bh_histograms(upper_env_hist, weighthist_dict["AUX_MUR10_MUF10"])
    central_hist = at.divide_bh_histograms(weighthist_dict["AUX_MUR10_MUF10"], weighthist_dict["AUX_MUR10_MUF10"])

# Multiply Hists
#     lower_env_hist = at.multiply_bh_histograms(lower_env_hist, dilepton_id_mass_kfactorbin_hist)
#     upper_env_hist = at.multiply_bh_histograms(upper_env_hist, dilepton_id_mass_kfactorbin_hist)

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

# Envelope Calculation
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
