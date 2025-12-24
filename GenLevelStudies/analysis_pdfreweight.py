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
DFDY_bins_list = [0, 33] + list(np.linspace(33, 100, 8))[1:] + [200, 2000]
bins_list = ttbar_bins_list 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_pdfreweight_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_pdfreweight_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
weight_name_list = tree["pdfReweight"].fields
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
    high_pT_lepton_mask
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
dilepton_id_mass_pdfreweight_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
dilepton_id_mass_pdfreweight_profilehist.fill(
    dilepton_vec.m / GeV, weight=scale_factor, sample=dilepton_vec.m / GeV
)
for weight_name in tree["pdfReweight"].fields:
    weighthist_dict[weight_name].fill(
        dilepton_vec.m / GeV, 
        weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
    )

# Print Statements:
print(f"Unweighted Events: {len(dilepton_vec)}")
print(f"Weighted Events: {len(dilepton_vec) * scale_factor}")
for weight in weight_name_list:
    print(f"{weight} Events: {weighthist_dict[weight_name].view().value.sum()}")

# Divide Hists
for weight_name in tree["pdfReweight"].fields:
    weighthist_idct[f"{weight_name}_Reweight"] = at.divide_bh_histograms(
        weighthist_idct[weight_name], dilepton_id_mass_pdfreweight_hist 
    )

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
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Linear pdfReweight Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Log pdfReweight Binning", yscale='log',
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_profilehist, "DiLepton Mass Linear K-Factor Binning Profile Hist",
                luminosity=luminosity)
# Special pdf plots
for weight_name in tree["pdfReweight"].fields:
    at.create_stacked_stair(
        [weighthist_dict[weight_name], dilepton_id_mass_pdfreweight_hist],
        f"DiLeptonMasspdfReweightBinning_{weight_name}",
        [weight_name, "CM09MTS"]
    )

# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    for weight_name in tree["pdfReweight"].fields:
        root_file[f"pdfReweight_{weight_name}"] =  weighthist_dict[f"{weight_name}_Reweight"]

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonpdfReweightFine": [dilepton_id_mass_pdfreweight_hist],
    "DileptonpdfReweightProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
# Store reweight histograms
for weight_name in tree["pdfReweight"].fields:
    f"pdfReweight_{weight_name}": [weighthist_dict[f"{weight_name}_Reweight"]],

with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
