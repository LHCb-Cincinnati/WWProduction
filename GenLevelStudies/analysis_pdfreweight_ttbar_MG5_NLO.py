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

# Functions
def plot_pdf_rms(hist_dict, pdf_str):
    """ hist should be in the order upper, lower, central
    """

    fig, axs = plt.subplots()
    plt.subplots_adjust(top=0.85)
    hist = hist_dict["PDFMember0Weight"]
    axs.stairs(
        hist.view().value, 
        edges=hist.axes[0].edges,
        label = pdf_str.upper() + " Central Value"
    )
    axs.bar(
        hist.axes[0].centers,
        hist.view().value,
        width = hist.axes[0].widths,
        yerr = hist.view().variance,
        fill = False,
        label = pdf_str.upper() + " RMS Envelope"
    )
    axs.set_xlim((0, 300))
    # axs.set_ylim((0, 35))
    axs.set_title("PDF Variations")
    axs.set_xlabel("$M_{e \\mu} (GeV)$")
    axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}}$")
    fig.savefig(f'DiLeptonMass_PDFVariations_{pdf_str}.png')
    plt.close()


# Parse Inputs
parser = at.Parser(sys.argv[1:])
args = parser.args
file_name_list =  [file.name for file in args.input_files]
cross_section = args.cross_section
ofile_name = args.output
luminosity = args.luminosity # Luminosity in fb^-1

# Get Tree
tree_iterator = uproot.iterate(
    [file_name + ":Tree"for file_name in file_name_list],
    step_size=45000
)

# Scale Factor
if args.cross_section == False:
    scale_factor = 1
else:
    scale_factor = args.cross_section

# Variables
unweighted_event_counter = 0
unweighted_event_counter = 0
num_nnpdf31lo_members = 101
num_msht20lo_members = 61
num_nnpdf31nlo_members = 101
num_ct18nlo_members = 59
num_msht20nlo_members = 65

# Binning Scheme
bins_list = list(np.linspace(0, 100, 4)) + [150, 200, 2000] 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_pdfreweight_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_pdfreweight_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
weight_name_list = ["nnpdf31lo", "ct18lo", "msht20lo", "nnpdf31nlo", "ct18nlo", "msht20nlo"]
weighthist_dict = {}
for weight_name in weight_name_list:
    weighthist_dict[weight_name] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
nnpdf31lo_hist_dict = {}
for i in range(num_nnpdf31lo_members):
    nnpdf31lo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
msht20lo_hist_dict = {}
for i in range(num_msht20lo_members):
    msht20lo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
nnpdf31nlo_hist_dict = {}
for i in range(num_nnpdf31nlo_members):
    nnpdf31nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
ct18nlo_hist_dict = {}
for i in range(num_ct18nlo_members):
    ct18nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
msht20nlo_hist_dict = {}
for i in range(num_msht20nlo_members):
    msht20nlo_hist_dict[f"PDFMember{i}Weight"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )

# Loop over files in small steps
for tree in tree_iterator:
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
        # gauss_high_pT_lepton_mask
        gauss_both_lepton_tight_acc_mask
        & gauss_mue_decay_mask
        & gauss_deltar_mask
    )
    high_pT_lepton_mask = ((lminus_vec.pt / GeV >20) & (lplus_vec.pt / GeV >20))
    low_pT_lepton_mask = ((lminus_vec.pt / GeV >5) & (lplus_vec.pt / GeV >5))
    tautau_invmass_cut = (ditau_vec.m  / GeV < 500)
    lepton_mask = (
        # both_lepton_tight_acc_mask
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
    for member_name in tree["NNPDF31LO_Members"].fields:
        nnpdf31lo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["NNPDF31LO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["MSHT20LO_Members"].fields:
        msht20lo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["MSHT20LO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["NNPDF31NLO_Members"].fields:
        nnpdf31nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["NNPDF31NLO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["CT18NLO_Members"].fields:
        ct18nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["CT18NLO_Members"][member_name][lepton_mask]*scale_factor
        )
    for member_name in tree["MSHT20NLO_Members"].fields:
        msht20nlo_hist_dict[member_name].fill(
            dilepton_vec.m / GeV, 
            weight=tree["MSHT20NLO_Members"][member_name][lepton_mask]*scale_factor
        )

    unweighted_event_counter += len(dilepton_vec)
    if args.debug:
        break

# Print Statements:
print(f"Unweighted Events: {unweighted_event_counter}")
print(f"Weighted Events: {unweighted_event_counter * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * unweighted_event_counter * scale_factor}")

if args.debug:
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
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Linear K-Factor Binning",
                luminosity=luminosity)
at.create_stair(dilepton_id_mass_pdfreweight_hist, "DiLepton Mass Log K-Factor Binning", yscale='log',
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
# NNPDF31LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(nnpdf31lo_hist_dict),
    "NNPD31LO"
)
# MSHT20LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(msht20lo_hist_dict), "MSHT20LO"
)
# NNPDF31NLO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(nnpdf31nlo_hist_dict),
    "NNPD31NLO"
)
# CT18NLO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(ct18nlo_hist_dict),
    "CT18NLO"
)
# MSHT20LO RMS Plots
plot_pdf_rms(
    at.calc_pdf_rms(msht20nlo_hist_dict),
    "MSHT20NLO"
)
# Many Family PDF Plots
fig, axs = plt.subplots()
plt.subplots_adjust(top=0.85)
# axs.stairs(
#     dilepton_id_mass_pdfreweight_hist.view().value, 
#     edges=dilepton_id_mass_pdfreweight_hist.axes[0].edges,
#     label="CT09MCS"
# )
for weight_name in tree["pdfReweight"].fields:
    if "nlo" in weight_name:
        axs.stairs(
            weighthist_dict[weight_name].view().value, 
            edges=weighthist_dict[weight_name].axes[0].edges,
            label=weight_name.upper()
        )
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}}$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations.png')
plt.close()


# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    # root_file["DileptonKFactorFine"] = dilepton_id_mass_pdfreweight_hist
    root_file["DileptonKFactorFine"] = msht20nlo_hist_dict["PDFMember0Weight"]

# Save histograms
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
pickle_dict = {
    "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
    "DileptonMassFine": [dilepton_id_mass_finebin_hist],
    "DileptonKFactorFine": [dilepton_id_mass_pdfreweight_hist],
    # "DileptonKFactorFine": [msht20nlo_hist_dict["PDFMember0Weight"]],
    "DileptonKFactorProfile": [dilepton_id_mass_pdfreweight_profilehist],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
