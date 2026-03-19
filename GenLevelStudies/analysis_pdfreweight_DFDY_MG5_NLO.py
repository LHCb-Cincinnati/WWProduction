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
bins_list = [0] + list(np.linspace(33, 100, 8)) + [200, 2000] 

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
    weighthist_dict[weight_name + "__kfactor__bin"] = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
    weighthist_dict[weight_name + "__rgh__bin"] = bh.Histogram(
        bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight()
    )
    weighthist_dict[weight_name + "__fine__bin"] = bh.Histogram(
        bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight()
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
    tautau_invmass_cut = (ditau_vec.m  / GeV > 500)
    lepton_mask = (
        both_lepton_tight_acc_mask
        & high_pT_lepton_mask
        & mue_decay_mask
        & deltar_mask
        & tautau_invmass_cut
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
        weighthist_dict[weight_name + "__kfactor__bin"].fill(
            dilepton_vec.m / GeV, 
            weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
        )
        weighthist_dict[weight_name + "__rgh__bin"].fill(
            dilepton_vec.m / GeV, 
            weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
        )
        weighthist_dict[weight_name + "__fine__bin"].fill(
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

# Mean PDF histogram calculations for rough and fine binned histogram
weighthist_pdf_dict = at.calc_pdf_mean(weighthist_dict, "__rgh__bin")
weighthist_pdf_dict = at.calc_pdf_mean(weighthist_dict, "__fine__bin")
# Special treatment for k-factor binning histogram w/ pdf uncertainties
# Create mean NLO PDF histogram from individual families
nnpdf31nlo_view = weighthist_dict["nnpdf31nlo__kfactor__bin"].view()
ct18nlo_view = weighthist_dict["ct18nlo__kfactor__bin"].view()
msht20nlo_view = weighthist_dict["msht20nlo__kfactor__bin"].view()
nlo_mean_hist = weighthist_dict["nnpdf31nlo__kfactor__bin"].copy()
nlo_mean_view = nlo_mean_hist.view()

# Per-bin mean of the three NLO histograms
values_stack = np.stack(
    [
        nnpdf31nlo_view.value,
        ct18nlo_view.value,
        msht20nlo_view.value,
    ],
    axis=0,
)
mean_values = np.mean(values_stack, axis=0)
nlo_mean_view.value = mean_values
sq_diff = (values_stack - mean_values) ** 2
rms_squared = np.mean(sq_diff, axis=0)
nlo_mean_view.variance = rms_squared
weighthist_dict["nlo_mean__kfactor__bin"] = nlo_mean_hist
# Create RMS NLO PDF histogram
nlo_rmslow_hist = weighthist_dict["nlo_mean__kfactor__bin"].copy()
nlo_rmslow_hist.view().value = nlo_rmslow_hist.view().value - nlo_rmslow_hist.view().variance
nlo_rmshigh_hist = weighthist_dict["nlo_mean__kfactor__bin"].copy()
nlo_rmshigh_hist.view().value = nlo_rmshigh_hist.view().value + nlo_rmshigh_hist.view().variance

# Print Statements:
print(f"Unweighted Events: {unweighted_event_counter}")
print(f"Weighted Events: {unweighted_event_counter * scale_factor}")
for weight_name in weight_name_list:
    print(f"{weight_name}: {weighthist_dict[weight_name + '__kfactor__bin'].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * unweighted_event_counter * scale_factor}")
print(f"Mean PDF: {weighthist_dict['nlo_mean__kfactor__bin'].view().value.sum() / dilepton_id_mass_pdfreweight_hist.view().value.sum() * unweighted_event_counter * scale_factor}")


if args.debug:
    exit()

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
at.create_stair(
    weighthist_dict["nlo_mean__rgh__bin"], 
    "DiLepton Mass Linear Rough Binning",
    luminosity=luminosity
)
at.create_stair(
    weighthist_dict["nlo_mean__rgh__bin"], 
    "DiLepton Mass Log Rough Binning",
    yscale="log",
    luminosity=luminosity
)
at.create_stair(
    weighthist_dict["nlo_mean__fine__bin"], 
    "DiLepton Mass Linear Fine Binning",
    luminosity=luminosity
)
at.create_stair(
    weighthist_dict["nlo_mean__fine__bin"], 
    "DiLepton Mass Log Fine Binning",
    yscale="log",
    luminosity=luminosity
)
at.create_stair(
    weighthist_dict["nlo_mean__kfactor__bin"], 
    "DiLepton Mass Linear K-Factor Binning",
    luminosity=luminosity
)
at.create_stair(
    weighthist_dict["nlo_mean__kfactor__bin"], 
    "DiLepton Mass Log K-Factor Binning",
    yscale="log",
    luminosity=luminosity
)

for weight_name in tree["pdfReweight"].fields:
    at.create_stacked_stair(
        [weighthist_dict[weight_name + "__kfactor__bin"], dilepton_id_mass_pdfreweight_hist],
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
# Plot NNPDF31NLO
axs.stairs(
    weighthist_dict["nnpdf31nlo__kfactor__bin"].view().value, 
    edges=weighthist_dict["nnpdf31nlo__kfactor__bin"].axes[0].edges,
    label="NNPDF31NLO",
    color="green"
)
# Plot CT18NLO
axs.stairs(
    weighthist_dict["ct18nlo__kfactor__bin"].view().value, 
    edges=weighthist_dict["ct18nlo__kfactor__bin"].axes[0].edges,
    label="CT18NLO",
    color="blue"
)
# Plot MSHT20NLO
axs.stairs(
    weighthist_dict["msht20nlo__kfactor__bin"].view().value, 
    edges=weighthist_dict["msht20nlo__kfactor__bin"].axes[0].edges,
    label="MSHT20NLO",
    color="red"
)
# Plot NLO Mean
axs.stairs(
    weighthist_dict["nlo_mean__kfactor__bin"].view().value, 
    edges=weighthist_dict["nlo_mean__kfactor__bin"].axes[0].edges,
    color="black"
)

axs.bar(
    x=weighthist_dict["nlo_mean__kfactor__bin"].axes[0].centers, 
    height= 2*np.sqrt(weighthist_dict["nlo_mean__kfactor__bin"].view().variance), 
    bottom=weighthist_dict["nlo_mean__kfactor__bin"].view().value - np.sqrt(weighthist_dict["nlo_mean__kfactor__bin"].view().variance), 
    width=weighthist_dict["nlo_mean__kfactor__bin"].axes[0].widths, 
    linewidth=0, 
    color="black", 
    alpha=0.25, 
    zorder=-1, 
    label="Mean w/ RMS Envelope"
)
axs.set_xlim((0, 300))
axs.set_title("PDF Variations")
axs.set_xlabel("$M_{e \\mu} (GeV)$")
axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
hist_handles, hist_labels = axs.get_legend_handles_labels()
axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
fig.savefig('DiLeptonMassPDFVariations.png')
plt.close()


# Save Histos in ROOT
os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
with uproot.recreate(ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = weighthist_dict["nlo_mean__kfactor__bin"]
    # root_file["DileptonKFactorFine"] = msht20nlo_hist_dict["PDFMember0Weight"]
# Save histograms
pickle_dict = {
    "DiLeptonMassRough": [weighthist_dict["nlo_mean__rgh__bin"]],
    "DileptonMassFine": [weighthist_dict["nlo_mean__fine__bin"]],
    "DileptonKFactorFine": [weighthist_dict["nlo_mean__kfactor__bin"]],
    # "DileptonKFactorFine": [msht20nlo_hist_dict["PDFMember0Weight"]],
}
with open(ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)

# Lower RMS Envelope Files
lower_rms_ofile_name = "DFDY_MG5_NLO_rwgt_mu10_LowerRMSMeanPDF_HighMass"
with uproot.recreate(lower_rms_ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = nlo_rmslow_hist
    # root_file["DileptonKFactorFine"] = msht20nlo_hist_dict["PDFMember0Weight"]
# Save histograms
pickle_dict = {
    "DiLeptonMassRough": [weighthist_dict["nlo_mean__rgh__bin"]],
    "DileptonMassFine": [weighthist_dict["nlo_mean__fine__bin"]],
    "DileptonKFactorFine": [nlo_rmslow_hist],
    # "DileptonKFactorFine": [msht20nlo_hist_dict["PDFMember0Weight"]],
}
with open(lower_rms_ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)

# Upper RMS Envelope Files
upper_rms_ofile_name = "DFDY_MG5_NLO_rwgt_mu10_UpperRMSMeanPDF_HighMass"
with uproot.recreate(upper_rms_ofile_name + ".root") as root_file:
    root_file["DileptonKFactorFine"] = nlo_rmshigh_hist
    # root_file["DileptonKFactorFine"] = msht20nlo_hist_dict["PDFMember0Weight"]
# Save histograms
pickle_dict = {
    "DiLeptonMassRough": [weighthist_dict["nlo_mean__rgh__bin"]],
    "DileptonMassFine": [weighthist_dict["nlo_mean__fine__bin"]],
    "DileptonKFactorFine": [nlo_rmshigh_hist],
    # "DileptonKFactorFine": [msht20nlo_hist_dict["PDFMember0Weight"]],
}
with open(upper_rms_ofile_name + ".pkl", "wb") as f:
    pickle.dump(pickle_dict, f)
