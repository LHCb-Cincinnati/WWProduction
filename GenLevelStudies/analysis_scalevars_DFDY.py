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
tree_iterator = uproot.iterate(
    [file_name + ":Tree"for file_name in file_name_list]
)

# Scale Factor
if args.cross_section == False:
    scale_factor = 1
else:
    scale_factor = args.cross_section

# Variables
weights_bool = True
unweighted_event_counter = 0

# Binning Scheme
bins_list = [0] + list(np.linspace(33, 100, 8)) + [200, 2000] 

# Define histograms
dilepton_id_mass_rghbin_hist = bh.Histogram(bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_finebin_hist = bh.Histogram(bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight())
dilepton_id_mass_kfactorbin_hist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.Weight()
)
dilepton_id_mass_kfactorbin_profilehist = bh.Histogram(
    bh.axis.Variable(bins_list), storage=bh.storage.WeightedMean()
)
eta_hist = bh.Histogram(
    bh.axis.Regular(12,-6, 6),
    bh.axis.Regular(12,-6, 6),
    storage=bh.storage.Weight()
)
if weights_bool:
    weight_name_list = ([
        "MUF=0.5_MUR=0.5_PDF=10770_MERGING=0.000",
        "MUF=0.5_MUR=1.0_PDF=10770_MERGING=0.000",
        "MUF=0.5_MUR=2.0_PDF=10770_MERGING=0.000",
        "MUF=1.0_MUR=0.5_PDF=10770_MERGING=0.000",
        "Weight",
        "MUF=1.0_MUR=2.0_PDF=10770_MERGING=0.000",
        "MUF=2.0_MUR=0.5_PDF=10770_MERGING=0.000",
        "MUF=2.0_MUR=1.0_PDF=10770_MERGING=0.000",
        "MUF=2.0_MUR=2.0_PDF=10770_MERGING=0.000",
    ])
    weight_name_list = ["AUX_MUR05_MUF05", "AUX_MUR05_MUF10", "AUX_MUR05_MUF20", "AUX_MUR10_MUF05", "AUX_MUR10_MUF10", "AUX_MUR10_MUF20", "AUX_MUR20_MUF05", "AUX_MUR20_MUF10", "AUX_MUR20_MUF20"]
    weighthist_scale_dict = {}
    lower_env_hist = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
    upper_env_hist = bh.Histogram(
        bh.axis.Variable(bins_list), storage=bh.storage.Weight()
    )
    for weight_name in weight_name_list:
        weighthist_scale_dict[weight_name] = bh.Histogram(
            bh.axis.Variable(bins_list), storage=bh.storage.Weight()
        )
    pdfweight_name_list = ["nnpdf31lo", "ct18lo", "msht20lo", "nnpdf31nlo", "ct18nlo", "msht20nlo"]
    weighthist_pdf_dict = {}
    for weight_name in pdfweight_name_list:
        weighthist_pdf_dict[weight_name + "__kfactor__bin"] = bh.Histogram(
            bh.axis.Variable(bins_list), storage=bh.storage.Weight()
        )
        weighthist_pdf_dict[weight_name + "__rgh__bin"] = bh.Histogram(
            bh.axis.Regular(26, 20, 306), storage=bh.storage.Weight()
        )
        weighthist_pdf_dict[weight_name + "__fine__bin"] = bh.Histogram(
            bh.axis.Regular(50, 20, 306), storage=bh.storage.Weight()
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
    both_lepton_tight_acc_mask = (
        (lminus_vec.eta>2.2)
        & (lminus_vec.eta<4.4)
        & (lplus_vec.eta>2.2)
        & (lplus_vec.eta<4.4)
    ) 
    mue_decay_mask = (((lminus_vec.pid==13) & (lplus_vec.pid==-11))
                        | ((lminus_vec.pid==11) & (lplus_vec.pid==-13)))
    invariant_mass_mask = (dilepton_vec.m / GeV >100)
    deltar_mask = (np.abs(lminus_vec.deltaR(lplus_vec)) > 0.1)
    high_pT_lepton_mask = ((lminus_vec.pt / GeV >20) & (lplus_vec.pt / GeV >20))
    tautau_invmass_cut = (ditau_vec.m  / GeV < 500)
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
    dilepton_id_mass_kfactorbin_hist.fill(dilepton_vec.m / GeV, weight=scale_factor)
    dilepton_id_mass_kfactorbin_profilehist.fill(
        dilepton_vec.m / GeV, weight=scale_factor, sample=dilepton_vec.m / GeV
    )
    eta_hist.fill(
        muon_vec.eta, electron_vec.eta, weight=scale_factor
    )
    if weights_bool:
        for weight_name in tree["Weights"].fields:
            weighthist_scale_dict[weight_name].fill(
                dilepton_vec.m / GeV, weight=tree["Weights"][weight_name][lepton_mask]
            )
        for weight_name in tree["pdfReweight"].fields:
            weighthist_pdf_dict[weight_name + "__kfactor__bin"].fill(
                dilepton_vec.m / GeV, 
                weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
            )
            weighthist_pdf_dict[weight_name + "__rgh__bin"].fill(
                dilepton_vec.m / GeV, 
                weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
            )
            weighthist_pdf_dict[weight_name + "__fine__bin"].fill(
                dilepton_vec.m / GeV, 
                weight=tree["pdfReweight"][weight_name][lepton_mask]*scale_factor
            )
    unweighted_event_counter += len(dilepton_vec)
    if args.debug:
        break

# Weight Hists
if weights_bool:
    # Scale Variation Weights 
    for index in range(len(dilepton_id_mass_kfactorbin_hist.view().value)):
        min_val = weighthist_scale_dict["AUX_MUR10_MUF10"].view()[index].value
        min_var = weighthist_scale_dict["AUX_MUR10_MUF10"].view()[index].variance 
        max_val = weighthist_scale_dict["AUX_MUR10_MUF10"].view()[index].value 
        max_var = weighthist_scale_dict["AUX_MUR10_MUF10"].view()[index].variance
        for weight_name in tree["Weights"].fields:
            hist_view = weighthist_scale_dict[weight_name].view()[index]
            if hist_view.value < min_val:
                min_val = hist_view.value
                min_var = hist_view.variance 
            if hist_view.value > max_val:
                max_val = hist_view.value
                max_var = hist_view.variance 
        lower_env_hist.view()[index] = [min_val, min_var]
        upper_env_hist.view()[index] = [max_val, max_var]
        
    # Mean PDF histogram calculations
    weighthist_pdf_dict = at.calc_pdf_mean(weighthist_pdf_dict, "__kfactor__bin")
    weighthist_pdf_dict = at.calc_pdf_mean(weighthist_pdf_dict, "__rgh__bin")
    weighthist_pdf_dict = at.calc_pdf_mean(weighthist_pdf_dict, "__fine__bin")
    pdf_scale_factor = (
        sum(weighthist_pdf_dict['nlo_mean__kfactor__bin'].view().value) 
        / sum(dilepton_id_mass_kfactorbin_hist.view().value) 
    )

# Print Statements:
print(f"Unweighted Events: {unweighted_event_counter}")
print(f"Weighted Events: {unweighted_event_counter * scale_factor}")
if weights_bool: 
    print(f"PDF-Scaled Weighted Events: {unweighted_event_counter * scale_factor * pdf_scale_factor}")
    print(f"Lower Bound: {sum(lower_env_hist.view().value) / sum(weighthist_scale_dict['AUX_MUR10_MUF10'].view().value) * unweighted_event_counter * scale_factor * pdf_scale_factor}")
    print(f"Upper Bound: {sum(upper_env_hist.view().value) / sum(weighthist_scale_dict['AUX_MUR10_MUF10'].view().value) * unweighted_event_counter * scale_factor * pdf_scale_factor}")

if weights_bool:
# Divide Hists
    lower_env_hist = at.divide_bh_histograms(lower_env_hist, weighthist_scale_dict["AUX_MUR10_MUF10"])
    upper_env_hist = at.divide_bh_histograms(upper_env_hist, weighthist_scale_dict["AUX_MUR10_MUF10"])
    central_hist = at.divide_bh_histograms(weighthist_scale_dict["AUX_MUR10_MUF10"], weighthist_scale_dict["AUX_MUR10_MUF10"])

# Multiply Hists
    lower_env_hist = at.multiply_bh_histograms(weighthist_pdf_dict["nlo_mean__kfactor__bin"], lower_env_hist)
    upper_env_hist = at.multiply_bh_histograms(upper_env_hist, weighthist_pdf_dict["nlo_mean__kfactor__bin"])
    central_hist = weighthist_pdf_dict["nlo_mean__kfactor__bin"]

if args.debug:
    exit()

# Save Figures
folder_path = at.create_folder_path(ofile_name, args.testing)
os.chdir(folder_path)
# Plot
if weights_bool:
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__rgh__bin"], 
        "DiLepton Mass Linear Rough Binning",
        luminosity=luminosity
    )
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__rgh__bin"], 
        "DiLepton Mass Log Rough Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__fine__bin"], 
        "DiLepton Mass Linear Fine Binning",
        luminosity=luminosity
    )
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__fine__bin"], 
        "DiLepton Mass Log Fine Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__kfactor__bin"], 
        "DiLepton Mass Linear K-Factor Binning",
        luminosity=luminosity
    )
    at.create_stair(
        weighthist_pdf_dict["nlo_mean__kfactor__bin"], 
        "DiLepton Mass Log K-Factor Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_kfactorbin_profilehist, 
        "DiLepton Mass Linear K-Factor Binning Profile Hist",
        luminosity=luminosity
    )
else:
    at.create_stair(
        dilepton_id_mass_rghbin_hist, 
        "DiLepton Mass Linear Rough Binning",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_rghbin_hist, 
        "DiLepton Mass Log Rough Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_finebin_hist, 
        "DiLepton Mass Linear Fine Binning",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_finebin_hist, 
        "DiLepton Mass Log Fine Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_kfactorbin_hist, 
        "DiLepton Mass Linear K-Factor Binning",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_kfactorbin_hist, 
        "DiLepton Mass Log K-Factor Binning",
        yscale="log",
        luminosity=luminosity
    )
    at.create_stair(
        dilepton_id_mass_kfactorbin_profilehist, 
        "DiLepton Mass Linear K-Factor Binning Profile Hist",
        luminosity=luminosity
    )
# Envelope Calculation
if weights_bool:
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
    axs.set_title("Scale Variations")
    axs.set_xlabel("$M_{e \\mu} (GeV)$")
    axs.set_ylabel("$ \\frac{d \\sigma}{d M_{e \\mu}} \\left( \\frac{\\mathrm{fb}}{\\mathrm{GeV}} \\right)$")
# hist_handles, hist_labels = axs.get_legend_handles_labels()
# axs.legend(hist_handles, hist_labels)
# Slightly fancy to remove whitespace
    fig.savefig('DiLeptonMassScaleVariations.png')
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
axs.set_ylabel("$\eta_e$")
axs.set_xlabel("$\eta_{\mu}$")
save_str = ''.join("EtaEtaYieldHist".split())
plt.savefig(folder_path + "/" + save_str + '.png')
plt.close()

# Save Histograms
if weights_bool:
# Mu=1.0 Histogram
    with uproot.recreate(ofile_name + ".root") as root_file:
        root_file["DileptonKFactorFine"] = weighthist_pdf_dict["nlo_mean__kfactor__bin"]
        # root_file["DileptonKFactorFine"] = upper_env_hist
        root_file["EtaEtaYield"] = eta_hist

    os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
    pickle_dict = {
        "DiLeptonMassRough": [weighthist_pdf_dict["nlo_mean__rgh__bin"]],
        "DileptonMassFine": [weighthist_pdf_dict["nlo_mean__fine__bin"]],
        "DileptonKFactorFine": [weighthist_pdf_dict["nlo_mean__kfactor__bin"]],
        # "DileptonKFactorFine": [upper_env_hist],
        "DileptonKFactorProfile": [dilepton_id_mass_kfactorbin_profilehist],
        "EtaEtaYield": [eta_hist]
    }
    with open(ofile_name + ".pkl", "wb") as f:
        pickle.dump(pickle_dict, f)
# Upper Envelope Histogram
    upper_env_ofile_name = "DFDY_MG5_NLO_rwgt_upperenv_CentralMeanPDF_LowMass"
    pickle_dict = {
        "DiLeptonMassRough": [weighthist_pdf_dict["nlo_mean__rgh__bin"]],
        "DileptonMassFine": [weighthist_pdf_dict["nlo_mean__fine__bin"]],
        "DileptonKFactorFine": [upper_env_hist],
        # "DileptonKFactorFine": [upper_env_hist],
        "DileptonKFactorProfile": [dilepton_id_mass_kfactorbin_profilehist],
        "EtaEtaYield": [eta_hist]
    }
    with open(upper_env_ofile_name + ".pkl", "wb") as f:
        pickle.dump(pickle_dict, f)
    with uproot.recreate(upper_env_ofile_name + ".root") as root_file:
        root_file["DileptonKFactorFine"] = upper_env_hist
        # root_file["DileptonKFactorFine"] = upper_env_hist
        root_file["EtaEtaYield"] = eta_hist

# Lower Envelope Histogram
    lower_env_ofile_name = "DFDY_MG5_NLO_rwgt_lowerenv_CentralMeanPDF_LowMass"
    pickle_dict = {
        "DiLeptonMassRough": [weighthist_pdf_dict["nlo_mean__rgh__bin"]],
        "DileptonMassFine": [weighthist_pdf_dict["nlo_mean__fine__bin"]],
        "DileptonKFactorFine": [lower_env_hist],
        # "DileptonKFactorFine": [upper_env_hist],
        "DileptonKFactorProfile": [dilepton_id_mass_kfactorbin_profilehist],
        "EtaEtaYield": [eta_hist]
    }
    with open(lower_env_ofile_name + ".pkl", "wb") as f:
        pickle.dump(pickle_dict, f)
    with uproot.recreate(lower_env_ofile_name + ".root") as root_file:
        root_file["DileptonKFactorFine"] = lower_env_hist
        # root_file["DileptonKFactorFine"] = upper_env_hist
        root_file["EtaEtaYield"] = eta_hist
# Mu=1.0 Histogram
else:
    with uproot.recreate(ofile_name + ".root") as root_file:
        root_file["DileptonKFactorFine"] = dilepton_id_mass_kfactorbin_hist
        # root_file["DileptonKFactorFine"] = upper_env_hist
        root_file["EtaEtaYield"] = eta_hist

    os.chdir(at.find_WW_path() + "/GenLevelStudies/Histograms")
    pickle_dict = {
        "DiLeptonMassRough": [dilepton_id_mass_rghbin_hist],
        "DileptonMassFine": [dilepton_id_mass_finebin_hist],
        "DileptonKFactorFine": [dilepton_id_mass_kfactorbin_hist],
        # "DileptonKFactorFine": [upper_env_hist],
        "DileptonKFactorProfile": [dilepton_id_mass_kfactorbin_profilehist],
        "EtaEtaYield": [eta_hist]
    }
    with open(ofile_name + ".pkl", "wb") as f:
        pickle.dump(pickle_dict, f)

