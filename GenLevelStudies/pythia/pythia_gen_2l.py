# Imports
# STL Packages
import sys
import pdb
# Scikit Packages
import numpy as np
# HEP Packages
import pythia8
import ROOT
# Personal Packages
sys.path.append(".") # Not great form.
import AnalysisTools as at
import PDFReweighter

# Classes
# Define a custom user hook class
class InAccHook(pythia8.UserHooks):
    def __init__(self, target_pid, min_eta, min_pT):
        super().__init__()
        self.target_pid = target_pid
        self.min_eta = min_eta
        self.min_pT = min_pT
        self.lepton_pid_array = np.array([11, 13])
        self.neutrino_pid_array = np.array([12, 14])

    # Allow process cross section to be modified...
    def canModifySigma(self): return True

    def multiplySigmaBy(self, sigmaProcessPtr, phaseSpacePtr, inEvent):
        # All events should be 2 -> 2, but kill them if not.
        if sigmaProcessPtr.nFinal() != 2: return 0.

        # Extract the pT for 2 -> 2 processes in the event generation chain
        # (inEvent = false for initialization).
        if inEvent:
            self.pTHat = phaseSpacePtr.pTHat()

        # Here we do not modify 2 -> 2 cross sections.
        return 1.

    # Allow a veto for the interleaved evolution in pT.
    # def canVetoPT(self): return True
    
    # Do the veto test at a pT scale of 5 GeV.
    # def scaleVetoPT(self): return 5.

    # The checkVetoEvent() method is called for each event before it is accepted
    def canVetoProcessLevel(self):
        return True  # tell PYTHIA we will sometimes veto whole events

    def doVetoProcessLevel(self, process):
        for index, particle in enumerate(pythia.process):
            if particle.id()==self.target_pid:
                itarget_particle = index
            elif particle.id()==-1*self.target_pid:
                itarget_antiparticle = index
        targetpart_mom = pythia.process[itarget_particle].p()
        targetantipart_mom = pythia.process[itarget_antiparticle].p()
        targetpart_in_acc = bool((targetpart_mom.eta() > self.min_eta) & (targetpart_mom.pT() > self.min_pT))
        targetantipart_in_acc = bool((targetantipart_mom.eta() > self.min_eta) & (targetantipart_mom.pT() > self.min_pT))
        # veto event if mass is below threshold
        if (not targetpart_in_acc) or (not targetantipart_in_acc):
            return True  # veto
        return False  # accept

# Functions
def deltaR(vec1, vec2):
    delta_r = np.sqrt(
        (vec1.phi() - vec2.phi())**2
        + (vec1.eta() - vec2.eta())**2
    )
    return(delta_r)

def fill_UE_array(array, pythia, lepton_index):
    lepton_vec = pythia.event[lepton_index].p()
    spT_sum = 0
    num_good_cones = 0
    isGoodCone_array = [0]*7
    spT_array = [0]*7
    for i in range(7):
        cone_spT = 0
        is_good_cone = True
        lepton_vec.rot(0, (2*np.pi) / 6)
        for particle in pythia.event:
            if ((deltaR(particle.p(), lepton_vec) < 0.50) and (particle.pT() < 20)):
                cone_spT += particle.pT()
            elif (particle.pT() > 20):
                is_good_cone = False
                isGoodCone_array[i] = False
                spT_array[i] = 0
                break
        isGoodCone_array[i] = True
        spT_array[i] = cone_spT
        num_good_cones +=1
        spT_sum += cone_spT
    array[:] = (
        spT_sum,
        num_good_cones,
        spT_array[0],
        isGoodCone_array[0],
        spT_array[1],
        isGoodCone_array[1],
        spT_array[2],
        isGoodCone_array[2],
        spT_array[3],
        isGoodCone_array[3],
        spT_array[4],
        isGoodCone_array[4],
        spT_array[5],
        isGoodCone_array[5]
    )
    return(array)

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/data/home/ganowak/WWProduction/pythia8311/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "ttbarProduction.cmnd"
ofile_name = "Test.root"
target_pid = -24
lepton_pid_array = np.array([11, 13])
neutrino_pid_array = np.array([12, 14])
# jet_array = np.array([-5,-3,-1])
jet_array = np.array([])

# NNPDF31 Reweighter
pdfrwgt_NNPDF31 = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="NNPDF31_lo_as_0130"
)
# CT18NLO Reweighter
pdfrwgt_CT18LO = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="CT18LO"
)
# MSHT20LO Reweighter
pdfrwgt_MSHT20LO = PDFReweighter.PDFReweighter(
    pdf_set_from="CT09MCS",
    pdf_set_to="MSHT20lo_as130"
)

# Read in Card File 
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/defaults.cmnd")
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Attach early veto hook
# hook = InAccHook(target_pid, 1, 20)
# pythia.setUserHooksPtr(hook)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:charge/F:status/F'
ue_str = "spTAve/F:NumGoodCones/F:Cone1spT/F:Cone1isGoodCone/F:Cone2spT/F:Cone2isGoodCone/F:Cone3spT/F:Cone3isGoodCone/F:Cone4spT/F:Cone4isGoodCone/F:Cone5spT/F:Cone5isGoodCone/F:Cone6spT/F:Cone6isGoodCone/F:"
initial_conditions_str = 'x1/F:id1/F:x2/F:id2/F'
pdfrwgt_str = 'nnpdf31lo/F:ct18lo/F:msht20lo/F'
evt_array = np.array([0], dtype=np.float32)
target_particle_array = np.array([0]*12, dtype=np.float32)
target_lepton_array = np.array([0]*12, dtype=np.float32)
target_lepton_UE_array = np.array([0]*14, dtype=np.float32)
target_neutrino_array = np.array([0]*12, dtype=np.float32)
target_antiparticle_array = np.array([0]*12, dtype=np.float32)
target_antilepton_array = np.array([0]*12, dtype=np.float32)
target_antilepton_UE_array = np.array([0]*14, dtype=np.float32)
target_antineutrino_array = np.array([0]*12, dtype=np.float32)
initial_conditions_array = np.array([0]*4, dtype=np.float32)
pdfrwgt_array = np.array([0]*3, dtype=np.float32)
if jet_array.any():
    target_jet_array = np.array([0]*12, dtype=np.float32)
    target_antijet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                        "RECREATE")
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F')
tree.Branch(
    'InitialConditions', initial_conditions_array, initial_conditions_str
)
tree.Branch('pdfReweight', pdfrwgt_array, pdfrwgt_str)
tree.Branch('TargetParticle', target_particle_array, var_str)
tree.Branch('TargetLepton', target_lepton_array, var_str)
tree.Branch('TargetLeptonUE', target_lepton_UE_array, ue_str)
tree.Branch('TargetNeutrino', target_neutrino_array, var_str)
tree.Branch('TargetAntiParticle', target_antiparticle_array, var_str)
tree.Branch('TargetAntiLepton', target_antilepton_array, var_str)
tree.Branch('TargetAntiLeptonUE', target_antilepton_UE_array, ue_str)
tree.Branch('TargetAntiNeutrino', target_antineutrino_array, var_str)
if jet_array.any():
    tree.Branch('TargetJet', target_jet_array, var_str)
    tree.Branch('TargetAntiJet', target_antijet_array, var_str)

# Testing
counter1 = 0
counter2 = 0

# Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent
    for index, particle in enumerate(pythia.event):
        if particle.id()==target_pid:
            itarget_particle = index
        elif particle.id()==-1*target_pid:
            itarget_antiparticle = index
    for idaughter in at.get_target_children(itarget_particle, pythia.event, child_index_list=[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in lepton_pid_array:
            ilepton_plus = idaughter
        elif daughter.id() in -1*neutrino_pid_array:
            ineutrino_plus = idaughter
        elif daughter.id() in jet_array:
            ijet = idaughter
    for idaughter in at.get_target_children(itarget_antiparticle, pythia.event, child_index_list=[]):
        daughter = pythia.event[idaughter]
        if daughter.id() in -1*lepton_pid_array:
            ilepton_minus = idaughter
        elif daughter.id() in neutrino_pid_array:
            ineutrino_minus = idaughter
        elif daughter.id() in -1*jet_array:
            iantijet = idaughter

    # Get Initial Conditions Info
    x1pdf = pythia.infoPython().x1pdf()
    x2pdf = pythia.infoPython().x2pdf()
    id1pdf = pythia.infoPython().id1pdf()
    id2pdf = pythia.infoPython().id2pdf()
    nnpdf_reweight = pdfrwgt_NNPDF31.calculate_reweighting_factor(
        x1pdf, x2pdf, id1pdf, id2pdf, sqrt_s=13e3
    )
    ct18lo_reweight = pdfrwgt_CT18LO.calculate_reweighting_factor(
        x1pdf, x2pdf, id1pdf, id2pdf, sqrt_s=13e3
    )
    msht20lo_reweight = pdfrwgt_MSHT20LO.calculate_reweighting_factor(
        x1pdf, x2pdf, id1pdf, id2pdf, sqrt_s=13e3
    )
    pdfrwgt_array[:] = (nnpdf_reweight, ct18lo_reweight, msht20lo_reweight)
    initial_conditions_array[:] = (x1pdf, id1pdf, x2pdf, id2pdf)
    target_particle_array = at.fill_array(target_particle_array, pythia.event,
                                        itarget_particle)
    target_lepton_array = at.fill_array(target_lepton_array, pythia.event,
                                        ilepton_plus)
    target_lepton_UE_array = fill_UE_array(target_lepton_UE_array, pythia, ilepton_minus)
    target_neutrino_array = at.fill_array(target_neutrino_array, pythia.event,
                                        ineutrino_plus)
    target_antiparticle_array = at.fill_array(target_antiparticle_array, pythia.event,
                                        itarget_antiparticle)
    target_antilepton_array = at.fill_array(target_antilepton_array, pythia.event,
                                        ilepton_minus)
    target_antilepton_UE_array = fill_UE_array(target_antilepton_UE_array, pythia, ilepton_plus)
    target_antineutrino_array = at.fill_array(target_antineutrino_array, pythia.event,
                                        ineutrino_minus)
    if jet_array.any():
        target_jet_array = at.fill_array(target_jet_array, pythia.event,
                                        ijet)
        target_antijet_array = at.fill_array(target_antijet_array, pythia.event,
                                        iantijet)
    tree.Fill()

print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write()
