""" Contains all the DaVinci sequences used in the WW production analysis.
"""

# Imports
from Configurables import CombineParticles
from PhysConf.Selections import Selection
from PhysConf.Selections import SelectionSequence
from PhysConf.Selections import CombineSelection
from Configurables import HltParticleFlow, HltJetBuilder
from Configurables import McParticleFlow, McJetBuilder

# Particle Collections
from StandardParticles import StdAllLooseMuons as Muons
from StandardParticles import StdAllLooseElectrons as Electrons
from StandardParticles import (StdLooseKsDD, StdLooseKsLL, StdLooseKsLD,
                               StdLooseLambdaDD, StdLooseLambdaLL, 
                               StdLooseLambdaLD)

def GetDiLeptonSequence():
    """ Creates a muon electron combine particles object.

    Creates a muon electron combine particles object.  Taken from the 
    StdAllLooseElectrons and StdAllLooseMuons collections.  Very few cuts are
    imposed on any aspect of the system to accept as many pairs as possible.
    A fake H_10 particle is used in the decay descriptor as a work around
    so that DaVinci can understand the dilepton object.

    Args:

    Returns:
        dilepton_seq (DaVinci Sequence): Muon Electron DaVinci sequence.
    """

    # dilepton Object Cuts
    # Cuts on the individual leptons, dilepton object, and parent particle.
    # All the cuts here are fake to accept almost any lepton pair
    # within the LHCb acceptance.
    dilepton_decay_products = {
        'mu+': '(PT > 5000*MeV) & (ETA < 5) & (ETA > 2)',
        'e-': '(PT > 5000*MeV) & (ETA < 5) & (ETA > 2)'
    }
    dilepton_comb = '(AM > 100*MeV)'
    dilepton_mother = '(MM > 100*MeV)'


    # dilepton Object Definition
    # Decay: H_10 -> mu+ mu-
    # Where H_10 is a fake parent particle created only to get information
    # on the leptons in the event.
    dilepton_sel = CombineSelection(
        'dilepton_sel',
        [Muons, Electrons],
        DecayDescriptor='[H_10 -> e- mu+]cc',
        DaughtersCuts=dilepton_decay_products,
        CombinationCut=dilepton_comb,
        MotherCut=dilepton_mother
    )
    dilepton_seq = SelectionSequence('dilepton_Seq', TopSelection=dilepton_sel)
    return(dilepton_seq)

def GetMCParticleFlow():
    """ Creates a MC particle flow object.

    Creates a MC particle flow object.  Taken from the MC/Particles
    collection.  Particles included in the flow are pions, kaons,
    lambdas, Chis, J/psi, Sigmas, and f0s.  Neutrinos are excluded.

    Args:

    Returns:
        mc_particle_flow (McParticleFlow): Particle flow of all MC particles
            in the event matching the inputs and exluding those banned.
    """

    mc_particle_flow = McParticleFlow('mcpf')
    mc_particle_flow.Inputs = [
        ['PID', 'ban', '12,-12,14,-14,16,-16'], #ban neutrinos
        ['PID', 'particle',
            '111,-111,' #pi0
            '321,-321,' #K+ and K-
            '130,-130,' #K0_L
            '310,-310,' #K0_S
            '211,-211,' #pi+ and pi-
            '3222,-3222,' #Sigma+
            '3122,-3122,' #Lambda
            '3112,-3112,' #Sigma-
            '3312,-3312,' #Xi-
            '3322,-3322,' #Xi0
            '30221,-30221,' #no particle assigned to this ID in PDG (?)
            '9010221,-9010221,' #f0(980)
            '443,-443' #J/psi
            ],
        ['MCParticle', 'daughters', 'MC/Particles']
        ]
    mc_particle_flow.Output = 'Phys/PF/MCParticles'
    return(mc_particle_flow)

def GetMCJetBuilder(mc_particle_flow):
    """ Creates a MC jet builder object.

    Creates a MC jet builder object.  Input is taken from the input
    mc_particle_flow.

    Args: 
        mc_particle_flow (McParticleFlow): Particle flow of all MC particles.

    Returns:
        mc_jet_builder (McJetBuilder): Jet of MC particles that shows up in
        Gaudi Python in the form of a MCParticle.  The PID will be set to 89.
    """

    mc_jet_builder = McJetBuilder('mcjb')
    mc_jet_builder.JetPtMin = 10000
    mc_jet_builder.JetR = 0.5
    mc_jet_builder.ChrVrt = True
    mc_jet_builder.NeuVrt = True
    mc_jet_builder.Inputs = [mc_particle_flow.Output]
    mc_jet_builder.Output = 'Phys/JB/MCParticles'
    return(mc_jet_builder)

def GetHLTParticleFlow(*args):
    """ Creates a HLT particle flow object.

    Creates a HLT particle flow object.  Taken from the ProtoParticle
    Charged and Neutral collections, as well as any user defined sequences.
    Particles included in the flow are protons, electrons, muons, pions, and
    kaons.

    Args:
        args (List[DaVinci Sequence]): A list of user defined DaVinci
            Sequences.

    Returns:
        hlt_particle_flow (HLTParticleFlow): Particle flow of all
            reconstructed particles in the event.
    """

    seq_list = ([['Particle', 'particle', arg.outputLocation()] for arg in args])
    hlt_particle_flow = HltParticleFlow('pf')
    hlt_particle_flow.Inputs = [*seq_list,
        ['ProtoParticle',  'best',     'Rec/ProtoP/Charged'],
        ['ProtoParticle',  'gamma',    'Rec/ProtoP/Neutrals']]
    hlt_particle_flow.Output = 'Phys/PF/Particles'
    # Particle names to use when assigning the 'best' PID to ProtoParticles.
    hlt_particle_flow.ProBestNames = ['mu+', 'e+', 'p+', 'K+', 'pi+']
    # ExtraInfo keys to check when assigning the 'best' PID to ProtoParticles.
    hlt_particle_flow.ProBestKeys  = [701,   700,  704,  703,  702]
    # Minimum required ExtraInfo values when assigning the 'best' PID to
    # ProtoParticles.
    hlt_particle_flow.ProBestMins  = [0.5,   0.5,  0.5,  0.5,  0.5]
    # Flag to only match the best ECAL CaloCluster for a Track.
    hlt_particle_flow.EcalBest = True 
    # Flag to perform neutral recovery using expected calorimeter response.
    hlt_particle_flow.SprRecover = False
    # Fix momentum to long Tracks with curvature / error below this.
    hlt_particle_flow.TrkLnErrMax = 10
    # Fix momentum for upstream Tracks with curvature / error below this.
    hlt_particle_flow.TrkUpErrMax = 10
    # Fix momentum for downstream Tracks with curvature / error below this.
    hlt_particle_flow.TrkDnErrMax = 10
    return(hlt_particle_flow)

def GetHLTJetBuilder(hlt_particle_flow):
    """
    """

    hlt_jet_builder = HltJetBuilder('jb')
    # If supplied, perform jet energy correction using histograms from this
    # path.
    hlt_jet_builder.JetEcPath = ''
    hlt_jet_builder.Inputs = [hlt_particle_flow.Output]
    hlt_jet_builder.Output = 'Phys/JB/Particles'
    hlt_jet_builder.JetPtMin = 10000
    return(hlt_jet_builder)