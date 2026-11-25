/**
 * @file pythia_gen_lhe_tautau.cc
 * @brief Shower events and select the leading opposite-sign tau pair.
 *
 * Compile from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe pythia_gen_lhe_tautau
 * @endcode
 * Run with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_tautau input.lhe output.root --seed 1
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include <string>
#include <vector>

namespace {

/**
 * Finds the highest-pT root tau of one charge.
 *
 * Same-ID descendants are excluded so radiation copies cannot compete with
 * their original tau. The pT comparison uses the earliest copy in each tau
 * history, before following that history to its decay.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   tauId: Requested tau PDG ID, either 15 or -15.
 *
 * Returns:
 *   Index of the highest-pT root tau, or -1 when that charge is absent.
 */
int highestPtRootTau(const Pythia8::Event& event, int tauId) {
  int bestIndex = -1;
  double bestPt = -1.0;
  for (int index = 1; index < event.size(); ++index) {
    if (event[index].id() != tauId
        || pythia_lhe::hasSameIdMother(event, index)) {
      continue;
    }
    if (event[index].pT() > bestPt) {
      bestIndex = index;
      bestPt = event[index].pT();
    }
  }
  return bestIndex;
}

/**
 * Finds the unique direct electron or muon from a decaying tau.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   firstTau: Earliest copy of the selected tau history.
 *   bornLepton: Output index of the direct charged-lepton daughter.
 *   bareLepton: Output index of the final same-ID lepton copy.
 *
 * Returns:
 *   True exactly when the tau has one charge-compatible e/mu daughter.
 */
bool tauDecayLepton(const Pythia8::Event& event,
                    int firstTau,
                    int& bornLepton,
                    int& bareLepton) {
  const int tau = pythia_lhe::lastSameIdCopy(event, firstTau);
  const int tauId = event[tau].id();
  const std::vector<int> daughters = event[tau].daughterList();
  std::vector<int> candidates;
  for (std::vector<int>::const_iterator daughter = daughters.begin();
       daughter != daughters.end(); ++daughter) {
    if (*daughter <= 0 || *daughter >= event.size()) continue;
    const int id = event[*daughter].id();
    if ((tauId == 15 && (id == 11 || id == 13))
        || (tauId == -15 && (id == -11 || id == -13))) {
      candidates.push_back(*daughter);
    }
  }
  if (candidates.size() != 1) return false;
  bornLepton = candidates.front();
  bareLepton = pythia_lhe::lastSameIdCopy(event, bornLepton);
  return true;
}

/** Selects charged leptons from the leading root tau of each charge. */
class LeadingTauPairSelector : public pythia_lhe::EventSelector {
 public:
  /** Returns the stable tau-derived dilepton particle-branch order. */
  std::vector<std::string> particleBranchNames() const override {
    return std::vector<std::string>{"BornLepton", "BornAntiLepton",
                                    "BareLepton", "BareAntiLepton"};
  }

  /**
   * Forces every tau to use decay modes containing an electron or muon.
   *
   * Args:
   *   pythia: Configured but not yet initialized Pythia instance.
   *
   * Returns:
   *   True when all three Pythia settings are accepted.
   */
  bool configure(Pythia8::Pythia& pythia) const override {
    return pythia.readString("15:mayDecay = on")
           && pythia.readString("15:onMode = off")
           && pythia.readString("15:onIfAny = 11 13");
  }

  /**
   * Selects charged daughters of the highest-pT opposite-sign root tau pair.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output born and bare lepton indices.
   *   failureReason: Explanation when the decay history is unsuitable.
   *
   * Returns:
   *   True when both tau charges and their forced leptonic decays are found.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::ParticleSelection& selected,
              std::string& failureReason) const override {
    const int tauMinus = highestPtRootTau(event, 15);
    const int tauPlus = highestPtRootTau(event, -15);
    if (tauMinus < 0 || tauPlus < 0) {
      failureReason = "could not find at least one root tau of each charge";
      return false;
    }
    int bornMinus = -1;
    int bareMinus = -1;
    int bornPlus = -1;
    int barePlus = -1;
    if (!tauDecayLepton(event, tauMinus, bornMinus, bareMinus)
        || !tauDecayLepton(event, tauPlus, bornPlus, barePlus)) {
      failureReason =
          "the highest-pT opposite-sign taus do not each have one direct e/mu daughter";
      return false;
    }
    selected.indices =
        std::vector<int>{bornMinus, bornPlus, bareMinus, barePlus};
    return true;
  }

  /** Returns the ROOT key used for the tau-decay interpretation note. */
  std::string metadataName() const override { return "TauDecayMetadata"; }

  /** Returns the tau forcing and normalization interpretation note. */
  std::string metadataText() const override {
    return "Both tau charges were forced to decay through modes containing an "
           "electron or muon (15:onMode = off; 15:onIfAny = 11 13).\n"
           "The highest-pT root tau- and highest-pT root tau+ are selected "
           "independently; they are not required to share a parent. "
           "BornLepton and BornAntiLepton are their direct charged-lepton "
           "daughters; Bare branches are the last same-PDG-ID copies after "
           "radiation.\nPythia renormalizes enabled "
           "tau branching fractions to unity; stored weights and the reported "
           "production cross section do not include the physical double-"
           "leptonic tau branching fraction.";
  }
};

}  // namespace

/**
 * Runs the forced leading-opposite-sign-tau LHE showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const LeadingTauPairSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
