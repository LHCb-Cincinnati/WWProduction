/**
 * @file pythia_gen_lhe_tautau.cc
 * @brief Shower Z-to-tau-pair events with forced leptonic tau decays.
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

#include <sstream>
#include <string>
#include <vector>

namespace {

/**
 * Reports whether a particle has a mother with the same PDG ID.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   index: Particle index to inspect.
 *
 * Returns:
 *   True when a valid same-ID mother exists.
 */
bool hasSameIdMother(const Pythia8::Event& event, int index) {
  if (index <= 0 || index >= event.size()) return false;
  const std::vector<int> mothers = event[index].motherList();
  for (std::vector<int>::const_iterator mother = mothers.begin();
       mother != mothers.end(); ++mother) {
    if (*mother > 0 && *mother < event.size()
        && event[*mother].id() == event[index].id()) {
      return true;
    }
  }
  return false;
}

/**
 * Finds a unique direct daughter with a requested PDG ID.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   parent: Parent particle index.
 *   daughterId: Requested daughter PDG ID.
 *   daughter: Output daughter index.
 *
 * Returns:
 *   True exactly when one matching direct daughter exists.
 */
bool uniqueDirectDaughter(const Pythia8::Event& event,
                          int parent,
                          int daughterId,
                          int& daughter) {
  daughter = -1;
  int matches = 0;
  const std::vector<int> daughters = event[parent].daughterList();
  for (std::vector<int>::const_iterator candidate = daughters.begin();
       candidate != daughters.end(); ++candidate) {
    if (*candidate > 0 && *candidate < event.size()
        && event[*candidate].id() == daughterId) {
      daughter = *candidate;
      ++matches;
    }
  }
  return matches == 1;
}

/**
 * Finds the unique direct electron or muon from a decaying tau.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   firstTau: Earliest tau copy descended directly from the Z.
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

/** Selects charged leptons from an explicit Z-to-tau-pair decay history. */
class ZTauTauSelector : public pythia_lhe::EventSelector {
 public:
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
   * Selects the charged daughters of the unique Z-sourced tau pair.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output born and bare lepton indices.
   *   failureReason: Explanation when the decay history is unsuitable.
   *
   * Returns:
   *   True when one unambiguous Z-to-tau-pair history is found.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::SelectedLeptons& selected,
              std::string& failureReason) const override {
    struct TauPair { int minus; int plus; };
    std::vector<TauPair> pairs;
    for (int index = 1; index < event.size(); ++index) {
      if (event[index].id() != 23 || hasSameIdMother(event, index)) continue;
      const int zDecay = pythia_lhe::lastSameIdCopy(event, index);
      TauPair pair = {-1, -1};
      if (uniqueDirectDaughter(event, zDecay, 15, pair.minus)
          && uniqueDirectDaughter(event, zDecay, -15, pair.plus)) {
        pairs.push_back(pair);
      }
    }

    if (pairs.size() != 1) {
      std::ostringstream reason;
      reason << "expected one explicit Z -> tau- tau+ history; found "
             << pairs.size();
      failureReason = reason.str();
      return false;
    }
    if (!tauDecayLepton(event, pairs[0].minus, selected.bornMinus,
                        selected.bareMinus)
        || !tauDecayLepton(event, pairs[0].plus, selected.bornPlus,
                           selected.barePlus)) {
      failureReason = "the Z-sourced taus do not each have one direct e/mu daughter";
      return false;
    }
    return true;
  }

  /** Returns the ROOT key used for the tau-decay interpretation note. */
  std::string metadataName() const override { return "TauDecayMetadata"; }

  /** Returns the tau forcing and normalization interpretation note. */
  std::string metadataText() const override {
    return "Both tau charges were forced to decay through modes containing an "
           "electron or muon (15:onMode = off; 15:onIfAny = 11 13).\n"
           "BornLepton and BornAntiLepton are direct daughters of the tau- and "
           "tau+ from an explicit Z decay; Bare branches are their last "
           "same-PDG-ID copies after radiation.\nPythia renormalizes enabled "
           "tau branching fractions to unity; stored weights and the reported "
           "production cross section do not include the physical double-"
           "leptonic tau branching fraction.";
  }
};

}  // namespace

/**
 * Runs the forced Z-to-tau-pair LHE showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const ZTauTauSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
