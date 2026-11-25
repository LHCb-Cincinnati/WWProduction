/**
 * @file pythia_gen_lhe_zz_lo.cc
 * @brief Shower unmatched LO ZZ LHE events and store direct Z daughters.
 *
 * Compile from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe pythia_gen_lhe_zz_lo
 * @endcode
 * Run with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_zz_lo input.lhe output.root --seed 1
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <vector>

namespace {

/** Direct charged-lepton daughters associated with one Z history. */
struct ChargedDaughters {
  int minus = -1;
  int plus = -1;
  int minusCount = 0;
  int plusCount = 0;
};

/** One qualifying leptonic Z history and its direct daughters. */
struct LeptonicZ {
  int z = -1;
  int minus = -1;
  int plus = -1;
};

/**
 * Finds root Z histories in the Pythia event record.
 *
 * Args:
 *   event: Complete Pythia event record.
 *
 * Returns:
 *   Z indices without a same-ID mother.
 */
std::vector<int> rootZHistories(const Pythia8::Event& event) {
  std::vector<int> roots;
  for (int index = 1; index < event.size(); ++index) {
    if (event[index].id() == 23
        && !pythia_lhe::hasSameIdMother(event, index)) {
      roots.push_back(index);
    }
  }
  return roots;
}

/**
 * Collects direct electron and muon daughters by electric-charge sign.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   parent: Index of the decaying Z particle.
 *
 * Returns:
 *   Candidate indices and multiplicities for each charge.
 */
ChargedDaughters directChargedDaughters(const Pythia8::Event& event,
                                        int parent) {
  ChargedDaughters result;
  const std::vector<int> daughters = event[parent].daughterList();
  for (std::vector<int>::const_iterator daughter = daughters.begin();
       daughter != daughters.end(); ++daughter) {
    if (*daughter <= 0 || *daughter >= event.size()) continue;
    const int id = event[*daughter].id();
    if (id == 11 || id == 13) {
      result.minus = *daughter;
      ++result.minusCount;
    } else if (id == -11 || id == -13) {
      result.plus = *daughter;
      ++result.plusCount;
    }
  }
  return result;
}

/** Selects direct charged-lepton daughters of two leptonic Z histories. */
class ZzLoSelector : public pythia_lhe::EventSelector {
 public:
  /** Returns the unmatched-LO command-file basename. */
  std::string defaultCommandFile() const override {
    return "MG5LOShower.cmd";
  }

  /** Enables analysis-compatible LO PDF-weight aliases. */
  bool useLoWeightAliases() const override { return true; }

  /** Returns the ZZ analysis particle-branch names in output order. */
  std::vector<std::string> particleBranchNames() const override {
    return std::vector<std::string>{"Z1_lminus", "Z1_lplus",
                                    "Z2_lminus", "Z2_lplus"};
  }

  /**
   * Selects and pT-orders exactly two direct leptonic Z decays.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output Z1-minus, Z1-plus, Z2-minus, and Z2-plus indices.
   *   failureReason: Explanation when two unique histories are unavailable.
   *
   * Returns:
   *   True when exactly two qualifying Z histories are found.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::ParticleSelection& selected,
              std::string& failureReason) const override {
    const std::vector<int> roots = rootZHistories(event);
    std::vector<LeptonicZ> candidates;
    for (std::vector<int>::const_iterator root = roots.begin();
         root != roots.end(); ++root) {
      const int z = pythia_lhe::lastSameIdCopy(event, *root);
      const ChargedDaughters daughters = directChargedDaughters(event, z);
      if (daughters.minusCount == 1 && daughters.plusCount == 1) {
        LeptonicZ candidate;
        candidate.z = z;
        candidate.minus = daughters.minus;
        candidate.plus = daughters.plus;
        candidates.push_back(candidate);
      }
    }

    if (candidates.size() != 2) {
      std::ostringstream reason;
      reason << "expected exactly two root Z histories with unique direct e/mu "
             << "pairs; found " << candidates.size()
             << " qualifying histories";
      failureReason = reason.str();
      return false;
    }

    std::sort(candidates.begin(), candidates.end(),
              [&event](const LeptonicZ& first, const LeptonicZ& second) {
                if (event[first.z].pT() != event[second.z].pT()) {
                  return event[first.z].pT() > event[second.z].pT();
                }
                return first.z < second.z;
              });
    selected.indices =
        std::vector<int>{candidates[0].minus, candidates[0].plus,
                         candidates[1].minus, candidates[1].plus};
    return true;
  }
};

}  // namespace

/**
 * Runs the unmatched LO ZZ showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const ZzLoSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
