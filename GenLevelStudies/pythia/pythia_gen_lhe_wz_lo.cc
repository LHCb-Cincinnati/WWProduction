/**
 * @file pythia_gen_lhe_wz_lo.cc
 * @brief Shower unmatched LO WZ LHE events and store direct boson daughters.
 *
 * Compile from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe pythia_gen_lhe_wz_lo
 * @endcode
 * Run with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_wz_lo input.lhe output.root --seed 1
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

/** Direct charged-lepton daughters associated with one boson history. */
struct ChargedDaughters {
  int minus = -1;
  int plus = -1;
  int minusCount = 0;
  int plusCount = 0;
};

/**
 * Finds root histories for one signed or unsigned boson PDG ID.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   id: Requested signed PDG ID.
 *   useAbsoluteId: Whether either particle charge is accepted.
 *
 * Returns:
 *   Indices without a same-ID mother.
 */
std::vector<int> rootHistories(const Pythia8::Event& event,
                               int id,
                               bool useAbsoluteId) {
  std::vector<int> roots;
  for (int index = 1; index < event.size(); ++index) {
    const bool idMatches = useAbsoluteId
                               ? std::abs(event[index].id()) == std::abs(id)
                               : event[index].id() == id;
    if (idMatches && !pythia_lhe::hasSameIdMother(event, index)) {
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
 *   parent: Index of the decaying parent particle.
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

/** Selects direct charged-lepton daughters of one W and one Z history. */
class WzLoSelector : public pythia_lhe::EventSelector {
 public:
  /** Returns the unmatched-LO command-file basename. */
  std::string defaultCommandFile() const override {
    return "MG5LOShower.cmd";
  }

  /** Enables analysis-compatible LO PDF-weight aliases. */
  bool useLoWeightAliases() const override { return true; }

  /** Returns the WZ analysis particle-branch names in output order. */
  std::vector<std::string> particleBranchNames() const override {
    return std::vector<std::string>{"Z_lminus", "Z_lplus", "W_lepton"};
  }

  /**
   * Selects unique direct e/mu daughters of one root W and one root Z.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output Z-minus, Z-plus, and W-lepton indices.
   *   failureReason: Explanation when histories or daughters are ambiguous.
   *
   * Returns:
   *   True when exactly one valid WZ decay history is found.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::ParticleSelection& selected,
              std::string& failureReason) const override {
    const std::vector<int> zRoots = rootHistories(event, 23, false);
    const std::vector<int> wRoots = rootHistories(event, 24, true);
    if (zRoots.size() != 1 || wRoots.size() != 1) {
      std::ostringstream reason;
      reason << "expected one root Z and one root W history; found "
             << zRoots.size() << " Z and " << wRoots.size() << " W histories";
      failureReason = reason.str();
      return false;
    }

    const int z = pythia_lhe::lastSameIdCopy(event, zRoots.front());
    const int w = pythia_lhe::lastSameIdCopy(event, wRoots.front());
    const ChargedDaughters zDaughters = directChargedDaughters(event, z);
    const ChargedDaughters wDaughters = directChargedDaughters(event, w);
    if (zDaughters.minusCount != 1 || zDaughters.plusCount != 1) {
      failureReason =
          "the decaying Z does not have one direct e/mu of each charge";
      return false;
    }
    const int wLeptonCount = wDaughters.minusCount + wDaughters.plusCount;
    if (wLeptonCount != 1) {
      failureReason =
          "the decaying W does not have exactly one direct e/mu daughter";
      return false;
    }

    const int wLepton = wDaughters.minusCount == 1
                            ? wDaughters.minus
                            : wDaughters.plus;
    selected.indices = std::vector<int>{zDaughters.minus, zDaughters.plus,
                                        wLepton};
    return true;
  }
};

}  // namespace

/**
 * Runs the unmatched LO WZ showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const WzLoSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
