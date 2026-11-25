/**
 * @file pythia_gen_lhe_nlo.cc
 * @brief Shower an NLO LHE file and store exactly one prompt dilepton pair.
 *
 * Compile from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe pythia_gen_lhe_nlo
 * @endcode
 * Run with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_nlo input.lhe output.root --seed 1
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include <sstream>
#include <string>
#include <vector>

namespace {

/** Selects events containing exactly one prompt e/mu of each charge. */
class PromptDileptonSelector : public pythia_lhe::EventSelector {
 public:
  /** Returns the stable NLO dilepton particle-branch order. */
  std::vector<std::string> particleBranchNames() const override {
    return std::vector<std::string>{"BornLepton", "BornAntiLepton",
                                    "BareLepton", "BareAntiLepton"};
  }

  /**
   * Selects the unique status-23 charged-lepton pair.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output born and bare lepton indices.
   *   failureReason: Explanation when the prompt pair is not unique.
   *
   * Returns:
   *   True when exactly one negative and one positive prompt lepton exist.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::ParticleSelection& selected,
              std::string& failureReason) const override {
    std::vector<int> negative;
    std::vector<int> positive;
    for (int index = 1; index < event.size(); ++index) {
      const Pythia8::Particle& particle = event[index];
      if (particle.statusAbs() != 23) continue;
      if (particle.id() == 11 || particle.id() == 13) {
        negative.push_back(index);
      } else if (particle.id() == -11 || particle.id() == -13) {
        positive.push_back(index);
      }
    }

    if (negative.size() != 1 || positive.size() != 1) {
      std::ostringstream reason;
      reason << "expected exactly one negative and one positive status-23 e/mu; found "
             << negative.size() << " negative and " << positive.size()
             << " positive candidates";
      failureReason = reason.str();
      return false;
    }

    selected.indices =
        std::vector<int>{negative.front(), positive.front(),
                         pythia_lhe::lastSameIdCopy(event, negative.front()),
                         pythia_lhe::lastSameIdCopy(event, positive.front())};
    return true;
  }
};

}  // namespace

/**
 * Runs the prompt-dilepton LHE showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const PromptDileptonSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
