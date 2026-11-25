/**
 * @file pythia_gen_lhe_tw_lo.cc
 * @brief Shower unmatched LO tW LHE events and store born/bare e-mu pairs.
 *
 * Compile from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe pythia_gen_lhe_tw_lo
 * @endcode
 * Run with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_tw_lo input.lhe output.root --seed 1
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

/** Selects the leading opposite-sign, different-flavour hard tW pair. */
class TwLoSelector : public pythia_lhe::EventSelector {
 public:
  /** Returns the unmatched-LO command-file basename. */
  std::string defaultCommandFile() const override {
    return "MG5LOShower.cmd";
  }

  /** Enables analysis-compatible LO PDF-weight aliases. */
  bool useLoWeightAliases() const override { return true; }

  /** Returns born and bare dilepton branch names in output order. */
  std::vector<std::string> particleBranchNames() const override {
    return std::vector<std::string>{"BornLepton", "BornAntiLepton",
                                    "BareLepton", "BareAntiLepton"};
  }

  /**
   * Selects the highest-summed-pT status-23 opposite-sign e-mu pair.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output born and last-copy particle indices.
   *   failureReason: Explanation when no suitable pair exists.
   *
   * Returns:
   *   True when a suitable hard-process pair is found.
   */
  bool select(const Pythia8::Event& event,
              pythia_lhe::ParticleSelection& selected,
              std::string& failureReason) const override {
    std::vector<int> negative;
    std::vector<int> positive;
    for (int index = 1; index < event.size(); ++index) {
      if (event[index].statusAbs() != 23) continue;
      const int id = event[index].id();
      if (id == 11 || id == 13) negative.push_back(index);
      if (id == -11 || id == -13) positive.push_back(index);
    }

    int bestNegative = -1;
    int bestPositive = -1;
    double bestPtSum = -1.0;
    for (std::vector<int>::const_iterator minus = negative.begin();
         minus != negative.end(); ++minus) {
      for (std::vector<int>::const_iterator plus = positive.begin();
           plus != positive.end(); ++plus) {
        if (std::abs(event[*minus].id()) == std::abs(event[*plus].id())) {
          continue;
        }
        const double ptSum = event[*minus].pT() + event[*plus].pT();
        const bool betterTie =
            ptSum == bestPtSum
            && (bestNegative < 0 || *minus < bestNegative
                || (*minus == bestNegative && *plus < bestPositive));
        if (ptSum > bestPtSum || betterTie) {
          bestPtSum = ptSum;
          bestNegative = *minus;
          bestPositive = *plus;
        }
      }
    }

    if (bestNegative < 0 || bestPositive < 0) {
      std::ostringstream reason;
      reason << "expected an opposite-sign, different-flavour status-23 e/mu "
             << "pair; found " << negative.size() << " negative and "
             << positive.size() << " positive hard e/mu candidates";
      failureReason = reason.str();
      return false;
    }

    selected.indices =
        std::vector<int>{bestNegative, bestPositive,
                         pythia_lhe::lastSameIdCopy(event, bestNegative),
                         pythia_lhe::lastSameIdCopy(event, bestPositive)};
    return true;
  }
};

}  // namespace

/**
 * Runs the unmatched LO tW showering executable.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *
 * Returns:
 *   Process exit status.
 */
int main(int argc, char* argv[]) {
  const TwLoSelector selector;
  return pythia_lhe::run(argc, argv, selector);
}
