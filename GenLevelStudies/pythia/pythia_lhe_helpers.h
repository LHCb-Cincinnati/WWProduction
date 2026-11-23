/**
 * @file pythia_lhe_helpers.h
 * @brief Shared interfaces for converting LHE events to ROOT tuples.
 *
 * Compile both showering programs from the repository root with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe all
 * @endcode
 * Run a program with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_nlo input.lhe output.root
 * @endcode
 */

#ifndef GENLEVELSTUDIES_PYTHIA_PYTHIA_LHE_HELPERS_H_
#define GENLEVELSTUDIES_PYTHIA_PYTHIA_LHE_HELPERS_H_

#include "Pythia8/Pythia.h"

#include <string>

namespace pythia_lhe {

/** Indices of the particles stored in the four lepton ROOT branches. */
struct SelectedLeptons {
  int bornMinus = -1;
  int bornPlus = -1;
  int bareMinus = -1;
  int barePlus = -1;
};

/** Strategy interface for process-specific Pythia setup and event selection. */
class EventSelector {
 public:
  /** Provides a polymorphic destructor for selector implementations. */
  virtual ~EventSelector();

  /**
   * Applies process-specific settings before Pythia initialization.
   *
   * Args:
   *   pythia: Pythia instance that has read the common command file.
   *
   * Returns:
   *   True when all requested settings were accepted.
   */
  virtual bool configure(Pythia8::Pythia& pythia) const;

  /**
   * Selects the charged leptons to store for one accepted Pythia event.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output indices for born and bare charged leptons.
   *   failureReason: Human-readable explanation when selection fails.
   *
   * Returns:
   *   True exactly when all four output indices are valid.
   */
  virtual bool select(const Pythia8::Event& event,
                      SelectedLeptons& selected,
                      std::string& failureReason) const = 0;

  /** Returns the ROOT key for optional selector-specific metadata. */
  virtual std::string metadataName() const;

  /** Returns optional selector-specific metadata text. */
  virtual std::string metadataText() const;
};

/**
 * Follows a particle through consecutive daughters with the same PDG ID.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   index: Index of the earliest particle copy.
 *
 * Returns:
 *   Index of the last same-ID copy, or the input index when no copy exists.
 */
int lastSameIdCopy(const Pythia8::Event& event, int index);

/**
 * Runs the common LHE showering and ROOT-output workflow.
 *
 * Args:
 *   argc: Command-line argument count.
 *   argv: Command-line arguments.
 *   selector: Process-specific event selector.
 *
 * Returns:
 *   Zero on success and one after a reported error.
 */
int run(int argc, char* argv[], const EventSelector& selector);

}  // namespace pythia_lhe

#endif  // GENLEVELSTUDIES_PYTHIA_PYTHIA_LHE_HELPERS_H_
