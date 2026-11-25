/**
 * @file pythia_lhe_helpers.h
 * @brief Shared interfaces for converting LHE events to ROOT tuples.
 *
 * Compile all showering programs from the repository root with:
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
#include <vector>

namespace pythia_lhe {

/** Ordered particle indices returned by a process-specific selector. */
struct ParticleSelection {
  std::vector<int> indices;
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
   * Returns the ordered ROOT particle-branch names produced by this selector.
   *
   * Returns:
   *   One unique branch name for every index returned by select().
   */
  virtual std::vector<std::string> particleBranchNames() const = 0;

  /**
   * Returns the command file used when the user does not pass --command-file.
   *
   * Returns:
   *   Command-file basename located beside the executable or in the source
   *   directory.
   */
  virtual std::string defaultCommandFile() const;

  /**
   * Reports whether legacy LO PDF-reweighting aliases should be written.
   *
   * Returns:
   *   True for LO selectors requiring the analysis-compatible weight schema.
   */
  virtual bool useLoWeightAliases() const;

  /**
   * Selects the charged leptons to store for one accepted Pythia event.
   *
   * Args:
   *   event: Complete Pythia event record.
   *   selected: Output indices ordered like particleBranchNames().
   *   failureReason: Human-readable explanation when selection fails.
   *
   * Returns:
   *   True exactly when all required output indices are valid.
   */
  virtual bool select(const Pythia8::Event& event,
                      ParticleSelection& selected,
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
 * Tests whether a particle has a valid mother with the same PDG ID.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   index: Particle index to inspect.
 *
 * Returns:
 *   True when the particle is a radiation or history copy of its mother.
 */
bool hasSameIdMother(const Pythia8::Event& event, int index);

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
