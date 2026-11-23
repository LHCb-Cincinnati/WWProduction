/**
 * @file pythia_lhe_helpers.cc
 * @brief Shared Pythia, JetMatching, LHE-weight, and ROOT-output implementation.
 *
 * Compile and link this helper through both executables with:
 * @code{.sh}
 * conda activate hep
 * make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe all
 * @endcode
 * Run a linked executable with:
 * @code{.sh}
 * GenLevelStudies/pythia/pythia_gen_lhe_nlo input.lhe output.root
 * @endcode
 */

#include "pythia_lhe_helpers.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/CombineMatchingInput.h"

#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pythia_lhe {

/** Parsed command-line options shared by all showering executables. */
struct Options {
  std::string lheFile;
  std::string outputFile;
  std::string commandFile = "GenLevelStudies/pythia/MG5NLOShower.cmd";
  int numberEvents = -1;
  int seed = -1;
  bool seedProvided = false;
  bool commandFileProvided = false;
  bool verbose = false;
};

/** Metadata for one named LHE variation weight. */
struct WeightInfo {
  std::string weightId;
  std::string leafName;
  std::string description;
};

/** Metadata for one LHE weight group and its ROOT branch. */
struct WeightGroup {
  std::string branchName;
  std::string originalName;
  std::vector<WeightInfo> weights;
};

/** Runtime storage backing one variable-size ROOT weight branch. */
struct WeightBranch {
  WeightGroup group;
  std::vector<float> values;
  std::string leafList;
};

/** Flat particle representation matching the historical ROOT leaf layout. */
struct ParticleRecord {
  float px = 0.0;
  float py = 0.0;
  float pz = 0.0;
  float pT = 0.0;
  float p = 0.0;
  float eta = 0.0;
  float e = 0.0;
  float phi = 0.0;
  float m0 = 0.0;
  float pid = 0.0;
  float charge = 0.0;
  float status = 0.0;
};

/** Initial parton flavours, momentum fractions, and factorization scale. */
struct InitialConditionsRecord {
  float x1 = 0.0;
  float id1 = 0.0;
  float x2 = 0.0;
  float id2 = 0.0;
  float Q2 = 0.0;
};

static_assert(sizeof(InitialConditionsRecord) == 5 * sizeof(float),
              "InitialConditionsRecord must match its ROOT leaf list");

/**
 * Prints command-line help.
 *
 * Args:
 *   program: Executable name displayed in the usage line.
 */
void printUsage(const char* program) {
  std::cerr
      << "Usage: " << program
      << " <lhe_file> <output_file> [-c command_file] [-n number_events]"
      << " [--seed seed] [--verbose]\n";
}

/**
 * Tests whether a path can be opened for reading.
 *
 * Args:
 *   path: Candidate input path.
 *
 * Returns:
 *   True when the file is readable.
 */
bool fileExists(const std::string& path) {
  std::ifstream input(path.c_str());
  return input.good();
}

/**
 * Extracts the directory component of a path.
 *
 * Args:
 *   path: Input filesystem path.
 *
 * Returns:
 *   Directory component or a suitable root/current-directory fallback.
 */
std::string dirname(const std::string& path) {
  const std::string::size_type slash = path.find_last_of('/');
  if (slash == std::string::npos) return ".";
  if (slash == 0) return "/";
  return path.substr(0, slash);
}

/**
 * Resolves the default production Pythia command file.
 *
 * Args:
 *   programPath: Path used to invoke the executable.
 *
 * Returns:
 *   Existing repository-relative or executable-relative command-file path.
 */
std::string resolveDefaultCommandFile(const std::string& programPath) {
  const std::string repoRootDefault = "GenLevelStudies/pythia/MG5NLOShower.cmd";
  if (fileExists(repoRootDefault)) return repoRootDefault;

  const std::string executableDirDefault = dirname(programPath) + "/MG5NLOShower.cmd";
  if (fileExists(executableDirDefault)) return executableDirDefault;

  return repoRootDefault;
}

/**
 * Parses a complete decimal integer option.
 *
 * Args:
 *   value: Text to parse.
 *   optionName: Option name used in error messages.
 *
 * Returns:
 *   Parsed integer.
 *
 * Throws:
 *   std::runtime_error: The value is not a complete integer.
 */
int parseInteger(const std::string& value, const std::string& optionName) {
  char* end = 0;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    throw std::runtime_error(optionName + " requires an integer value");
  }
  return static_cast<int>(parsed);
}

/**
 * Parses the common showering command line.
 *
 * Args:
 *   argc: Argument count.
 *   argv: Argument values.
 *
 * Returns:
 *   Validated shared options.
 */
Options parseArgs(int argc, char* argv[]) {
  Options options;
  std::vector<std::string> positional;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-c" || arg == "--command-file") {
      if (++i >= argc) throw std::runtime_error(arg + " requires a value");
      options.commandFile = argv[i];
      options.commandFileProvided = true;
    } else if (arg == "-n" || arg == "--number-events") {
      if (++i >= argc) throw std::runtime_error(arg + " requires a value");
      options.numberEvents = parseInteger(argv[i], arg);
    } else if (arg == "--seed") {
      if (++i >= argc) throw std::runtime_error(arg + " requires a value");
      options.seed = parseInteger(argv[i], arg);
      options.seedProvided = true;
    } else if (arg == "--verbose") {
      options.verbose = true;
    } else if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      std::exit(0);
    } else if (!arg.empty() && arg[0] == '-') {
      throw std::runtime_error("Unknown option: " + arg);
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.size() != 2) {
    printUsage(argv[0]);
    throw std::runtime_error("Expected exactly two positional arguments");
  }

  options.lheFile = positional[0];
  options.outputFile = positional[1];
  if (!options.commandFileProvided) {
    options.commandFile = resolveDefaultCommandFile(argv[0]);
  }
  return options;
}

/**
 * Removes leading and trailing whitespace.
 *
 * Args:
 *   input: Original text.
 *
 * Returns:
 *   Trimmed text.
 */
std::string trim(const std::string& input) {
  std::string::size_type first = 0;
  while (first < input.size()
         && std::isspace(static_cast<unsigned char>(input[first]))) {
    ++first;
  }

  std::string::size_type last = input.size();
  while (last > first
         && std::isspace(static_cast<unsigned char>(input[last - 1]))) {
    --last;
  }

  return input.substr(first, last - first);
}

/**
 * Replaces each whitespace run with a single space.
 *
 * Args:
 *   input: Original text.
 *
 * Returns:
 *   Normalized text.
 */
std::string collapseWhitespace(const std::string& input) {
  std::ostringstream output;
  bool inWhitespace = false;
  for (std::string::const_iterator it = input.begin(); it != input.end(); ++it) {
    if (std::isspace(static_cast<unsigned char>(*it))) {
      inWhitespace = true;
    } else {
      if (inWhitespace && output.tellp() > 0) output << ' ';
      output << *it;
      inWhitespace = false;
    }
  }
  return output.str();
}

/**
 * Converts ASCII text to lowercase.
 *
 * Args:
 *   input: Text copied for conversion.
 *
 * Returns:
 *   Lowercase text.
 */
std::string lowercase(std::string input) {
  std::transform(input.begin(), input.end(), input.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return input;
}

/**
 * Converts arbitrary text into a safe ROOT branch or leaf name.
 *
 * Args:
 *   name: Original name.
 *   fallback: Name used when no valid characters remain.
 *
 * Returns:
 *   Sanitized ROOT identifier.
 */
std::string sanitizeRootName(const std::string& name,
                             const std::string& fallback) {
  std::string clean;
  bool previousUnderscore = false;

  for (std::string::const_iterator it = name.begin(); it != name.end(); ++it) {
    const unsigned char c = static_cast<unsigned char>(*it);
    const bool keep = std::isalnum(c) || *it == '_';
    const char out = keep ? static_cast<char>(*it) : '_';

    if (out == '_') {
      if (!previousUnderscore && !clean.empty()) clean.push_back(out);
      previousUnderscore = true;
    } else {
      clean.push_back(out);
      previousUnderscore = false;
    }
  }

  while (!clean.empty() && clean[clean.size() - 1] == '_') clean.erase(clean.size() - 1);
  if (clean.empty()) clean = fallback;
  if (std::isdigit(static_cast<unsigned char>(clean[0]))) clean = "_" + clean;
  return clean;
}

/**
 * Makes a name unique within a mutable set.
 *
 * Args:
 *   name: Preferred name.
 *   used: Names already reserved; updated by this function.
 *
 * Returns:
 *   Reserved unique name.
 */
std::string uniqueName(const std::string& name, std::set<std::string>& used) {
  if (used.insert(name).second) return name;

  int index = 2;
  while (used.find(name + "_" + std::to_string(index)) != used.end()) ++index;

  const std::string unique = name + "_" + std::to_string(index);
  used.insert(unique);
  return unique;
}

/**
 * Reads a quoted XML-style attribute from tag text.
 *
 * Args:
 *   attributes: Raw opening-tag attributes.
 *   name: Attribute name.
 *   fallback: Value returned when the attribute is absent.
 *
 * Returns:
 *   Parsed value or fallback.
 */
std::string getAttribute(const std::string& attributes,
                         const std::string& name,
                         const std::string& fallback) {
  std::string::size_type pos = attributes.find(name);
  while (pos != std::string::npos) {
    const bool startsAtBoundary =
        pos == 0
        || !(std::isalnum(static_cast<unsigned char>(attributes[pos - 1]))
             || attributes[pos - 1] == '_');
    const std::string::size_type afterName = pos + name.size();
    const bool endsAtBoundary =
        afterName >= attributes.size()
        || !(std::isalnum(static_cast<unsigned char>(attributes[afterName]))
             || attributes[afterName] == '_');

    if (startsAtBoundary && endsAtBoundary) {
      std::string::size_type equal = afterName;
      while (equal < attributes.size()
             && std::isspace(static_cast<unsigned char>(attributes[equal]))) {
        ++equal;
      }
      if (equal < attributes.size() && attributes[equal] == '=') {
        ++equal;
        while (equal < attributes.size()
               && std::isspace(static_cast<unsigned char>(attributes[equal]))) {
          ++equal;
        }
        if (equal < attributes.size()
            && (attributes[equal] == '"' || attributes[equal] == '\'')) {
          const char quote = attributes[equal];
          const std::string::size_type valueStart = equal + 1;
          const std::string::size_type valueEnd =
              attributes.find(quote, valueStart);
          if (valueEnd != std::string::npos) {
            return attributes.substr(valueStart, valueEnd - valueStart);
          }
        }
      }
    }

    pos = attributes.find(name, pos + name.size());
  }
  return fallback;
}

/**
 * Searches text without case sensitivity.
 *
 * Args:
 *   input: Text to search.
 *   needle: Substring to locate.
 *
 * Returns:
 *   True when the substring occurs.
 */
bool containsCaseInsensitive(const std::string& input,
                             const std::string& needle) {
  return lowercase(input).find(lowercase(needle)) != std::string::npos;
}

/**
 * Extracts the first unsigned integer following a token.
 *
 * Args:
 *   input: Text to inspect.
 *   token: Case-insensitive marker preceding the integer.
 *   value: Output parsed integer.
 *
 * Returns:
 *   True when both token and integer are found.
 */
bool extractIntegerAfterToken(const std::string& input,
                              const std::string& token,
                              int& value) {
  const std::string lowerInput = lowercase(input);
  const std::string lowerToken = lowercase(token);
  const std::string::size_type tokenPos = lowerInput.find(lowerToken);
  if (tokenPos == std::string::npos) return false;

  std::string::size_type pos = tokenPos + lowerToken.size();
  while (pos < input.size()
         && !std::isdigit(static_cast<unsigned char>(input[pos]))) {
    ++pos;
  }
  if (pos == input.size()) return false;

  std::string::size_type end = pos;
  while (end < input.size()
         && std::isdigit(static_cast<unsigned char>(input[end]))) {
    ++end;
  }

  value = parseInteger(input.substr(pos, end - pos), token);
  return true;
}

/**
 * Derives a stable PDF-member leaf name when metadata identifies a PDF weight.
 *
 * Args:
 *   groupName: LHE weight-group name.
 *   description: Individual weight description.
 *   weightIndex: Position used when no member number is present.
 *
 * Returns:
 *   PDFMember name or an empty string for non-PDF weights.
 */
std::string pdfMemberLeafName(const std::string& groupName,
                              const std::string& description,
                              int weightIndex) {
  const bool looksLikePdfWeight =
      containsCaseInsensitive(groupName, "pdf")
      || containsCaseInsensitive(description, "pdf")
      || containsCaseInsensitive(description, "member");

  if (!looksLikePdfWeight) return "";

  int member = -1;
  if (!extractIntegerAfterToken(description, "member", member)) {
    member = weightIndex;
  }

  return "PDFMember" + std::to_string(member);
}

/**
 * Locates one XML-style opening tag and extracts its attributes.
 *
 * Args:
 *   text: Text containing tags.
 *   tag: Tag name without angle brackets.
 *   start: Search offset.
 *   tagEnd: Output offset of the closing angle bracket.
 *   attributes: Output raw attribute text.
 *
 * Returns:
 *   True when a complete opening tag is found.
 */
bool parseOpeningTag(const std::string& text,
                     const std::string& tag,
                     std::string::size_type start,
                     std::string::size_type& tagEnd,
                     std::string& attributes) {
  const std::string open = "<" + tag;
  start = text.find(open, start);
  if (start == std::string::npos) return false;

  const std::string::size_type afterOpen = start + open.size();
  if (afterOpen < text.size()
      && !std::isspace(static_cast<unsigned char>(text[afterOpen]))
      && text[afterOpen] != '>') {
    return parseOpeningTag(text, tag, afterOpen, tagEnd, attributes);
  }

  tagEnd = text.find('>', afterOpen);
  if (tagEnd == std::string::npos) return false;
  attributes = text.substr(afterOpen, tagEnd - afterOpen);
  return true;
}

/**
 * Parses initrwgt metadata from an LHE header.
 *
 * Args:
 *   lheFile: Input LHE path.
 *
 * Returns:
 *   Ordered weight groups with ROOT-safe names.
 */
std::vector<WeightGroup> parseLheWeightGroups(const std::string& lheFile) {
  std::ifstream input(lheFile.c_str());
  if (!input) throw std::runtime_error("Could not open LHE file: " + lheFile);

  std::string line;
  std::ostringstream initrwgt;
  bool inInitrwgt = false;

  while (std::getline(input, line)) {
    if (line.find("<initrwgt") != std::string::npos) inInitrwgt = true;
    if (inInitrwgt) initrwgt << line << '\n';
    if (line.find("</initrwgt>") != std::string::npos) break;
  }

  const std::string text = initrwgt.str();
  if (text.empty()) return std::vector<WeightGroup>();

  std::vector<WeightGroup> groups;
  std::set<std::string> usedBranches;

  int groupIndex = 0;
  std::string::size_type groupSearch = 0;

  while (true) {
    std::string groupAttributes;
    std::string::size_type groupOpenEnd = 0;
    if (!parseOpeningTag(text, "weightgroup", groupSearch, groupOpenEnd,
                         groupAttributes)) {
      break;
    }

    const std::string::size_type groupClose =
        text.find("</weightgroup>", groupOpenEnd + 1);
    if (groupClose == std::string::npos) break;

    WeightGroup group;
    group.originalName =
        getAttribute(groupAttributes, "name",
                     "weightgroup_" + std::to_string(groupIndex));
    group.branchName = uniqueName(
        sanitizeRootName(group.originalName, "WeightGroup" + std::to_string(groupIndex)),
        usedBranches);

    std::set<std::string> usedLeaves;
    const std::string groupBody =
        text.substr(groupOpenEnd + 1, groupClose - groupOpenEnd - 1);
    int weightIndex = 0;
    std::string::size_type weightSearch = 0;

    while (true) {
      std::string weightAttributes;
      std::string::size_type weightOpenEnd = 0;
      if (!parseOpeningTag(groupBody, "weight", weightSearch, weightOpenEnd,
                           weightAttributes)) {
        break;
      }

      const std::string::size_type weightClose =
          groupBody.find("</weight>", weightOpenEnd + 1);
      if (weightClose == std::string::npos) break;

      WeightInfo weight;
      weight.weightId =
          getAttribute(weightAttributes, "id", std::to_string(weightIndex));
      weight.description = collapseWhitespace(
          trim(groupBody.substr(weightOpenEnd + 1,
                                weightClose - weightOpenEnd - 1)));
      const std::string pdfLeafName =
          pdfMemberLeafName(group.originalName, weight.description, weightIndex);
      const std::string leafSeed =
          !pdfLeafName.empty()
              ? pdfLeafName
              : (weight.description.empty() ? "weight_" + weight.weightId
                                            : weight.description);
      weight.leafName = uniqueName(
          sanitizeRootName(leafSeed, "Weight" + std::to_string(weightIndex)),
          usedLeaves);
      group.weights.push_back(weight);

      weightSearch = weightClose + std::string("</weight>").size();
      ++weightIndex;
    }

    if (!group.weights.empty()) groups.push_back(group);
    groupSearch = groupClose + std::string("</weightgroup>").size();
    ++groupIndex;
  }

  return groups;
}

/**
 * Allocates runtime ROOT-branch buffers from parsed metadata.
 *
 * Args:
 *   groups: Parsed LHE weight groups.
 *
 * Returns:
 *   Branch buffers and leaf lists.
 */
std::vector<WeightBranch> buildWeightBranches(
    const std::vector<WeightGroup>& groups) {
  std::vector<WeightBranch> branches;
  branches.reserve(groups.size());

  for (std::vector<WeightGroup>::const_iterator group = groups.begin();
       group != groups.end(); ++group) {
    WeightBranch branch;
    branch.group = *group;
    branch.values.assign(group->weights.size(), 0.0);

    std::ostringstream leafList;
    for (std::size_t i = 0; i < group->weights.size(); ++i) {
      if (i != 0) leafList << ':';
      leafList << group->weights[i].leafName << "/F";
    }
    branch.leafList = leafList.str();
    branches.push_back(branch);
  }

  return branches;
}

/**
 * Copies Pythia's ordered weight vector into grouped ROOT buffers.
 *
 * Args:
 *   branches: Output branch buffers.
 *   pythiaWeights: Nominal weight followed by variation weights.
 */
void fillWeightBranches(std::vector<WeightBranch>& branches,
                        const std::vector<double>& pythiaWeights) {
  std::size_t offset = 1;
  for (std::vector<WeightBranch>::iterator branch = branches.begin();
       branch != branches.end(); ++branch) {
    std::fill(branch->values.begin(), branch->values.end(), 0.0);
    for (std::size_t i = 0; i < branch->values.size(); ++i) {
      const std::size_t pythiaIndex = offset + i;
      if (pythiaIndex < pythiaWeights.size()) {
        branch->values[i] = static_cast<float>(pythiaWeights[pythiaIndex]);
      }
    }
    offset += branch->values.size();
  }
}

/**
 * Writes the LHE-to-ROOT weight-name mapping.
 *
 * Args:
 *   rootFile: Destination ROOT file.
 *   branches: Weight metadata and buffers.
 */
void writeWeightMetadata(TFile& rootFile,
                         const std::vector<WeightBranch>& branches) {
  std::ostringstream metadata;

  for (std::vector<WeightBranch>::const_iterator branch = branches.begin();
       branch != branches.end(); ++branch) {
    metadata << "[" << branch->group.branchName << "] "
             << branch->group.originalName << '\n';
    for (std::vector<WeightInfo>::const_iterator weight =
             branch->group.weights.begin();
         weight != branch->group.weights.end(); ++weight) {
      metadata << weight->leafName << ": id=" << weight->weightId << "; "
               << weight->description << '\n';
    }
  }

  TObjString metadataObject(metadata.str().c_str());
  rootFile.cd();
  metadataObject.Write("LHEWeightMetadata");
}

/**
 * Converts one Pythia particle into the stable ROOT leaf layout.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   index: Valid particle index.
 *
 * Returns:
 *   Flat particle record.
 */
ParticleRecord makeParticleRecord(const Pythia8::Event& event, int index) {
  if (index <= 0 || index >= event.size()) {
    throw std::runtime_error("Invalid particle index selected for ROOT output");
  }

  const Pythia8::Particle& particle = event[index];
  ParticleRecord record;
  record.px = static_cast<float>(particle.px());
  record.py = static_cast<float>(particle.py());
  record.pz = static_cast<float>(particle.pz());
  record.pT = static_cast<float>(particle.pT());
  record.p = static_cast<float>(
      std::sqrt(particle.px() * particle.px() + particle.py() * particle.py()
                + particle.pz() * particle.pz()));
  record.eta = static_cast<float>(particle.eta());
  record.e = static_cast<float>(particle.e());
  record.phi = static_cast<float>(particle.phi());
  record.m0 = static_cast<float>(particle.m0());
  record.pid = static_cast<float>(particle.id());
  record.charge = static_cast<float>(particle.charge());
  record.status = static_cast<float>(particle.status());
  return record;
}

/**
 * Follows consecutive same-ID radiation copies.
 *
 * Args:
 *   event: Complete Pythia event record.
 *   index: Earliest particle-copy index.
 *
 * Returns:
 *   Last same-ID descendant index.
 */
int lastSameIdCopy(const Pythia8::Event& event, int index) {
  if (index <= 0 || index >= event.size()) return index;

  const Pythia8::Particle& particle = event[index];
  const std::vector<int> daughters = particle.daughterList();
  for (std::vector<int>::const_iterator daughter = daughters.begin();
       daughter != daughters.end(); ++daughter) {
    if (*daughter > 0 && *daughter < event.size()
        && event[*daughter].id() == particle.id()) {
      return lastSameIdCopy(event, *daughter);
    }
  }
  return index;
}

/** Provides the out-of-line polymorphic selector destructor. */
EventSelector::~EventSelector() {}

/** Returns success because the default selector adds no Pythia settings. */
bool EventSelector::configure(Pythia8::Pythia&) const { return true; }

/** Returns no metadata key by default. */
std::string EventSelector::metadataName() const { return ""; }

/** Returns no metadata text by default. */
std::string EventSelector::metadataText() const { return ""; }

/**
 * Executes common Pythia initialization, matching, selection, and ROOT output.
 *
 * Args:
 *   argc: Argument count.
 *   argv: Argument values.
 *   selector: Process-specific configuration and particle selection.
 *
 * Returns:
 *   Zero on success or one after a reported exception.
 */
int run(int argc, char* argv[], const EventSelector& selector) {
  try {
    const Options options = parseArgs(argc, argv);

    if (!fileExists(options.lheFile)) {
      throw std::runtime_error("LHE file does not exist: " + options.lheFile);
    }
    if (!fileExists(options.commandFile)) {
      throw std::runtime_error("Pythia command file does not exist: "
                               + options.commandFile);
    }

    if (options.verbose) {
      std::cerr << "Reading LHE weight metadata from " << options.lheFile
                << std::endl;
    }
    std::vector<WeightBranch> weightBranches =
        buildWeightBranches(parseLheWeightGroups(options.lheFile));

    if (options.verbose) {
      std::cerr << "Configuring Pythia with " << options.commandFile
                << std::endl;
    }
    Pythia8::Pythia pythia;
    if (!pythia.readFile(options.commandFile)) {
      throw std::runtime_error("Could not read Pythia command file: "
                               + options.commandFile);
    }

    pythia.readString("Beams:LHEF = " + options.lheFile);

    if (!selector.configure(pythia)) {
      throw std::runtime_error("Could not apply process-specific Pythia settings");
    }

    if (options.seedProvided) {
      if (options.seed <= 0) {
        throw std::runtime_error("--seed must be positive");
      }
      pythia.readString("Random:setSeed = on");
      pythia.readString("Random:seed = " + std::to_string(options.seed));
    }

    if (options.numberEvents > 0) {
      pythia.readString("Main:numberOfEvents = "
                        + std::to_string(options.numberEvents));
    }

    // The hook must stay alive until the event loop has completed.
    if (options.verbose) {
      std::cerr << "Installing CombineMatchingInput user hook" << std::endl;
    }
    Pythia8::CombineMatchingInput jetMatchingHook;
    jetMatchingHook.setHook(pythia);

    if (options.verbose) {
      std::cerr << "Initializing Pythia" << std::endl;
    }
    if (!pythia.init()) {
      throw std::runtime_error("Pythia failed to initialize");
    }

    const int nEvent = pythia.mode("Main:numberOfEvents");
    const bool runToEof = nEvent < 0;
    const int nAbort = pythia.mode("Main:timesAllowErrors");
    int iAbort = 0;
    int iEvent = 0;
    int rejectedEvents = 0;

    if (options.verbose) {
      std::cerr << "Creating ROOT output " << options.outputFile << std::endl;
    }
    TFile rootFile(options.outputFile.c_str(), "RECREATE");
    if (rootFile.IsZombie()) {
      throw std::runtime_error("Could not create ROOT output file: "
                               + options.outputFile);
    }

    TTree tree("Tree", "Tree");

    float eventValue = 0.0;
    float nominalWeight = 0.0;
    float pythiaWeight = 0.0;
    InitialConditionsRecord initialConditions;
    ParticleRecord bornLepton;
    ParticleRecord bornAntiLepton;
    ParticleRecord bareLepton;
    ParticleRecord bareAntiLepton;

    const char* particleLeafList =
        "px/F:py/F:pz/F:pT/F:p/F:eta/F:e/F:phi/F:m0/F:pid/F:charge/F:status/F";

    tree.Branch("Event", &eventValue, "Event/F");
    tree.Branch("InitialConditions", &initialConditions,
                "x1/F:id1/F:x2/F:id2/F:Q2/F");
    tree.Branch("BornLepton", &bornLepton, particleLeafList);
    tree.Branch("BornAntiLepton", &bornAntiLepton, particleLeafList);
    tree.Branch("BareLepton", &bareLepton, particleLeafList);
    tree.Branch("BareAntiLepton", &bareAntiLepton, particleLeafList);
    tree.Branch("NominalWeight", &nominalWeight, "NominalWeight/F");
    tree.Branch("PythiaWeight", &pythiaWeight, "PythiaWeight/F");

    for (std::vector<WeightBranch>::iterator branch = weightBranches.begin();
         branch != weightBranches.end(); ++branch) {
      tree.Branch(branch->group.branchName.c_str(),
                  branch->values.empty() ? 0 : &branch->values[0],
                  branch->leafList.c_str());
    }

    if (options.verbose) {
      std::cerr << "Starting event loop" << std::endl;
    }
    while (runToEof || iEvent < nEvent) {
      if (!pythia.next()) {
        if (pythia.info.atEndOfFile()) break;
        ++iAbort;
        if (iAbort > nAbort) {
          throw std::runtime_error("Pythia stopped after too many failed events");
        }
        continue;
      }

      std::vector<double> weightValues = pythia.info.weightValueVector();

      SelectedLeptons selected;
      std::string failureReason;
      if (!selector.select(pythia.event, selected, failureReason)) {
        ++rejectedEvents;
        std::cerr << "WARNING: event " << (iEvent + rejectedEvents - 1)
                  << " rejected: " << failureReason << '\n';
        continue;
      }

      eventValue = static_cast<float>(iEvent);
      nominalWeight = !weightValues.empty()
                          ? static_cast<float>(weightValues[0])
                          : static_cast<float>(pythia.info.weight());
      pythiaWeight = static_cast<float>(pythia.info.weightValueByIndex(0));

      // These are the flavours, momentum fractions, and factorization scale
      // used for the hard-process PDFs. They are the appropriate quantities
      // for evaluating PDF reweighting factors event by event.
      initialConditions.x1 = static_cast<float>(pythia.info.x1pdf());
      initialConditions.id1 = static_cast<float>(pythia.info.id1pdf());
      initialConditions.x2 = static_cast<float>(pythia.info.x2pdf());
      initialConditions.id2 = static_cast<float>(pythia.info.id2pdf());
      initialConditions.Q2 = static_cast<float>(pythia.info.Q2Fac());

      bornLepton = makeParticleRecord(pythia.event, selected.bornMinus);
      bornAntiLepton = makeParticleRecord(pythia.event, selected.bornPlus);
      bareLepton = makeParticleRecord(pythia.event, selected.bareMinus);
      bareAntiLepton = makeParticleRecord(pythia.event, selected.barePlus);
      fillWeightBranches(weightBranches, weightValues);

      tree.Fill();
      ++iEvent;
    }

    if (options.verbose) {
      std::cerr << "Writing metadata and finalizing ROOT file" << std::endl;
    }
    writeWeightMetadata(rootFile, weightBranches);
    if (!selector.metadataName().empty()) {
      TObjString metadata(selector.metadataText().c_str());
      rootFile.cd();
      metadata.Write(selector.metadataName().c_str());
    }
    std::cout << "Selection summary: wrote " << iEvent << " events; rejected "
              << rejectedEvents << " events.\n";
    pythia.stat();
    rootFile.Write();
    rootFile.Close();
    if (options.verbose) {
      std::cerr << "Done" << std::endl;
    }
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }

  return 0;
}

}  // namespace pythia_lhe
