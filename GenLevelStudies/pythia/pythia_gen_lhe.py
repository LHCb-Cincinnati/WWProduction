"""Shower an MG5_aMC@NLO LHE file with Pythia8 and write weights to ROOT."""

import argparse
import gzip
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import ROOT

DEFAULT_CMD = Path(__file__).with_name("MG5NLOShower.cmd")
# ROOT = None


@dataclass
class WeightInfo:
    """Description of one LHE weight from the initrwgt header."""

    weight_id: str
    leaf_name: str
    description: str


@dataclass
class WeightGroup:
    """One LHE weight group to become one ROOT branch."""

    branch_name: str
    original_name: str
    weights: list[WeightInfo]


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Shower an MG5_aMC@NLO LHE file with Pythia8."
    )
    parser.add_argument("lhe_file", help="Input LHE or LHE.gz file.")
    parser.add_argument("output_file", help="Output ROOT file.")
    parser.add_argument(
        "-c",
        "--command-file",
        default=str(DEFAULT_CMD),
        help="Pythia command file. Defaults to GenLevelStudies/pythia/MG5NLOShower.cmd.",
    )
    parser.add_argument(
        "-n",
        "--number-events",
        type=int,
        default=None,
        help="Override Main:numberOfEvents. Use a negative value to run to LHE EOF.",
    )
    return parser.parse_args(argv)


def open_text_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def sanitize_root_name(name, fallback):
    clean = re.sub(r"[^0-9A-Za-z_]", "_", name.strip())
    clean = re.sub(r"_+", "_", clean).strip("_")
    if not clean:
        clean = fallback
    if clean[0].isdigit():
        clean = "_" + clean
    return clean


def unique_name(name, used):
    if name not in used:
        used.add(name)
        return name
    index = 2
    while f"{name}_{index}" in used:
        index += 1
    unique = f"{name}_{index}"
    used.add(unique)
    return unique


def parse_lhe_weight_groups(lhe_file):
    """Parse LHE3 initrwgt metadata without loading the full event file."""

    with open_text_maybe_gzip(lhe_file) as handle:
        text = []
        in_initrwgt = False
        for line in handle:
            if "<initrwgt" in line:
                in_initrwgt = True
            if in_initrwgt:
                text.append(line)
            if "</initrwgt>" in line:
                break

    if not text:
        return []

    initrwgt = ET.fromstring("".join(text))
    used_branches = set()
    groups = []

    for igroup, group in enumerate(initrwgt.findall("weightgroup")):
        original_name = group.attrib.get("name", f"weightgroup_{igroup}")
        branch_name = unique_name(
            sanitize_root_name(original_name, f"WeightGroup{igroup}"),
            used_branches,
        )
        used_leaves = set()
        weights = []

        for iweight, weight in enumerate(group.findall("weight")):
            weight_id = weight.attrib.get("id", str(iweight))
            description = " ".join((weight.text or "").split())
            leaf_seed = description or f"weight_{weight_id}"
            leaf_name = unique_name(
                sanitize_root_name(leaf_seed, f"Weight{iweight}"),
                used_leaves,
            )
            weights.append(
                WeightInfo(
                    weight_id=weight_id,
                    leaf_name=leaf_name,
                    description=description,
                )
            )

        if weights:
            groups.append(
                WeightGroup(
                    branch_name=branch_name,
                    original_name=original_name,
                    weights=weights,
                )
            )

    return groups


def build_weight_arrays(weight_groups):
    arrays = {}
    leaflists = {}
    for group in weight_groups:
        arrays[group.branch_name] = np.zeros(len(group.weights), dtype=np.float32)
        leaves = [f"{weight.leaf_name}/F" for weight in group.weights]
        leaflists[group.branch_name] = ":".join(leaves)
    return arrays, leaflists


def fill_weight_arrays(weight_groups, arrays, pythia_weight_values):
    """Fill grouped arrays from Pythia weights.

    Pythia's first weight is the nominal "Weight" entry. The detailed
    LHE initrwgt values follow in the same order as the LHE header, even
    when Pythia rewrites their labels.
    """

    offset = 1
    for group in weight_groups:
        values = arrays[group.branch_name]
        values.fill(0.0)
        for index in range(len(group.weights)):
            pythia_index = offset + index
            if pythia_index < len(pythia_weight_values):
                values[index] = pythia_weight_values[pythia_index]
        offset += len(group.weights)


def fill_event_content(pythia, tree):
    """Hook for adding analysis-specific branches later."""

    return pythia, tree


def write_weight_metadata(root_file, weight_groups):
    metadata = []
    for group in weight_groups:
        metadata.append(f"[{group.branch_name}] {group.original_name}")
        for weight in group.weights:
            metadata.append(
                f"{weight.leaf_name}: id={weight.weight_id}; {weight.description}"
            )
    root_file.WriteObject(
        ROOT.TObjString("\n".join(metadata)),
        "LHEWeightMetadata",
    )


def configure_pythia(command_file, lhe_file, number_events):
    import pythia8

    pythia = pythia8.Pythia()
    pythia.readFile(str(command_file))
    pythia.readString(f"Beams:LHEF = {lhe_file}")
    if number_events is not None:
        pythia.readString(f"Main:numberOfEvents = {number_events}")
    return pythia


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    lhe_file = Path(args.lhe_file).expanduser()
    output_file = Path(args.output_file).expanduser()
    command_file = Path(args.command_file).expanduser()

    if not lhe_file.exists():
        raise FileNotFoundError(f"LHE file does not exist: {lhe_file}")
    if not command_file.exists():
        raise FileNotFoundError(f"Pythia command file does not exist: {command_file}")

    weight_groups = parse_lhe_weight_groups(lhe_file)
    weight_arrays, weight_leaflists = build_weight_arrays(weight_groups)

    pythia = configure_pythia(command_file, lhe_file, args.number_events)
    if not pythia.init():
        raise RuntimeError("Pythia failed to initialize.")

    n_event = pythia.mode("Main:numberOfEvents")
    run_to_eof = n_event < 0

    output_file.parent.mkdir(parents=True, exist_ok=True)
    root_file = ROOT.TFile.Open(str(output_file), "RECREATE")
    tree = ROOT.TTree("Tree", "Tree")

    event_array = np.zeros(1, dtype=np.float32)
    nominal_weight_array = np.zeros(1, dtype=np.float32)
    pythia_weight_array = np.zeros(1, dtype=np.float32)

    tree.Branch("Event", event_array, "Event/F")
    tree.Branch("NominalWeight", nominal_weight_array, "NominalWeight/F")
    tree.Branch("PythiaWeight", pythia_weight_array, "PythiaWeight/F")
    for group in weight_groups:
        tree.Branch(
            group.branch_name,
            weight_arrays[group.branch_name],
            weight_leaflists[group.branch_name],
        )

    i_event = 0
    i_abort = 0
    n_abort = pythia.mode("Main:timesAllowErrors")

    while run_to_eof or i_event < n_event:
        if not pythia.next():
            if pythia.infoPython().atEndOfFile():
                break
            i_abort += 1
            if i_abort > n_abort:
                raise RuntimeError("Pythia stopped after too many failed events.")
            continue

        info = pythia.infoPython()
        weight_values = list(info.weightValueVector())

        event_array[0] = i_event
        nominal_weight_array[0] = weight_values[0] if weight_values else info.weight()
        pythia_weight_array[0] = info.weightValueByIndex()
        fill_weight_arrays(weight_groups, weight_arrays, weight_values)
        fill_event_content(pythia, tree)
        tree.Fill()
        i_event += 1

    write_weight_metadata(root_file, weight_groups)
    pythia.stat()
    tree.Print()
    root_file.Write()
    root_file.Close()


if __name__ == "__main__":
    main()
