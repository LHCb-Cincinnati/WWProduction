#!/usr/bin/env python3
"""Run many independent LHE showering jobs with a selected C++ executable.

This wrapper intentionally does not merge the produced ROOT files.  Each input
LHE file is showered by one serial C++ process with its own deterministic seed.

Compile the showering executables:
    conda activate hep
    make -C GenLevelStudies/pythia -f Makefile.pythia_gen_lhe all

Run the NLO dilepton shower:
    python GenLevelStudies/pythia/run_lhe_shower_jobs.py \
        .agents/skills/showering/assets/ex_WW.lhe \
        --shower-script GenLevelStudies/pythia/pythia_gen_lhe_nlo \
        -o /tmp/showered_nlo -j 1

Run the forced tau-decay shower:
    python GenLevelStudies/pythia/run_lhe_shower_jobs.py \
        .agents/skills/showering/assets/Ex_Z_tautau_NLO.lhe \
        --shower-script GenLevelStudies/pythia/pythia_gen_lhe_tautau \
        -o /tmp/showered_tautau -j 1

Run an unmatched LO WZ shower using the executable's LO command-file default:
    python GenLevelStudies/pythia/run_lhe_shower_jobs.py \
        .agents/skills/showering/assets/Ex_WZ.lhe \
        --shower-script GenLevelStudies/pythia/pythia_gen_lhe_wz_lo \
        -o /tmp/showered_wz -j 1
"""

import argparse
import concurrent.futures
import os
import re
import subprocess
import sys
from pathlib import Path


def executable_path(value):
    """Validate and normalize a showering executable path for argparse.

    Args:
        value: User-provided filesystem path.

    Returns:
        Absolute path to an existing executable file.

    Raises:
        argparse.ArgumentTypeError: The path is absent, not a file, or not
            executable by the current user.
    """
    path = Path(value).expanduser()
    if not path.is_file():
        raise argparse.ArgumentTypeError(
            f"showering script does not exist or is not a file: {path}"
        )
    if not os.access(path, os.X_OK):
        raise argparse.ArgumentTypeError(
            f"showering script is not executable: {path}"
        )
    return str(path.resolve())


def parse_args(argv):
    """Parse launcher arguments.

    Args:
        argv: Command-line arguments excluding the program name.

    Returns:
        Parsed argparse namespace.
    """
    parser = argparse.ArgumentParser(
        description="Launch one serial, selectable C++ shower per LHE file."
    )
    parser.add_argument(
        "lhe_files",
        nargs="+",
        help="Input LHE files. Shell globs should be expanded by the shell.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory where per-LHE ROOT files will be written.",
    )
    parser.add_argument(
        "-s",
        "--shower-script",
        type=executable_path,
        default=str(Path(__file__).with_name("pythia_gen_lhe_nlo")),
        metavar="EXECUTABLE",
        help=(
            "Path to an existing executable showering script "
            "(default: pythia_gen_lhe_nlo beside this launcher)."
        ),
    )
    parser.add_argument(
        "-c",
        "--command-file",
        default=None,
        help=(
            "Optional Pythia command file passed to every job. When omitted, "
            "the selected executable uses its process-specific default."
        ),
    )
    parser.add_argument(
        "-n",
        "--number-events",
        type=int,
        default=-1,
        help="Number of accepted events per LHE file; negative means run to EOF.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=os.cpu_count() or 1,
        help="Maximum number of showering processes to run at once.",
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=12345,
        help="Seed for the first file. Later files use base_seed + file_index.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running them.",
    )
    return parser.parse_args(argv)


def sanitize_name(name):
    """Convert a filename component into a portable ROOT output stem.

    Args:
        name: Original filename component.

    Returns:
        Sanitized nonempty name.
    """
    clean = re.sub(r"[^0-9A-Za-z_]", "_", name)
    clean = re.sub(r"_+", "_", clean).strip("_")
    return clean or "lhe"


def output_name_for(lhe_file, used_names):
    """Build a collision-free ROOT filename for one LHE input.

    Args:
        lhe_file: Input LHE path.
        used_names: Set of output stems already assigned.

    Returns:
        Unique ROOT filename.
    """
    path = Path(lhe_file)
    stem = path.name
    if stem.endswith(".lhe"):
        stem = stem[:-4]

    parent = sanitize_name(path.parent.name)
    candidate = sanitize_name(f"{parent}_{stem}")
    if candidate in used_names:
        index = 2
        while f"{candidate}_{index}" in used_names:
            index += 1
        candidate = f"{candidate}_{index}"

    used_names.add(candidate)
    return candidate + ".root"


def build_jobs(args):
    """Construct deterministic per-file showering jobs.

    Args:
        args: Parsed launcher arguments.

    Returns:
        List of job tuples.
    """
    output_dir = Path(args.output_dir)
    used_names = set()
    jobs = []

    for index, lhe_file in enumerate(args.lhe_files):
        lhe_path = Path(lhe_file)
        if not lhe_path.exists():
            raise FileNotFoundError(f"LHE file does not exist: {lhe_path}")

        output_file = output_dir / output_name_for(lhe_path, used_names)
        seed = args.base_seed + index

        command = [
            args.shower_script,
            str(lhe_path),
            str(output_file),
            "-n",
            str(args.number_events),
            "--seed",
            str(seed),
        ]
        if args.command_file is not None:
            command[3:3] = ["-c", args.command_file]
        jobs.append((lhe_path, output_file, seed, command))

    return jobs


def run_job(job):
    """Run one shower process while inheriting the launcher's output streams.

    Args:
        job: Tuple containing input, output, seed, and command.

    Returns:
        Job metadata and the completed process.
    """
    lhe_path, output_file, seed, command = job
    output_file.parent.mkdir(parents=True, exist_ok=True)

    completed = subprocess.run(command, check=False)
    return lhe_path, output_file, seed, completed


def main(argv=None):
    """Run all requested showering jobs.

    Args:
        argv: Optional argument list used by tests.

    Returns:
        Zero when every job succeeds, otherwise one.
    """
    args = parse_args(sys.argv[1:] if argv is None else argv)
    jobs = build_jobs(args)

    if args.dry_run:
        for _, _, _, command in jobs:
            print(" ".join(command))
        return 0

    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_job = {executor.submit(run_job, job): job for job in jobs}
        for future in concurrent.futures.as_completed(future_to_job):
            lhe_path, output_file, seed, completed = future.result()
            if completed.returncode == 0:
                print(f"OK seed={seed} input={lhe_path} output={output_file}")
            else:
                failures += 1
                print(
                    f"FAILED seed={seed} input={lhe_path} output={output_file}",
                    file=sys.stderr,
                )

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
