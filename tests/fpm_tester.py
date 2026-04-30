#!/usr/bin/env python3

from __future__ import annotations

import argparse
from contextlib import contextmanager
from typing import TextIO
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STANDALONE_TEST_DIR = PROJECT_ROOT / "tests" / "standalone_tests"
STANDALONE_LOG_DIR = PROJECT_ROOT / "build" / "fpm_standalone_tests"
STANDALONE_THRESHOLD = 1e-4
DEFAULT_LOG_FILE = PROJECT_ROOT / "tests" / "fpm_test.log"
FILTERED_LINES = {"Project is up to date"}

TESTS = [
    ("ddx_core", []),
    ("ddx_operators", []),
    ("bessel", []),
    ("force", ["tests/Input_force.txt"]),
    (
        "ddx_driver",
        [
            "tests/data/ddpcm_force_fmm.in",
            "tests/data/ddpcm_force_fmm.out",
            "1E-12",
        ],
    ),
    (
        "ddx_driver",
        [
            "tests/data/ddcosmo_force_fmm.in",
            "tests/data/ddcosmo_force_fmm.out",
            "1E-12",
        ],
    ),
    ("force_ddlpb", ["tests/data/ddlpb_force.txt"]),
    ("ddlpb_esolv", ["tests/data/ddlpb_force.txt"]),
    ("matrix_derivatives", ["tests/data/ddlpb_force.txt"]),
    ("matrix_adjoint", ["tests/data/ddlpb_force.txt"]),
    ("matrix_solvers", ["tests/data/ddlpb_force.txt"]),
    ("m2l", []),
    ("multipolar_solutes", []),
    ("error", []),
    ("test_gradients", ["tests/Input_cosmo_small.txt"]),
]

STANDALONE_TESTS = [
    "cosmo",
    "cosmo_fmm",
    "cosmo_incore",
    "pcm",
    "pcm_fmm",
    "pcm_incore",
    "lpb",
    "lpb_fmm",
    "lpb_incore",
]


def format_test_name(name: str, args: list[str]) -> str:
    if not args:
        return name
    return f"{name} -- {' '.join(args)}"


def colorize(text: str, color: str) -> str:
    colors = {
        "green": "\033[38;2;0;170;80m",
        "red": "\033[38;2;220;40;40m",
        "cyan": "\033[38;2;0;160;200m",
        "reset": "\033[0m",
    }
    if not status_is_terminal():
        return text
    return f"{colors[color]}{text}{colors['reset']}"


def print_result(label: str, status: int, message: str = "") -> None:
    result = "PASS" if status == 0 else "FAIL"
    color = "green" if status == 0 else "red"
    detail = f": {message}" if message else ""
    print_status(f"{colorize(result, color)}: {label}{detail}")


def status_is_terminal() -> bool:
    if sys.stdout.isatty():
        return True
    try:
        with Path("/dev/tty").open("w", encoding="utf-8"):
            return True
    except OSError:
        return False


@contextmanager
def status_stream() -> TextIO:
    if sys.stdout.isatty():
        yield sys.stdout
        return

    try:
        with Path("/dev/tty").open("w", encoding="utf-8") as terminal:
            yield terminal
    except OSError:
        yield sys.stderr


def print_status(message: str = "") -> None:
    with status_stream() as stream:
        print(message, file=stream, flush=True)


def filter_output(output: str) -> str:
    lines = [
        line
        for line in output.splitlines()
        if line.strip() not in FILTERED_LINES
    ]
    if output.endswith("\n"):
        return "\n".join(lines) + ("\n" if lines else "")
    return "\n".join(lines)


def write_test_output(
    log_file,
    cmd: list[str],
    output: str,
    verbose: bool,
) -> None:
    header = "\n" + "=" * 80 + "\n" + " ".join(cmd) + "\n" + "=" * 80 + "\n"
    log_file.write(header)
    log_file.write(output)
    log_file.flush()

    if verbose:
        print(header, end="")
        print(output, end="")


def run_command(cmd: list[str], log_file, verbose: bool) -> int:
    result = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    write_test_output(log_file, cmd, filter_output(result.stdout), verbose)
    return result.returncode


def run_test(name: str, args: list[str], log_file, verbose: bool) -> int:
    cmd = ["fpm", "test", "--target", name]
    if args:
        cmd.extend(["--", *args])

    return run_command(cmd, log_file, verbose)


def read_driver_log(path: Path) -> tuple[float, list[list[float]]]:
    energy = 0.0
    forces = []
    section = ""

    with path.open("r", encoding="utf-8") as log:
        for line in log:
            if "Solvation energy (Hartree):" in line:
                tokens = line.split()
                energy = float(tokens[3])
            elif section == "forces":
                tokens = line.split()
                if len(tokens) >= 4:
                    forces.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
            elif "Full forces (kcal/mol/A)" in line:
                section = "forces"

    return energy, forces


def inf_norm(values: list[list[float]]) -> float:
    return max((abs(value) for row in values for value in row), default=0.0)


def relative_error(value: float, reference: float) -> float:
    scale = abs(reference) if reference != 0.0 else 1.0
    return abs(value - reference) / scale


def compare_driver_logs(output_file: Path, ref_file: Path) -> tuple[bool, str]:
    energy, forces = read_driver_log(output_file)
    ref_energy, ref_forces = read_driver_log(ref_file)

    energy_error = relative_error(energy, ref_energy)
    if energy_error >= STANDALONE_THRESHOLD:
        return (
            False,
            f"energy relative error {energy_error:.3e} exceeds "
            f"{STANDALONE_THRESHOLD:.3e}",
        )

    if len(forces) != len(ref_forces):
        return False, f"force count {len(forces)} differs from reference {len(ref_forces)}"

    force_diff = [
        [force_value - ref_value for force_value, ref_value in zip(force, ref_force)]
        for force, ref_force in zip(forces, ref_forces)
    ]
    force_error = inf_norm(force_diff) / max(inf_norm(ref_forces), 1.0)
    if force_error >= STANDALONE_THRESHOLD:
        return (
            False,
            f"force relative error {force_error:.3e} exceeds "
            f"{STANDALONE_THRESHOLD:.3e}",
        )

    return True, "passed"


def run_standalone_test(name: str, log_file, verbose: bool) -> tuple[int, str]:
    input_file = STANDALONE_TEST_DIR / f"{name}.txt"
    output_file = STANDALONE_LOG_DIR / f"{name}.log"
    ref_file = STANDALONE_TEST_DIR / f"{name}.ref"
    input_arg = str(input_file.relative_to(PROJECT_ROOT))
    cmd = ["fpm", "run", "--target", "ddx_driver", "--", input_arg]

    result = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    output = filter_output(result.stdout)
    write_test_output(log_file, cmd, output, verbose)

    STANDALONE_LOG_DIR.mkdir(parents=True, exist_ok=True)
    output_file.write_text(output, encoding="utf-8")

    if result.returncode != 0:
        return result.returncode, f"exited with status {result.returncode}"

    passed, message = compare_driver_logs(output_file, ref_file)
    return (0, message) if passed else (1, message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the complete ddX fpm test suite.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""\
Why this script exists:
  A plain 'fpm test' runs the test executables declared as [[test]] targets in
  fpm.toml, but it does not know which command-line input files/arguments some
  of the ddX tests require. This script runs each fpm test target with the
  expected arguments.

  The ddX project also has standalone driver tests, which run ddx_driver on 
  input fixtures and compare the resulting energies and forces against 
  reference output.

  This script provides one fpm-facing command that runs both groups:
    1. unit-style fpm test targets
    2. standalone ddx_driver fixture/reference tests

Output:
  By default, the terminal shows compact progress, PASS/FAIL results, and a
  final summary. Full test output is written to a log file.

Options:
  --verbose
      Also print the full test output to the terminal while still writing it
      to the log file.

  --log-file PATH
      Write full test output to PATH instead of the default log file.
      Default: {DEFAULT_LOG_FILE.relative_to(PROJECT_ROOT)}
""",
    )
    parser.add_argument(
        "--log-file",
        default=str(DEFAULT_LOG_FILE.relative_to(PROJECT_ROOT)),
        help="write full test output to this file (default: %(default)s)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="also print full test output to the terminal",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_path = Path(args.log_file)
    if not log_path.is_absolute():
        log_path = PROJECT_ROOT / log_path
    log_path.parent.mkdir(parents=True, exist_ok=True)

    results = []

    print_status(f"Full output is written to {log_path}")
    if args.verbose:
        print_status("Verbose mode enabled; full output is also printed below.")

    with log_path.open("w", encoding="utf-8") as log_file:
        for name, test_args in TESTS:
            label = format_test_name(name, test_args)
            print_status(f"{colorize('RUNNING', 'cyan')}: {label}")
            status = run_test(name, test_args, log_file, args.verbose)
            message = "passed" if status == 0 else f"exited with status {status}"
            print_result(label, status, message)
            results.append(("fpm", label, status, message))

        for name in STANDALONE_TESTS:
            label = f"standalone {name}"
            print_status(f"{colorize('RUNNING', 'cyan')}: {label}")
            status, message = run_standalone_test(name, log_file, args.verbose)
            print_result(label, status, message)
            results.append(("standalone", label, status, message))

    failures = [result for result in results if result[2] != 0]
    passed = len(results) - len(failures)

    print_status()
    print_status("Test summary:")
    print_status(f"  Total:      {len(results)}")
    print_status(f"  Passed:     {passed}")
    print_status(f"  Failed:     {len(failures)}")
    print_status(f"  Unit:       {sum(1 for kind, *_ in results if kind == 'fpm')}")
    print_status(
        f"  standalone: {sum(1 for kind, *_ in results if kind == 'standalone')}"
    )
    print_status(f"  log file:   {log_path}")

    if failures:
        print_status()
        print_status("Failed tests:")
        for _, label, _, message in failures:
            print_status(f"  {label}: {message}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
