from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
WOLFRAM_DIR = PROJECT_ROOT / "wolfram"
OUTPUT_DIR = WOLFRAM_DIR / "output"
RUN_MODEL_SCRIPT = WOLFRAM_DIR / "runners" / "RunModel.wl"

EFT_ORDER = 5
LOOP_ORDER = 1

# Historical benchmarks used to check the generic builder against validated
# singlet/doublet/triplet T3 models.
INTERESTING = [("A", 0), ("B", -1), ("C", -1), ("D", -2), ("E", 0)]
SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]
EXTENDED = [
    ("A", 0), ("A", -2), ("B", -1), ("C", -1),
    ("D", -2), ("E", 0), ("E", -2),
]


@dataclass
class RunRecord:
    name: str
    alpha: int
    return_code: int
    summary: dict
    output_dir: Path


def encode_alpha(alpha: int) -> str:
    """Encode signed alpha without a leading '-' for Wolfram CLI parsing."""
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def valid_t3_dimensions(d_s1: int, d_s2: int, d_f: int) -> bool:
    """Check the SU(2)-representation conditions required by T3.

    Each scalar must occur in 2⊗dF, so dS=dF±1. The external HH pair is in
    the symmetric triplet channel, so S1⊗S2 must contain total isospin J=1.
    """
    if min(d_s1, d_s2, d_f) < 1:
        return False
    if abs(d_s1 - d_f) != 1 or abs(d_s2 - d_f) != 1:
        return False

    # d=2j+1; apply the usual SU(2) angular-momentum addition rule for J=1.
    j1 = (d_s1 - 1) / 2
    j2 = (d_s2 - 1) / 2
    return abs(j1 - j2) <= 1 <= j1 + j2 and float(j1 + j2).is_integer()


def run_model(name: str, alpha: int, output_dir: Path, model_args: list[str], verbose: bool) -> RunRecord:
    """Run one Wolfram model and read the summary produced by RunModel.wl."""
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "wolframscript",
        "-file",
        str(RUN_MODEL_SCRIPT),
        str(output_dir),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        *model_args,
    ]
    process = subprocess.run(
        command,
        cwd=WOLFRAM_DIR,
        capture_output=True,
        text=True,
        check=False,
    )

    (output_dir / "wolfram_stdout.log").write_text(process.stdout, encoding="utf-8")
    (output_dir / "wolfram_stderr.log").write_text(process.stderr, encoding="utf-8")

    if verbose:
        if process.stdout:
            print(process.stdout)
        if process.stderr:
            print(process.stderr, file=sys.stderr)
    elif process.returncode != 0:
        # Keep normal scans compact while retaining enough context to diagnose a failure.
        if process.stdout:
            print("\n".join(process.stdout.splitlines()[-40:]))
        if process.stderr:
            print("\n".join(process.stderr.splitlines()[-20:]), file=sys.stderr)

    summary_path = output_dir / "comparison_summary.json"
    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {"BuildStatus": "ProcessFailed", "MatchingStatus": "NotRun"}
    )
    return RunRecord(name, alpha, process.returncode, summary, output_dir)


def run_class(model_class: str, alpha: int, verbose: bool = False) -> RunRecord:
    name = f"T3-{model_class}"
    output_dir = OUTPUT_DIR / f"T3_{model_class}_alpha_{encode_alpha(alpha)}"
    print(f"Running {name}, alpha={alpha} ...", flush=True)
    return run_model(name, alpha, output_dir, [model_class, encode_alpha(alpha)], verbose)


def run_dimensions(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    verbose: bool = False,
) -> RunRecord:
    if not valid_t3_dimensions(d_s1, d_s2, d_f):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) is not a valid T3 SU(2) assignment: "
            "each scalar must have dS=dF±1 and S1⊗S2 must contain the triplet."
        )

    name = f"T3-d{d_s1}-d{d_s2}-F{d_f}"
    output_dir = OUTPUT_DIR / (
        f"T3_d{d_s1}_d{d_s2}_F{d_f}_alpha_{encode_alpha(alpha)}"
    )
    model_args = ["DIMS", str(d_s1), str(d_s2), str(d_f), encode_alpha(alpha)]

    print(f"Running {name}, alpha={alpha} ...", flush=True)
    return run_model(name, alpha, output_dir, model_args, verbose)


def print_summary(records: list[RunRecord]) -> int:
    print("\n" + "=" * 72 + "\nT3 MODEL SUMMARY\n" + "=" * 72)
    successful = 0

    for record in records:
        summary = record.summary
        build_ok = summary.get("BuildStatus") == "Success"
        match_ok = summary.get("MatchingStatus") == "Success"
        successful += int(build_ok and match_ok)

        print(
            f"{record.name} alpha={record.alpha}: "
            f"build={summary.get('BuildStatus')}, match={summary.get('MatchingStatus')}, "
            f"T3={summary.get('T3IngredientsPresent')}, "
            f"Weinberg={summary.get('WeinbergOperatorPresent')}"
        )

        if not summary.get("WeinbergOperatorPresent"):
            continue

        extraction = summary.get("WeinbergExtractionStatus", "Unknown")
        if extraction == "Success":
            n_holo = summary.get("WeinbergHolomorphicTermCount", 0)
            n_hc = summary.get("WeinbergConjugateTermCount", 0)
            print(f"  C5 extraction: Success ({n_holo} holomorphic + {n_hc} HC terms)")
            if coefficient_file := summary.get("WeinbergCoefficientFile"):
                print(f"  C5: {record.output_dir / coefficient_file}")
        else:
            print(
                f"  Weinberg terms: {summary.get('WeinbergTermCount', 0)}; "
                f"C5: {extraction}"
            )

    aggregate = OUTPUT_DIR / "t3_model_comparison.json"
    aggregate.write_text(
        json.dumps([record.summary for record in records], indent=2),
        encoding="utf-8",
    )
    print(f"\n{successful}/{len(records)} completed build+matching.\nAggregate: {aggregate}")
    return 0 if successful == len(records) else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run T3 matching for legacy benchmarks or arbitrary valid SU(2) irreps."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--smoke", action="store_true", help="five historical T3 regression models")
    mode.add_argument("--extended", action="store_true", help="seven historical benchmark points")
    mode.add_argument(
        "--dims",
        nargs=3,
        type=int,
        metavar=("DS1", "DS2", "DF"),
        help="run one generic representation assignment, e.g. --dims 3 5 4",
    )
    parser.add_argument(
        "--alpha",
        type=int,
        default=0,
        help="T3 hypercharge parameter for --dims mode (default: 0)",
    )
    parser.add_argument("--verbose", action="store_true", help="show full Wolfram/Matchete output")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.dims:
        d_s1, d_s2, d_f = args.dims
        try:
            record = run_dimensions(d_s1, d_s2, d_f, args.alpha, verbose=args.verbose)
        except ValueError as exc:
            parser.error(str(exc))
        return print_summary([record])

    if args.smoke:
        mode_name, points = "smoke", SMOKE
    elif args.extended:
        mode_name, points = "extended", EXTENDED
    else:
        mode_name, points = "interesting", INTERESTING

    print(f"T3 scan mode: {mode_name}; {len(points)} model(s).")
    records = [run_class(model_class, alpha, args.verbose) for model_class, alpha in points]
    return print_summary(records)


if __name__ == "__main__":
    raise SystemExit(main())
