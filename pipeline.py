from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
W = ROOT / "wolfram"
OUT = W / "output"
RUN = W / "RunModel.wl"
EFT_ORDER = 5
LOOP_ORDER = 1

# Historical A-E regression points.  These remain useful because they compare
# the generic builder against the already validated singlet/doublet/triplet
# implementation.
INTERESTING = [("A", 0), ("B", -1), ("C", -1), ("D", -2), ("E", 0)]
SMOKE = [("B", -1), ("C", -1), ("A", 0), ("D", -2), ("E", 0)]
EXTENDED = [("A", 0), ("A", -2), ("B", -1), ("C", -1), ("D", -2), ("E", 0), ("E", -2)]


@dataclass
class Record:
    name: str
    alpha: int
    rc: int
    summary: dict
    path: Path


def alpha_arg(alpha: int) -> str:
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def class_label(cls: str, alpha: int) -> str:
    return f"T3_{cls}_alpha_{alpha_arg(alpha)}"


def dims_label(d1: int, d2: int, df: int, alpha: int) -> str:
    return f"T3_d{d1}_d{d2}_F{df}_alpha_{alpha_arg(alpha)}"


def _run(cmd: list[str], path: Path, verbose: bool) -> tuple[int, dict]:
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(cmd, cwd=W, capture_output=True, text=True, check=False)
    (path / "wolfram_stdout.log").write_text(proc.stdout, encoding="utf-8")
    (path / "wolfram_stderr.log").write_text(proc.stderr, encoding="utf-8")

    if verbose:
        if proc.stdout:
            print(proc.stdout)
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
    elif proc.returncode != 0:
        if proc.stdout:
            print("\n".join(proc.stdout.splitlines()[-40:]))
        if proc.stderr:
            print("\n".join(proc.stderr.splitlines()[-20:]), file=sys.stderr)

    summary_path = path / "comparison_summary.json"
    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {"BuildStatus": "ProcessFailed", "MatchingStatus": "NotRun"}
    )
    return proc.returncode, summary


def run_class(cls: str, alpha: int, verbose: bool = False) -> Record:
    path = OUT / class_label(cls, alpha)
    print(f"Running T3-{cls}, alpha={alpha} ...", flush=True)
    cmd = [
        "wolframscript",
        "-file",
        str(RUN),
        str(path),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        cls,
        alpha_arg(alpha),
    ]
    rc, summary = _run(cmd, path, verbose)
    return Record(f"T3-{cls}", alpha, rc, summary, path)


def valid_t3_dimensions(d1: int, d2: int, df: int) -> bool:
    """Representation-only T3 test.

    Each scalar must occur in 2 x dF, hence dS=dF±1.  The HH pair is in the
    symmetric triplet channel, so S1 x S2 must contain a triplet.
    """
    if min(d1, d2, df) < 1:
        return False
    if abs(d1 - df) != 1 or abs(d2 - df) != 1:
        return False

    # d=2j+1.  Check J=1 belongs to j1 x j2.
    j1 = (d1 - 1) / 2
    j2 = (d2 - 1) / 2
    return abs(j1 - j2) <= 1 <= j1 + j2 and float(j1 + j2).is_integer()


def run_dims(d1: int, d2: int, df: int, alpha: int, verbose: bool = False) -> Record:
    if not valid_t3_dimensions(d1, d2, df):
        raise ValueError(
            f"({d1}, {d2}, {df}) is not a valid T3 SU(2) assignment: "
            "each scalar must have dS=dF±1 and S1⊗S2 must contain the triplet."
        )

    path = OUT / dims_label(d1, d2, df, alpha)
    name = f"T3-d{d1}-d{d2}-F{df}"
    print(f"Running {name}, alpha={alpha} ...", flush=True)
    cmd = [
        "wolframscript",
        "-file",
        str(RUN),
        str(path),
        str(EFT_ORDER),
        str(LOOP_ORDER),
        "DIMS",
        str(d1),
        str(d2),
        str(df),
        alpha_arg(alpha),
    ]
    rc, summary = _run(cmd, path, verbose)
    return Record(name, alpha, rc, summary, path)


def print_summary(records: list[Record]) -> int:
    print("\n" + "=" * 72 + "\nT3 MODEL SUMMARY\n" + "=" * 72)
    ok = 0
    for record in records:
        s = record.summary
        success = s.get("BuildStatus") == "Success" and s.get("MatchingStatus") == "Success"
        ok += int(success)
        print(
            f"{record.name} alpha={record.alpha}: "
            f"build={s.get('BuildStatus')}, match={s.get('MatchingStatus')}, "
            f"T3={s.get('T3IngredientsPresent')}, Weinberg={s.get('WeinbergOperatorPresent')}"
        )
        if s.get("WeinbergOperatorPresent"):
            extraction = s.get("WeinbergExtractionStatus", "Unknown")
            hcount = s.get("WeinbergHolomorphicTermCount", 0)
            hccount = s.get("WeinbergConjugateTermCount", 0)
            if extraction == "Success":
                print(f"  C5 extraction: Success ({hcount} holomorphic + {hccount} HC terms)")
                c5_file = s.get("WeinbergCoefficientFile")
                if c5_file:
                    print(f"  C5: {record.path / c5_file}")
            else:
                print(f"  Weinberg terms: {s.get('WeinbergTermCount', 0)}; C5: {extraction}")

    aggregate = OUT / "t3_model_comparison.json"
    aggregate.write_text(json.dumps([r.summary for r in records], indent=2), encoding="utf-8")
    print(f"\n{ok}/{len(records)} completed build+matching.\nAggregate: {aggregate}")
    return 0 if ok == len(records) else 1


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

    OUT.mkdir(parents=True, exist_ok=True)

    if args.dims:
        d1, d2, df = args.dims
        try:
            record = run_dims(d1, d2, df, args.alpha, verbose=args.verbose)
        except ValueError as exc:
            parser.error(str(exc))
        return print_summary([record])

    points = SMOKE if args.smoke else EXTENDED if args.extended else INTERESTING
    mode_name = "smoke" if args.smoke else "extended" if args.extended else "interesting"
    print(f"T3 scan mode: {mode_name}; {len(points)} model(s).")
    records = [run_class(c, a, verbose=args.verbose) for c, a in points]
    return print_summary(records)


if __name__ == "__main__":
    raise SystemExit(main())
