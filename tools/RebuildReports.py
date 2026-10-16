"""Rebuild selected presentation reports from saved T3 numerical results.

Run from the repository root:
    python -m tools.RebuildReports
    python -m tools.RebuildReports --study staged-rg --study weinberg
    python -m tools.RebuildReports --study comparison --study correlations

Does not run RGEs, matching, fitting, or change saved numerical scan files.
Missing optional saved-study inputs are reported and skipped in the default run.
Explicitly requested missing inputs cause a nonzero exit status.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
STUDIES = ("comparison", "staged-rg", "weinberg", "correlations")


def _load_rows(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, list) or not all(isinstance(item, dict) for item in data):
        raise ValueError(f"Expected a list of result records in {path}")
    return data


def _run_study(study: str, *, scenario: str, output_root: Path,
               data_root: Path, correlation_source: Path) -> None:
    if study == "comparison":
        from Numerical.plotting.BuildDisplayResults import build_display
        build_display(data_root, PROJECT_ROOT / "Reports" / "output" / "full",
                      output_root / "comparison_display")
        return

    if study == "correlations":
        from Numerical.plotting.BuildInteractiveCorrelations import build
        target = PROJECT_ROOT / "configs" / "nufit_6_1_ic24_no.json"
        build(correlation_source, output_root / "interactive_correlations.html",
              target if target.is_file() else None)
        return

    source_name = {
        "staged-rg": "staged_rg_comparison",
        "weinberg": "weinberg_running_comparison",
    }[study]
    source = PROJECT_ROOT / "Reports" / "output" / "full" / source_name / scenario / "scan.json"
    rows = _load_rows(source)
    dest = output_root / "rebuilt_studies" / source_name / scenario
    dest.mkdir(parents=True, exist_ok=True)
    if study == "staged-rg":
        from Numerical.orchestration.RunFittedStagedRGComparison import write_plots
    else:
        from Numerical.orchestration.RunFittedWeinbergRunningComparison import write_plots
    write_plots(rows, dest)
    print(f"Rebuilt {study} plots from {source} into {dest}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=STUDIES, action="append",
                        help="Select one or more studies; default: all available studies")
    parser.add_argument("--scenario", default="smallY_smallL")
    parser.add_argument("--data-root", type=Path, default=PROJECT_ROOT / "output" / "full")
    parser.add_argument("--correlation-source", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output" / "full" / "fitted_sensitivity")
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output" / "display_stage5")
    args = parser.parse_args(argv)
    selected = tuple(dict.fromkeys(args.study or STUDIES))
    out = args.output.resolve()
    if out == PROJECT_ROOT or out == args.data_root.resolve():
        parser.error("--output cannot overwrite the project or numerical data root")
    failed = []
    for study in selected:
        try:
            _run_study(study, scenario=args.scenario, output_root=out,
                       data_root=args.data_root.resolve(),
                       correlation_source=args.correlation_source.resolve())
        except (FileNotFoundError, RuntimeError, ValueError) as exc:
            print(f"{study}: unavailable ({exc})")
            if args.study:
                failed.append(study)
        except Exception as exc:
            print(f"{study}: FAILED ({type(exc).__name__}: {exc})")
            failed.append(study)
    if failed:
        print("Failed studies:", ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
