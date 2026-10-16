"""Rebuild selected T3 study plots from saved scan.json, without numerical reruns.

Project path: tools/RebuildSavedStudyPlots.py
Run from repository root:
    python tools/RebuildSavedStudyPlots.py --study staged-rg
    python tools/RebuildSavedStudyPlots.py --study weinberg

The source scan.json is read only. Default figure outputs are under
Reports/output/display/rebuilt_studies/<study>/<scenario>/.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCENARIO = "smallY_smallL"
STUDIES = {
    "staged-rg": ("staged_rg_comparison", "RunFittedStagedRGComparison"),
    "weinberg": ("weinberg_running_comparison", "RunFittedWeinbergRunningComparison"),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, choices=tuple(STUDIES))
    parser.add_argument("--scenario", default=DEFAULT_SCENARIO)
    parser.add_argument("--source", type=Path, help="Override saved scan.json path")
    parser.add_argument("--output-dir", type=Path, help="Override plot destination")
    args = parser.parse_args(argv)

    folder, module = STUDIES[args.study]
    source = args.source or PROJECT_ROOT / "Reports" / "output" / "full" / folder / args.scenario / "scan.json"
    destination = args.output_dir or PROJECT_ROOT / "Reports" / "output" / "display" / "rebuilt_studies" / args.study / args.scenario
    source = source.resolve()
    destination = destination.resolve()
    if not source.is_file():
        parser.error(f"Saved scan missing: {source}. This tool never recomputes the scan.")
    with source.open(encoding="utf-8-sig") as handle:
        rows = json.load(handle)
    if not isinstance(rows, list):
        parser.error(f"Expected a JSON list of scan rows: {source}")
    usable = [row for row in rows if isinstance(row, dict) and row.get("status") == "Success"]
    if not usable:
        parser.error(f"No successful scan rows: {source}")
    if destination == source.parent:
        parser.error("Destination must differ from the saved scan directory")

    # Import plotting functions only after verifying that saved data exists.
    if module == "RunFittedStagedRGComparison":
        from Numerical.orchestration.RunFittedStagedRGComparison import write_plots
    else:
        from Numerical.orchestration.RunFittedWeinbergRunningComparison import write_plots
    destination.mkdir(parents=True, exist_ok=True)
    write_plots(rows, destination)
    print(f"Input (read only): {source}")
    print(f"Successful rows: {len(usable)}")
    print(f"Plot destination: {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
