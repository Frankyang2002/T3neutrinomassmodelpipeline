"""Run saved-artifact T3 UV sensitivity scans separately from benchmark diagnostics.

Default: inventory only. --run executes selected cases, --limit caps attempts.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

from Numerical.fitting.BenchmarkSensitivity import run_t3_sensitivity_scan
from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, CONFIG_ROOT, REQUIRED,
    discover_cases, config_for_case, inspect_case,
)

SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")


def _output_path(case: Path) -> Path:
    return case / "data" / "sensitivity_diagnostics.json"


def _target(ordering: str) -> Path:
    return PROJECT_ROOT / "configs" / f"nufit_6_1_ic24_{ordering.lower()}.json"


def run(*, execute: bool, limit: int, force: bool, scenario: str | None,
        model: str | None, gamma_min: float, gamma_max: float,
        points: int) -> dict:
    rows = []
    attempted = 0
    for sc, key, case in discover_cases(OUTPUT_ROOT):
        if scenario and sc != scenario:
            continue
        if model and key != model:
            continue
        inspected = inspect_case(sc, key, case, CONFIG_ROOT)
        destination = _output_path(case)
        row = {"scenario": sc, "model": key, "output": str(destination),
               "status": "", "message": ""}
        try:
            if inspected["missing_artifacts"]:
                row["status"] = "MissingArtifacts"
                row["message"] = ", ".join(inspected["missing_artifacts"])
            elif inspected["config_error"]:
                row["status"] = "InvalidConfig"
                row["message"] = inspected["config_error"]
            elif destination.is_file() and not force and json.loads(
                    destination.read_text(encoding="utf-8")).get("status") == "Success":
                row["status"] = "AlreadyComplete"
            elif not execute:
                row["status"] = "Ready"
            elif attempted >= limit:
                row["status"] = "Deferred"
            else:
                attempted += 1
                config_path = config_for_case(sc, key, CONFIG_ROOT)
                # This is a separate in-memory copy. Original configs are untouched.
                payload = json.loads(config_path.read_text(encoding="utf-8-sig"))
                ordering = str(payload.get("ordering", "NO")).upper()
                if ordering not in {"NO", "IO"}:
                    raise ValueError(f"Explicit NO/IO ordering required, got {ordering}")
                if payload.get("benchmark_search", {}).get("enabled", False):
                    raise ValueError("Refusing benchmark-search-enabled config")
                target = _target(ordering)
                if not target.is_file():
                    raise FileNotFoundError(target)
                payload["sensitivity"] = {
                    "enabled": True,
                    "gamma_min": gamma_min,
                    "gamma_max": gamma_max,
                    "n_points": points,
                }
                data = case / "data"
                result = run_t3_sensitivity_scan(
                    base_payload=payload,
                    uv_rgbeta_path=data / REQUIRED[0],
                    intermediate_rgbeta_path=data / REQUIRED[1],
                    final_weinberg_path=data / REQUIRED[2],
                    oscillation_target_path=target,
                )
                if result.get("status") != "Success":
                    raise ValueError(f"Sensitivity returned {result.get('status')!r}")
                result.update({
                    "scenario": sc, "model_key": key, "ordering": ordering,
                    "config_source": str(config_path),
                    "note": "Separate saved-artifact sensitivity scan; original running diagnostics unchanged.",
                })
                destination.parent.mkdir(parents=True, exist_ok=True)
                temp = destination.with_suffix(".json.tmp")
                temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                encoding="utf-8")
                temp.replace(destination)
                row["status"] = "Success"
        except Exception as exc:
            row["status"] = "Failed"
            row["message"] = f"{type(exc).__name__}: {exc}"
            print(traceback.format_exc(), flush=True)
        rows.append(row)
        print(f"{row['status']:18} {sc}/{key}", flush=True)
    summary = {"attempted": attempted, "cases": rows}
    report = PROJECT_ROOT / "Reports" / "output" / "full" / "observable_sensitivity"
    report.mkdir(parents=True, exist_ok=True)
    (report / "sensitivity_run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="Run sensitivity scans (default: inventory)")
    parser.add_argument("--force", action="store_true", help="Replace previously successful sensitivity outputs")
    parser.add_argument("--limit", type=int, default=1, help="Maximum number of attempted cases (default 1)")
    parser.add_argument("--scenario", choices=SCENARIOS)
    parser.add_argument("--model", help="Exact T3_dS1_... model key")
    parser.add_argument("--gamma-min", type=float, default=0.95)
    parser.add_argument("--gamma-max", type=float, default=1.05)
    parser.add_argument("--points", type=int, default=11)
    args = parser.parse_args(argv)
    if args.limit < 1 or args.points < 3 or not (0 < args.gamma_min < args.gamma_max):
        parser.error("Require --limit >= 1, --points >= 3, 0 < gamma-min < gamma-max")
    summary = run(execute=args.run, limit=args.limit, force=args.force,
                  scenario=args.scenario, model=args.model,
                  gamma_min=args.gamma_min, gamma_max=args.gamma_max, points=args.points)
    return 1 if any(r["status"] == "Failed" for r in summary["cases"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
