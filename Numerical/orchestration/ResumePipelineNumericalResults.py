"""Resume fixed-comparison numerical outputs without Matchete or RGBeta execution.

Run: python -m Numerical.orchestration.ResumePipelineNumericalResults [--run]
Default mode is a read-only inventory. Run mode invokes only the existing
pipeline numerical-results function for complete, unsolved cases.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

from common.RunRecords import RunRecord
from Numerical.orchestration.PipelineNumericalResults import (
    PIPELINE_NUMERICAL_CONFIG_KIND,
    run_pipeline_numerical_results,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "output" / "full" / "comparison"
CONFIG_ROOT = PROJECT_ROOT / "configs" / "generated_model_sets" / "comparison"
REPORT_ROOT = PROJECT_ROOT / "Reports" / "output" / "full"
SUMMARY_JSON = PROJECT_ROOT / "Reports" / "output" / "full" / "numerical_resume_summary.json"
SUMMARY_CSV = PROJECT_ROOT / "Reports" / "output" / "full" / "numerical_resume_summary.csv"

MODEL_RE = re.compile(r"^T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)$")
REQUIRED = ("uv_rgbeta_rge.json", "eft1_rgbeta_rge.json", "final_weinberg_coefficient.json")


def parse_model_key(key: str) -> tuple[int, int, int, int]:
    match = MODEL_RE.fullmatch(key)
    if match is None:
        raise ValueError(f"Unrecognized T3 model key: {key}")
    d1, d2, df, sign, magnitude = match.groups()
    return int(d1), int(d2), int(df), int(magnitude) * (-1 if sign == "m" else 1)


def discover_cases(output_root: Path = OUTPUT_ROOT) -> list[tuple[str, str, Path]]:
    """Find only expected comparison/model/physical-run directory layouts."""
    cases: list[tuple[str, str, Path]] = []
    if not output_root.is_dir():
        return cases
    for scenario_dir in sorted(p for p in output_root.iterdir() if p.is_dir()):
        for model_dir in sorted(p for p in scenario_dir.iterdir() if p.is_dir()):
            try:
                parse_model_key(model_dir.name)
            except ValueError:
                continue
            for run_dir in sorted(p for p in model_dir.iterdir() if p.is_dir()):
                if run_dir.name.startswith("T3_"):
                    cases.append((scenario_dir.name, model_dir.name, run_dir))
    return cases


def config_for_case(scenario: str, key: str, config_root: Path = CONFIG_ROOT) -> Path:
    """Prefer immutable scenario input; resolved configs may contain stale fits."""
    candidates = (
        config_root / scenario / "input" / f"{key}.json",
        config_root / scenario / "resolved" / f"{key}.json",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def inspect_case(scenario: str, key: str, case: Path, config_root: Path = CONFIG_ROOT) -> dict[str, Any]:
    d1, d2, df, alpha = parse_model_key(key)
    config = config_for_case(scenario, key, config_root)
    missing = [name for name in REQUIRED if not (case / "data" / name).is_file()]
    issue = None
    if not config.is_file():
        issue = f"Numerical configuration not found: {config}"
    else:
        try:
            payload = json.loads(config.read_text(encoding="utf-8-sig"))
            rep = payload["representation"]
            dims = (int(rep["d_s1"]), int(rep["d_s2"]), int(rep["d_f"]), int(rep["alpha"]))
            if payload.get("kind") != PIPELINE_NUMERICAL_CONFIG_KIND or dims != (d1, d2, df, alpha):
                issue = "Numerical configuration kind/representation mismatch"
            if payload.get("benchmark_search", {}).get("enabled", False):
                issue = "Numerical config enables benchmark search; refusing automatic resume"
        except (OSError, ValueError, KeyError, TypeError) as exc:
            issue = f"Unreadable numerical configuration: {exc}"
    diagnostics = case / "data" / "running_diagnostics.json"
    finished = False
    if diagnostics.is_file():
        try:
            finished = json.loads(diagnostics.read_text(encoding="utf-8-sig")).get("status") == "Success"
        except (OSError, ValueError):
            pass
    return {
        "scenario": scenario, "model_key": key, "model": case.name,
        "case_directory": str(case), "config": str(config), "resume_config": "",
        "missing_artifacts": missing, "config_error": issue,
        "existing_diagnostics_success": finished,
        "status": "NotRun", "message": "",
    }


def comparison_resume_config(source: Path, scenario: str, key: str,
                             config_root: Path = CONFIG_ROOT) -> Path:
    """Persist fixed-point resume settings; never modify scenario input/resolved configs.

    Sensitivity scans are outside the four-scenario fixed comparison. They can
    reference an incorrect NuFIT path and dramatically lengthen a resume.
    """
    payload = json.loads(source.read_text(encoding="utf-8-sig"))
    if payload.get("kind") != PIPELINE_NUMERICAL_CONFIG_KIND:
        raise ValueError(f"Wrong numerical config kind: {source}")
    if payload.get("benchmark_search", {}).get("enabled", False):
        raise ValueError("Fixed comparison resume must not trigger benchmark search")
    payload["sensitivity"] = {**payload.get("sensitivity", {}), "enabled": False}
    target = config_root / scenario / "numerical_resume" / f"{key}.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    temp = target.with_suffix(".json.tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temp.replace(target)
    return target


def resume(output_root: Path = OUTPUT_ROOT, config_root: Path = CONFIG_ROOT, *, execute: bool = False,
           force: bool = False, limit: int | None = None, runner=run_pipeline_numerical_results) -> dict[str, Any]:
    """Inventory first; optionally run available fixed points sequentially."""
    cases = discover_cases(output_root)
    results: list[dict[str, Any]] = []
    attempts = 0
    for scenario, key, case in cases:
        row = inspect_case(scenario, key, case, config_root)
        if row["missing_artifacts"]:
            row["status"] = "MissingArtifacts"
        elif row["config_error"]:
            row["status"] = "InvalidConfig"
        elif row["existing_diagnostics_success"] and not force:
            row["status"] = "AlreadyComplete"
        elif not execute:
            row["status"] = "Ready"
        elif limit is not None and attempts >= limit:
            row["status"] = "Deferred"
        else:
            attempts += 1
            d1, d2, df, alpha = parse_model_key(key)
            record = RunRecord(name=case.name.replace("_", "-", 1), alpha=alpha,
                               d_s1=d1, d_s2=d2, d_f=df, return_code=0,
                               summary={}, output_dir=case)
            try:
                effective_config = comparison_resume_config(
                    Path(row["config"]), scenario, key, config_root
                )
                row["resume_config"] = str(effective_config)
                success = runner(record, effective_config)
                row["status"] = "Success" if success else "Failed"
                row["message"] = str(record.summary.get("PipelineNumericalResultsError", ""))
            except Exception as exc:
                row["status"] = "Failed"
                row["message"] = f"{type(exc).__name__}: {exc}"
        results.append(row)
        print(f"{row['status']:18} {scenario}/{key}/{case.name}", flush=True)
    counts = {status: sum(row["status"] == status for row in results)
              for status in sorted({row["status"] for row in results})}
    return {"mode": "run" if execute else "dry-run", "discovered": len(results),
            "attempted": attempts, "counts": counts, "cases": results}


def write_summary(summary: dict[str, Any], json_path: Path = SUMMARY_JSON,
                  csv_path: Path = SUMMARY_CSV) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    fields = ("scenario", "model_key", "model", "status", "message", "case_directory", "config",
              "missing_artifacts", "config_error", "existing_diagnostics_success", "resume_config")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in summary["cases"]:
            writer.writerow({field: ", ".join(row[field]) if field == "missing_artifacts" else row[field]
                             for field in fields})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="Actually run numerical continuation")
    parser.add_argument("--force", action="store_true", help="Redo already completed numerical outputs")
    parser.add_argument("--limit", type=int, default=None, help="Maximum number of newly attempted cases")
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be >= 1")
    summary = resume(execute=args.run, force=args.force, limit=args.limit)
    write_summary(summary)
    print("Summary:", SUMMARY_JSON, "and", SUMMARY_CSV)
    print("Counts:", summary["counts"])
    return 1 if summary["counts"].get("Failed", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
