"""Fermion-first 16-symbolic / 64-numerical comparison runner.

Standalone entry point: python -m studies.FullT3ReuseStudy
Does not alter the existing pipeline.py --full behaviour.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from common.RunRecords import RunRecord
from studies.FullT3Study import (
    MODELS, COMPARISON_SCENARIOS, DEFAULT_TEMPLATE, PROJECT_ROOT,
    OUTPUT_ROOT, REPORT_ROOT, CONFIG_SET_ROOT, _model_key,
    _load_json, _write_json, _fixed_comparison_config, _dashboard,
)


def _symbolic_record(model_class: str, ds1: int, ds2: int, df: int, alpha: int,
                     *, rebuild: bool) -> RunRecord:
    key = _model_key(ds1, ds2, df, alpha)
    study = f"full/symbolic/{key}"
    study_dir = PROJECT_ROOT / "output" / study
    aggregate = study_dir / "t3_model_comparison.json"
    if rebuild or not aggregate.is_file():
        command = [
            sys.executable, str(PROJECT_ROOT / "pipeline.py"),
            "--dims", str(ds1), str(ds2), str(df), "--alpha", str(alpha),
            "--study", study, "--threshold", "F", "--threshold", "S1", "S2",
        ]
        print("\n> " + " ".join(command), flush=True)
        result = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Symbolic pipeline failed for {key}: exit {result.returncode}")
    if not aggregate.is_file():
        raise FileNotFoundError(f"Missing symbolic summary: {aggregate}")
    summaries = json.loads(aggregate.read_text(encoding="utf-8"))
    if len(summaries) != 1:
        raise ValueError(f"Expected one symbolic record for {key}, found {len(summaries)}")
    summary = summaries[0]
    if any(summary.get(k) != "Success" for k in
           ("BuildStatus", "MatchingStatus", "UVRGEStatus")):
        raise RuntimeError(f"Symbolic prerequisites failed for {key}")
    # The pipeline creates a model-specific subdirectory under the study directory.
    candidates = [p for p in study_dir.iterdir() if p.is_dir() and (p / "data").is_dir()]
    if len(candidates) != 1:
        raise RuntimeError(f"Cannot identify unique symbolic model output in {study_dir}: {candidates}")
    return RunRecord(
        name=candidates[0].name, alpha=alpha, d_s1=ds1, d_s2=ds2,
        d_f=df, return_code=0, summary=summary, output_dir=candidates[0],
        shared_scalar=False,
    )


def _scenario_record(source: RunRecord, scenario: str) -> RunRecord:
    destination = OUTPUT_ROOT / "comparison" / scenario / source.name
    destination.mkdir(parents=True, exist_ok=True)
    # Numerical evaluation requires the symbolic JSON artifacts; do not share
    # a mutable output directory across different numerical benchmarks.
    shutil.copytree(source.output_dir / "data", destination / "data", dirs_exist_ok=True)
    return RunRecord(
        name=source.name, alpha=source.alpha, d_s1=source.d_s1,
        d_s2=source.d_s2, d_f=source.d_f, return_code=0,
        summary=dict(source.summary), output_dir=destination,
        shared_scalar=False,
    )


def run(*, rebuild: bool = False, template_path: Path = DEFAULT_TEMPLATE) -> int:
    from Numerical.orchestration import PipelineNumericalResults as numerical_module
    from Numerical.orchestration.PipelineNumericalResults import run_pipeline_numerical_results

    # Work around the verified parents[1] bug without changing repository files.
    numerical_module.PROJECT_ROOT = PROJECT_ROOT
    numerical_module.RAW_OUTPUT_ROOT = PROJECT_ROOT / "output"
    numerical_module.REPORT_OUTPUT_ROOT = PROJECT_ROOT / "Reports" / "output"

    template = _load_json(template_path)
    started = time.monotonic()
    results = []
    failed = False
    for model_class, ds1, ds2, df, alpha in MODELS:
        key = _model_key(ds1, ds2, df, alpha)
        try:
            source = _symbolic_record(model_class, ds1, ds2, df, alpha, rebuild=rebuild)
        except Exception as exc:
            print(f"SYMBOLIC FAILED {key}: {exc}", flush=True)
            failed = True
            continue
        for scenario, (yukawa, scalar) in COMPARISON_SCENARIOS.items():
            config = _fixed_comparison_config(
                template, model_class, ds1, ds2, df, alpha, yukawa, scalar,
            )
            config_path = CONFIG_SET_ROOT / "comparison" / scenario / "input" / f"{key}.json"
            _write_json(config_path, config)
            record = _scenario_record(source, scenario)
            print(f"NUMERICAL {key} / {scenario}", flush=True)
            ok = run_pipeline_numerical_results(record, config_path)
            _write_json(record.output_dir / "numerical_run_summary.json", record.summary)
            failed |= not ok
            results.append({
                "model": key, "scenario": scenario,
                "status": "Success" if ok else "Failed",
                "error": record.summary.get("PipelineNumericalResultsError"),
                "lambdaT3_effective": config["comparison_point"]["lambdaT3_effective"],
            })
    for scenario, (yukawa, scalar) in COMPARISON_SCENARIOS.items():
        _dashboard(OUTPUT_ROOT / "comparison" / scenario,
                   REPORT_ROOT / "comparison" / scenario,
                   f"T3 normalized comparison: {scenario} (Y={yukawa:g}, lambda_ref={scalar:g})")
    _write_json(OUTPUT_ROOT / "full_numerical_summary.json", {
        "status": "Failed" if failed else "Success",
        "study_mode": "symbolic_reuse_comparison_only",
        "symbolic_model_count": len(MODELS),
        "numerical_case_count": len(results),
        "runtime_seconds": time.monotonic() - started,
        "results": results,
    })
    return int(failed)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Regenerate symbolic results for all models")
    parser.add_argument("--numerical", type=Path, default=DEFAULT_TEMPLATE)
    args = parser.parse_args()
    path = args.numerical if args.numerical.is_absolute() else PROJECT_ROOT / args.numerical
    return run(rebuild=args.rebuild, template_path=path)


if __name__ == "__main__":
    raise SystemExit(main())
