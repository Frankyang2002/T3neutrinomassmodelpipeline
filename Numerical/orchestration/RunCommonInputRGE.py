"""Numerical-only common-UV-input T3 RGE study; reuse existing symbolic artifacts.

From the repository root:
    python -m Numerical.orchestration.RunCommonInputRGE
    python -m Numerical.orchestration.RunCommonInputRGE --run --limit 2
    python -m Numerical.orchestration.RunCommonInputRGE --run
    python -m Numerical.orchestration.RunCommonInputRGE --run --force

All eligible smallY_smallL representation configurations receive the UV ordinary
parameters and scales of fitted A-m2. By default, E receives lambdaT3=0.05
rather than 0.10, preserving the existing two-charge-channel benchmark
normalisation. No fitting, Matchete, RGBeta or other symbolic code is called.
No existing fitted configurations, original diagnostics or symbolic artifacts
are modified. This is a controlled comparison, not a new phenomenological fit.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile

from common.RunRecords import RunRecord
from Numerical.orchestration.ResumePipelineNumericalResults import (
    PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases, parse_model_key,
)
from Numerical.orchestration.RunFittedSensitivityComparison import best_candidate
from Numerical.orchestration.RunFittedRGEDisplay import (
    CLASS_BY_DIMS, fixed_payload, short_label,
)
from Numerical.orchestration.PipelineNumericalResults import _run_and_build_payload
from Numerical.plotting.InteractiveModelComparison import write_dashboard
from Numerical.plotting.RunningResultFigures import generate_running_result_figures

SCENARIO = "smallY_smallL"
REFERENCE_MODEL = "T3_dS1_1_dS2_3_dF_2_alpha_m2"
DEST_ROOT = PROJECT_ROOT / "Reports" / "output" / "full" / "common_input_rge"
DISPLAY_ROOT = PROJECT_ROOT / "Reports" / "output" / "display"
COMMON_ORDINARY = ("y1", "y2", "lambdaT3", "MF", "mS1Sq", "mS2Sq")


def save_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def input_signature(payload: dict, artifacts: list[Path]) -> str:
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                      ensure_ascii=False).encode("utf-8"))
    for artifact in artifacts:
        stat = artifact.stat()
        digest.update(str(artifact.resolve()).encode("utf-8"))
        digest.update(f"{stat.st_size}:{stat.st_mtime_ns}".encode("ascii"))
    return digest.hexdigest()


def common_payload(model: str, source: dict, reference: dict, *, class_e_half: bool) -> dict:
    """Keep the target's representation; change only requested UV parameters/scales."""
    payload = deepcopy(source)
    source_ordinary = payload["base_state"]["ordinary"]
    reference_ordinary = reference["base_state"]["ordinary"]
    for name in COMMON_ORDINARY:
        if name not in source_ordinary or name not in reference_ordinary:
            raise KeyError(f"Missing common UV parameter {name!r} for {model}")
        source_ordinary[name] = deepcopy(reference_ordinary[name])
    payload["scales"] = deepcopy(reference["scales"])
    # Preserve target-model representation, model-specific nonselected scalar
    # interactions, and the numerical integration options.
    d1, d2, df, alpha = parse_model_key(model)
    rep = payload["representation"]
    if (int(rep["d_s1"]), int(rep["d_s2"]), int(rep["d_f"]), int(rep["alpha"])) != (d1, d2, df, alpha):
        raise ValueError(f"Representation mismatch for {model}")
    class_name = CLASS_BY_DIMS.get((d1, d2, df))
    if class_name is None:
        raise ValueError(f"Unexpected representation dimensions for {model}")
    if class_name == "E" and class_e_half:
        scalar = source_ordinary["lambdaT3"]
        if isinstance(scalar, dict) and isinstance(scalar.get("real"), (int, float)):
            scalar["real"] *= 0.5
            scalar["imag"] = float(scalar.get("imag", 0.0)) * 0.5
        elif isinstance(scalar, (int, float)):
            source_ordinary["lambdaT3"] = scalar * 0.5
        else:
            raise TypeError("Unexpected lambdaT3 encoding; cannot apply class-E convention")
    payload["benchmark_search"] = {**payload.get("benchmark_search", {}), "enabled": False}
    payload["sensitivity"] = {**payload.get("sensitivity", {}), "enabled": False}
    # Avoid accidentally using old fit/scan bindings when constructing the UV state.
    payload["scan"] = {"bindings": {}}
    return payload


def make_dashboard(records: list[dict], destination: Path) -> None:
    """Adapter for the existing read-only interactive dashboard generator."""
    with tempfile.TemporaryDirectory(prefix="t3_common_input_display_") as tmp:
        study = Path(tmp)
        count = 0
        for row in records:
            if row["status"] not in {"Success", "AlreadyComplete"}:
                continue
            original = Path(row["diagnostics"])
            diagnostic = json.loads(original.read_text(encoding="utf-8-sig"))
            if diagnostic.get("status") != "Success":
                continue
            benchmark = dict(diagnostic.get("benchmark_search") or {})
            representation = dict(benchmark.get("model") or {})
            representation["model_key"] = short_label(row["model"])
            benchmark["model"] = representation
            diagnostic["benchmark_search"] = benchmark
            path = study / row["model"] / "data" / "running_diagnostics.json"
            save_json(path, diagnostic)
            count += 1
        if not count:
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        write_dashboard(study, destination, f"Common-input T3 RGE: {SCENARIO} ({count} models)")
        print(f"Interactive dashboard ({count} models): {destination}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", default=SCENARIO, choices=(SCENARIO,))
    parser.add_argument("--run", action="store_true", help="Run numerical integration; default is inventory only")
    parser.add_argument("--limit", type=int, default=0, help="Maximum eligible models to inspect (0: all)")
    parser.add_argument("--force", action="store_true", help="Recalculate common-input results even when cached")
    parser.add_argument("--no-figures", action="store_true", help="Skip per-model PNG figures")
    parser.add_argument("--raw-class-e-lambda", action="store_true",
                        help="Use raw lambdaT3=0.1 also for E, rather than benchmark-normalised 0.05")
    args = parser.parse_args(argv)
    if args.limit < 0:
        parser.error("--limit must be nonnegative")
    if args.raw_class_e_lambda:
        # Avoid silently mixing two incompatible conventions in the same output folder.
        parser.error("The raw-class-E alternative requires a separate output convention; "
                     "use the default normalised benchmark for this study.")
    selected = best_candidate(args.scenario, REFERENCE_MODEL)
    if selected is None:
        parser.error("Reference A-m2 fitted candidate not found")
    _, _, ref_config, _ = selected
    reference = fixed_payload(Path(ref_config))
    cases = [(model, case) for scenario, model, case in discover_cases(OUTPUT_ROOT)
             if scenario == args.scenario]
    if len({model for model, _ in cases}) != len(cases):
        parser.error("Ambiguous duplicate physical-run directories")
    cases = sorted(cases)[:args.limit or None]
    output = DEST_ROOT / args.scenario
    records: list[dict] = []
    for model, case in cases:
        dest = output / model
        row = {"model": model, "label": short_label(model), "status": "Ready",
               "diagnostics": str(dest / "running_diagnostics.json"), "message": ""}
        try:
            chosen = best_candidate(args.scenario, model)
            if chosen is None:
                raise FileNotFoundError("Saved fitted model config missing")
            _, _, config, _ = chosen
            payload = common_payload(model, fixed_payload(Path(config)), reference, class_e_half=True)
            artifacts = [case / "data" / file for file in REQUIRED]
            missing = [str(p) for p in artifacts if not p.is_file()]
            if missing:
                raise FileNotFoundError("Missing precomputed numerical artifacts: " + ", ".join(missing))
            key = input_signature(payload, artifacts)
            diagnostic_file = Path(row["diagnostics"])
            if diagnostic_file.is_file() and not args.force:
                prev = json.loads(diagnostic_file.read_text(encoding="utf-8-sig"))
                if (prev.get("status") == "Success"
                        and prev.get("common_input_rge", {}).get("input_signature") == key):
                    row["status"] = "AlreadyComplete"
                    records.append(row)
                    print(f"{row['status']:18} {row['label']}", flush=True)
                    continue
            if not args.run:
                row["status"] = "NeedsRGE" if diagnostic_file.is_file() else "Ready"
            else:
                # The runner is numerical only: files in case/data contain already-generated
                # RGBeta RGEs and matching artifacts and are never written to here.
                save_json(dest / "common_input_config.json", payload)
                d1, d2, df, alpha = parse_model_key(model)
                record = RunRecord(name=case.name.replace("_", "-", 1), alpha=alpha,
                                   d_s1=d1, d_s2=d2, d_f=df, return_code=0,
                                   summary={}, output_dir=case)
                result = _run_and_build_payload(record=record, payload=payload,
                                                config_path=dest / "common_input_config.json")
                result["common_input_rge"] = {
                    "input_signature": key, "reference_model": REFERENCE_MODEL,
                    "reference_config": str(ref_config),
                    "target_model_fitted_config": str(config),
                    "lambda_convention": "class E: 0.5 * reference; other classes: reference",
                    "note": "Unrefitted common UV point; reuse existing symbolic artifacts",
                }
                save_json(diagnostic_file, result)
                if result.get("status") != "Success":
                    raise RuntimeError(f"Numerical diagnostics status: {result.get('status')}")
                if not args.no_figures:
                    generate_running_result_figures(diagnostic_file, dest / "figures")
                row["status"] = "Success"
        except Exception as exc:
            row["status"] = "Failed"
            row["message"] = f"{type(exc).__name__}: {exc}"
        records.append(row)
        print(f"{row['status']:18} {row['label']} {row['message']}", flush=True)
    if args.run:
        save_json(output / "common_input_summary.json", {
            "scenario": args.scenario, "reference_model": REFERENCE_MODEL,
            "note": "Fixed common numerical inputs, NOT independently fitted points",
            "results": records,
        })
        make_dashboard(records, DISPLAY_ROOT / f"common_input_rge_{args.scenario}.html")
    else:
        print("Inventory only. Add --run to evaluate existing numerical artifacts.")
    return int(any(row["status"] == "Failed" for row in records))


if __name__ == "__main__":
    raise SystemExit(main())
