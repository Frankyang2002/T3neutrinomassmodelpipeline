"""Run the improved exploratory UV/EFT1 benchmark for 14 ordinary T3 models.

Uses the existing GeneralRunningBenchmark orchestration and cached canonical
symbolic RGE exports. Does NOT perform C5 matching or low-energy evolution.
Run from the repository root with ``python -X utf8 -m scripts.RunGeneralT3ImprovedRunning``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from Numerical.orchestration.GeneralRunningBenchmark import (
    load_canonical_general_exports,
    run,
)

MODELS = (
    "A_-4", "A_-2", "A_+0",
    "B_-3", "B_+1", "C_-3", "C_+1",
    "D_-2", "D_+0", "D_+2",
    "E_-4", "E_-2", "E_+0", "E_+2",
)


def model_tag(model: str) -> str:
    family, signed_alpha = model.split("_", 1)
    alpha = int(signed_alpha)
    return f"{family}_{'m' if alpha < 0 else 'p'}{abs(alpha)}"


def complex_matrix_diagonal(entries: list[float]) -> list[list[float]]:
    return [[entries[i] if i == j else 0.0 for j in range(3)] for i in range(3)]


def improved_overrides(model: str) -> dict[str, Any]:
    overrides = {
        "g3": 0.75,
        "yu": complex_matrix_diagonal([1e-5, 0.003, 0.5]),
        "y1": complex_matrix_diagonal([0.0004] * 3),
        "y2": complex_matrix_diagonal([0.0004] * 3),
    }
    if model in ("A_-2", "D_-2", "E_-2"):
        overrides["y2Conj"] = complex_matrix_diagonal([0.0004] * 3)
    if model in ("A_+0", "D_+0", "E_+0"):
        overrides["y1Conj"] = complex_matrix_diagonal([0.0004] * 3)
    return overrides


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-run-id", default="full16_general_v1")
    parser.add_argument("--output", type=Path,
                        default=Path("output/general_numerical_all14_improved/running"))
    parser.add_argument("--mu-uv", type=float, default=1e7)
    parser.add_argument("--mu-f", type=float, default=1e5)
    parser.add_argument("--mu-s", type=float, default=1e3)
    parser.add_argument("--model", choices=MODELS, help="Run one model only")
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    statuses = []
    for model in ((args.model,) if args.model else MODELS):
        print(f"{model}: starting improved UV/EFT1 running", flush=True)
        try:
            uv, eft, provenance = load_canonical_general_exports(
                model, args.canonical_run_id)
            result = run(uv, eft, mu_uv=args.mu_uv, mu_f=args.mu_f,
                         mu_s=args.mu_s, overrides=improved_overrides(model))
            result["symbolic_provenance"] = provenance
            result["benchmark"] = "improved_exploratory_14"
            result["uv_benchmark_overrides"] = improved_overrides(model)
            result["physical_mass_parameters_available"] = False
            result["basis_mapping_required_for_physical_interpretation"] = True
            destination = args.output / f"{model_tag(model)}.json"
            destination.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                   encoding="utf-8")
            status = {"model": model, "status": "Success",
                      "running_json": str(destination),
                      "uv_solver": result.get("uv_solver"),
                      "eft1_solver": result.get("eft1_solver")}
            print(f"{model}: UV/EFT1 completed -> {destination}", flush=True)
        except Exception as exc:
            status = {"model": model, "status": "Failed",
                      "error_type": type(exc).__name__, "error": str(exc)}
            print(f"{model}: FAILED: {type(exc).__name__}: {exc}",
                  file=sys.stderr, flush=True)
        statuses.append(status)
        (args.output / "summary.json").write_text(
            json.dumps({"canonical_run_id": args.canonical_run_id,
                        "benchmark": "improved_exploratory_14",
                        "scope": "UV/EFT1 renormalizable running only; no C5 or MZ evolution",
                        "results": statuses}, indent=2, allow_nan=False) + "\n",
            encoding="utf-8")
    successful = sum(s["status"] == "Success" for s in statuses)
    print(f"UV/EFT1 running: {successful}/{len(statuses)} successful")
    return 0 if successful == len(statuses) else 1


if __name__ == "__main__":
    raise SystemExit(main())
