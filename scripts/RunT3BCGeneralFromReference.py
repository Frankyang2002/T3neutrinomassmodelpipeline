"""Exploratory B(-1)/C(-1) general-scalar fitting from a tracked fitted reference.

This runner uses the existing B/C general-quartic + mixed-bilinear RGEs,
existing one-scale Matchete C5 evaluator, Sobol sampler, and bounded least
squares. It does NOT require historical 1e-7/5e-9 fitted study outputs or
Stage23/26 evidence JSON; the existing conditional sign map is used as an
UNCERTIFIED assumption. No sequential EFT Wilson matching is claimed.

Examples (from repository root):
  python -X utf8 -m scripts.RunT3BCGeneralFromReference --model B_-1 --run-id b_m1_pilot --export-betas --sobol-points 128
  python -X utf8 -m scripts.RunT3BCGeneralFromReference --model C_-1 --run-id c_m1_pilot --export-betas --sobol-points 128
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from scripts import RunT3GeneralScalarBCExploratoryFit as bc
from scripts import RunT3GeneralScalarBCPilot as pilot
from scripts import RunT3ZeroExtraBoundaryBenchmark as base

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "Numerical" / "benchmarks" / "fitted_reference"
OUT = ROOT / "output" / "studies" / "general_scalar_extended" / "fresh_reference_bc"
MODELS = ("B_-1", "C_-1")


def _model_key(model: str) -> str:
    rep = {"B_-1": (2, 2, 1, -1), "C_-1": (2, 2, 3, -1)}[model]
    d1, d2, df, alpha = rep
    return f"T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_m{abs(alpha)}"


def _read(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required input: {path}")
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _save(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build_reference_source(model: str, *, lambda_t3: float | None = None) -> dict:
    """Use a tracked reference config, optionally rescaling Yukawas as a seed.

    Rescaling sqrt(lambda_old/lambda_new) is a starting-guess heuristic only.
    It is not a fitted or general scalar result.
    """
    from scripts import CoupledUVRunning as s31
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry

    if model not in MODELS:
        raise ValueError(f"Only {MODELS} currently supported")
    path = REF / (_model_key(model) + ".json")
    cfg = deepcopy(_read(path))
    registry = load_registry(base.REGISTRY)
    entry, new, reused = s31.inventories(registry, model)
    expected_rep = {
        "d_s1": int(entry["dimensions"][0]),
        "d_s2": int(entry["dimensions"][1]),
        "d_f": int(entry["dimensions"][2]),
        "alpha": int(entry["alpha"]),
    }
    expected_rep["shared_scalar"] = False
    if cfg.get("representation") != expected_rep:
        raise ValueError(f"Representation mismatch: {path}")
    if float(cfg["base_state"]["mu_gev"]) != 1e7:
        raise ValueError("Expected UV boundary at 1e7 GeV")
    if cfg.get("scales") != {"mu_fermion_threshold_gev": 1e5,
                              "mu_scalar_threshold_gev": 1e3,
                              "mu_low_gev": 1e2}:
        raise ValueError("The B/C exploratory evaluator currently fixes these thresholds")
    ordinary = cfg["base_state"]["ordinary"]
    old = ordinary["lambdaT3"]
    old_lam = complex(float(old["real"]), float(old["imag"]))
    if old_lam.imag != 0 or old_lam.real <= 0:
        raise ValueError("Reference scalar lambdaT3 must be positive real")
    new_lam = old_lam.real if lambda_t3 is None else float(lambda_t3)
    if not math.isfinite(new_lam) or new_lam <= 0:
        raise ValueError("--lambda-t3 must be positive and finite")
    if new_lam != old_lam.real:
        factor = math.sqrt(old_lam.real / new_lam)
        for name in ("y1", "y2"):
            spec = ordinary[name]
            spec["real"] = (np.asarray(spec["real"], dtype=float) * factor).tolist()
            spec["imag"] = (np.asarray(spec["imag"], dtype=float) * factor).tolist()
        ordinary["lambdaT3"] = {"real": new_lam, "imag": 0.0}
    for key in reused:
        value = ordinary.get(key)
        if value is not None and abs(base.complex_of(value)) > 1e-14:
            raise ValueError(f"Reference reused quartic {key} is already nonzero")
    return {
        "cfg": cfg, "cfg_path": path, "entry": entry, "new": new,
        "reused": reused, "lambdaT3": complex(new_lam),
        "uv_on": True, "registry": registry,
        "reference_lambda_t3": old_lam.real,
        "seed_rescaling": math.sqrt(old_lam.real / new_lam),
        "provenance": "tracked_fitted_reference; not a newly fitted general-scalar input",
    }


def ensure_betas(model: str, source: dict, *, export: bool, wolframscript: str) -> dict[str, Path]:
    """Use existing live RGBeta export helper without Stage33 historical files."""
    from scripts.RunGeneralT3QuadraticEFT1 import export_live
    from Numerical.core.State import T3Representation

    rep = source["cfg"]["representation"]
    obj = T3Representation(**rep).validated()
    paths = {stage: base.beta_path(model, stage) for stage in ("UV", "EFT1")}
    for stage, path in paths.items():
        if export or not path.is_file():
            if not export:
                raise FileNotFoundError(
                    f"Missing general beta export: {path}. Re-run with --export-betas "
                    "(requires installed Wolfram/RGBeta)."
                )
            print(f"Exporting {model} {stage} general-quadratic RGBeta: {path}", flush=True)
            export_live(model, obj, path, stage=stage, exe=wolframscript)
    return paths


def run(args: argparse.Namespace) -> Path:
    from Numerical.orchestration import RunSmallLambdaRefit
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.fitting.SobolBenchmarkSearch import (
        generate_sobol_parameter_points, refine_sobol_candidates,
    )
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE

    if not (args.sobol_points >= 2 and args.sobol_points & (args.sobol_points-1) == 0):
        raise ValueError("--sobol-points must be a power of two >= 2")
    if args.seed_count < 1 or args.seed_count > 8:
        raise ValueError("--seed-count must be 1..8")
    if not (0 < args.quartic_fraction <= 1 and 0 <= args.mu2_fraction <= 0.2):
        raise ValueError("Require 0 < --quartic-fraction <= 1 and 0 <= --mu2-fraction <= 0.2")
    run_id = pilot.parse_run_id(args.run_id)
    dest = OUT / run_id / model_tag(args.model)
    if dest.exists():
        raise FileExistsError(f"Refusing to replace existing run: {dest}; choose another --run-id")
    source = build_reference_source(args.model, lambda_t3=args.lambda_t3)
    betas = ensure_betas(args.model, source, export=args.export_betas,
                         wolframscript=args.wolframscript)
    coefficient = DEFAULT_SOURCE / (
        {"B_-1": "B_m1_final", "C_-1": "C_m1_final"}[args.model]
        + "_one_loop_probe_weinberg_ordered_coefficient.inputform.txt"
    )
    compiled = bc.compile_source(coefficient)
    bundle = bc.beta_bundle(args.model, source, source["registry"])
    target = OscillationFitTarget.from_json(RunSmallLambdaRefit.TARGET)
    evaluator = bc.CandidateEvaluator(
        "reference_scalar_general", args.model, source, source["registry"],
        bundle, compiled, target, quartic_fraction=args.quartic_fraction,
        mu2_fraction=args.mu2_fraction,
    )
    seed = bc.seed_vector(source)
    if max(abs(seed)) >= args.bound:
        raise ValueError("Rescaled Yukawa seed exceeds --bound; increase bound")
    defs, local = bc.local_parameters(seed, args.bound, args.sobol_fraction,
                                      args.sobol_floor)
    point = dict(zip(bc.NAMES, map(float, seed)))
    candidates, failures = [], []
    seed_chi2 = None
    try:
        first = evaluator(point)
        seed_chi2 = float(first.chi2)
        candidates.append(bc.result_record(first, point, "reference_seed"))
        print(f"Reference general-scalar seed chi2={first.chi2:.8g}", flush=True)
    except Exception as exc:
        failures.append({"source": "reference_seed", "message": f"{type(exc).__name__}: {exc}"})
        print(f"Reference general-scalar seed failed: {exc}", flush=True)
    for i, trial in enumerate(generate_sobol_parameter_points(local, args.sobol_points,
                                                               seed=args.seed), 1):
        try:
            result = evaluator(trial)
            candidates.append(bc.result_record(result, trial, f"sobol_{i}"))
            candidates.sort(key=lambda r: r["chi2"])
            candidates = candidates[:max(args.seed_count, 8)]
            if candidates[0]["source"] == f"sobol_{i}":
                print(f"Sobol {i}/{args.sobol_points}: best chi2={result.chi2:.8g}", flush=True)
        except Exception as exc:
            if len(failures) < 30:
                failures.append({"source": f"sobol_{i}", "message": f"{type(exc).__name__}: {exc}"})
    if not candidates:
        raise RuntimeError("All general-scalar numerical evaluations failed. "
                           f"First failures: {failures[:3]}")
    best = min(candidates, key=lambda r: r["chi2"])
    history = {"enabled": False, "reason": "--no-local-refinement"}
    if not args.no_local_refinement and best["chi2"] > args.target_chi2:
        try:
            better, history = refine_sobol_candidates(
                evaluator=evaluator, observable_names=tuple(target.observable_names),
                definitions=defs, candidates=candidates,
                target_chi2=args.target_chi2, seed_count=args.seed_count,
                active_counts=(5, 8, 12),
                max_nfev=(20, 30, args.max_nfev), sensitivity_step=0.02,
            )
            if better["chi2"] < best["chi2"]:
                best = better
        except Exception as exc:
            history = {"enabled": True, "status": "FailedRetainedBestSobol",
                       "error": f"{type(exc).__name__}: {exc}"}
    elif not args.no_local_refinement:
        history = {"enabled": True, "skipped": True,
                   "reason": "Best Sobol/reference point already passes target chi2"}
    final, details = evaluator.evaluate(best["parameters"], detail=True)
    fitted = bc.trial_record(source, np.array([best["parameters"][k] for k in bc.NAMES]))["cfg"]
    fitted["benchmark_search"] = {"enabled": False}
    fitted["sensitivity"] = {"enabled": False}
    fitted["exploratory_general_scalar_fit"] = {
        "model": args.model, "run_id": run_id, "chi2": float(final.chi2),
        "quartic_fraction": args.quartic_fraction, "mu2_fraction": args.mu2_fraction,
        "source": "tracked_fitted_reference", "physical_matching_certified": False,
        "sequential_eft_matching_performed": False,
    }
    inputs = [source["cfg_path"], base.REGISTRY, coefficient, *betas.values(),
              RunSmallLambdaRefit.TARGET]
    payload = {
        "model": args.model, "run_id": run_id,
        "status": "ExploratoryGeneralScalarOneScaleConditionalFit",
        "chi2_seed": seed_chi2,
        "chi2_final": float(final.chi2),
        "observable_names": list(final.observable_names),
        "prediction": dict(zip(final.observable_names, map(float, final.prediction))),
        "whitened_residual": np.asarray(final.whitened_residual, dtype=float).tolist(),
        "input_lambda_t3": float(source["lambdaT3"].real),
        "reference_lambda_t3": float(source["reference_lambda_t3"]),
        "reference_yukawa_seed_rescale_factor": float(source["seed_rescaling"]),
        "extra_quartic_fraction": args.quartic_fraction,
        "extra_mixed_scalar_bilinear_fraction": args.mu2_fraction,
        "best_parameters": best["parameters"], "best_source": best["source"],
        "sobol_attempts": args.sobol_points, "failures_recorded": failures,
        "local_refinement": history, "numerical_evaluator_calls": evaluator.calls,
        "details": details,
        "input_sha256": {path.relative_to(ROOT).as_posix(): _sha(path) for path in inputs},
        "physical_matchete_rgbeta_coupling_map_certified": False,
        "conditional_lambda_m_minus_lambda_rgbeta_assumed": True,
        "full_sequential_general_scalar_matching_performed": False,
        "general_extra_yukawas_included": False,
        "note": "B(-1)/C(-1) have no extra conjugate Yukawa in the 16-model catalogue. "
                "The saved B/C Matchete one-scale C5 and candidate general-scalar RGEs "
                "are evaluated; no physical full-sequential matching is claimed.",
    }
    dest.mkdir(parents=True, exist_ok=False)
    _save(dest / "fitted_config.json", fitted)
    _save(dest / "result.json", payload)
    with (dest / "observables.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["observable", "prediction", "target"])
        for k, value, target_value in zip(final.observable_names, final.prediction, final.target):
            writer.writerow([k, float(value), float(target_value)])
    print(f"General-scalar fit complete: chi2={float(final.chi2):.8g}", flush=True)
    print(f"Output: {dest}", flush=True)
    return dest


def model_tag(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--lambda-t3", type=float, default=None,
                        help="Fixed scalar coupling; default uses tracked reference (0.1)")
    parser.add_argument("--quartic-fraction", type=float, default=0.1)
    parser.add_argument("--mu2-fraction", type=float, default=0.01)
    parser.add_argument("--export-betas", action="store_true",
                        help="Regenerate general UV and EFT1 betas with Wolfram/RGBeta")
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--sobol-points", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--bound", type=float, default=3.0)
    parser.add_argument("--sobol-fraction", type=float, default=0.5)
    parser.add_argument("--sobol-floor", type=float, default=1e-6)
    parser.add_argument("--seed-count", type=int, default=3)
    parser.add_argument("--target-chi2", type=float, default=1e-8)
    parser.add_argument("--max-nfev", type=int, default=60)
    parser.add_argument("--no-local-refinement", action="store_true")
    args = parser.parse_args(argv)
    if not (args.bound > 0 and args.sobol_floor > 0 and
            args.sobol_fraction > 0 and args.max_nfev > 0 and args.target_chi2 > 0):
        parser.error("Invalid optimisation bounds/budgets")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
