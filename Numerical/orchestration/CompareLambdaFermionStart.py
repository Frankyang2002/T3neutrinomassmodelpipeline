"""Compare T3-A-m2 predictions with/without UV running at two lambda benchmarks.

Run from repository root:
    python -m Numerical.orchestration.CompareLambdaFermionStart

Inputs: saved individually fitted 5e-9 A(-2) configuration and the exact prepared
5e-10 A(-2) seed. No fitting, Sobol sampling, or symbolic regeneration.

'UV off' transports the given UV boundary state unchanged to mu_F, then performs
fermion matching, full EFT1 running, scalar matching and final Weinberg running.
This is an RGE ablation with couplings defined at mu_F, not a refit at mu_F.
"""
from __future__ import annotations

import argparse
import csv
import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration.studies.fitted_rge.RunFittedStagedRGComparison import (
    frozen_result,
    with_uv_eft_states,
)
from Numerical.orchestration.studies.fitted_rge.RunFittedWeinbergRunningComparison import predictions
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
from Numerical.running.T3Trajectory import (
    run_renormalisable_t3_trajectory,
    continue_with_weinberg_running,
)
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator

MODEL = "T3_dS1_1_dS2_3_dF_2_alpha_m2"
ROOT = base.ROOT
DEFAULT_OUTPUT = (ROOT / "output/studies/lambda_5e-10_refit/diagnostics/"
                  "a_m2_fermion_start_comparison")
INPUTS = {
    "5e-9_fitted": ROOT / "output/studies/lambda_5e-9_refit/fitted" / (MODEL + ".json"),
    "5e-10_exact_seed": ROOT / "output/studies/lambda_5e-10_refit/prepared" / (MODEL + ".json"),
}
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12", "sin2_theta13", "sin2_theta23")


def evaluate(config, uv_beta, eft_beta, c5_evaluator, target):
    scales = config["scales"]
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(scales["mu_low_gev"])
    numerical = config.get("numerical", {})
    rtol = float(numerical.get("rtol", 1e-8))
    atol = float(numerical.get("atol", 1e-11))
    vev = float(numerical.get("vev_gev", 246.22))

    payload = deepcopy(config)
    payload.setdefault("scan", {"bindings": {}})
    initial = build_uv_state_from_config(SimpleNamespace(raw=payload), {})
    mu_uv = float(initial.mu_gev)
    if not mu_uv > mu_f > mu_s > mu_low > 0:
        raise ValueError("Expected mu_UV > mu_F > mu_S > mu_low > 0")

    # 111: full running of all three RGE regions.
    full = run_renormalisable_t3_trajectory(
        initial, uv_beta, eft_beta,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        rtol=rtol, atol=atol,
    )
    full_c5 = np.asarray(c5_evaluator({}, full), dtype=complex)
    full_low = continue_with_weinberg_running(
        full, full_c5, mu_low, rtol=rtol, atol=atol,
    )
    common_ye = full_low.weinberg.ye
    with_uv = predictions(full_low.weinberg.K, common_ye, target, vev)

    # 011: UV parameters held fixed between mu_UV and mu_F; matching unchanged.
    no_uv = frozen_result(initial, mu_f, uv=True)
    eft_start = build_intermediate_scalar_boundary(
        no_uv.final_state, matching_scale_gev=mu_f,
    )
    eft_run = run_intermediate_scalar_segment(
        eft_start, eft_beta, mu_final_gev=mu_s, rtol=rtol, atol=atol,
    )
    renorm_no_uv = with_uv_eft_states(no_uv, eft_run, mu_uv, mu_f, mu_s)
    no_uv_c5 = np.asarray(c5_evaluator({}, renorm_no_uv), dtype=complex)
    no_uv_low = continue_with_weinberg_running(
        renorm_no_uv, no_uv_c5, mu_low, rtol=rtol, atol=atol,
    )
    # Compare in the same evolved low-energy charged-lepton basis, as in the
    # existing cumulative ablation. This isolates changes to the neutrino C5.
    without_uv = predictions(no_uv_low.weinberg.K, common_ye, target, vev)
    return {"111_full_uv": with_uv, "011_start_at_fermion": without_uv}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--artifact-map", type=Path,
                        help="Optional exact per-model symbolic artifact mapping")
    args = parser.parse_args()
    path_by_model = INPUTS
    mapping = base.load(args.artifact_map) if args.artifact_map else None
    artifact_paths = base.locate_artifacts(MODEL, mapping)
    uv_beta = load_rgbeta_payload(artifact_paths[0])
    eft_beta = load_rgbeta_payload(artifact_paths[1])
    c5 = FinalC5TrajectoryEvaluator(artifact_paths[2])
    target = OscillationFitTarget.from_json(base.TARGET)

    rows = []
    for label, path in path_by_model.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing source configuration: {path}")
        cfg = base.load(path)
        if cfg["representation"] != next(r["representation"] for r in base.original_rows()
                                          if r["model"] == MODEL):
            raise ValueError(f"Representation mismatch in {path}")
        result = evaluate(cfg, uv_beta, eft_beta, c5, target)
        lambda_t3 = float(cfg["base_state"]["ordinary"]["lambdaT3"]["real"])
        row = {"benchmark": label, "source": str(path), "lambdaT3": lambda_t3,
               "model": MODEL, "results": result}
        rows.append(row)
        print(f"\n{label}: lambdaT3={lambda_t3:.9g}")
        for state, prediction in result.items():
            print(f"  {state}: chi2={prediction['chi2']:.8g}")
            for obs in OBS:
                print(f"    {obs}: {prediction['prediction'][obs]:.10g}")
        print("  UV-on / UV-off observable ratios:")
        for obs in OBS:
            off = result["011_start_at_fermion"]["prediction"][obs]
            on = result["111_full_uv"]["prediction"][obs]
            print(f"    {obs}: {on/off:.8g}" if off else f"    {obs}: undefined")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "comparison.json").write_text(
        json.dumps({"description": "Same-source two-scale UV-off comparison; UV-off retains EFT1 and Weinberg running.",
                    "basis": "Shared full-evolution low-scale charged-lepton matrix per lambda point",
                    "data": rows}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["benchmark", "lambdaT3", "model", "stage", "chi2", *OBS]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for stage, result in row["results"].items():
                writer.writerow({"benchmark": row["benchmark"], "lambdaT3": row["lambdaT3"],
                                 "model": MODEL, "stage": stage, "chi2": result["chi2"],
                                 **result["prediction"]})
    print(f"\nSaved: {output / 'comparison.json'}")
    print(f"       {output / 'comparison.csv'}")


if __name__ == "__main__":
    main()
