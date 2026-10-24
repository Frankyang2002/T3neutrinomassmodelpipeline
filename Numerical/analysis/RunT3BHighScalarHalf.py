"""One-point T3-B-m1 scalar-coupling sensitivity test at lambdaT3=0.5.

Uses an existing fitted benchmark, with identical UV Yukawas and masses.
Creates a fresh output tree and runs the *unchanged* fixed-order production
pipeline plus the optional exploratory coupled C12/C5 sidecar.

Install: Numerical/analysis/RunT3BHighScalarHalf.py
Run: python -m Numerical.analysis.RunT3BHighScalarHalf --mode optimised
"""
from __future__ import annotations

import argparse
import copy
import json
import shutil
from pathlib import Path

from Numerical.orchestration.RunSmallLambdaRefit import ROOT, locate_artifacts
from Numerical.orchestration.T3BIntegratedCoupledPipeline import run

MODEL = "T3_dS1_2_dS2_2_dF_1_alpha_m1"
LAMBDA = 0.5
BASE = ROOT / "output/studies/lambda_5e-9_refit/comparison"


def load(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing required input: {path}")
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", default="optimised", choices=("optimised", "shared_a_m2"))
    ap.add_argument("--config", type=Path, default=None,
                    help="Override reference fitted config")
    ap.add_argument("--artifact-map", type=Path, default=None)
    args = ap.parse_args()

    source = (args.config if args.config is not None else
              BASE / args.mode / "models" / MODEL / "shared_config.json")
    config = load(source)
    rep = config["representation"]
    actual = (int(rep["d_s1"]), int(rep["d_s2"]), int(rep["d_f"]),
              int(rep["alpha"]), bool(rep.get("shared_scalar", False)))
    if actual != (2, 2, 1, -1, False):
        raise ValueError(f"Wrong model {actual}")

    previous = config["base_state"]["ordinary"]["lambdaT3"]
    modified = copy.deepcopy(config)
    modified["base_state"]["ordinary"]["lambdaT3"] = {
        "real": LAMBDA, "imag": 0.0
    }
    modified["benchmark_search"] = {"enabled": False}
    modified["sensitivity"] = {"enabled": False}
    modified.pop("comparison_point", None)

    directory = (ROOT / "output/studies/t3b_lambdaT3_0p5" /
                 args.mode / MODEL)
    config_out = directory / "benchmark_lambdaT3_0p5.json"
    report_out = directory / "comparison_summary.json"
    data_dir = directory / "data"
    if directory.exists():
        raise FileExistsError(f"Output already exists; refusing to overwrite: {directory}")

    paths = locate_artifacts(MODEL, load(args.artifact_map) if args.artifact_map else None)
    filenames = ("uv_rgbeta_rge.json", "eft1_rgbeta_rge.json",
                 "final_weinberg_coefficient.json")

    data_dir.mkdir(parents=True, exist_ok=False)
    for src, name in zip(paths, filenames):
        shutil.copy2(src, data_dir / name)
    config_out.write_text(json.dumps(modified, indent=2) + "\n", encoding="utf-8")

    result = run(config=config_out, output_dir=directory,
                 enable_coupled=True, max_order=4, samples=120)

    exploratory = load(data_dir / "t3b_coupled_ll_exploratory.json")
    comparison = {
        "status": result["status"],
        "model": MODEL,
        "mode": args.mode,
        "old_lambdaT3": previous,
        "new_lambdaT3": LAMBDA,
        "fixed_uv_yukawas_and_masses": True,
        "fixed_order_and_coupled_consistent_inputs": True,
        "exploratory_coupled": True,
        "relative_c5_change_scalar_threshold":
            exploratory["coupled_at_ms"]["relative_difference_from_fixed_order"],
        "relative_c5_change_low_energy":
            exploratory["low_energy"]["relative_c5_difference_at_low_scale"],
        "relative_reverse_feedback":
            exploratory["coupled_at_ms"]["relative_reverse_feedback"],
        "low_energy_observable_differences":
            exploratory["low_energy"]["differences"],
        "matching_checks": exploratory["matching_checks"],
        "perturbative_accuracy": exploratory["perturbative_accuracy"],
        "production_unchanged": result["production_default_unchanged"],
    }
    report_out.write_text(json.dumps(comparison, indent=2) + "\n",
                          encoding="utf-8")
    print("Saved:", report_out)
    print("Original lambdaT3:", previous)
    print("New lambdaT3:", LAMBDA)
    print("Relative C5 difference at low energy:",
          comparison["relative_c5_change_low_energy"])
    print("Relative reverse feedback:", comparison["relative_reverse_feedback"])
    print("Low energy observable differences:")
    print(json.dumps(comparison["low_energy_observable_differences"], indent=2))
    print("Exploratory result only; not a validated NLO prediction.")


if __name__ == "__main__":
    main()
