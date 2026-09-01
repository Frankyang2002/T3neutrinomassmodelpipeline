from __future__ import annotations

"""
Construct a T3 Yukawa point that reproduces a chosen low-scale neutrino target
after one-loop SM running.

Strategy
--------
1. Start from the SM couplings specified at the numerical matching scale.
2. Run those SM couplings down to the requested low scale with K=0.
3. Build the desired low-scale neutrino target C5.
4. Evolve that target C5 backwards to the matching scale using the low-scale
   SM couplings obtained in step 2.
5. Evaluate the actual Matchete T3 loop factors F_r for the chosen heavy masses.
6. Construct y1 and y2 so the matched T3 C5 equals the required high-scale C5.

This keeps the existing forward numerical pipeline unchanged. The output is
another numerical JSON config with only the T3 Yukawa matrices replaced.
"""

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from NumericalWeinbergRGE import SMInitialConditions, evolve_weinberg
from T3NeutrinoTarget import build_normal_ordering_target
from T3PhysicalYukawaPoint import evaluate_t3_loop_factors
from T3YukawaFit import fit_three_heavy_balanced, reconstruct_c5


def run_sm_to_low_scale(config: dict) -> SMInitialConditions:
    """Return the SM couplings at the numerical low scale."""

    sm = config["sm"]

    zero_k = np.zeros((3, 3), dtype=complex)

    initial = SMInitialConditions(
        gY=sm["gY"],
        g2=sm["g2"],
        g3=sm["g3"],
        lambdaH=sm["lambdaH"],
        ye=np.asarray(sm["ye"], dtype=float),
        yu=np.asarray(sm["yu"], dtype=float),
        yd=np.asarray(sm["yd"], dtype=float),
        K=zero_k,
    )

    result = evolve_weinberg(
        initial,
        config["mu_initial_gev"],
        config["mu_final_gev"],
    )

    return SMInitialConditions(
        gY=result.gY,
        g2=result.g2,
        g3=result.g3,
        lambdaH=result.lambdaH,
        ye=result.ye,
        yu=result.yu,
        yd=result.yd,
        K=zero_k,
    )


def required_high_scale_c5(
    config: dict,
    *,
    lightest_mass_ev: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (target low-scale C5, required high-scale C5)."""

    vev_gev = float(config.get("vev_gev", 246.22))

    target = build_normal_ordering_target(
        lightest_mass_ev=lightest_mass_ev,
        vev_gev=vev_gev,
    )

    low_sm = run_sm_to_low_scale(config)

    backward_initial = SMInitialConditions(
        gY=low_sm.gY,
        g2=low_sm.g2,
        g3=low_sm.g3,
        lambdaH=low_sm.lambdaH,
        ye=low_sm.ye,
        yu=low_sm.yu,
        yd=low_sm.yd,
        K=target.c5_matrix_gev_inv,
    )

    backward = evolve_weinberg(
        backward_initial,
        config["mu_final_gev"],
        config["mu_initial_gev"],
    )

    return target.c5_matrix_gev_inv, backward.K


def build_rge_corrected_config(
    c5_path: Path,
    config: dict,
    *,
    lightest_mass_ev: float = 0.01,
) -> tuple[dict, dict]:
    """Return a config whose T3 Yukawas reproduce the low-scale target."""

    low_target, high_required = required_high_scale_c5(
        config,
        lightest_mass_ev=lightest_mass_ev,
    )

    loop_factors = evaluate_t3_loop_factors(
        c5_path,
        config,
    )

    fit = fit_three_heavy_balanced(
        high_required,
        loop_factors,
    )

    reconstructed_high = reconstruct_c5(
        fit.y1,
        fit.y2,
        loop_factors,
    )

    updated = copy.deepcopy(config)
    model = updated["t3"]

    model["y1_real"] = fit.y1.real.tolist()
    model["y1_imag"] = fit.y1.imag.tolist()
    model["y2_real"] = fit.y2.real.tolist()
    model["y2_imag"] = fit.y2.imag.tolist()

    high_error = float(
        np.linalg.norm(reconstructed_high - high_required)
        / np.linalg.norm(high_required)
    )

    scale_ratio = float(
        np.linalg.norm(high_required) / np.linalg.norm(low_target)
    )

    diagnostics = {
        "Status": "Success",
        "Ordering": "NO",
        "LightestMassEV": lightest_mass_ev,
        "MatchingScaleGeV": float(config["mu_initial_gev"]),
        "LowScaleGeV": float(config["mu_final_gev"]),
        "RequiredHighToLowC5NormRatio": scale_ratio,
        "RelativeHighScaleC5Residual": high_error,
        "MaxAbsY1": fit.max_abs_y1,
        "MaxAbsY2": fit.max_abs_y2,
        "PerturbativeBySqrt4Pi": bool(
            max(fit.max_abs_y1, fit.max_abs_y2)
            < np.sqrt(4.0 * np.pi)
        ),
        "LoopFactorsGeVInv": [
            {"re": float(value.real), "im": float(value.imag)}
            for value in loop_factors
        ],
    }

    return updated, diagnostics


def write_rge_corrected_config(
    c5_path: Path,
    input_config_path: Path,
    output_config_path: Path,
    *,
    lightest_mass_ev: float = 0.01,
) -> dict:
    """Write an RGE-corrected fitted numerical configuration."""

    input_config_path = Path(input_config_path)
    output_config_path = Path(output_config_path)

    config = json.loads(
        input_config_path.read_text(encoding="utf-8")
    )

    updated, diagnostics = build_rge_corrected_config(
        c5_path,
        config,
        lightest_mass_ev=lightest_mass_ev,
    )

    output_config_path.write_text(
        json.dumps(updated, indent=2),
        encoding="utf-8",
    )

    summary_path = output_config_path.with_name(
        output_config_path.stem + "_fit_summary.json"
    )
    summary_path.write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    diagnostics["OutputConfig"] = str(output_config_path)
    diagnostics["FitSummary"] = str(summary_path)

    return diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Construct T3 Yukawas that reproduce a low-scale neutrino target "
            "after one-loop SM RGE running."
        )
    )
    parser.add_argument(
        "c5",
        type=Path,
        help="Path to Matchete c5_coefficient.txt",
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Existing numerical JSON config",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("t3_numerical_fitted_rge.json"),
    )
    parser.add_argument(
        "--m-lightest",
        type=float,
        default=0.01,
    )
    args = parser.parse_args()

    summary = write_rge_corrected_config(
        args.c5,
        args.config,
        args.output,
        lightest_mass_ev=args.m_lightest,
    )

    print("=" * 72)
    print("T3 RGE-CORRECTED YUKAWA FIT")
    print("=" * 72)
    print(
        "required ||C5(M)|| / ||C5(low)|| =",
        summary["RequiredHighToLowC5NormRatio"],
    )
    print("max |y1| =", summary["MaxAbsY1"])
    print("max |y2| =", summary["MaxAbsY2"])
    print(
        "relative high-scale C5 residual =",
        summary["RelativeHighScaleC5Residual"],
    )
    print(
        "perturbative (< sqrt(4pi)) =",
        summary["PerturbativeBySqrt4Pi"],
    )
    print("output config =", summary["OutputConfig"])
    print("fit summary =", summary["FitSummary"])


if __name__ == "__main__":
    main()
