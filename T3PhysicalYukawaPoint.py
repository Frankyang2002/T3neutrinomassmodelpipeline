from __future__ import annotations

"""
Construct a neutrino-data-matched T3 Yukawa point from the actual Matchete loop
kernel for a chosen model and heavy spectrum.

The one-generation Matchete result is factorized as

    kappa_1g = F_loop * conjugate(y1) * conjugate(y2).

For three heavy generations, the model-dependent loop factors are

    F_r = F_loop(MF -> MF_r),

evaluated using the same masses and scalar coupling that the numerical pipeline
already uses.

Those F_r are then passed to the constructive three-generation Yukawa solution
in T3YukawaFit.py so that the resulting flavor C5 reproduces a chosen neutrino
target matrix.
"""

import argparse
import copy
import json
from pathlib import Path

import numpy as np
import sympy as sp

from FlavorMatchedC5 import extract_t3_loop_kernel
from MatchedEFTRGE import parse_matchete_c5
from T3NeutrinoTarget import build_normal_ordering_target
from T3YukawaFit import fit_three_heavy_balanced, reconstruct_c5


def evaluate_t3_loop_factors(
    c5_path: Path,
    config: dict,
) -> np.ndarray:
    """Evaluate the actual Matchete loop kernel once for each heavy mass."""

    c5_path = Path(c5_path)

    kappa_1g = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )
    kernel = extract_t3_loop_kernel(kappa_1g)

    model = config["t3"]
    heavy_masses = np.asarray(model["MF"], dtype=float)

    if len(heavy_masses) != 3:
        raise ValueError(
            "The constructive physical-point builder currently requires "
            "exactly three heavy generations."
        )

    base_substitutions = {
        sp.Symbol("MS1"): model["MS1"],
        sp.Symbol("MS2"): model["MS2"],
        sp.Symbol("lambdaT3"): model["lambdaT3"],
        sp.Symbol("hbar"): model.get(
            "hbar",
            1.0 / (16.0 * np.pi**2),
        ),
    }

    MF = sp.Symbol("MF")
    values = np.empty(3, dtype=complex)

    for r, mass in enumerate(heavy_masses):
        kernel_r = kernel.subs(
            {
                **base_substitutions,
                MF: float(mass),
            }
        )
        values[r] = complex(sp.N(kernel_r, 18))

    return values


def build_matched_yukawa_config(
    c5_path: Path,
    config: dict,
    *,
    lightest_mass_ev: float = 0.01,
    vev_gev: float | None = None,
) -> tuple[dict, dict]:
    """Return a copy of the numerical config with neutrino-matched y1 and y2."""

    if vev_gev is None:
        vev_gev = float(config.get("vev_gev", 246.22))

    loop_factors = evaluate_t3_loop_factors(
        c5_path,
        config,
    )

    target = build_normal_ordering_target(
        lightest_mass_ev=lightest_mass_ev,
        vev_gev=vev_gev,
    )

    fit = fit_three_heavy_balanced(
        target.c5_matrix_gev_inv,
        loop_factors,
    )

    reconstructed = reconstruct_c5(
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

    relative_target_error = float(
        np.linalg.norm(reconstructed - target.c5_matrix_gev_inv)
        / np.linalg.norm(target.c5_matrix_gev_inv)
    )

    diagnostics = {
        "Status": "Success",
        "Ordering": "NO",
        "LightestMassEV": lightest_mass_ev,
        "VevGeV": vev_gev,
        "LoopFactorsGeVInv": [
            {"re": float(value.real), "im": float(value.imag)}
            for value in loop_factors
        ],
        "RelativeC5Residual": fit.relative_residual,
        "RelativeTargetError": relative_target_error,
        "MaxAbsY1": fit.max_abs_y1,
        "MaxAbsY2": fit.max_abs_y2,
        "PerturbativeBySqrt4Pi": bool(
            max(fit.max_abs_y1, fit.max_abs_y2)
            < np.sqrt(4.0 * np.pi)
        ),
    }

    return updated, diagnostics


def write_matched_yukawa_config(
    c5_path: Path,
    input_config_path: Path,
    output_config_path: Path,
    *,
    lightest_mass_ev: float = 0.01,
) -> dict:
    """Create a numerical config whose T3 Yukawas reproduce neutrino data."""

    input_config_path = Path(input_config_path)
    output_config_path = Path(output_config_path)

    config = json.loads(
        input_config_path.read_text(encoding="utf-8")
    )

    updated, diagnostics = build_matched_yukawa_config(
        c5_path,
        config,
        lightest_mass_ev=lightest_mass_ev,
    )

    output_config_path.write_text(
        json.dumps(updated, indent=2),
        encoding="utf-8",
    )

    diagnostics_path = output_config_path.with_name(
        output_config_path.stem + "_fit_summary.json"
    )
    diagnostics_path.write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    diagnostics["OutputConfig"] = str(output_config_path)
    diagnostics["FitSummary"] = str(diagnostics_path)

    return diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Use the actual Matchete T3 loop kernel to construct Yukawa "
            "matrices matched to a NuFIT normal-ordering target."
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
        help="Existing numerical pipeline JSON config",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("t3_numerical_fitted.json"),
        help="Output numerical config containing the fitted Yukawas",
    )
    parser.add_argument(
        "--m-lightest",
        type=float,
        default=0.01,
        help="Normal-ordering lightest neutrino mass m1 in eV",
    )
    args = parser.parse_args()

    summary = write_matched_yukawa_config(
        args.c5,
        args.config,
        args.output,
        lightest_mass_ev=args.m_lightest,
    )

    print("=" * 72)
    print("T3 MATCHEte-KERNEL YUKAWA FIT")
    print("=" * 72)

    print("loop factors [GeV^-1]:")
    for i, value in enumerate(summary["LoopFactorsGeVInv"], start=1):
        print(f"  F{i} = {value['re']:+.12e} {value['im']:+.12e}j")

    print()
    print("max |y1| =", summary["MaxAbsY1"])
    print("max |y2| =", summary["MaxAbsY2"])
    print("relative C5 residual =", summary["RelativeC5Residual"])
    print("perturbative (< sqrt(4pi)) =", summary["PerturbativeBySqrt4Pi"])
    print("output config =", summary["OutputConfig"])
    print("fit summary =", summary["FitSummary"])


if __name__ == "__main__":
    main()
