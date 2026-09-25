"""T3 running diagnostics for the fitted T3 benchmark.

This is a diagnostics layer only. It reuses the authoritative numerical pipeline
and records the scale-dependent quantities that are most useful for thesis plots:

- UV coupling measures before the F threshold;
- renormalisable scalar-only couplings between F and the scalar threshold;
- final-SM+Weinberg neutrino masses, mass splittings and mixing angles;
- final-SM+Weinberg C5 and charged-lepton-basis m_nu matrix elements.

The module deliberately does not invent a numerical LLSS Wilson trajectory. The
current numerical state stores only the renormalisable scalar-only EFT; the
higher-dimensional intermediate Wilson running remains in the RGE/matching layer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.BestFitDiagnostics import load_optimizer_best_parameters
from Numerical.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.RGBetaEvaluator import load_rgbeta_payload
from Numerical.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
from Numerical.ScanCLI import (
    build_uv_state_from_config,
    load_scan_config,
    validate_scan_inputs,
)
from Numerical.T3Trajectory import run_renormalisable_t3_trajectory
from Numerical.WeinbergTrajectory import run_weinberg_trajectory
from physics.NeutrinoObservables import calculate_neutrino_observables
from physics.NeutrinoTrajectory import (
    charged_lepton_mass_basis_matrix,
    neutrino_mass_from_c5,
    scale_dependent_neutrino_observables,
)


def _log_save_scales(
    mu_initial_gev: float,
    mu_final_gev: float,
    n_points: int,
) -> np.ndarray:
    count = int(n_points)
    if count < 2:
        raise ValueError("n_scale_points must be at least 2.")

    mu_initial = float(mu_initial_gev)
    mu_final = float(mu_final_gev)
    if (
        not np.isfinite(mu_initial)
        or not np.isfinite(mu_final)
        or mu_initial <= 0.0
        or mu_final <= 0.0
    ):
        raise ValueError("Trajectory scales must be finite and positive.")

    return np.geomspace(mu_initial, mu_final, count)


def mixing_angles_from_pmns_abs(
    pmns_abs: np.ndarray,
) -> tuple[float, float, float]:
    """Return sin^2(theta12), sin^2(theta13), sin^2(theta23)."""

    values = np.asarray(pmns_abs, dtype=float)
    if values.shape != (3, 3):
        raise ValueError("pmns_abs must be a 3x3 matrix.")

    sin2_theta13 = float(values[0, 2] ** 2)
    denominator = 1.0 - sin2_theta13
    if denominator <= 0.0:
        raise RuntimeError(
            "Cannot extract theta12/theta23 because 1-|Ue3|^2 <= 0."
        )

    return (
        float(values[0, 1] ** 2 / denominator),
        sin2_theta13,
        float(values[1, 2] ** 2 / denominator),
    )


def _complex_matrix_payload(matrix: np.ndarray) -> dict[str, Any]:
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError("Expected a 3x3 complex matrix.")
    return {
        "real": matrix.real.tolist(),
        "imag": matrix.imag.tolist(),
        "abs": np.abs(matrix).tolist(),
    }


def _uv_running_payload(result) -> dict[str, Any]:
    mu: list[float] = []
    y1_norm: list[float] = []
    y2_norm: list[float] = []
    lambda_t3_abs: list[float] = []
    mf_singular: list[list[float]] = []

    for index in range(result.n_points):
        state = result.state_at_index(index)
        if not hasattr(state, "y1") or not hasattr(state, "y2"):
            raise TypeError(
                "RunningDiagnostics currently supports the ordinary split-scalar scan branch."
            )

        mu.append(float(state.mu_gev))
        y1_norm.append(float(np.linalg.norm(state.y1, ord="fro")))
        y2_norm.append(float(np.linalg.norm(state.y2, ord="fro")))
        lambda_t3_abs.append(float(abs(state.lambdaT3)))
        mf_singular.append(
            np.sort(
                np.linalg.svd(
                    np.asarray(state.MF, dtype=complex),
                    compute_uv=False,
                )
            ).astype(float).tolist()
        )

    return {
        "mu_gev": mu,
        "y1_frobenius_norm": y1_norm,
        "y2_frobenius_norm": y2_norm,
        "lambdaT3_abs": lambda_t3_abs,
        "fermion_singular_masses_gev": mf_singular,
    }


def _intermediate_running_payload(result) -> dict[str, Any]:
    mu: list[float] = []
    lambda_h1: list[float] = []
    lambda_h2: list[float] = []
    lambda12: list[float] = []
    lambda_t3_abs: list[float] = []
    scalar_mass_sq: list[list[float]] = []

    for index in range(result.n_points):
        state = result.state_at_index(index)
        if not hasattr(state, "lambdaH1") or not hasattr(state, "lambdaH2"):
            raise TypeError(
                "RunningDiagnostics currently supports the ordinary split-scalar scan branch."
            )

        mu.append(float(state.mu_gev))
        lambda_h1.append(float(state.lambdaH1))
        lambda_h2.append(float(state.lambdaH2))
        lambda12.append(float(state.lambda12))
        lambda_t3_abs.append(float(abs(state.lambdaT3)))
        scalar_mass_sq.append([float(state.mS1Sq), float(state.mS2Sq)])

    return {
        "mu_gev": mu,
        "lambdaH1": lambda_h1,
        "lambdaH2": lambda_h2,
        "lambda12": lambda12,
        "lambdaT3_abs": lambda_t3_abs,
        "scalar_mass_squared_gev2": scalar_mass_sq,
    }


def run_running_diagnostics(
    *,
    scan_config_path: Path,
    optimizer_result_path: Path,
    output_path: Path,
    n_scale_points: int = 120,
) -> Path:
    """Rerun the best-fit point and save plotting-oriented running data."""

    config = load_scan_config(Path(scan_config_path))
    validate_scan_inputs(config)

    parameters, stored_chi2, stored_ordering = load_optimizer_best_parameters(
        Path(optimizer_result_path)
    )

    scales = config.raw["scales"]
    numerical = config.raw.get("numerical", {})
    if not isinstance(numerical, dict):
        raise ValueError("numerical must be an object when supplied.")

    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(scales["mu_low_gev"])
    vev = float(numerical.get("vev_gev", 246.22))
    rtol = float(numerical.get("rtol", 1.0e-8))
    atol = float(numerical.get("atol", 1.0e-11))

    uv_payload = load_rgbeta_payload(config.uv_rgbeta_path)
    intermediate_payload = load_rgbeta_payload(config.eft1_rgbeta_path)

    uv_state = build_uv_state_from_config(config, parameters).validated()

    renormalisable = run_renormalisable_t3_trajectory(
        uv_state,
        uv_payload,
        intermediate_payload,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        uv_save_scales_gev=_log_save_scales(
            uv_state.mu_gev,
            mu_f,
            n_scale_points,
        ),
        eft1_save_scales_gev=_log_save_scales(
            mu_f,
            mu_s,
            n_scale_points,
        ),
        rtol=rtol,
        atol=atol,
    )

    c5_evaluator = FinalC5TrajectoryEvaluator(config.final_weinberg_path)
    c5_threshold = np.asarray(
        c5_evaluator(parameters, renormalisable),
        dtype=complex,
    )

    initial = build_sm_weinberg_initial_conditions(
        renormalisable.final_sm_boundary,
        c5_threshold,
    )

    final_trajectory = run_weinberg_trajectory(
        initial,
        mu_s,
        mu_low,
        save_scales_gev=_log_save_scales(mu_s, mu_low, n_scale_points),
        rtol=rtol,
        atol=atol,
    )

    scale_points = scale_dependent_neutrino_observables(
        final_trajectory,
        vev_gev=vev,
        ordering=stored_ordering,
    )

    running_scales: list[float] = []
    masses: list[list[float]] = []
    c5_abs: list[list[list[float]]] = []
    mnu_flavor_abs_ev: list[list[list[float]]] = []
    dm21: list[float] = []
    dm3l: list[float] = []
    s12: list[float] = []
    s13: list[float] = []
    s23: list[float] = []

    for index, point in enumerate(scale_points):
        obs = point.observables
        running = final_trajectory.point_at_index(index)

        running_scales.append(float(point.mu_gev))
        masses.append(np.asarray(obs.masses_ev, dtype=float).tolist())
        c5_abs.append(np.abs(point.c5).tolist())
        dm21.append(float(obs.delta_m21_sq_ev2))
        dm3l.append(
            float(
                obs.delta_m31_sq_ev2
                if obs.ordering.upper() == "NO"
                else obs.delta_m32_sq_ev2
            )
        )
        angle12, angle13, angle23 = mixing_angles_from_pmns_abs(obs.pmns_abs)
        s12.append(angle12)
        s13.append(angle13)
        s23.append(angle23)

        mass_flavor = charged_lepton_mass_basis_matrix(
            point.mass_matrix_gev,
            running.ye,
        )
        mnu_flavor_abs_ev.append(
            (np.abs(mass_flavor) * 1.0e9).tolist()
        )

    low = final_trajectory.final_point
    mnu_low = charged_lepton_mass_basis_matrix(
        neutrino_mass_from_c5(low.c5, vev_gev=vev),
        low.ye,
    )
    low_observables = calculate_neutrino_observables(
        mnu_low,
        ordering=stored_ordering,
    )

    payload = {
        "status": "Success",
        "source": "RunningDiagnostics",
        "scan_config": str(config.source_path),
        "optimizer_result": str(Path(optimizer_result_path).resolve()),
        "stored_best_chi2": stored_chi2,
        "ordering": low_observables.ordering,
        "vev_gev": vev,
        "scales_gev": {
            "mu_uv": float(uv_state.mu_gev),
            "mu_fermion_threshold": mu_f,
            "mu_scalar_threshold": mu_s,
            "mu_low": mu_low,
        },
        "best_parameters": parameters,
        "uv_running": _uv_running_payload(renormalisable.uv),
        "intermediate_running": _intermediate_running_payload(
            renormalisable.intermediate
        ),
        "scalar_threshold": {
            "c5": _complex_matrix_payload(c5_threshold),
        },
        "final_running": {
            "mu_gev": running_scales,
            "masses_ev": masses,
            "c5_abs": c5_abs,
            "mnu_charged_lepton_abs_ev": mnu_flavor_abs_ev,
            "delta_m21_sq_ev2": dm21,
            "delta_m3l_sq_ev2": dm3l,
            "sin2_theta12": s12,
            "sin2_theta13": s13,
            "sin2_theta23": s23,
        },
        "low_energy": {
            "masses_ev": np.asarray(
                low_observables.masses_ev,
                dtype=float,
            ).tolist(),
            "delta_m21_sq_ev2": float(low_observables.delta_m21_sq_ev2),
            "delta_m31_sq_ev2": float(low_observables.delta_m31_sq_ev2),
            "delta_m32_sq_ev2": float(low_observables.delta_m32_sq_ev2),
            "pmns_abs": np.asarray(
                low_observables.pmns_abs,
                dtype=float,
            ).tolist(),
            "mnu_charged_lepton_basis_gev": _complex_matrix_payload(mnu_low),
            "c5": _complex_matrix_payload(low.c5),
            "takagi_residual": float(low_observables.takagi_residual),
        },
        "notes": {
            "intermediate_dimension_five_running": (
                "Not numerically sampled here. The scalar-only numerical state "
                "contains renormalisable parameters only; LLSS Wilson running "
                "remains in the RGE/matching pipeline."
            ),
            "uv_yukawa_display": (
                "Frobenius norms are stored for y1 and y2 to summarize their "
                "overall scale without selecting arbitrary matrix entries."
            ),
        },
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary.replace(output_path)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate T3 running diagnostics for a fitted T3 point."
    )
    parser.add_argument("scan_config", type=Path)
    parser.add_argument("--optimizer-result", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n-scale-points", type=int, default=120)
    args = parser.parse_args()

    output = run_running_diagnostics(
        scan_config_path=args.scan_config,
        optimizer_result_path=args.optimizer_result,
        output_path=args.output,
        n_scale_points=args.n_scale_points,
    )
    print(f"T3 running diagnostics written to: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
