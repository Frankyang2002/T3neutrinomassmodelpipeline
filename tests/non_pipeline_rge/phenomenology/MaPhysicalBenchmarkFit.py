from __future__ import annotations

"""Fit a high-scale Ma Yukawa matrix to low-energy neutrino data.

The scalar, gauge, SM Yukawa, and Majorana-mass inputs are held fixed.  Only
the complex 3x3 scotogenic Yukawa matrix h(mu_uv) is fitted.  Every objective
evaluation performs the complete calculation

    Ma UV running -> common-threshold matching -> SM Weinberg running.

To remove the complex-orthogonal redundancy of a quadratic Majorana matching
formula, the fit varies a complex symmetric right correction to an analytic
Takagi seed: h = h_seed @ (I + X), with X = X.T.
"""

import argparse
import copy
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from RGE.phenomenology.NeutrinoDataComparison import (
    compare_to_nufit6_io,
    compare_to_nufit6_no,
)
from RGE.phenomenology.NeutrinoObservables import calculate_neutrino_observables
from RGE.phenomenology.T3NeutrinoTarget import build_neutrino_target
from RGE.running.MaFullRunningComparison import (
    _threshold_state,
    compare_ma_running,
    state_from_config,
)
from RGE.running.MaUVRGE import (
    LOOP,
    MaUVInitialConditions,
    evolve_ma_uv,
    sample_ma_uv_trajectory,
    takagi_majorana,
)
from RGE.running.NumericalWeinbergRGE import evolve_weinberg, neutrino_mass_matrix


EV_TO_GEV = 1.0e-9
UPPER = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


@dataclass(frozen=True)
class MaPhysicalFitResult:
    fitted_config: dict[str, Any]
    summary: dict[str, Any]
    comparison: dict[str, Any]


def _target_c5(
    *,
    ordering: str,
    lightest_mass_ev: float,
    vev_gev: float,
    delta_cp: float,
    alpha21: float,
    alpha31: float,
):
    target = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=lightest_mass_ev,
        vev_gev=vev_gev,
        delta_cp=delta_cp,
        alpha21=alpha21,
        alpha31=alpha31,
    )
    # The project executes m_nu=-v^2 C5.  Keep the target matrix, including
    # its chosen phases, rather than relying only on the unphysical global sign.
    return target, -target.mass_matrix_ev * EV_TO_GEV / vev_gev**2


def _takagi_symmetric(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return U,s satisfying U.T @ matrix @ U = diag(s), s >= 0."""

    matrix = np.asarray(matrix, dtype=complex)
    _, singular_values, vh = np.linalg.svd(matrix)
    U = vh.conj().T
    candidate = U.T @ matrix @ U
    phases = np.exp(-0.5j * np.angle(np.diag(candidate)))
    U = U @ np.diag(phases)
    diagonalized = U.T @ matrix @ U
    singular_values = np.abs(np.diag(diagonalized))
    order = np.argsort(singular_values)
    U = U[:, order]
    singular_values = singular_values[order]
    residual = np.linalg.norm(U.T @ matrix @ U - np.diag(singular_values))
    residual /= max(1.0e-30, np.linalg.norm(matrix))
    if residual > 1.0e-8:
        raise RuntimeError(f"Target Takagi factorization failed: {residual:.3e}")
    return U, singular_values


def analytic_yukawa_seed(
    state: MaUVInitialConditions,
    target_c5: np.ndarray,
) -> np.ndarray:
    """Construct h satisfying the target at fixed parameters before running."""

    U_mass, masses = takagi_majorana(state.M)
    r = state.mEta2 / masses**2
    epsilon = r - 1.0
    f = np.empty(3, dtype=float)
    near = np.abs(epsilon) < 1.0e-5
    f[near] = -0.5 + epsilon[near] / 3.0 - epsilon[near] ** 2 / 4.0
    f[~near] = (
        np.log(r[~near]) / (1.0 - r[~near]) ** 2
        + 1.0 / (1.0 - r[~near])
    )
    weights = f / masses
    coefficients = state.lambda5 / LOOP * weights
    if np.any(np.abs(coefficients) < 1.0e-40):
        raise ValueError("The Ma matching coefficient is singular or zero.")

    U_c5, singular_values = _takagi_symmetric(target_c5)
    factor = np.diag(np.sqrt(singular_values)) @ U_c5.conj().T
    h_mass = np.diag(1.0 / np.sqrt(coefficients.astype(complex))) @ factor

    # h_mass = U_mass.T @ h in the Takagi convention used by MaUVRGE.
    h = U_mass.conj() @ h_mass
    return h


def _unpack_symmetric(parameters: np.ndarray) -> np.ndarray:
    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (12,):
        raise ValueError("The symmetric correction requires 12 real parameters.")
    result = np.zeros((3, 3), dtype=complex)
    for index, (i, j) in enumerate(UPPER):
        value = parameters[index] + 1j * parameters[index + len(UPPER)]
        result[i, j] = value
        result[j, i] = value
    return result


def _symmetric_residual(matrix: np.ndarray, scale: float) -> np.ndarray:
    """Pack a complex symmetric residual with Frobenius-equivalent weights."""

    values: list[complex] = []
    for i, j in UPPER:
        weight = np.sqrt(2.0) if i != j else 1.0
        values.append(weight * matrix[i, j] / scale)
    values = np.asarray(values, dtype=complex)
    return np.concatenate([values.real, values.imag])


def _full_low_c5(
    base_state: MaUVInitialConditions,
    h_high: np.ndarray,
    mu_uv: float,
    mu_match: float,
    mu_low: float,
) -> tuple[np.ndarray, MaUVInitialConditions]:
    state = MaUVInitialConditions(**{**base_state.__dict__, "h": h_high})
    uv = evolve_ma_uv(state, mu_uv, mu_match)
    sm, _metadata = _threshold_state(uv.state)
    eft = evolve_weinberg(sm, mu_match, mu_low)
    return eft.K, uv.state


def _observables_payload(c5: np.ndarray, vev_gev: float, ordering: str) -> dict[str, Any]:
    mass = neutrino_mass_matrix(c5, vev_gev=vev_gev)
    obs = calculate_neutrino_observables(mass, ordering=ordering)
    return {
        "Ordering": obs.ordering,
        "MassesEV": obs.masses_ev.tolist(),
        "DeltaM21SqEV2": obs.delta_m21_sq_ev2,
        "DeltaM31SqEV2": obs.delta_m31_sq_ev2,
        "DeltaM32SqEV2": obs.delta_m32_sq_ev2,
        "PMNSAbs": obs.pmns_abs.tolist(),
        "TakagiResidual": obs.takagi_residual,
    }


def _trajectory_checks(
    state: MaUVInitialConditions,
    mu_uv: float,
    mu_match: float,
    vev_gev: float,
) -> dict[str, Any]:
    """Check perturbativity, boundedness, and inert masses along the UV path."""

    trajectory = sample_ma_uv_trajectory(state, mu_uv, mu_match, samples=64)
    sqrt_4pi = float(np.sqrt(4.0 * np.pi))
    four_pi = float(4.0 * np.pi)
    minimum_bfb = {
        "lambda1": (np.inf, 0.0),
        "lambda2": (np.inf, 0.0),
        "lambda3_plus_sqrt": (np.inf, 0.0),
        "lambda345_plus_sqrt": (np.inf, 0.0),
    }
    minimum_inert_mass_sq = (np.inf, 0.0, "")
    maximum_yukawa = 0.0
    maximum_gauge = 0.0
    maximum_quartic = 0.0

    for scale, point in zip(trajectory.scales_gev, trajectory.states):
        root = np.sqrt(max(0.0, point.lambda1 * point.lambda2))
        margins = {
            "lambda1": point.lambda1,
            "lambda2": point.lambda2,
            "lambda3_plus_sqrt": point.lambda3 + root,
            "lambda345_plus_sqrt": (
                point.lambda3 + point.lambda4 - abs(point.lambda5) + root
            ),
        }
        for name, value in margins.items():
            if value < minimum_bfb[name][0]:
                minimum_bfb[name] = (float(value), float(scale))

        inert_masses = {
            "charged": point.mEta2 + 0.5 * vev_gev**2 * point.lambda3,
            "neutral_R": point.mEta2
            + 0.5 * vev_gev**2 * (point.lambda3 + point.lambda4 + point.lambda5),
            "neutral_I": point.mEta2
            + 0.5 * vev_gev**2 * (point.lambda3 + point.lambda4 - point.lambda5),
        }
        name, value = min(inert_masses.items(), key=lambda item: item[1])
        if value < minimum_inert_mass_sq[0]:
            minimum_inert_mass_sq = (float(value), float(scale), name)

        maximum_yukawa = max(
            maximum_yukawa,
            float(np.max(np.abs(point.Ye))),
            float(np.max(np.abs(point.Yu))),
            float(np.max(np.abs(point.Yd))),
            float(np.max(np.abs(point.h))),
        )
        maximum_gauge = max(maximum_gauge, abs(point.gY), abs(point.g2), abs(point.g3))
        maximum_quartic = max(
            maximum_quartic,
            abs(point.lambda1),
            abs(point.lambda2),
            abs(point.lambda3),
            abs(point.lambda4),
            abs(point.lambda5),
        )

    bfb_pass = all(value[0] > 0.0 for value in minimum_bfb.values())
    inert_pass = minimum_inert_mass_sq[0] > 0.0
    perturbative = (
        maximum_yukawa < sqrt_4pi
        and maximum_gauge < sqrt_4pi
        and maximum_quartic < four_pi
    )
    return {
        "samples": len(trajectory.states),
        "bounded_from_below_all_samples": bool(bfb_pass),
        "positive_inert_scalar_masses_all_samples": bool(inert_pass),
        "perturbative_all_samples": bool(perturbative),
        "all_checks_pass": bool(bfb_pass and inert_pass and perturbative),
        "minimum_bfb_margins": {
            name: {"value": value, "scale_gev": scale}
            for name, (value, scale) in minimum_bfb.items()
        },
        "minimum_inert_mass_sq_gev2": {
            "value": minimum_inert_mass_sq[0],
            "scale_gev": minimum_inert_mass_sq[1],
            "state": minimum_inert_mass_sq[2],
        },
        "maximum_absolute_yukawa": maximum_yukawa,
        "maximum_absolute_gauge_coupling": maximum_gauge,
        "maximum_absolute_quartic": maximum_quartic,
        "criteria": {
            "yukawa_and_gauge": "absolute value < sqrt(4 pi)",
            "quartic": "absolute value < 4 pi",
            "bounded_from_below": [
                "lambda1 > 0",
                "lambda2 > 0",
                "lambda3 + sqrt(lambda1 lambda2) > 0",
                "lambda3 + lambda4 - |lambda5| + sqrt(lambda1 lambda2) > 0",
            ],
        },
    }


def fit_ma_physical_benchmark(
    config: dict[str, Any],
    *,
    ordering: str = "NO",
    lightest_mass_ev: float = 0.001,
    delta_cp: float = 0.0,
    alpha21: float = 0.0,
    alpha31: float = 0.0,
    max_nfev: int = 300,
    target_tolerance: float = 1.0e-7,
) -> MaPhysicalFitResult:
    """Fit h(mu_uv) and return a new runnable configuration and diagnostics."""

    ordering = ordering.upper()
    mu_uv = float(config["mu_uv_gev"])
    mu_match = float(config["mu_matching_gev"])
    mu_low = float(config.get("mu_low_gev", 91.1876))
    vev_gev = float(config.get("vev_gev", 246.22))
    base_state = state_from_config(config["uv"])
    target, target_c5 = _target_c5(
        ordering=ordering,
        lightest_mass_ev=lightest_mass_ev,
        vev_gev=vev_gev,
        delta_cp=delta_cp,
        alpha21=alpha21,
        alpha31=alpha31,
    )
    target_norm = float(np.linalg.norm(target_c5))
    seed = analytic_yukawa_seed(base_state, target_c5)

    evaluations = 0

    def objective(parameters: np.ndarray) -> np.ndarray:
        nonlocal evaluations
        evaluations += 1
        correction = np.eye(3, dtype=complex) + _unpack_symmetric(parameters)
        h_high = seed @ correction
        try:
            low_c5, _threshold = _full_low_c5(
                base_state,
                h_high,
                mu_uv,
                mu_match,
                mu_low,
            )
            return _symmetric_residual(low_c5 - target_c5, target_norm)
        except (ValueError, RuntimeError, FloatingPointError):
            return np.full(12, 1.0e6, dtype=float)

    solution = least_squares(
        objective,
        np.zeros(12, dtype=float),
        method="trf",
        x_scale="jac",
        ftol=1.0e-11,
        xtol=1.0e-11,
        gtol=1.0e-11,
        max_nfev=max_nfev,
    )
    correction = np.eye(3, dtype=complex) + _unpack_symmetric(solution.x)
    fitted_h = seed @ correction
    final_c5, threshold_state = _full_low_c5(
        base_state,
        fitted_h,
        mu_uv,
        mu_match,
        mu_low,
    )
    relative_residual = float(np.linalg.norm(final_c5 - target_c5) / target_norm)
    if not solution.success or relative_residual > target_tolerance:
        raise RuntimeError(
            "Physical Ma fit did not converge: "
            f"solver={solution.message}; relative residual={relative_residual:.3e}."
        )

    fitted_config = copy.deepcopy(config)
    fitted_config["uv"]["h"] = {
        "real": fitted_h.real.tolist(),
        "imag": fitted_h.imag.tolist(),
    }
    fitted_config["fit_target"] = {
        "ordering": ordering,
        "lightest_mass_ev": float(lightest_mass_ev),
        "delta_cp_rad": float(delta_cp),
        "alpha21_rad": float(alpha21),
        "alpha31_rad": float(alpha31),
    }

    observables = _observables_payload(final_c5, vev_gev, ordering)
    data_comparison = (
        compare_to_nufit6_no(observables)
        if ordering == "NO"
        else compare_to_nufit6_io(observables)
    )
    comparison = compare_ma_running(fitted_config)
    fitted_high_state = MaUVInitialConditions(**{**base_state.__dict__, "h": fitted_h})
    theory_checks = _trajectory_checks(
        fitted_high_state,
        mu_uv,
        mu_match,
        vev_gev,
    )

    summary = {
        "status": "Success",
        "ordering": ordering,
        "target": {
            "lightest_mass_ev": float(lightest_mass_ev),
            "masses_ev": target.masses_ev.tolist(),
            "delta_cp_rad": float(delta_cp),
            "alpha21_rad": float(alpha21),
            "alpha31_rad": float(alpha31),
        },
        "fit": {
            "solver_success": bool(solution.success),
            "solver_message": str(solution.message),
            "function_evaluations": int(evaluations),
            "reported_nfev": int(solution.nfev),
            "relative_low_scale_c5_residual": relative_residual,
            "max_abs_h_uv": float(np.max(np.abs(fitted_h))),
            "max_abs_h_matching": float(np.max(np.abs(threshold_state.h))),
            "perturbative_below_sqrt_4pi": bool(
                max(np.max(np.abs(fitted_h)), np.max(np.abs(threshold_state.h)))
                < np.sqrt(4.0 * np.pi)
            ),
        },
        "low_energy_observables": observables,
        "nufit_comparison": data_comparison,
        "uv_theory_checks": theory_checks,
        "running_effects": comparison["effects"],
    }
    return MaPhysicalFitResult(
        fitted_config=fitted_config,
        summary=summary,
        comparison=comparison,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="unfitted Ma benchmark JSON")
    parser.add_argument("--output", type=Path, required=True, help="fitted config JSON")
    parser.add_argument("--summary", type=Path, default=None, help="fit summary JSON")
    parser.add_argument("--comparison", type=Path, default=None, help="four-way comparison JSON")
    parser.add_argument("--ordering", choices=("NO", "IO"), default="NO")
    parser.add_argument("--m-lightest", type=float, default=0.001)
    parser.add_argument("--delta-cp", type=float, default=0.0)
    parser.add_argument("--alpha21", type=float, default=0.0)
    parser.add_argument("--alpha31", type=float, default=0.0)
    parser.add_argument("--max-nfev", type=int, default=300)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = fit_ma_physical_benchmark(
        config,
        ordering=args.ordering,
        lightest_mass_ev=args.m_lightest,
        delta_cp=args.delta_cp,
        alpha21=args.alpha21,
        alpha31=args.alpha31,
        max_nfev=args.max_nfev,
    )
    summary_path = args.summary or args.output.with_name(args.output.stem + "_fit_summary.json")
    comparison_path = args.comparison or args.output.with_name(
        args.output.stem + "_comparison.json"
    )
    for path in (args.output, summary_path, comparison_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result.fitted_config, indent=2) + "\n", encoding="utf-8")
    summary_path.write_text(json.dumps(result.summary, indent=2) + "\n", encoding="utf-8")
    comparison_path.write_text(json.dumps(result.comparison, indent=2) + "\n", encoding="utf-8")

    fit = result.summary["fit"]
    obs = result.summary["low_energy_observables"]
    print(f"Wrote fitted config: {args.output}")
    print(f"Wrote fit summary: {summary_path}")
    print(f"Wrote comparison: {comparison_path}")
    print(f"relative low-scale C5 residual={fit['relative_low_scale_c5_residual']:.6g}")
    print(f"max |h(mu_uv)|={fit['max_abs_h_uv']:.6g}")
    print(f"masses [eV]={obs['MassesEV']}")
    print(f"all five NuFIT observables inside 3 sigma={result.summary['nufit_comparison']['AllFiveInside3Sigma']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
