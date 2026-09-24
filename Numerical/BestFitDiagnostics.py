"""
Physics diagnostics for a fitted T3 oscillation point.

This module reconstructs the best parameter point from an
OscillationOptimizer result and reruns the authoritative numerical pipeline
while retaining the final SM+Weinberg trajectory.

It produces:
- low-energy neutrino masses and their sum;
- the matched C5 matrix at the scalar threshold;
- the low-energy C5 and m_nu matrices;
- |PMNS|;
- scale-dependent neutrino masses;
- scale-dependent |C5_ij|;
- threshold mass diagnostics;
- standard thesis-ready plots.

No matching formula or oscillation observable definition is duplicated here.
The module reuses FinalC5Bridge, ScanCLI, T3Trajectory, WeinbergTrajectory and
RGE.phenomenology.NeutrinoObservables.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np


def _load_json_object(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return payload


def load_optimizer_best_parameters(
    optimizer_result_path: Path,
) -> tuple[dict[str, float], float, str]:
    """Return (parameters, chi2, ordering) from a successful optimiser result."""

    payload = _load_json_object(Path(optimizer_result_path))

    if payload.get("status") != "Success":
        raise ValueError(
            "Optimizer result must have status='Success'."
        )

    best = payload.get("best")
    if not isinstance(best, dict):
        raise ValueError(
            "Optimizer result must contain a 'best' object."
        )

    raw_parameters = best.get("parameters")
    if not isinstance(raw_parameters, dict) or not raw_parameters:
        raise ValueError(
            "Optimizer best point must contain a non-empty parameter mapping."
        )

    parameters = {
        str(name): float(value)
        for name, value in raw_parameters.items()
    }

    if not np.all(np.isfinite(list(parameters.values()))):
        raise ValueError(
            "Optimizer best parameters must all be finite."
        )

    chi2 = float(best["chi2"])
    if not np.isfinite(chi2) or chi2 < 0.0:
        raise ValueError("Optimizer best chi2 must be finite and non-negative.")

    ordering = str(best.get("ordering", "AUTO")).upper()
    if ordering not in {"NO", "IO", "AUTO"}:
        raise ValueError(
            "Optimizer best ordering must be NO, IO or AUTO."
        )

    return parameters, chi2, ordering


def _complex_matrix_payload(matrix: np.ndarray) -> dict[str, Any]:
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.ndim != 2:
        raise ValueError("Complex matrix payload requires a 2D matrix.")
    if not np.all(np.isfinite(matrix.real)) or not np.all(
        np.isfinite(matrix.imag)
    ):
        raise ValueError("Complex matrix contains non-finite entries.")

    return {
        "real": matrix.real.tolist(),
        "imag": matrix.imag.tolist(),
        "abs": np.abs(matrix).tolist(),
    }


def _save_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def _plot_matrix_abs(
    matrix: np.ndarray,
    *,
    title: str,
    output_path: Path,
    scientific: bool = True,
) -> None:
    values = np.abs(np.asarray(matrix, dtype=complex))

    plt.figure(figsize=(5.4, 4.8))
    image = plt.imshow(values)
    plt.colorbar(image, label="Absolute value")
    plt.xticks(range(3), ["1", "2", "3"])
    plt.yticks(range(3), ["1", "2", "3"])
    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.title(title)

    for i in range(3):
        for j in range(3):
            text = (
                f"{values[i, j]:.2e}"
                if scientific
                else f"{values[i, j]:.4g}"
            )
            plt.text(
                j,
                i,
                text,
                ha="center",
                va="center",
            )

    _save_figure(output_path)


def _plot_pmns_abs(
    pmns_abs: np.ndarray,
    output_path: Path,
) -> None:
    values = np.asarray(pmns_abs, dtype=float)

    plt.figure(figsize=(5.4, 4.8))
    image = plt.imshow(values, vmin=0.0, vmax=1.0)
    plt.colorbar(image, label=r"$|U_{\rm PMNS}|$")
    plt.xticks(range(3), ["1", "2", "3"])
    plt.yticks(range(3), [r"$e$", r"$\mu$", r"$\tau$"])
    plt.xlabel("Mass eigenstate")
    plt.ylabel("Charged-lepton flavour")
    plt.title(r"$|U_{\rm PMNS}|$")

    for i in range(3):
        for j in range(3):
            plt.text(
                j,
                i,
                f"{values[i, j]:.4f}",
                ha="center",
                va="center",
            )

    _save_figure(output_path)


def _plot_low_energy_masses(
    masses_ev: np.ndarray,
    output_path: Path,
) -> None:
    masses = np.asarray(masses_ev, dtype=float)

    plt.figure(figsize=(6.2, 4.6))
    x = np.arange(3)
    plt.bar(x, masses)
    plt.xticks(x, [r"$m_1$", r"$m_2$", r"$m_3$"])
    plt.ylabel("Mass [eV]")
    plt.title("Low-energy neutrino masses")
    plt.grid(True, axis="y", alpha=0.25)
    _save_figure(output_path)


def _plot_mass_running(
    scales_gev: np.ndarray,
    masses_ev: np.ndarray,
    output_path: Path,
) -> None:
    scales = np.asarray(scales_gev, dtype=float)
    masses = np.asarray(masses_ev, dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    for index, label in enumerate(
        [r"$m_1$", r"$m_2$", r"$m_3$"]
    ):
        plt.plot(
            scales,
            masses[:, index],
            label=label,
            linewidth=1.4,
        )

    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel("Neutrino mass [eV]")
    plt.title("Neutrino-mass running in the final SM+Weinberg EFT")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save_figure(output_path)


def _plot_c5_running(
    scales_gev: np.ndarray,
    c5_matrices: np.ndarray,
    output_path: Path,
) -> None:
    scales = np.asarray(scales_gev, dtype=float)
    c5 = np.asarray(c5_matrices, dtype=complex)

    plt.figure(figsize=(7.6, 5.2))

    for i in range(3):
        for j in range(i, 3):
            values = np.abs(c5[:, i, j])
            positive = values > 0.0

            if not np.any(positive):
                continue

            plt.plot(
                scales[positive],
                values[positive],
                label=rf"$|C_5^{{{i+1}{j+1}}}|$",
                linewidth=1.0,
            )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$|C_5^{ij}|$ [GeV$^{-1}$]")
    plt.title(r"Running Weinberg coefficient $C_5$")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save_figure(output_path)


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


def run_best_fit_diagnostics(
    *,
    scan_config_path: Path,
    optimizer_result_path: Path,
    output_dir: Path,
    n_scale_points: int = 120,
) -> Path:
    """Recompute and characterize the best fitted T3 point."""

    # Lazy imports keep the pure parsing/plot helpers independently testable.
    from Numerical.FinalC5Bridge import FinalC5TrajectoryBuilder
    from Numerical.FinalSMBoundary import build_weinberg_initial_conditions
    from Numerical.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.ScanCLI import (
        build_uv_state_from_config,
        load_scan_config,
        validate_scan_inputs,
    )
    from Numerical.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.WeinbergTrajectory import (
        charged_lepton_mass_basis_matrix,
        neutrino_mass_from_c5,
        run_weinberg_trajectory,
        scale_dependent_neutrino_observables,
    )
    from RGE.phenomenology.NeutrinoObservables import (
        calculate_neutrino_observables,
    )

    config = load_scan_config(Path(scan_config_path))
    validate_scan_inputs(config)

    parameters, stored_chi2, stored_ordering = (
        load_optimizer_best_parameters(
            Path(optimizer_result_path)
        )
    )

    raw = config.raw
    scales = raw["scales"]
    numerical = raw.get("numerical", {})
    if not isinstance(numerical, dict):
        raise ValueError("numerical must be an object when supplied.")

    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    mu_low = float(scales["mu_low_gev"])
    vev = float(numerical.get("vev_gev", 246.22))
    rtol = float(numerical.get("rtol", 1.0e-8))
    atol = float(numerical.get("atol", 1.0e-11))

    uv_payload = load_rgbeta_payload(config.uv_rgbeta_path)
    eft1_payload = load_rgbeta_payload(config.eft1_rgbeta_path)

    uv_state = build_uv_state_from_config(
        config,
        parameters,
    ).validated()

    renormalisable = run_renormalisable_t3_trajectory(
        uv_state,
        uv_payload,
        eft1_payload,
        mu_fermion_threshold_gev=mu_f,
        mu_scalar_threshold_gev=mu_s,
        rtol=rtol,
        atol=atol,
    )

    c5_builder = FinalC5TrajectoryBuilder(
        config.final_weinberg_path
    )
    c5_threshold = np.asarray(
        c5_builder(parameters, renormalisable),
        dtype=complex,
    )

    initial = build_weinberg_initial_conditions(
        renormalisable.final_sm_boundary,
        c5_threshold,
    )

    save_scales = _log_save_scales(
        mu_s,
        mu_low,
        n_scale_points,
    )

    final_trajectory = run_weinberg_trajectory(
        initial,
        mu_s,
        mu_low,
        save_scales_gev=save_scales,
        rtol=rtol,
        atol=atol,
    )

    ordering = stored_ordering
    if ordering == "AUTO":
        ordering = "AUTO"

    scale_points = scale_dependent_neutrino_observables(
        final_trajectory,
        vev_gev=vev,
        ordering=ordering,
    )

    low_running = final_trajectory.final_point

    mnu_low_weak = neutrino_mass_from_c5(
        low_running.c5,
        vev_gev=vev,
    )
    mnu_low_flavor = charged_lepton_mass_basis_matrix(
        mnu_low_weak,
        low_running.ye,
    )

    observables = calculate_neutrino_observables(
        mnu_low_flavor,
        ordering=ordering,
    )

    masses_ev = np.asarray(
        observables.masses_ev,
        dtype=float,
    )

    mass_running = np.asarray(
        [
            np.asarray(point.observables.masses_ev, dtype=float)
            for point in scale_points
        ],
        dtype=float,
    )
    c5_running = np.asarray(
        [point.c5 for point in scale_points],
        dtype=complex,
    )
    running_scales = np.asarray(
        [point.mu_gev for point in scale_points],
        dtype=float,
    )

    scalar_diag = renormalisable.scalar_threshold
    fermion_diag = renormalisable.fermion_threshold

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    physics_payload = {
        "status": "Success",
        "scan_config": str(config.source_path),
        "optimizer_result": str(Path(optimizer_result_path).resolve()),
        "stored_best_chi2": stored_chi2,
        "ordering": observables.ordering,
        "vev_gev": vev,
        "scales_gev": {
            "mu_uv": float(uv_state.mu_gev),
            "mu_fermion_threshold": mu_f,
            "mu_scalar_threshold": mu_s,
            "mu_low": mu_low,
        },
        "best_parameters": parameters,
        "threshold_diagnostics": {
            "fermion_singular_masses_gev": np.asarray(
                fermion_diag.masses_gev,
                dtype=float,
            ).tolist(),
            "scalar_running_masses_gev": list(
                scalar_diag.masses_gev
            ),
        },
        "low_energy": {
            "masses_ev": masses_ev.tolist(),
            "sum_masses_ev": float(np.sum(masses_ev)),
            "delta_m21_sq_ev2": float(
                observables.delta_m21_sq_ev2
            ),
            "delta_m31_sq_ev2": float(
                observables.delta_m31_sq_ev2
            ),
            "delta_m32_sq_ev2": float(
                observables.delta_m32_sq_ev2
            ),
            "pmns_abs": np.asarray(
                observables.pmns_abs,
                dtype=float,
            ).tolist(),
            "takagi_residual": float(
                observables.takagi_residual
            ),
            "c5": _complex_matrix_payload(
                low_running.c5
            ),
            "mnu_weak_basis_gev": _complex_matrix_payload(
                mnu_low_weak
            ),
            "mnu_charged_lepton_basis_gev": (
                _complex_matrix_payload(mnu_low_flavor)
            ),
        },
        "scalar_threshold": {
            "c5": _complex_matrix_payload(c5_threshold),
        },
        "running": {
            "mu_gev": running_scales.tolist(),
            "masses_ev": mass_running.tolist(),
            "c5_abs": np.abs(c5_running).tolist(),
        },
    }

    json_path = output_dir / "best_fit_physics.json"
    json_path.write_text(
        json.dumps(physics_payload, indent=2),
        encoding="utf-8",
    )

    _plot_low_energy_masses(
        masses_ev,
        output_dir / "neutrino_masses.png",
    )
    _plot_mass_running(
        running_scales,
        mass_running,
        output_dir / "neutrino_mass_running.png",
    )
    _plot_c5_running(
        running_scales,
        c5_running,
        output_dir / "c5_running.png",
    )
    _plot_matrix_abs(
        mnu_low_flavor * 1.0e9,
        title=r"$|m_\nu|$ in charged-lepton basis [eV]",
        output_path=output_dir / "mnu_matrix_abs.png",
    )
    _plot_matrix_abs(
        low_running.c5,
        title=r"Low-energy $|C_5|$ [GeV$^{-1}$]",
        output_path=output_dir / "c5_matrix_abs.png",
    )
    _plot_pmns_abs(
        observables.pmns_abs,
        output_dir / "pmns_abs.png",
    )

    summary_lines = [
        "T3 best-fit physics diagnostics",
        "",
        f"stored_best_chi2: {stored_chi2:.12g}",
        f"ordering: {observables.ordering}",
        f"mu_F [GeV]: {mu_f:.12g}",
        f"mu_S [GeV]: {mu_s:.12g}",
        f"mu_low [GeV]: {mu_low:.12g}",
        "",
        "Low-energy neutrino masses [eV]:",
        f"  m1 = {masses_ev[0]:.12g}",
        f"  m2 = {masses_ev[1]:.12g}",
        f"  m3 = {masses_ev[2]:.12g}",
        f"  sum = {np.sum(masses_ev):.12g}",
        "",
        "Mass splittings [eV^2]:",
        f"  Delta m21^2 = {observables.delta_m21_sq_ev2:.12g}",
        f"  Delta m31^2 = {observables.delta_m31_sq_ev2:.12g}",
        f"  Delta m32^2 = {observables.delta_m32_sq_ev2:.12g}",
        "",
        f"Takagi residual = {observables.takagi_residual:.12g}",
        "",
        "Heavy-fermion singular masses at mu_F [GeV]:",
    ]
    summary_lines.extend(
        f"  {value:.12g}"
        for value in fermion_diag.masses_gev
    )
    summary_lines.append("")
    summary_lines.append(
        "Running scalar masses at mu_S [GeV]:"
    )
    summary_lines.extend(
        f"  {value:.12g}"
        for value in scalar_diag.masses_gev
    )

    (output_dir / "summary.txt").write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print(f"Diagnostics written to: {output_dir}")
    print(
        "Low-energy masses [eV]: "
        + ", ".join(f"{x:.8g}" for x in masses_ev)
    )
    print(f"Sum m_nu [eV]: {np.sum(masses_ev):.8g}")

    return output_dir


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute and characterize the best fitted T3 oscillation point."
        )
    )
    parser.add_argument(
        "scan_config",
        type=Path,
        help=(
            "Scan config defining the representation, base state, bindings, "
            "matching files and scales."
        ),
    )
    parser.add_argument(
        "--optimizer-result",
        type=Path,
        required=True,
        help="Successful OscillationOptimizer result JSON.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for physics JSON, summary and figures.",
    )
    parser.add_argument(
        "--n-scale-points",
        type=int,
        default=120,
        help=(
            "Number of logarithmically spaced final-EFT trajectory points. "
            "Default: 120."
        ),
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    run_best_fit_diagnostics(
        scan_config_path=args.scan_config,
        optimizer_result_path=args.optimizer_result,
        output_dir=args.output_dir,
        n_scale_points=args.n_scale_points,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
