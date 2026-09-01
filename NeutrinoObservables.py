from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


GEV_TO_EV = 1.0e9


@dataclass(frozen=True)
class NeutrinoObservables:
    masses_ev: np.ndarray
    delta_m21_sq_ev2: float
    delta_m31_sq_ev2: float
    pmns_abs: np.ndarray
    takagi_residual: float


def takagi_factorization(
    mass_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Takagi-factorize a complex symmetric Majorana mass matrix.

    Returns U, masses, residual such that

        U.T @ M @ U ~= diag(masses),

    with non-negative masses sorted from smallest to largest.
    """

    M = np.asarray(mass_matrix, dtype=complex)

    if M.shape != (3, 3):
        raise ValueError("Neutrino mass matrix must be 3x3.")

    if not np.allclose(M, M.T, rtol=1e-9, atol=1e-18):
        raise ValueError("Neutrino mass matrix must be symmetric.")

    # For a Takagi decomposition, the right singular vectors are the
    # Takagi vectors up to column phases.
    _, singular_values, vh = np.linalg.svd(M)
    U = vh.conj().T

    diagonal_candidate = U.T @ M @ U

    phases = np.ones(3, dtype=complex)

    for i in range(3):
        value = diagonal_candidate[i, i]

        if abs(value) > 0:
            phases[i] = np.exp(-0.5j * np.angle(value))

    U = U @ np.diag(phases)
    diagonalized = U.T @ M @ U

    masses = np.abs(np.diag(diagonalized))

    order = np.argsort(masses)
    masses = masses[order]
    U = U[:, order]

    diagonalized = U.T @ M @ U
    target = np.diag(masses)

    scale = max(
        np.linalg.norm(M),
        1.0e-30,
    )
    residual = float(
        np.linalg.norm(diagonalized - target) / scale
    )

    return U, masses, residual


def calculate_neutrino_observables(
    mass_matrix_gev: np.ndarray,
) -> NeutrinoObservables:
    """Return masses, mass splittings and |PMNS| from m_nu."""

    U, masses_gev, residual = takagi_factorization(
        mass_matrix_gev
    )

    masses_ev = masses_gev * GEV_TO_EV

    delta_m21_sq = float(
        masses_ev[1] ** 2 - masses_ev[0] ** 2
    )
    delta_m31_sq = float(
        masses_ev[2] ** 2 - masses_ev[0] ** 2
    )

    return NeutrinoObservables(
        masses_ev=masses_ev,
        delta_m21_sq_ev2=delta_m21_sq,
        delta_m31_sq_ev2=delta_m31_sq,
        pmns_abs=np.abs(U),
        takagi_residual=residual,
    )


def _complex_matrix_from_file(path: Path) -> np.ndarray:
    payload = json.loads(
        Path(path).read_text(encoding="utf-8")
    )

    real = np.asarray(payload["real"], dtype=float)
    imag = np.asarray(payload["imag"], dtype=float)

    return real + 1j * imag


def run_neutrino_observables_stage(
    mass_matrix_path: Path,
    output_dir: Path,
) -> dict:
    """Read the low-scale mass matrix and save neutrino observables."""

    mass_matrix = _complex_matrix_from_file(
        mass_matrix_path
    )

    observables = calculate_neutrino_observables(
        mass_matrix
    )

    if observables.takagi_residual > 1.0e-7:
        raise RuntimeError(
            "Takagi decomposition residual is too large: "
            f"{observables.takagi_residual}"
        )

    output_dir = Path(output_dir)
    output_path = output_dir / "neutrino_observables.json"

    payload = {
        "NeutrinoObservableStatus": "Success",
        "MassesEV": observables.masses_ev.tolist(),
        "DeltaM21SqEV2": observables.delta_m21_sq_ev2,
        "DeltaM31SqEV2": observables.delta_m31_sq_ev2,
        "PMNSAbs": observables.pmns_abs.tolist(),
        "TakagiResidual": observables.takagi_residual,
    }

    output_path.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    return {
        "NeutrinoObservableStatus": "Success",
        "NeutrinoObservablesFile": output_path.name,
        **payload,
    }
