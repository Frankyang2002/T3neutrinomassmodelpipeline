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
    delta_m32_sq_ev2: float = 0.0
    ordering: str = "NO"


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


def _label_mass_eigenstates(
    U: np.ndarray,
    masses: np.ndarray,
    ordering: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Relabel Takagi states into the conventional (m1, m2, m3) order."""

    ordering = ordering.upper()

    if ordering == "NO":
        return U, masses

    if ordering == "IO":
        # Takagi returns ascending masses [m3, m1, m2] for IO.
        # Relabel both masses and PMNS columns to [m1, m2, m3].
        relabel = np.array([1, 2, 0], dtype=int)
        return U[:, relabel], masses[relabel]

    raise ValueError("ordering must be 'NO' or 'IO'")


def calculate_neutrino_observables(
    mass_matrix_gev: np.ndarray,
    ordering: str = "AUTO",
) -> NeutrinoObservables:
    """Return masses, mass splittings and |PMNS| from m_nu."""

    U, masses_gev, residual = takagi_factorization(
        mass_matrix_gev
    )

    if ordering.upper() == "AUTO":
        low_gap = masses_gev[1] ** 2 - masses_gev[0] ** 2
        high_gap = masses_gev[2] ** 2 - masses_gev[1] ** 2
        ordering = "NO" if low_gap < high_gap else "IO"

    U, masses_gev = _label_mass_eigenstates(
        U,
        masses_gev,
        ordering,
    )

    masses_ev = masses_gev * GEV_TO_EV

    delta_m21_sq = float(
        masses_ev[1] ** 2 - masses_ev[0] ** 2
    )
    delta_m31_sq = float(
        masses_ev[2] ** 2 - masses_ev[0] ** 2
    )
    delta_m32_sq = float(
        masses_ev[2] ** 2 - masses_ev[1] ** 2
    )

    return NeutrinoObservables(
        masses_ev=masses_ev,
        delta_m21_sq_ev2=delta_m21_sq,
        delta_m31_sq_ev2=delta_m31_sq,
        pmns_abs=np.abs(U),
        takagi_residual=residual,
        delta_m32_sq_ev2=delta_m32_sq,
        ordering=ordering.upper(),
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
    ordering: str = "AUTO",
) -> dict:
    """Read the low-scale mass matrix and save neutrino observables."""

    mass_matrix = _complex_matrix_from_file(
        mass_matrix_path
    )

    observables = calculate_neutrino_observables(
        mass_matrix,
        ordering=ordering,
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
        "Ordering": observables.ordering,
        "MassesEV": observables.masses_ev.tolist(),
        "DeltaM21SqEV2": observables.delta_m21_sq_ev2,
        "DeltaM31SqEV2": observables.delta_m31_sq_ev2,
        "DeltaM32SqEV2": observables.delta_m32_sq_ev2,
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
