from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from NeutrinoDataComparison import NUFIT6_NO


EV_TO_GEV = 1.0e-9


@dataclass(frozen=True)
class NeutrinoTarget:
    masses_ev: np.ndarray
    pmns: np.ndarray
    mass_matrix_ev: np.ndarray
    c5_matrix_gev_inv: np.ndarray


def pmns_matrix(
    sin2_theta12: float,
    sin2_theta23: float,
    sin2_theta13: float,
    delta_cp: float = 0.0,
    alpha21: float = 0.0,
    alpha31: float = 0.0,
) -> np.ndarray:
    """Return the standard three-neutrino PMNS matrix.

    Phases are in radians. Majorana phases multiply columns 2 and 3 by
    exp(i alpha21/2) and exp(i alpha31/2).
    """
    s12 = math.sqrt(sin2_theta12)
    s23 = math.sqrt(sin2_theta23)
    s13 = math.sqrt(sin2_theta13)
    c12 = math.sqrt(1.0 - sin2_theta12)
    c23 = math.sqrt(1.0 - sin2_theta23)
    c13 = math.sqrt(1.0 - sin2_theta13)

    eid = np.exp(1j * delta_cp)
    emid = np.exp(-1j * delta_cp)

    u = np.array(
        [
            [c12 * c13, s12 * c13, s13 * emid],
            [
                -s12 * c23 - c12 * s23 * s13 * eid,
                c12 * c23 - s12 * s23 * s13 * eid,
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13 * eid,
                -c12 * s23 - s12 * c23 * s13 * eid,
                c23 * c13,
            ],
        ],
        dtype=complex,
    )

    majorana = np.diag(
        [
            1.0 + 0.0j,
            np.exp(0.5j * alpha21),
            np.exp(0.5j * alpha31),
        ]
    )
    return u @ majorana


def build_normal_ordering_target(
    lightest_mass_ev: float = 0.01,
    vev_gev: float = 246.22,
    delta_cp: float = 0.0,
    alpha21: float = 0.0,
    alpha31: float = 0.0,
) -> NeutrinoTarget:
    """Build a low-scale normal-ordering target from NuFIT 6.0 best fits.

    The present target fixes delta_CP and both Majorana phases explicitly
    rather than fitting them. By default they are set to zero. The five
    observables already validated by NeutrinoDataComparison are therefore
    reproduced exactly at the target point.
    """
    if lightest_mass_ev < 0:
        raise ValueError("lightest_mass_ev must be non-negative")

    dm21 = NUFIT6_NO["delta_m21_sq_ev2"].central
    dm31 = NUFIT6_NO["delta_m31_sq_ev2"].central

    m1 = lightest_mass_ev
    m2 = math.sqrt(m1 * m1 + dm21)
    m3 = math.sqrt(m1 * m1 + dm31)
    masses = np.array([m1, m2, m3], dtype=float)

    u = pmns_matrix(
        NUFIT6_NO["sin2_theta12"].central,
        NUFIT6_NO["sin2_theta23"].central,
        NUFIT6_NO["sin2_theta13"].central,
        delta_cp=delta_cp,
        alpha21=alpha21,
        alpha31=alpha31,
    )

    # Majorana Takagi convention:
    #   U^T M_nu U = diag(m_i)
    # therefore
    #   M_nu = U^* diag(m_i) U^\dagger.
    mass_matrix_ev = u.conj() @ np.diag(masses) @ u.conj().T

    # The matched-operator normalization was directly checked in the
    # T3-B Matchete output: |M_nu| = v^2 |C5|.
    c5_matrix = mass_matrix_ev * EV_TO_GEV / (vev_gev * vev_gev)

    return NeutrinoTarget(
        masses_ev=masses,
        pmns=u,
        mass_matrix_ev=mass_matrix_ev,
        c5_matrix_gev_inv=c5_matrix,
    )


def _complex_matrix_to_json(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    return [
        [
            {"re": float(value.real), "im": float(value.imag)}
            for value in row
        ]
        for row in np.asarray(matrix, dtype=complex)
    ]


def target_to_json(target: NeutrinoTarget) -> dict:
    return {
        "Ordering": "NO",
        "MassesEV": target.masses_ev.tolist(),
        "PMNS": _complex_matrix_to_json(target.pmns),
        "MassMatrixEV": _complex_matrix_to_json(target.mass_matrix_ev),
        "C5MatrixGeVInv": _complex_matrix_to_json(target.c5_matrix_gev_inv),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a NuFIT 6.0 low-scale neutrino target matrix."
    )
    parser.add_argument(
        "--m-lightest",
        type=float,
        default=0.01,
        help="Lightest neutrino mass m1 in eV for normal ordering.",
    )
    parser.add_argument(
        "--vev",
        type=float,
        default=246.22,
        help="Electroweak Higgs VEV in GeV.",
    )
    parser.add_argument(
        "--delta-cp",
        type=float,
        default=0.0,
        help="Dirac CP phase in radians.",
    )
    parser.add_argument(
        "--alpha21",
        type=float,
        default=0.0,
        help="Majorana phase alpha21 in radians.",
    )
    parser.add_argument(
        "--alpha31",
        type=float,
        default=0.0,
        help="Majorana phase alpha31 in radians.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("neutrino_target.json"),
    )
    args = parser.parse_args()

    target = build_normal_ordering_target(
        lightest_mass_ev=args.m_lightest,
        vev_gev=args.vev,
        delta_cp=args.delta_cp,
        alpha21=args.alpha21,
        alpha31=args.alpha31,
    )
    args.output.write_text(
        json.dumps(target_to_json(target), indent=2),
        encoding="utf-8",
    )

    print("=" * 72)
    print("NORMAL-ORDERING NEUTRINO TARGET")
    print("=" * 72)
    print("masses [eV]:", target.masses_ev)
    print()
    print("M_nu [eV]:")
    print(target.mass_matrix_ev)
    print()
    print("C5 [GeV^-1]:")
    print(target.c5_matrix_gev_inv)
    print()
    print("wrote:", args.output)


if __name__ == "__main__":
    main()
