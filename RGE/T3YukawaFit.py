from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from RGE.T3NeutrinoTarget import build_normal_ordering_target


@dataclass(frozen=True)
class YukawaFit:
    y1: np.ndarray
    y2: np.ndarray
    reconstructed_c5: np.ndarray
    relative_residual: float
    max_abs_y1: float
    max_abs_y2: float


def reconstruct_c5(
    y1: np.ndarray,
    y2: np.ndarray,
    loop_factors: Sequence[complex],
) -> np.ndarray:
    """Reconstruct the flavor-lifted T3 C5 matrix.

    Convention:
        C5_pq = 1/2 sum_r F_r [
            y1_pr^* y2_qr^* + y2_pr^* y1_qr^*
        ].
    """
    y1 = np.asarray(y1, dtype=complex)
    y2 = np.asarray(y2, dtype=complex)
    f = np.asarray(loop_factors, dtype=complex)

    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape")
    if y1.shape[0] != 3:
        raise ValueError("This fitter currently expects three lepton flavors")
    if y1.shape[1] != len(f):
        raise ValueError("Number of heavy columns must match loop_factors")

    a = y1.conj()
    b = y2.conj()
    d = np.diag(f)

    return 0.5 * (a @ d @ b.T + b @ d @ a.T)


def fit_three_heavy_balanced(
    target_c5: np.ndarray,
    loop_factors: Sequence[complex],
    singular_tolerance: float = 1e-30,
) -> YukawaFit:
    """Construct an exact three-heavy-generation Yukawa solution.

    We use A = y1^* = s I and choose B = y2^* symmetric. Then

        C5_ij = s/2 (F_i + F_j) B_ij,

    so

        B_ij = 2 C5_ij / [s (F_i + F_j)].

    The real positive scale s is chosen to balance the largest entries of y1
    and y2. This is a constructive T3 analogue of using neutrino data as the
    input rather than randomly scanning Yukawa matrices.

    The construction requires three heavy generations and F_i + F_j != 0.
    """
    k = np.asarray(target_c5, dtype=complex)
    f = np.asarray(loop_factors, dtype=complex)

    if k.shape != (3, 3):
        raise ValueError("target_c5 must be a 3x3 matrix")
    if len(f) != 3:
        raise ValueError("This exact construction currently requires 3 heavy generations")
    if not np.allclose(k, k.T, rtol=1e-10, atol=1e-30):
        raise ValueError("target_c5 must be complex symmetric")

    b0 = np.zeros((3, 3), dtype=complex)

    for i in range(3):
        for j in range(3):
            denominator = f[i] + f[j]
            if abs(denominator) <= singular_tolerance:
                raise ValueError(
                    f"Loop-factor denominator F[{i}] + F[{j}] is singular"
                )
            b0[i, j] = 2.0 * k[i, j] / denominator

    # With s=1, max|y1|=1 and max|y2|=max|B0|.
    # Rescaling y1 -> s y1, y2 -> y2/s leaves C5 unchanged.
    # s=sqrt(max|B0|) equalizes the two maxima.
    max_b0 = float(np.max(np.abs(b0)))
    scale = np.sqrt(max_b0) if max_b0 > 0.0 else 1.0

    a = scale * np.eye(3, dtype=complex)
    b = b0 / scale

    y1 = a.conj()
    y2 = b.conj()

    reconstructed = reconstruct_c5(y1, y2, f)
    norm = np.linalg.norm(k)
    residual = float(
        np.linalg.norm(reconstructed - k) / norm
        if norm > 0
        else np.linalg.norm(reconstructed - k)
    )

    return YukawaFit(
        y1=y1,
        y2=y2,
        reconstructed_c5=reconstructed,
        relative_residual=residual,
        max_abs_y1=float(np.max(np.abs(y1))),
        max_abs_y2=float(np.max(np.abs(y2))),
    )


def _complex_matrix_to_json(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    return [
        [
            {"re": float(value.real), "im": float(value.imag)}
            for value in row
        ]
        for row in np.asarray(matrix, dtype=complex)
    ]


def _parse_complex(value: str) -> complex:
    return complex(value.replace("i", "j"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Construct three-generation T3 Yukawa matrices for a target "
            "normal-ordering neutrino mass matrix."
        )
    )
    parser.add_argument(
        "--loop-factors",
        nargs=3,
        required=True,
        help=(
            "Three loop factors F_r in GeV^-1. Complex values may use j, "
            "for example 1e-8+2e-9j."
        ),
    )
    parser.add_argument(
        "--m-lightest",
        type=float,
        default=0.01,
        help="Lightest normal-ordering neutrino mass in eV.",
    )
    parser.add_argument(
        "--vev",
        type=float,
        default=246.22,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("t3_yukawa_fit.json"),
    )
    args = parser.parse_args()

    loop_factors = [_parse_complex(value) for value in args.loop_factors]
    target = build_normal_ordering_target(
        lightest_mass_ev=args.m_lightest,
        vev_gev=args.vev,
    )
    fit = fit_three_heavy_balanced(
        target.c5_matrix_gev_inv,
        loop_factors,
    )

    payload = {
        "Status": "Success",
        "Ordering": "NO",
        "LightestMassEV": args.m_lightest,
        "LoopFactorsGeVInv": [
            {"re": float(x.real), "im": float(x.imag)}
            for x in loop_factors
        ],
        "Y1": _complex_matrix_to_json(fit.y1),
        "Y2": _complex_matrix_to_json(fit.y2),
        "ReconstructedC5GeVInv": _complex_matrix_to_json(
            fit.reconstructed_c5
        ),
        "RelativeResidual": fit.relative_residual,
        "MaxAbsY1": fit.max_abs_y1,
        "MaxAbsY2": fit.max_abs_y2,
        "Assumptions": [
            "three heavy generations",
            "diagonal heavy-fermion mass basis",
            "fixed loop factors F_r",
            "normal ordering",
            "delta_CP = 0",
            "Majorana phases = 0",
            "|M_nu| = v^2 |C5| normalization",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=" * 72)
    print("T3 CONSTRUCTIVE YUKAWA FIT")
    print("=" * 72)
    print("loop factors [GeV^-1]:", loop_factors)
    print("max |y1| =", fit.max_abs_y1)
    print("max |y2| =", fit.max_abs_y2)
    print("relative C5 residual =", fit.relative_residual)
    print("wrote:", args.output)


if __name__ == "__main__":
    main()
