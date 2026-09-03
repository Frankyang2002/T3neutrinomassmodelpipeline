from __future__ import annotations

"""Compare UV and low-energy RGE effects in the minimal Ma model.

This module performs a common-threshold calculation:

    Ma model at mu_uv -> Ma RGEs -> matching at mu_match
    -> SM + Weinberg operator RGEs -> mu_low.

Four calculations are returned to isolate the effects:

``frozen``
    No running.  Match the input UV parameters and keep C5 fixed.
``uv_only``
    Run the Ma parameters to the threshold and match, but do not run below it.
``eft_only``
    Match the frozen input parameters, then run only the low-energy EFT.
``full``
    Run above and below the common matching threshold.

The matching uses the small-lambda5 expression and the project convention
``m_nu = -v^2 C5``.  Split thresholds and finite threshold corrections are
deliberately outside this first benchmark.
"""

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    from .MaUVRGE import MaUVInitialConditions, evolve_ma_uv, match_scotogenic_c5
    from .NumericalWeinbergRGE import (
        SMInitialConditions,
        evolve_weinberg,
        neutrino_mass_matrix,
    )
except ImportError:  # Allow direct execution from RGE/running.
    from MaUVRGE import MaUVInitialConditions, evolve_ma_uv, match_scotogenic_c5
    from NumericalWeinbergRGE import (
        SMInitialConditions,
        evolve_weinberg,
        neutrino_mass_matrix,
    )


GEV_TO_EV = 1.0e9


def _complex_matrix(value: Any, name: str) -> np.ndarray:
    """Read a 3x3 real matrix or a {real, imag} JSON matrix."""

    if isinstance(value, dict):
        if "real" not in value:
            raise ValueError(f"{name} requires a 'real' matrix.")
        real = np.asarray(value["real"], dtype=float)
        imag = np.asarray(value.get("imag", np.zeros_like(real)), dtype=float)
        matrix = real + 1j * imag
    else:
        matrix = np.asarray(value, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix.")
    return matrix


def state_from_config(data: dict[str, Any]) -> MaUVInitialConditions:
    """Construct validated UV initial conditions from a JSON mapping."""

    required_scalars = (
        "gY",
        "g2",
        "g3",
        "lambda1",
        "lambda2",
        "lambda3",
        "lambda4",
        "lambda5",
        "mH2",
        "mEta2",
    )
    missing = [name for name in required_scalars if name not in data]
    missing += [name for name in ("Ye", "Yu", "Yd", "h", "M") if name not in data]
    if missing:
        raise ValueError("Missing UV inputs: " + ", ".join(missing))

    return MaUVInitialConditions(
        **{name: float(data[name]) for name in required_scalars},
        Ye=_complex_matrix(data["Ye"], "Ye"),
        Yu=_complex_matrix(data["Yu"], "Yu"),
        Yd=_complex_matrix(data["Yd"], "Yd"),
        h=_complex_matrix(data["h"], "h"),
        M=_complex_matrix(data["M"], "M"),
    ).validated()


def _charged_lepton_basis(
    Ye: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return ascending Yukawas, the lepton-doublet rotation, and residual.

    With ``Ye = U_R diag(ye) V_L^dagger``, the field convention is
    ``L_old = V_L L_mass``.  Consequently ``h_mass = h V_L`` and a symmetric
    Weinberg coefficient transforms as ``C5_mass = V_L.T C5 V_L``.
    """

    singlet_rotation, singular_values, left_dagger = np.linalg.svd(Ye)
    left_rotation = left_dagger.conj().T
    order = np.argsort(singular_values)
    singular_values = singular_values[order]
    singlet_rotation = singlet_rotation[:, order]
    left_rotation = left_rotation[:, order]

    # The SVD leaves one arbitrary phase per mass eigenstate.  Fix it so that
    # matrix comparisons between scenarios are not polluted by basis phases.
    for column in range(left_rotation.shape[1]):
        pivot = int(np.argmax(np.abs(left_rotation[:, column])))
        phase = np.angle(left_rotation[pivot, column])
        phase_factor = np.exp(-1j * phase)
        left_rotation[:, column] *= phase_factor
        singlet_rotation[:, column] *= phase_factor
        if left_rotation[pivot, column].real < 0.0:
            left_rotation[:, column] *= -1.0
            singlet_rotation[:, column] *= -1.0

    diagonalized = singlet_rotation.conj().T @ Ye @ left_rotation
    target = np.diag(singular_values)
    scale = max(1.0e-30, float(np.linalg.norm(Ye)))
    residual = float(np.linalg.norm(diagonalized - target) / scale)
    if residual > 1.0e-10:
        raise RuntimeError(
            "Charged-lepton singular-value decomposition failed: "
            f"relative residual={residual:.3e}."
        )

    return singular_values, left_rotation, residual


def _diagonal_yukawa(matrix: np.ndarray, name: str) -> np.ndarray:
    """Extract a diagonal Yukawa vector and reject unsupported mixing."""

    diagonal = np.diag(matrix)
    off_diagonal = matrix - np.diag(diagonal)
    scale = max(1.0e-16, float(np.linalg.norm(matrix)))
    if np.linalg.norm(off_diagonal) > 1.0e-8 * scale:
        raise ValueError(
            f"{name} developed appreciable off-diagonal entries. "
            "The present SMEFT runner assumes diagonal quark Yukawas."
        )
    return np.abs(diagonal)


def _threshold_state(state: MaUVInitialConditions) -> tuple[SMInitialConditions, dict]:
    """Match a Ma state and rotate C5 to the charged-lepton mass basis."""

    match = match_scotogenic_c5(state)
    ye, left_rotation, diagonalization_residual = _charged_lepton_basis(state.Ye)
    C5 = left_rotation.T @ match["C5"] @ left_rotation
    C5 = 0.5 * (C5 + C5.T)

    sm = SMInitialConditions(
        gY=state.gY,
        g2=state.g2,
        g3=state.g3,
        lambdaH=state.lambda1,
        ye=ye,
        yu=_diagonal_yukawa(state.Yu, "Yu"),
        yd=_diagonal_yukawa(state.Yd, "Yd"),
        K=C5,
    ).validated()
    metadata = {
        "heavy_masses_gev": match["heavy_masses"],
        "loop_weights_gev_inverse": match["loop_weights"],
        "charged_lepton_yukawas": ye,
        "charged_lepton_left_rotation": left_rotation,
        "charged_lepton_diagonalization_relative_residual": diagonalization_residual,
    }
    return sm, metadata


def _run_eft(initial: SMInitialConditions, mu_match: float, mu_low: float):
    result = evolve_weinberg(initial, mu_match, mu_low)
    return result.K, result


def _matrix_json(matrix: np.ndarray) -> dict[str, list]:
    matrix = np.asarray(matrix, dtype=complex)
    return {"real": matrix.real.tolist(), "imag": matrix.imag.tolist()}


def _array_json(array: np.ndarray) -> list:
    return np.asarray(array, dtype=float).tolist()


def _state_summary(state: MaUVInitialConditions) -> dict[str, Any]:
    return {
        "gY": state.gY,
        "g2": state.g2,
        "g3": state.g3,
        "lambda1": state.lambda1,
        "lambda2": state.lambda2,
        "lambda3": state.lambda3,
        "lambda4": state.lambda4,
        "lambda5": state.lambda5,
        "mH2_gev2": state.mH2,
        "mEta2_gev2": state.mEta2,
        "M_gev": _matrix_json(state.M),
        "h": _matrix_json(state.h),
    }


def _case_summary(C5: np.ndarray, vev_gev: float) -> dict[str, Any]:
    mass_ev = neutrino_mass_matrix(C5, vev_gev=vev_gev) * GEV_TO_EV
    masses_ev = np.sort(np.linalg.svd(mass_ev, compute_uv=False))
    return {
        "C5_gev_inverse": _matrix_json(C5),
        "neutrino_mass_matrix_ev": _matrix_json(mass_ev),
        "neutrino_masses_ev": _array_json(masses_ev),
        "C5_frobenius_norm_gev_inverse": float(np.linalg.norm(C5)),
        "mass_frobenius_norm_ev": float(np.linalg.norm(mass_ev)),
    }


def _relative_change(new: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    if denominator == 0.0:
        raise ValueError("Cannot compare RGE effects because the frozen C5 is zero.")
    return float(np.linalg.norm(new - reference) / denominator)


def compare_ma_running(config: dict[str, Any]) -> dict[str, Any]:
    """Run and compare the common-threshold UV/EFT scenarios."""

    mu_uv = float(config["mu_uv_gev"])
    mu_match = float(config["mu_matching_gev"])
    mu_low = float(config.get("mu_low_gev", 91.1876))
    vev_gev = float(config.get("vev_gev", 246.22))
    if not (mu_uv >= mu_match >= mu_low > 0.0):
        raise ValueError("Require mu_uv >= mu_matching >= mu_low > 0.")

    high_state = state_from_config(config["uv"])
    frozen_sm, frozen_match = _threshold_state(high_state)

    uv_result = evolve_ma_uv(high_state, mu_uv, mu_match)
    evolved_sm, evolved_match = _threshold_state(uv_result.state)

    C5_frozen = frozen_sm.K
    C5_uv_only = evolved_sm.K
    C5_eft_only, eft_result = _run_eft(frozen_sm, mu_match, mu_low)
    C5_full, full_eft_result = _run_eft(evolved_sm, mu_match, mu_low)

    interaction = C5_full - C5_uv_only - C5_eft_only + C5_frozen
    cases = {
        "frozen": _case_summary(C5_frozen, vev_gev),
        "uv_only": _case_summary(C5_uv_only, vev_gev),
        "eft_only": _case_summary(C5_eft_only, vev_gev),
        "full": _case_summary(C5_full, vev_gev),
    }
    effects = {
        "uv_only_relative_to_frozen": _relative_change(C5_uv_only, C5_frozen),
        "eft_only_relative_to_frozen": _relative_change(C5_eft_only, C5_frozen),
        "full_relative_to_frozen": _relative_change(C5_full, C5_frozen),
        "uv_eft_interaction_relative_to_frozen": float(
            np.linalg.norm(interaction) / np.linalg.norm(C5_frozen)
        ),
        "full_mass_norm_ratio_to_frozen": (
            cases["full"]["mass_frobenius_norm_ev"]
            / cases["frozen"]["mass_frobenius_norm_ev"]
        ),
    }

    return {
        "calculation": "Ma UV -> common-threshold C5 -> SM EFT",
        "fit_target": config.get("fit_target"),
        "conventions": {
            "loop_derivative": "16 pi^2 d/d ln(mu)",
            "neutrino_mass": "m_nu = -v^2 C5",
            "scalar_potential": "lambda1/2 (Hdag H)^2; lambdaT3 = -lambda5",
            "matchete_yukawa_bridge": "y1 = y2 = conjugate(h)",
            "limitations": [
                "small-lambda5 matching formula",
                "one common matching threshold",
                "diagonal quark Yukawas below the threshold",
                "one-loop running",
            ],
        },
        "scales_gev": {"uv": mu_uv, "matching": mu_match, "low": mu_low},
        "uv_running": {
            "solver_success": uv_result.solver_success,
            "solver_message": uv_result.solver_message,
            "function_evaluations": uv_result.nfev,
            "initial": _state_summary(high_state),
            "at_matching_scale": _state_summary(uv_result.state),
        },
        "matching": {
            "frozen_heavy_masses_gev": _array_json(frozen_match["heavy_masses_gev"]),
            "evolved_heavy_masses_gev": _array_json(evolved_match["heavy_masses_gev"]),
            "frozen_loop_weights_gev_inverse": _array_json(
                frozen_match["loop_weights_gev_inverse"]
            ),
            "evolved_loop_weights_gev_inverse": _array_json(
                evolved_match["loop_weights_gev_inverse"]
            ),
            "frozen_charged_lepton_yukawas": _array_json(
                frozen_match["charged_lepton_yukawas"]
            ),
            "evolved_charged_lepton_yukawas": _array_json(
                evolved_match["charged_lepton_yukawas"]
            ),
            "frozen_charged_lepton_left_rotation": _matrix_json(
                frozen_match["charged_lepton_left_rotation"]
            ),
            "evolved_charged_lepton_left_rotation": _matrix_json(
                evolved_match["charged_lepton_left_rotation"]
            ),
            "frozen_charged_lepton_diagonalization_relative_residual": (
                frozen_match["charged_lepton_diagonalization_relative_residual"]
            ),
            "evolved_charged_lepton_diagonalization_relative_residual": (
                evolved_match["charged_lepton_diagonalization_relative_residual"]
            ),
        },
        "eft_running": {
            "frozen_input_function_evaluations": eft_result.nfev,
            "evolved_input_function_evaluations": full_eft_result.nfev,
        },
        "cases": cases,
        "effects": effects,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="Ma benchmark JSON file")
    parser.add_argument("--output", type=Path, required=True, help="output JSON file")
    args = parser.parse_args()

    with args.config.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    result = compare_ma_running(config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    effects = result["effects"]
    print(f"Wrote {args.output}")
    print(
        "Relative C5 shifts: "
        f"UV={effects['uv_only_relative_to_frozen']:.6g}, "
        f"EFT={effects['eft_only_relative_to_frozen']:.6g}, "
        f"full={effects['full_relative_to_frozen']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
