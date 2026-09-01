from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


# NuFIT 6.0 (2024), normal ordering, including SK + IceCube atmospheric data.
# Source: https://www.nu-fit.org/?q=node/294
#
# We use the quoted one-sigma best-fit uncertainties for a simple diagnostic
# chi-square. This is NOT a replacement for the full correlated NuFIT
# likelihood; it is intended as a pipeline validation / first-pass fit target.
@dataclass(frozen=True)
class Measurement:
    central: float
    sigma_minus: float
    sigma_plus: float
    range_3sigma: tuple[float, float]


NUFIT6_NO = {
    "sin2_theta12": Measurement(
        central=0.308,
        sigma_minus=0.011,
        sigma_plus=0.012,
        range_3sigma=(0.275, 0.345),
    ),
    "sin2_theta23": Measurement(
        central=0.470,
        sigma_minus=0.013,
        sigma_plus=0.017,
        range_3sigma=(0.435, 0.585),
    ),
    "sin2_theta13": Measurement(
        central=0.02215,
        sigma_minus=0.00058,
        sigma_plus=0.00056,
        range_3sigma=(0.02030, 0.02388),
    ),
    "delta_m21_sq_ev2": Measurement(
        central=7.49e-5,
        sigma_minus=0.19e-5,
        sigma_plus=0.19e-5,
        range_3sigma=(6.92e-5, 8.05e-5),
    ),
    "delta_m31_sq_ev2": Measurement(
        central=2.513e-3,
        sigma_minus=0.019e-3,
        sigma_plus=0.021e-3,
        range_3sigma=(2.451e-3, 2.578e-3),
    ),
}


def mixing_angles_from_pmns_abs(pmns_abs: np.ndarray) -> dict[str, float]:
    """Extract the three standard mixing angles from |U_PMNS|.

    This assumes the charged-lepton mass basis and the standard three-neutrino
    parameterisation. The three angles can be obtained from magnitudes alone;
    delta_CP cannot.
    """
    u = np.asarray(pmns_abs, dtype=float)
    if u.shape != (3, 3):
        raise ValueError(f"PMNSAbs must be 3x3, got {u.shape}")

    ue3_sq = float(u[0, 2] ** 2)
    denom = 1.0 - ue3_sq
    if denom <= 0:
        raise ValueError("Invalid PMNS matrix: 1 - |U_e3|^2 <= 0")

    sin2_theta13 = ue3_sq
    sin2_theta12 = float(u[0, 1] ** 2 / denom)
    sin2_theta23 = float(u[1, 2] ** 2 / denom)

    return {
        "sin2_theta12": sin2_theta12,
        "sin2_theta23": sin2_theta23,
        "sin2_theta13": sin2_theta13,
        "theta12_deg": math.degrees(math.asin(math.sqrt(sin2_theta12))),
        "theta23_deg": math.degrees(math.asin(math.sqrt(sin2_theta23))),
        "theta13_deg": math.degrees(math.asin(math.sqrt(sin2_theta13))),
    }


def asymmetric_pull(value: float, measurement: Measurement) -> float:
    sigma = (
        measurement.sigma_plus
        if value >= measurement.central
        else measurement.sigma_minus
    )
    return (value - measurement.central) / sigma


def compare_to_nufit6_no(observables: dict[str, Any]) -> dict[str, Any]:
    """Compare a normal-ordering observable point with NuFIT 6.0.

    The input is the JSON structure written by NeutrinoObservables.py.
    """
    pmns_abs = np.asarray(observables["PMNSAbs"], dtype=float)
    angles = mixing_angles_from_pmns_abs(pmns_abs)

    predicted = {
        "sin2_theta12": angles["sin2_theta12"],
        "sin2_theta23": angles["sin2_theta23"],
        "sin2_theta13": angles["sin2_theta13"],
        "delta_m21_sq_ev2": float(observables["DeltaM21SqEV2"]),
        "delta_m31_sq_ev2": float(observables["DeltaM31SqEV2"]),
    }

    rows: dict[str, Any] = {}
    chi2 = 0.0

    for name, measurement in NUFIT6_NO.items():
        value = predicted[name]
        pull = asymmetric_pull(value, measurement)
        contribution = pull * pull
        chi2 += contribution

        low, high = measurement.range_3sigma
        rows[name] = {
            "Predicted": value,
            "BestFit": measurement.central,
            "SigmaMinus": measurement.sigma_minus,
            "SigmaPlus": measurement.sigma_plus,
            "Pull": pull,
            "Chi2Contribution": contribution,
            "Inside3Sigma": bool(low <= value <= high),
            "Range3Sigma": [low, high],
        }

    all_inside = all(row["Inside3Sigma"] for row in rows.values())

    return {
        "ComparisonStatus": "Success",
        "Ordering": "NO",
        "Dataset": "NuFIT 6.0 (2024), SK+IceCube atmospheric included",
        "Warning": (
            "Diagnostic independent-Gaussian chi-square only; "
            "NuFIT parameter correlations are not included."
        ),
        "PredictedMixingAngles": angles,
        "Observables": rows,
        "Chi2Diagnostic": chi2,
        "AllFiveInside3Sigma": all_inside,
        "DeltaCPCompared": False,
    }


def compare_file(
    observables_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    data = json.loads(observables_path.read_text(encoding="utf-8"))
    result = compare_to_nufit6_no(data)

    if output_path is None:
        output_path = observables_path.with_name("neutrino_data_comparison.json")

    output_path.write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare neutrino observables against NuFIT 6.0 normal-ordering data."
    )
    parser.add_argument(
        "observables",
        type=Path,
        help="Path to neutrino_observables.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output JSON path",
    )
    args = parser.parse_args()

    result = compare_file(args.observables, args.output)

    print("=" * 72)
    print("NEUTRINO DATA COMPARISON — NuFIT 6.0, NORMAL ORDERING")
    print("=" * 72)

    for name, row in result["Observables"].items():
        print(
            f"{name}: predicted={row['Predicted']:.8g}, "
            f"best-fit={row['BestFit']:.8g}, "
            f"pull={row['Pull']:.3g} sigma, "
            f"inside 3sigma={row['Inside3Sigma']}"
        )

    print()
    print(f"Diagnostic chi2 = {result['Chi2Diagnostic']:.8g}")
    print(f"All five inside 3sigma = {result['AllFiveInside3Sigma']}")
    print("delta_CP comparison = not yet available from |U_PMNS| alone")


if __name__ == "__main__":
    main()
