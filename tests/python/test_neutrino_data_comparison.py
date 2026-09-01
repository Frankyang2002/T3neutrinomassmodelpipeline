from __future__ import annotations

import numpy as np

from RGE.NeutrinoDataComparison import compare_to_nufit6_no, mixing_angles_from_pmns_abs


def pmns_abs_from_angles(theta12: float, theta23: float, theta13: float) -> np.ndarray:
    """CP-conserving PMNS magnitude matrix for a regression test."""
    s12, c12 = np.sin(theta12), np.cos(theta12)
    s23, c23 = np.sin(theta23), np.cos(theta23)
    s13, c13 = np.sin(theta13), np.cos(theta13)

    u = np.array(
        [
            [c12 * c13, s12 * c13, s13],
            [
                -s12 * c23 - c12 * s23 * s13,
                c12 * c23 - s12 * s23 * s13,
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13,
                -c12 * s23 - s12 * c23 * s13,
                c23 * c13,
            ],
        ],
        dtype=float,
    )
    return np.abs(u)


def main() -> None:
    target_s12 = 0.308
    target_s23 = 0.470
    target_s13 = 0.02215

    theta12 = np.arcsin(np.sqrt(target_s12))
    theta23 = np.arcsin(np.sqrt(target_s23))
    theta13 = np.arcsin(np.sqrt(target_s13))

    pmns_abs = pmns_abs_from_angles(theta12, theta23, theta13)
    extracted = mixing_angles_from_pmns_abs(pmns_abs)

    assert abs(extracted["sin2_theta12"] - target_s12) < 1e-12
    assert abs(extracted["sin2_theta23"] - target_s23) < 1e-12
    assert abs(extracted["sin2_theta13"] - target_s13) < 1e-12

    data = {
        "PMNSAbs": pmns_abs.tolist(),
        "DeltaM21SqEV2": 7.49e-5,
        "DeltaM31SqEV2": 2.513e-3,
    }
    result = compare_to_nufit6_no(data)

    assert result["ComparisonStatus"] == "Success"
    assert result["AllFiveInside3Sigma"] is True
    assert result["Chi2Diagnostic"] < 1e-20

    print("PASS: PMNS angle extraction")
    print("PASS: exact NuFIT best-fit point gives chi2 ~ 0")
    print("PASS: five-observable comparison")


if __name__ == "__main__":
    main()
