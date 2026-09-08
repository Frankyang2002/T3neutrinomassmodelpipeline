from __future__ import annotations

"""Consolidated neutrino-observable, ordering, Majorana-phase and NuFIT tests."""

import numpy as np

from RGE.phenomenology.NeutrinoDataComparison import (
    compare_to_nufit6,
    compare_to_nufit6_no,
    mixing_angles_from_pmns_abs,
)
from RGE.phenomenology.NeutrinoObservables import (
    calculate_neutrino_observables,
    takagi_factorization,
)
from RGE.phenomenology.T3NeutrinoTarget import build_neutrino_target

EV_TO_GEV = 1.0e-9


def observable_payload(target, ordering: str) -> dict:
    observables = calculate_neutrino_observables(
        target.mass_matrix_ev * EV_TO_GEV,
        ordering=ordering,
    )
    return {
        "PMNSAbs": observables.pmns_abs.tolist(),
        "DeltaM21SqEV2": observables.delta_m21_sq_ev2,
        "DeltaM31SqEV2": observables.delta_m31_sq_ev2,
        "DeltaM32SqEV2": observables.delta_m32_sq_ev2,
    }


def check_takagi_observables() -> None:
    masses_expected = np.array([0.01, 0.012, 0.05])
    masses_gev = masses_expected * EV_TO_GEV

    theta12, theta23, theta13 = 0.59, 0.84, 0.15
    c12, s12 = np.cos(theta12), np.sin(theta12)
    c23, s23 = np.cos(theta23), np.sin(theta23)
    c13, s13 = np.cos(theta13), np.sin(theta13)

    U = np.array(
        [
            [c12*c13, s12*c13, s13],
            [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
            [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13],
        ],
        dtype=complex,
    )
    U = U @ np.diag([1.0, np.exp(0.3j), np.exp(-0.4j)])

    M = U.conj() @ np.diag(masses_gev) @ U.conj().T
    U_fit, masses_fit, residual = takagi_factorization(M)
    obs = calculate_neutrino_observables(M)

    assert np.allclose(obs.masses_ev, masses_expected, rtol=1e-8, atol=1e-12)
    assert residual <= 1e-8
    assert np.allclose(np.abs(U_fit), np.abs(U), rtol=1e-7, atol=1e-9)


def check_ordering(ordering: str) -> None:
    target = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=0.01,
    )
    observables = calculate_neutrino_observables(
        target.mass_matrix_ev * EV_TO_GEV,
        ordering=ordering,
    )

    result = compare_to_nufit6(
        {
            "PMNSAbs": observables.pmns_abs.tolist(),
            "DeltaM21SqEV2": observables.delta_m21_sq_ev2,
            "DeltaM31SqEV2": observables.delta_m31_sq_ev2,
            "DeltaM32SqEV2": observables.delta_m32_sq_ev2,
        },
        ordering=ordering,
    )

    assert result["AllFiveInside3Sigma"] is True
    assert result["Chi2Diagnostic"] < 1e-12

    if ordering == "NO":
        assert observables.masses_ev[0] < observables.masses_ev[1] < observables.masses_ev[2]
        assert observables.delta_m31_sq_ev2 > 0
    else:
        assert observables.masses_ev[2] < observables.masses_ev[0] < observables.masses_ev[1]
        assert observables.delta_m32_sq_ev2 < 0


def check_majorana_phases(ordering: str) -> None:
    zero = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=0.01,
        alpha21=0.0,
        alpha31=0.0,
    )
    phased = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=0.01,
        alpha21=0.7,
        alpha31=1.3,
    )

    assert not np.allclose(
        zero.mass_matrix_ev,
        phased.mass_matrix_ev,
        rtol=1e-12,
        atol=1e-14,
    )

    zero_result = compare_to_nufit6(observable_payload(zero, ordering), ordering=ordering)
    phased_result = compare_to_nufit6(observable_payload(phased, ordering), ordering=ordering)

    assert zero_result["Chi2Diagnostic"] < 1e-10
    assert phased_result["Chi2Diagnostic"] < 1e-10
    assert phased_result["AllFiveInside3Sigma"] is True


def check_angle_extraction() -> None:
    target_s12, target_s23, target_s13 = 0.308, 0.470, 0.02215
    t12 = np.arcsin(np.sqrt(target_s12))
    t23 = np.arcsin(np.sqrt(target_s23))
    t13 = np.arcsin(np.sqrt(target_s13))

    s12, c12 = np.sin(t12), np.cos(t12)
    s23, c23 = np.sin(t23), np.cos(t23)
    s13, c13 = np.sin(t13), np.cos(t13)

    U = np.abs(
        np.array(
            [
                [c12*c13, s12*c13, s13],
                [-s12*c23-c12*s23*s13, c12*c23-s12*s23*s13, s23*c13],
                [s12*s23-c12*c23*s13, -c12*s23-s12*c23*s13, c23*c13],
            ]
        )
    )

    extracted = mixing_angles_from_pmns_abs(U)
    assert abs(extracted["sin2_theta12"] - target_s12) < 1e-12
    assert abs(extracted["sin2_theta23"] - target_s23) < 1e-12
    assert abs(extracted["sin2_theta13"] - target_s13) < 1e-12

    result = compare_to_nufit6_no(
        {
            "PMNSAbs": U.tolist(),
            "DeltaM21SqEV2": 7.49e-5,
            "DeltaM31SqEV2": 2.513e-3,
        }
    )
    assert result["ComparisonStatus"] == "Success"
    assert result["AllFiveInside3Sigma"] is True
    assert result["Chi2Diagnostic"] < 1e-20


def main() -> None:
    check_takagi_observables()
    check_angle_extraction()

    for ordering in ("NO", "IO"):
        check_ordering(ordering)
        check_majorana_phases(ordering)

    print("PASS: Takagi decomposition and observables")
    print("PASS: PMNS/NuFIT comparison")
    print("PASS: normal and inverted ordering")
    print("PASS: Majorana phases")


if __name__ == "__main__":
    main()
