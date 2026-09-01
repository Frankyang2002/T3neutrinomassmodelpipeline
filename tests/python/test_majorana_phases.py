from __future__ import annotations

import numpy as np

from RGE.NeutrinoDataComparison import compare_to_nufit6
from RGE.NeutrinoObservables import calculate_neutrino_observables
from RGE.T3NeutrinoTarget import build_neutrino_target


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


def check_majorana_phases(ordering: str) -> None:
    zero_phase = build_neutrino_target(
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

    if np.allclose(
        zero_phase.mass_matrix_ev,
        phased.mass_matrix_ev,
        rtol=1e-12,
        atol=1e-14,
    ):
        raise AssertionError(
            f"{ordering}: Majorana phases did not change M_nu."
        )

    zero_result = compare_to_nufit6(
        observable_payload(zero_phase, ordering),
        ordering=ordering,
    )
    phased_result = compare_to_nufit6(
        observable_payload(phased, ordering),
        ordering=ordering,
    )

    if zero_result["Chi2Diagnostic"] > 1e-10:
        raise AssertionError(
            f"{ordering}: zero-phase target no longer reproduces NuFIT."
        )

    if phased_result["Chi2Diagnostic"] > 1e-10:
        raise AssertionError(
            f"{ordering}: Majorana phases changed oscillation observables."
        )

    if not phased_result["AllFiveInside3Sigma"]:
        raise AssertionError(
            f"{ordering}: phased point failed the oscillation comparison."
        )

    print(
        f"PASS: {ordering} Majorana phases change M_nu "
        f"but not oscillation observables"
    )


def main() -> None:
    check_majorana_phases("NO")
    check_majorana_phases("IO")


if __name__ == "__main__":
    main()
