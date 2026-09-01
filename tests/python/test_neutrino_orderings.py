from __future__ import annotations

import numpy as np

from RGE.NeutrinoDataComparison import compare_to_nufit6
from RGE.NeutrinoObservables import calculate_neutrino_observables
from RGE.T3NeutrinoTarget import build_neutrino_target


EV_TO_GEV = 1.0e-9


def check_ordering(ordering: str) -> None:
    target = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=0.01,
    )

    observables = calculate_neutrino_observables(
        target.mass_matrix_ev * EV_TO_GEV,
        ordering=ordering,
    )

    payload = {
        "PMNSAbs": observables.pmns_abs.tolist(),
        "DeltaM21SqEV2": observables.delta_m21_sq_ev2,
        "DeltaM31SqEV2": observables.delta_m31_sq_ev2,
        "DeltaM32SqEV2": observables.delta_m32_sq_ev2,
    }

    result = compare_to_nufit6(payload, ordering=ordering)

    assert result["AllFiveInside3Sigma"] is True
    assert result["Chi2Diagnostic"] < 1.0e-12

    if ordering == "NO":
        assert observables.masses_ev[0] < observables.masses_ev[1]
        assert observables.masses_ev[1] < observables.masses_ev[2]
        assert observables.delta_m31_sq_ev2 > 0.0
    else:
        assert observables.masses_ev[2] < observables.masses_ev[0]
        assert observables.masses_ev[0] < observables.masses_ev[1]
        assert observables.delta_m32_sq_ev2 < 0.0

    print(
        f"PASS: {ordering} target -> Takagi relabelling -> "
        f"NuFIT comparison"
    )


def main() -> None:
    check_ordering("NO")
    check_ordering("IO")


if __name__ == "__main__":
    main()
