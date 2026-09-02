from __future__ import annotations

import numpy as np

from RGE.phenomenology.NeutrinoObservables import (
    calculate_neutrino_observables,
    takagi_factorization,
)


def main() -> int:
    masses_ev_expected = np.array([
        0.01,
        0.012,
        0.05,
    ])

    masses_gev = masses_ev_expected / 1.0e9

    theta12 = 0.59
    theta23 = 0.84
    theta13 = 0.15

    c12, s12 = np.cos(theta12), np.sin(theta12)
    c23, s23 = np.cos(theta23), np.sin(theta23)
    c13, s13 = np.cos(theta13), np.sin(theta13)

    U = np.array([
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
    ], dtype=complex)

    majorana_phases = np.diag([
        1.0,
        np.exp(0.3j),
        np.exp(-0.4j),
    ])

    U = U @ majorana_phases

    M = U.conj() @ np.diag(masses_gev) @ U.conj().T

    U_fit, masses_fit_gev, residual = takagi_factorization(M)
    obs = calculate_neutrino_observables(M)

    print("=" * 72)
    print("NEUTRINO OBSERVABLE REGRESSION")
    print("=" * 72)
    print()
    print("Recovered masses [eV]:")
    print(obs.masses_ev)
    print()
    print("Expected masses [eV]:")
    print(masses_ev_expected)
    print()
    print("Delta m21^2 [eV^2]:")
    print(obs.delta_m21_sq_ev2)
    print()
    print("Delta m31^2 [eV^2]:")
    print(obs.delta_m31_sq_ev2)
    print()
    print("|PMNS|:")
    print(obs.pmns_abs)
    print()
    print("Takagi residual:")
    print(residual)
    print()

    if not np.allclose(
        obs.masses_ev,
        masses_ev_expected,
        rtol=1e-8,
        atol=1e-12,
    ):
        print("FAIL: masses")
        return 1

    if residual > 1e-8:
        print("FAIL: Takagi residual")
        return 1

    if not np.allclose(
        np.abs(U_fit),
        np.abs(U),
        rtol=1e-7,
        atol=1e-9,
    ):
        print("FAIL: mixing matrix")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
