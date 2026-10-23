"""T3-B LLSS C12 diagonal flavour action from the generic one-loop tensor RGE.

Source: RGE/general/WilsonTensorRGE.py and AnomalousDimensions.py.

For the *LLSS/LLHH two-operator projection*, Higgs-only SM Yukawa
interactions yield the lepton-collinear term

    16 pi^2 beta_C12 |_{Ye} = (E @ C12 + C12 @ E.T)/2,
    E = Ye @ Ye.conj().T.

The gauge/quartic coefficient and the LLHH -> LLSS mixing coefficient are
from the already-audited one-generation T3-B tensor projection. This module
does NOT establish closure of the two-operator projection under all one-loop
generated operators or the sequential heavy-threshold matching scheme.
"""
from __future__ import annotations

import argparse
import json
import numpy as np

REVERSE_COEFFICIENT = 8.0 * np.sqrt(3.0) / 3.0


def _matrix(name: str, value: np.ndarray) -> np.ndarray:
    array = np.asarray(value, dtype=np.complex128)
    if array.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 complex matrix")
    if not np.isfinite(array.real).all() or not np.isfinite(array.imag).all():
        raise ValueError(f"{name} must have finite entries")
    return array


def c12_beta_numerator(
    c12: np.ndarray,
    c5_ordered: np.ndarray,
    *,
    ye: np.ndarray,
    g2: complex,
    lambda12: complex,
    lambda12Adj: complex,
    lambdaT3: complex,
    reverse: bool = True,
) -> np.ndarray:
    """Return the 16*pi^2 numerator, in the existing ordered-C5 convention.

    A symmetric C12/C5 input yields a symmetric beta. The first scalar term
    comes from the exact T3-B scalar-quartic basis map; the lepton-flavour
    term is obtained from the fermion collinear anomalous dimension of the
    full master RGE. Reverse mixing is the one-generation scalar projection.
    """
    c12 = _matrix("c12", c12)
    c5 = _matrix("c5_ordered", c5_ordered)
    ye = _matrix("ye", ye)
    for name, value in (("c12", c12), ("c5_ordered", c5)):
        if not np.allclose(value, value.T, rtol=1e-10, atol=1e-24):
            raise ValueError(f"{name} must be symmetric")
    for name, val in (("g2", g2), ("lambda12", lambda12),
                      ("lambda12Adj", lambda12Adj), ("lambdaT3", lambdaT3)):
        if not np.isfinite(complex(val)):
            raise ValueError(f"{name} must be finite")

    e = ye @ ye.conj().T
    scalar = -3.0 * g2**2 + 2.0 * lambda12 - 0.5 * lambda12Adj
    numerator = scalar * c12 + 0.5 * (e @ c12 + c12 @ e.T)
    if reverse:
        numerator = numerator + REVERSE_COEFFICIENT * np.conj(lambdaT3) * c5
    return 0.5 * (numerator + numerator.T)


def verify_lepton_flavour_action_against_master() -> dict:
    """Evaluate a 2-generation complex-Yukawa regression using the *real* master.

    The gauge and quartic inputs are zero so the result isolates the charged-
    lepton Yukawa contribution. The LLSS test tensor has two independent lepton
    flavours and is symmetric in both fermions and scalars. This tests the
    flavour *action*, not gauge/scalar terms or operator-space closure.

    Intended for execution inside the full T3 repository (SymPy installed).
    """
    import sympy as sp
    from RGE.general.FermionBasis import FermionBasis, WeylFermion
    from RGE.general.ScalarBasis import ScalarBasis, ComplexScalar
    from RGE.general.WilsonTensorRGE import (
        WilsonRGEInputs, calculate_complete_wilson_tensor_rge,
    )
    from RGE.matching.WeinbergTensorAdapter import SparseWilsonLookup

    model = ScalarBasis((
        ComplexScalar("H", 2, sp.Rational(1, 2)),
        ComplexScalar("S1", 2, sp.Rational(1, 2)),
        ComplexScalar("S2", 2, -sp.Rational(1, 2)),
    ))
    fb = FermionBasis((
        WeylFermion("L", 2, -sp.Rational(1, 2), multiplicity=2),
        WeylFermion("eC", 1, sp.Integer(1), multiplicity=2),
    ))

    # Complex, non-diagonal and non-symmetric Yukawa matrix: transpose vs
    # Hermitian-conjugate mistakes cannot accidentally pass this regression.
    Y = sp.Matrix([
        [sp.Rational(1, 3) + sp.I / 7, sp.Rational(2, 5) - sp.I / 11],
        [sp.Rational(-1, 4) + sp.I / 13, sp.Rational(3, 8) + sp.I / 9],
    ])
    K = sp.Matrix([
        [sp.Rational(2, 7), sp.Rational(1, 3) + sp.I / 5],
        [sp.Rational(1, 3) + sp.I / 5, sp.Rational(-3, 11)],
    ])
    h = model.block("H")
    hpr, hpi, h0r, h0i = (h.local_to_global(i) for i in range(1, 5))
    nus = [fb.global_index("L", 1, copy=p+1) for p in range(2)]
    electrons = [fb.global_index("L", 2, copy=p+1) for p in range(2)]
    ecs = [fb.global_index("eC", 1, copy=r+1) for r in range(2)]

    yukawas = {}
    for p in range(2):
        for r in range(2):
            for lepton, hr, hi in ((nus[p], hpr, hpi),
                                    (electrons[p], h0r, h0i)):
                for scalar, multiplier in ((hr, 1/sp.sqrt(2)),
                                           (hi, -sp.I/sp.sqrt(2))):
                    value = multiplier * Y[p, r]
                    yukawas[lepton, ecs[r], scalar] = value
                    yukawas[ecs[r], lepton, scalar] = value

    def yukawa(i, j, a):
        return yukawas.get((i, j, a), sp.S.Zero)

    s1r = model.block("S1").local_to_global(1)
    s2r = model.block("S2").local_to_global(1)
    components = {}
    for p in range(2):
        for q in range(2):
            for si, sj in ((s1r, s2r), (s2r, s1r)):
                components[nus[p], nus[q], si, sj] = K[p, q]
    coefficient = SparseWilsonLookup(components)
    inputs = WilsonRGEInputs(
        fermion_dimension=fb.dimension,
        yukawa=yukawa,
        quartic=lambda a, b, c, d: sp.S.Zero,
        gauge_sectors=(),
    )
    E = Y * Y.conjugate().T
    expected = (E*K + K*E.T)/2

    residuals = {}
    other_term_nonzero = {}
    for p, q in ((0, 0), (0, 1), (1, 1)):
        contributions = calculate_complete_wilson_tensor_rge(
            model, inputs, (nus[p], nus[q], s1r, s2r),
            coefficient=coefficient, simplify_each=False,
        )
        label = f"{p+1},{q+1}"
        residual = sp.simplify(contributions["total"] - expected[p, q])
        residuals[label] = str(residual)
        for name, expression in contributions.items():
            if name not in {"total", "fermion_anomalous_dimension"}:
                if sp.simplify(expression) != 0:
                    other_term_nonzero[f"{label}:{name}"] = str(expression)
        if sp.simplify(contributions["fermion_anomalous_dimension"] - expected[p, q]) != 0:
            other_term_nonzero[f"{label}:fermion_anomalous_dimension"] = "incorrect matrix action"
    passed = all(value == "0" for value in residuals.values()) and not other_term_nonzero
    return {
        "status": "C12LeptonFlavorSourceVerified" if passed else "C12LeptonFlavorMismatch",
        "source": "RGE/general/WilsonTensorRGE.py and RGE/general/AnomalousDimensions.py",
        "scope": "two-generation complex non-diagonal Ye, LLSS Yukawa action; gauge/quartic zero",
        "matrix_formula": "(Ye Ye_dagger C12 + C12 (Ye Ye_dagger)^T)/2",
        "master_minus_matrix_residuals": residuals,
        "unexpected_terms": other_term_nonzero,
        "operator_closure_verified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-master", action="store_true", required=True)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    result = verify_lepton_flavour_action_against_master()
    if args.output:
        from pathlib import Path
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if result["status"] != "C12LeptonFlavorSourceVerified":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
