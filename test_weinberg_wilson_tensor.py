from __future__ import annotations

import sympy as sp

from GeneralWeinbergRGEGenerator import RGEModel
from T3RGETensors import FermionBasis, WeylFermion
from WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
    validate_weinberg_tensor_symmetry,
)


def main() -> int:
    kappa = sp.Symbol("kappa")

    model = RGEModel.t3(
        d_s1=2,
        y_s1=-sp.Rational(1, 2),
        d_s2=2,
        y_s2=sp.Rational(1, 2),
    )

    fermions = FermionBasis(
        fermions=(
            WeylFermion(
                name="L",
                su2_dimension=2,
                hypercharge=-sp.Rational(1, 2),
            ),
            WeylFermion(
                name="F",
                su2_dimension=1,
                hypercharge=0,
            ),
        )
    )

    tensor = build_weinberg_wilson_tensor(model, fermions, kappa)
    validate_weinberg_tensor_symmetry(tensor)

    print("Nonzero C_ijab components:", len(tensor))
    print()

    for key in sorted(tensor):
        print(f"  C{key} = {sp.simplify(tensor[key])}")

    print()
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
