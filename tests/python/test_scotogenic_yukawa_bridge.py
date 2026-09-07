from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from tests.python.Weinberg.T3YukawaAdapter import (
    T3WeylConvention,
    real_yukawa_tensor_from_exchange,
)


PROJECT_ROOT = Path(__file__).resolve().parent
EXCHANGE = (
    PROJECT_ROOT
    / "wolfram"
    / "output"
    / "T3_B_alpha_m1"
    / "rge_tensor_exchange.json"
)


def main() -> int:
    if not EXCHANGE.exists():
        raise FileNotFoundError(
            "Missing exchange JSON. First run:\n"
            "  wolframscript -file "
            r".\wolfram\tests\TestT3RGEExchangeExport.wl"
        )

    data = json.loads(EXCHANGE.read_text(encoding="utf-8"))

    scalar_model, fermion_basis, tensor = real_yukawa_tensor_from_exchange(
        data,
        T3WeylConvention(
            lepton_multiplicity=1,
            heavy_fermion_multiplicity=1,
            symmetrize_fermion_indices=True,
        ),
    )

    if not tensor:
        raise RuntimeError("No y_ija components were produced.")

    print("Model:", data.get("model_name"))
    print("Scalar real dimension:", scalar_model.total_real_scalar_dimension)
    print("Fermion Weyl dimension:", fermion_basis.dimension)
    print("Nonzero y_ija components:", len(tensor))
    print()

    print("Fermion basis:")
    for block in fermion_basis.blocks:
        print(
            f"  {block.fermion.name}[copy={block.copy}] "
            f"= {block.first}..{block.last}"
        )

    print()
    print("Nonzero y_ija:")

    for key in sorted(tensor):
        print(f"  y{key} = {sp.simplify(tensor[key])}")

    print()
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
