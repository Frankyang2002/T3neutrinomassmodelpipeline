from __future__ import annotations

import argparse
from pathlib import Path

from tests.python.Weinberg.T3YukawaAdapter import (
    T3WeylConvention,
    real_yukawa_tensor_from_exchange,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inspect the y_ija tensor generated from a T3 RGE exchange JSON."
    )
    parser.add_argument("exchange_json", type=Path)
    args = parser.parse_args()

    import json

    data = json.loads(args.exchange_json.read_text(encoding="utf-8"))

    scalar_model, fermion_basis, tensor = real_yukawa_tensor_from_exchange(
        data,
        T3WeylConvention(
            lepton_multiplicity=1,
            heavy_fermion_multiplicity=1,
            symmetrize_fermion_indices=True,
        ),
    )

    print("Scalar real dimension:", scalar_model.total_real_scalar_dimension)
    print("Fermion Weyl dimension:", fermion_basis.dimension)
    print("Nonzero y_ija components:", len(tensor))

    print()
    print("Fermion blocks:")
    for block in fermion_basis.blocks:
        print(
            f"  {block.fermion.name}[copy={block.copy}] "
            f"-> {block.first}..{block.last}"
        )

    print()
    print("First nonzero y_ija components:")

    for key in sorted(tensor)[:40]:
        print(f"  y{key} = {tensor[key]}")

    if not tensor:
        raise RuntimeError("No y_ija components were generated.")

    print()
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
