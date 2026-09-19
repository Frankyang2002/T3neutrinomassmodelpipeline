from __future__ import annotations

"""Shared parser for Matchete Clebsch--Gordan registry entries."""

from dataclasses import dataclass

import sympy as sp

from RGE.running.MatcheteParsing import parse_reps, parse_sparse_array


@dataclass(frozen=True)
class CGTensor:
    name: str
    reps: tuple[str, ...]
    tensor: sp.MutableDenseNDimArray


def load_cg_registry(seed: dict) -> dict[str, CGTensor]:
    """Parse a Matchete ``CGRegistry`` payload into exact tensors."""
    result: dict[str, CGTensor] = {}

    for item in seed.get("CGRegistry", []):
        result[item["Name"]] = CGTensor(
            name=item["Name"],
            reps=parse_reps(item["RepsInputForm"]),
            tensor=parse_sparse_array(item["TensorInputForm"]),
        )

    return result
