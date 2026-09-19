from __future__ import annotations

"""Sparse Wilson-coefficient lookup helpers used by production RGE stages."""

from typing import Mapping

import sympy as sp


class SparseWilsonLookup:
    """Sparse C_ijab lookup with zero for components absent from the mapping."""

    def __init__(
        self,
        components: Mapping[tuple[int, int, int, int], sp.Expr],
    ):
        self.components = components

    def __getitem__(
        self,
        key: tuple[int, int, int, int],
    ) -> sp.Expr:
        return sp.sympify(
            self.components.get(key, sp.S.Zero)
        )

    def __call__(
        self,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> sp.Expr:
        return self[i, j, a, b]


def wilson_component_function(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
) -> SparseWilsonLookup:
    """Return a sparse Wilson-coefficient lookup callable."""

    return SparseWilsonLookup(components)
