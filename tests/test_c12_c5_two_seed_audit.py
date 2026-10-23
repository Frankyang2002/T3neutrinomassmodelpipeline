"""Algebra-only tests for the independent C12/C5 audit helpers."""
import sympy as sp
import pytest

from RGE.running.intermediate.AuditC12C5Mixing import (
    canonical, complete_symmetric_entries, mixed_yukawa_part
)


def test_canonical_components():
    assert canonical((2, 1, 6, 5)) == (1, 2, 5, 6)


def test_symmetric_expansion():
    entries = complete_symmetric_entries({(1, 2, 5, 6): sp.Symbol("a")})
    assert len(entries) == 4
    assert entries[(2, 1, 6, 5)] == sp.Symbol("a")


def test_conflicting_symmetry_rejected():
    with pytest.raises(ValueError):
        complete_symmetric_entries({(1, 2, 5, 6): 1, (2, 1, 6, 5): 2})


def test_mixed_yukawa_extraction():
    y1, y2 = sp.symbols("y1 y2")
    e = 3*y1**2 + 5*y2**2 + 7*y1*y2 + 11*sp.conjugate(y1)*y2
    assert sp.simplify(
        mixed_yukawa_part(e) - (7*y1*y2 + 11*sp.conjugate(y1)*y2)
    ) == 0
