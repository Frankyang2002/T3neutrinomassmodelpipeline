"""Tests for the component-support audit's conservative counting boundary."""
import pytest
import sympy as sp
from Numerical.diagnostics.AuditT3ComponentMultiplicity import expand_tensor


def test_singlet_axes_restored():
    t = sp.MutableDenseNDimArray([1, 0, sp.sqrt(2), 0], (2, 2))
    rows = list(expand_tensor(t, (2, 1, 2)))
    assert rows == [((1, 1, 1), 1), ((2, 1, 1), sp.sqrt(2))]


def test_mismatched_tensor_dimensions_fail():
    t = sp.MutableDenseNDimArray([1, 2], (2,))
    with pytest.raises(ValueError, match="do not agree"):
        list(expand_tensor(t, (2, 2, 1)))
