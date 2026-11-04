"""Regression tests for the saved 16-model fermion-threshold scan fixes.

No physics integrations, Matchete calls, Yukawa fits, or saved output writes.
The existing study-contract tests remain unchanged.
"""
from __future__ import annotations

import ast
import inspect
import math
import sys
import types

import pytest

from scripts import RunT3GeneralFermionThresholdScan as study


def _calls(function) -> set[str]:
    """Ensure the repaired helpers are actually used by production branches."""
    tree = ast.parse(inspect.getsource(function))
    return {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }


def test_nonmixed_branch_uses_mass_basis_helper():
    assert "_nonmixed_heavy_mass_basis" in _calls(study._nonmixed)


def test_nonmixed_mass_basis_passes_all_required_tolerances(monkeypatch):
    called = []
    state = object()
    result = object()

    def required_signature(
        arg, *, offdiagonal_atol, imaginary_atol,
        takagi_residual_rtol, takagi_residual_atol,
    ):
        called.append((arg, offdiagonal_atol, imaginary_atol,
                       takagi_residual_rtol, takagi_residual_atol))
        return result

    # The test does not import the full numerical pipeline. The strict stub
    # reproduces the keyword-only function signature in FinalC5TrajectoryAdapter.
    base = types.ModuleType("Numerical")
    base.__path__ = []
    running = types.ModuleType("Numerical.running")
    running.__path__ = []
    adapter = types.ModuleType("Numerical.running.FinalC5TrajectoryAdapter")
    adapter._heavy_mass_basis_data = required_signature
    monkeypatch.setitem(sys.modules, "Numerical", base)
    monkeypatch.setitem(sys.modules, "Numerical.running", running)
    monkeypatch.setitem(sys.modules, adapter.__name__, adapter)

    assert study._nonmixed_heavy_mass_basis(state) is result
    assert called == [(state, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10)]


def test_mixed_branch_uses_provenance_aware_fractions():
    assert "_mixed_fitted_scalar_fractions" in _calls(study._mixed)


def test_original_b_minus_one_direct_fractions():
    original = {
        "extra_quartic_fraction": 0.1,
        "extra_mixed_scalar_bilinear_fraction": 0.01,
    }
    assert study._mixed_fitted_scalar_fractions(original) == (0.1, 0.01)


def test_continued_c_minus_one_provenance_fractions():
    continued = {
        "matching_and_rge_provenance": {
            "quartic_fraction": 0.1,
            "mu2_fraction": 0.01,
            "modified_scalar_boundary": False,
        }
    }
    assert study._mixed_fitted_scalar_fractions(continued) == (0.1, 0.01)


def test_consistent_duplicate_fraction_storage_accepted():
    both = {
        "extra_quartic_fraction": 0.1,
        "extra_mixed_scalar_bilinear_fraction": 0.01,
        "matching_and_rge_provenance": {
            "quartic_fraction": 0.1,
            "mu2_fraction": 0.01,
            "modified_scalar_boundary": False,
        },
    }
    assert study._mixed_fitted_scalar_fractions(both) == (0.1, 0.01)


def test_conflicting_saved_boundary_is_rejected():
    conflicting = {
        "extra_quartic_fraction": 0.1,
        "extra_mixed_scalar_bilinear_fraction": 0.01,
        "matching_and_rge_provenance": {
            "quartic_fraction": 0.2,
            "mu2_fraction": 0.01,
        },
    }
    with pytest.raises(ValueError, match="Conflicting mixed-scalar boundary"):
        study._mixed_fitted_scalar_fractions(conflicting)


def test_changed_continuation_boundary_is_rejected():
    changed = {
        "matching_and_rge_provenance": {
            "quartic_fraction": 0.1,
            "mu2_fraction": 0.01,
            "modified_scalar_boundary": True,
        },
    }
    with pytest.raises(ValueError, match="Continuation changed"):
        study._mixed_fitted_scalar_fractions(changed)


@pytest.mark.parametrize("bad", [
    {},
    {"matching_and_rge_provenance": "invalid"},
    {"extra_quartic_fraction": float("nan"),
     "extra_mixed_scalar_bilinear_fraction": 0.01},
    {"extra_quartic_fraction": 0.0,
     "extra_mixed_scalar_bilinear_fraction": 0.01},
    {"extra_quartic_fraction": 1.01,
     "extra_mixed_scalar_bilinear_fraction": 0.01},
    {"extra_quartic_fraction": 0.1,
     "extra_mixed_scalar_bilinear_fraction": 0.21},
])
def test_invalid_or_missing_fractions_fail_closed(bad):
    with pytest.raises((ValueError, KeyError, TypeError)):
        study._mixed_fitted_scalar_fractions(bad)


def test_scale_sequence_and_mass_scan_unchanged():
    assert study.FACTORS == (0.5, 1.0, 2.0)
    assert study.MIXED == frozenset({"B_-1", "C_-1"})
    assert len(study.MODELS) == 16
    assert study.PINNED_HEAD == "1b4258b19c8da91afe7f352c2df1b1665de61f26"
