"""Regression for public Matchete scalar parser API used by the coupled B/C runner."""
from __future__ import annotations

import inspect
import sys
from types import ModuleType

import numpy as np
import pytest

from Numerical.experimental import T3BCCoupledWilsonTensor as tensor


def test_import_uses_public_parser_name():
    source = inspect.getsource(tensor._matchete_complex)
    assert "import parse_mathematica_scalar" in source
    assert "import _parse_mathematica_scalar" not in source


def test_parse_matchete_complex_with_public_api(monkeypatch):
    # A minimal stand-in exercises the import through the published name even
    # when tests run outside the full public-repository checkout.
    parser = ModuleType("RGE.running.intermediate.IntermediateMatcheteParsing")
    import sympy as sp

    def parse_mathematica_scalar(raw):
        cleaned = raw.replace("Sqrt[", "sqrt(").replace("]", ")").replace("^", "**")
        return sp.sympify(cleaned, locals={"I": sp.I, "sqrt": sp.sqrt})

    parser.parse_mathematica_scalar = parse_mathematica_scalar
    monkeypatch.setitem(sys.modules, parser.__name__, parser)
    assert np.isclose(tensor._matchete_complex("(-I)/Sqrt[2]"), -1j/np.sqrt(2))
    assert np.isclose(tensor._matchete_complex("3^(-1/4)"), 3**(-.25))
    assert tensor._matchete_complex("0") == 0j


def test_real_repository_parser_if_present():
    try:
        from RGE.running.intermediate.IntermediateMatcheteParsing import parse_mathematica_scalar
    except ModuleNotFoundError:
        pytest.skip("Repository parser unavailable in isolated environment")
    assert np.isclose(complex(parse_mathematica_scalar("(-I)/Sqrt[2]")), -1j/np.sqrt(2))
