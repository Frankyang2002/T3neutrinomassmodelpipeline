"""Offline tests for the additive three-stage T3 report adapter.

Run inside the full project as:
    python -X utf8 -m pytest -q tests/python/test_t3_three_stage_rge.py
These do not test the live numerical integration; its saved-endpoint parity gate
runs when the report is built from actual per-model centre data.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

MODULE = Path(__file__).resolve().parents[2] / 'scripts' / 'BuildT3SharedA2ThreeStageRGE.py'
spec = importlib.util.spec_from_file_location('t3_three_stage_report', MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_two_saved_matrix_encodings():
    mat = np.array([[1. + 0j, 2. - 1j, 0j], [2. - 1j, 3. + 0j, 0j], [0j, 0j, 4j]])
    compact = {'real': mat.real.tolist(), 'imag': mat.imag.tolist()}
    per_entry = [[{'real': float(z.real), 'imag': float(z.imag)} for z in row] for row in mat]
    np.testing.assert_allclose(mod.matrix(compact), mat)
    np.testing.assert_allclose(mod.matrix(per_entry), mat)
    with pytest.raises(ValueError):
        mod.matrix({'real': [[1]], 'imag': [[0]]})


def test_scalar_validation():
    assert mod.real({'real': 0.3, 'imag': 0.}) == pytest.approx(0.3)
    with pytest.raises(ValueError):
        mod.real({'real': 1, 'imag': 1})


def test_combine_never_fabricates_eft1_c5():
    old = {'status': 'Success',
           'uv': {'mu_gev': [1e7, 1e5], 'series': {'g2': [.5, .6]}},
           'eft1': {'mu_gev': [1e5, 1e3], 'series': {'g2': [.6, .7]}}}
    new = {'status': 'Success', 'mu_gev': [1e3, 1e2], 'series': {'m3_ev': [.06, .05]}}
    combined = mod.combine(old, new)
    assert set(combined) == {'uv', 'eft1', 'smeft'}
    assert 'c5' not in combined['eft1']['series']
    output = mod.render_html({'A_m2': combined})
    assert 'SMEFT' in output
    assert 'EFT1 C5 is not plotted' in output
    assert 'm3_ev' in output
    assert 'script src=' not in output  # must work offline
