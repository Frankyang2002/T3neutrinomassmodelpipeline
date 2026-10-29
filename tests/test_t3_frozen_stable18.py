"""Offline syntax and file-structure regression tests for the stable fitter."""
from pathlib import Path
import ast

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'RefineT3Frozen18Stable.py'

def test_syntax():
    ast.parse(SOURCE.read_text(encoding='utf-8'))

def test_uses_18d_scipy_residual_fit():
    code = SOURCE.read_text(encoding='utf-8')
    assert 'least_squares(residual, seed_z' in code
    assert "jac='3-point'" in code
    assert "x_scale='jac'" in code

def test_no_matching_or_physical_threshold_override():
    code = SOURCE.read_text(encoding='utf-8')
    assert 'bulk.Evaluator(' in code
    assert 'bulk.matching(' in code
    assert 'bulk.update_yukawas(' in code

def test_has_termination_diagnostics():
    code = SOURCE.read_text(encoding='utf-8')
    assert 'optimizer_nfev' in code
    assert 'optimizer_message' in code
    assert 'optimality' in code
