"""Isolated schedule checks; live runner requires the project's physics modules."""
import ast
from pathlib import Path
import pytest

p = Path(__file__).resolve().parents[1] / 'scripts' / 'RunT3FrozenEstablishedSobol.py'
module = ast.parse(p.read_text(encoding='utf-8'))
fn = next(node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == 'points_for_round')
ns = {}
exec(compile(ast.Module(body=[fn], type_ignores=[]), str(p), 'exec'), ns)
points_for_round = ns['points_for_round']

def test_escalation():
    assert [points_for_round(128, i, 4096) for i in range(1,9)] == [128,256,512,1024,2048,4096,4096,4096]

def test_invalid_counts():
    with pytest.raises(ValueError): points_for_round(100,1,4096)
    with pytest.raises(ValueError): points_for_round(128,1,64)
