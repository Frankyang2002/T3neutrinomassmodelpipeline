from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'scripts' / 'SweepT3Frozen18Radius.py'


def test_syntax_and_stage_names():
    tree = ast.parse(SRC.read_text(encoding='utf-8'))
    assert any(isinstance(x, ast.FunctionDef) and x.name == 'stage_id' for x in tree.body)
    assert 'frozen' in SRC.read_text(encoding='utf-8').lower()


def test_radius_schedule():
    text = SRC.read_text(encoding='utf-8')
    assert 'RADII = (0.005, 0.01, 0.02)' in text
    assert "if result.returncode not in (0, 2)" in text
    assert 'imports.insert(0, previous)' in text


def test_checkpoint_checks():
    text = SRC.read_text(encoding='utf-8')
    assert 'checked_completed(state, args.model, args.case)' in text
    assert 'raise FileExistsError' in text
