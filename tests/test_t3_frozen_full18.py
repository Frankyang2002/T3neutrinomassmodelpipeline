import ast
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'RunT3FrozenEstablishedSobol.py'

def test_syntax():
    ast.parse(SOURCE.read_text(encoding='utf-8'))

def test_full_18_stage_present():
    text = SOURCE.read_text(encoding='utf-8')
    assert 'active_counts=(5,8,12,18)' in text
    assert 'active_counts=(18,)' in text
    assert "--resume-from-run" in text
    assert "--refine-only" in text

def test_previous_run_preserved():
    text=SOURCE.read_text(encoding='utf-8')
    assert "old_dir.glob('round_*.json')" in text
    assert "output / 'full18_refinement.json'" in text
