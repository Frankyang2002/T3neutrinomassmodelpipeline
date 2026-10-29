from pathlib import Path
import ast

SOURCE = Path(__file__).resolve().parents[1] / 'scripts' / 'RunT3FrozenEstablishedSobol.py'

def test_source_compiles():
    compile(SOURCE.read_text(encoding='utf-8'), str(SOURCE), 'exec')

def test_uses_established_refinement():
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    imports=[n for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
    assert any(n.module=='Numerical.fitting.SobolBenchmarkSearch' and
        {'generate_sobol_parameter_points','refine_sobol_candidates'} <= {x.name for x in n.names}
        for n in imports)

def test_frozen_only():
    s=SOURCE.read_text(encoding='utf-8')
    assert "CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')" in s
    assert "if history['uv_on']:" in s

def test_preserves_existing():
    s=SOURCE.read_text(encoding='utf-8')
    assert "not args.refit_accepted" in s
    assert 'previous_best' in s

def test_staged_optimisation():
    s=SOURCE.read_text(encoding='utf-8')
    assert 'active_counts=(5,8,12)' in s
    assert 'seed_count=args.seed_count' in s
