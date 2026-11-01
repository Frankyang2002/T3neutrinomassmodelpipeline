"""Presentation-only tests; no Wolfram or physics calculations."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / 'Reports' / 'GeneralSymbolicPresentation.py'
ADAPTER = ROOT / 'Reports' / 'GeneralSymbolicAdapter.py'


def _functions():
    import re
    src = PRESENTATION.read_text(encoding='utf-8')
    tree = ast.parse(src)
    wanted = {'notation_rows', '_notation_symbol', 'field_content'}
    subset = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in wanted]
    module = ast.fix_missing_locations(ast.Module(body=subset, type_ignores=[]))
    scope = {
        're': re,
        'ORIGINAL_DEFINITIONS': {'gY','lambdaH1'},
        'NAMED_QUARTICS': {
            'lambdaS1bar2S2bar2': (r'(S_1^\dagger)^2(S_2^\dagger)^2', r'\lambda_9'),
        },
        'MASS': re.compile(r'(?:mu|m)T3[A-Za-z0-9_]+'),
        'latex_escape_text': lambda x: x.replace('_', r'\_'),
        'generated_symbol': lambda name: r'\lambda_{X_1}' if 'T3' in name else None,
        '_registry': lambda: {'models': {'A_-2': {
            'new_rgbeta_symbols': [{
                'name': 'lambdaT3GeneralQ110101I1',
                'counts': [1,1,0,1,0,1],
                'invariant': 'del[SU2L, #1, #2]',
                'matchete_basis_name':'lambdaT3ExtraQ110101I#',
                'matchete_count':1,
            }], 'preexisting_rgbeta_symbols': []}}},
    }
    from types import SimpleNamespace
    # Avoid imports from the repository by stubbing the source notation map.
    import sys
    old = sys.modules.get('Reports.ReportGeneration')
    sys.modules['Reports.ReportGeneration'] = SimpleNamespace(PAPER_SYMBOL_LATEX={})
    try:
        exec(compile(module, str(PRESENTATION), 'exec'), scope)
        return scope
    finally:
        if old is None:
            sys.modules.pop('Reports.ReportGeneration',None)
        else:
            sys.modules['Reports.ReportGeneration'] = old


def test_every_encountered_coupling_has_index_row():
    import sys
    from types import SimpleNamespace
    scope = _functions()
    old = sys.modules.get('Reports.ReportGeneration')
    sys.modules['Reports.ReportGeneration'] = SimpleNamespace(PAPER_SYMBOL_LATEX={})
    try:
        rows = scope['notation_rows'](['A_-2'], {
            'lambdaS1bar2S2bar2','y1Conj','muT3GeneralS1Sq',
            'lambdaT3ExtraQ110101I1','unregisteredCoupling','gY'})
    finally:
        if old is None:sys.modules.pop('Reports.ReportGeneration', None)
        else:sys.modules['Reports.ReportGeneration']=old
    rendered = ' '.join(' '.join(row) for row in rows)
    assert r'\lambda_9' in rendered
    assert 'y1Conj' in rendered
    assert 'muT3GeneralS1Sq' in rendered
    assert 'lambdaT3ExtraQ110101I1' in rendered
    assert 'unregisteredCoupling' in rendered
    assert 'gY' not in rendered
    assert 'del[SU2L' in rendered


def test_original_writers_still_used_with_dynamic_inventory():
    source = ADAPTER.read_text(encoding='utf-8')
    ast.parse(source)
    assert 'lagrangian_names.update' in source
    assert "all_names['UV'] | all_names['EFT1'] | lagrangian_names" in source
    for name in ('write_bsm_uv_field_table(', 'write_rge_comparison(',
                 'write_eft1_rge_comparison(', '_saved_rge_beta_sections('):
        assert name in source
