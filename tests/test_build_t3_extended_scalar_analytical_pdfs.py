"""Static rendering helpers only; live 16-model source validation uses --check."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'BuildT3ExtendedScalarAnalyticalPDFs.py'
spec = importlib.util.spec_from_file_location('extended_scalar_pdf_builder', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_field_content_conjugate_powers():
    assert module.field_content([0,2,1,0,1,0]) == (
        r'\left(H^\dagger\right)^{2}\, S_1\, S_2')
    assert module.field_content([1,1,0,0,0,0]) == r'H\, H^\dagger'


def test_tex_names():
    assert module.model_tex('A_+0') == r'\mathrm{T3-A}\,(\alpha=0)'
    assert module.model_tex('C_-1') == r'\mathrm{T3-C}\,(\alpha=-1)'
    assert module.tag('C_-1') == 'C_m1'


def test_literal_source_preserved():
    source = '4*Bar[lambdaT3]*(Tr[y1 . Trans[Bar[y2]]])'
    assert source in module.listing(source)
    assert r'\begin{lstlisting}' in module.listing(source)


def test_document_has_explicit_candidate_flags():
    content = module.document('Test', 'Overview', ['Result'])
    assert 'not certified' in content
    assert 'paper coupling notation' in content
    assert content.rstrip().endswith(r'\end{document}')


def test_scalars_are_not_original_lambda_t3():
    registry={'models':{m:{'dimensions':[1,3,2],'alpha':int(m.split('_')[1]),
         'total_additional_quartics':0,'quadratics_not_implemented':0} for m in module.MODEL_ORDER}}
    table='\n'.join(module.inventory_table(registry))
    assert 'extra $\\lambda$' in table
    assert len(registry['models'])==16


def test_append_to_original_does_not_mutate_source(tmp_path):
    old_root = module.ROOT
    try:
        module.ROOT = tmp_path
        folder=tmp_path/'Reports'/'output'/'analytical'/'RGE'
        folder.mkdir(parents=True)
        source=folder/'UV.tex'
        original=r"\documentclass{article}\begin{document}Original\end{document}"
        source.write_text(original,encoding='utf-8')
        merged=module.append_to_existing_report('RGE',appendix=[r'\section*{New}Results'])
        assert 'Original' in merged and r'\section*{New}' in merged
        assert r'\usepackage{listings}' in merged
        assert source.read_text(encoding='utf-8')==original
    finally:
        module.ROOT = old_root


def test_notation_preserves_old_coupling_labels():
    assert module.BASE_PAPER_NAMES['lambdaT3'] == r'\lambda_{5}'
    assert module.REUSED_PAPER_NAMES['lambdaHHdagS1barS1bar'] == r'\lambda_{7}'
    assert module.REUSED_PAPER_NAMES['lambdaHHdagS1barS2barCross'] == r'\lambda_{12}'


def test_new_couplings_have_distinct_global_labels():
    reg={'models':{m:{'new_rgbeta_symbols':[], 'preexisting_rgbeta_symbols':[]}
                    for m in module.MODEL_ORDER}}
    a='lambdaT3GeneralQ021010I1'
    b='lambdaT3GeneralQ000211I1'
    c='lambdaT3GeneralQ000211I2'
    reg['models']['A_+0']['new_rgbeta_symbols']=[{'name':a,'counts':[0,2,1,0,1,0]},
        {'name':b,'counts':[0,0,0,2,1,1]}]
    reg['models']['E_+0']['new_rgbeta_symbols']=[{'name':c,'counts':[0,0,0,2,1,1]},
        {'name':b,'counts':[0,0,0,2,1,1]}]
    reg['models']['B_-1']['preexisting_rgbeta_symbols']=[
        {'existing_name':'lambdaHHdagS2S2','counts':[0,2,0,0,2,0]}]
    mapping=module.notation_map(reg)
    assert mapping[b]==r'\lambda_{13}'
    assert mapping[c]==r'\lambda_{14}'
    assert mapping[a]==r'\lambda_{15}'
    assert mapping['lambdaHHdagS2S2']==r'\lambda_{8}'
    assert len({mapping[a],mapping[b],mapping[c]})==3


def test_latex_math_listing_hides_internal_source_names_without_changing_coefficients():
    m=dict(module.BASE_PAPER_NAMES, **module.REUSED_PAPER_NAMES,
       lambdaT3GeneralQ021010I1=r'\lambda_{13}')
    raw='4*Bar[lambdaT3GeneralQ021010I1]*lambdaT3 + 2*lambdaH1 - 3*Bar[y1]'
    rendered=module.listing(raw,m)
    assert 'lambdaT3GeneralQ' not in rendered
    assert r'$\lambda_{13}^{*}$' in rendered
    assert r'$\lambda_{5}$' in rendered
    assert r'$\lambda_{HS_1}^{(1)}$' in rendered
    assert r'$y_1^{*}$' in rendered
    assert '4*' in rendered and ' + 2*' in rendered
    assert 'mathescape=true' in rendered


def test_reject_unmapped_extra_scalar_source():
    import pytest
    with pytest.raises(ValueError, match='unmapped source symbol'):
        module.listing('lambdaT3GeneralQ999999I1',module.BASE_PAPER_NAMES)


def test_conjugated_quadratic_and_decorated_couplings_are_tex_safe():
    """Regression: TeX rejects ``\\mu^2^{*}`` as a double superscript."""
    m = module.BASE_PAPER_NAMES
    raw = ('12*Bar[muT3GeneralS2Sq] + 4*Bar[muT3GeneralS1Sq]'
           ' + Bar[muT3GeneralS1S2] + 3*Bar[lambdaH1]'
           ' + 2*Bar[lambda12Cross] + Bar[lambdaT3]')
    converted = module.listing(raw, m)
    assert r'$\left(\mu_{S_2}^2\right)^{*}$' in converted
    assert r'$\left(\mu_{S_1}^2\right)^{*}$' in converted
    assert r'$\left(\mu_{12}^2\right)^{*}$' in converted
    assert r'$\left(\lambda_{HS_1}^{(1)}\right)^{*}$' in converted
    assert r'$\left(\lambda_{12}^{(\times)}\right)^{*}$' in converted
    assert r'$\lambda_{5}^{*}$' in converted
    for source in ('muT3GeneralS1Sq', 'muT3GeneralS2Sq',
                   'muT3GeneralS1S2', 'lambdaH1', 'lambda12Cross'):
        assert 'Bar[' + source + ']' not in converted
    assert r'\mu_{S_2}^2^{*}' not in converted


def test_all_known_paper_symbols_can_be_conjugated():
    """Every mapped source symbol must be safe when it occurs in Bar[...]."""
    notation = dict(module.BASE_PAPER_NAMES, **module.REUSED_PAPER_NAMES)
    notation['lambdaT3GeneralQ021010I1'] = r'\lambda_{13}'
    for name, symbol in notation.items():
        result = module.listing('Bar[' + name + ']', notation)
        expected = module.conjugated_symbol_tex(symbol)
        assert '$' + expected + '$' in result
        assert not (symbol + r'^{*}' in result and '^' in symbol)


def test_render_actual_mathematical_beta_not_inputform():
    import importlib.util
    from pathlib import Path
    file = Path(__file__).resolve().parents[1] / 'scripts' / 'BuildT3ExtendedScalarAnalyticalPDFs.py'
    spec = importlib.util.spec_from_file_location('t3_math_display', file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    paper = dict(module.BASE_PAPER_NAMES)
    paper['lambdaT3GeneralQ021010I1'] = r'\lambda_{14}'
    expression = r'4*Bar[lambdaT3GeneralQ021010I1]*lambdaT3 + Tr[y1 . Trans[Bar[y2]]] - 12*g2^2*lambdaT3'
    result = module.beta_equation(paper['lambdaT3'], expression, paper)
    assert r'\operatorname{Tr}' in result
    assert r'\lambda_{14}^{*}' in result
    assert r'g_2^{2}' in result
    assert 'Trans[' not in result and 'Tr[' not in result and '*lambda' not in result


def test_render_mass_square_conjugate_without_double_superscript():
    import importlib.util
    from pathlib import Path
    file = Path(__file__).resolve().parents[1] / 'scripts' / 'BuildT3ExtendedScalarAnalyticalPDFs.py'
    spec = importlib.util.spec_from_file_location('t3_mass_display', file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    expression = '12*Bar[muT3GeneralS2Sq]*lambdaT3'
    result = module.beta_equation(r'\lambda_5', expression, module.BASE_PAPER_NAMES)
    assert r'\left(\mu_{S_2}^2\right)^{*}' in result
    assert '^2^{*}' not in result
