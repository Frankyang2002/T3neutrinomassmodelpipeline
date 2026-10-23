"""Regressions for local T3-B GroupFactors HTML companions."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

try:
    from Reports.GroupFactorCoupledHTML import (
        render_group_factor_html, write_group_factor_html_companions, main,
    )
except ModuleNotFoundError:
    # Allows tests on the standalone ZIP before extraction into the project.
    module_path = Path(__file__).resolve().parents[1] / 'Reports' / 'GroupFactorCoupledHTML.py'
    spec = importlib.util.spec_from_file_location('t3b_groupfactor_html', module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    render_group_factor_html = module.render_group_factor_html
    write_group_factor_html_companions = module.write_group_factor_html_companions
    main = module.main


def _record(*, d_s1=2, d_s2=2, d_f=1, alpha=-1, label='EFT_1_after_F',
            active=('S1', 'S2'), integrated=('F',), name='T3_B_alpha_m1'):
    return SimpleNamespace(
        name=name, d_s1=d_s1, d_s2=d_s2, d_f=d_f, alpha=alpha,
        shared_scalar=False,
        summary={'EFTStages': [dict(Label=label, ActiveHeavyFields=list(active),
                                    IntegratedFields=list(integrated))]},
    )


def _tex(tmp_path, stem='GF_EFT_1_after_F'):
    folder = tmp_path / 'GroupFactors'
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / (stem + '.tex')
    path.write_text('\\documentclass{article}', encoding='utf-8')
    return path


def test_t3b_f_first_contains_two_wilson_betas_and_factors(tmp_path):
    source = _tex(tmp_path)
    out, = write_group_factor_html_companions([_record()], [source])
    html = out.read_text(encoding='utf-8')
    assert 'Coupled dimension-five RGEs: T3-B' in html
    assert 'dC<sub>12</sub>/dln&mu;' in html
    assert 'dC<sub>5</sub>/dln&mu;' in html
    assert '(8&radic;3/3)&lambda;<sub>T3</sub><sup>*</sup>' in html
    assert '(4&radic;3/3)&lambda;<sub>T3</sub>' in html
    assert 'E<sup>T</sup>' in html
    assert 'GF_EFT_1_after_F.pdf' in html
    assert 'ordered Wilson-tensor convention' in html


def test_non_t3b_model_cannot_inherit_t3b_group_factors(tmp_path):
    path = _tex(tmp_path)
    out, = write_group_factor_html_companions(
        [_record(d_s1=1, d_s2=3, d_f=2, alpha=-2)], [path]
    )
    html = out.read_text(encoding='utf-8')
    assert 'T3-B-specific' in html
    assert '(8&radic;3/3)' not in html


def test_wrong_threshold_order_cannot_claim_f_first_equation(tmp_path):
    path = _tex(tmp_path)
    out, = write_group_factor_html_companions(
        [_record(active=('F', 'S2'), integrated=('S1',))], [path]
    )
    assert 'dC<sub>12</sub>/dln' not in out.read_text(encoding='utf-8')


def test_fully_decoupled_stage_shows_only_smeft_weinberg(tmp_path):
    path = _tex(tmp_path, 'GF_EFT_2_after_S1_S2')
    rec = _record(label='EFT_2_after_S1_S2', active=(), integrated=('S1', 'S2'))
    out, = write_group_factor_html_companions([rec], [path])
    html = out.read_text(encoding='utf-8')
    assert 'Final SMEFT: Weinberg operator running' in html
    assert 'dC<sub>5</sub>/dln' in html
    assert 'mixing no longer occurs' in html
    assert '8&radic;3/3' not in html


def test_uv_stage_excludes_scalar_only_wilson_equation(tmp_path):
    path = _tex(tmp_path, 'GF_UV')
    out, = write_group_factor_html_companions([_record()], [path])
    html = out.read_text(encoding='utf-8')
    assert 'Heavy fermions are still active' in html
    assert 'dC<sub>12</sub>/dln' not in html


def test_explicit_cli_mode_rebuilds_saved_tex_without_pipeline(tmp_path):
    path = _tex(tmp_path)
    assert main(['--report-root', str(tmp_path), '--t3b']) == 0
    generated = path.with_suffix('.html')
    assert generated.is_file()
    assert '8&radic;3/3' in generated.read_text(encoding='utf-8')


def test_cli_without_explicit_model_does_not_guess(tmp_path):
    path = _tex(tmp_path)
    assert main(['--report-root', str(tmp_path)]) == 0
    assert '8&radic;3/3' not in path.with_suffix('.html').read_text(encoding='utf-8')


def test_html_escapes_untrusted_record_names(tmp_path):
    path = _tex(tmp_path)
    name = '<script>alert(1)</script>'
    out, = write_group_factor_html_companions([_record(name=name)], [path])
    html = out.read_text(encoding='utf-8')
    assert name not in html
    assert '&lt;script&gt;' in html


def test_missing_source_tex_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        write_group_factor_html_companions([], [tmp_path / 'GF_EFT_1_after_F.tex'],
                                           explicit_t3b=True)


def test_html_is_complete_and_offline_readable():
    html = render_group_factor_html(
        stage_label='EFT_1_after_F', pdf_name='test.pdf', tex_name='test.tex',
        t3b_names=['T3-B'],
    )
    assert html.startswith('<!doctype html>')
    assert html.rstrip().endswith('</html>')
    assert 'MathJax' not in html
    assert '<table>' in html
    assert 'Copyable LaTeX equations' in html


def test_existing_with_c12_variant_gets_rge_when_t3b_selected(tmp_path):
    path = _tex(tmp_path, 'GF_EFT_1_after_F_with_C12')
    assert main(['--report-root', str(tmp_path), '--t3b']) == 0
    assert '8&radic;3/3' in path.with_suffix('.html').read_text(encoding='utf-8')


def test_integrated_opt_in_regenerates_saved_groupfactor_html(monkeypatch, tmp_path):
    """Exercise the automatic report hook without running expensive physics."""
    import json
    import sys
    import types
    import importlib.util

    module_path = (Path(__file__).resolve().parents[1] / 'Numerical'
                   / 'orchestration' / 'T3BIntegratedCoupledPipeline.py')
    spec = importlib.util.spec_from_file_location('t3b_html_hook_test', module_path)
    coupled = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(coupled)

    # Relocate its report search root to the disposable test project.
    project = tmp_path / 'project'
    monkeypatch.setattr(coupled, '__file__', str(project / 'Numerical'
                       / 'orchestration' / 'T3BIntegratedCoupledPipeline.py'))
    saved_tex = project / 'Reports' / 'output' / 'analytical' / 'GroupFactors' / 'GF_EFT_1_after_F.tex'
    saved_tex.parent.mkdir(parents=True)
    saved_tex.write_text('\\documentclass{article}', encoding='utf-8')

    config = project / 'benchmark.json'
    config.write_text(json.dumps({'representation': {
        'd_s1': 2, 'd_s2': 2, 'd_f': 1, 'alpha': -1},
        'benchmark_search': {'enabled': False}}), encoding='utf-8')
    output_dir = project / 'output'
    production_file = output_dir / 'data' / 'running_diagnostics.json'
    expected_production = '{"status":"Success","authoritative":true}\n'

    def production(_record, _config):
        production_file.parent.mkdir(parents=True)
        production_file.write_text(expected_production, encoding='utf-8')
        return True

    numerical_pkg = types.ModuleType('Numerical')
    numerical_pkg.__path__ = []
    orchestr_pkg = types.ModuleType('Numerical.orchestration')
    orchestr_pkg.__path__ = []
    production_module = types.ModuleType('Numerical.orchestration.PipelineNumericalResults')
    production_module.run_pipeline_numerical_results = production
    monkeypatch.setitem(sys.modules, 'Numerical', numerical_pkg)
    monkeypatch.setitem(sys.modules, 'Numerical.orchestration', orchestr_pkg)
    monkeypatch.setitem(sys.modules, 'Numerical.orchestration.PipelineNumericalResults',
                        production_module)
    common_pkg = types.ModuleType('common')
    common_pkg.__path__ = []
    records_module = types.ModuleType('common.RunRecords')
    class Record:
        def __init__(self, **attrs):
            self.__dict__.update(attrs)
    records_module.RunRecord = Record
    monkeypatch.setitem(sys.modules, 'common', common_pkg)
    monkeypatch.setitem(sys.modules, 'common.RunRecords', records_module)

    monkeypatch.setattr(coupled, 'compute_opt_in_coupled', lambda **kwargs: {
        'coupled_at_ms': {'relative_difference_from_fixed_order': 0.025,
                          'relative_reverse_feedback': 1.5e-5},
        'low_energy': {'relative_c5_difference_at_low_scale': 0.025},
        'perturbative_accuracy': {'bookkeeping': {'terms': {
            'order_2': {'norm_relative_to_order1': .033},
            'order_3': {'reverse_contribution_relative_to_order1': .000015},
        }}},
    })
    result = coupled.run(config=config, output_dir=output_dir,
                         enable_coupled=True)
    rendered = saved_tex.with_suffix('.html')
    assert rendered.exists()
    assert str(rendered) in result['group_factor_html_files']
    assert 'dC<sub>12</sub>/dln' in rendered.read_text(encoding='utf-8')
    assert production_file.read_text(encoding='utf-8') == expected_production
