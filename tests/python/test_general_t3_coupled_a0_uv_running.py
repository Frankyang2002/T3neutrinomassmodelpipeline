"""Stage30 local contracts; live RGBeta export is a separate Wolfram run."""
from copy import deepcopy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.CoupledA0UVRunning import (
    BETA_CONVENTION, CONFIG, MODEL, START, STOP, file_hash, integrate_vector,
    quartic_names, uv_quartics, validate_export, verify_stage29_registry_provenance,
)


def fixture():
    names = tuple(f'lambdaT3GeneralQ{i:06d}I1' for i in range(5))
    reg = {'totals':{'rgbeta_new_vertices':38,'matchete_extra_quartics':50},
           'models': {MODEL: {
               'total_additional_quartics':5,
               'preexisting_rgbeta_symbols':[],
               'new_rgbeta_symbols':[{'name':n} for n in names]}}}
    s29={'status':'NumericalBoundaryCandidatesOnlyNoRGEIntegration',
         'models_checked':16,'slots_checked':50,'blocks_checked':46,
         'benchmark_scale_gev':START,'physical_coupling_mapping_certified':False,
         'rge_evolution_performed':False,
         'models':{MODEL:{'rgbeta_quartics':{
             n:{'registry_class':'new','exact':'-1/100', 'float':-0.01}
             for n in names}}}}
    bases=('gY','g2','lambdaH')
    export={'format':'stage30-A-p0-coupled-uv-rgbeta-betas-v1',
            'status':'Success','model':MODEL,
            'physical_coupling_mapping_certified':False,
            'metadata':{'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':0,
                        'ReportBetaConvention':BETA_CONVENTION},
            'base_beta_symbols':list(bases),
            'new_quartic_symbols':list(names),
            'report_betas':{x:'0' for x in bases + names}}
    return reg,s29,export,names,bases


def test_source_fixture_exact_five():
    reg,*_=fixture()
    assert len(quartic_names(reg))==5


def test_uv_boundary_exact_conventional_value():
    reg,s29,_,names,_=fixture()
    assert uv_quartics(s29,reg)=={name:-0.01+0j for name in names}


def test_eight_beta_names_all_accounted():
    _,_,export,names,bases=fixture()
    assert set(validate_export(export,names,bases))==set(bases)|set(names)

@pytest.mark.parametrize('change',[
    lambda x: x['models'][MODEL]['rgbeta_quartics'].pop(next(iter(x['models'][MODEL]['rgbeta_quartics']))),
    lambda x: x.update(physical_coupling_mapping_certified=True),
    lambda x: x.update(rge_evolution_performed=True),
    lambda x: x.update(benchmark_scale_gev=1e6),
])
def test_reject_bad_stage29(change):
    reg,s29,_,_,_=fixture()
    change(s29)
    with pytest.raises(ValueError,match='Stage30 fail-closed'):
        uv_quartics(s29,reg)


def test_reject_mismatching_boundary_decimal():
    reg,s29,_,names,_=fixture()
    s29['models'][MODEL]['rgbeta_quartics'][names[0]]['float']=-0.02
    with pytest.raises(ValueError,match='exact and float'):
        uv_quartics(s29,reg)


def test_reject_symbolic_boundary():
    reg,s29,_,names,_=fixture()
    s29['models'][MODEL]['rgbeta_quartics'][names[0]]['exact']='x'
    with pytest.raises(ValueError,match='not an exact real'):
        uv_quartics(s29,reg)

@pytest.mark.parametrize('alter',[
    lambda x: x['metadata'].update(dS2=2),
    lambda x: x['metadata'].update(ReportBetaConvention='native gauge dg^2'),
    lambda x: x['base_beta_symbols'].pop(),
    lambda x: x['new_quartic_symbols'].pop(),
    lambda x: x['report_betas'].pop('gY'),
    lambda x: x.update(physical_coupling_mapping_certified=True),
    lambda x: x['report_betas'].update(gY='BetaTerm[gY,1]'),
    lambda x: x['report_betas'].update(gY=''),
])
def test_reject_modified_or_incomplete_live_beta_export(alter):
    _,_,export,names,bases=fixture()
    alter(export)
    with pytest.raises(ValueError,match='Stage30 fail-closed'):
        validate_export(export,names,bases)


def test_backward_integrator_matches_exact_exponential():
    vec=np.array([0.1,0.3,2.0])
    k=0.123
    sol=integrate_vector(vec,lambda t,y:k*y,n_points=13,rtol=1e-11,atol=1e-13)
    assert sol.success
    assert sol.y.shape==(3,13)
    assert np.allclose(sol.y[:,-1],vec*(STOP/START)**k,rtol=3e-11,atol=1e-12)
    assert np.allclose(sol.y[:,0],vec,rtol=0,atol=1e-14)


def test_negative_scale_ratio_and_no_hidden_loop_factor():
    y0=np.array([1.])
    b=16*math.pi**2
    result=integrate_vector(y0,lambda t,y:np.array([b/(16*math.pi**2)]),n_points=4)
    assert abs(result.y[0,-1]-(1+math.log(STOP/START)))<1e-8

@pytest.mark.parametrize('y0,n',[(np.array([np.nan]),5),(np.array([1.,2.]),2)])
def test_bad_numerical_settings_rejected(y0,n):
    with pytest.raises(ValueError,match='Stage30 fail-closed'):
        integrate_vector(y0,lambda t,y:y,n_points=n)


def test_json_export_fixture_roundtrips():
    _,_,export,names,bases=fixture()
    assert set(validate_export(json.loads(json.dumps(export)),names,bases))==set(names)|set(bases)


def test_no_production_code_replacements_in_stage30_package():
    root=Path(__file__).resolve().parents[2]
    script=(root/'scripts/CoupledA0UVRunning.py').read_text(encoding='utf-8')
    wl=(root/'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    assert 'T3RGBetaOneLoopBetas[]' in wl and 'BetaTerm[Quartic,1]' in wl
    assert 'beta_values_to_derivative' in script
    assert 'state_environment' in script
    assert 'DOP853' in script
    assert 'lambdaT3GeneralQ' in wl


def test_stage30_fixture_is_explicitly_uncertified():
    _,_,e,_,_=fixture()
    assert e['physical_coupling_mapping_certified'] is False


def test_coupled_rhs_includes_cross_feedback_and_complex_parts(monkeypatch):
    """Probe the actual sidecar packing (mock only the external core modules)."""
    import sys
    import types
    from scripts.CoupledA0UVRunning import combined_rhs
    names=('lambdaT3GeneralQ000211I1', 'lambdaT3GeneralQ021010I1')
    class Layout:
        size=2
        names=('gY','lambdaH')
    def unpack(vec,layout,mu_gev):
        return types.SimpleNamespace(gY=vec[0],lambdaH=vec[1],mu_gev=mu_gev)
    def env(s):
        return {'gY':s.gY,'lambdaH':s.lambdaH}
    def parser(expr,environment):
        return environment[expr]
    def base_pack(values,layout,divide_by_loop_factor):
        assert divide_by_loop_factor is True
        return np.array([values[k].real for k in layout.names]) / (16*math.pi**2)
    imports={
        'Numerical':types.ModuleType('Numerical'),
        'Numerical.core':types.ModuleType('Numerical.core'),
        'Numerical.core.StateVector':types.ModuleType('Numerical.core.StateVector'),
        'Numerical.core.RGBetaEvaluator':types.ModuleType('Numerical.core.RGBetaEvaluator'),
        'Numerical.core.BetaVector':types.ModuleType('Numerical.core.BetaVector'),
    }
    imports['Numerical.core.StateVector'].unpack_uv_state=unpack
    imports['Numerical.core.RGBetaEvaluator'].state_environment=env
    imports['Numerical.core.RGBetaEvaluator'].evaluate_inputform_expression=parser
    imports['Numerical.core.BetaVector'].beta_values_to_derivative=base_pack
    for key,val in imports.items():
        monkeypatch.setitem(sys.modules,key,val)
    betas={
        'gY':names[0],  # standard beta depends on the extra quartic
        'lambdaH':names[1],
        names[0]:'gY',  # new beta depends on the standard gauge coupling
        names[1]:'lambdaH',
    }
    initial=np.array([0.4,0.2, -0.01,0.003, -0.02,-0.004])
    derivative=combined_rhs(math.log(START),initial,layout=Layout(),names=names,betas=betas)
    expected=np.array([-0.01,-0.02,0.4,0,0.2,0])/(16*math.pi**2)
    assert np.allclose(derivative,expected,rtol=0,atol=1e-15)
    assert derivative.shape==initial.shape


def test_new_coupling_beta_must_be_a_scalar(monkeypatch):
    import sys
    import types
    from scripts.CoupledA0UVRunning import combined_rhs
    n='lambdaT3GeneralQ000211I1'
    layout=types.SimpleNamespace(size=1,names=('gY',))
    modules={
        'Numerical':types.ModuleType('Numerical'),
        'Numerical.core':types.ModuleType('Numerical.core'),
        'Numerical.core.StateVector':types.ModuleType('Numerical.core.StateVector'),
        'Numerical.core.RGBetaEvaluator':types.ModuleType('Numerical.core.RGBetaEvaluator'),
        'Numerical.core.BetaVector':types.ModuleType('Numerical.core.BetaVector'),
    }
    modules['Numerical.core.StateVector'].unpack_uv_state=lambda y,lay,mu_gev:object()
    modules['Numerical.core.RGBetaEvaluator'].state_environment=lambda s:{'gY':0.1}
    modules['Numerical.core.RGBetaEvaluator'].evaluate_inputform_expression=lambda expr,env:(
        np.eye(2) if expr=='MATRIX' else 0.0)
    modules['Numerical.core.BetaVector'].beta_values_to_derivative=lambda d,l,divide_by_loop_factor:np.zeros(1)
    for key,val in modules.items():monkeypatch.setitem(sys.modules,key,val)
    with pytest.raises(ValueError,match='matrix-valued'):
        combined_rhs(math.log(START),np.array([0.4,0.01,0.0]),layout=layout,
                     names=(n,),betas={'gY':'0','lambdaT3GeneralQ000211I1':'MATRIX'})


def test_rgbeta_loaded_before_beta_symbol_is_parsed():
    """Do not pre-intern Global`BetaTerm before RGBeta declares its own symbol."""
    root = Path(__file__).resolve().parents[2]
    wl = (root / 'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    needs_position = wl.index('Needs["RGBeta`"]')
    exporter_definition_position = wl.index('T3S30Export[] :=')
    beta_term_position = wl.index('BetaTerm[Quartic,1]')
    assert needs_position < exporter_definition_position < beta_term_position
    assert wl.count('Needs["RGBeta`"]') == 1
    assert 'CheckAbort[Quiet@Check[Needs["RGBeta`"]' in wl


def test_all_quartic_betas_use_native_full_class_projection():
    """No single-quartic BetaTerm for new slots: RGBeta mixes invariants."""
    root = Path(__file__).resolve().parents[2]
    wl = (root / 'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    assert 'BetaTerm[Quartic,1]' in wl
    assert 'BetaTerm[Symbol[name],1]' not in wl
    assert 'projectedBase=Association@Table' in wl
    assert 'extended=Association@Table' in wl
    assert 'beta=Join[beta,projectedBase]' in wl
    assert 'Sort[expectedQ]=!=actualQ' in wl
    assert 'nativeResult=BetaTerm[Quartic,1]' in wl
    assert 'quartic=CheckAbort[' in wl
    assert 'Quiet@Check[BetaTerm[Quartic,1]' not in wl


def test_projected_quartic_failure_remains_fail_closed():
    root = Path(__file__).resolve().parents[2]
    wl = (root / 'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    assert 'Projected new quartic beta missing or unevaluated:' in wl
    assert 'Projected base quartic beta missing or unevaluated:' in wl
    assert 'class-level BetaTerm[Quartic,1] returned' in wl
    assert 'Missing["AbsentQuarticBeta",name]' in wl


def test_native_projector_diagnostics_not_masked_by_check():
    """A Wolfram message may accompany a valid Association; Check is too strict."""
    root = Path(__file__).resolve().parents[2]
    wl = (root/'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    assert 'Block[{$MessageList={}}' in wl
    assert 'nativeResult=BetaTerm[Quartic,1]' in wl
    assert 'CheckAbort[' in wl
    assert 'ProjectionToUnmixedBetas::singular' in wl
    assert 'a projector singularity cannot be bypassed' in wl
    assert 'Check[BetaTerm[Quartic,1]' not in wl


def test_outer_export_must_not_use_check_on_native_monitor_warning():
    """FrontEndObject::notavail from Monitor must not discard good beta output.

    Fail-closed guards remain: native abort, explicit Throw, result-schema checks.
    """
    wl=(Path(__file__).resolve().parents[2] /
        'scripts/ExportGeneralT3A0UVOneLoop.wl').read_text(encoding='utf-8')
    assert 'result=CheckAbort[Catch[T3S30Export[],"T3S30"],$Aborted];' in wl
    assert 'Check[Catch[T3S30Export[]' not in wl
    assert 'If[!AssociationQ[result]||Lookup[result,"status",""]=!="Success"' in wl
    assert 'FrontEndObject::notavail' in wl
    assert 'quartic=CheckAbort[' in wl


@pytest.mark.parametrize('path_key', [
    'RGE/general/generated/general_rgbeta_quartic_registry.json',
    r'RGE\general\generated\general_rgbeta_quartic_registry.json',
])
def test_stage29_registry_provenance_accepts_windows_and_posix_keys(tmp_path, path_key):
    registry = tmp_path / 'general_rgbeta_quartic_registry.json'
    registry.write_text('{"verified": "example"}', encoding='utf-8')
    s29 = {'input_sha256': {path_key: file_hash(registry)}}
    verify_stage29_registry_provenance(s29, registry)


def test_stage29_registry_provenance_rejects_stale_hash(tmp_path):
    registry = tmp_path / 'general_rgbeta_quartic_registry.json'
    registry.write_text('old', encoding='utf-8')
    s29 = {'input_sha256': {
        r'RGE\general\generated\general_rgbeta_quartic_registry.json': file_hash(registry)}}
    registry.write_text('new', encoding='utf-8')
    with pytest.raises(ValueError, match='source changed'):
        verify_stage29_registry_provenance(s29, registry)


@pytest.mark.parametrize('mapping', [
    {},
    {'RGE/general/generated/general_rgbeta_quartic_registry.json': 'a',
     r'RGE\general\generated\general_rgbeta_quartic_registry.json': 'a'},
])
def test_stage29_registry_provenance_rejects_missing_or_duplicate_key(tmp_path, mapping):
    registry = tmp_path / 'general_rgbeta_quartic_registry.json'
    registry.write_text('example', encoding='utf-8')
    with pytest.raises(ValueError, match='absent or ambiguous'):
        verify_stage29_registry_provenance({'input_sha256': mapping}, registry)
