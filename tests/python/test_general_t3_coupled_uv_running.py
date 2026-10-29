"""Stage31 exact inventory, reuse, projection-contract and ODE regression tests.

Native Wolfram export is tested on the user's RGBeta installation, not mocked
as a successful physical result by this test module.
"""
from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import types

import numpy as np
import pytest

from scripts.CoupledUVRunning import (
    MODELS, START, STOP, BETA_CONVENTION, COMPLEX_REUSED,
    inventories, stage29_values, config_path, prepare_state,
    validate_export, summarize,
)
from scripts.CoupledA0UVRunning import integrate_vector

COUNTS = {'A_-2': (5, 0), 'A_-4': (0, 0), 'B_+1': (1, 0),
          'B_-1': (1, 6), 'B_-3': (1, 0), 'C_+1': (1, 0),
          'C_-1': (1, 6), 'C_-3': (1, 0), 'D_+0': (5, 0),
          'D_+2': (0, 0), 'D_-2': (5, 0), 'E_+0': (6, 0),
          'E_+2': (0, 0), 'E_-2': (6, 0), 'E_-4': (0, 0)}
REUSED = tuple(sorted(COMPLEX_REUSED))


def make_registry():
    reps = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
            'D': (3, 1, 2), 'E': (3, 3, 2)}
    models = {}
    for model, (n, r) in COUNTS.items():
        cls, alpha = model.split('_')
        symbols = [f'lambdaT3GeneralQ{i:06d}I1' for i in range(n)]
        models[model] = {'model': model, 'dimensions': list(reps[cls]),
                         'alpha': int(alpha),
                         'total_additional_quartics': n+r,
                         'new_rgbeta_symbols': [{'name': x} for x in symbols],
                         'preexisting_rgbeta_symbols': [
                             {'existing_name': s} for s in REUSED[:r]]}
    models['A_+0'] = {'model': 'A_+0', 'dimensions': [1,3,2],
                       'alpha': 0, 'total_additional_quartics': 5,
                       'new_rgbeta_symbols': [{'name': f'lambdaT3GeneralQ{i:06d}I1'} for i in range(5)],
                       'preexisting_rgbeta_symbols': []}
    return {'totals': {'matchete_extra_quartics': 50,
                       'rgbeta_new_vertices': 38,
                       'rgbeta_preexisting_structural': 12},
            'models': models}


def make_s29(model, new, reused):
    values = {}
    for sym in new + reused:
        values[sym] = {'registry_class': 'new' if sym in new else 'reused',
                       'exact': '-1/100', 'float': -0.01}
    return {'status': 'NumericalBoundaryCandidatesOnlyNoRGEIntegration',
            'models_checked': 16, 'slots_checked': 50, 'blocks_checked': 46,
            'benchmark_scale_gev': START,
            'physical_coupling_mapping_certified': False,
            'rge_evolution_performed': False,
            'models': {model: {'rgbeta_quartics': values}}}


def make_export(entry, new, reused, base_names):
    d1,d2,df = entry['dimensions']
    return {'format': 'stage31-multi-model-coupled-uv-rgbeta-betas-v1',
            'status': 'Success', 'model': entry['model'],
            'physical_coupling_mapping_certified': False,
            'metadata': {'SharedScalar': False,
                         'dS1': d1, 'dS2': d2, 'dF': df,
                         'alpha': entry['alpha'],
                         'ReportBetaConvention': BETA_CONVENTION},
            'new_quartic_symbols': list(new),
            'reused_quartic_symbols': list(reused),
            'base_beta_symbols': list(base_names),
            'report_betas': {k: '0' for k in base_names + new}}


def test_inventory_complete_and_disjoint():
    reg = make_registry()
    assert len(MODELS) == 15
    counts = [inventories(reg, name) for name in MODELS]
    assert sum(len(new) for _, new, _ in counts) == 33
    assert sum(len(reused) for _, _, reused in counts) == 12
    assert sum(e['total_additional_quartics'] for e,_,_ in counts) == 45
    assert sum(not new and not reused for _,new,reused in counts) == 4


@pytest.mark.parametrize('model', MODELS)
def test_each_registered_model_inventory_and_config_path(model):
    entry, new, reused = inventories(make_registry(), model)
    assert (len(new), len(reused)) == COUNTS[model]
    path = config_path(entry).name
    assert 'T3_dS1_' in path and path.endswith('.json')
    assert f"alpha_{'m' if entry['alpha'] < 0 else 'p'}{abs(entry['alpha'])}" in path


@pytest.mark.parametrize('model', ['A_-4', 'D_+2', 'E_+2', 'E_-4'])
def test_zero_extra_models_require_real_base_running(model):
    entry, new, reused = inventories(make_registry(),model)
    vals = stage29_values(make_s29(model,new,reused),entry,new,reused)
    assert vals == {} and new == reused == ()
    assert validate_export(make_export(entry,new,reused,('gY','lambdaH')),
                           entry,new,reused,('gY','lambdaH')) == {'gY':'0','lambdaH':'0'}


@pytest.mark.parametrize('model', ['B_-1','C_-1'])
def test_reused_quartics_are_separate_from_one_new_slot(model):
    entry,new,reused = inventories(make_registry(),model)
    values = stage29_values(make_s29(model,new,reused),entry,new,reused)
    assert len(values)==7 and len(new)==1 and len(reused)==6
    assert all(values[k] == -0.01+0j for k in values)
    bases=('gY','lambdaH') + reused
    assert set(validate_export(make_export(entry,new,reused,bases),entry,new,reused,bases)) == set(bases+new)


@pytest.mark.parametrize('mutator,pattern', [
    (lambda x: x['models'].pop('A_-4'), 'model selection'),
    (lambda x: x['totals'].update(rgbeta_new_vertices=37), 'totals'),
    (lambda x: x['models']['B_-1']['new_rgbeta_symbols'].append({'name':'lambdaX'}), 'inventory'),
    (lambda x: x['models']['B_-1']['preexisting_rgbeta_symbols'][0].update(existing_name='lambdaUnknown'), 'inventory'),
])
def test_reject_tampered_registry(mutator,pattern):
    reg=make_registry();mutator(reg)
    with pytest.raises(ValueError,match=pattern): inventories(reg,'A_-4' if pattern=='model selection' else 'B_-1')


def test_reject_bad_stage29_exact_or_certification():
    entry,new,reused=inventories(make_registry(),'E_+0')
    stage=make_s29('E_+0',new,reused)
    stage['models']['E_+0']['rgbeta_quartics'][new[0]]['exact']='some_symbol'
    with pytest.raises(ValueError,match='exactly real'):
        stage29_values(stage,entry,new,reused)
    stage=make_s29('E_+0',new,reused)
    stage['physical_coupling_mapping_certified']=True
    with pytest.raises(ValueError,match='Stage29'):
        stage29_values(stage,entry,new,reused)


@pytest.mark.parametrize('field,change', [
    ('metadata',lambda x: x['metadata'].update(dF=3)),
    ('source',lambda x: x['report_betas'].update(gY='BetaTerm[gY,1]')),
    ('inventory',lambda x: x['new_quartic_symbols'].clear()),
    ('certification',lambda x: x.update(physical_coupling_mapping_certified=True)),
    ('reused',lambda x: x['reused_quartic_symbols'].pop()),
])
def test_reject_corrupt_live_rgbeta_export(field,change):
    entry,new,reused=inventories(make_registry(),'B_-1')
    base=('gY','lambdaH')+reused
    exp=make_export(entry,new,reused,base)
    change(exp)
    with pytest.raises(ValueError,match='Stage31 fail-closed'):
        validate_export(exp,entry,new,reused,base)


def test_both_reused_benchmark_values_are_overwritten_in_memory(monkeypatch):
    import scripts.CoupledUVRunning as module
    entry,new,reused=inventories(make_registry(),'B_-1')
    original={'representation':{'d_s1':2,'d_s2':2,'d_f':1,'alpha':-1,'shared_scalar':False},
              'base_state':{'mu_gev':START,'ordinary':{
                  name:{'real':0.0,'imag':0.0} for name in reused}}}
    values={name:complex(-0.015,0) for name in reused}
    captured=[]
    fake_scan=types.ModuleType('Numerical.fitting.ScanCLI')
    fake_vector=types.ModuleType('Numerical.core.StateVector')
    def make_state(ns,params):
        captured.append(ns.raw)
        return types.SimpleNamespace(mu_gev=START,validated=lambda:types.SimpleNamespace(mu_gev=START))
    fake_scan.build_uv_state_from_config=make_state
    fake_vector.pack_uv_state=lambda st: (np.array([0.]),types.SimpleNamespace(names=reused))
    monkeypatch.setitem(sys.modules,'Numerical',types.ModuleType('Numerical'))
    monkeypatch.setitem(sys.modules,'Numerical.fitting',types.ModuleType('Numerical.fitting'))
    monkeypatch.setitem(sys.modules,'Numerical.core',types.ModuleType('Numerical.core'))
    monkeypatch.setitem(sys.modules,'Numerical.fitting.ScanCLI',fake_scan)
    monkeypatch.setitem(sys.modules,'Numerical.core.StateVector',fake_vector)
    module.prepare_state(original,entry,reused,values)
    assert captured[0]['base_state']['ordinary'][reused[0]]=={'real':-0.015,'imag':0.0}
    assert original['base_state']['ordinary'][reused[0]]=={'real':0.0,'imag':0.0}


def test_backward_integrator_for_zero_and_six_new_complex_slots():
    for n in (0,6):
        vec=np.array([0.35,0.5] + [-0.01,0.0]*n)
        solution=integrate_vector(vec,lambda t,y:0.18*y,n_points=8)
        assert solution.success
        np.testing.assert_allclose(solution.y[:,-1],vec*(STOP/START)**0.18,rtol=1e-8,atol=1e-10)


def test_wolfram_export_uses_projected_native_quartic_betas_without_headless_check():
    path=Path(__file__).resolve().parents[2]/'scripts/ExportGeneralT3AllModelsUVOneLoop.wl'
    s=path.read_text(encoding='utf-8')
    assert s.index('Needs["RGBeta`"]') < s.index('T3S31Export[] :=')
    assert 'BetaTerm[Quartic,1]' in s
    assert 'projectedBase=Association@Table' in s
    assert 'T3RGBetaOneLoopBetas[]' in s
    assert 'Check[Catch[T3S31Export[]' not in s
    assert 'T3GeneralQuarticReusedCount' in s


def test_runner_only_targets_15_and_summary_does_not_fake_completion():
    root=Path(__file__).resolve().parents[2]
    s=(root/'scripts/RunGeneralT3CoupledUVRunning.ps1').read_text(encoding='utf-8')
    assert 'ExportGeneralT3AllModelsUVOneLoop.wl' in s
    assert 'CoupledUVRunning --model' in s
    assert 'CoupledUVRunning --summary' in s
    assert '$LASTEXITCODE -ne 0' in s


def test_summary_fails_closed_without_saved_results(tmp_path,monkeypatch):
    from scripts import CoupledUVRunning as module
    reg=make_registry()
    fake_registry=tmp_path/'registry.json'
    fake_registry.write_text(json.dumps(reg),encoding='utf-8')
    monkeypatch.setattr(module,'REGISTRY',fake_registry)
    with pytest.raises(ValueError, match='Missing input'):
        summarize(folder=tmp_path)


def test_stage30_user_uploaded_beta_export_matches_input_convention_if_available():
    # Optional local fixture; tests cannot substitute for live native RGBeta.
    path=Path('/mnt/data/stage30_A_p0_full_uv_betas.json')
    if not path.is_file(): pytest.skip('No user-provided Stage30 native beta export')
    d=json.loads(path.read_text(encoding='utf-8'))
    assert d['status']=='Success'
    assert len(d['new_quartic_symbols'])==5
    assert len(d['base_beta_symbols'])==20
    assert d['metadata']['ReportBetaConvention']==BETA_CONVENTION


def test_wolfram_registry_empty_symbol_lists_preserved():
    """Wolfram Lookup[{},key] yields Missing, not the expected empty list.

    Class A_-2 has five new couplings but no reused couplings; A_-4 has
    neither. Both symbol lists must be mapped record-by-record.
    """
    path = Path(__file__).resolve().parents[2] / 'scripts/ExportGeneralT3AllModelsUVOneLoop.wl'
    source = path.read_text(encoding='utf-8')
    assert 'newRows=Lookup[entry,"new_rgbeta_symbols",Missing["new_rgbeta_symbols"]]' in source
    assert 'reusedRows=Lookup[entry,"preexisting_rgbeta_symbols",Missing["preexisting_rgbeta_symbols"]]' in source
    assert 'names=(Lookup[#,"name",Missing["name"]]& /@ newRows)' in source
    assert 'reused=(Lookup[#,"existing_name",Missing["existing_name"]]& /@ reusedRows)' in source
    assert '!ListQ[newRows]||!ListQ[reusedRows]' in source
    assert '!AllTrue[Join[newRows,reusedRows],AssociationQ]' in source
    assert '!AllTrue[Join[names,reused],StringQ]' in source
    assert 'Lookup[Lookup[entry,"new_rgbeta_symbols",{}],"name"]' not in source
    assert 'Lookup[Lookup[entry,"preexisting_rgbeta_symbols",{}],"existing_name"]' not in source


def test_stage31_empty_inventory_models_are_in_real_registry():
    """Read the actual project registry when present; no mock of Wolfram runs."""
    source = Path(__file__).resolve().parents[2] / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    if not source.exists():
        pytest.skip('Current repository registry is not available in this test tree')
    models = json.loads(source.read_text(encoding='utf-8'))['models']
    for name in ('A_-2', 'A_-4', 'D_+2', 'E_-4'):
        assert models[name]['preexisting_rgbeta_symbols'] == []
    for name in ('A_-4', 'D_+2', 'E_-4'):
        assert models[name]['new_rgbeta_symbols'] == []
