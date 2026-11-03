"""Small independent algebraic tests; run in the full local repository."""
from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.running.GeneralT3SequentialLeadingLog import (
    MODELS, anchor_to_one_scale, decode, encode, geometry_contribution,
    general_extra_channels, number,
)


def test_all_16_models_explicit():
    assert len(MODELS) == len(set(MODELS)) == 16
    assert {'B_-1','C_-1','E_+2','A_-4'} <= set(MODELS)


def test_real_imag_fourth_derivative_independent_contributions():
    rows_r=[{'H1_H2_S1_S2_0based':[0,1,1,2], 'value_exact':'-sqrt(2)'}]
    rows_i=[{'H1_H2_S1_S2_0based':[0,1,1,2], 'value_exact':'1/2'}]
    q=geometry_contribution(n=8,real_rows=rows_r,imag_rows=rows_i,
                            z=2+3j,label='test')
    expected=-2*np.sqrt(2)+1.5
    assert q.shape==(4,4,8,8)
    assert np.isclose(q[0,1,1,2],expected)
    assert np.isclose(q[1,0,2,1],expected)
    assert np.isclose(q[0,1,2,1],expected)
    assert np.count_nonzero(q)==4


def test_full_sector_extras_with_real_complex_and_exclusion():
    z={'counts_H_BarH_S1_BarS1_S2_BarS2':[1,1,1,1,0,0],
       'source_registry_symbols':[{'rgbeta_symbol':'lambdaT3GeneralQ1'},
                                  {'rgbeta_symbol':'lambdaT3GeneralQ2'}],
       'exact_reference_real_component_tensor_basis':[
          {'id':'r1','Re_c_unit_tensor_components':[
              {'H1_H2_S1_S2_0based':[0,0,0,0],'value_exact':'4'}],
           'Im_c_unit_tensor_components':[]},
          {'id':'r2','Re_c_unit_tensor_components':[],
           'Im_c_unit_tensor_components':[
              {'H1_H2_S1_S2_0based':[1,1,1,1],'value_exact':'-3'}]},
       ]}
    geometry={'monomial_sectors':[z], 'registry_HHSS_extra_symbols':2}
    registry={'models':{'B_-1':{'total_additional_quartics':2}}}
    state={'lambdaT3GeneralQ1':2,
           'lambdaT3GeneralQ2':{'real':0,'imag':5}}
    names,parts,info=general_extra_channels('B_-1',registry,state,geometry=geometry)
    assert names==['lambdaT3GeneralQ1','lambdaT3GeneralQ2']
    assert parts[0][0,0,0,0]==8
    assert parts[1][1,1,1,1]==-15
    names,parts,_=general_extra_channels('B_-1',registry,state,
        excluded={'lambdaT3GeneralQ2'},geometry=geometry)
    assert names==['lambdaT3GeneralQ1']
    with pytest.raises(ValueError,match='missing GENERAL scalar coupling'):
        general_extra_channels('B_-1',registry,{},geometry=geometry)


def test_one_scale_anchor_no_double_counting():
    full=np.array([[3e-12,1j*1e-12,0],[1j*1e-12,2e-12,0],[0,0,1e-12]])
    log=np.eye(3,dtype=complex)*0.3e-12
    on=anchor_to_one_scale(full,log,wilson_running_on=True)
    off=anchor_to_one_scale(full,log,wilson_running_on=False)
    np.testing.assert_allclose(decode(on['c5_at_scalar_threshold']),full)
    np.testing.assert_allclose(decode(off['c5_at_scalar_threshold']),full-log)
    np.testing.assert_allclose(decode(on['c5_residual_matching']),full-log)
    assert on['source'].startswith('one_scale_minus_log')


def test_decode_and_number_validation():
    a=np.eye(3,dtype=complex)*(0.3+0.7j)
    assert np.allclose(decode(encode(a)),a)
    assert number({'real':-2.,'imag':0.5},'x')==-2+0.5j
    with pytest.raises(ValueError,match='nonfinite'):
        number(float('nan'),'x')


def test_malformed_general_tensor_rejected():
    with pytest.raises(ValueError,match='bad derivative component'):
        geometry_contribution(n=8,real_rows=[
            {'H1_H2_S1_S2_0based':[0,3,8,0],'value_exact':'1'}],
            imag_rows=[],z=1,label='invalid')


def test_one_scale_wilson_ablation_exact_on_and_off(monkeypatch):
    from scripts import RunT3GeneralSequentialFitAblation as ab
    target=np.array([[4,0.1j,0],[0.1j,3,0],[0,0,2]],complex)*1e-12
    leading=np.eye(3,dtype=complex)*0.2e-12
    bubble=np.eye(3,dtype=complex)*0.4e-12
    fake_result={
       'c5_eft1_leading_log':encode(leading),
       'c5_scalar_finite_bubble':encode(bubble),
       'per_channel':[], 'quartic_provenance':{},
       'matching_scale_cancellation_residual':1e-14,
       'weinberg_beta_projection_residual':0.,
       'weinberg_bubble_projection_residual':0.,
    }
    from Numerical.running import GeneralT3SequentialLeadingLog as module
    monkeypatch.setattr(ab, 'response',lambda *a,**k:fake_result)
    monkeypatch.setattr(ab, 'low_energy',lambda state,c5,**k:{
        'c5_at_low_scale':encode(c5),
        'masses_ev_ascending':[1.,2.,3.],
        'prediction':{'sin2_theta13':0.02},'chi2':0.})
    fake=SimpleNamespace(couplings={'mS1Sq':1e6,'mS2Sq':1e6},
        provenance={'general_all':True})
    cases=ab._final_cases('A_-4',target,SimpleNamespace(),fake,None,
                          uv_on=True,eft_on=False,uv_nfev=20,eft_nfev=0,
                          spinor_factor=1.,vev=246.22)
    assert set(cases)=={'1000','1001','1010','1011'}
    for stage in ('0','1'):
        np.testing.assert_allclose(decode(cases['101'+stage]['c5_at_scalar_threshold']),target)
        np.testing.assert_allclose(decode(cases['100'+stage]['c5_at_scalar_threshold']),target-leading)
        np.testing.assert_allclose(decode(cases['101'+stage]['c5_bookkeeping_residual_matching']),
                                   target-leading-bubble)


def test_candidate_quartic_mapping_keeps_extra(monkeypatch):
    from Numerical.running.GeneralT3SequentialLeadingLog import native_aliases
    import types
    # The live registry generator provides this exact shape; fixture avoids
    # requesting local Matchete files in the pure algebra test.
    parent=types.ModuleType('Numerical.lo_certification')
    coupling=types.ModuleType('Numerical.lo_certification.couplings')
    child=types.ModuleType('Numerical.lo_certification.couplings.extra_quartic_source_v6')
    child.model_registry_quartics=lambda model,reg:[
        {'Matchete_generator_candidate_symbol':'lambdaT3ExtraQ110200I1',
         'rgbeta_symbol':'lambdaT3GeneralQ110200I1'}]
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification',parent)
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification.couplings',coupling)
    monkeypatch.setitem(sys.modules,'Numerical.lo_certification.couplings.extra_quartic_source_v6',child)
    names=['lambdaT3','lambdaT3ExtraQ110200I1']
    vals,alias=native_aliases('A_-2',names,{'lambdaT3':0.1,
            'lambdaT3GeneralQ110200I1':0.25}, {})
    assert vals['lambdaT3ExtraQ110200I1']==0.25
    assert alias['lambdaT3ExtraQ110200I1']=='lambdaT3GeneralQ110200I1'


def test_dynamic_tree_contracts_all_flavour_axes_without_dropping_terms(monkeypatch):
    import types
    from collections import defaultdict
    from Numerical.running.GeneralT3DynamicFermionWilson import _dense_operator
    p=types.ModuleType('Numerical.running.GeneralT3WilsonPolynomialBatch')
    def orbit(k):
        i,a,j,b,x,y=k
        return sorted({(i,a,j,b,x,y),(j,b,i,a,x,y),
                       (i,a,j,b,y,x),(j,b,i,a,y,x)})
    def polynomial(raw):
        grouped=defaultdict(complex)
        for key,z in raw.items():
            k=min(orbit(key))
            grouped[k]+=z
        return {k:v/len(orbit(k)) for k,v in grouped.items()},dict(grouped)
    p.source_term_alignment=lambda op,row:{'flavour_axis_permutation_for_lepton_order':[1,0]}
    p.canonical_polynomial=polynomial
    p.distinct_orbit=orbit
    b=types.ModuleType('Numerical.running.GeneralT3OrderedTensor')
    b.scalar_basis=lambda model:SimpleNamespace(total_dimension=12)
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3WilsonPolynomialBatch',p)
    monkeypatch.setitem(sys.modules,'Numerical.running.GeneralT3OrderedTensor',b)
    m=np.zeros((3,3),complex);m[0,1]=2j
    op={'source_index':7,'components':[{'lepton_su2':[1,2],
         'scalar_real':[5,6],
         'ordered_flavour_coefficient_gev_inv':encode(m)}]}
    by,indices,meta=_dense_operator([op],[{'source_index':7}],model='A_-4')
    assert indices==(7,)
    assert by.shape==(1,3,3,2,2,8,8)
    # Matchete lepton axes [1,0] map input flavour [p,q]=[1,2] -> [2,1].
    assert by[0,1,0,0,1,0,1]==0.5j
    assert by[0,0,1,1,0,1,0]==0.5j
    assert np.count_nonzero(by)==4
    np.testing.assert_allclose(by[0],by[0].transpose(1,0,3,2,4,5))
    np.testing.assert_allclose(by[0],by[0].swapaxes(4,5))
    assert meta['operator_symmetrisation'].startswith('exact')


def test_mixed_reference_flavour_lift_uses_matcher_axis_order(monkeypatch):
    """Regression for B/C: reference tensors have (flav,flav,SU2,SU2,S,S)."""
    import types
    from Numerical.running import GeneralT3DynamicFermionWilson as dyn

    fake_bridge = types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')

    class SourceBoundary:
        def __init__(self, model, path, hash, mu_f_gev, raw, by_source,
                     source_indices, labels, normalization_certified, source_status):
            self.model = model
            self.path = path
            self.hash = hash
            self.mu_f_gev = mu_f_gev
            self.raw = raw
            self.by_source = by_source
            self.source_indices = source_indices
            self.labels = labels
            self.normalization_certified = normalization_certified
            self.source_status = source_status

    fake_bridge.SourceBoundary = SourceBoundary
    monkeypatch.setitem(sys.modules, 'Numerical.joint_scalar_matching.fermion_boundary_bridge', fake_bridge)
    monkeypatch.setattr(dyn, '_heavy_bilinear', lambda state: state.bilinear)

    reference = np.array([[1., 2j, 0.], [2j, 3., 0.], [0., 0., 4.]], complex)
    current = np.array([[3., 1j, 0.], [1j, 4., 0.], [0., 0., 5.]], complex)
    geometry = np.zeros((2, 2, 8, 8), complex)
    geometry[0, 0, 0, 0] = 1.2
    geometry[1, 1, 1, 1] = -0.3
    ref_tensor = np.einsum('pq,ijab->pqijab', reference, geometry)
    old = SourceBoundary('B_-1', Path('reference.json'), 'reference_hash', 1e5,
                         ref_tensor, np.stack([ref_tensor]), (7,),
                         tuple('S' for _ in range(8)), False, 'reference')
    result, details = dyn.mixed_reference_boundary(
        'B_-1', old, SimpleNamespace(bilinear=reference),
        SimpleNamespace(bilinear=current))
    expected = np.einsum('pq,ijab->pqijab', current, geometry)
    assert result.raw.shape == (3, 3, 2, 2, 8, 8)
    np.testing.assert_allclose(result.raw, expected, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(result.raw, result.raw.transpose(1, 0, 3, 2, 4, 5))
    np.testing.assert_allclose(result.raw, result.raw.swapaxes(4, 5))
    assert details['reference_bilinear_reconstruction_max_relative_residual'] < 1e-12


def test_dynamic_tensor_is_accepted_by_scalar_response_interface(monkeypatch):
    """The exact array layout must produce (3,3,2,2,4,4) after HHSS contraction."""
    import types
    from Numerical.running.GeneralT3SequentialLeadingLog import TensorBundle, response
    from Numerical.running import GeneralT3DynamicFermionWilson as dyn
    from collections import defaultdict
    polynomial = types.ModuleType('Numerical.running.GeneralT3WilsonPolynomialBatch')

    def orbit(k):
        i, p, j, q, a, b = k
        return sorted({(i, p, j, q, a, b), (j, q, i, p, a, b),
                       (i, p, j, q, b, a), (j, q, i, p, b, a)})

    def symmetrise(raw):
        sums = defaultdict(complex)
        for key, value in raw.items():
            sums[min(orbit(key))] += value
        return {key: z / len(orbit(key)) for key, z in sums.items()}, dict(sums)

    polynomial.source_term_alignment = lambda *args: {
        'flavour_axis_permutation_for_lepton_order': [0, 1]}
    polynomial.canonical_polynomial = symmetrise
    polynomial.distinct_orbit = orbit
    basis = types.ModuleType('Numerical.running.GeneralT3OrderedTensor')
    basis.scalar_basis = lambda model: SimpleNamespace(total_dimension=12)
    monkeypatch.setitem(sys.modules, 'Numerical.running.GeneralT3WilsonPolynomialBatch', polynomial)
    monkeypatch.setitem(sys.modules, 'Numerical.running.GeneralT3OrderedTensor', basis)

    m = np.eye(3, dtype=complex) * 1e-5
    op = {'source_index': 3, 'components': [
        {'lepton_su2': [1, 1], 'scalar_real': [5, 5],
         'ordered_flavour_coefficient_gev_inv': encode(m)}]}
    by_source, indexes, _ = dyn._dense_operator([op], [{'source_index': 3}], 'A_-4')
    raw = by_source.sum(axis=0)
    fake = SimpleNamespace(model='A_-4', raw=raw, by_source=by_source,
                           source_indices=indexes, mu_f_gev=1e5)
    quartic = np.zeros((4, 4, 8, 8))
    quartic[0, 0, 0, 0] = 0.5
    bundle = TensorBundle('A_-4', np.eye(8) * 1e6, quartic,
                          quartic[None, ...], ('quartic',), {})
    large = types.ModuleType('Numerical.joint_scalar_matching.large_batch')

    def scalar_matcher(C, q, mass2, mu_f, mu_s):
        assert C.shape == (3, 3, 2, 2, 8, 8)
        amp = np.einsum('fgijab,ABab->fgijAB', C, q)
        assert amp.shape == (3, 3, 2, 2, 4, 4)
        projected = amp[:, :, 0, 0, 0, 0]
        c5 = 0.5 * (projected + projected.T)
        return {'beta': c5, 'running': c5 * 0.1,
                'bubble': c5 * 0.2, 'sum': c5 * 0.3,
                'beta_residual': 0., 'bubble_residual': 0.}

    large._response = scalar_matcher
    monkeypatch.setitem(sys.modules, 'Numerical.joint_scalar_matching.large_batch', large)
    result = response(fake, bundle, mu_f=1e5, mu_s=1e3, spinor_factor=1.)
    assert result['model'] == 'A_-4'
    assert len(result['per_channel']) == 2  # quartic + fermion-tree source
    assert decode(result['c5_eft1_leading_log']).shape == (3, 3)


def test_mixed_boundary_fractions_original_and_continuation():
    from scripts.RunT3GeneralSequentialFitAblation import mixed_fit_fractions
    original={'extra_quartic_fraction':0.1,
              'extra_mixed_scalar_bilinear_fraction':0.01}
    continuation={'matching_and_rge_provenance':{
                  'method':'mixed-scalar-conditional-one-scale',
                  'quartic_fraction':0.1,'mu2_fraction':0.01,
                  'modified_scalar_boundary':False}}
    assert mixed_fit_fractions(original)==(0.1,0.01)
    assert mixed_fit_fractions(continuation)==(0.1,0.01)
    assert mixed_fit_fractions({**original,**continuation})==(0.1,0.01)


def test_mixed_boundary_fractions_never_silently_invented():
    from scripts.RunT3GeneralSequentialFitAblation import mixed_fit_fractions
    with pytest.raises(ValueError,match='missing mixed-fit boundary fraction'):
        mixed_fit_fractions({'matching_and_rge_provenance':{}})
    with pytest.raises(ValueError,match='conflicting saved mixed-fit boundary fractions'):
        mixed_fit_fractions({'extra_quartic_fraction':0.2,
            'extra_mixed_scalar_bilinear_fraction':0.01,
            'matching_and_rge_provenance':{'quartic_fraction':0.1,'mu2_fraction':0.01}})
    with pytest.raises(ValueError,match='continuation changed scalar boundary'):
        mixed_fit_fractions({'matching_and_rge_provenance':{
            'quartic_fraction':0.1,'mu2_fraction':0.01,'modified_scalar_boundary':True}})


def test_resume_flag_exposed():
    from scripts.RunT3GeneralSequentialFitAblation import parse_args
    args=parse_args(['--source-run','thesis16_merged_v1','--run-id','general_seq_05',
                     '--spinor-factor','1','--resume'])
    assert args.resume is True
    assert args.source_run=='thesis16_merged_v1'


def test_resume_reuses_15_saved_models_and_retries_only_mixed_failure(monkeypatch, tmp_path):
    """Exercise the existing-run path without invoking symbolic or numerical physics."""
    import json
    import types
    from scripts import RunT3GeneralSequentialFitAblation as ab
    from scripts.RunT3GeneralSequentialLeadingLog import tag

    output = tmp_path / 'ablation'
    saved_run = output / 'general_seq_05'
    saved_run.mkdir(parents=True)
    fitroot = tmp_path / 'fits'
    (fitroot/'thesis16_merged_v1').mkdir(parents=True)
    (fitroot/'thesis16_merged_v1'/'settings.json').write_text(
        json.dumps({'canonical_run_id':'full16_general_v1'}),encoding='utf-8')
    registry = tmp_path / 'registry.json'
    registry.write_text('{}',encoding='utf-8')
    previous={'format':'t3-full-general-16_U_E_W_S_ablation-v1',
              'requested_models':list(MODELS),'spinor_factor_assumed':1.,
              'computed':15}
    ab.save(saved_run/'summary.json',previous)
    cases={f'{i:04b}':{} for i in range(16)}
    prior={}
    for model in MODELS:
        if model=='C_-1':
            continue
        record={'model':model,'status':'Computed','cases':cases,
                'provenance':{'test':'unchanged'}}
        path=saved_run/tag(model)/'cases.json'
        ab.save(path,record)
        prior[model]=path.read_bytes()
    ab.save(saved_run/tag('C_-1')/'error.json',
            {'model':'C_-1','status':'Blocked','error':'KeyError: extra_quartic_fraction'})

    general=types.ModuleType('scripts.RunT3GeneralAll16FromReference')
    general.OUTPUT=fitroot
    bridge=types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')
    bridge.load_release_manifest=lambda _: {m:tmp_path/m for m in MODELS}
    bridge.load_source=lambda _,expected_model:object()
    q=types.ModuleType('Numerical.joint_scalar_matching.source_quartic_v4')
    q.source_seed_path=lambda *_:None
    monkeypatch.setitem(sys.modules,general.__name__,general)
    monkeypatch.setitem(sys.modules,bridge.__name__,bridge)
    monkeypatch.setitem(sys.modules,q.__name__,q)
    monkeypatch.setattr(ab,'OUT',output)
    monkeypatch.setattr(ab,'preflight',lambda *args,**kwargs:{'ready':True,'missing':{},'issues':{}})
    monkeypatch.setattr(ab,'_load_channels',lambda *args,**kwargs:object())
    called=[]
    def fake_mixed(model,*args,**kwargs):
        called.append(model)
        return cases, {'test':'new C(−1) result'}
    monkeypatch.setattr(ab,'_mixed',fake_mixed)
    monkeypatch.setattr(ab,'write_csv',lambda *args,**kwargs:None)
    rc=ab.main(['--source-run','thesis16_merged_v1','--run-id','general_seq_05',
                '--spinor-factor','1','--registry',str(registry),'--resume'])
    assert rc==0
    assert called==['C_-1']
    assert ab.load(saved_run/'summary.json')['computed']==16
    assert ab.load(saved_run/tag('C_-1')/'cases.json')['status']=='Computed'
    assert not (saved_run/tag('C_-1')/'error.json').exists()
    for model,content in prior.items():
        assert (saved_run/tag(model)/'cases.json').read_bytes()==content
