"""Fail-closed contract tests for direct one-scale C5 -> SM+Weinberg endpoint."""
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.RunGeneralT3LowScaleOneScale import (
    STATUS, _complex_entry_matrix, _check_inputs, _masses, evaluate, fail,
)

class State:
    def __init__(self, *, gY, g2, g3, lambdaH, ye, yu, yd, K):
        self.gY=gY;self.g2=g2;self.g3=g3;self.lambdaH=lambdaH
        self.ye=ye;self.yu=yu;self.yd=yd;self.K=K
    def validated(self):
        assert np.asarray(self.K).shape == (3,3)
        return self

def test_entry_roundtrip_and_takagi():
    z=[{'real':float(i==j)*1e-13,'imag':0.} for i in range(3) for j in range(3)]
    k=_complex_entry_matrix([z[i:i+3] for i in (0,3,6)],'k')
    assert np.allclose(k,1e-13*np.eye(3))
    assert np.allclose(_masses(k,246.22),[3.03121442]*3,rtol=1e-7)

def test_entry_rejects_malformed():
    with pytest.raises(ValueError,match='expected 3x3'):_complex_entry_matrix([],'k')
    with pytest.raises(ValueError,match='malformed'):
        _complex_entry_matrix([[{'real':0.,'imag':0.}]*3]*2+[[{'real':0.,'imag':'x'}]*3],'k')

def fixture_files(tmp_path, model='B_-1'):
    slug='T3_dS1_2_dS2_2_dF_'+('1' if model=='B_-1' else '3')+'_alpha_m1'
    frozen=tmp_path/(slug+'.json')
    m=np.eye(3)*1e-4
    mat=lambda a:{'real':np.asarray(a).real.tolist(),'imag':np.asarray(a).imag.tolist()}
    sm={'gY':0.39,'g2':0.6,'g3':0.8,'lambdaH':0.25,
        'ye':mat(m),'yu':mat(m),'yd':mat(m)}
    config={'representation':{'d_s1':2,'d_s2':2,'d_f':1 if model=='B_-1' else 3,
                              'alpha':-1,'shared_scalar':False},
       'numerical':{'vev_gev':246.22,'rtol':1e-8,'atol':1e-11},
       'scales':{'mu_fermion_threshold_gev':100000.,'mu_scalar_threshold_gev':1000.,
                 'mu_low_gev':100.},
       'base_state':{'mu_gev':1e7,'sm':sm},
       'benchmark_search':{'enabled':False},'sensitivity':{'enabled':False}}
    frozen.write_text(json.dumps(config),encoding='utf-8')
    digest=hashlib.sha256(frozen.read_bytes()).hexdigest()
    rel=f'Numerical/benchmarks/fitted_reference/{slug}.json'
    manifest={'models':[{'model':slug,'config_path':rel,'sha256':digest}]+[
        {'model':f'dummy-{i}'} for i in range(15)]}
    mp=tmp_path/'benchmark_manifest.json'; mp.write_text(json.dumps(manifest),encoding='utf-8')
    coeff=1e-13*np.eye(3)
    source={'model':model,'status':'FrozenScalarDiagonalDirectOneScaleC5NotPhysicalEFT1',
            'frozen_scalar_aligned_one_scale_scenario':True,
            'mass_basis_compatible_with_frozen_reference':True,
            'source_mass_parameters_used_as_numeric_scalar_input':False,
            'n_fermion_flavours':3,'sequential_EFT1_performed':False,
            'physical_c5_certified':False,'scalar_masses_gev':[1000.,1142.857],
            'benchmark_config':{
                'format':'t3-frozen-scalar-diagonal-direct-C5-input-v1',
                'physical_input_certified':False,'sequential_EFT1':False,
                'vev_gev':246.22,'reference_provenance':{
                    'frozen_config_path':rel,'frozen_sha256':digest,
                    'quadratic_mixing_mu2_gev2':{'real':0.,'imag':0.},
                    'quadratic_mixing_source':'zero_by_explicit_default_scenario_not_in_fitted_reference',
                    'physical_EFT1_matching_performed':False}},
            'symmetrized_C5_gev_inv':[[{'real':float(coeff[i,j].real), 'imag':0.}
                                       for j in range(3)] for i in range(3)],
            'largest_takagi_mass_eV':_masses(coeff,246.22)[0]}
    sf=tmp_path/'source.json';sf.write_text(json.dumps(source),encoding='utf-8')
    return sf,frozen,mp,source,config

def fake_evolve(initial, mu_initial, mu_final, *,rtol,atol):
    assert mu_initial>mu_final and rtol==1e-8 and atol==1e-11
    factor=(mu_final/mu_initial)**0.01
    return SimpleNamespace(solver_success=True,gY=initial.gY,g2=initial.g2,g3=initial.g3,
        lambdaH=initial.lambdaH,ye=initial.ye,yu=initial.yu,yd=initial.yd,
        K=np.asarray(initial.K)*factor)

@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_endpoint_injects_one_scale_C5_after_background(tmp_path,model):
    source,frozen,manifest,_,_=fixture_files(tmp_path,model)
    record=evaluate(model,source,frozen,manifest,evolution=fake_evolve,state_cls=State)
    assert record['status']==STATUS
    assert record['fermion_first_EFT1_matching_performed'] is False
    assert record['renormalizable_BSM_RGE_included_above_scalar_threshold'] is False
    assert record['physical_C5_prediction_certified'] is False
    assert record['masses_eV_at_low_scale_descending'][0]<record['masses_eV_at_scalar_boundary_descending'][0]
    assert record['relative_largest_mass_change']==pytest.approx(10**(-0.01)-1)
    assert record['mu_scalar_boundary_gev']==1000.
    assert record['mu_low_gev']==100.

def test_manifest_hash_change_is_rejected(tmp_path):
    source,frozen,manifest,_,_=fixture_files(tmp_path)
    f=json.loads(frozen.read_text());f['numerical']['vev_gev']=123
    frozen.write_text(json.dumps(f),encoding='utf-8')
    with pytest.raises(ValueError,match='SHA256'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)

def test_source_status_fail_closed(tmp_path):
    source,frozen,manifest,orig,_=fixture_files(tmp_path)
    orig['physical_c5_certified']=True
    source.write_text(json.dumps(orig),encoding='utf-8')
    with pytest.raises(ValueError,match='source is not'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)

def test_source_reference_sha_fail_closed(tmp_path):
    source,frozen,manifest,orig,_=fixture_files(tmp_path)
    orig['benchmark_config']['reference_provenance']['frozen_sha256']='0'*64
    source.write_text(json.dumps(orig),encoding='utf-8')
    with pytest.raises(ValueError,match='provenance'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)

def test_nonsymmetric_C5_rejected(tmp_path):
    source,frozen,manifest,orig,_=fixture_files(tmp_path)
    orig['symmetrized_C5_gev_inv'][0][1]['real']=2e-13
    source.write_text(json.dumps(orig),encoding='utf-8')
    with pytest.raises(ValueError,match='not symmetric'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)

def test_invalid_scales_rejected(tmp_path):
    source,frozen,manifest,orig,config=fixture_files(tmp_path)
    config['scales']['mu_low_gev']=2000.
    frozen.write_text(json.dumps(config),encoding='utf-8')
    digest=hashlib.sha256(frozen.read_bytes()).hexdigest()
    man=json.loads(manifest.read_text()); man['models'][0]['sha256']=digest
    manifest.write_text(json.dumps(man),encoding='utf-8')
    orig['benchmark_config']['reference_provenance']['frozen_sha256']=digest
    source.write_text(json.dumps(orig),encoding='utf-8')
    with pytest.raises(ValueError,match='threshold ordering'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)

def test_background_must_preserve_zero_C5(tmp_path):
    source,frozen,manifest,_,_=fixture_files(tmp_path)
    def buggy(initial,a,b,**kwargs):
        obj=fake_evolve(initial,a,b,**kwargs)
        obj.K=np.asarray(obj.K)+1e-12
        return obj
    with pytest.raises(ValueError,match='unexpectedly generated C5'):
        evaluate('B_-1',source,frozen,manifest,evolution=buggy,state_cls=State)

def test_source_Takagi_must_match_coeff(tmp_path):
    source,frozen,manifest,orig,_=fixture_files(tmp_path)
    orig['largest_takagi_mass_eV']=99.
    source.write_text(json.dumps(orig),encoding='utf-8')
    with pytest.raises(ValueError,match='internally inconsistent'):
        evaluate('B_-1',source,frozen,manifest,evolution=fake_evolve,state_cls=State)
