"""Regression tests for the experimental hard-boundary and RGBeta pilots."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pytest

from Numerical.experimental.GeneralT3C5Pilot import _safe_kernel, load_kernels, evaluate_response
from Numerical.experimental.GeneralT3RGSourcePilot import one_step_source

L1 = "lambdaT3ExtraQ020200I1"
L2 = "lambdaT3ExtraQ020020I1"


def mock_atom_file(tmp_path: Path, model='B_m1') -> Path:
    atoms = []
    for name, yy in ((L1, 'y1'), (L2, 'y2')):
        c = f"Bar[Coupling[{yy}, {{Index[d$$1, Flavor], Index[d$$1, NFlavor]}}, 0]]"
        expr = f"2*hbar*{c}^2*(1+Log[μbar2/Coupling[MF,{{}},0]^2])/Coupling[MF,{{}},0]"
        atoms.append({"coupling_atom":f"Bar[Coupling[{name},{{}},0]]",
                      "linear_coefficient_inputform":expr,
                      "coefficient_formally_zero":False,
                      "depends_on_other_extra_atoms":False})
    p = tmp_path/'atoms.json'
    p.write_text(json.dumps({"model":model,"atoms":atoms,
                             "response_reconstruction_residual_proven_zero":True,
                             "response_at_zero_proven_zero":True,
                             "polynomial_in_extra_atoms":True}))
    return p


def mock_config():
    return {"MF":[1000.,2000.,3000.],
            "y1_real":(np.eye(3)*1e-5).tolist(), "y1_imag":np.zeros((3,3)).tolist(),
            "y2_real":(np.eye(3)*2e-5).tolist(), "y2_imag":np.zeros((3,3)).tolist(),
            "matchete_extra_couplings":{
                L1:{"real":0.01,"imag":0}, L2:{"real":0.02,"imag":0}}}


def test_kernels_load_and_match_scale(tmp_path):
    kernels=load_kernels(mock_atom_file(tmp_path),'B_m1')
    assert len(kernels)==2
    assert {k.yukawa for k in kernels}=={'y1','y2'}
    assert all(str(k.expression)=='2*hbar/MF' for k in kernels)


def test_hard_response_matrix(tmp_path):
    kernels=load_kernels(mock_atom_file(tmp_path),'C_m1') if False else load_kernels(mock_atom_file(tmp_path),'B_m1')
    result=evaluate_response(kernels,mock_config())
    m=np.asarray(result['delta_c5']['real'])
    assert m.shape==(3,3)
    assert np.allclose(m,np.diag(np.diag(m)))
    assert m[0,0]>0
    expected=(2/(16*np.pi**2)/1000)*(0.01*(1e-5)**2+0.02*(2e-5)**2)
    assert np.isclose(m[0,0],expected)


def test_zero_extra_responds_zero(tmp_path):
    pars=mock_config()
    for value in pars['matchete_extra_couplings'].values(): value.update(real=0.,imag=0.)
    result=evaluate_response(load_kernels(mock_atom_file(tmp_path),'B_m1'),pars)
    assert result['delta_c5']['frobenius_norm']==0


def test_reject_poles():
    with pytest.raises(ValueError,match='pole'):
        _safe_kernel('hbar*Bar[Coupling[y1,{},0]]^2/(ϵ*Coupling[MF,{},0])','y1')


def test_reject_mixed_yukawa():
    with pytest.raises(ValueError,match='other Yukawa'):
        _safe_kernel('hbar*Bar[Coupling[y1,{},0]]*Bar[Coupling[y2,{},0]]/Coupling[MF,{},0]','y1')


def test_reject_unrecognised_matchete_symbols():
    with pytest.raises(ValueError,match='Unsupported symbols'):
        _safe_kernel('hbar*Bar[Coupling[y1,{},0]]^2*z/Coupling[MF,{},0]','y1')


def test_reject_unproven_artifact(tmp_path):
    path=mock_atom_file(tmp_path)
    data=json.loads(path.read_text())
    data['polynomial_in_extra_atoms']=False
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='polynomial_in_extra_atoms'):
        load_kernels(path,'B_m1')


def mock_beta_proof(tmp_path):
    p=tmp_path/'rgproof.json'
    rows=[]
    for name in ('B_m1','C_m1'):
        rows.append({'model':name,'status':'Analyzed','couplings':[
            {'coupling':c,'status':'Analyzed', 'nonzero_formal_source':True,
             'source_independent_of_special_quartics':True,
             'source_at_zero_special_quartics_inputform':'Bar[lambdaT3]*Tr[y1 . Trans[Bar[y2]]]'}
            for c in ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar')]})
    p.write_text(json.dumps({'status':'OneLoopSourcesVerified','results':rows}))
    return p


def test_rgbeta_source_numerical(tmp_path):
    pars=mock_config()
    pars['rgbeta_lambdaT3']={'real':.08,'imag':.02}
    pars['scale_uv_gev']=10000
    pars['scale_ir_gev']=1000
    result=one_step_source(mock_beta_proof(tmp_path),'B_m1',pars)
    assert result['status']=='ExperimentalFrozenBackgroundLeadingLog'
    val=result['results']['lambdaHHdagS2S2']['delta_lambda_leading_log']
    assert val['real']<0
    assert val['imag']>0


def test_rgbeta_scale_equal(tmp_path):
    pars=mock_config()
    pars.update(rgbeta_lambdaT3={'real':.08,'imag':0.0},scale_uv_gev=1000,scale_ir_gev=1000)
    result=one_step_source(mock_beta_proof(tmp_path),'C_m1',pars)
    assert all(row['delta_lambda_leading_log']['real']==0 for row in result['results'].values())


def test_rgbeta_reject_unverified(tmp_path):
    p=mock_beta_proof(tmp_path)
    d=json.loads(p.read_text());d['status']='ReviewRequired';p.write_text(json.dumps(d))
    pars=mock_config();pars.update(rgbeta_lambdaT3={'real':.08,'imag':0},scale_uv_gev=1000,scale_ir_gev=900)
    with pytest.raises(ValueError,match='not verified'):
        one_step_source(p,'B_m1',pars)


def test_physical_scalar_log_is_not_erased():
    # The renormalisation log vanishes at mu=MF but the mass-ratio log survives.
    raw = (
        'hbar*Bar[Coupling[y1,{},0]]^2*'
        '(1+Log[μbar2/Coupling[MF,{},0]^2]'
        '+Log[Coupling[MF,{},0]^2/Coupling[MS1,{},0]^2])'
        '/Coupling[MF,{},0]'
    )
    kernel = _safe_kernel(raw, 'y1')
    assert kernel.has(__import__('sympy').log)
    assert 'MS1' in str(kernel)
    assert 'μbar2' not in str(kernel)


def test_scalar_mass_log_evaluation(tmp_path):
    path = mock_atom_file(tmp_path)
    obj = json.loads(path.read_text())
    first = obj['atoms'][0]
    first['linear_coefficient_inputform'] = (
        first['linear_coefficient_inputform'].replace(
            '1+Log[μbar2/Coupling[MF,{},0]^2]',
            '1+Log[μbar2/Coupling[MF,{},0]^2]+'
            'Log[Coupling[MF,{},0]^2/Coupling[MS1,{},0]^2]'
        )
    )
    path.write_text(json.dumps(obj))
    pars = mock_config()
    with pytest.raises(ValueError, match='provide physical MS1'):
        evaluate_response(load_kernels(path,'B_m1'), pars)
    pars['MS1'] = 800.0
    r = evaluate_response(load_kernels(path,'B_m1'), pars)
    assert r['finite_mass_ratio_logs_preserved'] is True
    assert r['delta_c5']['frobenius_norm'] > 0


def test_reject_dimensional_log():
    with pytest.raises(ValueError, match='not demonstrably dimensionless'):
        _safe_kernel('hbar*Bar[Coupling[y1,{},0]]^2*Log[Coupling[MF,{},0]]/Coupling[MF,{},0]','y1')


def test_reject_unrecognised_log_argument():
    with pytest.raises(ValueError, match='Unsupported symbols'):
        _safe_kernel('hbar*Bar[Coupling[y1,{},0]]^2*Log[z/Coupling[MF,{},0]^2]/Coupling[MF,{},0]','y1')


def test_scalar_mass_square_consistency():
    from Numerical.experimental.GeneralT3C5Pilot import _mass_substitutions
    import sympy as sp
    v=_mass_substitutions({'MS1':1200.0,'MS2':1500.0})
    assert v[sp.Symbol('mS1Sq')]==1200.0**2
    assert v[sp.Symbol('mS2Sq')]==1500.0**2
    with pytest.raises(ValueError,match='Inconsistent'):
        _mass_substitutions({'MS1':1200.0,'mS1Sq':1.0})
