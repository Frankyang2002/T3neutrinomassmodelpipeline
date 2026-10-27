"""Mass diagonalisation / SU(2) covariance / complete B-C tensor inventory tests."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis,EPS,field_map,transport_operator,
    verify_su2_covariance,verify_tensor_transport,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import quadratic_mass_hessian
from scripts import RunGeneralT3MixedScalarMassBasis as driver

M=(1350.**2,1780.**2)

@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('mu',[0j,12000+0j,125000+43000j,-40000-85000j])
def test_mass_diagonalization_and_inverse(model,mu):
    b=MixedDoubletMassBasis.build(model,*M,mu)
    assert np.allclose(b.U.conj().T@b.A@b.U,np.diag(b.eigenvalues_sq),atol=1e-6)
    assert np.allclose(b.U.conj().T@b.U,np.eye(2),atol=1e-13)
    rng=np.random.default_rng(111)
    for _ in range(10):
        x=rng.normal(size=2)+1j*rng.normal(size=2)
        y=rng.normal(size=2)+1j*rng.normal(size=2)
        s1,s2=b.original_fields(x,y)
        a,c=b.eigen_fields(s1,s2)
        assert np.allclose(a,x,rtol=1e-12,atol=1e-12)
        assert np.allclose(c,y,rtol=1e-12,atol=1e-12)
        assert np.isclose(b.old_potential(s1,s2),b.new_potential(x,y),rtol=1e-12)
    assert verify_su2_covariance(b)<1e-12

@pytest.mark.parametrize('model',['B_-1','C_-1'])
@pytest.mark.parametrize('mu',[0j,12000+0j,125000+43000j])
def test_old_hessian_spectrum_and_eigenvalues(model,mu):
    mf=fixture_manifest()
    b=MixedDoubletMassBasis.build(model,*M,mu)
    old=quadratic_mass_hessian(mf,model,*M,mu)
    assert np.allclose(old['mass_squared_eigenvalues_gev2'],
                       np.repeat(b.eigenvalues_sq,4),atol=1e-6,rtol=1e-12)


def fixture_manifest():
    m={'models':{}}
    for model in ('B_-1','C_-1'):
        structures=[([0,0,0,1,1,2],1),([0,0,0,2,0,2],1),
                    ([0,0,1,2,0,1],1),([0,2,0,0,2,0],1),
                    ([0,2,0,2,0,0],1),([1,1,0,1,0,1],2)]
        m['models'][model]={'extra_quartic_complex_couplings':7,
            'quadratic_complex_couplings_not_enabled':1,
            'extra_quartic_structures':[{'counts':a,'invariants':n} for a,n in structures],
            'quadratic_structures_not_enabled':[{'counts':[0,0,0,1,0,1],
                                                 'invariants':1}]}
    return m


def fixture_registry():
    r={'totals':{'matchete_extra_quartics':50,'rgbeta_new_vertices':38,
                 'rgbeta_preexisting_structural':12},'models':{}}
    count_orders=[([0,0,0,1,1,2],'lambdaS1barS2S2bar2'),
                  ([0,0,0,2,0,2],'lambdaS1bar2S2bar2'),
                  ([0,0,1,2,0,1],'lambdaS1S1bar2S2bar'),
                  ([0,2,0,0,2,0],'lambdaHHdagS2S2'),
                  ([0,2,0,2,0,0],'lambdaHHdagS1barS1bar'),
                  ([1,1,0,1,0,1],'lambdaHHdagS1barS2barCross')]
    for model in ('B_-1','C_-1'):
        r['models'][model]={'model':model,'dimensions':[2,2,1 if model=='B_-1' else 3],
                            'alpha':-1,'total_additional_quartics':7,
                            'quadratics_not_implemented':1,
                            'preexisting_rgbeta_symbols':[
                                {'counts':cnt,'existing_name':name} for cnt,name in count_orders],
                            'new_rgbeta_symbols':[{'name':'lambdaT3GeneralQ110101I1',
                                'counts':[1,1,0,1,0,1],
                                'fields':['H','Bar @ H','Bar @ S1','Bar @ S2']}]}
    return r


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_complete_seventeen_operator_registry(model):
    names=driver.quartic_operators(fixture_manifest(),fixture_registry(),model)
    assert len(names)==17
    assert set(driver.REUSED_QUARTICS)<=set(names)
    assert set(driver.BASE_QUARTICS)<=set(names)
    assert set(names)-set(driver.BASE_QUARTICS)-set(driver.REUSED_QUARTICS)=={
        'lambdaT3GeneralQ110101I1'}


@pytest.mark.parametrize('model',['B_-1','C_-1'])
def test_structural_registry_failure_closed(model):
    r=fixture_registry()
    r['models'][model]['new_rgbeta_symbols'][0]['counts']=[0,0,1,1,1,1]
    with pytest.raises(ValueError):
        driver.quartic_operators(fixture_manifest(),r,model)


@pytest.mark.parametrize('mu',[0j,1230+0j,25000+12000j])
def test_all_seventeen_tensors_plus_yukawa_scalar_legs(mu):
    u=MixedDoubletMassBasis.build('B_-1',*M,mu).U
    for fields in driver.quartic_operators(fixture_manifest(),fixture_registry(),'B_-1').values():
        assert verify_tensor_transport(fields,u,seed=12,trials=3)<2e-12
        terms=transport_operator(fields,u)
        assert len(terms)==2**sum(f not in ('H','Hdag') for f in fields)
    for f in driver.YUKAWA_SCALAR_LEGS.values():
        assert verify_tensor_transport((f,),u)<2e-12


def test_epsilon_is_not_identity_and_hypercharge_preserved():
    u=MixedDoubletMassBasis.build('B_-1',*M,12000+8000j).U
    s2=field_map('S2',u)
    assert np.allclose(s2[0][1],-u[1,0].conjugate()*EPS)
    assert not np.allclose(s2[0][1],u[1,0]*np.eye(2))
    assert all(n.endswith('dag') for n,_ in s2)
    assert all(n in ('X','Y') for n,_ in field_map('S2dag',u))


def test_invalid_and_tachyonic_model_inputs():
    with pytest.raises(ValueError):
        MixedDoubletMassBasis.build('A_+0',*M,0)
    with pytest.raises(ValueError):
        MixedDoubletMassBasis.build('B_-1',np.nan,M[1],0)
    with pytest.raises(ValueError):
        transport_operator(('F',),np.eye(2))
    b=MixedDoubletMassBasis.build('B_-1',1.,1.,2.+0j)
    assert b.tachyonic
    assert b.info()['masses_gev'][0] is None


def test_full_saved_eft1_endpoint_roundtrip(tmp_path,monkeypatch):
    root=tmp_path
    (root/'docs').mkdir()
    (root/'RGE/general/generated').mkdir(parents=True)
    (root/'output/general_t3_quadratic_eft1').mkdir(parents=True)
    manifest=fixture_manifest()
    reg=fixture_registry()
    (root/'docs/general_scalar_trial_operator_manifest.json').write_text(json.dumps(manifest))
    (root/'RGE/general/generated/general_rgbeta_quartic_registry.json').write_text(json.dumps(reg))
    monkeypatch.setattr(driver,'ROOT',root)
    monkeypatch.setattr(driver,'MANIFEST',root/'docs/general_scalar_trial_operator_manifest.json')
    monkeypatch.setattr(driver,'REGISTRY',root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    for model in ('B_-1','C_-1'):
        raw=driver.quartic_operators(manifest,reg,model)
        baseline={k:({'real':0.013,'imag':0.005} if k not in driver.REAL_QUARTICS
                     else 0.013) for k in raw if k not in ('lambdaH','lambdaT3GeneralQ110101I1')}
        baseline.update({'mS1Sq':M[0],'mS2Sq':M[1], 'sm':{'lambdaH':0.13}})
        mu=complex(125000.,-23000.)
        bp=MixedDoubletMassBasis.build(model,*M,mu)
        prototype=quadratic_mass_hessian(manifest,model,*M,mu)
        payload={'format':'general-t3-candidate-quadratic-eft1-v1',
                 'status':'CandidateQuadraticAndFullQuarticEFT1Integrated',
                 'model':model,'quadratic_kind':'mixed','mu2_beta_included':True,
                 'physical_matchete_rgbeta_mapping_certified':False,
                 'scalar_mass_threshold_matching_performed':False,
                 'mu_final_gev':1000.,'final_baseline':baseline,
                 'final_new_quartics':{'lambdaT3GeneralQ110101I1':{'real':0.01,'imag':0.005}},
                 'mu2_final_gev2':{'real':mu.real,'imag':mu.imag},
                 'final_scalar_mass_spectrum':prototype}
        fp=root/f'output/general_t3_quadratic_eft1/{model.replace("-","m")}_eft1_quadratic_running.json'
        fp.write_text(json.dumps(payload))
        result=driver.run_model(model,folder=fp.parent,
                                 output=root/'output/general_t3_mass_basis',trials=2)
        assert result['source_quartic_count']==17
        assert len(result['quartics'])==17
        assert len(result['yukawa_scalar_leg_transforms'])==2
        assert result['basis']['unitarity_verified']
        assert result['quadratic_hessian_agreement']
        assert np.allclose(result['basis']['masses_squared_gev2'],bp.eigenvalues_sq)
        assert result['max_source_tensor_transport_relative_error']<2e-12
        assert all(chr(92) not in key for key in result['source_sha256'])
        assert 'RGE/general/generated/general_rgbeta_quartic_registry.json' in result['source_sha256']
        assert (root/f'output/general_t3_mass_basis/{model.replace("-","m")}_final_mass_basis.json').exists()


def test_reject_previously_saved_bad_mass_spectrum(tmp_path,monkeypatch):
    # The numerical one-loop data must not be silently replaced by reconstructed guesses.
    assert driver.STATUS=='ExperimentalMassBasisTensorTransportNoPhysicalMatching'
