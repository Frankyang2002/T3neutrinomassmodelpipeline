"""Offline mathematical/regression checks for the optional quadratic batch.

Tests the standalone scalar Hessian without requiring a Wolfram front end;
when installed in the repository it also checks the exact runner CLI contract.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
PATH=ROOT/'Numerical/experimental/GeneralT3ScalarQuadraticMasses.py'
spec=importlib.util.spec_from_file_location('t3_candidate_quadratic',PATH)
assert spec and spec.loader
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

KINDS={
 'A_+0':'self_s1','A_-2':'self_s2',
 'B_-1':'mixed','C_-1':'mixed',
 'D_+0':'self_s1','D_-2':'self_s2',
 'E_+0':'self_s1','E_-2':'self_s2',
}
ALL=['A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
     'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
     'E_-4','E_-2','E_+0','E_+2']


def fixture_manifest():
    models={}
    patterns={'self_s1':[0,0,0,2,0,0],
              'self_s2':[0,0,0,0,0,2],
              'mixed':[0,0,0,1,0,1]}
    for name in ALL:
        kind=KINDS.get(name)
        entries=[{'counts':patterns[kind], 'invariants':1}] if kind else []
        models[name]={'quadratic_structures_not_enabled':entries,
                      'quadratic_complex_couplings_not_enabled':len(entries)}
    return {'model_count':16,'models':models}


@pytest.mark.parametrize('model',ALL)
def test_exact_16_model_quadratic_inventory(model):
    m=fixture_manifest()
    assert module.quadratic_kind(m,model)==KINDS.get(model,'none')
    result=module.quadratic_mass_hessian(m,model,2000000.0,3000000.0)
    expected=2*sum(module.REP[model[0]][:2])
    assert len(result['mass_squared_eigenvalues_gev2'])==expected
    assert result['analytic_spectrum_verified'] is True
    assert result['matchete_mass_basis_verified'] is False
    assert result['included_in_threshold_matching'] is False
    assert np.allclose(np.asarray(result['real_scalar_hessian_gev2']),
                       np.diag([2000000.0]*module.REP[model[0]][0] +
                               [3000000.0]*module.REP[model[0]][1] +
                               [2000000.0]*module.REP[model[0]][0] +
                               [3000000.0]*module.REP[model[0]][1]))


@pytest.mark.parametrize('model',sorted(KINDS))
@pytest.mark.parametrize('mu',[complex(82000,0),complex(0,71000),complex(-41000,27000)])
def test_all_eight_complex_quadratic_spectra_against_analytic_and_potential(model,mu):
    manifest=fixture_manifest()
    d1,d2,_=module.REP[model[0]]
    result=module.quadratic_mass_hessian(manifest,model,1.21e6,2.45e6,mu)
    H=np.asarray(result['real_scalar_hessian_gev2'])
    eigen=np.linalg.eigvalsh(H)
    analytic=module.analytic_mass_squared(manifest,model,1.21e6,2.45e6,mu)
    assert np.allclose(eigen,analytic,atol=1e-8,rtol=1e-13)
    assert np.allclose(H,H.T)
    # Independent direct field-space contraction, not a second Hessian formula.
    rng=np.random.default_rng(2026)
    x=rng.normal(size=d1+d2)
    y=rng.normal(size=d1+d2)
    z=(x+1j*y)/np.sqrt(2)
    s1,s2=z[:d1],z[d1:]
    V=1.21e6*np.vdot(s1,s1).real+2.45e6*np.vdot(s2,s2).real
    kind=KINDS[model]
    if kind=='self_s1':
        extra=mu*np.conj(s1).T@module._metric(d1)@np.conj(s1)
    elif kind=='self_s2':
        extra=mu*np.conj(s2).T@module._metric(d2)@np.conj(s2)
    else:
        extra=mu*np.conj(s1).T@module._metric(d1)@np.conj(s2)
    V+=2*extra.real
    R=np.concatenate([x,y])
    from_hessian=0.5*R@H@R
    assert np.isclose(V,from_hessian,atol=1e-8,rtol=1e-12)


def test_zero_splitting_and_mass_phase_invariance():
    manifest=fixture_manifest()
    for model in KINDS:
        a=module.analytic_mass_squared(manifest,model,4e6,5e6,1e5+0j)
        b=module.analytic_mass_squared(manifest,model,4e6,5e6,0+1e5j)
        assert np.allclose(a,b)


def test_nonzero_quadratic_forbidden_on_other_models():
    manifest=fixture_manifest()
    for model in ALL:
        if model not in KINDS:
            with pytest.raises(ValueError):
                module.quadratic_mass_hessian(manifest,model,4e6,5e6,1+0j)


def test_self_bilinear_tachyon_bound():
    manifest=fixture_manifest()
    val=module.quadratic_mass_hessian(manifest,'D_+0',1e6,2e6,0.6e6)
    assert val['tachyonic'] is True
    assert val['minimum_mass_squared_gev2']<0
    no=module.quadratic_mass_hessian(manifest,'D_+0',1e6,2e6,0.1e6)
    assert no['tachyonic'] is False


def test_mixed_bilinear_analytical_masses():
    manifest=fixture_manifest()
    m1,m2,mu=1.2e6,2e6,2e5+3e5j
    eig=module.analytic_mass_squared(manifest,'B_-1',m1,m2,mu)
    delta=np.hypot((m1-m2)/2,abs(mu))
    assert np.allclose(eig[:4],(m1+m2)/2-delta)
    assert np.allclose(eig[4:],(m1+m2)/2+delta)


def test_tampered_singlet_manifest_rejected():
    manifest=fixture_manifest()
    manifest['models']['B_-1']['quadratic_structures_not_enabled'][0]['invariants']=2
    with pytest.raises(ValueError):
        module.quadratic_kind(manifest,'B_-1')


def test_repository_manifest_if_present():
    file=ROOT/'docs/general_scalar_trial_operator_manifest.json'
    if not file.exists():
        pytest.skip('Local full repository manifest not present in packaging environment')
    manifest=module.load_manifest(file)
    assert {m for m in manifest['models'] if module.quadratic_kind(manifest,m)!='none'}==set(KINDS)


def test_quadratic_cli_command_models_in_repository():
    # The full CLI is meaningful only in the user's complete repository.
    src=ROOT/'scripts/RunGeneralT3QuadraticEFT1.py'
    content=src.read_text(encoding='utf-8')
    assert "alpha_token = f'm{-rep.alpha}' if rep.alpha < 0 else str(rep.alpha)" in content
    assert "stage='EFT1'" in content
    assert "'physical_matchete_rgbeta_mapping_certified':False" in content
    assert "'wilson_coefficient_matching_performed':False" in content


def test_independent_quadratic_enumeration_over_sixteen_models(tmp_path):
    sys.path.insert(0,str(ROOT))
    from scripts.AuditGeneralT3QuadraticBasis import audit
    path=tmp_path/'manifest.json'
    path.write_text(json.dumps(fixture_manifest()),encoding='utf-8')
    report=audit(manifest_path=path,output=tmp_path/'audit.json')
    assert report['models_checked']==16
    assert report['quadratic_slots']==8
    assert report['neutral_self_bilinear_slots']==6
    assert report['opposite_hypercharge_doublet_mixing_slots']==2
    assert report['claim_of_physical_basis_redundancy'] is False
    assert (tmp_path/'audit.json').exists()


def test_renormalisable_mass_beta_export_tracks_all_required_terms():
    text=(ROOT/'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl').read_text()
    for coupling in ('muT3GeneralS1Sq','muT3GeneralS2Sq','muT3GeneralS1S2'):
        assert f'AddScalarMass[{coupling}' in text
    assert 'BetaTerm[Quartic,1]' in text
    assert 'BetaTerm[mS1Sq,1]' in text
    assert 'BetaTerm[mS2Sq,1]' in text
    assert 'BetaTerm[extraMassSymbol,1]' in text
    assert 'SelfConjugate -> False' in text
    assert 'PhysicalMatcheteRGBetaMassBasisCertified' in text
