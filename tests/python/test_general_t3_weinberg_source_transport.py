"""Stage19 source-selected Weinberg projection tests (no production code changed)."""
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from Numerical.experimental import GeneralT3WeinbergProjection as mod


def checked_normalization():
    return {'status':'ProductionDirectC5NormalizationVerified'}


def checked_certificate(model='B_m1'):
    return {'status':'FullExtraWeinbergTensorProjectionVerified',
            'model':model,'independent_components_verified':30,'all_component_residuals_zero':True}


def source_payload(model='B_m1',scale=1e-26):
    k=np.array([[1+2j,2+0.5j,0.7-2j],[2+0.5j,-1j,3+2j],[.7-2j,3+2j,0.5j]])*scale
    t=mod._analytic_mixed_template()
    return {'status':mod.STAGE18_STATUS,'model':model,
            'extra_source_only':True,'eft1_wilson_seed_frozen_at_mu_F':True,
            'full_eft1_wilson_anomalous_dimension_implemented':False,
            'higgs_real_components':[
                {'higgs_real_indices':list(pair),
                 'integrated_delta_C_ijab_gev_inverse':{'real':(z*k).real.tolist(),
                                                      'imag':(z*k).imag.tolist()}}
                 for pair,z in t.items()]}, k


def test_template_matches_production_polynomial_and_derivative_factor():
    r=sp.symbols('r1:5',real=True)
    nu,e=sp.symbols('nu e')
    hp=(r[0]+sp.I*r[1])/sp.sqrt(2)
    h0=(r[2]+sp.I*r[3])/sp.sqrt(2)
    poly=(nu*h0-e*hp)**2
    for (a,b),ratio in mod._analytic_mixed_template().items():
        got=sp.diff(poly,nu,e,r[a-1],r[b-1])/4
        assert sp.simplify(got-sp.sympify(ratio))==0


def test_integration_projects_correct_complex_symmetric_matrix():
    payload,expected=source_payload()
    got,report=mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())
    assert np.allclose(got,expected,atol=1e-40,rtol=1e-14)
    assert report['projection_relative_residual']<1e-13
    assert report['complete_two_loop_prediction'] is False


def test_projection_zero_source_is_zero():
    payload,expected=source_payload(scale=0)
    got,report=mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())
    assert np.count_nonzero(got)==0
    assert report['projection_relative_residual']==0


def test_integration_rejects_wrong_relative_phase():
    payload,_=source_payload()
    c=payload['higgs_real_components'][1]['integrated_delta_C_ijab_gev_inverse']
    c['real'],c['imag']=[[-v for v in row] for row in c['real']],[[-v for v in row] for row in c['imag']]
    with pytest.raises(ValueError,match='phase/sign'):
        mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())


def test_integration_rejects_missing_complex_component():
    payload,_=source_payload()
    payload['higgs_real_components'].pop()
    with pytest.raises(ValueError,match='four'):
        mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())


def test_integration_rejects_unsupported_component():
    payload,_=source_payload()
    payload['higgs_real_components'][0]['higgs_real_indices']=[1,1]
    with pytest.raises(ValueError,match='support'):
        mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())


def test_integration_rejects_wrong_model_and_certificate():
    payload,_=source_payload()
    with pytest.raises(ValueError,match='Stage18'):
        mod.project_stage18_integral(payload,'C_m1',checked_certificate('C_m1'),checked_normalization())
    with pytest.raises(ValueError,match='certificates'):
        mod.project_stage18_integral(payload,'B_m1',checked_certificate('C_m1'),checked_normalization())


def test_integration_rejects_unrecognized_source_ode():
    payload,_=source_payload()
    payload['extra_source_only']=False
    with pytest.raises(ValueError,match='single-insertion'):
        mod.project_stage18_integral(payload,'B_m1',checked_certificate(),checked_normalization())


def test_production_normalization_positive_and_negative_b_c():
    for sign in (1,-1):
        prefactor=f'{sign*4}*sqrt(3)*lambdaT3/3' if sign==1 else '-4*lambdaT3/3'
        b=f'{prefactor}*conjugate(y1)*conjugate(y2)/MF'
        r={'status':'Success','weinberg_subspace_validation':{'matches_weinberg_subspace':True,'beta_kappa_16pi2':b}}
        f={'status':'Success','combined':{'ready_for_physical_majorana_numerics':True},
           'direct_running_contribution':{'manifestly_symmetric_under_pq':True,
                                          'raw_expression':f'({prefactor}*log(MS/MF))*C12[p,q]'}}
        assert mod.verify_production_c5_normalization(r,f)['direct_running_log_coefficient_equals_beta_kappa']


def test_production_normalization_rejects_factor_or_sign_error():
    r={'status':'Success','weinberg_subspace_validation':{'matches_weinberg_subspace':True,'beta_kappa_16pi2':'a*conjugate(y1)*conjugate(y2)/MF'}}
    f={'status':'Success','combined':{'ready_for_physical_majorana_numerics':True},
       'direct_running_contribution':{'manifestly_symmetric_under_pq':True,
                                      'raw_expression':'(-a*log(MS/MF))*C12[p,q]'}}
    with pytest.raises(ValueError,match='mismatch'):
        mod.verify_production_c5_normalization(r,f)


def test_production_normalization_rejects_unsuccessful_final_matching():
    with pytest.raises(ValueError,match='certificate'):
        mod.verify_production_c5_normalization({}, {})


def test_full_symbolic_certification_uses_all_thirty_components(monkeypatch):
    # A fast synthetic module-level test; the real Matchete seed is tested on Windows.
    from itertools import combinations_with_replacement
    template={}
    k=sp.Symbol('k')
    for i,j in combinations_with_replacement((1,2),2):
        for a,b in combinations_with_replacement((1,2,3,4),2):
            # Exact template is nonzero for precisely six (nu,nu)/(e,e) plus
            # four (nu,e) entries. Synthetic has only mixed support.
            if (i,j)==(1,2) and (a,b) in mod._analytic_mixed_template():
                template[(i,j,a,b)]=sp.nsimplify(mod._analytic_mixed_template()[(a,b)])
    monkeypatch.setattr(mod,'_production_template',lambda:template)
    monkeypatch.setattr(mod,'real_quartic_entries',lambda:{})
    monkeypatch.setattr(mod,'extra_scalar_pair_component',lambda quartic,w,i,j,a,b:sp.simplify(k*template.get((i,j,a,b),0)))
    import types
    fake=types.ModuleType('RGE.running.intermediate.ScalarOnlyTensorAdapters')
    fake.build_eft1_wilson_tensor=lambda *a,**kw:object()
    monkeypatch.setitem(sys.modules,'RGE.running.intermediate.ScalarOnlyTensorAdapters',fake)
    source={'status':'ExtraEFT1ScalarPairMixingDerived','model':'B_m1',
            'components':[{'higgs_real_indices':list(pair),
                           'beta_16pi2_expression':sp.sstr(k*val)} for pair,val in mod._analytic_mixed_template().items()]}
    cert=mod.verify_full_seed_weinberg_projection('B_m1',{},source)
    assert cert['independent_components_verified']==30
    assert cert['all_component_residuals_zero'] is True
    source['components'][0]['beta_16pi2_expression']='3*k'
    with pytest.raises(ValueError,match='Stage17 tensor component'):
        mod.verify_full_seed_weinberg_projection('B_m1',{},source)


def test_full_symbolic_rejects_non_weinberg_component(monkeypatch):
    from itertools import combinations_with_replacement
    import types
    k=sp.Symbol('k')
    template={(1,2,*pair):sp.nsimplify(value)
              for pair,value in mod._analytic_mixed_template().items()}
    monkeypatch.setattr(mod,'_production_template',lambda:template)
    monkeypatch.setattr(mod,'real_quartic_entries',lambda:{})
    # Deliberately introduce one extra nu-nu HH tensor entry unsupported by the
    # certified template; checking only four (nu,e) components would miss it.
    def beta(quartic,w,i,j,a,b):
        return sp.simplify(k*template.get((i,j,a,b),0) +
                           (k if (i,j,a,b)==(1,1,3,3) else 0))
    monkeypatch.setattr(mod,'extra_scalar_pair_component',beta)
    fake=types.ModuleType('RGE.running.intermediate.ScalarOnlyTensorAdapters')
    fake.build_eft1_wilson_tensor=lambda *a,**kw:object()
    monkeypatch.setitem(sys.modules,'RGE.running.intermediate.ScalarOnlyTensorAdapters',fake)
    source={'status':'ExtraEFT1ScalarPairMixingDerived','model':'B_m1',
            'components':[{'higgs_real_indices':list(pair),
                           'beta_16pi2_expression':sp.sstr(k*val)}
                          for pair,val in mod._analytic_mixed_template().items()]}
    with pytest.raises(ValueError,match='NOT pure Weinberg'):
        mod.verify_full_seed_weinberg_projection('B_m1',{},source)


def test_local_stage18_if_available():
    import json
    p=ROOT/'output/general_t3_stage18/B_m1_extra_eft1_wilson_source_integral.json'
    if not p.exists():pytest.skip('Actual Stage18 reports only in user local repository')
    j=json.loads(p.read_text(encoding='utf-8-sig'))
    assert len(j['higgs_real_components'])==4
    assert j['status']==mod.STAGE18_STATUS
