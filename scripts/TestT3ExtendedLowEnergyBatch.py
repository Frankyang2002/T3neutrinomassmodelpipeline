#!/usr/bin/env python3
"""Regression for conditional low-energy batch; independent fake-RGE controls.

The mocked SM evolution is explicitly a software test. On the user's full
repository, a separate REAL SM+Weinberg evolution smoke test is also run.
"""
from __future__ import annotations
import math
from pathlib import Path
import json
import sys
import tempfile
import types
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3ExtendedLowEnergyStudy as LE
from Numerical.experimental import T3ExtendedConventionLedger as LEDGER
from Numerical.experimental import T3RGAssistedExtendedC5 as BR
from scripts.TestT3RGAssistedExtendedC5 import fixture_scenario, synthetic_uv


def fake_evolution(state,mu_initial,mu_final,*,rtol,atol):
    assert mu_initial>mu_final>0
    assert rtol>0 and atol>0
    factor=math.exp(-0.01*math.log(mu_initial/mu_final))
    return types.SimpleNamespace(K=state.K*factor,ye=state.ye,
                                 solver_success=True,nfev=7)


class FakeSM:
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)
    def validated(self):
        assert np.shape(self.K)==(3,3)
        return self


def model_state(model,*,extra_zero=False):
    x=synthetic_uv(model,extra_zero=extra_zero)
    x.update({'gY':0.36,'g2':0.63,'g3':1.01,'lambdaH':0.13,
              'yu':np.diag([.00001,.003,.8]).astype(complex),
              'yd':np.diag([.00002,.0004,.018]).astype(complex)})
    return x


def real_sm_weinberg_smoke()->None:
    """The real production solver is tested on Windows, never replaced by fake data."""
    try:
        from Numerical.running.SMWeinbergEvolution import (
            SMWeinbergInitialConditions,evolve_sm_weinberg)
    except ModuleNotFoundError:
        print('SKIP: isolated environment lacks production SM+Weinberg solver; live smoke deferred')
        return
    val=SMWeinbergInitialConditions(gY=.36,g2=.63,g3=1.0,lambdaH=.13,
         ye=np.diag([3e-6,.0006,.01]).astype(complex),
         yu=np.diag([1e-5,.004,.83]).astype(complex),
         yd=np.diag([2e-5,.0004,.018]).astype(complex),
         K=np.array([[1e-15,.1e-15j,0],[.1e-15j,2e-15,0],[0,0,3e-15]],complex))
    data=evolve_sm_weinberg(val,3000.,91.1876,rtol=3e-8,atol=1e-24)
    assert data.solver_success and np.isfinite(data.K).all()
    assert np.allclose(data.K,data.K.T,rtol=1e-12,atol=1e-27)
    assert data.nfev>0
    print('PASS: real production 76-variable SM+Weinberg ODE smoke, complex symmetric C5')



def test_source_ledger() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        src={
            'Lagrangian/interactions/T3Topology.wl':
                '"Expression" -> PlusHc[ lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels] ]\n',
            'Lagrangian/interactions/GeneralScalarOperators.wl':
                'contribution = coupling[] product If[cg === None, 1, CG[cg, labels]]; "Expression" -> PlusHc[contribution]\n',
            'RGE/running/rgbeta/T3RGBetaModel.wl':
                'T3RGBetaAddMixingQuartic[ x ] := AddQuartic[lambdaT3, {H, H, S1, Bar @ S2}]; T3RGBetaYukawaInvariant[x];\n',
        }
        for name,body in src.items():
            p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(body)
        for model,(mat,rg) in BR.EXTRA_QUARTIC.items():
            alpha=0 if model.endswith('p0') else -2
            word='{Bar @ H, Bar @ H, S1, S2}' if alpha==0 else '{Bar @ H, Bar @ H, Bar @ S1, Bar @ S2}'
            path=root/'RGE/general/generated'/f'{model[0]}_{alpha:+d}_extra_quartics.wl'
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(f'AddQuartic[{rg},\n {word},\n GroupInvariant -> (1 &), SelfConjugate -> False];')
        def proof(_container,model):
            mat=BR.EXTRA_QUARTIC[model][0]
            return {'extra_quartic_matchete':mat,'extra_yukawa':BR.UV.EXTRAS[model],
                    'source_sha256':{'full_C5':'synthetic_only'}}
        with patch.object(LEDGER.C5,'verify_model',side_effect=proof):
            report=LEDGER.audit(root)
            assert report['n_native_matchete_reductions']==6
            assert report['physical_coupling_dictionary_certified'] is False
            altered=root/'RGE/general/generated/E_+0_extra_quartics.wl'
            altered.write_text(altered.read_text().replace('Bar @ H, Bar @ H','H, Bar @ H'))
            try:LEDGER.audit(root)
            except ValueError:pass
            else:raise AssertionError('Corrupt source monomial accepted by convention ledger')
        print('PASS: six-model source ledger, provenance and wrong-field-word negative control')

def self_test()->None:
    LE.self_test()
    test_source_ledger()
    scenario=fixture_scenario()
    study=LE.load_study({'format':'t3-extended-c5-conditional-low-energy-v1',
                         'physical_dictionary_certified':False,
                         'hierarchical_threshold_matching_certified':False,
                         'allow_unverified_sudden_decoupling':True,
                         'endpoint_buffers':[1.15,1.5],
                         'mu_low_gev':91.1876,'atol':1e-24})
    with patch.dict(sys.modules,{'Numerical.running.SMWeinbergEvolution':types.SimpleNamespace(
                SMWeinbergInitialConditions=FakeSM)}):
        for model in LE.MODELS:
            end=model_state(model)
            zero=model_state(model,extra_zero=True)
            outcome=LE.evaluate_endpoint(model,end,end,zero,4e5,scenario,study,evolve=fake_evolution)
            assert outcome['status']==LE.STATUS
            assert len(outcome['phase_scenarios'])==2
            assert outcome['threshold_diagnostics']['fermion_first_hierarchy_in_mass_proxies']
            for row in outcome['phase_scenarios']:
                assert row['c5_linearity_residual']<1e-11
                assert len(row['low_energy_observables']['masses_ev'])==3
                assert row['direct_extra_fraction_at_low']>0
                assert row['seeded_minus_zero_seed_fraction_at_low']>0
                assert not outcome['low_energy_predictions_physical_certified']
        # Nonlinear C5 evolution would invalidate the decomposition and MUST fail.
        def nonlinear(state,*args,**kwargs):
            r=fake_evolution(state,*args,**kwargs)
            return types.SimpleNamespace(**{**vars(r),'K':r.K+1e-17*np.eye(3,dtype=complex)})
        try:LE.evaluate_endpoint('A_p0',model_state('A_p0'),model_state('A_p0'),
                                  model_state('A_p0',extra_zero=True),4e5,scenario,study,evolve=nonlinear)
        except ValueError:pass
        else:raise AssertionError('Nonlinear Weinberg evolution passed exact linearity guard')
    # End-to-end 2 matching scales, 2 phases with TRUE SciPy UV integration and
    # mocked zero-beta functions; the SM EFT mock is clearly labelled above.
    model='A_p0'
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)
        rg=BR.EXTRA_QUARTIC[model][1]
        p=root/'RGE/general/generated/A_+0_extra_quartics.wl'
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text('AddQuartic['+rg+', {Bar @ H, Bar @ H, S1, S2}, GroupInvariant -> (1 &), SelfConjugate -> False];')
        initial=model_state(model)
        payload={'report_betas':{name:'unused' for name in initial}}
        previous=BR.UV.integrate
        def integrate_zero(*args,**kw):
            return previous(*args,**kw,evaluator=lambda eq,state:{name:0*state[name] for name in eq})
        proof={'status':'ExactNativeOneStepOrderedCoefficientReduced',
               'exact_reduction_of_all_three_sectors_after_identification':True,
               'source_sha256':{'full_C5':'fake-source'}}
        with patch.object(LE.C5,'verify_model',return_value=proof),\
             patch.object(LE.UV,'load_case',return_value=(payload,{})),\
             patch.object(LE.UV,'prepare_state',return_value=initial),\
             patch.object(LE.UV,'endpoint_for',return_value=(1e6,4e5)),\
             patch.object(LE.UV,'integrate',side_effect=integrate_zero),\
             patch.dict(sys.modules,{'Numerical.running.SMWeinbergEvolution':types.SimpleNamespace(
                 SMWeinbergInitialConditions=FakeSM)}),\
             patch.object(LE,'transport',wraps=LE.transport) as spy:
            # Use explicit mock only for synthetic ODE, not real physics.
            original=LE.evaluate_endpoint
            def fake_evaluate(*args,**kwargs):
                return original(*args,**kwargs,evolve=fake_evolution)
            with patch.object(LE,'evaluate_endpoint',side_effect=fake_evaluate):
                output=LE.run_model(model,scenario,study,root=root)
            assert len(output['matching_scale_scan'])==2
            assert output['matching_scale_scan'][0]['full_RGBeta_beta_blocks']==len(initial)
            assert output['final_SM_Weinberg_RGE_performed_under_sudden_decoupling'] is True
            assert output['physical_low_energy_prediction_certified'] is False
        p.write_text(p.read_text().replace('Bar @ H, Bar @ H','H, Bar @ H'))
        try:BR.check_rgbeta_quartic(model,root,payload)
        except ValueError:pass
        else:raise AssertionError('Tampered quartic word passed')
    print('PASS: 6/6 conditional C5 low-energy mock, 2-endpoint 2-phase SciPy UV, nonlinear countercheck, source tamper')
    print('SYNTHETIC low-energy evolution in primary regression; real RGE smoke is separate')
    real_sm_weinberg_smoke()


if __name__=='__main__':self_test()
