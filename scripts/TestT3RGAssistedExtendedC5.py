#!/usr/bin/env python3
"""Synthetic regression tests for the declared UV->Matchete hypothesis bridge.

Tests are NOT evidence for cross-package physical action normalization.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3RGAssistedExtendedC5 as B
from Numerical.experimental import T3ExtendedC5ExpressionReduction as C


def fixture_scenario():
    return B.read_scenario({
        'format':B.SCHEMA,'basis_hypothesis':'complex-linear-unverified-source-components',
        'physical_coupling_dictionary_certified':False,
        'hierarchical_threshold_matching_certified':False,
        'matchete_from_rgbeta':{name:{'real':-1.,'imag':0.}
                                for name in ('y1','y2','extra_yukawa','lambdaT3','extra_quartic')},
        'extra_quartic_relative_phase_radians':[0.,math.pi],
        'endpoint_buffer':1.5,'seed_factor':0.3,'n_save':5,
    })


def synthetic_uv(model:str, *,extra_zero:bool=False):
    extra=B.UV.EXTRAS[model];rg=B.EXTRA_QUARTIC[model][1]
    y1=np.array([[.004,.002,-.001],[.001,.002,.003],[.002,-.001,.006]],complex)
    y2=np.array([[.003,-.001,.001],[.001,.004,.002],[.002,.002,.003]],complex)
    u={
        'MF':np.array([[1.3e5,2e3,0],[0,1.65e5,3e3],[0,0,2.1e5]],complex),
        'mS1Sq':float(1.3e4**2),'mS2Sq':float(1.5e4**2),
        'y1':y1,'y2':y2,'lambdaT3':complex(.001,.0002),
        rg:complex(.0005,-.0001),
        'ye':np.diag([2.9e-6,.0006,.01]).astype(complex),
        extra:np.zeros((3,3),complex) if extra_zero else .12*y1*(1+.2j),
    }
    return u


def self_test()->None:
    B.self_test()
    scenario=fixture_scenario()
    for model in B.EXTRA_QUARTIC:
        a=synthetic_uv(model)
        b=synthetic_uv(model,extra_zero=True)
        converted=B.to_matchete_inputs(model,a,scenario)
        assert converted['parameters'][B.EXTRA_QUARTIC[model][0]]['real'] != []
        basis=B.vectorlike_basis(a['MF'])
        assert basis.residual<1e-12
        assert np.allclose(basis.left.conj().T@a['MF']@basis.right,np.diag(basis.masses))
        # Independent formula regression only. The two evaluators have
        # different provenance contracts; a Matchete-only object is constructed
        # here expressly for algebra comparison, not for live RGBeta transport.
        old=dict(converted,basis_declaration='explicit-matchete-source-couplings-no-rgbeta-transport')
        reference=C.evaluate_conditional_c5(model,old)
        hypothetical=B.evaluate_hypothetically_transported_c5(model,converted)
        assert np.allclose(B.from_json(reference['C5_total_gev_inverse']),
                           B.from_json(hypothetical['C5_total_gev_inverse']),rtol=2e-12,atol=1e-30)
        output=B.compare_solution(model,a,b,scenario,1e6)
        assert output['n_phases']==2
        rows=output['scenarios']
        assert abs(rows[0]['extra_relative_phase_radians'])<1e-9
        for row in rows:
            masses=row['conditional_observables']['masses_ev']
            assert len(masses)==3 and all(np.isfinite(masses))
            assert row['frobenius_C5_total_gev_inverse']>0
            assert row['frobenius_C5_extra_gev_inverse']>0
        assert abs(rows[0]['frobenius_C5_total_gev_inverse'] -
                   rows[1]['frobenius_C5_total_gev_inverse'])>1e-24
        # Turning off extra coupling suppresses the direct term exactly.
        z=B.to_matchete_inputs(model,b,scenario)
        direct=B.evaluate_hypothetically_transported_c5(model,z)['C5_extra_gev_inverse']
        assert np.max(np.abs(B.from_json(direct)))==0.
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)
        for model in B.EXTRA_QUARTIC:
            _,rg=B.EXTRA_QUARTIC[model]
            alpha=B.UV.SPECS[model][3]
            src=root/'RGE/general/generated'/f'{model[0]}_{alpha:+d}_extra_quartics.wl'
            src.parent.mkdir(parents=True,exist_ok=True)
            fields='{Bar @ H, Bar @ H, S1, S2}' if alpha==0 else '{Bar @ H, Bar @ H, Bar @ S1, Bar @ S2}'
            src.write_text(f'AddQuartic[{rg},\n  {fields},\n GroupInvariant -> (1 &),\n SelfConjugate -> False];')
            payload={'report_betas':{key:'x' for key in ('y1','y2',B.UV.EXTRAS[model],
                           'MF','mS1Sq','mS2Sq','lambdaT3',rg)}}
            result=B.check_rgbeta_quartic(model,root,payload)
            assert result['rgbeta_name']==rg
            src.write_text(src.read_text().replace('Bar @ H, Bar @ H','H, Bar @ H'))
            try:B.check_rgbeta_quartic(model,root,payload)
            except ValueError:pass
            else:raise AssertionError('Wrong RGBeta monomial accepted')
    # End-to-end orchestration with the REAL SciPy solver but a zero-beta
    # fake payload. No synthetic result is confused with the native exports.
    real_integrate=B.UV.integrate
    with tempfile.TemporaryDirectory() as t:
        root=Path(t)
        model='A_p0'
        _,rgb=B.EXTRA_QUARTIC[model]
        path=root/'RGE/general/generated/A_+0_extra_quartics.wl'
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(f'AddQuartic[{rgb}, {{Bar @ H, Bar @ H, S1, S2}}, '
                        'GroupInvariant -> (1 &), SelfConjugate -> False];')
        state=synthetic_uv(model)
        payload={'report_betas':{k:'synthetic_unused' for k in state}}
        def fake_integrate(source,initial,start,end,*,n_save,**kwargs):
            return real_integrate(source,initial,start,end,n_save=n_save,
                                  evaluator=lambda eq,st: {name:0*st[name] for name in eq})
        evidence={'status':'ExactNativeOneStepOrderedCoefficientReduced',
                  'exact_reduction_of_all_three_sectors_after_identification':True,
                  'source_sha256':{'full_C5':'synthetic_fake'}}
        with patch.object(B.C5,'verify_model',return_value=evidence), \
             patch.object(B.UV,'load_case',return_value=(payload,{})), \
             patch.object(B.UV,'prepare_state',return_value=state), \
             patch.object(B.UV,'endpoint_for',return_value=(1e6,4e5)), \
             patch.object(B.UV,'integrate',side_effect=fake_integrate):
            outcome=B.run_model(model,scenario,root=root)
            assert outcome['status']==B.STATUS
            assert outcome['uv_rhs_evaluations']>0
            assert outcome['uv_control_rhs_evaluations']>0
            assert outcome['phase_scan']['n_phases']==2
            assert outcome['absolute_physical_matchete_rgbeta_dictionary_certified'] is False
            assert outcome['sequential_threshold_matching_performed'] is False
    try:B.validate_endpoint('A_p0',synthetic_uv('A_p0'),1e5)
    except ValueError:pass
    else:raise AssertionError('Premature heavy threshold crossing accepted')
    print('PASS: six-model synthetic UV-to-one-step C5, SVD, phase branches, controls and source-monomial tamper')
    print('SYNTHETIC ONLY: native coupling sign/phase dictionary and hierarchical matching remain uncertified')


if __name__=='__main__':self_test()
