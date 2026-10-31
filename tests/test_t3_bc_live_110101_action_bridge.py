"""Regression tests for exact exported-basis B/C action bridge.

Synthetic live JSON is only a schema fixture.  Live Wolfram/Matchete evidence
must be generated independently by the user's installed packages.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import sympy as sp
import numpy as np

from scripts.AuditT3BC110101MassBasisMap import (
    MODELS, NATIVE, RGB_NATIVE, RGB_LIVE,
    exact_map, original_matchete_records, exact_number, check_live, load,
)
from scripts.AuditT3BCLive110101ActionBridge import (
    EXPECTED_A, exact_action_identity, run, run_one,
)
from scripts.VerifyT3BCNativeLam110101Coefficients import (
    FORMAT as LAM_FORMAT, STATUS as LAM_STATUS, SYMBOLS as LAM_SYMBOLS,
)
from scripts.VerifyT3BCMatchete110101Action import (
    FORMAT as MAT_FORMAT, SOURCE as MAT_SOURCE, EXPECTED as MAT_EXPECTED,
)

ROOT = Path(__file__).resolve().parents[1]


def synthetic_matchete() -> dict:
    return {
        'format': MAT_FORMAT,
        'status': 'NativeMatcheteInteractionConstructionDiagnostic',
        'source': MAT_SOURCE,
        'matchete_scalar_off_diagonal_mass_included': False,
        'cross_tool_action_equality_verified': False,
        'physical_coupling_mapping_certified': False,
        'physical_matching_certified': False,
        'models': {model: {
            'model':model,'status':'LiveMatchete110101ActionTermsValidated',
            'matchete_action_validated':True,'native_model_built':True,
            'scalar_extra_quartics_count':7,
            'selected_110101_couplings_count':2,
            'combined_check_lagrangian_passed':True,
            'combined_check_lagrangian_return_inputform':'True',
            'source_operator_uses_plus_hc':'verified separately from GeneralScalarOperators.wl',
            'cross_tool_action_equality_verified':False,
            'physical_coupling_mapping_certified':False,
            'physical_matching_certified':False,
            'selected_terms':[{
                'name':term,'coupling':coupling,
                'coupling_present_in_native_expression':True,
                'check_lagrangian_passed':True,
                'check_lagrangian_return_inputform':'True',
                'native_expression_inputform':f'PlusHc[{coupling}[] tensor]',
            } for term,coupling in MAT_EXPECTED],
        } for model in MODELS},
    }


def synthetic_lam() -> dict:
    return {
        'format': LAM_FORMAT,
        'status':'LiveLam110101SourceCoefficientComponentAudit',
        'source':'installed RGBeta`PackageScope`Lam',
        'action_level_matchete_independently_evaluated':False,
        'physical_coupling_mapping_certified':False,
        'physical_matching_certified':False,
        'models': {model: {
            'model':model, 'status':LAM_STATUS,
            'live_native_lam_called':True, 'symbol':'RGBeta`PackageScope`Lam',
            'source_symbols':list(LAM_SYMBOLS),
            'source_scalar_field_positions':[[1,2,4,6],[2,1,3,5]],
            'registered_complex_quartic_normalization':'6*Matrix[...]',
            'unique_arrangements_per_coupling':24,
            'component_tests':64,'expected_component_tests':64,
            'expected_lam_over_registered_invariant':'1/4',
            'coefficient_extraction':'symbolic Coefficient on evaluated native Lam',
            'native_evaluated_coupling_forms':[{ 
                'orientation':orientation,'coupling':coupling,
                'evaluated_coupling_form':(
                    'Bar['+coupling+']' if orientation=='conjugate' else coupling),
            } for orientation in ('direct','conjugate') for coupling in LAM_SYMBOLS],
            'failed_component_count':0,'first_mismatches':[],
            'maximum_absolute_component_residual':0.0,
            'matchete_runtime_evaluated':False,
            'physical_coupling_mapping_certified':False,
            'physical_matching_certified':False,
        } for model in MODELS},
    }


def native_input():
    return load(NATIVE),load(RGB_NATIVE),load(RGB_LIVE)


def tensor_matrices(model):
    native,inv,cov=native_input()
    mats=original_matchete_records(native,model)
    M=sp.Matrix.hstack(*(sp.Matrix([exact_number(s) for s in row]) for row in mats))
    rgb=check_live(cov,inv,model)
    R=sp.Matrix.hstack(*(sp.Matrix([sp.Integer(int(round(float(np.real(z))))) for z in t.flat]) for t in rgb))
    return M,R,exact_map(native,model,rgb)


@pytest.mark.parametrize('model',MODELS)
def test_exact_direct_and_hc_all_components(model):
    M,R,A=tensor_matrices(model)
    report=exact_action_identity(M,R,A)
    assert A==EXPECTED_A
    assert report['exact_nonzero_component_residuals']==0
    assert report['direct_component_coefficients_compared']==32
    assert report['hermitian_conjugate_component_coefficients_compared']==32
    assert report['spurious_extra_quarter_rejected']


@pytest.mark.parametrize('model',MODELS)
def test_wrong_sign_and_quarter_both_fail(model):
    M,R,A=tensor_matrices(model)
    with pytest.raises(ValueError,match='basis map'):
        exact_action_identity(M,R,-A)
    with pytest.raises(ValueError,match='exact signed action identity'):
        exact_action_identity(M,R/4,A)


def test_source_matchete_and_rgbeta_live_schema_gates(tmp_path):
    n,inv,cov=native_input()
    paths=[tmp_path/x for x in ('mat.json','lam.json','native.json','inv.json','cov.json')]
    import json
    for p,v in zip(paths,(synthetic_matchete(),synthetic_lam(),n,inv,cov)):
        p.write_text(json.dumps(v),encoding='utf-8')
    source=ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl'
    result=run(*paths,source)
    assert result['total_direct_and_hc_coefficients_checked']==128
    assert result['live_matchete_interactions_checked']==4
    assert result['direct_cross_tool_action_subtraction_in_live_matchete_verified'] is False
    assert result['physical_coupling_mapping_certified'] is False


@pytest.mark.parametrize('mutation',(
    lambda x:x['models']['B_-1'].update(combined_check_lagrangian_passed=False),
    lambda x:x['models']['C_-1']['selected_terms'][1].update(check_lagrangian_passed=False),
    lambda x:x.update(cross_tool_action_equality_verified=True),
))
def test_matchete_live_fail_closed(tmp_path,mutation):
    import json
    n,inv,cov=native_input()
    mat=synthetic_matchete();mutation(mat)
    paths=[tmp_path/x for x in ('mat.json','lam.json','native.json','inv.json','cov.json')]
    for p,v in zip(paths,(mat,synthetic_lam(),n,inv,cov)):
        p.write_text(json.dumps(v),encoding='utf-8')
    with pytest.raises(ValueError):run(*paths,ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl')


def test_rgbeta_live_fail_closed(tmp_path):
    import json
    n,inv,cov=native_input()
    lam=synthetic_lam();lam['models']['B_-1']['failed_component_count']=1
    paths=[tmp_path/x for x in ('mat.json','lam.json','native.json','inv.json','cov.json')]
    for p,v in zip(paths,(synthetic_matchete(),lam,n,inv,cov)):
        p.write_text(json.dumps(v),encoding='utf-8')
    with pytest.raises(ValueError):run(*paths,ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl')


def test_matchete_native_tensor_mutation_rejected():
    native,inv,cov=native_input()
    fake=deepcopy(native)
    tensor=next(x for x in fake['models']['B_-1']['tensors']
                if x['symbol']=='lambdaT3ExtraQ110101I1')
    index=next(i for i,x in enumerate(tensor['components']) if x=='1')
    tensor['components'][index]='2'
    with pytest.raises(ValueError):run_one('B_-1',fake,inv,cov)


def test_no_matching_claim():
    native,inv,cov=native_input()
    result=run_one('C_-1',native,inv,cov)
    assert result['physical_coupling_mapping_certified'] is False
    assert result['direct_live_cross_tool_action_subtraction_evaluated'] is False
    assert result['mass_eigenstate_matchete_matching_certified'] is False
