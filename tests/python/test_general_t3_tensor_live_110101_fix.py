"""Post-fix live Wolfram audit must pass for both SU(2) 110101 invariants."""
from pathlib import Path
import copy
import json
import pytest
from scripts.VerifyGeneralT3SU2Tensor110101Correction import check

INPUT=Path('output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json')


def load():
    if not INPUT.is_file():
        pytest.skip('Live Wolfram output not yet generated')
    return json.loads(INPUT.read_text(encoding='utf-8-sig'))


def test_known_old_native_evidence_is_rejected():
    data=load()
    old=copy.deepcopy(data)
    # A pre-fix output MUST be rejected. The current live file may itself be old.
    old['models']['B_-1']['native_expression_inputform'][1]='eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
    with pytest.raises(ValueError,match='corrected native expressions'):
        check(old)


def test_current_live_native_is_postfix_when_present():
    data=load()
    if 'eps[SU2L[fund], t3i1, t3i3]' in data['models']['B_-1']['native_expression_inputform'][1]:
        pytest.skip('Pre-fix live evidence is historical; rerun Wolfram exporter after installing fix')
    assert check(data)=={'models_verified':2,'registered_symbols_verified':4,
                         'all_native_generator_residuals_zero':True,
                         'physical_coupling_mapping_certified':False}


def test_synthetic_postfix_live_wolfram_component_and_generator_crosscheck():
    from scripts.RGBeta110101Covariance import PROPOSED,components,INDEX_GRID
    data=copy.deepcopy(load())
    for model in ('B_-1','C_-1'):
        r=data['models'][model]
        r['native_expression_inputform'][1]=PROPOSED
        r['native_covariance_audits'][1]['su2_covariant']=True
        r['native_covariance_audits'][1]['nonzero_residual_components_by_generator']=[0,0,0]
        r['native_covariance_audits'][1]['component_tensor_flat_inputform']=[str(components(PROPOSED)[idx]) for idx in INDEX_GRID]
    assert check(data)['registered_symbols_verified']==4
    data['models']['B_-1']['native_covariance_audits'][1]['nonzero_residual_components_by_generator']=[1,0,0]
    with pytest.raises(ValueError,match=r'not SU\(2\)-covariant'):
        check(data)
