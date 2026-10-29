"""Stage 13: independent tensor identity and safe coupled-RG orchestration."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.experimental.GeneralT3TensorMap import (
    symmetric_doublet_projector, epsilon_doublet,
    all_bar_tensor_from_projector, build_tensor_mapping,
)
from Numerical.experimental.GeneralT3CoupledUV import (
    MODELS, ALL_COMPLEX, _complex, _matrix,
)


def _tensor_report():
    return {"status":"MatcheteTensorConventionVerified","results":[
        {"model":m,"tensors":[
            {"coupling":c,"status":"CanonicalTensorProportional",
             "matchete_over_canonical_factor_inputform":"1",
             "all_component_residuals_zero":True}
            for c in ("lambdaT3ExtraQ020020I1","lambdaT3ExtraQ020200I1")
        ]} for m in ("B_m1","C_m1")
    ]}


def _source():
    return '''AddQuartic[
            lambdaHHdagS2S2,
            {Bar @ H, Bar @ H, S2, S2},
            GroupInvariant -> (
                delS2[SU2L, a, #1, #2] *
                delS2[SU2L, a, #3, #4] &
            ), SelfConjugate -> False];
        AddQuartic[ lambdaHHdagS1barS1bar,
            {Bar @ H, Bar @ H, Bar @ S1, Bar @ S1},
            GroupInvariant -> (delS2[SU2L, a, #1, #2] * delS2[SU2L, a, #3, #4] &),
            SelfConjugate -> False];'''


def test_projector_is_idempotent_and_rank_three():
    p = symmetric_doublet_projector().reshape(4,4)
    assert np.array_equal(p@p,p)
    assert np.trace(p) == 3


def test_epsilon_intertwiner_matches_allbar_tensor_exactly():
    eps = epsilon_doublet()
    expected = (np.einsum('ac,bd->abcd',eps,eps)+
                np.einsum('ad,bc->abcd',eps,eps))/2.0
    np.testing.assert_array_equal(all_bar_tensor_from_projector(),expected)


def test_tensor_mapping_four_factors_but_not_full_couplings(tmp_path):
    report=tmp_path/'matchete.json'; report.write_text(json.dumps(_tensor_report()))
    rgb=tmp_path/'RGBeta.wl'; rgb.write_text(_source())
    result=build_tensor_mapping(report,rgb)
    assert result['status']=='DoubletGroupTensorRatiosVerified'
    assert result['full_lagrangian_coupling_identification_verified'] is False
    assert {r['rgbeta_over_matchete_tensor_factor'] for r in result['mappings']}=={1.0}
    assert len(result['mappings'])==4


def test_tensor_mapping_missing_matchete_result_rejected(tmp_path):
    payload=_tensor_report();payload['results'][0]['tensors'].pop()
    (tmp_path/'mat.json').write_text(json.dumps(payload))
    (tmp_path/'rgb.wl').write_text(_source())
    with pytest.raises(ValueError, match='exactly the four'):
        build_tensor_mapping(tmp_path/'mat.json',tmp_path/'rgb.wl')


def test_tensor_mapping_changed_rgbeta_field_content_rejected(tmp_path):
    (tmp_path/'mat.json').write_text(json.dumps(_tensor_report()))
    (tmp_path/'rgb.wl').write_text(_source().replace('Bar @ S1, Bar @ S1','S1, S1'))
    with pytest.raises(ValueError, match='certified delS2 structure'):
        build_tensor_mapping(tmp_path/'mat.json',tmp_path/'rgb.wl')


def test_complex_and_yukawa_inputs_validate():
    assert _complex({'real':.1,'imag':-.2},'z')==complex(.1,-.2)
    with pytest.raises(ValueError,match='real/imag'):
        _complex({'real':1},'z')
    pars={'y1_real':np.eye(3).tolist(),'y1_imag':np.zeros((3,3)).tolist()}
    np.testing.assert_array_equal(_matrix(pars,'y1'),np.eye(3))
    pars['y1_imag']=[[0,0]]
    with pytest.raises(ValueError,match='3x3'):
        _matrix(pars,'y1')


def test_source_models_and_special_state_are_explicit():
    assert MODELS=={'B_m1':1,'C_m1':3}
    assert len(ALL_COMPLEX)==6
    assert {'lambdaHHdagS2S2','lambdaHHdagS1barS1bar'}.issubset(ALL_COMPLEX)


def test_stage13_scripts_exist_and_use_existing_uv_runner():
    source=(ROOT/'Numerical/experimental/GeneralT3CoupledUV.py').read_text()
    assert 'from Numerical.running.UVRunner import run_uv_segment' in source
    assert 'from Numerical.core.RGBetaEvaluator import load_rgbeta_payload' in source
    assert 'full_matchete_rgbeta_coupling_mapping_verified' in source
    assert 'endpoint > 1.001*np.max' in source
    assert (ROOT/'scripts/RunGeneralT3CoupledUVRGE.ps1').exists()
