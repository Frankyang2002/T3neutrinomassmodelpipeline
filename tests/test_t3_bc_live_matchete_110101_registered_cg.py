"""The fixtures below simulate JSON; they do NOT claim live Matchete execution."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from scripts.AuditT3BC110101MassBasisMap import (
    MODELS, MAT, NATIVE, RGB_NATIVE, RGB_LIVE,
    load, original_matchete_records,
)
from scripts.VerifyT3BCMatchete110101RegisteredCG import (
    CG_NAMES, FORMAT, SOURCE, STATUS,
    validate_probe, compare_one, verify_bridge,
)


def synthetic_probe() -> dict:
    source = load(NATIVE)
    rows = {}
    for model in MODELS:
        originals = original_matchete_records(source, model)
        records = []
        for i in range(2):
            records.append({
                'symbol': CG_NAMES[i],
                'matchete_coupling': MAT[i],
                'term_name': f'T3ExtraQ_1_1_0_1_0_1_{i+1}',
                'cg_present_in_native_plus_hc_interaction': True,
                'tensor_accessor': 'Matchete`GetCGTensor',
                'returned_head': 'List',
                'dimensions': [2, 2, 2, 2],
                'components': originals[i],
                'component_count': 16,
            })
        rows[model] = {
            'model': model,
            'status': 'LiveRegisteredMatcheteCGComponentsCaptured',
            'extra_quartic_inventory': 7,
            'registered_tensors': records,
            'registered_cg_components_extracted': 32,
            'live_registered_cg_components_verified': False,
            'direct_live_cross_tool_action_identity_verified': False,
            'physical_matching_certified': False,
        }
    return {
        'format': FORMAT, 'status': STATUS, 'source': SOURCE,
        'models': rows,
        'source_basis_tensors_substituted': False,
        'direct_live_cross_tool_action_identity_verified': False,
        'physical_coupling_mapping_certified': False,
        'mass_basis_matching_certified': False,
        'sequential_eft_matching_certified': False,
        'production_modified': False,
    }


@pytest.mark.parametrize('model', MODELS)
def test_registered_components_equal_stage22_exact_and_action(model):
    probe = validate_probe(synthetic_probe())
    result = compare_one(model, probe, load(NATIVE), load(RGB_NATIVE), load(RGB_LIVE))
    assert result['live_registered_components_equal_native_stage22_exact'] == 32
    assert result['live_matchete_registered_cg_direct_and_hc_coefficients_checked'] == 64
    assert result['exact_nonzero_component_residuals'] == 0
    assert result['signed_cross_tool_component_identity_verified'] is True
    assert result['physical_coupling_mapping_certified'] is False


@pytest.mark.parametrize('model', MODELS)
def test_modified_live_component_rejected(model):
    d = synthetic_probe()
    d['models'][model]['registered_tensors'][0]['components'][0] = '1001'
    probe = validate_probe(d)
    with pytest.raises(ValueError, match='differs from Stage22'):
        compare_one(model, probe, load(NATIVE), load(RGB_NATIVE), load(RGB_LIVE))


def test_swapping_registered_contractions_rejected():
    d = synthetic_probe()
    a,b = d['models']['B_-1']['registered_tensors']
    a['components'],b['components'] = b['components'],a['components']
    probe = validate_probe(d)
    with pytest.raises(ValueError, match='differs from Stage22'):
        compare_one('B_-1', probe, load(NATIVE), load(RGB_NATIVE), load(RGB_LIVE))


@pytest.mark.parametrize('change', [
    ('source_basis_tensors_substituted', True),
    ('physical_coupling_mapping_certified', True),
    ('mass_basis_matching_certified', True),
    ('sequential_eft_matching_certified', True),
    ('production_modified', True),
])
def test_false_physical_certification_rejected(change):
    d = synthetic_probe()
    d[change[0]] = change[1]
    with pytest.raises(ValueError, match='unwarranted'):
        validate_probe(d)


@pytest.mark.parametrize('tamper', ['missing_tensor', 'wrong_dims', 'missing_from_action',
                                    'missing_component', 'wrong_cg_name', 'inexact_symbol'])
def test_malformed_registration_rejected(tamper):
    d = synthetic_probe()
    a = d['models']['C_-1']['registered_tensors'][0]
    if tamper == 'missing_tensor':
        d['models']['C_-1']['registered_tensors'].pop()
    if tamper == 'wrong_dims':
        a['dimensions'] = [2, 2, 2]
    if tamper == 'missing_from_action':
        a['cg_present_in_native_plus_hc_interaction'] = False
    if tamper == 'missing_component':
        a['components'].pop()
    if tamper == 'wrong_cg_name':
        a['symbol'] = 'unrelatedCG'
    if tamper == 'inexact_symbol':
        a['components'][0] = 'unknownParam'
    with pytest.raises(ValueError):
        validate_probe(d)


def test_blocked_model_fails_without_false_success():
    d = synthetic_probe()
    d['models']['C_-1'] = {'model':'C_-1','status':'LiveRegisteredMatcheteCGBlocked',
                           'error':'GetCGTensor returned unexpected dimensions'}
    with pytest.raises(ValueError, match='blocked or incomplete'):
        validate_probe(d)


def test_bridge_requires_matching_source_hashes(tmp_path):
    source = tmp_path/'source.json'
    source.write_text('one', encoding='utf-8')
    from scripts.VerifyT3BCMatchete110101RegisteredCG import hash_bytes
    source_hash = hash_bytes(source)
    doc = {
        'format':'t3-bc-live-evidence-gated-110101-action-identity-v1',
        'status':'ExactExportedComponentActionIdentityWithSeparateLiveGates',
        'total_direct_and_hc_coefficients_checked':128,
        'live_matchete_interactions_checked':4,
        'live_rgbeta_lam_coefficients_checked':128,
        'direct_live_matchete_cg_component_equality_verified':False,
        'physical_coupling_mapping_certified':False,
        'source_sha256':{str(source.resolve()):source_hash},
        'models':[{'model':model,
                   'exact_action_identity_in_exported_su2_component_basis':True,
                   'matchete_live_cg_numerical_components_independently_extracted':False}
                  for model in MODELS],
    }
    verify_bridge(doc, {'source':source})
    source.write_text('two', encoding='utf-8')
    with pytest.raises(ValueError, match='stale'):
        verify_bridge(doc, {'source':source})


def test_wolfram_probe_uses_registered_api_and_frontend():
    script = Path(__file__).resolve().parents[1]/'scripts/ProbeT3BCMatchete110101RegisteredCG.wl'
    raw = script.read_text(encoding='utf-8')
    assert 'GetCGTensor[cg]' in raw
    assert 'UsingFrontEnd[T3BCGRun[]]' in raw
    assert 'T3ExtendedBuildExtraQuarticCandidates[model]' in raw
    assert 'CG[cg,_List]' in raw
    assert 'Exit[7]' in raw
    assert 'T3ExtendedBosonicTensorBasis[model' not in raw
