"""Stage24 four-quartic post-fix source and live-verifier regression coverage."""
from copy import deepcopy
from itertools import product
from pathlib import Path

import pytest

from scripts.Stage24BCRemainderCovariance import PROPOSED, FIELDS, MODELS, covariance, components
from scripts.VerifyGeneralT3Stage24BCPostFix import FORMAT, verify

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'RGE/running/rgbeta/T3RGBetaModel.wl'


@pytest.mark.parametrize('symbol',sorted(PROPOSED))
def test_postfix_expression_su2_covariance(symbol):
    assert covariance(PROPOSED[symbol],FIELDS[symbol])['nonzero_residuals_T1_T2_T3']==[0,0,0]


def test_four_full_source_replacements_and_existing_110101_fix():
    assert SOURCE.is_file()
    s=SOURCE.read_text(encoding='utf-8')
    for sym in PROPOSED:
        begin = s.index('AddQuartic[\n            '+sym+',')
        end = s.index('SelfConjugate -> False',begin)
        block = s[begin:end]
        assert 'delS2[' not in block
        assert ('eps[SU2L @ fund' in block)
        assert '/ 2 &' in block
    begin=s.index('AddQuartic[\n            lambdaHHdagS1barS2barCross,')
    end=s.index('SelfConjugate -> False',begin)
    assert 'del[SU2L @ fund, #1, #3]' in s[begin:end]


def _synthetic():
    """Generate in-memory fixture; this is NOT live RGBeta evidence."""
    rows={}
    mt={}
    for model in MODELS:
        rows[model]={}
        mt[model]={'tensors':[]}
        for sym in PROPOSED:
            v=[str(components(PROPOSED[sym])[i]) for i in product(range(2),repeat=4)]
            audit={'su2_covariant':True,
                   'nonzero_residual_components_by_generator':[0,0,0],
                   'component_tensor_flat_inputform':v}
            rows[model][sym]={'native':deepcopy(audit),'proposed':deepcopy(audit)}
            from scripts.Stage22RGBetaExplicitBasis import SOURCE_ORDER
            ct=[FIELDS[sym].count(f) for f in SOURCE_ORDER]
            mt[model]['tensors'].append({'symbol':'synthetic_'+sym,'species':FIELDS[sym],
                'counts':ct,'dimensions':[2,2,2,2],'components':v})
    w={'format':FORMAT,'native_non_covariant_slots':0,'native_covariant_slots':8,
        'models':rows,'physical_coupling_mapping_certified':False}
    return w,{'models':mt}


def test_synthetic_verifier_accepts_exact_components():
    w,m=_synthetic()
    summary=verify(w,m)
    assert summary['native_su2_pass']==8
    assert summary['native_exact_matchete_tensor_matches']==8
    assert summary['physical_coupling_mapping_certified'] is False


def test_synthetic_verifier_rejects_nonzero_gauge_residual():
    w,m=_synthetic()
    w['models']['B_-1'][next(iter(PROPOSED))]['native']['nonzero_residual_components_by_generator']=[1,0,0]
    with pytest.raises(ValueError,match='noncovariant'):
        verify(w,m)


def test_synthetic_verifier_rejects_wrong_native_tensor_component():
    w,m=_synthetic()
    w['models']['C_-1'][next(iter(PROPOSED))]['native']['component_tensor_flat_inputform'][0]='100'
    with pytest.raises(ValueError,match='incorrect native'):
        verify(w,m)


def test_synthetic_verifier_rejects_wrong_matchete_tensor():
    w,m=_synthetic()
    m['models']['B_-1']['tensors'][0]['components'][0]='100'
    with pytest.raises(Exception):
        verify(w,m)


def test_synthetic_verifier_rejects_historical_export_format():
    w,m=_synthetic()
    w['format']='stage24-live-rgbeta-bc-other-quartic-covariance-v1'
    with pytest.raises(ValueError,match='source phase'):
        verify(w,m)


def test_live_integration_when_postfix_export_available():
    import json
    folder=ROOT/'output/general_t3_stage22'
    p=folder/'stage24_bc_postfix_wolfram_native_covariance.json'
    m=folder/'matchete_tensors.json'
    if not (p.is_file() and m.is_file()):
        pytest.skip('Need fresh live Stage24 post-fix Wolfram export')
    status=verify(json.loads(p.read_text(encoding='utf-8-sig')),
                  json.loads(m.read_text(encoding='utf-8-sig')))
    assert status['native_exact_matchete_tensor_matches']==8
