"""Offline checks for a *saved-input* quadratic RG closure audit."""
import json
from pathlib import Path

import pytest

from scripts import AuditT3QuadraticRGClosure as a


def test_model_list_and_paths():
    assert len(a.MODELS) == len(set(a.MODELS)) == 16
    assert len(a.EXPECTED) == 8
    assert a.tag('A_+0') == 'A_p0'
    assert a.tag('C_-1') == 'C_m1'


def test_uv_trace_removed_when_mu_zero():
    expr='6*lambdaS1S1bar2S2bar*mS1Sq + muT3GeneralS1S2*(Tr[y1 . Trans[Bar[y1]]] + Tr[y2 . Trans[Bar[y2]]])'
    src,zero,sources=a.source_at_zero(expr,'muT3GeneralS1S2',{'lambdaS1S1bar2S2bar'})
    assert not zero
    assert '__TR' not in a.exprtext(src)
    assert sources==['lambdaS1S1bar2S2bar']
    assert a.evaluate(src,{'mS1Sq':1e6,'lambdaS1S1bar2S2bar':-0.01})==-60000


def test_source_from_sum_inside_parentheses():
    expr='(6*lambdaT3GeneralQ000211I1 + 2*lambdaT3GeneralQ000211I2)*mS2Sq + 12*muT3GeneralS1Sq*lambdaS1'
    src,zero,sources=a.source_at_zero(expr,'muT3GeneralS1Sq',{
        'lambdaT3GeneralQ000211I1','lambdaT3GeneralQ000211I2'})
    assert sources==['lambdaT3GeneralQ000211I1','lambdaT3GeneralQ000211I2']
    val=a.evaluate(src,{'mS2Sq':1e6,'lambdaT3GeneralQ000211I1':0.01,
                        'lambdaT3GeneralQ000211I2':-0.01})
    assert val==pytest.approx(40000)


def test_zero_all_extra_exactly():
    src,zero,symbols=a.source_at_zero('6*lam*mS1Sq+2*mu*base+Bar[mu]*other', 'mu',{'lam'})
    assert not zero
    assert 'lam' in a.exprtext(src)


def test_detect_baseline_only_source():
    with pytest.raises(ValueError,match='nonquartic source|does not close'):
        a.source_at_zero('3*lambdaT3*mS1Sq+2*mu*g2', 'mu',{'extra'})


def test_detect_non_closed_source_without_extra():
    with pytest.raises(ValueError,match='does not close'):
        a.source_at_zero('3*mS1Sq+2*mu*g2', 'mu',{'extra'})


def test_invalid_division_by_mu():
    with pytest.raises(ValueError,match='denominator'):
        a.source_at_zero('mS1Sq/mu','mu',{'lam'})


def test_unsupported_calls_rejected():
    with pytest.raises(ValueError):
        a.source_at_zero('Sum[lam,1,3]*mS1Sq','mu',{'lam'})


def test_complex_source_evaluation():
    expr='6*lam*mS1Sq + 6*mu*gY'
    src,_,_=a.source_at_zero(expr,'mu',{'lam'})
    assert a.evaluate(src,{'lam':0.01+0.02j,'mS1Sq':1e6})==60000+120000j


def test_read_json_fail_closed(tmp_path):
    with pytest.raises(ValueError,match='missing required'):
        a.read_json(tmp_path/'missing.json')


def create_fixture(tmp_path:Path):
    for model in a.MODELS:
        kind=a.EXPECTED.get(model,'none')
        stem=a.tag(model)
        sym=a.MASS_SYMBOLS[kind]
        expr='6*lambdaT3GeneralQ000211I1*mS2Sq+2*'+sym+'*lambdaS1' if sym else None
        for stage in ('UV','EFT1'):
            payload={
                'format':'general-t3-quadratic-rgbeta-experimental-v1',
                'status':'SuccessCandidateMassBeta',
                'metadata':{
                    'GeneralScalarModel':model,'GeneralScalarStage':stage,
                    'QuadraticKind':kind,'QuadraticMassSymbol':sym,
                    'PhysicalMatcheteRGBetaMassBasisCertified':False,
                    'QuadraticPrefactorConvention':'V2 = mu2 * invariant + h.c.; no 1/2',
                    'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
                    'GeneralScalarNewBetaKeys':['lambdaT3GeneralQ000211I1'] if kind!='none' else [],
                },
                'report_betas':{sym:expr} if sym else {'mS1Sq':'0'},
            }
            (tmp_path/f'{stem}_{stage.lower()}_quadratic_live_rgbeta.json').write_text(json.dumps(payload),encoding='utf-8')
        run={
            'model':model,'format':'general-t3-candidate-quadratic-eft1-v1',
            'quadratic_kind':kind,'status':'CandidateQuadraticAndFullQuarticEFT1Integrated',
            'physical_matchete_rgbeta_mapping_certified':False,
            'physical_neutrino_mass_prediction':False,
            'mu2_beta_included':kind!='none',
            'mu_initial_gev':1e5,'mu_final_gev':1e3,
            'initial_baseline':{'mS1Sq':1e6,'mS2Sq':2e6,'lambdaS1':0.1},
            'initial_new_quartics':{'lambdaT3GeneralQ000211I1':{'real':0.01,'imag':0}} if kind!='none' else {},
        }
        (tmp_path/f'{stem}_eft1_quadratic_running.json').write_text(json.dumps(run),encoding='utf-8')


def test_end_to_end_synthetic_16(tmp_path):
    create_fixture(tmp_path)
    report=a.build(tmp_path)
    assert report['model_count']==16
    assert report['quadratic_slots']==8
    assert len(report['models_with_nonzero_saved_boundary_source'])==8
    assert len(report['source_sha256'])==48
    b=report['models']['B_-1']
    assert abs(b['eft1_zero_mu2_source_at_saved_nonzero_quartic_boundary']['real']-120000)<1e-10
    assert 0<b['ratio_abs_generated_mu2_over_min_input_m2']<1e-2
    a.write(report,tmp_path/'out')
    assert (tmp_path/'out'/'quadratic_rg_zero_boundary_audit.json').is_file()
    assert (tmp_path/'out'/'quadratic_rg_zero_boundary_summary.md').is_file()
    assert len((tmp_path/'out'/'quadratic_rg_zero_boundary_audit.csv').read_text().splitlines())==17


def test_missing_stage_fails(tmp_path):
    create_fixture(tmp_path)
    (tmp_path/'B_m1_uv_quadratic_live_rgbeta.json').unlink()
    with pytest.raises(ValueError,match='missing required input'):
        a.build(tmp_path)


def test_stage_mismatch_fails(tmp_path):
    create_fixture(tmp_path)
    p=tmp_path/'B_m1_eft1_quadratic_live_rgbeta.json'
    o=json.loads(p.read_text())
    o['metadata']['GeneralScalarStage']='UV'
    p.write_text(json.dumps(o))
    with pytest.raises(ValueError,match='mismatched GeneralScalarStage'):
        a.build(tmp_path)


def test_source_proof_fails_for_extra_mass_source(tmp_path):
    create_fixture(tmp_path)
    p=tmp_path/'B_m1_eft1_quadratic_live_rgbeta.json'
    o=json.loads(p.read_text())
    o['report_betas']['muT3GeneralS1S2']='mS2Sq + muT3GeneralS1S2*lambdaS1'
    p.write_text(json.dumps(o))
    with pytest.raises(ValueError,match='does not close'):
        a.build(tmp_path)


def test_report_rejects_false_physical_claim(tmp_path):
    create_fixture(tmp_path)
    p=tmp_path/'C_m1_eft1_quadratic_running.json'
    o=json.loads(p.read_text())
    o['physical_neutrino_mass_prediction']=True
    p.write_text(json.dumps(o))
    with pytest.raises(ValueError,match='unexpected claims'):
        a.build(tmp_path)
