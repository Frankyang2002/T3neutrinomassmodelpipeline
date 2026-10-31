"""Independent synthetic controls for opt-in B/C Higgs mass runner.

Synthetic inputs do not certify installed RGBeta; run the live Wolfram
export and VerifyT3BCFullScalarMassClosure to do that.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import VerifyT3BCFullScalarMassClosure as audit


def sample(model='B_-1',stage='EFT1'):
    return {
        'format':audit.FORMAT, 'status':audit.STATUS,
        'model':model,'stage':stage,
        'registered_mass_symbols':list(audit.NAMES),
        'physical_mapping_certified':False,
        'mass_eigenstate_threshold_matching_performed':False,
        'metadata':{
            'SharedScalar':False,'dS1':2,'dS2':2,
            'dF':1 if model=='B_-1' else 3,'alpha':-1,
            'GeneralScalarModel':model,'GeneralScalarStage':stage,
            'HiggsScalarMassSymbol':'mHSq','RegisteredQuarticCount':7,
            'MassBetaParameters':list(audit.NAMES),
            'ReportBetaConvention':audit.CONVENTION,
            'PhysicalMatcheteRGBetaMassBasisCertified':False,
            'DimensionlessCouplingsRunInThisReport':False,
            'FullRenormalizableMassSectorClosed':'PendingIndependentVerification',
        },
        'report_betas':{
            'mHSq':'mHSq + mS1Sq + mS2Sq',
            'mS1Sq':'mS1Sq + mHSq',
            'mS2Sq':'mS2Sq + mHSq',
            'muT3GeneralS1S2':'muT3GeneralS1S2',
        }
    }


def toy_evaluator(expr,env):
    g2=env['g2']; gy=env['gY']; yh=env['yu']; yd=env['yd']; ye=env['ye']
    hcoef=(6*env['lambdaH']-4.5*g2*g2-1.5*gy*gy
           +6*np.real(np.trace(yh@yh.conj().T))
           +6*np.real(np.trace(yd@yd.conj().T))
           +2*np.real(np.trace(ye@ye.conj().T)))
    mu=env[audit.NAMES[-1]]
    q=env['lambdaT3GeneralQ110101I1']
    h=hcoef*env['mHSq']+4*env['lambdaH1']*env['mS1Sq']+4*env['lambdaH2']*env['mS2Sq']
    h+=2*(q*mu.conjugate()).real
    if expr=='h':return h
    if expr=='s1':
        return (4*env['lambdaH1']*env['mHSq']
                +(-4.5*g2*g2-1.5*gy*gy+6*env['lambdaS1'])*env['mS1Sq']
                +4*env['lambda12']*env['mS2Sq']
                +12*(mu*env['lambdaS1S1bar2S2bar'].conjugate()).real)
    if expr=='s2':
        return (4*env['lambdaH2']*env['mHSq']
                +(-4.5*g2*g2-1.5*gy*gy+6*env['lambdaS2'])*env['mS2Sq']
                +4*env['lambda12']*env['mS1Sq']
                +12*(mu*env['lambdaS1barS2S2bar2'].conjugate()).real)
    if expr=='mu':
        return (6*env['lambdaS1S1bar2S2bar']*env['mS1Sq']
                +6*env['lambdaS1barS2S2bar2']*env['mS2Sq']
                +2*q*env['mHSq']
                +(-4.5*g2*g2-1.5*gy*gy+2*env['lambda12']
                  -1.5*env['lambda12Adj'])*mu
                +12*env['lambdaS1bar2S2bar2']*mu.conjugate())
    raise AssertionError('bad synthetic expression')


def synthetic_env(model='B_-1',stage='EFT1'):
    e=audit.benchmark_environment(model,stage)
    return e


def synthetic_betas():
    return dict(zip(audit.NAMES,['h','s1','s2','mu']))


def test_native_runner_is_opt_in_and_reports_all_mass_betas():
    src=(Path(__file__).resolve().parents[1]/'RGE/running/rgbeta/RunT3BCFullScalarMassRGBeta.wl').read_text()
    assert 'Length[$ScriptCommandLine] =!= 4' in src
    assert 'AddScalarMass[mHSq,{Bar@H,H}' in src
    assert 'AddScalarMass[muT3GeneralS1S2,{Bar@S1,Bar@S2}' in src
    assert 'T3RGBetaBuildEFT1[2,2,dF,-1]' in src
    assert 'T3RGBetaBuild[2,2,dF,-1]' in src
    assert 'Get[extension]' in src
    assert 'BetaTerm[sym,1]' in src
    assert 'T3RGBetaModel.wl' in src


@pytest.mark.parametrize('model',audit.MODELS)
@pytest.mark.parametrize('stage',audit.STAGES)
def test_four_live_report_contracts(model,stage):
    data=sample(model,stage)
    assert set(audit.check_payload(data,model,stage))==set(audit.NAMES)


def test_higgs_mass_absence_rejected():
    d=sample(); del d['report_betas']['mHSq']
    with pytest.raises(ValueError,match='incomplete mass beta'):
        audit.check_payload(d,'B_-1','EFT1')


def test_wrong_stage_rejected():
    d=sample(stage='UV')
    with pytest.raises(ValueError,match='model or stage'):
        audit.check_payload(d,'B_-1','EFT1')


def test_false_physical_certificate_rejected():
    d=sample(); d['physical_mapping_certified']=True
    with pytest.raises(ValueError,match='physical certification'):
        audit.check_payload(d,'B_-1','EFT1')


def test_portal_missing_rejected():
    d=sample(); d['report_betas']['mS1Sq']='mS1Sq+lambdaH1'
    with pytest.raises(ValueError,match='back-reaction'):
        audit.check_payload(d,'B_-1','EFT1')


@pytest.mark.parametrize('model',audit.MODELS)
@pytest.mark.parametrize('stage',audit.STAGES)
def test_five_real_mass_closure(model,stage):
    env=synthetic_env(model,stage)
    result=audit.affine_mass_closure(synthetic_betas(),env,toy_evaluator)
    assert len(result['affine_beta_matrix_5_by_5'])==5
    assert max(abs(x) for x in result['beta_mass_vector_gev2_per_16pi2'])>0
    assert result['affine_residual_max']<1e-5
    assert result['mass_only_frozen_dimensionless_step']['dimensionless_running_included'] is False


def test_incorrect_higgs_portal_rejected():
    def broken(expr,env):
        z=toy_evaluator(expr,env)
        if expr=='h':z+= env['mS1Sq'] * env['lambdaH1']
        return z
    with pytest.raises(ValueError,match='Higgs-doublet portal coefficient mismatch'):
        audit.affine_mass_closure(synthetic_betas(),synthetic_env(),broken)


def test_missing_standard_model_higgs_mass_beta_rejected():
    def broken(expr,env):
        z=toy_evaluator(expr,env)
        if expr=='h':z+=0.5*env['mHSq']
        return z
    with pytest.raises(ValueError,match=r'SM\+two-doublet Higgs beta fails'):
        audit.affine_mass_closure(synthetic_betas(),synthetic_env(),broken)


def test_nonaffine_mass_flow_rejected():
    def broken(expr,env):
        z=toy_evaluator(expr,env)
        if expr=='h':z+=env['mHSq']**2/1e4
        return z
    with pytest.raises(ValueError,match='non-affine'):
        audit.affine_mass_closure(synthetic_betas(),synthetic_env(),broken)


def test_nonreal_diagonal_mass_beta_rejected():
    def broken(expr,env):
        z=toy_evaluator(expr,env)
        return z+1j*100 if expr=='h' else z
    with pytest.raises(ValueError,match='nonreal'):
        audit.affine_mass_closure(synthetic_betas(),synthetic_env(),broken)


def test_missing_output_rejected(tmp_path):
    with pytest.raises(ValueError,match='missing live beta export'):
        audit.audit(tmp_path,toy_evaluator)


def test_end_to_end_synthetic_reports(tmp_path):
    for model in audit.MODELS:
        for stage in audit.STAGES:
            p=tmp_path/f"{model.replace('_-','_m')}_{stage.lower()}_full_mass_v1.json"
            d=sample(model,stage)
            d['report_betas']=synthetic_betas()
            # The live contract requires strings with mass symbols; use a
            # separate accepted key phrase, and have synthetic evaluator map it.
            d['report_betas']={
              'mHSq':'h mHSq mS1Sq mS2Sq',
              'mS1Sq':'s1 mHSq mS1Sq',
              'mS2Sq':'s2 mHSq mS2Sq',
              'muT3GeneralS1S2':'mu',
            }
            p.write_text(json.dumps(d),encoding='utf8')
    def eval_with_prefix(expr,env):
        return toy_evaluator(expr.split(' ')[0],env)
    result=audit.audit(tmp_path,eval_with_prefix)
    assert set(result['models'])==set(audit.MODELS)
    assert result['five_real_mass_components'] is True
    assert result['physical_matchete_rgbeta_phase_certified'] is False
