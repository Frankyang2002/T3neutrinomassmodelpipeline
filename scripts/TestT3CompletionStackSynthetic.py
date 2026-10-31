#!/usr/bin/env python3
"""End-to-end 16-model synthetic integration test in a temporary directory.

Uses shape-correct, U(1)-covariant analytic toy beta functions. This exercises
all state, metadata, zero-control and ODE orchestration without the user's
Wolfram JSON artifacts. It is NOT a test of real RGBeta beta expressions.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.experimental.T3GeneralYukawaCoupledUV import (
    SPECS,EXTRAS,REAL,MATRIX,BETA_CONVENTION,run_model,require,
)
from Numerical.experimental.T3GeneralYukawaPhaseAudit import test_rephase_covariance


def fake_beta(expressions,environment):
    # Scalar/matrix beta with the same U(1) weight as its corresponding field.
    # Rephasing covariance, dimension and integration are real checked invariants.
    return {k:0.03*environment[k] for k in expressions}


def fixture(model:str, report_dir:Path, config_dir:Path):
    d1,d2,df,alpha=SPECS[model]
    report_dir.mkdir(parents=True,exist_ok=True)
    config_dir.mkdir(parents=True,exist_ok=True)
    keys=['gY','g2','g3','yu','yd','ye','lambdaH','y1','y2','MF',
          'mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12',
          'lambdaT3','lambdaT3GeneralQsynthetic']
    if model in EXTRAS:keys.append(EXTRAS[model])
    formulas={k:f'DUMMY[{k}]' for k in keys}
    report={
        'status':'Success',
        'metadata':{
            'dS1':d1,'dS2':d2,'dF':df,'alpha':alpha,'SharedScalar':False,
            'MajoranaFermion':alpha==-1,'GeneralScalarStage':'UV',
            'GeneralScalarQuarticBetaBasis':'RGBetaFullQuarticProjection',
            'ReportBetaConvention':BETA_CONVENTION},
        'report_betas':formulas,'betas':formulas,
    }
    (report_dir/(model+'.json')).write_text(json.dumps(report))
    diag={'real':np.eye(3,dtype=float).tolist(),
          'imag':np.zeros((3,3)).tolist()}
    heavy={'real':(1e5*np.eye(3,dtype=float)).tolist(),
           'imag':np.zeros((3,3)).tolist()}
    state={
        'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':alpha,'shared_scalar':False},
        'base_state':{
            'mu_gev':1e7,
            'sm':{'gY':.3,'g2':.6,'g3':.8,'lambdaH':.2,
                  'yu':diag,'yd':diag,'ye':diag},
            'ordinary':{'y1':diag,'y2':diag,'MF':heavy,'mS1Sq':1e6,
                'mS2Sq':1e6,'lambdaS1':.1,'lambdaS2':.12,
                'lambdaH1':.01,'lambdaH2':.01,'lambda12':.04,
                'lambdaT3':{'real':.1,'imag':0.}},
        },
        'scales':{'mu_fermion_threshold_gev':1e5},
    }
    short=f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_'+(
        f'm{-alpha}' if alpha<0 else f'p{alpha}')
    (config_dir/(short+'.json')).write_text(json.dumps(state))
    return report


def main()->int:
    with tempfile.TemporaryDirectory(prefix='t3_completion_synthetic_') as tmp:
        work=Path(tmp)
        rd=work/'native'
        cd=work/'reference'
        out=work/'trajectories'
        passed=[];phase_passed=[]
        for name in SPECS:
            report=fixture(name,rd,cd)
            output=out/(name+'.json')
            run=run_model(name,n_save=5,report_dir=rd,config_dir=cd,
                          output=output,evaluator=fake_beta)
            assert output.is_file()
            assert run['status']=='CoupledGeneralYukawaRGBetaOnlyUVPassed'
            assert run['extra_yukawa']==EXTRAS.get(name)
            assert run['endpoint_above_evolved_heavy_mass']
            if name in EXTRAS:
                assert run['seeded_vs_zero_control'] is not None
                assert run['seeded_vs_zero_control']['extra_final_zero_seed_frobenius'] < 1e-10
            else:assert run['seeded_vs_zero_control'] is None
            if SPECS[name][3]!=-1:
                p=test_rephase_covariance(report['report_betas'],
                     {k:(np.eye(3,dtype=complex) if k in MATRIX else 1.)
                      for k in report['report_betas']},evaluator=fake_beta)
                assert p['status']=='RGBetaVectorlikeFieldRephaseCovariant'
                phase_passed.append(name)
            passed.append(name)
        require(len(passed)==16 and len([x for x in passed if x in EXTRAS])==6,
                'Sixteen-model synthetic integration incomplete')
        require(len(phase_passed)==14,'Expected 14 charged Dirac models, excluding B/C(-1)')
        print('PASS: 16/16 synthetic full coupled UV models, 6/6 seed-vs-zero controls, 14/14 charged-Dirac phase tests')
        print('SYNTHETIC ONLY: live Wolfram RGBeta report integration requires a separate local run')
    return 0

if __name__=='__main__':raise SystemExit(main())
