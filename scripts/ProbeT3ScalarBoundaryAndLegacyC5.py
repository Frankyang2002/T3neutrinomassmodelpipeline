"""Read-only scalar-boundary sensitivity and legacy C5 parity investigation.

Candidate RGBeta general-quartic calculation; NOT a refit, or a certified EFT.
Uses existing campaign RGE and matching kernels unchanged. Outputs only new files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts.DiagnoseT3GeneralScalarScreens import screen_details, legacy_parity

DEFAULT_FACTORS = (0.25, 0.5, 0.75, 1.0, 1.25)
DEFAULT_CASE = 'small_uv_running_5e-9'
LEGACY_CASE = 'large_uv_frozen_0p5'


def mass_rescale(cfg, factor):
    """Scale the two UV scalar mass *squares* by factor**2; leave all else unchanged."""
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError('Scalar mass factor must be finite and positive')
    result = deepcopy(cfg)
    ordinary = result['base_state']['ordinary']
    for name in ('mS1Sq', 'mS2Sq'):
        old = float(ordinary[name])
        if old <= 0 or not math.isfinite(old):
            raise ValueError('Nonpositive original mass square '+name)
        ordinary[name] = old * factor**2
    return result


def complex_matrix(raw):
    return np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)


def comparison_metrics(ref, fast):
    """Diagnostic hypotheses only: no automatic parity correction is applied."""
    ref = np.asarray(ref, dtype=complex)
    fast = np.asarray(fast, dtype=complex)
    if ref.shape != (3, 3) or fast.shape != (3, 3):
        raise ValueError('Expected 3x3 C5 matrices')
    denom = max(float(np.max(np.abs(ref))), 1e-30)
    delta = fast-ref
    inner = np.vdot(fast, ref)
    norm_fast_sq = float(np.vdot(fast, fast).real)
    best_alpha = complex(inner/norm_fast_sq) if norm_fast_sq>0 else None
    ratios = {}
    for name, candidate in [('fast',fast), ('minus_fast',-fast),
                            ('conjugate_fast',fast.conj()),
                            ('transpose_fast',fast.T),
                            ('conjugate_transpose_fast',fast.conj().T)]:
        ratios[name] = float(np.max(np.abs(ref-candidate))/denom)
    if best_alpha is not None:
        ratios['best_complex_scale_fast'] = float(np.max(np.abs(ref-best_alpha*fast))/denom)
    return {
        'reference_frobenius_norm': float(np.linalg.norm(ref)),
        'fast_frobenius_norm': float(np.linalg.norm(fast)),
        'relative_max_error': float(np.max(np.abs(delta))/denom),
        'max_reference_abs': float(np.max(np.abs(ref))),
        'max_fast_abs': float(np.max(np.abs(fast))),
        'max_difference_abs': float(np.max(np.abs(delta))),
        'entry_max_difference': [int(i) for i in np.unravel_index(np.argmax(np.abs(delta)),delta.shape)],
        'hypothesis_relative_errors': ratios,
        'least_squares_scale_reference_over_fast': (None if best_alpha is None else
                     {'real':float(best_alpha.real),'imag':float(best_alpha.imag)}),
        'note':'Complex scaling/sign/conjugation tests are diagnostics, not modifications to the kernel.',
    }


def load_eft_beta(model, beta_campaign):
    p = bulk.OUT/beta_campaign/'betas'/f'{bulk.safe(model)}.json'
    return bulk.read(p),p


def scalar_probe(model, case, *, beta_campaign, factors, extra_real):
    registry = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(registry,model)
    if reused:
        raise ValueError('Mixed-quadratic B/C cases require separate evaluator')
    history = bulk.historical(case,model,entry)
    payload,beta_path=load_eft_beta(model,beta_campaign)
    q0=tuple(complex(extra_real) for _ in new)
    result=[]
    for factor in factors:
        row={'model':model,'case':case,'scalar_mass_factor':factor,
             'uv_mS1Sq_gev2':None,'uv_mS2Sq_gev2':None, 'status':'NotRun'}
        try:
            cfg=mass_rescale(history['cfg'],factor)
            row['uv_mS1Sq_gev2']=cfg['base_state']['ordinary']['mS1Sq']
            row['uv_mS2Sq_gev2']=cfg['base_state']['ordinary']['mS2Sq']
            vec,layout,unpack=bulk.uv_base(cfg,entry)
            if history['uv_on']:
                betas,uv_path=bulk.uv_betas(model,entry,new,layout)
                uv,qf,nuv=bulk.run_uv(vec,layout,betas,new,q0)
                row['uv_beta_path']=str(uv_path.relative_to(bulk.ROOT))
            else:
                uv=unpack(vec,layout,mu_gev=1e5)
                qf=q0;nuv=0
            eft,masses,qs,neft=bulk.run_eft(uv,qf,payload,model,new,())
            screen=screen_details(uv,eft,masses,qf,qs)
            original=bulk.testscreen(uv,eft,masses,qf,qs)
            if screen['passed']!=original['screen']:
                raise RuntimeError('Screen disagreement with bulk campaign')
            row.update({'status':'ScreenPassed' if screen['passed'] else 'ScreenRejected',
                        'screen':screen,'uv_nfev':int(nuv),'eft_nfev':int(neft),
                        'muF_mS1Sq_gev2':float(uv.mS1Sq),
                        'muF_mS2Sq_gev2':float(uv.mS2Sq),
                        'muS_mS1Sq_gev2':float(eft.mS1Sq),
                        'muS_mS2Sq_gev2':float(eft.mS2Sq)})
        except (ValueError,RuntimeError,ArithmeticError,OverflowError,TypeError,
                FloatingPointError,np.linalg.LinAlgError) as exc:
            row.update(status='EvaluationError',error=repr(exc))
        result.append(row)
        print(model,case,'factor',factor,row['status'],
              (row.get('screen') or {}).get('scalar_max_tev'),flush=True)
    return {'model':model,'case':case,'historical_config':str(history['cfg_path'].relative_to(bulk.ROOT)),
            'eft_beta':str(beta_path.relative_to(bulk.ROOT)),
            'eft_beta_sha256':sha256(beta_path.read_bytes()).hexdigest(),
            'rows':result,'physics_status':'candidate general-quartic RG, no quadratic extension or complete EFT matching'}


def legacy_probe(model='A_-4', case=LEGACY_CASE, *, beta_campaign='overnight_six_v1_bulk48',extra_real=-0.01):
    registry=bulk.read(bulk.REGISTRY)
    entry,new,reused=bulk.inventory(registry,model)
    if reused:raise ValueError('Not a legacy unsupported-quadratic case')
    history=bulk.historical(case,model,entry)
    if history['uv_on']:raise ValueError('Legacy diagnostic initially requires frozen UV')
    payload,beta_path=load_eft_beta(model,beta_campaign)
    vec,layout,unpack=bulk.uv_base(history['cfg'],entry)
    uv=unpack(vec,layout,mu_gev=1e5)
    q0=tuple(complex(extra_real) for _ in new)
    eft,masses,q,neft=bulk.run_eft(uv,q0,payload,model,new,())
    screen=screen_details(uv,eft,masses,q0,q)
    kernel,matching_path,kind=bulk.matching(model)
    if kind!='legacy-common-scale':raise ValueError('Matching kernel is not legacy')
    result={'model':model,'case':case,'status':'ScreenRejected',
            'matching_source':str(matching_path.relative_to(bulk.ROOT)),
            'matching_sha256':sha256(matching_path.read_bytes()).hexdigest(),
            'screen':screen,'eft_nfev':int(neft)}
    if screen['passed']:
        parity=legacy_parity(kernel,uv,eft)
        result['parity']=parity
        result['analysis']=comparison_metrics(complex_matrix(parity['reference_c5']),
                                             complex_matrix(parity['fast_c5']))
        result['status']='ParityPassed' if parity['passes'] else 'ParityFailed'
    return result


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--beta-campaign',default='overnight_six_v1_bulk48')
    ap.add_argument('--models',nargs='+',default=['A_-2'],help='Full-UV models to probe')
    ap.add_argument('--case',default=DEFAULT_CASE,choices=tuple(bulk.STUDIES))
    ap.add_argument('--factors',nargs='+',type=float,default=list(DEFAULT_FACTORS))
    ap.add_argument('--extra-real',type=float,default=-0.01)
    ap.add_argument('--skip-legacy',action='store_true')
    args=ap.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):ap.error('Invalid run ID')
    if len(set(args.factors))!=len(args.factors) or any(not math.isfinite(f) or f<=0 for f in args.factors):
        ap.error('Factors must be unique finite positive values')
    out=bulk.ROOT/'output/studies/general_scalar_extended/diagnostics'/args.run_id
    if out.exists():ap.error('Output exists: choose a new --run-id')
    records=[]
    for model in args.models:
        records.append(scalar_probe(model,args.case,beta_campaign=args.beta_campaign,
                                    factors=args.factors,extra_real=args.extra_real))
    legacy=None
    if not args.skip_legacy:
        legacy=legacy_probe(beta_campaign=args.beta_campaign,extra_real=args.extra_real)
        print('A_-4 parity',legacy['status'],
              legacy.get('analysis',{}).get('relative_max_error'),flush=True)
    out.mkdir(parents=True,exist_ok=False)
    artifact={'format':'t3-scalar-boundary-and-legacy-parity-probe-v1',
              'scalar_probes':records,'legacy_parity_probe':legacy,
              'read_only_physics_inputs':True,
              'no_changes_to_screen_thresholds_or_matching':True,
              'not_refitted':True}
    (out/'probe.json').write_text(json.dumps(artifact,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    cols=['model','case','scalar_mass_factor','status','uv_mS1Sq_gev2','uv_mS2Sq_gev2',
          'scalar_min_tev','scalar_max_tev','rhoF','yukawa_sigma','quartic_max',
          'failed_checks','uv_nfev','eft_nfev','error']
    with (out/'scalar_boundary.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=cols);writer.writeheader()
        for probe in records:
            for row in probe['rows']:
                s=row.get('screen') or {}
                writer.writerow({**{k:row.get(k) for k in cols},
                         **{k:s.get(k) for k in ('scalar_min_tev','scalar_max_tev','rhoF','yukawa_sigma','quartic_max')},
                         'failed_checks':';'.join(s.get('failed_checks',[]))})
    print('Saved',out/'probe.json',flush=True)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
