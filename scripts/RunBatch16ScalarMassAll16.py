"""Audit all 16 Batch15 scalar thresholds with explicit real mass Hessians.

Project path: scripts/RunBatch16ScalarMassAll16.py
Run: python scripts/RunBatch16ScalarMassAll16.py [--model A_+0]
The spectra are convention-dependent candidate physical inputs, NOT certified
Matchete matching inputs. No Wilson matching or SMEFT prediction is made.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.experimental.Batch16ScalarMassHessians import audit_mass_spectrum

MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
          'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
          'E_+0','E_+2','E_-2','E_-4')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',choices=MODELS)
    ap.add_argument('--input',type=Path,default=ROOT/'output/batch15/scalar_threshold')
    ap.add_argument('--out',type=Path,default=ROOT/'output/batch16/scalar_mass_audit')
    args=ap.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    summary={
        'format':'t3-batch16-scalar-mass-summary-v1',
        'status':'CandidateMassHessianAuditOnly',
        'physical_c5_certified':False,
        'models':{},
    }
    for model in (args.model,) if args.model else MODELS:
        try:
            source=args.input/f'{model}_scalar_threshold_handoff.json'
            report=json.loads(source.read_text(encoding='utf-8-sig'))
            if report.get('model')!=model:
                raise ValueError('Input model does not match filename')
            result=audit_mass_spectrum(report)
            (args.out/f'{model}_mass_hessian.json').write_text(
                json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            r=result['conventions']['rgbeta_potential_plus']
            summary['models'][model]={
                'passed':True,
                'quadratic_kind':result['quadratic_kind'],
                'tachyon_free_candidate':r['tachyon_free_candidate'],
                'eigenmass_hierarchy_candidate':r['eigenmass_threshold_hierarchy_candidate'],
                'max_proxy_eigenvalue_shift_gev2':result['maximum_shift_from_diagonal_proxy_gev2'],
            }
            print(f'[{model}] AUDITED | {result["quadratic_kind"]} | '
                  f'tachyon-free={r["tachyon_free_candidate"]} | '
                  f'eigenmass hierarchy={r["eigenmass_threshold_hierarchy_candidate"]} | '
                  f'C5=BLOCKED',flush=True)
        except Exception as exc:
            summary['models'][model]={'passed':False,'error':str(exc)}
            (args.out/f'{model}_error.log').write_text(traceback.format_exc(),encoding='utf-8')
            print(f'[{model}] FAILED: {exc}',flush=True)
    summary['audited']=sum(v['passed'] for v in summary['models'].values())
    summary['failed']=len(summary['models'])-summary['audited']
    summary['candidate_hierarchy_failures']=[k for k,v in summary['models'].items()
        if v['passed'] and not v['eigenmass_hierarchy_candidate']]
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(f'Summary: {args.out/"summary.json"} | audited={summary["audited"]} '
          f'failed={summary["failed"]} | candidate hierarchy failures='
          f'{len(summary["candidate_hierarchy_failures"])}')
    return 0 if not summary['failed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
