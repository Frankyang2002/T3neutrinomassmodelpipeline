#!/usr/bin/env python3
"""Read-only source-guarded all-16 first-fermion-threshold LLSS Gaussian bridge.

Run from the repository root. All outputs are written under the existing owned
lo_certification_direct_uv directory, never to production UV/EFT inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.lo_certification.direct_uv.gaussian_amplitude import batch,AmplitudeError


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=ROOT)
    parser.add_argument('--input-dir',type=Path,default=None)
    parser.add_argument('--output-dir',type=Path,default=None)
    parser.add_argument('--summary-only',action='store_true',
                        help='Omit full source-pair 3x3xnsxns tensors, keeping all norms')
    args=parser.parse_args()
    root=args.project_root.resolve()
    allowed=(root/'output/general_observable/lo_certification_direct_uv').resolve()
    out=(args.output_dir or allowed/'gaussian_amplitude_01').resolve()
    try:
        if not out.is_relative_to(allowed):
            raise AmplitudeError('Output must stay in owned lo_certification_direct_uv directory')
        summary=out/'all_16_gaussian_amplitude_status.json'
        if summary.exists() or any((out/'models').glob('*.json')) or (out/'all_16_gaussian_amplitude_summary.csv').exists():
            raise AmplitudeError('Output already exists: use a NEW --output-dir, never overwrite')
        result=batch(root,args.input_dir,include_tensors=not args.summary_only)
        modelout=out/'models';modelout.mkdir(parents=True,exist_ok=True)
        rows=[]
        for r in result['models']:
            entry=dict(r)
            output=entry.pop('output',None)
            if output is not None:
                p=modelout/(entry['model']+'.json')
                p.write_text(json.dumps({'model':entry['model'],'benchmark_sha256':entry.get('frozen_benchmark_sha256'),
                    'source_candidate_only':True,**output},indent=2,allow_nan=False)+'\n',encoding='utf-8')
                entry['amplitude_file']=str(p.relative_to(root))
                norm=output['total_charge_allowed_coefficient_norm_gev_minus1']
            else:norm=None
            rows.append({'model':r['model'],'status':r['status'],
                         'fermion_type':r.get('fermion_type',''),
                         'heavy_fermion_states':output['heavy_fermion_state_count'] if output else '',
                         'LLSS_total_tensor_norm_gev_minus1':norm if norm is not None else '',
                         'blocked_reason':r.get('reason',''),
                         'physically_certified':False})
            result['models'][len(rows)-1]=entry
        summary.parent.mkdir(parents=True,exist_ok=True)
        summary.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        with (out/'all_16_gaussian_amplitude_summary.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    except (AmplitudeError,ValueError,KeyError,TypeError,OSError,RuntimeError) as exc:
        print('LLSS Gaussian amplitude bridge BLOCKED: '+str(exc),file=sys.stderr)
        return 2
    print('Formal all-heavy-fermion LLSS at INPUT scale (not muF matched): {}/16, blocked {}/16'.format(
        result['tree_amplitude_computed_models'],result['blocked_models']))
    print('Physical source-to-Weyl dictionary certified: 0/16; hard one-loop matching: NOT INCLUDED')
    print('Report: '+str(summary))
    return 0 if not result['blocked_models'] else 2


if __name__=='__main__':raise SystemExit(main())
