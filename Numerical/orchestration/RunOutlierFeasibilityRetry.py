"""Retry only first-pass T3 cases with chi2 above a threshold; preserve original results."""
from __future__ import annotations
import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, discover_cases


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario',default='smallY_smallL')
    p.add_argument('--threshold',type=float,default=100.0)
    p.add_argument('--points',type=int,default=1024)
    p.add_argument('--local-evals',type=int,default=100)
    p.add_argument('--local-starts',type=int,default=5)
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--rerun',action='store_true')
    p.add_argument('--run',action='store_true')
    a=p.parse_args(argv)
    old=PROJECT_ROOT/'Reports'/'output'/'full'/'fit_feasibility'/a.scenario
    out=PROJECT_ROOT/'Reports'/'output'/'full'/'fit_feasibility_retry'/a.scenario
    cases=[]
    for scenario,key,_ in sorted(discover_cases(OUTPUT_ROOT)):
        if scenario!=a.scenario: continue
        f=old/key/'feasibility_summary.json'
        if not f.is_file(): continue
        try:
            chi=float(json.loads(f.read_text(encoding='utf-8'))['best']['chi2'])
        except (ValueError,KeyError,TypeError): continue
        if math.isfinite(chi) and chi>a.threshold:
            cases.append((key,chi))
    if a.limit: cases=cases[:a.limit]
    print(f'{len(cases)} cases above chi2={a.threshold:g}; output: {out}',flush=True)
    for key,chi in cases: print(f'  {key}: {chi:.6g}',flush=True)
    if not a.run:
        print('Dry run; add --run to execute.')
        return 0
    rows=[]
    for index,(key,old_chi) in enumerate(cases,1):
        path=out/key/'feasibility_summary.json'
        status='AlreadyComplete'
        if a.rerun or not path.is_file():
            print(f'[{index}/{len(cases)}] {key}',flush=True)
            command=[sys.executable,'-m','Numerical.orchestration.RunFitFeasibilityRetry',
                     '--scenario',a.scenario,'--model',key,'--points',str(a.points),
                     '--local-evals',str(a.local_evals),'--local-starts',str(a.local_starts),
                     '--seed','20261012','--run']
            proc=subprocess.run(command,cwd=PROJECT_ROOT,check=False)
            status='Success' if proc.returncode==0 and path.is_file() else f'Failed(exit={proc.returncode})'
        row={'scenario':a.scenario,'model':key,'original_chi2':old_chi,'retry_chi2':'',
             'improvement_factor':'','status':status}
        if path.is_file() and not status.startswith('Failed'):
            try:
                summary=json.loads(path.read_text(encoding='utf-8'))
                new=float(summary['best']['chi2'])
                row['retry_chi2']=new
                row['improvement_factor']=old_chi/new if new>0 else math.inf
                row.update(summary['best'].get('prediction',{}))
            except (ValueError,KeyError,TypeError) as exc:
                row['status']=f'InvalidSummary({exc})'
        rows.append(row)
        out.mkdir(parents=True,exist_ok=True)
        fields=list(dict.fromkeys(k for r in rows for k in r))
        with (out/'retry_summary.csv').open('w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=fields);w.writeheader();w.writerows(rows)
        (out/'retry_summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        print(f'  {row["status"]}: {old_chi:.5g} -> {row["retry_chi2"]}',flush=True)
    print('Saved:',out/'retry_summary.csv')
    return int(any(str(r['status']).startswith(('Failed','Invalid')) for r in rows))

if __name__=='__main__':
    raise SystemExit(main())
