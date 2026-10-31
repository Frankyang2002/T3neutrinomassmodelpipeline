"""One-point end-to-end operational smoke test of existing ordinary-T3 ScanCLI.

Does not validate physical accuracy or change production physics code.
Run at repository root: python scripts/maintenance/SmokeTestFullPipeline.py
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output' / 'audit' / 'end_to_end_smoke'
SKIP_DIRS = {'.git','.venv','venv','__pycache__','node_modules','.pytest_cache','output','Reports','archive','archives','legacy'}


def candidates():
    """Find existing scan configurations, without searching generated reports."""
    for current, dirs, files in os.walk(ROOT):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
        rel = Path(current).relative_to(ROOT)
        if len(rel.parts) > 5:
            dirs[:] = []
            continue
        for name in sorted(files):
            if not name.lower().endswith('.json'):
                continue
            p = Path(current) / name
            try:
                obj = json.loads(p.read_text(encoding='utf-8-sig'))
            except (OSError,UnicodeError,ValueError):
                continue
            if isinstance(obj,dict) and all(k in obj for k in ('representation','base_state','scan','files','output','scales')) and isinstance(obj.get('files'),dict) and {'uv_rgbeta','eft1_rgbeta','final_weinberg','oscillation_target'} <= set(obj['files']):
                yield p


def one_point(raw):
    result=deepcopy(raw)
    scan=result['scan']
    mode=str(scan.get('mode','')).lower()
    if mode == 'grid':
        values=scan['values']
        result['scan']={**scan,'mode':'grid','values':{k:[v[0]] for k,v in values.items()},'fail_fast':True}
    elif mode == 'random':
        result['scan']={**scan,'mode':'random','n_points':1,'seed':int(scan.get('seed') or 42),'fail_fast':True}
    else:
        raise ValueError('Unsupported scan mode '+mode)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,help='Existing production ScanCLI JSON config (optional; auto-detect otherwise)')
    args=parser.parse_args()
    os.chdir(ROOT)
    sys.path.insert(0,str(ROOT))
    OUT.mkdir(parents=True,exist_ok=True)
    result={'status':'NOT_RUN','started_utc':datetime.now(timezone.utc).isoformat(),'scope':'Operational one-point UV -> EFT1 -> C5 -> SM+Weinberg -> oscillation observables; not physical certification','configuration':None,'discovery':[]}
    path=OUT/'smoke_result.json'
    try:
        from Numerical.fitting import ScanCLI
        configs=[args.config.resolve()] if args.config else list(candidates())
        if not configs:
            raise FileNotFoundError('No existing ScanCLI config found. Supply --config path/to/your_scan_config.json')
        selected=None
        for p in configs:
            try:
                original=ScanCLI.load_scan_config(p)
                candidate=replace(original,raw=one_point(original.raw))
                # Check all existing native input files and their metadata before running ODEs.
                ScanCLI.validate_scan_inputs(candidate)
                selected=candidate
                break
            except Exception as exc:
                result['discovery'].append({'config':str(p),'reason':f'{type(exc).__name__}: {exc}'[:900]})
        if selected is None:
            raise RuntimeError('No eligible production scan config passed validation; inspect discovery in smoke_result.json. For a known config rerun with --config.')
        result['configuration']=str(selected.source_path)
        selected=replace(selected,output_path=OUT/'single_point_scan.json',checkpoint_path=OUT/'single_point_scan.checkpoint.json')
        print('Selected existing config:',selected.source_path,flush=True)
        print('Running exactly one benchmark point through production ScanCLI...',flush=True)
        scan=ScanCLI.execute_scan(selected)
        if len(scan.points)!=1:
            raise AssertionError(f'Expected one benchmark point, got {len(scan.points)}')
        point=scan.points[0]
        result['scan_point']={'status':point.status,'error':point.error,'parameters':point.parameters,'chi2':point.chi2,'ordering':point.ordering,'prediction':point.prediction}
        if not point.succeeded:
            raise RuntimeError('Numerical benchmark failed: '+str(point.error))
        if not isinstance(point.prediction,dict) or not point.prediction:
            raise AssertionError('No observable prediction dictionary returned')
        if not all(isinstance(v,(int,float)) and math.isfinite(float(v)) for v in point.prediction.values()):
            raise AssertionError('Nonfinite or nonnumeric observable')
        if not (OUT/'single_point_scan.json').is_file():
            raise FileNotFoundError('Production scan did not write its output JSON')
        result['status']='PASS_END_TO_END_OPERATIONAL'
        print('PASS: UV -> EFT1 -> matched C5 -> low-energy SM/Weinberg -> numerical observables',flush=True)
        return 0
    except Exception as exc:
        result['status']='FAIL_END_TO_END'
        result['error']=f'{type(exc).__name__}: {exc}'
        result['traceback']=traceback.format_exc()
        print(result['error'],file=sys.stderr,flush=True)
        return 1
    finally:
        result['finished_utc']=datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(result,indent=2,default=str)+'\n',encoding='utf-8')
        print('Report:',path,flush=True)

if __name__=='__main__':
    sys.exit(main())
