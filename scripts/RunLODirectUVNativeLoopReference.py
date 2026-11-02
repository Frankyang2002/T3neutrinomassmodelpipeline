#!/usr/bin/env python3
"""All-16 independent exact full-mass UV loop reference, source fail-closed.

Writes only to a NEW Workstream 4 output directory. Reads existing provisional
UV tensors and current source. No benchmark, source, EFT or production output
is modified. Comparison with existing numerical engine is optional and does
not promote the physical normalisation status.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

PROJECT=Path(__file__).resolve().parents[1]
if str(PROJECT) not in sys.path:sys.path.insert(0,str(PROJECT))
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import run_all
from Numerical.lo_certification.direct_uv.native_amplitude_contract import ContractError


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=PROJECT)
    p.add_argument('--input-dir',type=Path)
    p.add_argument('--native-export',type=Path,help='Optional existing live Matchete A_+0/B_-1 operator export')
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--rephase',action='store_true',help='Recalculate independent reference after genuine complex basis changes')
    p.add_argument('--no-production-comparison',action='store_true')
    a=p.parse_args()
    root=a.project_root.resolve()
    owned=(root/'output/general_observable/lo_certification_direct_uv').resolve()
    out=(a.output_dir or owned/'native_uv_reference_01').resolve()
    try:
        if not out.is_relative_to(owned):
            raise ContractError('Output must stay in owned Workstream 4 output tree')
        if out.exists() and any(out.iterdir()):
            raise ContractError('Output directory is not empty; select a NEW output directory')
        report=run_all(root,rephase=a.rephase,compare_production=not a.no_production_comparison,
                       input_dir=a.input_dir,native_export=a.native_export)
        entries=report.pop('reports_by_model')
        out.mkdir(parents=True,exist_ok=True)
        md=out/'models';md.mkdir(exist_ok=True)
        for label,entry in entries.items():
            (md/(label+'.json')).write_text(json.dumps(entry,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        for entry in report['models']:
            entry['report_path']=str((md/(entry['model']+'.json')).relative_to(root))
        summary=out/'all_16_native_uv_reference_status.json'
        summary.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        with (out/'all_16_native_uv_reference_summary.csv').open('w',newline='',encoding='utf-8') as s:
            columns=['model','status','independent_masses_ev','uv_pole_relative',
                     'cancellation_factor','relative_production_difference',
                     'physical_source_normalisation_certified']
            writer=csv.DictWriter(s,fieldnames=columns)
            writer.writeheader()
            for entry in report['models']:
                writer.writerow({key:entry.get(key) for key in columns})
    except (ContractError,ValueError,TypeError,OSError,KeyError,RuntimeError) as exc:
        print('WS4 NATIVE UV REFERENCE BLOCKED: '+str(exc),file=sys.stderr)
        return 2
    print('All-16 independent full-mass UV references computed conditionally: 16/16')
    print('Method: Hermitian mass-matrix spectral function, no Takagi/heavy-mass expansion')
    print('Physical native Matchete/RGBeta normalisation certified: 0/16')
    print('Native scalar portal bare->Inv mapping: unresolved, explicitly traced')
    print('Output: '+str(summary))
    return 0

if __name__=='__main__':raise SystemExit(main())
