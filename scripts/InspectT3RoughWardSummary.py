"""Flatten existing all-16 LLHH Ward diagnostics into a readable CSV/table.

Read-only to Workstreams 1/2; creates only a C-owned CSV.  No thresholds
or new RGEs are evaluated.  Run after RunT3RoughAll16.py.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.threshold.models import ALL_MODELS


def flatten_ward(doc: dict) -> dict:
    model=doc['model']
    ward=doc['ward_SU2xU1']
    alternative=doc['alternative_coordinate_hypotheses_DO_NOT_APPLY'][0]
    mf=doc.get('source_MF_orderzero_witness',{})
    c5=doc.get('stream2_separately_declared_C5',{})
    values={
        'model':model,
        'source_LLHH_norm_GeV_inv':doc['source_LLHH_norm_gev_inverse'],
        'Weinberg_residual':doc['relative_Weinberg_purity_residual'],
        'SU2_T1':ward['SU2_T1']['relative_to_LLHH'],
        'SU2_T2':ward['SU2_T2']['relative_to_LLHH'],
        'SU2_T3':ward['SU2_T3']['relative_to_LLHH'],
        'U1_Y':ward['U1_Y']['relative_to_LLHH'],
        'neutral_vs_LS_C5_relative':doc['C5_component_vs_LS_relative_to_larger'],
        'separate_C5_vs_LS_relative':c5.get('relative_to_C_projector'),
        'separate_C5_vs_neutral_relative':c5.get('relative_to_neutral_component'),
        'best_template_diagnostic_only':alternative['hypothesis_only'],
        'best_template_residual':alternative['relative_residual'],
        'MF_witness_kind':mf.get('availability','not_present_in_older_report'),
        'MF_orderzero_Weinberg_residual':mf.get('relative_Weinberg_purity_residual'),
        'first_investigation_hint':doc.get('diagnostic_summary',{}).get(
            'first_investigation_hint_NOT_A_REPAIR','older_diagnostic_no_classification'),
    }
    values['max_Ward']=max(values[k] for k in ('SU2_T1','SU2_T2','SU2_T3','U1_Y'))
    return values


def main(argv=None):
    p=argparse.ArgumentParser(description='View genuine A-source Ward numerics without nested PowerShell objects')
    p.add_argument('--root',type=Path,default=Path('output/parallel_t3/rough_results/rough_fermion_first_v1'))
    p.add_argument('--csv',type=Path,default=None)
    args=p.parse_args(argv)
    folder=args.root/'_diagnostics'
    records=[]
    for model in ALL_MODELS:
        slug=model.replace('_-','_m').replace('_+','_p')
        f=folder/f'{slug}_llhh_projection.json'
        if not f.is_file():
            print(model, 'DIAGNOSTIC_NOT_AVAILABLE',f)
            continue
        doc=json.loads(f.read_text(encoding='utf-8-sig'))
        if doc.get('model')!=model or doc.get('status')!='SOURCE_DIAGNOSTIC_NOT_MATCHING':
            raise ValueError(f'{f}: source-model/status mismatch')
        record=flatten_ward(doc)
        records.append(record)
        print(f'{model:<5} purity={record["Weinberg_residual"]:.6g} '
              f'SU2=( {record["SU2_T1"]:.4g}, {record["SU2_T2"]:.4g}, '
              f'{record["SU2_T3"]:.4g} )  Y={record["U1_Y"]:.4g} '
              f'best={record["best_template_diagnostic_only"]} '
              f'({record["best_template_residual"]:.5g}) '
              f'MF={record["MF_witness_kind"]}')
    if not records:
        raise SystemExit(f'No LLHH diagnostics in {folder}')
    output=args.csv or (folder/'ward_summary_flat.csv')
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',encoding='utf-8',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    print('FLAT WARD SUMMARY',output)
    return 0

if __name__=='__main__':
    sys.exit(main())
