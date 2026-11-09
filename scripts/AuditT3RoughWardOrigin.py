"""All-16 source origin audit of already-adapted Stream-2 LLHH tensors.

Reads ONLY C-owned `_adapted_inputs` from the existing all-16 runner.
No UV, MF, EFT1, scalar threshold, or SMEFT RGE is recomputed.
This is a diagnostic and NEVER produces genuine neutrino masses.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.threshold.models import ALL_MODELS, normalize_model
from Numerical.parallel_t3.threshold.rough_origin import FORMAT, audit_pair

DEFAULT=Path('output/parallel_t3/rough_results/rough_fermion_first_v1')


def flatten(record:dict)->dict:
    name=record['model']
    row={'model':name,'status':record['status'],
         'source_wilson_sha256':record.get('source_sha256',{}).get('npz'),
         'MF_origin_interpretation':record.get('origin_interpretation_NOT_proof')}
    if row['status']!='SourceWardOriginDiagnosticOnly':
        row['blocker']=record.get('blocker');return row
    x=record['MS_full_evolved_LLHH']; mf=record['MF_orderzero_LLHH']
    rows=x['charge_sector_fractions']
    ward=x['allowed_HH_SU2_Ward_relative']
    row.update({'MS_charge_category':x['status'],
                'MS_Weinberg_residual':record['MS_Weinberg_purity_relative'],
                'MS_fraction_HH':rows['HH'],
                'MS_fraction_H_Hbar':rows['H_Hbar_plus_Hbar_H'],
                'MS_fraction_Hbar_Hbar':rows['Hbar_Hbar'],
                'MS_HH_Singlet_residual':x['allowed_HH_singlet_residual_relative'],
                'MS_allowed_HH_SU2_T1':ward['T1'],
                'MS_allowed_HH_SU2_T2':ward['T2'],
                'MS_allowed_HH_SU2_T3':ward['T3'],
                'MS_U1_Y_Ward':x['direct_U1_Ward_relative'],
                'MS_U1_Y_sector_reconstruction':x['independently_reconstructed_U1_Ward_relative'],
                'MF_orderzero_availability':mf['availability'],
                'MF_orderzero_charge_category':mf.get('status'),
                'MF_orderzero_HH_singlet_residual':mf.get('allowed_HH_singlet_residual_relative'),
                'MF_orderzero_U1_Y_Ward':mf.get('direct_U1_Ward_relative'),
                'MF_orderzero_LLHH_norm':mf.get('full_real_LLHH_norm_gev_inverse')})
    return row


def parse_args(argv=None):
    p=argparse.ArgumentParser(description='Trace T3 Weinberg Ward failure to SU2 and hypercharge sectors')
    p.add_argument('--root',type=Path,default=DEFAULT,
        help='Output root from RunT3RoughAll16.py (contains _adapted_inputs)')
    p.add_argument('--models',nargs='+',default=None,help='Model labels, default all 16')
    p.add_argument('--output',type=Path,default=None,
        help='Optional results directory; default <root>/_diagnostics/origin_v5')
    return p.parse_args(argv)


def main(argv=None)->int:
    args=parse_args(argv)
    models=[normalize_model(v) for v in (args.models or ALL_MODELS)]
    if len(models)!=len(set(models)):raise SystemExit('Duplicate models in --models')
    target=args.output or args.root/'_diagnostics'/'origin_v5'
    target.mkdir(parents=True,exist_ok=True)
    rows=[]
    for model in models:
        slug=model.replace('_-','_m').replace('_+','_p')
        source=args.root/'_adapted_inputs'/slug/'wilson.json'
        try:
            result=audit_pair(source,model=model)
            x=result['MS_full_evolved_LLHH']
            mf=result['MF_orderzero_LLHH']
            print(f'{model:<5} {x["status"]:<45} '
                  f'HH={x["charge_sector_fractions"]["HH"]:.4g} '
                  f'HHbar={x["charge_sector_fractions"]["H_Hbar_plus_Hbar_H"]:.4g} '
                  f'HbarHbar={x["charge_sector_fractions"]["Hbar_Hbar"]:.4g} '
                  f'HH_SU2={max(x["allowed_HH_SU2_Ward_relative"].values()):.4g} '
                  f'MF={mf.get("status",mf["availability"])}',flush=True)
        except Exception as exc:
            result={'format':FORMAT,'model':model,
                    'status':'SOURCE_AUDIT_BLOCKED',
                    'source_manifest':str(source),
                    'blocker':f'{type(exc).__name__}: {exc}'}
            print(model,'SOURCE_AUDIT_BLOCKED',result['blocker'],flush=True)
        file=target/f'{slug}_origin.json'
        file.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        rows.append(flatten(result))
    keys=list(dict.fromkeys(key for r in rows for key in r))
    with (target/'origin_summary.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=keys,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    (target/'origin_summary.json').write_text(json.dumps({
        'format':'t3-rough-llhh-origin-summary-v5',
        'diagnostic_only':True,'models':rows,
        'count':len(rows),
        'missing_or_invalid_input':sum(r['status']=='SOURCE_AUDIT_BLOCKED' for r in rows),
        'source_A_B_files_never_modified':True,
        'no_source_tensor_repair_or_threshold_matching':True,
    },indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('WARD ORIGIN SUMMARY:',target/'origin_summary.csv',flush=True)
    return 1 if any(r['status']=='SOURCE_AUDIT_BLOCKED' for r in rows) else 0


if __name__=='__main__':sys.exit(main())
