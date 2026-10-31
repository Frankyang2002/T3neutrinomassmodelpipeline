"""Recompute B/C conditional Matchete restricted C5 under both Yukawa maps.

Invoke from the repository root: python -m scripts.RunT3BCConventionMatchingBatch
Requires the previously generated running UV and native probe reports.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import RunT3BCRunningOneScaleTriangle as runner
from Numerical.experimental import T3BCConventionMatchingBatch as core

OUT=ROOT/'output/studies/t3_bc_convention_matching_batch_v1'
NATIVE=ROOT/'output/studies/t3_bc_majorana_chiral_rotation_v1/native_matchete.json'
ACTION=ROOT/'output/studies/t3_bc_majorana_action_batch_v1/summary.json'


def main():
    os.chdir(ROOT)
    native=runner.read(NATIVE);action=runner.read(ACTION)
    v22=runner.import_source(runner.DIRECT,'t3_batch_v22')
    v18=runner.import_source(runner.V18,'t3_batch_v18')
    evidence=runner.read(runner.V22)
    summary={'format':'t3-bc-convention-matching-batch-v1',
      'status':'TwoConditionalRotationsComputedCrossPackageDictionaryPending',
      'models':{},'physical_C5_certified':False,
      'production_files_modified':False,
      'all_results_conditional_not_physical':True,
      'complete_sequential_EFT_threshold_matching':False}
    OUT.mkdir(parents=True,exist_ok=True)
    for model in ('B_-1','C_-1'):
        tag=runner.model_tag(model)
        uv_path=runner.BC_PARENT/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json'
        archived_path=runner.ROOT/'output/studies/t3_bc_running_one_scale_triangle_v1'/f'{tag}_running_one_scale.json'
        source,boundary=runner.load_bc(model,uv_path)
        archived=runner.read(archived_path)
        runner.require(archived.get('running_bc_source_sha256')==runner.digest(uv_path),
                'archived running endpoint provenance mismatch')
        check=core.check_gate(source,archived,native,action)
        evaluator,reference=runner.archived_source(model,v22,v18,evidence)
        runner.require(archived['frozen_source_parity']['native_direct_source_sha256']==reference['native_direct_source_sha256'],
                'archived direct source digest differs')
        row=core.evaluate(boundary,evaluator,reference['electroweak_vev_reference_gev'],archived)
        row['source_gates']=check
        row['running_input_sha256']=runner.digest(uv_path)
        row['archived_triangle_sha256']=runner.digest(archived_path)
        row['direct_native_source_sha256']=reference['native_direct_source_sha256']
        (OUT/f'{tag}_conventions.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        summary['models'][model]={'status':row['status'],
            'relative_C5_difference':row['relative_C5_same_U_minus_production'],
            'relative_spectrum_difference':row['relative_neutrino_spectrum_same_U_minus_production'],
            'production_masses_ev_NOT_PHYSICAL':row['production_baseline']['neutrino_masses_descending_ev_NOT_PHYSICAL'],
            'same_U_masses_ev_NOT_PHYSICAL':row['same_U_weyl_candidate']['neutrino_masses_descending_ev_NOT_PHYSICAL'],
            'baseline_reproduction_relative':row['recomputed_production_vs_archived_C5_relative'],
            'report':str((OUT/f'{tag}_conventions.json').relative_to(ROOT))}
        print(model,'C5 difference',row['relative_C5_same_U_minus_production'],
              'masses same-U',row['same_U_weyl_candidate']['neutrino_masses_descending_ev_NOT_PHYSICAL'],flush=True)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(summary['status']);print('Output:',OUT)
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr)
        raise
