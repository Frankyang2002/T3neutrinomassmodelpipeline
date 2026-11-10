#!/usr/bin/env python3
"""One-command, source-gated A(-4)/D(-2) consolidated numerical handoff.

Reuses existing complete renormalisable RG and raw/Laurent Matchete outputs.
Does NOT relabel unsubtracted formal expressions as finite EFT1 Wilsons.
When source-derived finite C5(MS) is present and passes the full matching
contract, runs the existing one-loop SMEFT solver and neutrino observables.
Otherwise emits ONE durable blocking report, not fictitious predictions.

Example:
  python -X utf8 -m scripts.RunT3ThesisFullGeneralConsolidated --all
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import re
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.parallel_t3.complete.full_general_completion import (
    MODELS, STAGES, REPORT_FORMAT, CompletionBlocked, existing_evidence,
    verify_completed_c5_document, finish_smeft, digest, put_json,
)
from Numerical.parallel_t3.background.general_batch1_api import read_json,encode_native
from Numerical.parallel_t3.complete.hard_llhh_source import (
    numerical_unsubtracted_source,HardSourceError,
)

B1=ROOT/'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1'
B2=ROOT/'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1'
LAURENT=ROOT/'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2'
C5=ROOT/'output/parallel_t3/full_general_complete_input'
OUT=ROOT/'output/parallel_t3/full_general_consolidated'
CONFIG=ROOT/'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json'


def run_model(model:str, args:argparse.Namespace, target:dict, vev:float, output:Path)->dict:
    statuses={stage:'NotYetCalculated' for stage in STAGES}
    report={'format':REPORT_FORMAT,'model':model,'created_utc':datetime.now(timezone.utc).isoformat(),
            'status':'Blocked','stages':statuses,
            'one_loop_physical_C5_provided':False,'full_general_Yukawa_fit_performed':False,
            'analytic_triangle_fallback_used':False,
            'scalar_bubble_only_fallback_used':False,
            'existing_symbolic_matches_reexecuted':False}
    evidence=None
    try:
        evidence=existing_evidence(model,args.batch1/model,args.raw/model,args.laurent/model)
        statuses.update(evidence.stage_statuses)
        report['inputs_sha256']=evidence.sources_sha256
        report['scale_hierarchy_GeV']={'UV':1e7,'F':evidence.mf_gev,'S':evidence.ms_gev,'low':100}
        report['all_running_coupling_names_at_S']=sorted(evidence.rg_at_ms)
        report['scalar_candidate_mass_squared_eigenvalues_GeV2']=\
            evidence.scalar_at_ms['candidate_spherical_scalar_quadratic_spectrum']['mass_squared_eigenvalues_gev2']
        report['scalar_handover_native_at_S']=\
            encode_native({k:v for k,v in evidence.rg_at_ms.items() if
                           k.startswith('lambda') or k.startswith('mS') or k.startswith('muT3')})
        # This is a strictly unrenormalised Matchete hard-source diagnostic,
        # NOT a substitute for the required full finite C5.
        raw_ledger=read_json(args.raw/model/'one_loop_operator_ledger.json')
        at_f=read_json(args.batch1/model/'fermion_matching_input.json')
        formal=numerical_unsubtracted_source(
            model,raw_ledger,at_f['native_uv_rgbeta_couplings_at_muF'],
            mu_f=evidence.mf_gev)
        put_json(output/model/'formal_F_LLHH_not_a_Wilson.json',formal)
        report['formal_raw_LLHH_diagnostic_file']='formal_F_LLHH_not_a_Wilson.json'
        report['formal_raw_LLHH_diagnostic_status']=formal['status']
        path=args.c5_root/model/'renormalized_c5_at_MS.json'
        if not path.is_file():
            raise CompletionBlocked(
                'No source-derived finite full-general C5 at mu_S. '
                'Saved Laurent WXF contains UV/EFT1 unsubtracted terms, not '
                'numerical Wilsons; full fermion matching subtraction, EFT1 '
                'Wilson RG, and joint scalar matching remain unimplemented.')
        boundary, provenance=verify_completed_c5_document(path,evidence)
        statuses.update({
            'F_finite_one_loop_Wilson':'ProvidedSourceDerivedFiniteBoundary',
            'EFT1_full_Wilson_RGE':'ProvidedSourceDerivedCompleteEvolution',
            'S_joint_full_one_loop':'ProvidedSourceDerivedJointScalarMatching'})
        results=finish_smeft(evidence,boundary,vev=vev,target=target)
        statuses['SMEFT_RG']='ComputedOneLoopSMWeinberg'
        statuses['nu_observables']='ComputedTakagiPMNSAndFiveGaussianResiduals'
        report['C5_MS_input_sha256']=digest(path)
        report['results']=results
        report['one_loop_physical_C5_provided']=True
        report['status']='BenchmarkNeutrinoObservablesComputed_FitNotPerformed'
    except (CompletionBlocked, FileNotFoundError, ValueError, OSError,
            ImportError, RuntimeError, KeyError, TypeError, HardSourceError) as error:
        report['status']='Blocked_MissingFullGeneralMatchingOrNumericHandover'
        report['reason']=str(error)
        if args.debug:
            report['traceback_tail']=traceback.format_exc(limit=8)
    put_json(output/model/'completion_report.json',report)
    print(f'{model}: {report["status"]}')
    if report['status'].startswith('Blocked'):
        print('  '+report['reason'])
    return report


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=MODELS)
    parser.add_argument('--batch1',type=Path,default=B1)
    parser.add_argument('--raw',type=Path,default=B2)
    parser.add_argument('--laurent',type=Path,default=LAURENT)
    parser.add_argument('--c5-root',type=Path,default=C5)
    parser.add_argument('--output-root',type=Path,default=OUT)
    parser.add_argument('--run-id',default='thesis_full_general_completion_v1')
    parser.add_argument('--config',type=Path,default=CONFIG)
    parser.add_argument('--debug',action='store_true')
    args=parser.parse_args(argv)
    if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is None:
        parser.error('Invalid run ID')
    root=args.output_root/args.run_id
    if root.exists():
        parser.error('Output run-id already exists: use a fresh --run-id to retain prior results')
    try:
        cfg=read_json(args.config)
        if cfg.get('schema')!='t3-thesis-benchmark-two-model-pilot-v1':
            raise ValueError('Wrong thesis configuration')
        if cfg.get('scales_GeV')!={'UV':1e7,'F':1e5,'S':1e3,'low':100}:
            raise ValueError('The fixed thesis scale hierarchy changed')
        vev=float(cfg['vev_GeV'])
        if not 0<vev<1000:
            raise ValueError('Invalid Higgs vev')
        target_path=args.config.parent.parent.parent/cfg['fitting']['target']
        # Absolute paths and project-root paths are supported.
        if not target_path.is_file():
            target_path=ROOT/cfg['fitting']['target']
        target=read_json(target_path)
    except (OSError,ValueError,KeyError,TypeError) as e:
        parser.error('Could not load thesis benchmark and NuFIT inputs: '+str(e))
    rows=[]
    for model in MODELS if args.all else (args.model,):
        rows.append(run_model(model,args,target,vev,root))
    put_json(root/'summary.json',{
        'format':REPORT_FORMAT+'-summary',
        'run_id':args.run_id,
        'total_models':len(rows),
        'benchmark_predictions_computed':sum(x['one_loop_physical_C5_provided'] for x in rows),
        'all_models_full_general_completed':False,
        '18_real_Yukawa_fit_performed':False,
        'models':[{k:x.get(k) for k in ('model','status','reason')} for x in rows],
        'source_principle':'No triangle-only or bubble-only substitution for missing complete matching',
    })
    print('Summary: '+str(root/'summary.json'))
    if any(x['status'].startswith('Blocked') for x in rows):
        print('FULL GENERAL COMPLETION: BLOCKED at source-derived finite Wilson/scalar matching.')
        return 2
    print('Both benchmark predictions available; 18-real-Yukawa refit still required.')
    return 3


if __name__=='__main__':
    raise SystemExit(main())
