#!/usr/bin/env python3
"""Six-model, conditional extended-C5 -> SM/Weinberg low-energy batch.

Use --self-test, --preflight, --pilot or --all; never modifies production.
Failures are written per-model and the process returns nonzero. The matching
is a one-step sudden-decoupling approximation, NOT certified hierarchical EFT.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import csv
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.experimental import T3ExtendedLowEnergyStudy as LE
from Numerical.experimental import T3ExtendedConventionLedger as LEDGER
from Numerical.experimental import T3RGAssistedExtendedC5 as BR

DEFAULT_SCENARIO=ROOT/'Numerical/configs/t3_rg_assisted_c5_convention_hypothesis.json'
DEFAULT_STUDY=ROOT/'Numerical/configs/t3_extended_low_energy_study.json'


def write_json(path:Path,obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(path)


def run(args:argparse.Namespace)->int:
    if args.self_test:
        BR.self_test();LE.self_test()
        from scripts.TestT3ExtendedLowEnergyBatch import self_test
        self_test()
        return 0
    config=json.loads(args.scenario.read_text(encoding='utf-8-sig'))
    scenario=BR.read_scenario(config)
    study=LE.load_study(json.loads(args.study.read_text(encoding='utf-8-sig')))
    if not study['allow_unverified_sudden_decoupling']:
        raise ValueError('Explicit allow_unverified_sudden_decoupling=true required; cannot silently treat one-step matching as hierarchical')
    native=args.native_source
    ledger=LEDGER.audit(ROOT,native_source=native)
    out=args.output_dir.resolve()
    write_json(out/'physical_dictionary_gate.json',ledger)
    print('PRECHECK: six native C5 reductions and UV/Matchete scalar field words PASS; absolute dictionary BLOCKED',flush=True)
    if args.preflight:
        print('Provenance:',out/'physical_dictionary_gate.json')
        return 0
    models=['A_p0'] if args.pilot else list(LE.MODELS)
    buffers=[study['buffers'][0]] if args.pilot else study['buffers']
    rows=[]
    summary={}
    for i,model in enumerate(models,1):
        try:
            result=LE.run_model(model,scenario,study,root=ROOT,buffers=buffers,native_source=native)
            path=out/'models'/f'{model}.json'
            write_json(path,result)
            for ep in result['matching_scale_scan']:
                for phase in ep['phase_scenarios']:
                    obs=phase['low_energy_observables']
                    row={
                        'model':model,'status':'CONDITIONAL_PASS',
                        'mu_match_gev':ep['mu_match_gev'],
                        'endpoint_buffer':ep['endpoint_buffer'],
                        'phase_radians':phase['extra_relative_phase_radians'],
                        'mu_low_gev':ep['mu_low_gev'],
                        'm1_ev':obs['masses_ev'][0],'m2_ev':obs['masses_ev'][1],
                        'm3_ev':obs['masses_ev'][2],
                        'delta_m21_sq_ev2':obs['delta_m21_sq_ev2'],
                        'delta_m31_sq_ev2':obs['delta_m31_sq_ev2'],
                        'ordering':obs['ordering'],
                        'direct_extra_fraction_at_low':phase['direct_extra_fraction_at_low'],
                        'seeded_minus_zero_seed_fraction_at_low':phase['seeded_minus_zero_seed_fraction_at_low'],
                        'c5_linearity_residual':phase['c5_linearity_residual'],
                        'fermion_first_mass_hierarchy':ep['threshold_diagnostics']['fermion_first_hierarchy_in_mass_proxies'],
                        'report':str(path),
                    }
                    rows.append(row)
            print(f'[{i}/{len(models)}] CONDITIONAL PASS {model}: {len(result["matching_scale_scan"])} matching scales; {len(result["matching_scale_scan"][0]["phase_scenarios"])} phases',flush=True)
        except Exception as exc:
            row={'model':model,'status':'BLOCKED','error':type(exc).__name__+': '+str(exc)}
            rows.append(row)
            print(f'[{i}/{len(models)}] BLOCKED {model}: {row["error"]}',file=sys.stderr,flush=True)
        summary={'status':'ConditionalLowEnergyPartialOrBlocked','models_requested':models,
                 'n_models_completed':len({r['model'] for r in rows if r['status']=='CONDITIONAL_PASS'}),
                 'n_models_blocked':len({r['model'] for r in rows if r['status']=='BLOCKED'}),
                 'model_results':rows,'physical_dictionary_gate':str(out/'physical_dictionary_gate.json'),
                 'physical_matchete_rgbeta_dictionary_certified':False,
                 'hierarchical_matching_certified':False,
                 'final_neutrino_mass_predictions_certified':False,
                 'final_SM_Weinberg_RGE_performed_with_one_step_matching':True,
                 'formula':'One-step native reduced C5 evaluated with explicitly unverified RGBeta/Matchete conversion at a common scale above all heavy mass proxies, then SM+Weinberg RGEs to mu_low',
                 'study':study['raw'],'scenario':config,
                 'generated_at_utc':datetime.now(timezone.utc).isoformat()}
        if summary['n_models_completed']==len(models) and summary['n_models_blocked']==0:
            summary['status']='AllRequestedConditionalLowEnergyStudiesPassed_NotHierarchical'
        write_json(out/'summary.json',summary)
    csvpath=out/'low_energy_scan.csv'
    good=[r for r in rows if r['status']=='CONDITIONAL_PASS']
    if good:
        with csvpath.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(good[0]))
            w.writeheader();w.writerows(good)
    print('Summary:',out/'summary.json')
    print('CSV:',csvpath if good else 'Not generated')
    print('PHYSICAL STATUS: all numerical mass/mixing values are conditional. EFT1/hierarchical matching and absolute coupling conversions are NOT CERTIFIED.')
    return 0 if summary.get('n_models_blocked',len(models))==0 and summary.get('n_models_completed')==len(models) else 4


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    action=ap.add_mutually_exclusive_group(required=True)
    for flag in ('self-test','preflight','pilot','all'):
        action.add_argument('--'+flag,action='store_true')
    ap.add_argument('--scenario',type=Path,default=DEFAULT_SCENARIO)
    ap.add_argument('--study',type=Path,default=DEFAULT_STUDY)
    ap.add_argument('--native-source',type=Path,default=None)
    ap.add_argument('--output-dir',type=Path,default=LE.OUT)
    return run(ap.parse_args())

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        print('BLOCKED:',type(exc).__name__+':',exc,file=sys.stderr)
        raise SystemExit(3)
