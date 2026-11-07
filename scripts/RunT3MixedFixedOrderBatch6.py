"""Batch 6: source-guarded three-heavy-flavour hard scalar C5; pole gate explicit.

ONLY B_-1/C_-1. Native scalar hard Match must have succeeded first. This
program does NOT manufacture an opposite inherited pole or declare a complete
fixed-order neutrino mass without source-specific counterterms.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any

from Numerical.running.T3MixedBatch6Numerics import evaluate_hard
from scripts.RunT3MixedScalarHardBatch6 import (
    MODELS, OUTPUT as HARD_OUTPUT, ROOT, RUNNING_ROOT, SYNC,
    preflight, paths, tag,
)

OUTPUT = ROOT/'output/general_sequential_pilot/mixed_fixed_order_batch6'


def load(path: Path) -> dict[str, Any]:
    result=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(result,dict):
        raise ValueError(f'Expected JSON object: {path}')
    return result


def save(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def verify_hard(model: str, hard_run: str) -> tuple[Path,dict]:
    folder=HARD_OUTPUT/hard_run/tag(model)
    report_file=folder/'hard_scalar_report.json'
    source=folder/'hard_c5_candidate.txt'
    audit=folder/'hard_scalar_input_audit.json'
    if not all(p.is_file() for p in (report_file,source,audit)):
        raise ValueError('Native scalar-only Match not exported for '+model+': '+str(folder))
    r=load(report_file); a=load(audit)
    if (r.get('status')!='MixedScalarHardCandidateExported' or
        r.get('model')!=model or r.get('native_scalar_one_loop_matched') is not True or
        r.get('canonical_EOMEFT_input') is not True or
        r.get('single_heavy_eigenvalue_template') is not True or
        r.get('fermion_matching_recomputed') is not False or
        r.get('full_pole_cancellation_certified') is not False or
        a.get('Model')!=model or
        a.get('OriginalFermionSourceUnmodified') is not True or
        a.get('All17NativeQuarticsFromOriginalBuilder') is not True or
        a.get('MFSourceMode') not in ('SavedUnindexedSingleMassTemplate',
                                      'IndexedMFConvertedToSingleMassTemplate') or
        r.get('mf_source_mode') != a.get('MFSourceMode') or
        r.get('independent_three_mass_native_matching') is not False or
        a.get('IndependentThreeMassNativeMatching') is not False):
        raise ValueError(model+': native scalar matching certificate not supported')
    if len(source.read_text(encoding='utf-8-sig'))<100:
        raise ValueError(model+': missing candidate coefficients')
    return source,{'hard_report_sha256':sha256(report_file.read_bytes()).hexdigest(),
                   'hard_candidate_sha256':sha256(source.read_bytes()).hexdigest(),
                   'hard_audit_sha256':sha256(audit.read_bytes()).hexdigest()}


def main(argv: list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id',default='mixed_fixed_01')
    parser.add_argument('--hard-run',default='mixed_hard_01')
    parser.add_argument('--models',default='B_-1,C_-1')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--conditional-quartic-map',action='store_true',help='UNCERTIFIED Matchete=-RGBeta map for two HH quartics')
    parser.add_argument('--require-pole-cancellation',action='store_true',help='Fail closed; independent renormalised inherited/soft source not available')
    parser.add_argument('--muS-gev',type=float,default=None)
    a=parser.parse_args(argv)
    names=tuple(x.strip() for x in a.models.split(','))
    if not names or len(set(names))!=len(names) or any(m not in MODELS for m in names):
        parser.error('Only unique B_-1,C_-1 are allowed')
    for label in (a.run_id,a.hard_run):
        if not re.fullmatch('[A-Za-z0-9_-]+',label):
            parser.error('Run identifiers must contain only letters/digits/_/-')
    sources=[]
    for model in names:
        try:
            ready=preflight(model)
            src,fps=verify_hard(model,a.hard_run)
            sources.append((model,ready,src,fps))
            print(model+': READY, native X/Y scalar hard C5, F threshold reused')
        except (ValueError,KeyError,OSError,TypeError) as exc:
            print(model+': BLOCKED: '+str(exc),file=sys.stderr)
    if len(sources)!=len(names):
        return 2
    # A genuinely inherited, source-projected and scheme-renormalised Weinberg
    # pole is not present in the saved B/C F-only tree. Archived raw F-loop
    # kernels are NOT renormalised soft counterterms and need independent work.
    if a.require_pole_cancellation:
        print('BLOCKED: B/C inherited first-F Weinberg/soft counterterm pole '+
              'is not independently normalised/projected in the current source. '+
              'No pole cancellation certificate can be issued.',file=sys.stderr)
        return 2
    if a.check_only:
        print('READY: native hard source available; inherited/soft pole cancellation '+
              'NOT VERIFIED. No physics or prior ordinary models recomputed.')
        return 0
    folder=OUTPUT/a.run_id
    config={'models':list(names),'hard_run':a.hard_run,
            'conditional_quartic_map':a.conditional_quartic_map,
            'muS_gev':a.muS_gev,
            'input_hashes':{model:dict(fps,source_fingerprints=ready['source_fingerprints'])
               for model,ready,_,fps in sources}}
    if folder.exists():
        if not a.resume or not (folder/'settings.json').is_file() or load(folder/'settings.json')!=config:
            raise ValueError('Output exists or resume settings/hash mismatch; choose a fresh run ID')
    else:
        if a.resume: raise ValueError('Cannot resume a nonexistent run')
        folder.mkdir(parents=True)
        save(folder/'settings.json',config)
    succeeded=0
    for model,ready,source,fps in sources:
        dest=folder/tag(model)
        result_path=dest/'result.json'
        if a.resume and result_path.is_file():
            prior=load(result_path)
            if prior.get('status')!='NativeHardThreeFlavorPoleSeparated_NotInheritedPoleCancelled':
                raise ValueError(model+': prior result unexpected, refusing overwrite')
            print(model+': reused completed Batch 6 hard-only numerical assembly')
            succeeded+=1
            continue
        try:
            ps=paths(model)
            result=evaluate_hard(model,source,load(ps['running']),load(ps['parameters']),
                conditional_quartic_map=a.conditional_quartic_map,muS=a.muS_gev)
            result['source_fingerprints']=dict(fps,source_fingerprints=ready['source_fingerprints'])
            result['inherited_pole_source']=None
            result['pole_cancellation_certificate']=None
            save(result_path,result)
            succeeded+=1
            print(model+': 3-flavour hard C5 computed; hard pole separated; inherited pole NOT CANCELLED')
            print('  hard-only conditional singular values [eV]: '+
                  str(result['hard_finite_only_conditional_takagi_masses_ev']))
        except Exception as exc:
            save(dest/'error.json',{'status':'Blocked','model':model,'reason':str(exc)})
            print(model+': BLOCKED: '+str(exc),file=sys.stderr)
    save(folder/'summary.json',{
        'status':'NativeHardNumericsOnly_InheritedPoleUnverified' if succeeded==len(names) else 'PartialOrBlocked',
        'models':list(names),'native_hard_three_flavour_successes':succeeded,
        'ordinary_models_rerun':False,'fermion_matching_rerun':False,
        'hard_pole_extracted':succeeded==len(names),
        'inherited_pole_cancellation_certified':False,
        'complete_fixed_order_C5_assembled':False,'physical_C5_certified':False,
        'conditional_quartic_map':a.conditional_quartic_map,
        'notes':'No inherited first-threshold renormalised counterterm projection available; never infer cancellation from hard pole alone.'})
    return 0 if succeeded==len(names) else 2


if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,KeyError,OSError,TypeError) as exc:
        print('BATCH6 ERROR: '+str(exc),file=sys.stderr)
        raise SystemExit(2)
