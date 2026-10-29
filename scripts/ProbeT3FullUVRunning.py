"""Read-only physics preflight for full-UV T3 general-scalar benchmarks.

Check frozen-successful UV Yukawa seeds under the actual full UV+EFT1 RGEs.
This is NOT a fit and never changes screening or mass boundaries.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from scripts import RunT3GeneralScalar48Campaign as bulk

CASE = 'small_uv_running_5e-9'
BLOCKED_LEGACY = {'A_-4','B_-3','C_-3','E_-4'}
REUSED = {'B_-1','C_-1'}
FORMAT = 't3-fulluv-mass-screen-preflight-v1'


def finite_or_none(v):
    try:
        v = float(v)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError, OverflowError):
        return None


def candidate_vectors(model, frozen_run, hist):
    """Try historical running seed and optionally saved frozen-fit Yukawas.

    Frozen accepted fits are candidates only, not purported UV-running fits.
    """
    yield 'historical_fulluv', hist['seed']
    if not frozen_run:
        return
    for case in ('small_uv_frozen_5e-9',):
        p = bulk.OUT / frozen_run / case / bulk.safe(model) / 'state.json'
        if not p.is_file():
            continue
        doc = bulk.read(p)
        best = doc.get('best') or {}
        params = best.get('parameters') or {}
        if not all(n in params for n in bulk.NAMES):
            continue
        x = np.asarray([params[n] for n in bulk.NAMES], dtype=float)
        if x.shape == (18,) and np.isfinite(x).all():
            yield f'frozen_fit:{frozen_run}', x


def probe(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        raise ValueError('Invalid run-id')
    registry = bulk.read(bulk.REGISTRY)
    models = list(registry['models'])
    if args.models:
        selected = [x.strip() for x in args.models.split(',')]
        unknown = set(selected) - set(models)
        if unknown:
            raise ValueError(f'Unknown models: {sorted(unknown)}')
        models = [m for m in models if m in selected]
    base = bulk.OUT / args.run_id
    if base.exists() and not args.resume:
        raise FileExistsError(f'Output already exists: {base}; use --resume')
    base.mkdir(parents=True, exist_ok=True)
    report = base / 'fulluv_preflight.jsonl'
    completed = set()
    if args.resume and report.is_file():
        for line in report.read_text(encoding='utf-8').splitlines():
            if line.strip():
                row=json.loads(line)
                completed.add(row['model'])
    for model in models:
        if model in completed:
            print('SKIP already checked',model,flush=True)
            continue
        row={'format':FORMAT,'model':model,'case':CASE,'status':'NotRun','seeds':[]}
        try:
            if model in BLOCKED_LEGACY:
                row['status']='BlockedLegacyC5Parity';row['reason']='Legacy fast/authoritative C5 mismatch; no fit attempted'
            elif model in REUSED:
                row['status']='SpecialMixedScalar';row['reason']='Requires separate verified B/C implementation'
            else:
                entry,new,reused=bulk.inventory(registry,model)
                if reused:
                    raise ValueError('Unsupported reused extra quartics')
                hist=bulk.historical(CASE,model,entry)
                if not hist['uv_on']:
                    raise ValueError('Expected full-UV historical scenario')
                path=bulk.OUT / args.source_campaign / 'betas' / f'{bulk.safe(model)}.json'
                payload=bulk.read(path)
                # Important: compute actual UV trajectory and EFT1 for each seed.
                vec,layout,_=bulk.uv_base(hist['cfg'],entry)
                uv_beta,uv_path=bulk.uv_betas(model,entry,new,layout)
                q0=tuple(complex(args.extra_real) for _ in new)
                row['uv_beta_file']=str(uv_path.relative_to(bulk.ROOT))
                row['eft_beta_file']=str(path.relative_to(bulk.ROOT))
                row['status']='AllSeedsMassUnsafe'
                for label,x in candidate_vectors(model,args.frozen_run,hist):
                    record={'seed':label,'status':'NotRun'}
                    try:
                        cfg=bulk.update_yukawas(hist,x)
                        vec,layout,_=bulk.uv_base(cfg,entry)
                        uv,qf,nuv=bulk.run_uv(vec,layout,uv_beta,new,q0)
                        eft,masses,qs,neft=bulk.run_eft(uv,qf,payload,model,new,())
                        chk=bulk.testscreen(uv,eft,masses,qf,qs)
                        record.update(chk)
                        record.update({'uv_nfev':nuv,'eft1_nfev':neft,
                                       'status':'MassSafe' if chk['screen'] else 'MassUnsafe'})
                        if chk['screen']:
                            row['status']='AtLeastOneMassSafeSeed'
                    except Exception as e:
                        record['status']='NumericalFailure'
                        record['reason']=f'{type(e).__name__}: {e}'[:500]
                    row['seeds'].append(record)
                    print(f'{model} {label}: {record["status"]} masses_TeV={record.get("scalar_masses_muS_tev")} rhoF={record.get("rhoF")}',flush=True)
        except Exception as e:
            row['status']='InputOrRuntimeFailure';row['reason']=f'{type(e).__name__}: {e}'[:500]
        with report.open('a',encoding='utf-8') as f:
            f.write(json.dumps(row,allow_nan=False)+'\n')
        print('CHECKPOINT',model,row['status'],flush=True)
    print('COMPLETE',report,flush=True)
    return 0


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--source-campaign',default='overnight_six_v1_bulk48')
    p.add_argument('--models',default='')
    p.add_argument('--frozen-run',default='',help='Optional run-id with saved small-frozen fitted states')
    p.add_argument('--extra-real',type=float,default=-.01)
    p.add_argument('--resume',action='store_true')
    return probe(p.parse_args(argv))

if __name__=='__main__':
    raise SystemExit(main())
