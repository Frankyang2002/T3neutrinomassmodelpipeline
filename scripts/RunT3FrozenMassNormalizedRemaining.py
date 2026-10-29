"""Sequential, checkpointed mass-normalized refinement of unresolved frozen T3 cases.

Does not modify matching, RGEs, existing fits, or GitHub. Runs only frozen cases.
Stops on first nonaccepted fit, error or manual interruption; --resume continues it.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts import SearchT3FrozenPersistent as search

CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
ACCEPTED = {'CandidateConstrainedAccepted', 'ExistingConditionalAccepted'}
# Authoritative C5 parity for these extreme-hypercharge kernels was not established.
PARITY_BLOCKED = {'A_-4', 'B_-3', 'C_-3', 'E_-4'}


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def accepted_chi2(doc, target):
    if not isinstance(doc, dict):
        return None
    chi = doc.get('chi2')
    if chi is None and isinstance(doc.get('best'), dict):
        chi = doc['best'].get('chi2')
    try:
        value = float(chi)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and 0 <= value <= target else None


def existing_success(root, case, model, target, previous_runs):
    """Return evidence of successful fit, without trusting an unchecked campaign status."""
    directory = root / case / bulk.safe(model)
    for run in previous_runs:
        state_path = root / run / case / bulk.safe(model) / 'state.json'
        if state_path.is_file():
            state = load_json(state_path)
            if state.get('model') == model and state.get('case') == case and state.get('status') == 'Accepted':
                chi = accepted_chi2(state.get('best'), target)
                if chi is not None:
                    return {'source': str(state_path), 'chi2': chi}
    # Historical campaign outcomes are separate from the new runner.
    if directory.joinpath('fit.json').is_file():
        data = load_json(directory / 'fit.json')
        if data.get('model') == model and data.get('case') == case and data.get('status') in ACCEPTED:
            chi = accepted_chi2(data.get('best_mass_safe'), target)
            if chi is not None:
                return {'source': str(directory / 'fit.json'), 'chi2': chi}
    return None


def plan(rows, cases, models, target, skipped_success, include_parity_blocked=False):
    """Return (pending, skipped). Treat absent or malformed rows as errors."""
    lookup = {(r['model'], r['study']): r for r in rows}
    pending, skipped = [], []
    for model in models:
        for case in cases:
            key = (model, case)
            if key not in lookup:
                raise ValueError(f'Missing original campaign row {key}')
            row = lookup[key]
            if key in skipped_success:
                skipped.append((key, 'previously accepted: ' + skipped_success[key]['source']))
            elif row['status'] in ACCEPTED and accepted_chi2(row, target) is not None:
                skipped.append((key, 'already accepted in baseline campaign'))
            elif model in bulk.REUSED:
                skipped.append((key, 'mixed-scalar special evaluator required'))
            elif model in PARITY_BLOCKED and not include_parity_blocked:
                skipped.append((key, 'legacy C5 parity unresolved; blocked'))
            else:
                pending.append(key)
    return pending, skipped


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline-campaign', default='overnight_six_v1_bulk48', help='Campaign containing campaign.json and betas/')
    ap.add_argument('--run-prefix', default='frozen_massnorm_remaining_v1')
    ap.add_argument('--previous-run', action='append', default=[], help='Runs containing accepted state.json files, repeatable')
    ap.add_argument('--from-run', action='append', default=[], help='Seed sources imported by the underlying optimizer, repeatable')
    ap.add_argument('--cases', default=','.join(CASES))
    ap.add_argument('--models', default='', help='Comma-separated subset, otherwise all 16 models')
    ap.add_argument('--points', type=int, default=256)
    ap.add_argument('--seeds', type=int, default=3)
    ap.add_argument('--nfev', type=int, default=300)
    ap.add_argument('--radii', default='0.02,0.05,0.10')
    ap.add_argument('--min-seed-distance', type=float, default=.08)
    ap.add_argument('--diff-step', type=float, default=1e-6)
    ap.add_argument('--mass-scaling-rtol', type=float, default=.03)
    ap.add_argument('--extra-real', type=float, default=-.01)
    ap.add_argument('--seed', type=int, default=20261030)
    ap.add_argument('--target', type=float, default=5.86)
    ap.add_argument('--max-cycles', type=int, default=0, help='0: run each model until accepted or interrupted')
    ap.add_argument('--resume', action='store_true', help='Resume checkpointed model currently in progress')
    ap.add_argument('--dry-run', action='store_true')
    args=ap.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,55}', args.run_prefix):
        ap.error('Invalid --run-prefix')
    cases=tuple(x.strip() for x in args.cases.split(','))
    if not cases or len(set(cases)) != len(cases) or any(c not in CASES for c in cases):
        ap.error('Only the two distinct frozen cases are supported')
    reg=bulk.read(bulk.REGISTRY)
    all_models=tuple(reg['models'])
    models=tuple(x.strip() for x in args.models.split(',')) if args.models else all_models
    if len(set(models)) != len(models) or any(m not in reg['models'] for m in models):
        ap.error('Invalid model list')
    baseline=bulk.OUT / args.baseline_campaign
    report=load_json(baseline/'campaign.json')
    if report.get('model_count') != 16 or report.get('study_count') != 3:
        raise ValueError('Not a 16 x 3 baseline campaign report')
    previous_runs=list(dict.fromkeys(args.previous_run + [args.run_prefix]))
    done={}
    for model in models:
        for case in cases:
            proof=existing_success(bulk.OUT,case,model,args.target,previous_runs)
            if proof is not None:
                done[(model,case)]=proof
    pending,skipped=plan(report['rows'],cases,models,args.target,done)
    print(f'PLAN: {len(pending)} pending; {len(skipped)} skipped',flush=True)
    for key, reason in skipped:
        print(f'SKIP {key[0]} {key[1]}: {reason}',flush=True)
    for model,case in pending:
        run_id=f'{args.run_prefix}_{bulk.safe(model)}_{"s" if case == CASES[0] else "l"}'
        state_path=bulk.OUT/run_id/case/bulk.safe(model)/'state.json'
        if state_path.exists() and not args.resume:
            raise FileExistsError(f'Checkpoint exists: {state_path}. Rerun with --resume.')
        print(f'PENDING {model} {case} -> {run_id}',flush=True)
    if args.dry_run:
        return 0
    for index,(model,case) in enumerate(pending,1):
        run_id=f'{args.run_prefix}_{bulk.safe(model)}_{"s" if case == CASES[0] else "l"}'
        state_path=bulk.OUT/run_id/case/bulk.safe(model)/'state.json'
        opts=['--run-id',run_id,'--model',model,'--case',case,
              '--source-campaign',args.baseline_campaign,
              '--points',str(args.points),'--seeds',str(args.seeds),
              '--nfev',str(args.nfev),'--radii',args.radii,
              '--min-seed-distance',str(args.min_seed_distance),
              '--diff-step',str(args.diff_step),
              '--mass-scaling-rtol',str(args.mass_scaling_rtol),
              '--extra-real',str(args.extra_real),
              '--seed',str(args.seed),'--target',str(args.target),
              '--max-cycles',str(args.max_cycles)]
        for source in dict.fromkeys(args.from_run + args.previous_run):
            source_folder=bulk.OUT/source/case/bulk.safe(model)
            # The underlying load_seeds() raises for absent per-model folders.
            # Previous A(+0) checkpoints should not prevent searches for D/E/etc.
            if source_folder.is_dir():
                opts.extend(['--from-run',source])
        if state_path.exists():
            opts.append('--resume')
        print(f'RUN {index}/{len(pending)} {model} {case}',flush=True)
        search.main(opts)
        if not state_path.is_file():
            raise RuntimeError(f'Missing checkpoint for {model} {case}')
        state=load_json(state_path)
        chi=accepted_chi2(state.get('best'),args.target)
        if state.get('status') != 'Accepted' or chi is None:
            print(f'HALTED {model} {case}: not accepted. Resume same command later.',flush=True)
            return 2
        print(f'ACCEPTED {model} {case} chi2={chi:.12g}',flush=True)
    print('COMPLETE: all selected supported frozen cases have accepted fits.',flush=True)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
