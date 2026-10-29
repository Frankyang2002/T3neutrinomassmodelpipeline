"""Sequential frozen-UV 18D local-refinement radius sweep.

Runs the existing RefineT3Frozen18Stable module with increasing trust windows.
Each stage imports the previous stage's best fit and all original seed sources.
Numerical status 2 means a valid but poor fit; other nonzero statuses abort.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

from scripts import RunT3GeneralScalar48Campaign as bulk

RADII = (0.005, 0.01, 0.02)


def stage_id(run_id, index):
    return f'{run_id}_r{index:02d}'


def saved_state(run_id, case, model):
    path = bulk.OUT / run_id / case / bulk.safe(model) / 'state.json'
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding='utf-8-sig'))


def checked_completed(state, model, case):
    if not isinstance(state, dict) or state.get('model') != model or state.get('case') != case:
        raise ValueError('Saved stage identity mismatch')
    if state.get('status') not in ('CandidateConstrainedAccepted', 'CandidateConstrainedPoorFit', 'NoMassSafeNeutrinoFit'):
        raise ValueError('Saved stage incomplete; cannot auto-resume an interrupted stage')
    return state


def run(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,49}', args.run_id):
        raise ValueError('Invalid run ID')
    if args.case not in ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5'):
        raise ValueError('Frozen scenarios only')
    if not 0 < args.target <= 100 or args.seeds < 1 or args.nfev < 20:
        raise ValueError('Invalid refinement settings')
    if not args.from_run:
        raise ValueError('Supply at least one existing --from-run')
    base = bulk.OUT / args.run_id
    base.mkdir(parents=True, exist_ok=True)
    summary_path = base / 'radius_sweep.json'
    previous = None
    history = []
    for i, radius in enumerate(RADII, 1):
        sid = stage_id(args.run_id, i)
        state = saved_state(sid, args.case, args.model)
        if state is None:
            imports = list(args.from_run)
            if previous:
                imports.insert(0, previous)
            cmd = [sys.executable, '-u', '-X', 'utf8', '-m', 'scripts.RefineT3Frozen18Stable',
                   '--run-id', sid, '--model', args.model, '--case', args.case,
                   '--source-campaign', args.source_campaign,
                   '--seeds', str(args.seeds), '--nfev', str(args.nfev),
                   '--radius', str(radius), '--diff-step', str(args.diff_step),
                   '--extra-real', str(args.extra_real), '--target', str(args.target)]
            for name in imports:
                cmd += ['--from-run', name]
            print(f'STAGE {i}/{len(RADII)} radius={radius} run-id={sid}', flush=True)
            result = subprocess.run(cmd, cwd=bulk.ROOT, check=False)
            if result.returncode not in (0, 2):
                raise RuntimeError(f'Refinement subprocess failed at stage {i}: exit={result.returncode}')
            state = saved_state(sid, args.case, args.model)
            if state is None:
                raise RuntimeError(f'Refinement stage {i} exited but wrote no completed state')
        elif not args.resume:
            raise FileExistsError(f'Stage {sid} already exists. Use --resume.')
        checked_completed(state, args.model, args.case)
        best = state.get('best') or {}
        chi2 = best.get('chi2')
        if chi2 is not None and not math.isfinite(float(chi2)):
            raise ValueError('Nonfinite saved chi2')
        history.append({'run_id': sid, 'radius': radius, 'status': state['status'], 'chi2': chi2})
        bulk.OUT.joinpath(args.run_id).mkdir(parents=True, exist_ok=True)
        temporary = summary_path.with_suffix('.tmp')
        temporary.write_text(json.dumps({'model': args.model, 'case': args.case,
            'target_chi2': args.target, 'stages': history}, indent=2) + '\n', encoding='utf-8')
        temporary.replace(summary_path)
        print(f'RESULT radius={radius} status={state["status"]} chi2={chi2}', flush=True)
        previous = sid
        if chi2 is not None and float(chi2) <= args.target:
            print('ACCEPTED: stopping radius sweep', flush=True)
            return 0
    print('UNRESOLVED: all radius stages completed; best remains above target', flush=True)
    return 2


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=True)
    p.add_argument('--model', default='A_+0')
    p.add_argument('--case', default='small_uv_frozen_5e-9')
    p.add_argument('--source-campaign', default='overnight_six_v1_bulk48')
    p.add_argument('--from-run', action='append', default=[])
    p.add_argument('--seeds', type=int, default=3)
    p.add_argument('--nfev', type=int, default=600)
    p.add_argument('--diff-step', type=float, default=1e-6)
    p.add_argument('--extra-real', type=float, default=-0.01)
    p.add_argument('--target', type=float, default=5.86)
    p.add_argument('--resume', action='store_true')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
