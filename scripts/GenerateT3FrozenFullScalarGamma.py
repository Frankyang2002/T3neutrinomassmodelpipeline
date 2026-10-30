"""T1/T2/T3 genuine saved response scans for frozen full-scalar candidates.

The Yukawa perturbation is made at the frozen 100-TeV boundary. The scalar EFT1
background is held fixed (heavy Yukawas have decoupled); neutrino matching and
SM+Weinberg running are RE-EVALUATED for each gamma. This is not a new fit.
Physics classifications and limitations inherit from RunT3FullScalarNeutrino.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
from scripts import T3Frozen64Sources as policy
from scripts import RunT3FullScalarNeutrino as nu
from scripts import RunT3FullScalarContinuous as full

ROOT = full.ROOT
BASE = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64'


def parse_grid(raw: str) -> tuple[float,...]:
    nums = tuple(float(x.strip()) for x in raw.split(',') if x.strip())
    policy.require(len(nums) >= 3 and len(nums) == len(set(nums)), 'invalid gamma grid')
    policy.require(all(math.isfinite(x) and x > 0 for x in nums), 'nonpositive gamma')
    policy.require(any(abs(x-1) < 1e-12 for x in nums), 'gamma grid requires gamma=1')
    return tuple(sorted(nums))


def perturb(uv, direction: str, gamma: float):
    policy.require(direction in policy.DIRECTIONS, 'unknown direction')
    a = np.array(uv.y1, dtype=complex, copy=True)
    b = np.array(uv.y2, dtype=complex, copy=True)
    policy.require(a.shape == (3,3) and b.shape == (3,3), 'Yukawa shape changed')
    if direction == 'T1':
        a[0, :] *= gamma
    elif direction == 'T2':
        b[1, :] *= gamma
    else:
        a[2, :] *= gamma
        b[2, :] *= gamma
    return replace(uv, y1=a, y2=b).validated()


def compute(cell: dict, grid: tuple[float,...], source_root: Path, path: Path) -> tuple[int, int]:
    model = cell['model']
    traj = policy.read_json(policy.validate_source(ROOT, cell['trajectory'],cell['trajectory_sha256']))
    central = policy.read_json(policy.validate_source(ROOT, cell['neutrino'],cell['neutrino_sha256']))
    policy.require(traj.get('boundary_mode') == 'FROZEN_UV' and traj.get('model') == model,
                   'not frozen full-general-scalar running')
    policy.require(traj.get('settings',{}).get('benchmark') == cell['benchmark'],
                   'trajectory benchmark mismatch')
    policy.require(central.get('model') == model and central.get('source_continuous_sha256') == cell['trajectory_sha256'],
                   'central prediction does not derive from trajectory')
    uv = nu.decode_uv(traj['uv_fermion_threshold_state'])
    eft = nu.decode_eft(traj['eft1_scalar_endpoint_state'])
    expected = central['oscillation_observables_normal_ordering']
    rows = []
    errors = 0
    for direction in policy.DIRECTIONS:
        for gamma in grid:
            row = {'model':model,'benchmark':cell['benchmark'], 'case':cell['case'],
                   'yukawa_mode':cell['yukawa_mode'], 'direction':direction,
                   'gamma':gamma, 'status':'BLOCKED',
                   'matching_status':central['status'],
                   'trajectory_sha256':cell['trajectory_sha256'],
                   'neutrino_sha256':cell['neutrino_sha256'],
                   'error':'', **{k:None for k in policy.FIELDS}}
            try:
                changed = perturb(uv,direction,gamma)
                if model in nu.MIXED:
                    c5,_,_,_ = nu.matched_bc_c5(model,changed,source_root)
                else:
                    c5,_,_ = nu.matched_nonmixed_c5(model,changed,eft)
                actual = nu.low_energy(model,c5,eft)['oscillation_observables_normal_ordering']
                for field in policy.FIELDS:
                    value = float(actual[field])
                    policy.require(math.isfinite(value), 'nonfinite response: '+field)
                    row[field] = value
                    if math.isclose(gamma,1.0,rel_tol=0,abs_tol=1e-12):
                        policy.require(math.isclose(value,float(expected[field]),rel_tol=1e-6,abs_tol=1e-10),
                            'gamma=1 not equal to saved central prediction')
                row['status'] = 'NeutrinoEvaluated'
            except (ValueError, KeyError, RuntimeError, TypeError, OSError, ArithmeticError) as error:
                errors += 1
                row['error'] = f'{type(error).__name__}: {error}'[:450]
            rows.append(row)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), errors


def run(args) -> int:
    base = BASE / args.run_id
    manifest = policy.validate_manifest(policy.read_json(base / 'manifest.json'))
    gamma = parse_grid(args.gamma_grid)
    targets = [c for c in manifest['cells'] if
        (args.benchmark == 'all' or c['benchmark'] == args.benchmark) and
        (args.models == 'all' or c['model'] in args.models.split(','))]
    errors = 0
    built = []
    out = base / args.scan_id
    policy.require(bool(policy.RUN_ID.fullmatch(args.scan_id)), 'unsafe scan ID')
    for c in targets:
        if c['status'] != 'NEUTRINO_COMPLETE':
            print('SKIP missing central prediction:',c['benchmark'],c['model'])
            continue
        model, benchmark = c['model'],c['benchmark']
        source_folder = base / benchmark / full.model_tag(model)
        path = out / benchmark / (full.model_tag(model)+'.csv')
        meta_path = path.with_suffix('.json')
        expected = {'model':model,'benchmark':benchmark,'gamma':gamma,
                    'trajectory_sha256':c['trajectory_sha256'],
                    'neutrino_sha256':c['neutrino_sha256'],
                    'generator_sha256':policy.digest(Path(__file__))}
        if path.exists() or meta_path.exists():
            policy.require(args.resume and path.is_file() and meta_path.is_file(),
                           'partial/old scan; choose new --scan-id: '+str(path))
            saved = policy.read_json(meta_path)
            policy.require(saved.get('source') == dict(expected,gamma=list(gamma)) and
                           saved.get('csv_sha256') == policy.digest(path),
                           'source or grid changed; choose new --scan-id')
            built.append(saved)
            continue
        total, failures = compute(c,gamma,source_folder,path)
        errors += failures
        summary = {'source':dict(expected,gamma=list(gamma)), 'csv_sha256':policy.digest(path),
                   'rows':total,'failures':failures,
                   'physics_status':'candidate full-scalar; self-bilinear neutrino controls are restricted'}
        full.write_json_atomic(meta_path,summary)
        built.append(summary)
        print(f'{benchmark}/{model}: {total} scan points, {failures} failures',flush=True)
    out.mkdir(parents=True, exist_ok=True)
    full.write_json_atomic(out/'summary.json',{
        'scan_id':args.scan_id,'run_id':args.run_id,'gamma_grid':list(gamma),
        'models_scanned':len(built), 'requested':len(targets),
        'failed_points_current_run':errors})
    return 2 if errors else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='frozen_general64_v1')
    p.add_argument('--scan-id', default='gamma_scans_v1')
    p.add_argument('--benchmark', choices=('all',*policy.BENCHMARKS),default='all')
    p.add_argument('--models',default='all')
    p.add_argument('--gamma-grid',default='0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5')
    p.add_argument('--resume',action='store_true')
    return run(p.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
