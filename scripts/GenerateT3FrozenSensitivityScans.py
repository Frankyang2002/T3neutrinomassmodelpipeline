from __future__ import annotations

"""Generate real frozen-model T1/T2/T3 gamma scan points from accepted fits.

Implemented directions:
- T1: scale the electron-flavour row of Y1 by gamma.
- T2: scale the muon-flavour row of Y2 by gamma.
- T3: scale the tau-flavour row of both Yukawas by gamma.

This uses the existing frozen evaluator. No search is rerun.
"""

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scripts import RunT3GeneralScalar48Campaign as bulk
from Numerical.orchestration import RunSmallLambdaRefit as refit
from Numerical.fitting.OscillationFit import OscillationFitTarget

OUT = Path('output/studies/general_scalar_extended/frozen_display')
CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')


def accepted_states(run_ids: list[str]):
    base = bulk.OUT
    for case in CASES:
        for model in bulk.read(bulk.REGISTRY)["models"]:
            for run_id in run_ids:
                matches = [base / run_id / case / bulk.safe(model) / 'state.json']
                matches += sorted(base.glob(f'{run_id}_*/{case}/{bulk.safe(model)}/state.json'))
                for path in matches:
                    if not path.is_file():
                        continue
                    doc = bulk.read(path)
                    best = doc.get('best') or {}
                    if doc.get('status') not in {'Accepted', 'CandidateConstrainedAccepted', 'ExistingConditionalAccepted'}:
                        continue
                    params = best.get('parameters') or {}
                    if not all(name in params for name in bulk.NAMES):
                        continue
                    yield model, case, path.parents[2].name, path, best
                    break
                else:
                    continue
                break
            continue



def build_evaluator(model: str, case: str, source_campaign: str, extra_real: float):
    reg = bulk.read(bulk.REGISTRY)
    entry, new, reused = bulk.inventory(reg, model)
    if reused:
        raise ValueError(f'{model}: mixed-scalar evaluator is not supported by this frozen scan script')
    hist = bulk.historical(case, model, entry)
    if hist['uv_on']:
        raise ValueError(f'{model} {case}: UV-running case not allowed here')
    beta_path = bulk.OUT / source_campaign / 'betas' / (bulk.safe(model) + '.json')
    payload = bulk.read(beta_path)
    kernel, _matching_path, _matching_kind = bulk.matching(model)
    target = OscillationFitTarget.from_json(refit.TARGET)
    return bulk.Evaluator(case, model, entry, new, hist, payload, kernel, target, extra_real)


def parse_gamma_grid(value: str) -> list[float]:
    vals = []
    for part in value.split(','):
        part = part.strip()
        if not part:
            continue
        vals.append(float(part))
    if not vals or not all(math.isfinite(v) and v > 0 for v in vals):
        raise ValueError('Invalid gamma grid')
    vals = sorted(set(vals))
    if 1.0 not in vals:
        vals.append(1.0)
        vals.sort()
    return vals


def apply_direction(vec: np.ndarray, direction: str, gamma: float) -> np.ndarray:
    out = vec.copy()
    if direction == 'T1':
        out[:3] *= gamma
    elif direction == 'T2':
        out[12:15] *= gamma
    elif direction == 'T3':
        out[6:9] *= gamma
        out[15:18] *= gamma
    else:
        raise ValueError(direction)
    return out


def fro_stats(vec: np.ndarray) -> tuple[float, float, float]:
    y1 = vec[:9].reshape(3, 3)
    y2 = vec[9:].reshape(3, 3)
    return float(np.linalg.norm(y1)), float(np.linalg.norm(y2)), float(np.max(np.abs(vec)))


def run(args) -> int:
    gamma_grid = parse_gamma_grid(args.gamma_grid)
    out_dir = OUT / args.run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / 'frozen_gamma_points.csv'
    run_ids = [x.strip() for x in args.run_ids.split(',') if x.strip()]
    fields = ['model', 'case', 'run_id', 'source_path', 'direction', 'gamma', 'status', 'chi2', 'y1_fro', 'y2_fro', 'max_abs_y',
              'mass_hi_ev', 'mass_mid_ev', 'mass_lo_ev', 'delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23']
    with out_csv.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for model, case, run_id, path, best in accepted_states(run_ids):
            evaluator = build_evaluator(model, case, args.source_campaign, args.extra_real)
            base = np.array([float(best['parameters'][n]) for n in bulk.NAMES], dtype=float)
            print(f'SCAN {model} {case}', flush=True)
            for direction in ('T1', 'T2', 'T3'):
                for gamma in gamma_grid:
                    vec = apply_direction(base, direction, gamma)
                    row = evaluator.evaluate(vec, f'{direction}_gamma_{gamma:g}')
                    pred = row.get('prediction') or {}
                    y1f, y2f, ymax = fro_stats(vec)
                    writer.writerow({
                        'model': model, 'case': case, 'run_id': run_id, 'source_path': str(path), 'direction': direction, 'gamma': gamma,
                        'status': row.get('status') or '', 'chi2': row.get('chi2'), 'y1_fro': y1f, 'y2_fro': y2f, 'max_abs_y': ymax,
                        'mass_hi_ev': pred.get('mass_hi_ev') or pred.get('largest_neutrino_mass_ev'),
                        'mass_mid_ev': pred.get('mass_mid_ev') or pred.get('middle_neutrino_mass_ev'),
                        'mass_lo_ev': pred.get('mass_lo_ev') or pred.get('smallest_neutrino_mass_ev'),
                        'delta_m21_sq_ev2': pred.get('delta_m21_sq_ev2'), 'delta_m3l_sq_ev2': pred.get('delta_m3l_sq_ev2'),
                        'sin2_theta12': pred.get('sin2_theta12'), 'sin2_theta13': pred.get('sin2_theta13'), 'sin2_theta23': pred.get('sin2_theta23'),
                    })
    print(f'Wrote {out_csv}')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', default='sensitivity_v1')
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--source-campaign', default='overnight_six_v1_bulk48')
    ap.add_argument('--gamma-grid', default='0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5')
    ap.add_argument('--extra-real', type=float, default=-0.01)
    return run(ap.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
