from __future__ import annotations

"""Generate saved frozen-model trajectory JSONs for display.

For frozen UV benchmarks, the input Yukawas are displayed as scale-independent
curves across the UV interval by construction. The low-energy segment stores the
accepted neutrino observables as a saved endpoint proxy across the EFT interval.

This is sufficient for the Lambda-style report requested here without rerunning
any fit/search procedure.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scripts import RunT3GeneralScalar48Campaign as bulk

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



def log_grid(lo: float, hi: float, n: int) -> list[float]:
    return [10 ** x for x in np.linspace(math.log10(lo), math.log10(hi), n)]


def fro_stats(vec: np.ndarray) -> tuple[float, float, float]:
    y1 = vec[:9].reshape(3, 3)
    y2 = vec[9:].reshape(3, 3)
    return float(np.linalg.norm(y1)), float(np.linalg.norm(y2)), float(np.max(np.abs(vec)))


def effective_weinberg_proxy(pred: dict) -> float | None:
    dm31 = pred.get('delta_m3l_sq_ev2')
    if dm31 is None:
        return None
    try:
        return float(abs(float(dm31)) ** 0.5)
    except Exception:
        return None


def compute_eft_trajectory(model, case, source_campaign, extra_real, points):
    """Integrate the implemented EFT1 RGEs; never interpolate endpoint observables."""
    from scipy.integrate import solve_ivp
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state, unpack_intermediate_scalar_state
    from scripts import RunT3GeneralScalar48Campaign as b
    reg = b.read(b.REGISTRY)
    entry, new, reused = b.inventory(reg, model)
    if reused:
        raise ValueError(f'{model}: special mixed-scalar matching not supported')
    hist = b.historical(case, model, entry)
    if hist['uv_on']:
        raise ValueError('Expected frozen-UV scenario')
    vec, layout, unpack = b.uv_base(hist['cfg'], entry)
    uv = unpack(vec, layout, mu_gev=1e5)
    initial_state = build_intermediate_scalar_boundary(uv, matching_scale_gev=1e5)
    payload = b.read(b.OUT / source_campaign / 'betas' / (b.safe(model) + '.json'))
    betas, scalar_layout = b.validate_eft_betas(payload, initial_state, model, new, ())
    initial, _ = pack_intermediate_scalar_state(initial_state)
    q0 = [v for _ in new for v in (float(extra_real), 0.0)]
    z0 = np.r_[initial, q0]
    t = np.linspace(math.log(1e5), math.log(1e3), points)
    sol = solve_ivp(lambda logmu, zz: b.eft_rhs(logmu, zz, scalar_layout, betas, new),
                    (t[0], t[-1]), z0, t_eval=t, method='DOP853', rtol=1e-8, atol=1e-11)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(f'{model}: EFT1 evolution failed: {sol.message}')
    rows=[]
    for i, ti in enumerate(sol.t):
        state = unpack_intermediate_scalar_state(sol.y[:scalar_layout.size, i], scalar_layout,
                                                   mu_gev=float(math.exp(ti)))
        extras = {name: {'real':float(sol.y[scalar_layout.size+2*j,i]),
                         'imag':float(sol.y[scalar_layout.size+2*j+1,i])}
                  for j,name in enumerate(new)}
        rows.append({'mu_gev':float(math.exp(ti)),
                     'mS1Sq_gev2':float(state.mS1Sq), 'mS2Sq_gev2':float(state.mS2Sq),
                     'extra_quartics':extras})
    return rows


def run(args) -> int:
    out_dir = OUT / args.run_id / 'trajectories_v1'
    out_dir.mkdir(parents=True, exist_ok=True)
    run_ids = [x.strip() for x in args.run_ids.split(',') if x.strip()]
    for model, case, run_id, path, best in accepted_states(run_ids):
        if model in {'B_-1','C_-1','A_-4','B_-3','C_-3','E_-4'}:
            print(f'SKIP {model}: special/legacy matching', flush=True)
            continue
        rows = compute_eft_trajectory(model, case, args.source_campaign, args.extra_real, args.eft_points)
        # UV segment is frozen by definition: use true UV boundary scalar masses.
        reg=bulk.read(bulk.REGISTRY)
        entry,_,_=bulk.inventory(reg,model)
        hist=bulk.historical(case,model,entry)
        uv0=hist['cfg']['base_state']['ordinary']
        uv=[]
        for scale in log_grid(args.uv_start_gev,args.fermion_threshold_gev,args.uv_points):
            uv.append({'mu_gev':float(scale),
                       'mS1Sq_gev2':float(uv0['mS1Sq']),
                       'mS2Sq_gev2':float(uv0['mS2Sq']),
                       'extra_quartics':{}})
        doc={'model':model,'case':case,'run_id':run_id,'source_path':str(path),
             'note':'UV state frozen; EFT1 scalar RGEs actually integrated. No synthetic Weinberg running.',
             'settings':{'uv_start_gev':args.uv_start_gev,
                         'fermion_threshold_gev':args.fermion_threshold_gev,
                         'scalar_endpoint_gev':args.scalar_endpoint_gev},
             'uv_trajectory':uv,'eft1_trajectory':rows}
        result=out_dir / f'{bulk.safe(model)}__{case}.json'
        result.write_text(json.dumps(doc,indent=2),encoding='utf-8')
        print(f'RGE {model} {case} -> {result}',flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', default='rge_v1')
    ap.add_argument('--run-ids', default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    ap.add_argument('--source-campaign', default='overnight_six_v1_bulk48')
    ap.add_argument('--extra-real', type=float, default=-0.01)
    ap.add_argument('--uv-start-gev', type=float, default=1e7)
    ap.add_argument('--fermion-threshold-gev', type=float, default=1e5)
    ap.add_argument('--scalar-endpoint-gev', type=float, default=1e3)
    ap.add_argument('--uv-points', type=int, default=24)
    ap.add_argument('--eft-points', type=int, default=24)
    return run(ap.parse_args(argv))


if __name__ == '__main__':
    raise SystemExit(main())
