"""Read-only A(-4) fixed-index vs scalar-only-diagnostic beta decomposition.

The earlier five-monomial diagnostic removed scalar-only differences from an
old RGBeta result. The permanent fresh-index fix regenerates *full* one-loop
betas, which can differ at nonzero gauge/Yukawa background. This script tests
that distinction at the same coupling points, without solving another ODE.

Run from project root:
  python scripts/InspectAminus4FixedFullBetaDifferences.py \
    --output output/a_minus4_next_batch/fresh_index_production/full_beta_diagnostic.json

No production source, beta payload, or previous numerical output is modified.
This does not certify the Matchete/RGBeta scalar dictionary or Wilson boundary.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.running import Aminus4EFT1Couplings as eft

FIXED = Path('output/a_minus4_next_batch/fresh_index_production/A_m4_fixed_17beta_compatible.json')
DIAG = Path('output/a_minus4_next_batch/quartic_five_monomial/scalar_only_reference_DIAGNOSTIC_beta.json')
OLD = Path('output/a_minus4_next_batch/quartic_five_monomial/native_global_projection_compatible.json')
FIXED_TRAJECTORY = Path('output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensor_basis.json')
DIAG_TRAJECTORY = Path('output/a_minus4_next_batch/conditional_quartic_ode/reference_diagnostic/A_m4_eft1_tensor_basis.json')
OLD_TRAJECTORY = Path('output/a_minus4_next_batch/conditional_quartic_ode/native/A_m4_eft1_tensor_basis.json')
QUARTIC_NAMES = tuple(sorted({'lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
                              'lambdaH1', 'lambdaH2', 'lambdaH2Adj',
                              'lambda12', 'lambdaT3'}))


def read(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f'Required previous result missing: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return obj


def evaluated(beta: dict, env: dict, name: str) -> complex:
    z = evaluate_inputform_expression(beta['report_betas'][name], env)
    value = np.asarray(z, dtype=complex)
    if value.ndim != 0 or not np.all(np.isfinite(value)):
        raise ValueError(f'{name} beta must be finite scalar')
    return complex(value.item())


def input_at(point: dict, *, gauge: bool, yukawa: bool) -> dict:
    p = eft._point(point)
    if not gauge:
        for k in ('gY', 'g2', 'g3'):
            p[k] = 0.0
    if not yukawa:
        for k in ('yu', 'yd', 'ye'):
            p[k] = np.zeros((3, 3), dtype=complex)
    return p


def decomposed_difference(fixed: dict, diagnostic: dict, point: dict,
                          name: str) -> dict[str, complex]:
    """Exact inclusion-exclusion at a fixed coupling background.

    scalar_only is evaluated with all gauge and SM Yukawas switched off;
    gauge_dependent includes pure-gauge and lambda*g^2 terms;
    yukawa_dependent includes the Yukawa terms; joint_interference is the
    remaining mixed gauge-Yukawa dependence (typically zero at one loop).
    """
    def difference(g: bool, y: bool) -> complex:
        env = input_at(point, gauge=g, yukawa=y)
        return evaluated(fixed, env, name) - evaluated(diagnostic, env, name)
    f11, f10, f01, f00 = difference(True, True), difference(True, False), \
                          difference(False, True), difference(False, False)
    return {
        'full': f11,
        'scalar_only': f00,
        'gauge_dependent': f10 - f00,
        'yukawa_dependent': f01 - f00,
        'joint_gauge_yukawa_interference': f11 - f10 - f01 + f00,
    }


def encode(z: complex) -> dict[str, float]:
    return {'real': float(z.real), 'imag': float(z.imag), 'absolute': float(abs(z))}


def point_coordinate(point: dict, name: str) -> float:
    p = eft._point(point)
    if name.endswith('.real') and name.startswith('lambdaT3'):
        return float(complex(p['lambdaT3']).real)
    if name.endswith('.imag') and name.startswith('lambdaT3'):
        return float(complex(p['lambdaT3']).imag)
    return float(p[name])


def run(root: Path) -> dict:
    fix, diag, old = (read(root / p) for p in (FIXED, DIAG, OLD))
    for p in (fix, diag, old):
        eft.validate_beta_payload(p)
    if set(fix['report_betas']) != set(diag['report_betas']):
        raise ValueError('Beta inventories differ')
    if diag['metadata'].get('Aminus4DoNotTreatAsProductionOrPhysical') is not True:
        raise ValueError('Previous reference does not carry nonphysical-diagnostic status')
    trajectories = [read(root / p) for p in (FIXED_TRAJECTORY, DIAG_TRAJECTORY, OLD_TRAJECTORY)]
    for t in trajectories:
        if t.get('model') != 'A_-4' or len(t.get('running_couplings', [])) != 7:
            raise ValueError('Expected seven-scale A(-4) trajectory')
    fixed_running, diag_running, old_running = trajectories
    for t in (diag_running, old_running):
        if not np.allclose(fixed_running['sampled_scales_gev'],t['sampled_scales_gev'],
                           rtol=1e-12,atol=1e-8):
            raise ValueError('Comparison trajectory scales mismatch')
    # The fixed and diagnostic formulas are evaluated at identical running
    # coupling values. This isolates beta-expression changes from ODE feedback.
    records = []
    peak = {kind: {'absolute': 0.0, 'name': None, 'mu_gev': None}
            for kind in ('full', 'scalar_only', 'gauge_dependent',
                         'yukawa_dependent', 'joint_gauge_yukawa_interference')}
    for idx in (0, 3, 6):
        point = diag_running['running_couplings'][idx]
        mu = float(diag_running['sampled_scales_gev'][idx])
        differences = {}
        for name in QUARTIC_NAMES:
            ds = decomposed_difference(fix, diag, point, name)
            differences[name] = {kind: encode(val) for kind,val in ds.items()}
            for kind,val in ds.items():
                if abs(val) > peak[kind]['absolute']:
                    peak[kind] = {'absolute': float(abs(val)),
                                  'name': name, 'mu_gev': mu}
        records.append({'mu_gev': mu, 'beta_difference_16pi2': differences})
    # Full nonquartic beta set was independently checked in the original
    # source-fix batch; we verify agreement with the diagnostic again at all
    # three present backgrounds, not merely that the name inventory matches.
    max_nonquartic = 0.0
    for row in records:
        k = int(np.argmin(np.abs(np.asarray(diag_running['sampled_scales_gev']) - row['mu_gev'])))
        env = input_at(diag_running['running_couplings'][k], gauge=True,yukawa=True)
        for name in set(eft.BETA_NAMES) - set(QUARTIC_NAMES):
            a = np.asarray(evaluate_inputform_expression(fix['report_betas'][name],env), complex)
            b = np.asarray(evaluate_inputform_expression(diag['report_betas'][name],env), complex)
            if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
                raise ValueError(f'Nonquartic beta malformed: {name}')
            max_nonquartic = max(max_nonquartic, float(np.max(np.abs(a-b))))
    if max_nonquartic > 2e-9:
        raise ValueError('Unexpected change to a nonquartic beta channel')
    if peak['scalar_only']['absolute'] > 2e-8:
        raise ValueError('Corrected native and five-monomial diagnostic scalar-only betas disagree')
    # Compare *completed* trajectories separately; they generally differ
    # from fixed-point beta differences due to nonlinear RG feedback.
    endpoints = {}
    for name in ('lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaS2Adj',
                 'lambdaH1', 'lambdaH2', 'lambdaH2Adj', 'lambda12',
                 'lambdaT3.real', 'lambdaT3.imag'):
        x = point_coordinate(fixed_running['running_couplings'][-1], name)
        y = point_coordinate(diag_running['running_couplings'][-1], name)
        z = point_coordinate(old_running['running_couplings'][-1], name)
        endpoints[name] = {
            'fixed': x, 'old_scalar_only_diagnostic': y, 'original_native': z,
            'fixed_minus_diagnostic': x-y,
            'fixed_minus_original_native': x-z,
            'diagnostic_minus_original_native': y-z,
        }
    return {
        'format': 'a-minus4-fixed-index-full-beta-gauge-yukawa-decomposition-v1',
        'model': 'A_-4',
        'status': 'READ_ONLY_FORMULA_COMPARISON_NOT_PHYSICAL_CERTIFICATION',
        'reference_points': 'previous diagnostic trajectory, identical background for both beta evaluations',
        'beta_normalization': '16pi2*dX/dln(mu)',
        'beta_differences': records,
        'maximum_beta_difference_by_origin': peak,
        'nonquartic_beta_max_absolute_difference': max_nonquartic,
        'endpoint_running_comparison': endpoints,
        'physical_matchete_quartic_dictionary_certified': False,
        'physical_Wilson_boundary_certified': False,
        'Wolfram_rerun_required': False,
        'new_ODE_rerun_required': False,
        'EFT_order': 'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars simultaneously -> SMEFT',
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output',type=Path,
                    default=Path('output/a_minus4_next_batch/fresh_index_production/full_beta_diagnostic.json'))
    args = ap.parse_args(argv)
    root = args.root.resolve()
    result = run(root)
    destination = args.output if args.output.is_absolute() else root / args.output
    if destination.exists():
        raise FileExistsError(f'Refusing to overwrite previous diagnostic: {destination}')
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print('Fixed vs scalar-only diagnostic beta origins (max |16pi² Δβ|):')
    for kind, record in result['maximum_beta_difference_by_origin'].items():
        print(f'  {kind}: {record["absolute"]:.8g} ({record["name"]})')
    print('Endpoint lambdaH2Adj fixed minus previous diagnostic:',
          result['endpoint_running_comparison']['lambdaH2Adj']['fixed_minus_diagnostic'])
    print('Matchete/Wilson physical certification: NO')
    print('Read-only report:',destination)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
