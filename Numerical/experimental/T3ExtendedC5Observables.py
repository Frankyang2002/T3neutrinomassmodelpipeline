"""Strict, source-labelled numerical C5 -> neutrino observables adapter.

Does NOT infer an extended C5 from RGBeta trajectories or a symbolic
Matchete matching expression. Inputs are independently supplied numerical
C5 matrices with an explicit source label. No physical certification is inferred
from the self-asserted metadata in the input file.

Usage:
  python Numerical/experimental/T3ExtendedC5Observables.py --self-test
  python Numerical/experimental/T3ExtendedC5Observables.py --input file.json --output file.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from physics.NeutrinoObservables import calculate_neutrino_observables
from physics.NeutrinoTrajectory import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix

ALLOWED = {"A_m2", "A_p0", "D_m2", "D_p0", "E_m2", "E_p0"}
KIND = 't3_external_extended_c5_numeric_v1'


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def complex_matrix(raw: Any, label: str) -> np.ndarray:
    need(isinstance(raw, dict) and set(raw) == {'real', 'imag'},
         f'{label}: real/imag fields required')
    re = np.asarray(raw['real'], dtype=float)
    im = np.asarray(raw['imag'], dtype=float)
    need(re.shape == im.shape == (3, 3), f'{label}: must be complex 3x3')
    result = re + 1j * im
    need(bool(np.isfinite(result).all()), f'{label}: non-finite matrix')
    return result


def encoded(a: np.ndarray) -> dict:
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def from_object(obj: dict) -> dict:
    need(isinstance(obj, dict) and obj.get('kind') == KIND, 'Incorrect input schema')
    model = obj.get('model')
    need(model in ALLOWED, 'Only six extended-Yukawa models supported')
    mu = float(obj.get('scale_gev', float('nan')))
    vev = float(obj.get('vev_gev', 246.22))
    need(math.isfinite(mu) and mu > 0 and math.isfinite(vev) and vev > 0,
         'Invalid positive scale or VEV')
    ordering = str(obj.get('ordering', 'NO')).upper()
    need(ordering in {'NO', 'IO'}, 'Explicit NO/IO ordering required')
    c5 = complex_matrix(obj.get('c5_gev_inverse'), 'c5_gev_inverse')
    asym=float(np.linalg.norm(c5-c5.T))
    scale=max(float(np.linalg.norm(c5)), 1e-300)
    need(asym/scale < 1e-10,
         f'C5 is not complex symmetric (relative Frobenius {asym/scale:.4e})')
    c5 = .5 * (c5 + c5.T)
    # Do not silently assume the charged-lepton mass basis. Require ye.
    ye = complex_matrix(obj.get('charged_lepton_yukawa'), 'charged_lepton_yukawa')
    m = neutrino_mass_from_c5(c5, vev_gev=vev)
    basis = charged_lepton_mass_basis_matrix(m, ye)
    obs = calculate_neutrino_observables(basis, ordering=ordering)
    need(obs.takagi_residual <= 1e-7, f'Takagi residual {obs.takagi_residual} exceeds tolerance')
    u = obs.pmns_abs
    s13 = float(u[0, 2]**2)
    den = 1-s13
    need(den > 1e-12, 'PMNS conventional mixing angle denominator vanishes')
    s12 = float(u[0, 1]**2 / den)
    s23 = float(u[1, 2]**2 / den)
    provenance = obj.get('provenance')
    need(isinstance(provenance, dict) and isinstance(provenance.get('description'), str)
         and len(provenance['description'].strip()) >= 10,
         'Explicit source provenance.description is required')
    return {
        'format': 't3-external-numeric-c5-neutrino-observables-v1',
        'status': 'NeutrinoObservablesFromSuppliedC5',
        'model': model, 'mu_gev': mu, 'vev_gev': vev,
        'ordering': obs.ordering, 'mass_matrix_gev': encoded(m),
        'mass_matrix_charged_lepton_basis_gev': encoded(basis),
        'masses_ev': obs.masses_ev.tolist(),
        'delta_m21_sq_ev2': obs.delta_m21_sq_ev2,
        'delta_m31_sq_ev2': obs.delta_m31_sq_ev2,
        'delta_m32_sq_ev2': obs.delta_m32_sq_ev2,
        'pmns_abs': u.tolist(),
        'sin2_theta12': s12, 'sin2_theta13': s13, 'sin2_theta23': s23,
        'takagi_residual': obs.takagi_residual,
        'source_provenance': provenance,
        'source_c5_user_supplied': True,
        'extended_matchete_rgbeta_dictionary_verified_by_this_script': False,
        'two_threshold_extended_matching_verified_by_this_script': False,
        'physical_prediction_certified_by_this_script': False,
        'note': ('This computes the exact repository C5->mnu->Takagi observables map '
                 'for a supplied symmetric 3x3 C5 and ye. The script has not '
                 'verified the UV coupling dictionary or the origin of C5.'),
    }


def save(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def self_test() -> None:
    masses = np.array([.009, .012, .052])  # eV, independent illustrative diagonal test
    vev = 246.22
    mgev = np.diag(masses) / 1e9
    c5 = -2 * mgev / (vev * vev)
    inp = {'kind': KIND, 'model': 'A_p0', 'scale_gev': 100., 'vev_gev': vev,
           'ordering': 'NO', 'c5_gev_inverse': encoded(c5),
           'charged_lepton_yukawa': encoded(np.diag([1e-6, .0006, .01]).astype(complex)),
           'provenance': {'description': 'Synthetic diagonal C5 exact-value test'}}
    out = from_object(inp)
    need(np.allclose(out['masses_ev'], masses, rtol=1e-8, atol=1e-10),
         'Wrong reconstructed neutrino masses')
    need(abs(out['sin2_theta13']) < 1e-8 and out['takagi_residual'] < 1e-8,
         'Unexpected diagonal mixing/Takagi residual')
    bad = json.loads(json.dumps(inp))
    bad['c5_gev_inverse']['real'][0][1] = 1e-11
    try:
        from_object(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Non-symmetric C5 accepted')
    bad = json.loads(json.dumps(inp))
    bad.pop('charged_lepton_yukawa')
    try:
        from_object(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Missing charged lepton basis accepted')
    bad = json.loads(json.dumps(inp))
    bad['provenance'] = {}
    try:
        from_object(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('Unattributed C5 accepted')
    print('PASS: source-labelled symmetric C5 -> mnu, charged-lepton rotation, Takagi observables, negative controls')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--self-test', action='store_true')
    g.add_argument('--input', type=Path)
    p.add_argument('--output', type=Path,
                   default=ROOT / 'Reports/output/workstream_a/extended_c5_native/external_c5_observables.json')
    a = p.parse_args()
    if a.self_test:
        self_test(); return 0
    try:
        raw = a.input.read_bytes()
        payload = from_object(json.loads(raw.decode('utf-8-sig')))
        payload['input_sha256'] = sha256(raw).hexdigest()
        save(a.output, payload)
        print('PASS: numerically evaluated externally provided C5; physical origin NOT certified')
        print('Report:', a.output)
        return 0
    except Exception as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 4


if __name__ == '__main__':
    raise SystemExit(main())
