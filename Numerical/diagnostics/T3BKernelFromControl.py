"""Quantify T3-B one-loop off-diagonal Wilson mixing on an actual EFT1 trajectory.

This is a kernel diagnostic, not a full-flavour C12/C5 solution or threshold
matching. It uses ONLY the independently audited off-diagonal coefficients.

    16 pi^2 ds/dln(mu)     += 8 sqrt(3)/3 conj(lambdaT3) kappa
    16 pi^2 dkappa/dln(mu) += 4 sqrt(3)/3 lambdaT3 s

Second-order feedback kernel integrates both directions in their actual time
ordering. This does not include diagonal gamma, lepton Yukawa wavefunctions,
or the full flavour dependence of the Wilson coefficient.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

LOOP = 16.0 * np.pi**2


def _complex_pair(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError('Expected complex coupling as [real, imaginary].')
    z = complex(float(value[0]), float(value[1]))
    if not np.isfinite(z.real) or not np.isfinite(z.imag):
        raise ValueError('Nonfinite lambdaT3')
    return z


def pair(value):
    z = complex(value)
    return [float(z.real), float(z.imag)]


def evaluate(control: dict):
    if control.get('status') != 'FixedOrderControlComplete':
        raise ValueError('Input must be the successful T3-B fixed-order control.')
    if control.get('model') != 'T3-B alpha=-1':
        raise ValueError('Unexpected model')
    data = control.get('intermediate_running', {})
    samples = data.get('samples')
    if not isinstance(samples, list) or len(samples) < 2:
        raise ValueError('Expected at least two EFT1 running samples')
    # Original trajectory descends from MF to MS; preserve that orientation.
    mu = np.asarray([float(p['mu_gev']) for p in samples], dtype=float)
    if not np.all(np.isfinite(mu)) or np.any(mu <= 0) or not np.all(np.diff(mu) < 0):
        raise ValueError('Samples must be positive, finite, and ordered MF -> MS')
    t = np.log(mu)
    lam = np.asarray([_complex_pair(p['lambdaT3']) for p in samples])
    forward = 4 * np.sqrt(3) / 3 * lam / LOOP  # s -> kappa
    backward = 8 * np.sqrt(3) / 3 * np.conjugate(lam) / LOOP  # kappa -> s
    # Integrals are with decreasing t. This sign is physical for downward running.
    cumul_f = cumulative_trapezoid(forward, t, initial=0)
    cumul_b = cumulative_trapezoid(backward, t, initial=0)
    # Dyson expansion: s -> kappa -> s, and vice versa, with time ordering.
    cumul_fb = cumulative_trapezoid(backward * cumul_f, t, initial=0)
    cumul_bf = cumulative_trapezoid(forward * cumul_b, t, initial=0)
    norm = control.get('fixed_order_c5_at_ms', {}).get('frobenius_norm')
    if norm is None or not np.isfinite(float(norm)):
        raise ValueError('Missing fixed-order full-flavour C5 Frobenius norm')
    return {
        'status': 'MixingKernelComplete',
        'model': 'T3-B alpha=-1',
        'scope': 'off-diagonal one-generation Wilson kernel on full-flavour renormalisable background',
        'n_points': len(samples),
        'mu_f_gev': float(mu[0]),
        'mu_s_gev': float(mu[-1]),
        'integrated_s_to_kappa': pair(cumul_f[-1]),
        'integrated_kappa_to_s': pair(cumul_b[-1]),
        'second_order_s_to_kappa_to_s': pair(cumul_fb[-1]),
        'second_order_kappa_to_s_to_kappa': pair(cumul_bf[-1]),
        'fixed_order_c5_frobenius_norm': float(norm),
        'lambdaT3_modulus_min': float(np.min(np.abs(lam))),
        'lambdaT3_modulus_max': float(np.max(np.abs(lam))),
        'warning': 'Kernels are dimensionless transfer coefficients; these are NOT neutrino masses, full-flavour beta functions, or a corrected C5. No Wilson boundary or diagonal running included.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = evaluate(json.loads(args.control.read_text(encoding='utf-8')))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
