"""Check absolute T3-B direct-Weinberg normalisation with no fitted factor.

Only validates production's fixed-order K_running against independently
constructed heavy-eigenstate-resolved LL tensor. Does not calculate physical
resummed corrections or claim complete three-generation RGEs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from Numerical.diagnostics.T3BHeavyLogAlignment import (
    flavor_kernel_by_heavy,
    matrix_to_json,
    run_production,
)


def evaluate(direct, y1, y2, heavy_masses, ms, lambda_t3, *, rtol=1e-9):
    """Apply symbolic gamma_5,12 and production Majorana conversion factor.

    One-generation beta: gamma_5,12 = (4 sqrt(3)/3) lambdaT3.
    Production's physical-Majorana direct coefficient is twice the ordered
    tensor. Its hbar is 1/(16*pi^2). Each heavy eigenstate uses ln(MS/Mr).
    """
    masses = np.asarray(heavy_masses, dtype=float)
    if masses.ndim != 1 or np.any(masses <= ms) or ms <= 0:
        raise ValueError('Expected heavy eigenmasses strictly larger than MS')
    tensors = flavor_kernel_by_heavy(y1, y2, masses)
    logs = np.log(ms / masses)
    resolved = np.einsum('r,rij->ij', logs, tensors)
    coefficient = (2.0 * (4.0 * np.sqrt(3.0) / 3.0)
                   * complex(lambda_t3) / (16.0 * np.pi**2))
    predicted = coefficient * resolved
    actual = np.asarray(direct, complex)
    if actual.shape != (3, 3) or not np.all(np.isfinite(actual)):
        raise ValueError('Production direct C5 must be finite 3x3')
    norm = float(np.linalg.norm(actual, 'fro'))
    if norm == 0:
        raise ValueError('Production direct C5 is zero; relative check undefined')
    discrepancy = float(np.linalg.norm(actual - predicted, 'fro') / norm)
    return {
        'status': 'UnfittedNormalizationPassed' if discrepancy <= rtol else 'UnfittedNormalizationMismatch',
        'relative_frobenius_residual': discrepancy,
        'acceptance_rtol': rtol,
        'symbolic_gamma_5_12': '(4*sqrt(3)/3)*lambdaT3',
        'physical_majorana_conversion_factor': 2,
        'loop_factor': '1/(16*pi**2)',
        'heavy_resolved_logs': logs.tolist(),
        'lambda_t3_at_ms': [complex(lambda_t3).real, complex(lambda_t3).imag],
        'unfitted_complex_coefficient': [coefficient.real, coefficient.imag],
        'production_direct_norm': norm,
        'production_direct_c5': matrix_to_json(actual),
        'predicted_direct_c5': matrix_to_json(predicted),
        'scope': 'unfitted absolute direct one-loop leading-log normalisation; not resummed running',
        'neutrino_mass_prediction_available': False,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('config', 'uv-rgbeta', 'eft1-rgbeta', 'final-c5', 'output'):
        ap.add_argument('--'+key, required=True, type=Path)
    ap.add_argument('--rtol', type=float, default=1e-9)
    args = ap.parse_args()
    if not 0 < args.rtol < 1:
        raise ValueError('--rtol must lie in (0,1)')
    direct, inputs, _ = run_production(args.config, args.uv_rgbeta,
                                       args.eft1_rgbeta, args.final_c5)
    outcome = evaluate(direct, inputs.y1, inputs.y2,
                       inputs.heavy_masses_gev,
                       inputs.grouped_scalar_matching_scale_gev,
                       inputs.lambda_t3, rtol=args.rtol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(outcome, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: outcome[key] for key in (
        'status','relative_frobenius_residual','acceptance_rtol',
        'unfitted_complex_coefficient','production_direct_norm',
        'neutrino_mass_prediction_available')}, indent=2))
    return 0 if outcome['status'] == 'UnfittedNormalizationPassed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
