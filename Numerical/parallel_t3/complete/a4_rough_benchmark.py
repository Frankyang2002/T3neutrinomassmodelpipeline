"""A(-4) executable, APPROXIMATE two-threshold neutrino benchmark.

This is a one-channel calculation intended to produce useful rough numerical
results, not the source-complete physical matching calculation.

Conventions/assumptions:
 * three heavy fermions removed together at F=1e5 GeV;
 * both heavy scalar multiplets removed together at S=1e3 GeV;
 * use the finite ordered A(-4) direct F-hard + S-soft Laurent remainder;
 * apply C5[p,q] = A[p,q] + A[q,p] from RunMatching.wl;
 * unit signed GammaCC -> Weyl operator normalisation (adjustable);
 * identify Matchete/RGBeta mixing-coupling coordinate without an
   independently derived full quartic dictionary;
 * keep one-loop couplings at F in the direct channel and scalar masses at S;
 * leave other first-threshold Wilsons, counterterms, and finite one-loop
   scalar operators OUT of this estimate.

The separate EFT1 master leading-log C5 is NOT added: the exact direct
F-hard+S-soft logarithm already includes the shared F/S logarithmic channel
in the common-scale expression. Adding a separate RGE-generated LLHH term
without a consistent two-scale hard/soft split could double-count that log.

SMEFT running and Takagi/PMNS are the repository's real numerical solvers.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import decode_native, encode_native
from Numerical.parallel_t3.complete.a4_direct_hard_soft import formal_numeric_direct_remainder


class RoughBenchmarkError(ValueError):
    pass


@dataclass(frozen=True)
class RoughBoundary:
    ordered: np.ndarray
    c5_at_ms: np.ndarray
    operator_factor: float
    hard_soft: dict[str, Any]


def _matrix(value: Any, shape: tuple[int, ...] = (3, 3)) -> np.ndarray:
    z = np.asarray(decode_native(value), complex)
    if z.shape != shape or not np.isfinite(z).all():
        raise RoughBenchmarkError(f'Expected finite {shape} complex array')
    return z


def estimate_boundary(f_input: Mapping[str, Any], s_input: Mapping[str, Any],
                      *, operator_factor: float = 1.0) -> RoughBoundary:
    """Convert the native ordered finite direct channel to a tentative C5.

    `operator_factor` represents the UNKNOWN signed GammaCC/Weyl
    normalisation. Default 1 is an approximation, not a certification.
    No extra loop factor is applied: the direct source already includes
    1/(16 pi^2) from the hbar convention.
    """
    if isinstance(operator_factor, bool) or not np.isfinite(operator_factor) or operator_factor == 0:
        raise RoughBenchmarkError('operator_factor must be finite nonzero real')
    formal = formal_numeric_direct_remainder(dict(f_input), dict(s_input))
    ordered = _matrix(formal['ordered_formal_finite_remainder_GeVinv'])
    c5 = float(operator_factor) * (ordered + ordered.T)
    if not np.allclose(c5, c5.T, atol=1e-30, rtol=1e-12):
        raise RoughBenchmarkError('Candidate C5 must be complex symmetric')
    if not np.isfinite(c5).all():
        raise RoughBenchmarkError('Nonfinite provisional C5')
    return RoughBoundary(ordered, c5, float(operator_factor), formal)


def run_smeft_observables(boundary: RoughBoundary, s_input: Mapping[str, Any],
                          target: Mapping[str, Any], *, vev_gev: float = 246.22,
                          mu_low_gev: float = 100.) -> tuple[dict, dict[str, np.ndarray]]:
    """Use existing numerical SM + Weinberg RGE, Takagi and NuFIT adapters.

    The input C5 is APPROXIMATE. Solver success does not upgrade the physics
    status of the matching input.
    """
    if not (np.isfinite(vev_gev) and 0 < vev_gev < 1000):
        raise RoughBenchmarkError('Invalid Higgs vev')
    if not (np.isfinite(mu_low_gev) and 0 < mu_low_gev < 1e3):
        raise RoughBenchmarkError('Expected mu_low below scalar threshold')
    if s_input.get('model') != 'A_-4' or float(s_input.get('scale_GeV', -1)) != 1000:
        raise RoughBenchmarkError('Expected A(-4) scalar endpoint at 1000 GeV')
    if s_input.get('both_heavy_scalar_multiplets_to_integrate_out_together') is not True:
        raise RoughBenchmarkError('Scalar threshold does not integrate both multiplets')

    values = s_input['native_eft1_rgbeta_couplings_at_muS']
    required = ('gY', 'g2', 'g3', 'lambdaH', 'ye', 'yu', 'yd')
    if any(k not in values for k in required):
        raise RoughBenchmarkError('Missing required SM + Weinberg running data')

    def real_scalar(name: str) -> float:
        z = complex(np.asarray(decode_native(values[name])))
        if not np.isfinite(z) or abs(z.imag) > 1e-10:
            raise RoughBenchmarkError(f'{name} must be real')
        return float(z.real)

    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg,
    )
    from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables
    from Numerical.parallel_t3.complete.full_general_completion import score_nufit

    initial = SMWeinbergInitialConditions(
        gY=real_scalar('gY'), g2=real_scalar('g2'), g3=real_scalar('g3'),
        lambdaH=real_scalar('lambdaH'),
        ye=_matrix(values['ye']), yu=_matrix(values['yu']), yd=_matrix(values['yd']),
        K=boundary.c5_at_ms,
    )
    # K is O(1e-15) GeV^-1, far below the stock integrator's 1e-11
    # absolute tolerance. Use a tolerance small enough to resolve K itself.
    low = evolve_sm_weinberg(initial, 1000., float(mu_low_gev),
                             rtol=1e-8, atol=1e-24)
    if not low.solver_success:
        raise RoughBenchmarkError('SM + Weinberg numerical RGE failed: '+low.solver_message)
    mnu_gev = -0.5 * vev_gev**2 * np.asarray(low.K, complex)
    obs, matrices = neutrino_observables(mnu_gev, low.ye)
    score = None
    if obs.get('mixing_defined') is True:
        score = score_nufit(obs, target)

    summary = {
        'status': 'ApproximateNumericalBenchmark_NotFullPhysicalMatching',
        'model': 'A_-4',
        'mu_F_gev': 1e5, 'mu_S_gev': 1e3,
        'mu_low_gev': float(mu_low_gev), 'vev_gev': float(vev_gev),
        'operator_factor_assumed': boundary.operator_factor,
        'C5_MS_approx_GeVinv': encode_native(boundary.c5_at_ms),
        'C5_MS_approx_frobenius_GeVinv': float(np.linalg.norm(boundary.c5_at_ms)),
        'C5_low_approx_GeVinv': encode_native(np.asarray(low.K, complex)),
        'mnu_low_approx_eV': encode_native(mnu_gev * 1e9),
        'observables': obs,
        'nufit_gaussian_approximation': score,
        'SMEFT_one_loop_solver_success': True,
        'SMEFT_one_loop_nfev': int(low.nfev),
        'neutrino_observable_algorithm': 'repository numerical Takagi/PMNS routine',
        'full_physical_C5_verified': False,
        'independent_certification_performed': False,
        'full_F_hard_operator_basis_and_counterterms_included': False,
        'full_EFT1_Wilson_RG_consistently_included': False,
        'full_joint_scalar_finite_matching_included': False,
        'double_counted_separate_formal_RG_term_added': False,
        'approximation': 'source-derived finite direct hard/soft channel only; unverified overall operator normalisation',
    }
    arrays = {'C5_MS_approx_GeVinv': boundary.c5_at_ms,
              'C5_low_approx_GeVinv': np.asarray(low.K, complex),
              'mnu_low_approx_GeV': mnu_gev, 'ye_low': low.ye}
    arrays.update(matrices)
    return summary, arrays
