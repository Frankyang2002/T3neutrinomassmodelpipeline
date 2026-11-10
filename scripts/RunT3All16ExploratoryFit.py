#!/usr/bin/env python3
"""Exploratory 16-model real-Yukawa fit, with a corrected B/C doublet portal.

For doublet S1 and S2 the source lambdaT3 term is H H S1 S2^dagger.
The old CG product treats the conjugate SU(2) leg as a fundamental and
produces a tensor orthogonal to the Weinberg operator for B(-3).
This wrapper changes only that portal component expansion, and uses the
existing consistent L2 projection for rough, explicitly conditional fits.
It does not rewrite the repository physics modules or archived outputs.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.rough_inputs import builder as builder_module
from Numerical.parallel_t3.threshold import preflight as preflight_module
from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
from Numerical.parallel_t3.threshold.matching import project_master
from Numerical.parallel_t3.wilson import workflow, rough_stream2
from scripts import RunT3All16RoughYukawaFit as runner_module
from Numerical.parallel_t3.complete import all16_rough_yukawa_fit as fit_module

ORIGINAL_QUARTIC_UNITS = builder_module.quartic_unit_tensors
ORIGINAL_PREFLIGHT = preflight_module.preflight_pair
ORIGINAL_ROUGH_CASE = numerics_module.run_rough_case



def corrected_quartic_unit_tensors(model: str):
    """Correct only the H H S1 S2* portal for the B/C doublet models.

    SU(2) contraction for H_i H_j S1_k (S2*)_l:
        (epsilon[i,k] delta[j,l] + epsilon[j,k] delta[i,l]) / 2.
    The 2nd S leg is anti-fundamental, NOT a second fundamental. The
    RGBeta lambdaT3 trajectory is left unchanged; only its real tensor map
    is corrected. All other quartic channels retain existing conventions.
    """
    units = ORIGINAL_QUARTIC_UNITS(model)
    if model.split('_', 1)[0] not in ('B', 'C'):
        return units
    from Numerical.parallel_t3.wilson.basis import build_basis
    basis = build_basis(model)
    if basis.dims[:2] != (2, 2):
        raise RuntimeError('B/C scalar dimensions changed: expected two doublets')
    ns = basis.ns
    epsilon = np.array([[0., 1.], [-1., 0.]])
    portal = np.zeros((ns, ns, ns, ns), complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for ell in range(2):
                    coefficient = (epsilon[i, k] * (j == ell)
                                   + epsilon[j, k] * (i == ell)) / 2.
                    if coefficient:
                        portal += coefficient * builder_module._monomial(
                            builder_module._u(ns, 0, i),
                            builder_module._u(ns, 0, j),
                            builder_module._u(ns, 4, k),
                            builder_module._u(ns, 8, ell, True))
    if not np.isfinite(portal).all() or not np.allclose(
            portal, portal.transpose(1, 0, 2, 3), rtol=1e-12, atol=1e-12):
        raise RuntimeError('Corrected portal tensor is invalid')
    if not np.allclose(portal, portal.swapaxes(0, 2)):
        raise RuntimeError('Corrected portal tensor is not fully symmetric')
    units['lambdaT3'] = portal
    return units


def consistent_project_c5(full, basis):
    """Orthogonal Weinberg-template projection of full LLHH; absolute residual.

    Must agree with the downstream threshold matcher. Does NOT discard the
    non-Weinberg remainder from the evolved full Wilson tensor.
    """
    c = np.asarray(full, dtype=np.complex128)
    lepton = [i for doublet in basis.leptons for i in doublet]
    higgs = list(basis.higgs)
    block = c[np.ix_(lepton, lepton, higgs, higgs)]
    # [p,i,q,j,A,B] -> [p,q,i,j,A,B]
    block = block.reshape(3, 2, 3, 2, 4, 4).transpose(0, 2, 1, 3, 4, 5)
    c5, relative = project_master(block, enforce_pure=False)
    absolute = relative * np.linalg.norm(block)
    return c5, float(absolute)


def exploratory_preflight(*args, **kwargs):
    kwargs['allow_projection_residual'] = True
    return ORIGINAL_PREFLIGHT(*args, **kwargs)


def exploratory_rough_case(*args, **kwargs):
    kwargs['allow_projection_residual'] = True
    result = ORIGINAL_ROUGH_CASE(*args, **kwargs)
    if result.get('status') == numerics_module.PROJECTED_STATUS:
        result = dict(result)
        result['status'] = 'RoughNumericalBenchmark'
        result['exploratory_nonweinberg_projection'] = True
    return result


def run_with_warm_starts(argv=None):
    """Use prior best UV Yukawas as the first seed of each independent fit.

    Existing optimizer, physics and original benchmark second seed are unchanged.
    The prior result is read-only; the caller must choose a new --output.
    """
    control = argparse.ArgumentParser(add_help=False)
    control.add_argument('--continue-from', type=Path, required=True,
                         help='Previous all-16 results directory (read-only)')
    args, remaining = control.parse_known_args(argv)
    prior_root = args.continue_from
    if not prior_root.is_absolute():
        prior_root = ROOT / prior_root
    prior_root = prior_root.resolve()
    if not prior_root.is_dir():
        control.error('Continuation results directory not found: ' + str(prior_root))

    original_seed_pair = fit_module.seed_pair
    original_optimize_one = runner_module.optimize_one

    def resumed_optimize_one(evaluate, cfg, *, model, **kwargs):
        path = prior_root / model / 'result.json'
        if not path.is_file():
            raise fit_module.PipelineBlocked('Prior result missing for ' + model + ': ' + str(path))
        doc = json.loads(path.read_text(encoding='utf-8-sig'))
        if doc.get('model') != model or not isinstance(doc.get('best'), dict):
            raise fit_module.PipelineBlocked('Prior result is not a valid fitted point: ' + str(path))
        best = doc['best']
        if 'x_scaled' in best:
            x = np.asarray(best['x_scaled'], dtype=float)
        else:
            y1 = np.asarray(best['UV_y1'], dtype=float)
            y2 = np.asarray(best['UV_y2'], dtype=float)
            if y1.shape != (3, 3) or y2.shape != (3, 3):
                raise fit_module.PipelineBlocked('Prior Yukawa matrices are not 3x3: ' + str(path))
            x = np.r_[y1.ravel(), y2.ravel()] / fit_module.YUKAWA_SCALE
        if x.shape != (18,) or not np.isfinite(x).all():
            raise fit_module.PipelineBlocked('Invalid 18 real Yukawa coordinates: ' + str(path))
        cap = float(kwargs['yukawa_bound']) / fit_module.YUKAWA_SCALE
        if np.max(np.abs(x)) >= cap:
            raise fit_module.PipelineBlocked(
                f'{model}: prior best max abs Yukawa {np.max(np.abs(x))*fit_module.YUKAWA_SCALE:.6g} '
                f'is not strictly below --yukawa-bound {kwargs["yukawa_bound"]:.6g}')
        old_chi = float(best['chi2'])
        if not np.isfinite(old_chi) or old_chi < 0:
            raise fit_module.PipelineBlocked('Invalid prior best chi2: ' + str(path))
        print(f'WARM START {model}: prior chi2={old_chi:.8g}, '
              f'max |Y|={np.max(np.abs(x))*fit_module.YUKAWA_SCALE:.5g}; '
              'first seed=previous best, second seed=original D_-2 benchmark', flush=True)

        def seeded(cfg_inner, seed):
            if seed == 'A_-4':
                return x.copy()
            return original_seed_pair(cfg_inner, seed)

        # Only one model is optimised at a time in this runner; restore the
        # seed function after each fit, including on exceptions.
        with patch.object(fit_module, 'seed_pair', seeded):
            result = original_optimize_one(evaluate, cfg, model=model, **kwargs)
        result['continuation'] = {
            'source_result': str(path), 'prior_best_chi2': old_chi,
            'first_seed_is_prior_best': True,
            'second_seed_is_original_D_minus2': True,
            'previous_best_retained_in_old_result': True,
        }
        return result

    with patch.object(runner_module, 'optimize_one', resumed_optimize_one):
        return runner_module.main(remaining)


if __name__ == '__main__':
    print('EXPLORATORY v3 CONTINUATION: B/C portal + L2 C5; previous best Yukawas as seeds.', flush=True)
    with patch.object(builder_module, 'quartic_unit_tensors', corrected_quartic_unit_tensors), \
         patch.object(workflow, 'project_c5', consistent_project_c5), \
         patch.object(rough_stream2, 'project_c5', consistent_project_c5), \
         patch.object(preflight_module, 'preflight_pair', exploratory_preflight), \
         patch.object(numerics_module, 'run_rough_case', exploratory_rough_case):
        raise SystemExit(run_with_warm_starts())
