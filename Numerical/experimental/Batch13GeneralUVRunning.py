"""Batch 13: full candidate RGBeta UV evolution, all 16 split-scalar T3 models.

Intended path: Numerical/experimental/Batch13GeneralUVRunning.py

Uses canonical T3UVState/StateVector and the existing RGBeta parser. This is a
renormalisable UV-only calculation, NOT a matched neutrino-mass prediction.
The native RGBeta <-> Matchete scalar normalisation remains uncertified.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import json
import math

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.core.GeneralCouplingLayout import layout_from_rgbeta
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
from Numerical.core.StateVector import pack_uv_state, unpack_uv_state

EXPECTED_FORMAT = 'general-t3-quadratic-rgbeta-batch13-restricted-yukawa-v2'
EXPECTED_STATUS = 'SuccessCandidateMassBeta'
EXPECTED_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
REPRESENTATIONS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
                   'D': (3, 1, 2), 'E': (3, 3, 2)}
FORBIDDEN_YUKAWAS = frozenset({'y1Conj', 'y2Conj'})
OPTIONAL_REAL = ('lambdaH1Adj', 'lambdaH2Adj', 'lambdaS1Adj',
                 'lambdaS2Adj', 'lambda12Adj', 'lambda12Cross')
OPTIONAL_COMPLEX = ('lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
                    'lambdaS1bar2S2bar2', 'lambdaS1barS2S2bar2',
                    'lambdaS1S1bar2S2bar', 'lambdaHHdagS1barS2barCross')


def validate_uv_export(payload: Mapping[str, Any], model: str) -> tuple[str, ...]:
    """Reject mismatched, incomplete or unrecognised candidate UV exports."""
    if payload.get('format') != EXPECTED_FORMAT or payload.get('status') != EXPECTED_STATUS:
        raise ValueError(
            f"{model}: unsupported candidate UV export format/status: "
            f"format={payload.get('format')!r}, status={payload.get('status')!r}; "
            f"expected format={EXPECTED_FORMAT!r}, status={EXPECTED_STATUS!r}"
        )
    meta = payload.get('metadata')
    if not isinstance(meta, Mapping):
        raise ValueError(f'{model}: missing metadata')
    prefix, separator, signed_alpha = model.partition('_')
    if not separator or prefix not in REPRESENTATIONS:
        raise ValueError(f'Unknown ordinary T3 model: {model}')
    alpha = int(signed_alpha)
    d1, d2, df = REPRESENTATIONS[prefix]
    required = {'GeneralScalarModel': model, 'GeneralScalarStage': 'UV',
                'dS1': d1, 'dS2': d2, 'dF': df, 'alpha': alpha,
                'ReportBetaConvention': EXPECTED_CONVENTION,
                'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete'}
    for key, wanted in required.items():
        if meta.get(key) != wanted:
            raise ValueError(f'{model}: metadata {key}: {meta.get(key)!r} != {wanted!r}')
    if payload.get('physical_mapping_certified') is not False:
        raise ValueError(f'{model}: unexpected physical mapping certification state')
    specs = layout_from_rgbeta(payload)
    names = tuple(spec.name for spec in specs)
    if FORBIDDEN_YUKAWAS.intersection(names):
        raise ValueError(f'{model}: forbidden independent conjugate Yukawa beta')
    required_names = {'gY', 'g2', 'g3', 'yu', 'yd', 'ye', 'lambdaH',
                      'y1', 'y2', 'MF', 'mS1Sq', 'mS2Sq',
                      'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2',
                      'lambda12', 'lambdaT3'}
    missing = required_names.difference(names)
    if missing:
        raise ValueError(f'{model}: missing required UV beta keys {sorted(missing)}')
    declared = set(meta.get('GeneralScalarNewBetaKeys', []))
    if declared - set(names):
        raise ValueError(f'{model}: declared general quartic beta absent from export')
    if meta.get('QuadraticKind') != 'none':
        symbol = meta.get('QuadraticMassSymbol')
        if not isinstance(symbol, str) or symbol not in names:
            raise ValueError(f'{model}: quadratic beta missing')
    return names


def _diag(values: tuple[float, float, float]) -> np.ndarray:
    return np.diag(np.asarray(values, dtype=complex))


def illustrative_initial_state(payload: Mapping[str, Any], model: str,
                               *, mu_gev: float = 1e7) -> T3UVState:
    """Make an explicitly *illustrative* native-RGBeta benchmark.

    New quartics and non-diagonal quadratic parameters are seeded at zero;
    they are nevertheless included in the ODE and may be radiatively sourced.
    All scalar couplings here are RGBeta-native candidate parameters, NOT
    matched Matchete coefficients. Replace these seeds for physical studies.
    """
    names = set(validate_uv_export(payload, model))
    family, alpha_string = model.split('_', 1)
    d1, d2, df = REPRESENTATIONS[family]
    rep = T3Representation(d1, d2, df, int(alpha_string))
    sm = SMNumericalState(
        gY=0.38, g2=0.62, g3=0.92, lambdaH=0.13,
        yu=_diag((1.5e-5, 7e-3, 0.80)),
        yd=_diag((3e-5, 5e-4, 2e-2)),
        ye=_diag((3e-6, 6e-4, 1e-2)),
    )
    base = dict(mu_gev=float(mu_gev), representation=rep, sm=sm,
                y1=_diag((0.01, 0.01, 0.01)),
                y2=_diag((0.01, 0.01, 0.01)),
                MF=_diag((3e5, 3.5e5, 4e5)),
                mS1Sq=(3e3)**2, mS2Sq=(3.2e3)**2,
                lambdaS1=0.10, lambdaS2=0.10, lambdaH1=0.05,
                lambdaH2=0.05, lambda12=0.04, lambdaT3=0.02+0j)
    for name in OPTIONAL_REAL:
        base[name] = 0.01 if name in names else None
    for name in OPTIONAL_COMPLEX:
        base[name] = 0.01 + 0j if name in names else None
    canonical = {'gY','g2','g3','yu','yd','ye','lambdaH',
                 'y1','y2','MF','mS1Sq','mS2Sq','lambdaS1','lambdaS2',
                 'lambdaH1','lambdaH2','lambda12','lambdaT3'}
    canonical.update(OPTIONAL_REAL)
    canonical.update(OPTIONAL_COMPLEX)
    new_names = names - canonical
    base['general_couplings'] = {name: 0j for name in new_names}
    state = T3UVState(**base).validated()
    _, layout = pack_uv_state(state)
    if set(layout.names) != names:
        raise ValueError(f'{model}: canonical state / beta inventory mismatch: '
                         f'missing={sorted(names-set(layout.names))}, '
                         f'extra={sorted(set(layout.names)-names)}')
    return state


def run_uv(payload: Mapping[str, Any], model: str, initial: T3UVState,
           *, mu_final_gev: float = 1e5, rtol: float = 1e-7,
           atol: float = 1e-10, max_nfev: int = 50000) -> dict[str, Any]:
    """Run the complete native-RGBeta candidate UV ODE down to mu_F.

    No F threshold matching is performed here. Fail closed on unknown beta
    syntax/shape and missing couplings; preserve Majorana MF symmetry.
    """
    names = set(validate_uv_export(payload, model))
    initial = initial.validated()
    if initial.representation.shared_scalar:
        raise ValueError('This flow requires both ordinary T3 scalars')
    if (initial.representation.d_s1, initial.representation.d_s2,
        initial.representation.d_f) != REPRESENTATIONS[model[0]] or \
            initial.representation.alpha != int(model.split('_')[1]):
        raise ValueError(f'{model}: initial state representation mismatch')
    y0, layout = pack_uv_state(initial)
    if set(layout.names) != names:
        raise ValueError(f'{model}: UV beta/state mismatch: '
                         f'missing={sorted(names-set(layout.names))}, '
                         f'extra={sorted(set(layout.names)-names)}')
    mu0, mu1 = float(initial.mu_gev), float(mu_final_gev)
    if not (np.isfinite(mu1) and mu0 > mu1 > 0):
        raise ValueError('Require mu_UV > mu_F > 0')
    betas = payload['report_betas']
    majorana = initial.representation.majorana_fermion
    mass_block = layout.block('MF') if majorana else None

    def symmetric_mass(v: np.ndarray) -> None:
        if mass_block is None:
            return
        n = int(np.prod(mass_block.shape))
        for start in (mass_block.start, mass_block.start + n):
            arr = v[start:start+n].reshape(mass_block.shape)
            arr[:] = (arr + arr.T) / 2

    nfev = 0
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        nonlocal nfev
        nfev += 1
        if nfev > max_nfev:
            raise RuntimeError(f'{model}: exceeded max_nfev={max_nfev}')
        vec = np.array(y, dtype=float, copy=True)
        symmetric_mass(vec)
        state = unpack_uv_state(vec, layout, mu_gev=math.exp(t))
        env = state_environment(state)
        evaluated = {}
        for name, expression in betas.items():
            try:
                evaluated[name] = evaluate_inputform_expression(expression, env)
            except Exception as error:
                raise ValueError(f'{model}: failed evaluating beta {name}: {error}') from error
        dy = beta_values_to_derivative(evaluated, layout, divide_by_loop_factor=True)
        symmetric_mass(dy)
        return dy

    # Validate the full RHS before invoking the adaptive solver.
    rhs(math.log(mu0), y0)
    solution = solve_ivp(rhs, (math.log(mu0), math.log(mu1)), y0,
                         method='DOP853', rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(f'{model}: UV integration failed: {solution.message}')
    final_vec = np.array(solution.y[:, -1], copy=True)
    symmetric_mass(final_vec)
    final = unpack_uv_state(final_vec, layout, mu_gev=mu1)
    return {'model': model, 'state': final, 'nfev': int(solution.nfev),
            'mu_initial_gev': mu0, 'mu_final_gev': mu1,
            'n_parameters_real': int(layout.size),
            'n_couplings': len(layout.names),
            'status': 'CandidateUVIntegrated',
            'matchete_rgbeta_normalisation_certified': False,
            'mass_eigenstate_threshold_matching_performed': False,
            'restricted_yukawa': True}


def serialise_state(state: T3UVState) -> dict[str, Any]:
    """JSON-compatible, lossless real/imag encoding for all ODE parameters."""
    from Numerical.core.RGBetaEvaluator import state_environment
    env = state_environment(state)
    out = {}
    for name, value in env.items():
        arr = np.asarray(value, dtype=complex)
        if arr.shape == ():
            z = complex(arr)
            out[name] = {'real': float(z.real), 'imag': float(z.imag)}
        else:
            out[name] = {'real': arr.real.tolist(), 'imag': arr.imag.tolist()}
    return out
