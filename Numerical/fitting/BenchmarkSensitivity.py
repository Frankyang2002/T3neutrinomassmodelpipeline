"""Controlled UV-parameter sensitivity trajectories for the T3 benchmark.

The three scaling directions preserve lambdaT3*(y1*y2+y2*y1) at the
input scale for positive real gamma.
T1: lambdaT3 -> gamma lambdaT3; y1,y2 -> gamma^(-1/2) y1,y2
T2: lambdaT3 -> gamma lambdaT3; y2 -> gamma^(-1) y2
T3: y1 -> gamma y1; y2 -> gamma^(-1) y2

Each point uses the existing numerical T3 trajectory and final-C5 evaluator.
Failed points are recorded, not interpolated or included in sensitivity ovals.
"""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping
import numpy as np
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.fitting.ScanCLI import build_uv_state_from_config


def _scale_complex_matrix_spec(spec: dict[str, Any], factor: float) -> None:
    real = np.asarray(spec['real'], dtype=float)
    imag = np.asarray(spec.get('imag', np.zeros_like(real)), dtype=float)
    spec['real'] = (real * float(factor)).tolist()
    spec['imag'] = (imag * float(factor)).tolist()


def _scale_complex_scalar_spec(spec: dict[str, Any], factor: float) -> None:
    spec['real'] = float(spec.get('real', 0.0)) * float(factor)
    spec['imag'] = float(spec.get('imag', 0.0)) * float(factor)


def scaled_config_payload(payload: Mapping[str, Any], direction: str, gamma: float) -> dict[str, Any]:
    gamma = float(gamma)
    if not np.isfinite(gamma) or gamma <= 0:
        raise ValueError('gamma must be finite and positive.')
    direction = str(direction).upper()
    if direction not in {'T1', 'T2', 'T3'}:
        raise ValueError('direction must be T1, T2 or T3.')
    result = deepcopy(dict(payload))
    ordinary = result['base_state']['ordinary']
    if direction == 'T1':
        _scale_complex_scalar_spec(ordinary['lambdaT3'], gamma)
        factor = gamma ** -0.5
        _scale_complex_matrix_spec(ordinary['y1'], factor)
        _scale_complex_matrix_spec(ordinary['y2'], factor)
    elif direction == 'T2':
        _scale_complex_scalar_spec(ordinary['lambdaT3'], gamma)
        _scale_complex_matrix_spec(ordinary['y2'], gamma ** -1)
    else:
        _scale_complex_matrix_spec(ordinary['y1'], gamma)
        _scale_complex_matrix_spec(ordinary['y2'], gamma ** -1)
    return result


def leading_prefactor_scale(direction: str, gamma: float) -> float:
    gamma = float(gamma)
    direction = str(direction).upper()
    if direction == 'T1': return gamma * gamma ** -0.5 * gamma ** -0.5
    if direction in {'T2', 'T3'}: return gamma * gamma ** -1
    raise ValueError('direction must be T1, T2 or T3.')


def _config_view(payload: Mapping[str, Any]) -> SimpleNamespace:
    raw = deepcopy(dict(payload))
    raw.setdefault('scan', {'bindings': {}})
    return SimpleNamespace(raw=raw)


def gamma_grid(gamma_min: float, gamma_max: float, n_points: int, spacing: str = 'log') -> np.ndarray:
    """Return odd-length gamma grid with gamma=1 represented exactly at centre."""
    lo, hi, n = float(gamma_min), float(gamma_max), int(n_points)
    if not (np.isfinite(lo) and np.isfinite(hi) and 0 < lo < 1 < hi and n >= 3 and n % 2):
        raise ValueError('Require 0 < gamma_min < 1 < gamma_max and odd n_points >= 3.')
    if spacing == 'log':
        # Allow nonsymmetric limits and still anchor gamma=1 at the middle.
        half = n // 2
        values = np.r_[np.geomspace(lo, 1.0, half + 1), np.geomspace(1.0, hi, half + 1)[1:]]
    elif spacing == 'linear':
        half = n // 2
        values = np.r_[np.linspace(lo, 1.0, half + 1), np.linspace(1.0, hi, half + 1)[1:]]
    else:
        raise ValueError('spacing must be log or linear.')
    values[n // 2] = 1.0
    return values


def run_t3_sensitivity_scan(*, base_payload: Mapping[str, Any], uv_rgbeta_path: Path,
                            intermediate_rgbeta_path: Path, final_weinberg_path: Path,
                            oscillation_target_path: Path) -> dict[str, Any]:
    sensitivity = base_payload.get('sensitivity', {})
    if not isinstance(sensitivity, Mapping):
        raise ValueError('sensitivity must be an object.')
    if not bool(sensitivity.get('enabled', True)):
        return {'status': 'Disabled'}
    gamma_min = float(sensitivity.get('gamma_min', 0.95))
    gamma_max = float(sensitivity.get('gamma_max', 1.05))
    n_points = int(sensitivity.get('n_points', 11))
    spacing = str(sensitivity.get('spacing', 'log')).lower()
    gammas = gamma_grid(gamma_min, gamma_max, n_points, spacing)
    uv_rgbeta = load_rgbeta_payload(Path(uv_rgbeta_path))
    intermediate_rgbeta = load_rgbeta_payload(Path(intermediate_rgbeta_path))
    target = OscillationFitTarget.from_json(Path(oscillation_target_path))
    scales = base_payload['scales']
    numerical = base_payload.get('numerical', {})
    mu_f = float(scales['mu_fermion_threshold_gev'])
    mu_s = float(scales['mu_scalar_threshold_gev'])
    mu_low = float(scales['mu_low_gev'])
    vev = float(numerical.get('vev_gev', 246.22))
    rtol = float(numerical.get('rtol', 1.0e-8))
    atol = float(numerical.get('atol', 1.0e-11))
    c5_evaluator = FinalC5TrajectoryEvaluator(Path(final_weinberg_path))
    results: dict[str, list[dict[str, Any]]] = {}
    for direction in ('T1', 'T2', 'T3'):
        points = []
        for gamma in gammas:
            entry: dict[str, Any] = {'gamma': float(gamma)}
            try:
                deformed = scaled_config_payload(base_payload, direction, float(gamma))
                view = _config_view(deformed)
                evaluator = make_t3_point_evaluator(
                    uv_state_builder=lambda _params, view=view: build_uv_state_from_config(view, {}),
                    c5_builder=c5_evaluator, uv_rgbeta_payload=uv_rgbeta,
                    eft1_rgbeta_payload=intermediate_rgbeta, target=target,
                    mu_fermion_threshold_gev=mu_f, mu_scalar_threshold_gev=mu_s,
                    mu_low_gev=mu_low, vev_gev=vev, rtol=rtol, atol=atol)
                fit = evaluator({})
                prediction = {name: float(value) for name, value in zip(fit.observable_names, fit.prediction)}
                if not np.isfinite(float(fit.chi2)) or any(not np.isfinite(v) for v in prediction.values()):
                    raise ValueError('Non-finite chi2 or observable')
                entry.update(status='Success', chi2=float(fit.chi2), prediction=prediction)
            except Exception as exc:
                entry.update(status='Failed', message=f'{type(exc).__name__}: {exc}')
            points.append(entry)
        results[direction] = points
    definitions = {
        'T1': 'lambdaT3 -> gamma lambdaT3; y1,y2 -> gamma^(-1/2) y1,y2',
        'T2': 'lambdaT3 -> gamma lambdaT3; y2 -> gamma^(-1) y2',
        'T3': 'y1 -> gamma y1; y2 -> gamma^(-1) y2',
    }
    failures = sum(p['status'] != 'Success' for pts in results.values() for p in pts)
    return {
        'status': 'Success',  # Completed scan; individual failed points remain explicitly tagged
        'gamma_min': gamma_min, 'gamma_max': gamma_max, 'n_points': n_points,
        'spacing': spacing, 'n_failed_points': failures,
        'directions': {d: {'definition': definitions[d], 'points': results[d]} for d in results},
        'target': {'ordering': target.ordering, 'observable_names': list(target.observable_names),
                   'central_values': target.central_values.tolist(), 'covariance': target.covariance.tolist(),
                   'source': target.source,
                   'note': 'Gaussian approximation; 1sigma/3sigma regions are not official Delta-chi2 contours.'},
    }
