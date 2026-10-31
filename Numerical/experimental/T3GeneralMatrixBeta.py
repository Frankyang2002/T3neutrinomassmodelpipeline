"""Additive full-matrix UV RGBeta integrator for the general-scalar T3 basis.

Uses the repository's existing RGBeta InputForm parser, packed baseline state,
and Stage29 extra quartics. The optional Yukawa matrix is an independent
3x3 complex state, *not* a spurion proportional to y1 or y2. No new physics
beta function is introduced here. This is a candidate RGBeta-coordinate flow;
Matchete/RGBeta physical component equivalence remains an external gate.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

LOOP = 16.0 * math.pi**2
EXTRA_MODELS = {
    'A_-2': 'y2Conj', 'A_+0': 'y1Conj',
    'D_-2': 'y2Conj', 'D_+0': 'y1Conj',
    'E_-2': 'y2Conj', 'E_+0': 'y1Conj',
}


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError('T3 general UV: ' + msg)


def finite_complex(value: Any, label: str) -> complex:
    z = complex(value)
    require(math.isfinite(z.real) and math.isfinite(z.imag), 'nonfinite ' + label)
    return z


def finite_matrix(value: Any, label: str) -> np.ndarray:
    a = np.asarray(value, dtype=np.complex128)
    require(a.shape == (3, 3) and bool(np.isfinite(a).all()),
            label + ' must be a finite complex 3x3 matrix')
    return a


def encode_matrix(a: Any) -> dict:
    z = finite_matrix(a, 'encoded matrix')
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def decode_matrix(raw: Any, label: str) -> np.ndarray:
    require(isinstance(raw, dict) and set(raw) == {'real', 'imag'},
            label + ' must have real/imag keys')
    return finite_matrix(np.asarray(raw['real'], dtype=float) +
                         1j*np.asarray(raw['imag'], dtype=float), label)


def trial_extra_seed(y1: np.ndarray, y2: np.ndarray, fraction: float) -> np.ndarray:
    """Nonphysical RGBeta-coordinate benchmark: an *independent* matrix.

    Zero is the only non-arbitrary seed; nonzero fraction is explicitly a
    numerical sensitivity control, never a derived physical UV coupling.
    """
    f = float(fraction)
    require(math.isfinite(f) and 0 <= f <= 1, 'invalid extra Yukawa seed fraction')
    a = finite_matrix(y1, 'y1'); b = finite_matrix(y2, 'y2')
    scale = max(float(np.max(abs(a))), float(np.max(abs(b))), 1e-6)
    template = np.array([[1, .2+.15j, -.1], [.3, .8, .2j], [-.05j, .1, .6]], complex)
    return f * scale * template


@dataclass(frozen=True)
class ExtendedUVEquationSystem:
    """A source-validated 1-loop real ODE in the original baseline layout.

    Vector = [production baseline blocks; 2 per Stage29 new quartic;
              18 per independent extra Yukawa matrix (when applicable)].
    """
    model: str
    base_layout: Any
    extra_quartics: tuple[str, ...]
    extra_yukawa: str | None
    expressions: Mapping[str, str]
    source_metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        require(self.extra_yukawa == EXTRA_MODELS.get(self.model),
                'extra Yukawa inventory differs from model catalogue')
        names = self.all_names
        require(len(set(names)) == len(names), 'duplicate coupling blocks')
        require(set(self.expressions) == set(names),
                'RGBeta equations mismatch packed UV state: missing=' +
                str(sorted(set(names)-set(self.expressions))) + ' unexpected=' +
                str(sorted(set(self.expressions)-set(names))))
        require(all(isinstance(v, str) and v.strip() and
                    all(token not in v for token in ('$Failed', '$Aborted', 'Missing[', 'BetaTerm[', 'Indeterminate'))
                    for v in self.expressions.values()), 'invalid/unresolved native beta expression')
        meta = self.source_metadata
        expect = {'GeneralScalarModel': self.model,
                  'GeneralScalarStage': 'UV',
                  'ReportBetaConvention': '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
                  'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection'}
        require(all(meta.get(k) == v for k, v in expect.items()),
                'native RGBeta source metadata mismatch')
        require(set(meta.get('GeneralScalarNewBetaKeys', ())) == set(self.extra_quartics),
                'new quartic inventory disagrees with native source')
        require(meta.get('GeneralScalarQuadraticMassMixing') == 'NotIncluded',
                'expected UV source without extra quadratic betas')

    @property
    def all_names(self) -> tuple[str, ...]:
        return tuple(self.base_layout.names) + self.extra_quartics + (
            (self.extra_yukawa,) if self.extra_yukawa else ())

    @property
    def size(self) -> int:
        return self.base_layout.size + 2*len(self.extra_quartics) + (18 if self.extra_yukawa else 0)

    def initial_vector(self, baseline_vector: np.ndarray,
                       quartics: Mapping[str, complex],
                       extra_matrix: np.ndarray | None = None) -> np.ndarray:
        baseline = np.asarray(baseline_vector, dtype=float)
        require(baseline.shape == (self.base_layout.size,) and np.isfinite(baseline).all(),
                'baseline vector has wrong size/values')
        require(set(quartics) == set(self.extra_quartics), 'quartic seed inventory mismatch')
        side = [v for name in self.extra_quartics
                for v in (finite_complex(quartics[name], name).real,
                          finite_complex(quartics[name], name).imag)]
        if self.extra_yukawa:
            require(extra_matrix is not None, 'missing extra Yukawa boundary matrix')
            m = finite_matrix(extra_matrix, self.extra_yukawa)
            side.extend(m.real.reshape(-1)); side.extend(m.imag.reshape(-1))
        else:
            require(extra_matrix is None, 'unexpected extra Yukawa boundary')
        v = np.concatenate((baseline, np.asarray(side, dtype=float)))
        require(v.shape == (self.size,), 'packed state size mismatch')
        return v

    def split(self, vector: np.ndarray) -> tuple[np.ndarray, dict[str, complex], np.ndarray | None]:
        v = np.asarray(vector, dtype=float)
        require(v.shape == (self.size,) and np.isfinite(v).all(), 'invalid ODE vector')
        n = self.base_layout.size
        qs = {key: complex(v[n+2*i], v[n+2*i+1]) for i, key in enumerate(self.extra_quartics)}
        n += 2*len(qs)
        extra = (v[n:n+9].reshape(3, 3) + 1j*v[n+9:n+18].reshape(3, 3)
                 if self.extra_yukawa else None)
        return v[:self.base_layout.size], qs, extra

    def environment(self, t: float, vector: np.ndarray) -> dict[str, Any]:
        from Numerical.core.RGBetaEvaluator import state_environment
        from Numerical.core.StateVector import unpack_uv_state
        base, quartics, extra = self.split(vector)
        state = unpack_uv_state(base, self.base_layout, mu_gev=float(math.exp(t)))
        env = state_environment(state)
        env.update(quartics)
        if self.extra_yukawa:
            env[self.extra_yukawa] = extra
        return env

    def derivative(self, t: float, vector: np.ndarray) -> np.ndarray:
        from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
        from Numerical.core.BetaVector import beta_values_to_derivative
        env = self.environment(t, vector)
        values = {}
        for key, expr in self.expressions.items():
            values[key] = evaluate_inputform_expression(expr, env)
        result = list(beta_values_to_derivative(
            {k: values[k] for k in self.base_layout.names}, self.base_layout,
            divide_by_loop_factor=True))
        for name in self.extra_quartics:
            z = finite_complex(values[name], name+' beta')
            require(np.ndim(values[name]) == 0, name+' beta must be a scalar')
            result.extend((z.real/LOOP, z.imag/LOOP))
        if self.extra_yukawa:
            m = finite_matrix(values[self.extra_yukawa], self.extra_yukawa+' beta')
            result.extend(m.real.reshape(-1)/LOOP)
            result.extend(m.imag.reshape(-1)/LOOP)
        dy = np.asarray(result, dtype=float)
        require(dy.shape == (self.size,) and np.isfinite(dy).all(),
                'nonfinite/wrong-size UV beta vector')
        return dy

    def summary(self, vector: np.ndarray) -> dict:
        _, quartics, extra = self.split(vector)
        return {
            'extra_quartic_values': {k: {'real': v.real, 'imag': v.imag} for k, v in quartics.items()},
            'extra_yukawa_matrix': encode_matrix(extra) if extra is not None else None,
            'extra_yukawa_frobenius_norm': float(np.linalg.norm(extra)) if extra is not None else None,
        }
