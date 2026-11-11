"""Cached real-scalar V4 fourth-derivative tensor response (A_-4, D_-2).

For an affine-linear scalar potential in the real/imaginary coupling coordinates,
construct Q(c)=Q_ref + sum_i (c_i-c_ref_i) dQ_i, with Q_abcd = d^4 V/dphi^4.
The constructor accepts the existing repository's source-aware Q evaluator;
this module neither modifies the RGBeta basis nor certifies its mapping.

All 12 real scalar components are retained: four Higgs and eight BSM components.
Complex coupling directions are independently probed. A user of this module
MUST keep the heavy-field real basis fixed while evaluating Q on a running
trajectory, and rotate at matching when appropriate.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping
import math
import numpy as np

SHAPE = (12, 12, 12, 12)


def scalar(value: Any) -> complex:
    if isinstance(value, Mapping):
        if 'real' not in value or 'imag' not in value:
            raise ValueError('Complex mapping requires real and imag')
        x = complex(float(value['real']), float(value['imag']))
    else:
        x = complex(np.asarray(value).item())
    if not (math.isfinite(x.real) and math.isfinite(x.imag)):
        raise ValueError('Non-finite scalar coupling')
    return x


def as_tensor(value: Any) -> np.ndarray:
    q = np.asarray(value)
    if q.shape != SHAPE or not np.isfinite(q).all():
        raise ValueError('Expected finite real-scalar tensor with shape (12,)*4')
    if np.iscomplexobj(q):
        if np.max(np.abs(q.imag)) > 1e-10 * max(float(np.max(np.abs(q.real))), 1e-100):
            raise ValueError('Quartic tensor has complex part')
        q = q.real
    return np.asarray(q, dtype=float)


def symmetry_residual(q: np.ndarray) -> float:
    q = as_tensor(q)
    norm = max(float(np.linalg.norm(q)), 1e-100)
    return max(float(np.linalg.norm(q - q.transpose(p))) / norm for p in (
        (1, 0, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2)))


def _set(x: Mapping[str, Any], name: str, value: complex, is_complex: bool) -> dict:
    z = dict(x)
    z[name] = {'real': float(value.real), 'imag': float(value.imag)} if is_complex else float(value.real)
    return z


@dataclass
class ScalarQuarticLinearResponse:
    """Fixed-basis affine response with exact direct-reconstruction tests."""
    reference: dict[str, Any]
    reference_tensor: np.ndarray
    directions: tuple[tuple[str, str], ...]
    columns: np.ndarray
    direct: Callable[[Mapping[str, Any]], np.ndarray]

    @classmethod
    def construct(cls, reference: Mapping[str, Any],
                  direct: Callable[[Mapping[str, Any]], np.ndarray],
                  names: list[str] | tuple[str, ...],
                  complex_names: set[str] | frozenset[str], *,
                  tolerance: float = 2e-9) -> 'ScalarQuarticLinearResponse':
        ref = dict(reference)
        names = tuple(names)
        if len(names) != len(set(names)) or not names:
            raise ValueError('Coupling names must be unique and nonempty')
        if not set(complex_names) <= set(names):
            raise ValueError('Complex names not found in scalar inventory')
        q0 = as_tensor(direct(ref))
        if symmetry_residual(q0) > tolerance:
            raise ValueError('Direct reference tensor lacks bosonic index symmetry')
        columns = []
        dirs = []
        for name in names:
            if name not in ref:
                raise ValueError('Missing scalar coupling: ' + name)
            z = scalar(ref[name])
            if name not in complex_names and abs(z.imag) > 1e-12:
                raise ValueError('Expected real scalar coupling: ' + name)
            for axis, step in [('real', 1.), ('imag', 1j)] if name in complex_names else [('real', 1.)]:
                test = _set(ref, name, z + step, name in complex_names)
                col = as_tensor(direct(test)) - q0
                if symmetry_residual(col) > tolerance:
                    raise ValueError(f'Quartic derivative for {name}.{axis} not symmetric')
                columns.append(col)
                dirs.append((name, axis))
        result = cls(ref, q0, tuple(dirs), np.stack(columns), direct)
        # Independent multi-coupling, non-unit interpolation test. One-hot
        # probes alone cannot reveal nonlinear or cross-coupling terms.
        trial = dict(ref)
        for j, name in enumerate(names):
            z = scalar(ref[name])
            dr = 0.137 * math.sin(j + 1)
            di = 0.083 * math.cos(j + 1) if name in complex_names else 0.
            trial = _set(trial, name, z + complex(dr, di), name in complex_names)
        comparison = result.compare(trial)
        if comparison > tolerance:
            raise ValueError(f'Running quartic basis not affine-linear: {comparison:.4g}')
        return result

    def tensor(self, values: Mapping[str, Any]) -> np.ndarray:
        increments = []
        for name, axis in self.directions:
            if name not in values:
                raise ValueError('Running state missing scalar coupling: ' + name)
            delta = scalar(values[name]) - scalar(self.reference[name])
            increments.append(delta.imag if axis == 'imag' else delta.real)
        return self.reference_tensor + np.tensordot(np.asarray(increments, float), self.columns, axes=(0, 0))

    def compare(self, values: Mapping[str, Any]) -> float:
        exact = as_tensor(self.direct(values))
        approx = self.tensor(values)
        return float(np.linalg.norm(exact - approx) / max(np.linalg.norm(exact), 1e-100))

    def metadata(self) -> dict:
        return {'tensor_shape': list(SHAPE), 'directions': [list(x) for x in self.directions],
                'reference_norm': float(np.linalg.norm(self.reference_tensor)),
                'direction_norms': [float(np.linalg.norm(c)) for c in self.columns],
                'symmetry_residual': symmetry_residual(self.reference_tensor),
                'physical_normalisation_certified': False,
                'scheme': 'cached source-aware affine real-quartic tensor; fixed scalar basis'}
