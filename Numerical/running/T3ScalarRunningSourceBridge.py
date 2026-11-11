"""Worker B source-aware A_-4/D_-2 running scalar background adapter.

Uses pre-existing cached Batch3B quartic response (never regenerates native
invariants) and the actual EFT1 running snapshot. Source HHSS must be supplied
by the caller, not replaced by candidate values. No production solver changes.

Conventions: V2=1/2 phi.T M2 phi; phi=(H4,S1,S2), complex S components
are represented as (sqrt(2) Re S, sqrt(2) Im S) interleaved.
For D_-2, the SOURCE native L_int=+mu (Sdag)^2+h.c. gives
V2=m2 Sdag S-[mu (Sdag)^2+h.c.], hence [[m2-2Re(mu),-2Im(mu)],
[-2Im(mu),m2+2Re(mu)]] for each component.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Mapping
import json
import math
import numpy as np

from Numerical.running.T3GeneralScalarBackground import (
    AffineTensor, ScalarBackground, build_scalar_background,
)

SUPPORTED = ('A_-4', 'D_-2')
MASS_KEYS = ('mS1Sq', 'mS2Sq')


def _number(value: Any, name: str, *, real: bool = False) -> complex:
    if isinstance(value, Mapping):
        if set(value) != {'real', 'imag'}:
            raise ValueError(f'{name}: expected real,imag mapping')
        z = complex(float(value['real']), float(value['imag']))
    else:
        z = complex(value)
    if not math.isfinite(z.real) or not math.isfinite(z.imag):
        raise ValueError(f'{name}: nonfinite')
    if real and abs(z.imag) > 1e-12:
        raise ValueError(f'{name}: expected real')
    return z


def source_quadratic_tensor(model: str, state: Mapping[str, Any]) -> np.ndarray:
    """Unbroken real-scalar Hessian; Higgs block is unspecified (zero)."""
    if model not in SUPPORTED:
        raise ValueError('Source mass convention only established here for A_-4/D_-2')
    for key in MASS_KEYS:
        if key not in state:
            raise ValueError(f'missing running mass {key}')
    a, b = (_number(state[k], k, real=True).real for k in MASS_KEYS)
    m = np.zeros((12, 12), float)
    d1 = 1 if model == 'A_-4' else 3
    m[4:4+2*d1,4:4+2*d1] = a*np.eye(2*d1)
    m[4+2*d1:,4+2*d1:] = b*np.eye(8-2*d1)
    if model == 'D_-2':
        name = 'muT3GeneralS2Sq'
        if name not in state:
            raise ValueError(f'D_-2: additional quadratic direction {name} missing')
        mu = _number(state[name], name)
        # Native Matchete S2 complex component is one SU2 singlet, r,i.
        i = 4+2*d1
        m[i:i+2,i:i+2] = [[b-2*mu.real, -2*mu.imag],
                         [-2*mu.imag, b+2*mu.real]]
    return m


def mass_affine(model: str, reference: Mapping[str, Any]) -> AffineTensor:
    """Exact mass response; no native scalar parameter is silently omitted."""
    names = MASS_KEYS + (('muT3GeneralS2Sq',) if model == 'D_-2' else ())
    for name in names:
        if name not in reference:
            raise ValueError(f'reference missing quadratic {name}')
    base = source_quadratic_tensor(model, reference)
    cols = []
    dirs = []
    for name in names:
        z = _number(reference[name], name)
        for axis in (('real','imag') if name.startswith('mu') else ('real',)):
            sample = dict(reference)
            dz = 1j if axis == 'imag' else 1
            new = z + dz
            sample[name] = {'real':new.real,'imag':new.imag} if name.startswith('mu') else new.real
            cols.append(source_quadratic_tensor(model,sample)-base)
            dirs.append((name,axis))
    return AffineTensor(dict(reference), tuple(dirs), base, np.stack(cols), 2)


def cached_quartic_affine(context: Any) -> AffineTensor:
    """Adapt checked Batch3B response WITHOUT asserting physical normalization."""
    q = np.asarray(context.reference_tensor, float)
    cols = np.asarray(context.columns, float)
    return AffineTensor(dict(context.reference), tuple(map(tuple, context.directions)), q, cols, 4)


def build_from_cached_context(model: str, scale_gev: float,
                              state: Mapping[str, Any], context: Any,
                              *, hhss_source: np.ndarray) -> ScalarBackground:
    """Build at a supplied running EFT1 scale with source HHSS injection.

    A context is an existing CandidateQuarticContext loaded with its source
    provenance checks; the caller must supply independently built source HHSS
    at the *same* scale and in the same fixed 12-real-coordinate basis.
    """
    if model not in SUPPORTED or context.model != model:
        raise ValueError('model/context mismatch or unsupported scalar source')
    if hhss_source is None:
        raise ValueError('source HHSS mandatory for this conditional bridge')
    q = cached_quartic_affine(context)
    m = mass_affine(model, context.reference)
    return build_scalar_background(
        model, scale_gev, state, q, m,
        hhss_source=np.asarray(hhss_source),
        expected_quartic_names={name for name, _ in q.directions},
        expected_quadratic_names=set(MASS_KEYS) |
          ({'muT3GeneralS2Sq'} if model == 'D_-2' else set()),
        certificate='Batch3B candidate non-HHSS + source HHSS; unbroken masses; not certified physical potential',
    )


def load_saved_running_state(project: Path, model: str) -> tuple[float, dict]:
    """Read existing renormalisable EFT1@S benchmark without modifying it."""
    if model not in SUPPORTED:
        raise ValueError('model not supported')
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import load_running_snapshot
    filename = Path(project)/'output/general_observable/all16_v1'/model/'running.json'
    snapshot = load_running_snapshot(filename, model)
    return float(snapshot.mu_s), dict(snapshot.couplings)


def compare_masses(mass: np.ndarray, expected: Any) -> dict:
    """Numerical endpoint comparison against supplied historical eigenvalues."""
    computed = np.sort(np.linalg.eigvalsh(np.asarray(mass, float)[4:,4:]))
    reference = np.sort(np.asarray(expected, float))
    if reference.shape != (8,) or not np.isfinite(reference).all():
        raise ValueError('expected eight real scalar mass-squared eigenvalues')
    return {'eigenvalues_gev2': computed.tolist(),
            'reference_gev2': reference.tolist(),
            'maximum_absolute_difference_gev2': float(np.max(abs(computed-reference))),
            'maximum_relative_difference': float(np.max(abs(computed-reference))/max(float(np.max(abs(reference))),1e-100))}
