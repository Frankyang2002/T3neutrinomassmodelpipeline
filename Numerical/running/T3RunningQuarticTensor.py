"""Running full real-scalar quartic tensor from an explicit response dictionary.

This module does not infer RGBeta/Matchete normalisations. Each supplied tensor
is the coefficient of +1 in the indicated RGBeta coupling component, expressed
in the same V4 = lambda_abcd phi_a phi_b phi_c phi_d / 4! convention used by
the seven-term Wilson RHS.

Only use dictionaries generated/checked against the actual local native source.
No implicit zero fallback for missing model-specific quartics.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping
import numpy as np

DIMS = {'A': (1, 3), 'B': (2, 2), 'C': (2, 2),
        'D': (3, 1), 'E': (3, 3)}
MIXED = frozenset({'B_-1', 'C_-1'})


class QuarticEvidenceError(ValueError):
    pass


def _check_model(model: str) -> int:
    if model not in {
        'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
        'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
        'E_+0', 'E_+2', 'E_-2', 'E_-4',
    }:
        raise QuarticEvidenceError(f'Unsupported model: {model}')
    return 4 + 2 * sum(DIMS[model[0]])


def _complex(value, name):
    if isinstance(value, Mapping):
        if set(value) != {'real', 'imag'}:
            raise QuarticEvidenceError(f'{name}: malformed complex object')
        value = complex(float(value['real']), float(value['imag']))
    else:
        value = complex(value)
    if not np.isfinite(value.real) or not np.isfinite(value.imag):
        raise QuarticEvidenceError(f'{name}: nonfinite coupling')
    return value


def _tensor(value, n: int, label: str):
    tensor = np.asarray(value, dtype=float)
    if tensor.shape != (n, n, n, n):
        raise QuarticEvidenceError(
            f'{label}: expected {(n,)*4}, got {tensor.shape}')
    if not np.isfinite(tensor).all():
        raise QuarticEvidenceError(f'{label}: nonfinite coefficients')
    for axis in (1, 2, 3):
        if not np.allclose(tensor, tensor.swapaxes(0, axis),
                           rtol=1e-11, atol=1e-12):
            raise QuarticEvidenceError(f'{label}: Bose symmetry failure')
    return tensor


@dataclass(frozen=True)
class RunningQuarticBasis:
    model: str
    symbols: tuple[str, ...]
    directions: tuple[tuple[str, str], ...]
    responses: np.ndarray
    constant: np.ndarray
    source: str

    @classmethod
    def from_payload(cls, payload: Mapping, expected_symbols: set[str]):
        """Fail-closed import of an explicitly sourced real-Cartesian basis.

        Format:
          model: one of the 16
          convention: 'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!'
          basis: 'fixed-original-Cartesian'
          source: nonempty evidence identifier
          constant: full n^4 tensor (typically zeros)
          responses: [{symbol,component:'real'|'imag',tensor:n^4},...]
        """
        if payload.get('format') != 't3-running-real-quartic-response-v1':
            raise QuarticEvidenceError('Unknown tensor response file format')
        model = payload.get('model')
        n = _check_model(model)
        if payload.get('convention') != 'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!':
            raise QuarticEvidenceError('Quartic potential normalization unspecified')
        if payload.get('basis') != 'fixed-original-Cartesian':
            raise QuarticEvidenceError('Moving-basis scalar tensors are forbidden')
        source = payload.get('source')
        if not isinstance(source, str) or not source.strip():
            raise QuarticEvidenceError('Missing native tensor provenance')
        rows = payload.get('responses')
        if not isinstance(rows, list):
            raise QuarticEvidenceError('Missing explicit coupling responses')
        directions, matrices = [], []
        for row in rows:
            symbol, part = row.get('symbol'), row.get('component')
            if not isinstance(symbol, str) or part not in ('real', 'imag'):
                raise QuarticEvidenceError('Invalid coupling direction')
            direction = (symbol, part)
            if direction in directions:
                raise QuarticEvidenceError(f'Duplicate direction {direction}')
            directions.append(direction)
            matrices.append(_tensor(row.get('tensor'), n, str(direction)))
        symbols = {s for s, _ in directions}
        if symbols != set(expected_symbols):
            missing = sorted(set(expected_symbols) - symbols)
            excess = sorted(symbols - set(expected_symbols))
            raise QuarticEvidenceError(
                f'Model coupling mismatch: missing={missing}, excess={excess}')
        if not all((s, 'real') in directions for s in expected_symbols):
            raise QuarticEvidenceError('Missing real response of required quartic')
        offset = _tensor(payload.get('constant'), n, 'constant')
        return cls(model, tuple(sorted(symbols)), tuple(directions),
                   np.stack(matrices) if matrices else np.empty((0,)*1 + (n,)*4),
                   offset, source)

    def evaluate(self, state: Mapping) -> np.ndarray:
        missing = set(self.symbols) - set(state)
        if missing:
            raise QuarticEvidenceError(f'Missing running quartics: {sorted(missing)}')
        values = {symbol: _complex(state[symbol], symbol)
                  for symbol in self.symbols}
        coords = np.asarray([
            getattr(values[symbol], part) for symbol, part in self.directions
        ], dtype=float)
        out = self.constant + np.einsum(
            'k,kabcd->abcd', coords, self.responses, optimize=True)
        if not np.isfinite(out).all():
            raise QuarticEvidenceError('Nonfinite running quartic tensor')
        return out

    def validate_hhss(self, state: Mapping, expected: np.ndarray,
                      *, rtol: float = 1e-9):
        """Check pre-existing HHSS implementation at this exact running state."""
        n = _check_model(self.model)
        nheavy = n - 4
        q = np.asarray(expected, dtype=float)
        if q.shape != (4, 4, nheavy, nheavy):
            raise QuarticEvidenceError('HHSS reference has wrong shape')
        got = self.evaluate(state)[:4, :4, 4:, 4:]
        rel = float(np.linalg.norm(got - q) /
                    max(float(np.linalg.norm(q)), 1e-300))
        if rel > rtol:
            raise QuarticEvidenceError(
                f'HHSS source parity failed: {rel:.3e}')
        return rel


def expected_scalar_symbols(registry: Mapping, model: str,
                            baseline_symbols: set[str]) -> set[str]:
    """Demand all registry extra and preexisting couplings, plus baseline."""
    _check_model(model)
    m = registry['models'][model]
    new = {item['name'] for item in m['new_rgbeta_symbols']}
    old = {item['existing_name'] for item in m['preexisting_rgbeta_symbols']}
    if len(new) + len(old) != m['total_additional_quartics']:
        raise QuarticEvidenceError('Registry mismatch')
    if not old <= baseline_symbols:
        raise QuarticEvidenceError(
            f'Missing preexisting RGBeta couplings: {sorted(old-baseline_symbols)}')
    return set(baseline_symbols) | new
