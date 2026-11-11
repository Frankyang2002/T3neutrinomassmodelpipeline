"""Exploratory non-HHSS quartic overlay for the EXISTING seven-term EFT1 solver.

The native Matchete coefficient for each invariant I_k is provisionally identified
with the RGBeta symbol in the same model/field-content block and invariant order.
This identification is NOT checked; calculations using it must say 'exploratory'.

Existing HHSS, Higgs and baseline scalar self-interactions continue to come
from the current general_quartic(model, running_state, HHSS) implementation.
Only independent, non-HHSS model-specific quartics are superimposed, to avoid
adding the HHSS channels twice. No scalar quadratic or field mixing is changed.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from collections import defaultdict
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp

from scripts.RGBetaTensorAlgebra import exact


class ExploratoryQuarticError(ValueError):
    pass


def _is_hhss(counts):
    if len(counts) != 6 or sum(counts) != 4:
        raise ExploratoryQuarticError('Invalid six-species field counts')
    return sum(counts[:2]) == 2 and sum(counts[2:]) == 2


def _read_json(path):
    b = Path(path).read_bytes()
    return json.loads(b), hashlib.sha256(b).hexdigest()


def _number(x):
    a = exact(x)
    if sp.simplify(sp.im(a)) != 0:
        raise ExploratoryQuarticError(f'Complex Cartesian response: {x}')
    y = float(sp.N(sp.re(a), 16))
    if not np.isfinite(y):
        raise ExploratoryQuarticError('Nonfinite Cartesian response')
    return y


def _expected_registry_by_block(entry):
    groups = defaultdict(list)
    for row in entry['new_rgbeta_symbols']:
        groups[tuple(row['counts'])].append(row['name'])
    for row in entry['preexisting_rgbeta_symbols']:
        groups[tuple(row['counts'])].append(row['existing_name'])
    if sum(len(v) for v in groups.values()) != entry['total_additional_quartics']:
        raise ExploratoryQuarticError('Registry slot count differs')
    return groups


def _matchete_index(symbol, counts):
    prefix = 'lambdaT3ExtraQ' + ''.join(map(str, counts)) + 'I'
    if not symbol.startswith(prefix):
        raise ExploratoryQuarticError(f'Matchete symbol {symbol} has unexpected block')
    tail = symbol[len(prefix):]
    if not tail.isdigit() or int(tail) <= 0:
        raise ExploratoryQuarticError('Matchete invariant index malformed')
    return int(tail)


def _response_tensor(record, part, n):
    payload = next((r for r in record['responses'] if r['component'] == part), None)
    if payload is None:
        raise ExploratoryQuarticError(f'Missing {part} response for {record["symbol"]}')
    out = np.zeros((n, n, n, n), dtype=float)
    seen = set()
    for row in payload['entries']:
        idx = tuple(row['indices'])
        if len(idx) != 4 or tuple(sorted(idx)) != idx or any(not isinstance(k,int) or k < 0 or k >= n for k in idx):
            raise ExploratoryQuarticError(f'Invalid sparse Cartesian indices: {idx}')
        if idx in seen:
            raise ExploratoryQuarticError(f'Duplicate sparse tensor index {idx}')
        seen.add(idx)
        component = _number(row['lambda_abcd'])
        for perm in set(permutations(idx)):
            out[perm] = component
    return out


@dataclass(frozen=True)
class ExploratoryQuarticOverlay:
    model: str
    symbols: tuple[str, ...]
    source: str
    directions: tuple[tuple[str,str], ...]
    responses: tuple[np.ndarray, ...]
    original_baseline: object
    hhss_provider: object
    n: int
    mapping_assumption: str = 'Matchete I_k equals RGBeta block-ordered coupling k'

    def evaluate(self, state):
        hhss = self.hhss_provider(state)
        base = np.asarray(self.original_baseline(self.model, state, hhss), dtype=float)
        if base.shape != (self.n,)*4:
            raise ExploratoryQuarticError(f'{self.model}: incorrect baseline tensor shape')
        out = base.copy()
        for (symbol, part), tensor in zip(self.directions, self.responses):
            if symbol not in state:
                raise ExploratoryQuarticError(f'{self.model}: missing running {symbol}')
            raw = state[symbol]
            # The native mixed RGBeta snapshot represents complex couplings
            # as {'real': ..., 'imag': ...}; accept ordinary complex/real too.
            if isinstance(raw, dict):
                if not {'real', 'imag'} <= raw.keys():
                    raise ExploratoryQuarticError(
                        f'{self.model}/{symbol}: malformed complex scalar dictionary')
                try:
                    c = complex(float(raw['real']), float(raw['imag']))
                except (TypeError, ValueError, OverflowError) as exc:
                    raise ExploratoryQuarticError(
                        f'{self.model}/{symbol}: invalid complex scalar values') from exc
            else:
                try:
                    c = complex(raw)
                except (TypeError, ValueError, OverflowError) as exc:
                    raise ExploratoryQuarticError(
                        f'{self.model}/{symbol}: invalid scalar coupling') from exc
            if not (np.isfinite(c.real) and np.isfinite(c.imag)):
                raise ExploratoryQuarticError(f'{symbol}: nonfinite running coupling')
            out += (c.real if part == 'real' else c.imag) * tensor
        if not np.isfinite(out).all():
            raise ExploratoryQuarticError('Nonfinite exploratory tensor')
        return out

    def validate_hhss(self, state, expected, *, rtol=1e-9):
        got = self.evaluate(state)[:4,:4,4:,4:]
        expect = np.asarray(expected, dtype=float)
        if got.shape != expect.shape:
            raise ExploratoryQuarticError('HHSS shape mismatch')
        rel = np.linalg.norm(got-expect)/max(np.linalg.norm(expect),1e-300)
        if rel > rtol:
            raise ExploratoryQuarticError(f'HHSS parity failed: {rel:.3e}')
        return float(rel)


def build_overlay(model, candidate_path, registry_path, original_baseline,
                  hhss_provider, *, required_quartics=()):
    """Build an explicit *provisional* model basis from current native exports.

    No false promise: this reads the existing Stage 22 Matchete candidate file;
    it does not pretend to export independent RGBeta tensors. The resulting
    basis is compatible with Batch 4 integrate_extended_eft1(), provided its
    adapter supplies the full extended running state at every RHS evaluation.
    """
    candidates, digest = _read_json(candidate_path)
    registry, _ = _read_json(registry_path)
    if candidates.get('format') != 't3-matchete-candidate-cartesian-quartics-v1' or candidates.get('rgbeta_mapping_certified') is not False:
        raise ExploratoryQuarticError('Expected explicitly uncertified Matchete candidate export')
    if model not in registry['models'] or model not in candidates['models']:
        raise ExploratoryQuarticError(f'Model missing from source files: {model}')
    entry = registry['models'][model]
    grouped = _expected_registry_by_block(entry)
    records = candidates['models'][model]['records']
    if len(records) != entry['total_additional_quartics']:
        raise ExploratoryQuarticError('Candidate count mismatch')
    # A record's six-species content is encoded in its stable Matchete symbol.
    dirs, tensors, all_symbols = [], [], set()
    d1,d2 = {'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}[model[0]]
    n=4+2*(d1+d2)
    matched = set()
    for counts, names in grouped.items():
        key = ''.join(map(str, counts))
        block = [r for r in records if r['symbol'].startswith('lambdaT3ExtraQ'+key+'I')]
        ordered = sorted(block, key=lambda r: _matchete_index(r['symbol'], counts))
        if len(ordered) != len(names) or [_matchete_index(r['symbol'],counts) for r in ordered] != list(range(1,len(names)+1)):
            raise ExploratoryQuarticError(f'{model}/{key}: incomplete or ambiguous Matchete basis')
        for symbol, record in zip(names, ordered):
            matched.add(record['symbol'])
            all_symbols.add(symbol)
            if record['cartesian_dimensions'] != n:
                raise ExploratoryQuarticError('Cartesian dimensions changed')
            if _is_hhss(counts):
                # Existing q_basis supplies HHSS, including these extra HHSS channels.
                continue
            for part in ('real','imag'):
                t = _response_tensor(record, part, n)
                # By construction a non-HHSS tensor must never touch HHSS.
                if np.any(t[:4,:4,4:,4:]):
                    raise ExploratoryQuarticError('Non-HHSS tensor unexpectedly contributes to HHSS')
                dirs.append((symbol,part))
                tensors.append(t)
    if matched != {r['symbol'] for r in records}:
        raise ExploratoryQuarticError('Not all candidate records matched')
    required = tuple(sorted(set(required_quartics))) if required_quartics else tuple(sorted(all_symbols))
    missing = set(all_symbols) - set(required)
    if missing:
        raise ExploratoryQuarticError(f'Missing model quartics in running state: {sorted(missing)}')
    return ExploratoryQuarticOverlay(model, required,
        'EXPLORATORY-Matchete-block-order-'+digest[:16],
        tuple(dirs), tuple(tensors), original_baseline, hhss_provider, n)


def run_with_adapter(adapter, ren0, c0, overlay, *, eft_on=True):
    """Call Batch 4 joint DOP853 with an exploratory overlay (not certified)."""
    from Numerical.running.T3ExtendedSimultaneousEFT1 import integrate_extended_eft1
    result = integrate_extended_eft1(adapter, ren0, c0, overlay, eft_on=eft_on)
    # Batch 4's generic flag is *not* physical completeness certification.
    result['complete_quartic_response_used'] = False
    result['exploratory_full_scalar_inventory'] = True
    result['physical_rgbeta_mapping_certified'] = False
    result['quartic_approximation'] = overlay.mapping_assumption
    return result
