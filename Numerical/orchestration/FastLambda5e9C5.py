"""Optional numerical C5 assembly for the lambda=5e-9 refit only.

Uses the authoritative hard/direct kernels and the same flavour tensor as
WeinbergFlavorMatching.build_majorana_c5_flavor_matrix. Does not modify repo physics.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as sp

from RGE.matching.FinalWeinbergAdapter import _load_payload, _parse_direct_c12_prefactor
from RGE.matching.MatcheteC5Parsing import parse_matchete_c5
from RGE.matching.WeinbergFlavorMatching import extract_one_generation_loop_kernel


@lru_cache(maxsize=32)
def _kernels(resolved_path: str):
    payload = _load_payload(Path(resolved_path))
    hard = payload.get('renormalized_hard_threshold', {})
    if hard.get('status') != 'Success' or not hard.get('expression'):
        raise ValueError('Missing successful hard C5 threshold')
    direct = payload.get('direct_running_contribution', {})
    if not direct.get('manifestly_symmetric_under_pq', False):
        raise ValueError('Direct C5 contribution not marked symmetric')
    hard_kernel = extract_one_generation_loop_kernel(parse_matchete_c5(hard['expression']))
    direct_prefactor = _parse_direct_c12_prefactor(direct.get('raw_expression', ''))
    running_kernel = sp.simplify(sp.Symbol('hbar') * direct_prefactor / sp.Symbol('MF'))
    return hard_kernel, running_kernel


def evaluate_fast(c5_path, config):
    """Numerically assemble exactly the two authoritative symmetric C5 kernels."""
    model = config['t3']
    masses = np.asarray(model['MF'], dtype=float)
    y1 = np.asarray(model['y1_real'], dtype=float) + 1j*np.asarray(model['y1_imag'], dtype=float)
    y2 = np.asarray(model['y2_real'], dtype=float) + 1j*np.asarray(model['y2_imag'], dtype=float)
    if y1.shape != y2.shape or y1.shape != (3, len(masses)):
        raise ValueError('Inconsistent heavy flavour or Yukawa dimensions')
    hard, running = _kernels(str(Path(c5_path).resolve()))
    ms = model.get('MS')
    if ms is None:
        if not np.isclose(model['MS1'], model['MS2']):
            raise ValueError('Missing grouped scalar threshold MS')
        ms = model['MS1']
    substitutions = {
        sp.Symbol('MS1'): float(model['MS1']),
        sp.Symbol('MS2'): float(model['MS2']),
        sp.Symbol('MS'): float(ms),
        sp.Symbol('lambdaT3'): complex(model['lambdaT3']),
        sp.Symbol('hbar'): float(model.get('hbar', 1/(16*np.pi**2))),
    }
    weights = np.zeros(len(masses), dtype=complex)
    for r, mass in enumerate(masses):
        values = {**substitutions, sp.Symbol('MF'): float(mass)}
        combined = sp.N(hard.subs(values), 18) + sp.N(running.subs(values), 18)
        if combined.free_symbols:
            raise ValueError(f'Unresolved C5 kernel symbols: {combined.free_symbols}')
        weights[r] = complex(combined)
    # Same conjugation and normalization as build_majorana_c5_flavor_matrix.
    c5 = (y1.conj()*weights[np.newaxis, :]) @ y2.conj().T
    c5 += (y2.conj()*weights[np.newaxis, :]) @ y1.conj().T
    if not np.all(np.isfinite(c5)) or not np.allclose(c5, c5.T, atol=1e-14, rtol=1e-10):
        raise ValueError('Non-finite or asymmetric C5 result')
    return c5


def install():
    """Patch this Python process only, including the already-imported adapter alias."""
    import Numerical.running.SMWeinbergStage as sm
    import Numerical.running.FinalC5TrajectoryAdapter as adapter
    sm.evaluate_final_weinberg_json = evaluate_fast
    adapter.evaluate_final_weinberg_json = evaluate_fast
