"""Read-only, validated access to saved full-comparison running diagnostics.

Both historical long and new short model folder names are supported.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import re
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.orchestration.ModelNames import is_model, shorten

# Retained for consumers that import the historical pattern.
MODEL_PATTERN = re.compile(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_[mp]\d+')

SCENARIOS = (
    'smallY_smallL', 'smallY_largeL', 'largeY_smallL', 'largeY_largeL',
)

QUANTITIES: dict[str, tuple[str, str, str, str]] = {
    'y1': ('uv_running', 'y1_frobenius_norm', 'scalar', r'$\Vert y_1\Vert_F$'),
    'y2': ('uv_running', 'y2_frobenius_norm', 'scalar', r'$\Vert y_2\Vert_F$'),
    'lambda': ('uv_running', 'lambdaT3_abs', 'scalar', r'$|\lambda_{T3}|$'),
    'direct': ('intermediate_direct_weinberg', 'delta_c5_abs', 'matrix_norm', r'$\Vert\Delta C_5^{\rm direct}\Vert_F$ [GeV$^{-1}$]'),
    'c5': ('final_running', 'c5_abs', 'matrix_norm', r'$\Vert C_5\Vert_F$ [GeV$^{-1}$]'),
    'dm21': ('final_running', 'delta_m21_sq_ev2', 'scalar', r'$\Delta m_{21}^2$ [eV$^2$]'),
    'dm3l': ('final_running', 'delta_m3l_sq_ev2', 'scalar_abs', r'$|\Delta m_{3\ell}^2|$ [eV$^2$]'),
    'm1': ('final_running', 'masses_ev', 'vector0', r'$m_1$ [eV]'),
    'm2': ('final_running', 'masses_ev', 'vector1', r'$m_2$ [eV]'),
    'm3': ('final_running', 'masses_ev', 'vector2', r'$m_3$ [eV]'),
}


@dataclass(frozen=True)
class SavedRun:
    scenario: str
    model_key: str
    path: Path
    payload: dict[str, Any]


def load_diagnostics(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(payload, dict) or payload.get('status') != 'Success':
        raise ValueError(f'Not a successful diagnostics object: {path}')
    for key in ('scales_gev', 'uv_running', 'intermediate_direct_weinberg', 'final_running'):
        if not isinstance(payload.get(key), dict):
            raise ValueError(f'Missing diagnostics block {key}: {path}')
    return payload


def extract_trajectory(payload: dict[str, Any], quantity: str) -> tuple[np.ndarray, np.ndarray]:
    block_name, data_key, selector, _ = QUANTITIES[quantity]
    block = payload[block_name]
    mu = np.asarray(block['mu_gev'], dtype=float)
    values = np.asarray(block[data_key], dtype=float)
    if mu.ndim != 1 or mu.size < 2 or np.any(~np.isfinite(mu)) or np.any(mu <= 0):
        raise ValueError(f'Invalid scale array for {quantity}')
    if selector == 'matrix_norm':
        if values.shape != (mu.size, 3, 3):
            raise ValueError(f'Expected (n,3,3) matrix data for {quantity}')
        values = np.linalg.norm(values.reshape(mu.size, 9), axis=1)
    elif selector.startswith('vector'):
        if values.shape != (mu.size, 3):
            raise ValueError(f'Expected (n,3) vector data for {quantity}')
        values = values[:, int(selector[-1])]
    elif values.shape != (mu.size,):
        raise ValueError(f'Expected length-n scalar data for {quantity}')
    if selector == 'scalar_abs':
        values = np.abs(values)
    if np.any(~np.isfinite(values)):
        raise ValueError(f'Nonfinite values for {quantity}')
    return mu, values


def discover_saved_runs(root: Path) -> dict[str, dict[str, SavedRun]]:
    """Discover both long and short model directories; expose short report keys."""
    root = Path(root)
    found: dict[str, dict[str, SavedRun]] = {}
    for scenario in SCENARIOS:
        scenario_dir = root / scenario
        if not scenario_dir.is_dir():
            continue
        for model_dir in sorted(scenario_dir.iterdir()):
            if not model_dir.is_dir() or not is_model(model_dir.name):
                continue
            key = shorten(model_dir.name)
            paths = sorted(model_dir.glob('T3_*/data/running_diagnostics.json'))
            if len(paths) > 1:
                raise ValueError(f'Multiple physical run diagnostics: {model_dir}')
            if not paths:
                continue
            if scenario in found.get(key, {}):
                raise ValueError(f'Duplicate long/short folders for {scenario}/{key}')
            path = paths[0]
            found.setdefault(key, {})[scenario] = SavedRun(
                scenario, key, path, load_diagnostics(path),
            )
    return found
