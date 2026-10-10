"""Canonical short T3 model labels; representation data is not changed.

Old long IDs and short IDs are both accepted by output discovery.
"""
from __future__ import annotations

import re

CLASS_DIMS = {
    'A': (1, 3, 2),
    'B': (2, 2, 1),
    'C': (2, 2, 3),
    'D': (3, 1, 2),
    'E': (3, 3, 2),
}
LONG = re.compile(r'^T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)$')
SHORT = re.compile(r'^T3-([A-E])-([mp])(\d+)$')


def model_label(ds1: int, ds2: int, df: int, alpha: int) -> str:
    dims = (int(ds1), int(ds2), int(df))
    candidates = [label for label, signature in CLASS_DIMS.items() if dims == signature]
    if len(candidates) != 1:
        raise ValueError(f'Unrecognised or ambiguous T3 representation: {dims}')
    sign = 'm' if int(alpha) < 0 else 'p'
    return f'T3-{candidates[0]}-{sign}{abs(int(alpha))}'


def shorten(name: str) -> str:
    match = SHORT.fullmatch(name)
    if match:
        return name
    match = LONG.fullmatch(name)
    if not match:
        raise ValueError(f'Not a supported T3 model key: {name}')
    ds1, ds2, df = map(int, match.group(1, 2, 3))
    alpha = int(match.group(5)) * (-1 if match.group(4) == 'm' else 1)
    return model_label(ds1, ds2, df, alpha)


def long_name(name: str) -> str:
    match = LONG.fullmatch(name)
    if match:
        shorten(name)  # validate dimensions
        return name
    match = SHORT.fullmatch(name)
    if not match:
        raise ValueError(f'Not a supported T3 model key: {name}')
    label, sign, magnitude = match.groups()
    ds1, ds2, df = CLASS_DIMS[label]
    return f'T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_{sign}{int(magnitude)}'


def is_model(name: str) -> bool:
    try:
        shorten(name)
        return True
    except ValueError:
        return False
