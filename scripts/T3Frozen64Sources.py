"""Strict, deterministic IDs and provenance for the frozen full-scalar T3 study.

No wildcard discovery. The manifest never treats historical/reduced results as
full-general-scalar evaluations. No repository mutations on import.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import re

MODELS = (
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
)
BENCHMARKS = {
    'optimised_low': ('optimised', 'small_uv_frozen_5e-9'),
    'optimised_high': ('optimised', 'large_uv_frozen_0p5'),
    'shared_a_m2_low': ('shared_a_m2', 'small_uv_frozen_5e-9'),
    'shared_a_m2_high': ('shared_a_m2', 'large_uv_frozen_0p5'),
}
FIELDS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
          'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
DIRECTIONS = ('T1', 'T2', 'T3')
FORMAT = 't3-full-general-scalar-frozen-64-manifest-v1'
RUN_ID = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$')


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError('T3 frozen64: ' + message)


def tag(model: str) -> str:
    require(model in MODELS, 'invalid model: ' + repr(model))
    return model.replace('+', 'p').replace('-', 'm')


def digest(path: Path) -> str:
    require(path.is_file(), 'missing source: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing source: ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'not a JSON object: ' + str(path))
    return data


def relative(root: Path, path: Path) -> str:
    root, path = root.resolve(), path.resolve()
    require(path.is_relative_to(root), 'source outside project: ' + str(path))
    return path.relative_to(root).as_posix()


def safe_manifest_path(root: Path, value: str) -> Path:
    require(isinstance(value, str) and value and '\\' not in value,
            'source path is not a project-relative POSIX path')
    p = Path(value)
    require(not p.is_absolute() and '..' not in p.parts, 'unsafe path: ' + value)
    return (root / p).resolve()


def expected_keys():
    return {(benchmark, model) for benchmark in BENCHMARKS for model in MODELS}


def validate_manifest(data: dict) -> dict:
    require(data.get('format') == FORMAT, 'wrong manifest format')
    cells = data.get('cells')
    require(isinstance(cells, list) and len(cells) == 64, 'expected exactly 64 cells')
    seen = set()
    for cell in cells:
        require(isinstance(cell, dict), 'non-object manifest cell')
        b, m = cell.get('benchmark'), cell.get('model')
        require((b, m) in expected_keys(), f'unknown cell {b}/{m}')
        require((b, m) not in seen, f'duplicate cell {b}/{m}')
        seen.add((b, m))
        require(cell.get('yukawa_mode') == BENCHMARKS[b][0], 'Yukawa mode mismatch')
        require(cell.get('case') == BENCHMARKS[b][1], 'scalar case mismatch')
        require(cell.get('status') in {'NOT_RUN', 'SCALAR_COMPLETE', 'NEUTRINO_COMPLETE', 'BLOCKED'},
                'unknown numerical status')
    require(seen == expected_keys(), 'missing model/benchmark cells')
    return data


def validate_source(root: Path, path_str: str, sha: str) -> Path:
    p = safe_manifest_path(root, path_str)
    require(digest(p) == sha, 'source hash mismatch: ' + path_str)
    return p
