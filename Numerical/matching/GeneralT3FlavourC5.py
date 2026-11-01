"""Exploratory full-general Matchete C5 flavour assembly for NON-mixed T3 models.

Uses the existing Lagrangian/EvaluateGeneralC5Channels.wl to evaluate every
quadratic Yukawa monomial of the *full-general* one-generation Matchete C5.
It does not infer missing matching terms or substitute restricted C5 kernels.
The UV Yukawa matrices are transformed into the existing heavy-fermion mass
basis; channel flavour tensors are summed before final Majorana symmetrisation.

This is one-scale provisional C5, not fermion-first sequential EFT matching.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Mapping

import numpy as np

from Numerical.matching.GeneralT3MultiModelNeutrinos import (
    COUPLING_PATTERN, YUK, kernel_parameters,
)

ROOT = Path(__file__).resolve().parents[2]
EVALUATOR_WL = ROOT / 'Lagrangian' / 'EvaluateGeneralC5Channels.wl'
MIXED = frozenset(('B_-1', 'C_-1'))
MARKER = re.compile(r'^T3YUK(y1Conj|y2Conj|y1|y2)(Bar|Plain)$')


def _require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError('Full-general flavour C5: ' + reason)


def _raw_matrix(raw: Any, name: str) -> np.ndarray:
    value = np.asarray(raw, dtype=complex)
    _require(value.shape == (3, 3) and np.isfinite(value.real).all()
             and np.isfinite(value.imag).all(), name + ' must be finite 3x3')
    return value


def prepare_heavy_flavour(uv_state: Any, required: set[str]):
    """Reuse the original Takagi/SVD transformation for all allowed channels.

    For an extra conjugate-scalar Yukawa the same heavy chiral leg as its
    parent y1 or y2 is assumed. This extra-vertex assignment is provisional.
    """
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data

    _require({'y1', 'y2'} <= required, 'source must contain both topology Yukawas')
    basis = _heavy_mass_basis_data(
        uv_state, offdiagonal_atol=1e-7, imaginary_atol=1e-9,
        takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
    )
    masses = np.asarray(basis.masses_gev, dtype=float)
    _require(masses.shape == (3,) and np.isfinite(masses).all() and (masses > 0).all(),
             'invalid heavy-fermion eigenmasses')
    yukawas = {'y1': _raw_matrix(basis.y1, 'y1'),
               'y2': _raw_matrix(basis.y2, 'y2')}
    extra = uv_state.general_couplings
    for key in sorted(required - {'y1', 'y2'}):
        _require(key in extra, 'matching requires missing UV Yukawa ' + key)
        raw = _raw_matrix(extra[key], key)
        if uv_state.representation.majorana_fermion:
            rotation = basis.rotation
        elif key == 'y1Conj':
            rotation = basis.rotation.conj()
        else:
            _require(basis.right_rotation is not None,
                     'missing vectorlike right rotation')
            rotation = basis.right_rotation
        yukawas[key] = raw @ rotation
    return masses, yukawas, {'heavy_basis_residual': float(max(
        basis.takagi_residual, basis.basis_residual)),
        'extra_yukawa_rotation': 'same-heavy-leg-as-parent-assumption'
        if required - {'y1', 'y2'} else 'not-applicable'}


def flavour_c5_from_channels(rows: list[dict], yukawas: Mapping[str, np.ndarray],
                             channels: list[str] | None = None) -> np.ndarray:
    """Keep different Yukawa pairings and complex conjugations independent."""
    _require(isinstance(rows, list) and len(rows) > 0, 'no Matchete channels')
    result = np.zeros((3, 3), dtype=complex)
    seen: set[tuple[int, str]] = set()
    present = set()
    for row in rows:
        _require(isinstance(row, dict) and row.get('ok') is True,
                 'unevaluated or malformed Matchete channel')
        idx = int(row['index'])
        _require(idx in (0, 1, 2), 'heavy flavour index out of range')
        channel = str(row['channel'])
        _require((idx, channel) not in seen, 'duplicated heavy-flavour channel')
        seen.add((idx, channel)); present.add(channel)
        factors = row.get('factors')
        _require(isinstance(factors, list) and len(factors) == 2,
                 'channel does not contain exactly two Yukawa factors')
        vectors = []
        for factor in factors:
            match = MARKER.fullmatch(str(factor))
            _require(match is not None, 'unrecognised Yukawa marker ' + str(factor))
            name, conjugation = match.groups()
            _require(name in yukawas, 'missing Yukawa matrix ' + name)
            matrix = _raw_matrix(yukawas[name], name)
            v = matrix[:, idx]
            vectors.append(np.conjugate(v) if conjugation == 'Bar' else v)
        coefficient = complex(float(row['real']), float(row['imag']))
        _require(math.isfinite(coefficient.real) and math.isfinite(coefficient.imag),
                 'non-finite numeric kernel')
        result += coefficient * np.outer(*vectors)
    if channels is not None:
        _require(set(channels) == present, 'Matchete channel inventory not fully evaluated')
        expected = {(r, channel) for r in range(3) for channel in channels}
        _require(seen == expected, 'Matchete channel/heavy-flavour grid incomplete')
    result = (result + result.T) / 2  # Majorana, NOT Hermitian conjugation.
    _require(np.isfinite(result.real).all() and np.isfinite(result.imag).all(),
             'nonfinite assembled C5')
    return result


@dataclass
class GeneralT3FlavourC5:
    """Callable with the same evaluate(uv_state, eft_state) API as old C5 kernels."""

    model: str
    coefficient_path: Path
    wolframscript: str = 'wolframscript'
    hbar_value: float = 1 / (16 * math.pi**2)
    timeout_sec: int = 240

    def __post_init__(self):
        self.coefficient_path = Path(self.coefficient_path)
        _require(self.model not in MIXED, 'mixed B/C requires its mass-eigenstate matching')
        _require(self.coefficient_path.is_file(), 'missing full-general Matchete expression ' +
                 str(self.coefficient_path))
        _require(EVALUATOR_WL.is_file(), 'missing Wolfram channel evaluator')
        _require(self.hbar_value > 0 and math.isfinite(self.hbar_value), 'invalid hbar')
        source = self.coefficient_path.read_text(encoding='utf-8-sig')
        self.symbols = set(COUPLING_PATTERN.findall(source))
        self.required_yukawas = self.symbols & YUK
        _require({'y1', 'y2'} <= self.required_yukawas,
                 'full-general matched coefficient lacks y1 or y2')
        self.last_diagnostics: dict = {}

    def evaluate(self, uv_state: Any, eft_state: Any) -> np.ndarray:
        from Numerical.orchestration.GeneralRunningBenchmark import _state_values

        masses, yukawas, basis_diagnostics = prepare_heavy_flavour(
            uv_state, self.required_yukawas)
        uv = _state_values(uv_state)
        eft = _state_values(eft_state)
        needed_splitting = any(s in self.symbols for s in
                               ('T31PlusMass','T31MinusMass',
                                'T32PlusMass','T32MinusMass','T3NeutralPhase'))
        kappa = None
        if needed_splitting:
            source = ('muT3GeneralS2Sq' if self.model.endswith('-2')
                      else 'muT3GeneralS1Sq')
            _require(source in eft, 'no running scalar splitting for ' + source)
            kappa = complex(eft[source]['real'], eft[source]['imag']) if isinstance(
                eft[source], dict) else complex(eft[source])
            if self.model in {'A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'}:
                _require(abs(kappa) > 0.0, 'self-bilinear mass basis requires nonzero mu2')
                _require(kappa.real > 0.0 and
                         abs(kappa.imag) <= 1e-10 * max(1.0, abs(kappa)),
                         'current real self-bilinear Matchete rotation only supports '
                         'positive real mu2; complex-phase rotation not yet derived')
        points = []
        mappings: list[str] = []
        for i, mf in enumerate(masses):
            scalars, notes = kernel_parameters(
                self.symbols, self.model, uv, eft, mf,
                kappa_override=kappa, strict_mapping=False)
            points.append({'index': i, 'couplings': scalars})
            mappings.extend(notes)
        with tempfile.TemporaryDirectory(prefix='t3_general_c5_') as temp:
            inputs = Path(temp) / 'points.json'
            outputs = Path(temp) / 'channels.json'
            inputs.write_text(json.dumps({'points': points,
                                          'hbar_value': self.hbar_value}),
                              encoding='utf-8')
            command = [self.wolframscript, '-file', str(EVALUATOR_WL),
                       str(self.coefficient_path), str(inputs), str(outputs)]
            proc = subprocess.run(command, cwd=ROOT, capture_output=True,
                                  text=True, timeout=self.timeout_sec, check=False)
            if proc.returncode or not outputs.is_file():
                raise RuntimeError('Full-general Matchete C5 channel evaluation failed '
                                   f'({proc.returncode}): ' +
                                   (proc.stdout + proc.stderr)[-3000:])
            payload = json.loads(outputs.read_text(encoding='utf-8-sig'))
        raw_rows = payload.get('results')
        channels = payload.get('channels')
        _require(isinstance(channels, list) and bool(channels),
                 'no channel inventory from Wolfram')
        c5 = flavour_c5_from_channels(raw_rows, yukawas, channels)
        self.last_diagnostics = {
            'matching_kind': 'full-general-Matchete-one-scale-Yukawa-channels',
            'self_bilinear_real_mass_basis_candidate': self.model in {
                'A_-2','A_+0','D_-2','D_+0','E_-2','E_+0'},
            'matched_yukawas': sorted(self.required_yukawas),
            'n_yukawa_channels': len(channels),
            'channels': channels,
            'quartic_basis_mapping_assumptions': sorted(set(mappings)),
            'hbar_value': self.hbar_value,
            **basis_diagnostics,
            'sequential_eft_wilson_matching_performed': False,
            'physical_general_c5_certified': False,
        }
        return c5
