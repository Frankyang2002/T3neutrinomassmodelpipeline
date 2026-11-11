"""Conditional full-quartic context for nonmixed T3 A_-4/D_-2 coupled IVPs.

The cached Batch3B tensor replaces the old *non-HHSS* potential ansatz only.
Source-resolved HHSS is retained exactly, since the two calculations have no
certified native/RGBeta basis and absolute normalization identification.

All twelve scalar coordinates stay in a fixed real basis over EFT1; the D_-2
saved mass-basis rotation is explicitly identified as candidate-only.
This is not an independent physical quartic or threshold matching certificate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping
import json
import math
import numpy as np

from Numerical.running.T3ScalarQuarticLinearResponse import scalar, symmetry_residual

SHAPE = (12, 12, 12, 12)
MODELS = ('A_-4', 'D_-2')


def digest(path: Path) -> str:
    h = sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda: f.read(1 << 20), b''):
            h.update(part)
    return h.hexdigest()


def relative(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-100))


def _put_hhss(q: np.ndarray, hhss: np.ndarray) -> np.ndarray:
    """Replace each H,H,S,S orbit ONCE, no repeated-orbit accumulation."""
    h = np.asarray(hhss, float)
    if h.shape != (4, 4, 8, 8) or not np.isfinite(h).all():
        raise ValueError('Expected finite 4x4x8x8 HHSS source tensor')
    if not np.allclose(h, h.swapaxes(0, 1), rtol=1e-10, atol=1e-12):
        raise ValueError('Source HHSS Higgs exchange symmetry failed')
    if not np.allclose(h, h.swapaxes(2, 3), rtol=1e-10, atol=1e-12):
        raise ValueError('Source HHSS heavy scalar exchange symmetry failed')
    # Exactly six placements for the two Higgs vs two heavy indices.
    # This writes each placement directly, with no costly permutation loops.
    q[:4, :4, 4:, 4:] = h
    q[:4, 4:, :4, 4:] = h.transpose(0, 2, 1, 3)
    q[:4, 4:, 4:, :4] = h.transpose(0, 2, 3, 1)
    q[4:, :4, :4, 4:] = h.transpose(2, 0, 1, 3)
    q[4:, :4, 4:, :4] = h.transpose(2, 0, 3, 1)
    q[4:, 4:, :4, :4] = h.transpose(2, 3, 0, 1)
    return q


@dataclass
class CandidateQuarticContext:
    model: str
    reference: dict[str, Any]
    directions: tuple[tuple[str, str], ...]
    reference_tensor: np.ndarray
    columns: np.ndarray
    input_hashes: dict[str, str]
    source_directory: str
    number_of_calls: int = 0
    max_source_hhss_relative: float = 0.0
    max_final_symmetry: float = 0.0
    max_cached_tensor_norm: float = 0.0
    sample_points: list[dict[str, float]] = field(default_factory=list)

    @classmethod
    def load(cls, project: Path, response_root: Path, model: str):
        if model not in MODELS:
            raise ValueError('Unsupported model; no generalized cached scalar basis')
        folder = Path(response_root) / model
        report_file, tensor_file = folder / 'report.json', folder / 'quartic_response.npz'
        report = json.loads(report_file.read_text(encoding='utf-8-sig'))
        if report.get('model') != model or report.get('schema') != 't3-batch3b-full-real-quartic-response-v1':
            raise ValueError(f'{model}: wrong Batch3B report/schema')
        if report.get('status') != 'FullReferenceQuarticAffineResponseChecked_PhysicalNormalizationOpen':
            raise ValueError(f'{model}: Batch3B tensor not validated')
        if report.get('output_npz_sha256') != digest(tensor_file):
            raise ValueError(f'{model}: saved tensor hash differs from report')
        hashes = report.get('input_sha256', {})
        registry_file = Path(project) / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
        if hashes.get('rgbeta_registry') != digest(registry_file):
            raise ValueError(f'{model}: RGBeta registry changed since cached tensor build')
        # Batch3B reference coefficients were sourced from this precise native MF snapshot.
        from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
        snap = load_snapshot(Path(project), model)
        if hashes.get('quartic_seed') != digest(snap.quartic_seed):
            raise ValueError(f'{model}: source quartic changed since cached tensor build')
        reference = dict(snap.trajectory.at_fermion_threshold())
        with np.load(tensor_file, allow_pickle=False) as data:
            ref_q = np.asarray(data['reference_tensor'], dtype=float).copy()
            columns = np.asarray(data['response_columns'], dtype=float).copy()
            native = np.asarray(data['native_rotation'], dtype=float)
            mass = np.asarray(data['heavy_mass_rotation'], dtype=float)
        dirs = tuple(tuple(d) for d in report['tensor']['directions'])
        if ref_q.shape != SHAPE or columns.shape != (len(dirs), *SHAPE):
            raise ValueError(f'{model}: response tensor dimensions changed')
        for name, axis in dirs:
            if name not in reference or axis not in ('real','imag'):
                raise ValueError(f'{model}: bad cached coupling coordinate {name}/{axis}')
        for x in (ref_q, columns, native, mass):
            if not np.isfinite(x).all():
                raise ValueError('Nonfinite response or mass-basis matrix')
        if native.shape != (8, 8) or mass.shape != (8, 8) or not np.allclose(mass.T @ mass, np.eye(8), atol=1e-9):
            raise ValueError(f'{model}: invalid fixed scalar rotation')
        if symmetry_residual(ref_q) > 1e-9:
            raise ValueError(f'{model}: tensor exchange symmetry failed')
        return cls(model, reference, dirs, ref_q, columns, dict(hashes), str(folder))

    def full_reference(self, state: Mapping[str, Any]) -> np.ndarray:
        coordinates = []
        for name, axis in self.directions:
            if name not in state or state[name] is None:
                raise ValueError(f'{self.model}: missing running coupling {name}')
            delta = scalar(state[name]) - scalar(self.reference[name])
            coordinates.append(delta.real if axis == 'real' else delta.imag)
        return self.reference_tensor + np.tensordot(np.asarray(coordinates), self.columns, axes=(0,0))

    def tensor(self, state: Mapping[str, Any], hhss: np.ndarray) -> np.ndarray:
        raw = self.full_reference(state)
        existing = np.asarray(raw[:4, :4, 4:, 4:])
        source = np.asarray(hhss, float)
        discrepancy = relative(existing, source)
        q = _put_hhss(raw, source)
        symmetric = symmetry_residual(q)
        if symmetric > 1e-9:
            raise ValueError(f'{self.model}: hybrid quartic violates Bose symmetry: {symmetric}')
        if not np.isfinite(q).all():
            raise ValueError('Non-finite running quartic')
        self.number_of_calls += 1
        self.max_source_hhss_relative = max(self.max_source_hhss_relative, discrepancy)
        self.max_final_symmetry = max(self.max_final_symmetry, symmetric)
        norm = float(np.linalg.norm(q))
        self.max_cached_tensor_norm = max(self.max_cached_tensor_norm, norm)
        if self.number_of_calls in (1, 2, 10, 40, 100, 200, 500):
            self.sample_points.append({'call': self.number_of_calls,
                                       'source_vs_reference_hhss_rel': discrepancy,
                                       'full_quartic_norm': norm})
        return q

    def diagnostics(self) -> dict:
        return {
            'model': self.model,
            'status': 'ExperimentalCachedNonHHSS_OriginalSourceHHSSRetained',
            'number_of_context_calls': self.number_of_calls,
            'max_reference_vs_source_hhss_relative': self.max_source_hhss_relative,
            'max_quartic_symmetry_residual': self.max_final_symmetry,
            'max_full_tensor_norm': self.max_cached_tensor_norm,
            'sample_points': list(self.sample_points),
            'response_source_directory': self.source_directory,
            'quartic_response_directions': [list(x) for x in self.directions],
            'input_hashes': self.input_hashes,
            'full_physical_scalar_potential_certified': False,
            'full_physical_fermion_matching_certified': False,
            'joint_scalar_one_loop_matching_certified': False,
        }


def self_test():
    h = np.zeros((4,4,8,8)); h[1,2,3,4] = h[2,1,3,4] = h[1,2,4,3] = h[2,1,4,3] = 0.7
    q = np.zeros(SHAPE)
    _put_hhss(q, h)
    assert symmetry_residual(q) < 1e-15
    assert np.allclose(q[:4,:4,4:,4:], h)
    baseline = q.copy()
    baseline[5,5,5,5] = 2.0
    assert np.isclose(baseline[5,5,5,5], 2.)
    _put_hhss(baseline, h)
    assert baseline[5,5,5,5] == 2.
    assert np.allclose(baseline[:4,:4,4:,4:], h)
    assert np.array_equal(q[:4,4:,:4,4:], h.transpose(0,2,1,3))
    print('PASS: HHSS six-orbit replacement, non-HHSS preservation, Bose symmetry')
