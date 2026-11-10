"""Read-only native-RGBeta trajectory and threshold-state interface for T3 Batch 2.

This is a coordinate interface, not a Matchete <-> RGBeta basis conversion.
Every represented running parameter can be retrieved at the exact endpoints or
interpolated in log(mu); no Wilson coefficient is generated or interpreted here.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

FORMAT = 't3-thesis-full-general-batch1-rg-v1'


def _sha(path: Path) -> str:
    h = sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError('JSON root must be an object: ' + str(path))
    return obj


def encode_native(value: Any) -> Any:
    """JSON representation retaining all complex phases and matrix ordering."""
    if isinstance(value, np.ndarray):
        return encode_native(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        z = complex(value)
        if not np.isfinite([z.real, z.imag]).all():
            raise ValueError('Nonfinite complex coupling')
        return {'real': float(z.real), 'imag': float(z.imag)}
    if isinstance(value, (float, int, np.integer, np.floating)):
        v = float(value)
        if not np.isfinite(v):
            raise ValueError('Nonfinite scalar')
        return v
    if isinstance(value, (list, tuple)):
        return [encode_native(x) for x in value]
    if isinstance(value, dict):
        return {str(k): encode_native(v) for k, v in value.items()}
    raise TypeError('Cannot encode native coupling of type ' + type(value).__name__)


def decode_native(value: Any) -> np.ndarray:
    """Decode a scalar or nested matrix from the documented native JSON form."""
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        return np.asarray(complex(float(value['real']), float(value['imag'])))
    if isinstance(value, (int, float)):
        return np.asarray(value, dtype=np.complex128)
    if isinstance(value, list):
        # np.array of decoded scalar values; only rank-2 matrices are expected
        return np.asarray([decode_native(v) for v in value], dtype=np.complex128)
    raise ValueError('Invalid native coupling representation')


@dataclass
class NativeBatch1Trajectory:
    folder: Path
    manifest: dict[str, Any]
    arrays: dict[str, np.ndarray]

    @classmethod
    def load(cls, folder: str | Path) -> 'NativeBatch1Trajectory':
        folder = Path(folder)
        meta = read_json(folder / 'manifest.json')
        if meta.get('format') != FORMAT or meta.get('status') != 'CandidateGeneralRenormalisableRGIntegrated':
            raise ValueError('Not a completed T3 Batch 1 renormalisable trajectory')
        npz = folder / 'trajectories.npz'
        if _sha(npz) != meta.get('trajectory_sha256'):
            raise ValueError('Trajectory checksum mismatch')
        with np.load(npz, allow_pickle=False) as src:
            arrays = {key: np.array(src[key], copy=True) for key in src.files}
        expected = {'uv__mu_gev', 'eft1__mu_gev'}
        for stage in ('uv', 'eft1'):
            names = meta['coupling_order'][stage]
            if len(names) != len(set(names)):
                raise ValueError('Duplicate coupling names')
            expected.update(f'{stage}__{n}' for n in names)
            scales = arrays.get(f'{stage}__mu_gev')
            if scales is None or scales.ndim != 1 or scales.size < 3:
                raise ValueError('Missing scale trajectory: ' + stage)
            if not np.all(np.isfinite(scales)) or not np.all(np.diff(scales) < 0):
                raise ValueError('Scales must be finite and strictly decreasing')
            for name in names:
                arr = arrays.get(f'{stage}__{name}')
                shape = tuple(meta['coupling_shapes'][stage][name])
                if arr is None or arr.dtype != np.complex128 or arr.shape != (len(scales), *shape):
                    raise ValueError('Missing or wrongly shaped coupling: ' + stage + '/' + name)
                if not np.isfinite(arr).all():
                    raise ValueError('Nonfinite coupling: ' + stage + '/' + name)
        if set(arrays) != expected:
            raise ValueError('Extra or missing native trajectory arrays')
        m = meta['scales_gev']
        if not (np.isclose(arrays['uv__mu_gev'][0], m['UV'], atol=1e-6)
                and np.isclose(arrays['uv__mu_gev'][-1], m['F'], atol=1e-6)
                and np.isclose(arrays['eft1__mu_gev'][0], m['F'], atol=1e-6)
                and np.isclose(arrays['eft1__mu_gev'][-1], m['S'], atol=1e-6)):
            raise ValueError('Trajectory scale endpoints inconsistent with manifest')
        uv_names = set(meta['coupling_order']['uv'])
        eft_names = set(meta['coupling_order']['eft1'])
        forbidden = {'MF', 'y1', 'y2', 'y1Conj', 'y2Conj'}
        if forbidden & eft_names or not eft_names <= uv_names:
            raise ValueError('Heavy fermion not decoupled or lost EFT1 running coordinates')
        shifts = meta.get('eft1_boundary_deltas', {})
        if not isinstance(shifts, dict) or not set(shifts) <= eft_names:
            raise ValueError('Unknown EFT1 boundary shifts in the manifest')
        for name in eft_names:
            lhs = arrays[f'uv__{name}'][-1] + (decode_native(shifts[name]) if name in shifts else 0)
            rhs = arrays[f'eft1__{name}'][0]
            if not np.allclose(lhs, rhs, rtol=1e-10, atol=1e-9):
                raise ValueError('Discontinuous fermion boundary after declared shift: ' + name)
        return cls(folder, meta, arrays)

    @property
    def model(self) -> str:
        return str(self.manifest['model'])

    def at(self, stage: str, mu_gev: float) -> dict[str, np.ndarray]:
        """Logarithmic interpolation over a *saved* RG segment.

        For exact matching use the saved threshold inputs, not interpolation.
        Interpolating complex coordinates is a numerical convenience, and is
        not equivalent to re-solving the RGEs at the requested scale.
        """
        if stage not in ('uv', 'eft1'):
            raise ValueError('stage must be uv or eft1')
        scale = float(mu_gev)
        s = self.arrays[f'{stage}__mu_gev']
        if not np.isfinite(scale) or not s[-1] <= scale <= s[0]:
            raise ValueError('Requested scale is outside the saved trajectory')
        x = np.log(s[::-1])
        out: dict[str, np.ndarray] = {}
        for name in self.manifest['coupling_order'][stage]:
            data = self.arrays[f'{stage}__{name}'][::-1]
            shape = data.shape[1:]
            flat = data.reshape((len(s), -1))
            interp = np.array([np.interp(np.log(scale), x, flat[:, i].real)
                               + 1j * np.interp(np.log(scale), x, flat[:, i].imag)
                               for i in range(flat.shape[1])], dtype=np.complex128)
            out[name] = interp.reshape(shape)
        return out

    def at_fermion_threshold(self) -> dict[str, np.ndarray]:
        return {n: self.arrays['uv__' + n][-1].copy()
                for n in self.manifest['coupling_order']['uv']}

    def at_scalar_threshold(self) -> dict[str, np.ndarray]:
        return {n: self.arrays['eft1__' + n][-1].copy()
                for n in self.manifest['coupling_order']['eft1']}

    def fermion_matching_input(self) -> dict[str, Any]:
        if _sha(self.folder / 'fermion_matching_input.json') != self.manifest['fermion_matching_input_sha256']:
            raise ValueError('Fermion matching input checksum mismatch')
        return read_json(self.folder / 'fermion_matching_input.json')

    def scalar_matching_input(self) -> dict[str, Any]:
        if _sha(self.folder / 'scalar_matching_input.json') != self.manifest['scalar_matching_input_sha256']:
            raise ValueError('Scalar matching input checksum mismatch')
        return read_json(self.folder / 'scalar_matching_input.json')


# Alias with a short, stable name for the Batch 2 integration code.
load_batch1 = NativeBatch1Trajectory.load
