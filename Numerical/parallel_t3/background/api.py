"""Read-only coordinate-aware consumers for Workstreams A and C.

A may interpolate native RGBeta EFT1 couplings in ln(mu). This is NOT an
RGBeta->Matchete scalar tensor conversion. C may read scalar coordinates at
MS, but B never supplies an invented component mass/quartic tensor.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from .contracts import ContractError, HEAVY, check_finite, hash_file, read_json


@dataclass
class Background:
    manifest: dict
    arrays: dict[str, np.ndarray]

    @classmethod
    def load(cls, folder: Path) -> 'Background':
        folder = Path(folder)
        meta = read_json(folder/'background_manifest.json')
        path = folder/'background.npz'
        if meta.get('schema') != 't3-parallel-B-background-v1':
            raise ContractError('Wrong B background schema')
        if hash_file(path) != meta.get('background_npz_sha256'):
            raise ContractError('B background NPZ differs from manifest SHA256')
        with np.load(path, allow_pickle=False) as data:
            arrays = {n:np.asarray(data[n]) for n in data.files}
        uv_names = set(meta['beta_symbol_order_uv'])
        eft_names = set(meta['beta_symbol_order_eft1'])
        if HEAVY & eft_names or uv_names - HEAVY != eft_names:
            raise ContractError('B EFT1 inventory still contains heavy fermion coordinates')
        expected = ({'mu_uv_gev','mu_eft1_gev','t_uv','t_eft1'} |
                    {f'uv__{n}' for n in uv_names} |
                    {f'eft1__{n}' for n in eft_names})
        if set(arrays) != expected:
            raise ContractError('Unexpected/missing B NPZ arrays')
        for stage, names, scales in [('uv',uv_names,arrays['mu_uv_gev']),
                                      ('eft1',eft_names,arrays['mu_eft1_gev'])]:
            if scales.ndim != 1 or len(scales) < 3 or not np.all(np.diff(scales)<0):
                raise ContractError(f'{stage}: expected descending log-scale samples')
            for n in names:
                v = arrays[f'{stage}__{n}']
                if v.shape[0] != scales.size or v.dtype != np.complex128:
                    raise ContractError(f'{stage}/{n}: array shape or dtype incompatible')
                check_finite(v)
        for n in eft_names:
            if not np.allclose(arrays[f'uv__{n}'][-1],
                               arrays[f'eft1__{n}'][0], atol=1e-9, rtol=1e-11):
                raise ContractError('Mismatch across MF for '+n)
        return cls(meta, arrays)

    def native_eft1_at(self, mu_gev: float) -> dict[str, np.ndarray]:
        """Log-linear interpolation of saved trajectory, NOT re-integration.

        Use the original solver / more samples for higher precision. This
        provides A with a numerical background, not A's Wilson normalization.
        """
        mu = float(mu_gev)
        scales = self.arrays['mu_eft1_gev']
        if not np.isfinite(mu) or not scales[-1] <= mu <= scales[0]:
            raise ContractError('Interpolation scale outside [MS,MF]')
        t = np.log(scales[::-1])
        out = {}
        for n in self.manifest['beta_symbol_order_eft1']:
            v = self.arrays['eft1__'+n][::-1]
            # Interpolate arbitrary 3x3 / scalar complex coordinates.
            flat = v.reshape((v.shape[0], -1))
            x = np.log(mu)
            rv = np.array([np.interp(x,t,flat[:,j].real) for j in range(flat.shape[1])])
            iv = np.array([np.interp(x,t,flat[:,j].imag) for j in range(flat.shape[1])])
            out[n] = (rv+1j*iv).reshape(v.shape[1:])
        return out

    def native_scalar_ms(self) -> dict[str, np.ndarray]:
        keys = self.manifest['beta_symbol_order_eft1']
        return {n:self.arrays['eft1__'+n][-1].copy() for n in keys
                if n.startswith(('lambda','mS','muT3'))}


def load_formal_boundary(folder: Path) -> tuple[dict[str, np.ndarray], dict]:
    folder = Path(folder)
    meta = read_json(folder/'fermion_boundary_manifest.json')
    if meta.get('schema') != 't3-parallel-B-fermion-boundary-v1':
        raise ContractError('Wrong B F-threshold schema')
    if not meta.get('formal_LLsS_tree_available'):
        raise ContractError('This model has no evaluated formal F-tree LLSS boundary')
    path = folder/'fermion_boundary_formal.npz'
    if hash_file(path) != meta.get('fermion_boundary_npz_sha256'):
        raise ContractError('F-threshold NPZ fingerprint mismatch')
    with np.load(path, allow_pickle=False) as data:
        arrays = {n:np.asarray(data[n]) for n in data.files}
    c = arrays.get('C_LLSS_formal')
    if c is None or c.ndim != 6 or c.shape[:4] != (3,3,2,2):
        raise ContractError('Expected formal C_LLSS axes p,q,i,j,a,b')
    if c.shape[-1] != c.shape[-2] or c.dtype != np.complex128:
        raise ContractError('Formal scalar real axis dimensions/dtype wrong')
    check_finite(c)
    return arrays,meta
