"""Explicit T3 fermion-threshold LLSS -> master-Weyl coordinate bridge.

The input is the *formal* polynomial in the project's SourceBoundary layout:
    (flavour_p, flavour_q, L-SU2_i, L-SU2_j, heavy_real_a, heavy_real_b)
This is NOT a certified Matchete GammaCC -> two-component Weyl coefficient.
No scalar-basis identification or spinor factor is silently inserted.

Use an external, independently justified real basis change and signed spinor
factor before treating the result as a physical Wilson coefficient.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Mapping, Sequence

import numpy as np


class SourceBridgeError(ValueError):
    pass


@dataclass(frozen=True)
class FormalSource:
    model: str
    raw: np.ndarray
    mu_f_gev: float
    metadata: dict
    npz_path: Path


def load_formal_source(path: str | Path, model: str | None = None) -> FormalSource:
    """Read exact `write_provisional_source` files, without certifying them."""
    p = Path(path)
    meta_path = p.with_suffix('.json')
    metadata = json.loads(meta_path.read_text(encoding='utf-8-sig'))
    if metadata.get('status') != 'FormalLLSSExported_PhysicalNormalizationUncertified':
        raise SourceBridgeError('Require formal original SourceBoundary export')
    if metadata.get('gamma_cc_weyl_normalization_certified') is not False:
        raise SourceBridgeError('Source normalization status is not the expected unresolved status')
    found = metadata.get('model')
    if model is not None and found != model:
        raise SourceBridgeError(f'Model mismatch: expected {model}, source is {found}')
    with np.load(p, allow_pickle=False) as data:
        if 'C_LLSS_formal_muF' not in data:
            raise SourceBridgeError('Missing formal source tensor')
        raw = np.array(data['C_LLSS_formal_muF'], dtype=np.complex128, copy=True)
    if raw.ndim != 6 or raw.shape[:4] != (3, 3, 2, 2):
        raise SourceBridgeError(f'Invalid source axes: {raw.shape}')
    n = raw.shape[-1]
    if raw.shape[-2] != n or metadata.get('n_heavy_real') != n:
        raise SourceBridgeError('Heavy scalar axis or metadata mismatch')
    if len(metadata.get('source_scalar_real_labels', [])) != n:
        raise SourceBridgeError('Source labels do not cover the heavy real scalar axes')
    mu = float(metadata['mu_f_gev'])
    if not np.isfinite(mu) or mu <= 0 or not np.isfinite(raw).all():
        raise SourceBridgeError('Non-finite source tensor or invalid matching scale')
    _validate_source_orbits(raw)
    return FormalSource(str(found), raw, mu, metadata, p)


def _validate_source_orbits(raw: np.ndarray, *, rtol: float = 1e-10,
                            atol: float = 1e-28) -> None:
    if not np.allclose(raw, raw.transpose(1, 0, 3, 2, 4, 5), rtol=rtol, atol=atol):
        raise SourceBridgeError('Lepton-field interchange symmetry fails')
    if not np.allclose(raw, raw.swapaxes(4, 5), rtol=rtol, atol=atol):
        raise SourceBridgeError('Scalar interchange symmetry fails')


def embed_formal_source(
    raw: np.ndarray, *, model: str, fermion_labels: Sequence[str],
    scalar_labels: Sequence[str], heavy_source_from_target: np.ndarray | None = None,
    spinor_signed_factor: float | None = None,
) -> tuple[np.ndarray, dict]:
    """Return C[nf,nf,ns,ns] in the EXPLICIT target coordinates.

    `heavy_source_from_target[a,A]`: phi_source[a] = R[a,A] phi_target[A].
    If R is omitted, the identity is used only as a *formal coordinate test*.
    A supplied factor does not constitute independent normalization evidence.

    For source coefficient P = sum raw[p,q,i,j,a,b] L[p,i] L[q,j]
    phi_source[a] phi_source[b], map coordinates without combinatoric factors.
    """
    raw = np.asarray(raw, dtype=np.complex128)
    if raw.ndim != 6 or raw.shape[:4] != (3, 3, 2, 2):
        raise SourceBridgeError('Expected source axes (3,3,2,2,n,n)')
    n = raw.shape[-1]
    if raw.shape[-2] != n or not np.isfinite(raw).all():
        raise SourceBridgeError('Invalid/nonfinite source')
    _validate_source_orbits(raw)
    from Numerical.running.GeneralT3OrderedTensor import scalar_basis
    basis = scalar_basis(model)
    ns = int(basis.total_dimension)
    if ns != 4+n or len(scalar_labels) != ns:
        raise SourceBridgeError('Real-scalar dimension differs from source or target basis')
    if len(set(scalar_labels)) != ns or len(set(fermion_labels)) != len(fermion_labels):
        raise SourceBridgeError('Duplicate target basis labels')
    higgs = ('H[1,R]', 'H[1,I]', 'H[2,R]', 'H[2,I]')
    if tuple(scalar_labels[:4]) != higgs:
        raise SourceBridgeError('Target Higgs real-coordinate ordering differs')
    fmap = {label: i for i, label in enumerate(fermion_labels)}
    required = [f'L[{p},{i}]' for p in (1, 2, 3) for i in (1, 2)]
    if any(key not in fmap for key in required):
        raise SourceBridgeError('Missing explicit three-family SU(2) lepton components')
    leptons = [fmap[k] for k in required]
    if len(set(leptons)) != 6:
        raise SourceBridgeError('Repeated target lepton component')

    if heavy_source_from_target is None:
        R = np.eye(n)
        map_status = 'IDENTITY_ASSUMED_FORMAL_ONLY'
    else:
        R = np.asarray(heavy_source_from_target, dtype=float)
        if R.shape != (n, n) or not np.isfinite(R).all():
            raise SourceBridgeError('Heavy scalar basis map must be finite real n x n')
        if not np.allclose(R @ R.T, np.eye(n), atol=1e-10, rtol=1e-10):
            raise SourceBridgeError('Heavy scalar basis map is not orthogonal')
        map_status = 'EXTERNAL_MATRIX_UNCERTIFIED'
    if spinor_signed_factor is None:
        factor = 1.0
        spin_status = 'UNIT_FORMAL_POLYNOMIAL_NOT_WEYL_NORMALIZATION'
    else:
        if isinstance(spinor_signed_factor, (bool, np.bool_)) or not np.isfinite(spinor_signed_factor) or spinor_signed_factor == 0:
            raise SourceBridgeError('Invalid signed spinor factor')
        factor = float(spinor_signed_factor)
        spin_status = 'EXTERNAL_FACTOR_UNCERTIFIED'

    # Change REAL scalar coordinates, not complex multiplet components.
    transformed = np.einsum('pqijab,aA,bB->pqijAB', raw, R, R, optimize=True)
    ordered = transformed.transpose(0, 2, 1, 3, 4, 5).reshape(6, 6, n, n)
    result = np.zeros((len(fermion_labels), len(fermion_labels), ns, ns),
                      dtype=np.complex128)
    heavy = np.arange(4, ns)
    result[np.ix_(leptons, leptons, heavy, heavy)] = factor * ordered
    if not np.allclose(result, result.swapaxes(0, 1), rtol=1e-10, atol=1e-28):
        raise SourceBridgeError('Embedded master fermion pair not symmetric')
    if not np.allclose(result, result.swapaxes(2, 3), rtol=1e-10, atol=1e-28):
        raise SourceBridgeError('Embedded master scalar pair not symmetric')
    details = {
        'model': model, 'dimensions': list(result.shape),
        'fermion_external_label_order': list(fermion_labels),
        'scalar_external_label_order': list(scalar_labels),
        'lepton_indices_zero_based': leptons,
        'heavy_scalar_indices_zero_based': heavy.tolist(),
        'scalar_coordinate_map_status': map_status,
        'spinor_normalization_status': spin_status,
        'claimed_physical_master_boundary': False,
        'note': 'No double-counted symmetry/orbit factor. No certificate supplied by this conversion.',
    }
    return result, details
