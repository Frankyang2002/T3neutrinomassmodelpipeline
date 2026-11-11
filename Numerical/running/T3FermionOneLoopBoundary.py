"""Worker A: safe, separately inspectable EFT1 fermion matching boundary.

This module does NOT calculate UV-minus-EFT1 matching. A finite one-loop
coefficient must be supplied from an independently renormalised, soft-subtracted
source. Raw Matchete Laurent-regular expressions are explicitly forbidden.

No production solver is changed by importing this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import numpy as np


class UncertifiedMatchingError(ValueError):
    """Finite one-loop boundary has not been established."""


@dataclass(frozen=True)
class FermionBoundaryParts:
    model: str
    matching_scale_gev: float
    tree: np.ndarray
    finite_one_loop: np.ndarray | None
    provenance: dict[str, Any]

    @property
    def complete_one_loop(self) -> bool:
        return self.finite_one_loop is not None

    @property
    def total(self) -> np.ndarray:
        if self.finite_one_loop is None:
            raise UncertifiedMatchingError(
                'One-loop boundary unavailable: do not substitute zero or raw Laurent regular terms.')
        return self.tree + self.finite_one_loop


def _tensor(raw: Any, name: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    arr = np.asarray(raw, dtype=np.complex128)
    if arr.ndim != 6 or arr.shape[:4] != (3, 3, 2, 2):
        raise ValueError(f'{name}: required axis order (3,3,2,2,n_scalar,n_scalar)')
    if arr.shape[4] != arr.shape[5] or arr.shape[4] <= 0:
        raise ValueError(f'{name}: scalar dimensions must be equal and positive')
    if shape is not None and arr.shape != shape:
        raise ValueError(f'{name}: mismatching tree/loop tensor shapes')
    if not np.all(np.isfinite(arr)):
        raise ValueError(f'{name}: contains non-finite numbers')
    if not np.allclose(arr, arr.transpose(1, 0, 3, 2, 4, 5), rtol=1e-12, atol=1e-30):
        raise ValueError(f'{name}: failed simultaneous lepton exchange')
    if not np.allclose(arr, arr.swapaxes(4, 5), rtol=1e-12, atol=1e-30):
        raise ValueError(f'{name}: failed scalar exchange')
    frozen = arr.copy()
    frozen.setflags(write=False)
    return frozen


def assemble_fermion_boundary(
    model: str,
    tree: Any,
    *,
    matching_scale_gev: float = 1e5,
    finite_one_loop: Any | None = None,
    provenance: Mapping[str, Any] | None = None,
) -> FermionBoundaryParts:
    """Store tree C5 and separately certified finite 1L C5, never add raw terms.

    This is the same six-axis convention as GeneralT3DynamicFermionWilson:
    [lepton generation, lepton generation, SU2, SU2, real scalar, real scalar].
    No field-basis conversion or GammaCC normalisation is performed here.

    Certification metadata is a *caller assertion*, not independent proof.
    Requiring it prevents accidental promotion of known raw results, but a
    human/source audit remains essential before production use.
    """
    if not model or not isinstance(model, str):
        raise ValueError('model must be a non-empty label')
    mu = float(matching_scale_gev)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError('matching scale must be finite and positive')
    a = _tensor(tree, 'tree')
    info = dict(provenance or {})
    b = None
    if finite_one_loop is not None:
        required = {
            'status': 'FiniteRenormalisedUVMinusEFT1',
            'all_heavy_fermions_integrated_together': True,
            'both_scalar_multiplets_retained': True,
            'eft1_soft_subtraction_complete': True,
            'counterterms_included': True,
            'operator_normalization_verified': True,
            'restricted_y1_y2_only': True,
            'conjugate_scalar_yukawas_absent': True,
        }
        problems = [key for key, expected in required.items() if info.get(key) != expected]
        if not isinstance(info.get('source_sha256'), str) or len(info['source_sha256']) != 64:
            problems.append('source_sha256')
        if problems:
            raise UncertifiedMatchingError('Finite matching provenance missing/invalid: '
                                            + ', '.join(problems))
        if info.get('source_type') in {'raw_matchete_laurent', 'unsubtracted_regular',
                                       'candidate_hard_source'}:
            raise UncertifiedMatchingError('Raw source is not a finite matching coefficient')
        b = _tensor(finite_one_loop, 'finite_one_loop', a.shape)
    return FermionBoundaryParts(model, mu, a, b, info)
