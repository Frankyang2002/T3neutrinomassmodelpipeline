"""Running heavy-fermion mass crossings on an UNDECOUPLED T3 UV trajectory.

Uses the same mass definition as FermionThresholdBoundary: ascending singular
values of the complete running complex 3x3 MF (Takagi masses for symmetric
Majorana MF; bi-unitary singular masses for vectorlike MF). This locates
candidate crossings mu = M_i(mu) only. It neither rotates Yukawas nor performs
matching, removes a heavy field, or evolves the intervening EFTs.

The caller supplies the existing UV ODE dense solution. Each scan and root
re-diagonalizes MF itself; no linear interpolation of heavy masses or assumed
fixed UV mass hierarchy is used. Labels fermion:0,1,2 mean ascending MASS
RANK at each scale, not a conserved flavour/eigenvector identity.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class FermionCrossing:
    particle: str
    kind: str
    scale_gev: float
    mass_gev: float
    relative_mass_scale_residual: float
    mass_order_at_crossing_gev: tuple[float, ...]
    method: str = 'dense full-UV RGBeta ODE; re-SVD of running MF; brentq in ln(mu)'
    status: str = 'candidate_crossing_matching_required_not_applied'


def singular_masses(mf: np.ndarray, *, majorana: bool, gap_rtol: float = 1e-7) -> np.ndarray:
    """Return sorted positive singular masses, rejecting ambiguous degeneracy.

    For Majorana MF require complex symmetry as in T3UVState.validated().
    SVD singular values are the Takagi masses, independent of basis phases.
    No Takagi rotation is needed for this threshold-only diagnostic.
    """
    matrix = np.asarray(mf, dtype=complex)
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
        raise ValueError('MF must be finite and 3x3')
    if majorana and not np.allclose(matrix, matrix.T, rtol=1e-10, atol=1e-14):
        raise ValueError('Majorana running MF is not complex symmetric')
    masses = np.sort(np.asarray(np.linalg.svd(matrix, compute_uv=False), dtype=float))
    if not np.all(np.isfinite(masses)) or masses[0] <= 0:
        raise ValueError('Running heavy MF has a nonpositive/nonfinite singular mass')
    if np.any(np.diff(masses) <= gap_rtol * np.maximum(masses[1:], 1.0)):
        raise ValueError('Heavy singular masses degenerate: mass-rank tracking is ambiguous')
    return masses


def mass_matrix_from_vector(vector: np.ndarray, layout: object, *,
                            trailing_components: int = 0) -> np.ndarray:
    """Read MF from a complete UV ODE vector in the existing real/imag layout.

    Most models have a vector of length ``layout.size``. The eight ordinary
    T3 models with a holomorphic scalar bilinear have precisely two additional
    real ODE coordinates (Re(mu2), Im(mu2)) *after* the RGBeta state blocks.
    They are NOT part of StateVectorLayout. The caller must explicitly specify
    ``trailing_components=2`` for such a UV ODE; arbitrary trailing entries
    are rejected rather than silently ignored.
    """
    block = layout.block('MF')
    n = int(np.prod(block.shape))
    vec = np.asarray(vector, dtype=float)
    if trailing_components not in (0, 2) or \
       block.kind != 'complex_matrix' or block.shape != (3, 3) or \
       block.stop-block.start != 2*n or \
       vec.shape != (layout.size + trailing_components,) or \
       not np.isfinite(vec).all():
        raise ValueError('Invalid packed RGBeta MF/layout')
    return (vec[block.start:block.start+n].reshape((3,3))
            + 1j*vec[block.start+n:block.stop].reshape((3,3)))


def locate_fermion_crossings(matrix_at_logmu: Callable[[float], np.ndarray],
                              mu_start_gev: float, mu_stop_gev: float,
                              *, majorana: bool, scan_points: int = 65) -> dict:
    """Find descending mu=M_i(mu) crossings on ONE still-full UV solution.

    No extrapolation outside [stop,start]. An unbracketed branch is reported
    as such, not silently assigned an endpoint matching scale. An inverse or
    repeated crossing fails closed. For distinct masses the rank order remains
    usable; physical eigenvector tracking across degeneracy is not attempted.
    """
    hi, lo = float(mu_start_gev), float(mu_stop_gev)
    if not (math.isfinite(hi) and math.isfinite(lo) and hi > lo > 0):
        raise ValueError('UV bounds must be finite, positive, descending')
    if not isinstance(scan_points, int) or scan_points < 9 or scan_points > 2001:
        raise ValueError('scan_points must be an integer in [9,2001]')
    logs = np.linspace(math.log(hi), math.log(lo), scan_points)

    def masses(t: float) -> np.ndarray:
        return singular_masses(matrix_at_logmu(float(t)), majorana=majorana)

    samples = np.asarray([masses(float(t)) for t in logs])
    ratios = samples/np.exp(logs)[:, None]-1.0
    events = []
    unresolved = []
    for k in range(3):
        roots = []
        for j in range(len(logs)-1):
            g0, g1 = float(ratios[j,k]), float(ratios[j+1,k])
            if g0 > 0 and g1 < 0:
                raise ValueError(f'fermion:{k}: inverse crossing on undecoupled UV flow')
            if g0 <= 0 <= g1:
                if g0 == 0:
                    root = float(logs[j])
                elif g1 == 0:
                    root = float(logs[j+1])
                else:
                    root = float(brentq(lambda t: masses(t)[k]/math.exp(t)-1,
                                        float(logs[j+1]), float(logs[j]), xtol=1e-13))
                if not any(abs(root-old) < 1e-10 for old in roots):
                    roots.append(root)
        if len(roots) > 1:
            raise ValueError(f'fermion:{k}: multiple crossings in still-full UV')
        if roots:
            t=roots[0]
            eigen=masses(t)
            scale=math.exp(t)
            err=abs(float(eigen[k]/scale-1))
            if err > 1e-9:
                raise RuntimeError(f'fermion:{k}: root residual {err}')
            events.append(FermionCrossing(
                f'fermion:{k}', 'fermion', scale, float(eigen[k]), err,
                tuple(map(float,eigen))))
        else:
            start_ratio=float(ratios[0,k])
            end_ratio=float(ratios[-1,k])
            if end_ratio < 0:
                reason='below_saved_UV_stop_undetermined'
            elif start_ratio > 0:
                reason='above_saved_UV_start_undetermined'
            else:
                reason='unbracketed_no_crossing'
            unresolved.append({'particle':f'fermion:{k}', 'status':reason,
                               'mass_at_UV_start_gev':float(samples[0,k]),
                               'mass_at_grouped_F_gev':float(samples[-1,k]),
                               'mass_over_scale_minus_one_at_grouped_F':end_ratio})
    return {
        'crossings':[vars(v) for v in sorted(events,key=lambda c:-c.scale_gev)],
        'unresolved_crossings':unresolved,
        'fermion_mass_method':('Majorana complex-symmetric MF -> singular/Takagi masses'
                               if majorana else 'vectorlike MF -> biunitary singular masses'),
        'scan_points':scan_points,
        'mu_uv_start_gev':hi,
        'mu_grouped_fermion_endpoint_gev':lo,
        'masses_at_uv_start_gev':samples[0].tolist(),
        'masses_at_grouped_F_gev':samples[-1].tolist(),
        'all_three_crossings_bracketed':len(events)==3,
        'sequential_heavy_field_decoupling_performed':False,
        'physical_matching_applied':False,
    }
